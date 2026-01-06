import openmm as mm
from openmm import app, unit
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.units.openmm import to_openmm
from openff.interchange import Interchange
from rdkit import Chem
import numpy as np
import os
import time

os.environ["OPENFF_FORCE_FIELD_PARALLEL"]="14"

n_carp=50           
n_water=50       
box_size=10.0       
timestep=4.0       

smiles="c12-c3c4c(c5cc(cc(-c6cccc(-c7c8c(c9cc(cc(-c(ccc1)n2)c9[nH]8)c1ccccc1)cc(c7)c1ccccc1)n6)c5[nH]4)c1ccccc1)cc(c3)c1ccccc1"

print(f"{n_carp} assembler + {n_water} water in {box_size}nm box")


carp_mol=Molecule.from_pdb_and_smiles('carpyridine.pdb', smiles,allow_undefined_stereo=True)
carp_coords=to_openmm(carp_mol.conformers[0])

water=Chem.AddHs(Chem.MolFromSmiles('O'))
water_mol=Molecule.from_rdkit(water)

off_topology=Topology.from_molecules([carp_mol] * n_carp + [water_mol] * n_water)

def rotate():
    u=np.random.uniform(0, 1, 3)
    q=np.array([np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])])
    r=np.array([[1 - 2*q[1]**2 - 2*q[2]**2, 2*q[0]*q[1] - 2*q[3]*q[2], 2*q[0]*q[2] + 2*q[3]*q[1]],[2*q[0]*q[1] + 2*q[3]*q[2], 1 - 2*q[0]**2 - 2*q[2]**2, 2*q[1]*q[2] - 2*q[3]*q[0]],[2*q[0]*q[2] - 2*q[3]*q[1], 2*q[1]*q[2] + 2*q[3]*q[0], 1 - 2*q[0]**2 - 2*q[1]**2]])
    return r

positions=[]
np.random.seed(42)
limit=box_size 
placed_centers=[]

carp_coords_nm=carp_coords.value_in_unit(unit.nanometer)
com=np.mean(carp_coords_nm, axis=0)
relative_coords=carp_coords_nm - com

print("placing assemblers")
for i in range(n_carp):
    while True:
        offset=np.random.uniform(0.5, limit - 0.5, 3)
        too_close=False
        for center in placed_centers:
            dist=np.linalg.norm(offset - center)
            if dist < 1.5:
                too_close=True
                break
        if not too_close:
            placed_centers.append(offset)
            rot_matrix=rotate()
            rotated_coords=np.dot(relative_coords, rot_matrix.T)
            final_coords=(rotated_coords + offset) * unit.nanometer
            for atom_pos in final_coords:
                positions.append(atom_pos)
            break

print(f"placing water")

n_buffer=int(n_water * 1.2)
all_offsets=np.random.uniform(0.1, limit - 0.1, (n_buffer, 3))
centers_array=np.array(placed_centers)
diff=all_offsets[:, np.newaxis, :] - centers_array[np.newaxis, :, :]
dist_sq=np.sum(diff**2, axis=2)
min_dists=np.sqrt(np.min(dist_sq, axis=1))
valid_offsets=all_offsets[min_dists > 0.6]

h1_vec=np.array([0.0958, 0, 0])
h2_vec=np.array([-0.024, 0.092, 0])

for i in range(n_water):
    off=valid_offsets[i]
    positions.append(off * unit.nanometer)
    positions.append((off + h1_vec) * unit.nanometer)
    positions.append((off + h2_vec) * unit.nanometer)

print(f"waters placed")

off_ff=ForceField('openff-2.1.0.offxml')

print("assembler topology")
carp_topo=Topology.from_molecules([carp_mol] * n_carp)
carp_interchange=Interchange.from_smirnoff(force_field=off_ff, topology=carp_topo)

print("water topology")
water_topo=Topology.from_molecules([water_mol] * n_water)
water_interchange=Interchange.from_smirnoff(force_field=off_ff, topology=water_topo)

print("combined topology")
interchange=carp_interchange + water_interchange

print("total system")
system=interchange.to_openmm(hydrogen_mass=4.0, combine_nonbonded_forces=True)
matched_topology=interchange.to_openmm_topology()

for force in system.getForces():
    if isinstance(force, mm.NonbondedForce):
        force.setNonbondedMethod(mm.NonbondedForce.PME)
        force.setCutoffDistance(0.9 * unit.nanometer)
        force.setEwaldErrorTolerance(0.001)
        force.setReactionFieldDielectric(78.3) 

integrator=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, timestep*unit.femtoseconds)
platform=mm.Platform.getPlatformByName('CPU')
properties={'Threads': '14'} 

simulation=app.Simulation(matched_topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)
size_nm=box_size * unit.nanometer
simulation.context.setPeriodicBoxVectors([size_nm, 0, 0], [0, size_nm, 0], [0, 0, size_nm])

print("minimising energy")
simulation.minimizeEnergy(tolerance=1.0*unit.kilojoule_per_mole / unit.nanometer, maxIterations=2000)

with open('topology.pdb', 'w') as f:
    state=simulation.context.getState(getPositions=True)
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

runs=1000
run_ns=0.1
steps=int((run_ns * 1000000) / timestep)

direc="checkpoints"
if not os.path.exists(direc):
    os.makedirs(direc)

simulation.reporters.append(app.DCDReporter('assembly.dcd', 50000))
simulation.reporters.append(app.StateDataReporter('assembly.log', 50000, step=True, potentialEnergy=True, temperature=True, density=True))

print(f"starting sim")
start=time.time()
for i in range(1, runs + 1):
    segment_start=time.time()
    print(f"segment {i}/{runs}...")
    simulation.step(steps)
    path=os.path.join(direc, f'checkpoint_{i}.chk')
    state_path=os.path.join(direc, f'state_{i}.xml')
    simulation.saveCheckpoint(path)
    simulation.saveState(state_path)
    segment_end=time.time()
    segment_duration=(segment_end - segment_start) / 60.0 
    print(f"segment {i} completed in {segment_duration} minutes")

end=time.time()
seconds= end-start
hours=seconds/3600
mins=hours/60

print(f"done")
print(f"time : {int(total_hours)} hrs, {total_mins:.1f} mins")

