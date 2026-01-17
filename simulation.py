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
n_water=1000      
box_size=10.0       
timestep=4.0       

smiles_dict={'2H-CAR-H':'c12c3c(-c4nc(-c5c6c(c7c(c(-c8cccc(-c9c(c1ccc9)[nH]3)n8)ccc7)[nH]6)ccc5)ccc4)ccc2',
             '2H-CAR-C1':'c12c3c(-c4nc(-c5c6c(c7c(c(-c8cccc(-c9c(c1cc(c9)C)[nH]3)n8)cc(C)c7)[nH]6)cc(C)c5)ccc4)cc(c2)C',
             '2H-CAR-C4':'c12c3c(-c4nc(-c5c6c(c7c(c(-c8cccc(-c9c(c1cc(c9)CCCC)[nH]3)n8)cc(CCCC)c7)[nH]6)cc(CCCC)c5)ccc4)cc(c2)CCCC',
             '2H-CAR-C6':'c12c3c(-c4cccc(-c5c6c(c7cc(CCCCCC)cc(-c8cccc(-c9c(c1cc(c9)CCCCCC)[nH]3)n8)c7[nH]6)cc(CCCCCC)c5)n4)cc(CCCCCC)c2',
             '2H-CAR-C8':'c12c3c(-c4nc(-c5c6c(c7c(c(-c8cccc(-c9c(c1cc(c9)CCCCCCCC)[nH]3)n8)cc(CCCCCCCC)c7)[nH]6)cc(CCCCCCCC)c5)ccc4)cc(c2)CCCCCCCC',
             '2H-CAR-Ph':"c12-c3c4c(c5cc(cc(-c6cccc(-c7c8c(c9cc(cc(-c(ccc1)n2)c9[nH]8)c1ccccc1)cc(c7)c1ccccc1)n6)c5[nH]4)c1ccccc1)cc(c3)c1ccccc1",
             '2H-CAR-Tol':'[C-]12=[CH+][CH-]=[CH+][C-](=[N-]2)[C+]2=[CH-][C+](=[CH-][C+]3=[C-]2N[C+]2=[C-]3[CH+]=[C-]([CH+][C-]2[C+]2=N[C+](=[CH-][CH+]=[CH-]2)[C-]2=[CH+][C-](=[CH+][C-]3=[C+]2N[C-]2=[C+]1[CH-]=[C+]([CH-]=[C+]32)[C-]1=[CH+][CH-]=[C+](C)[CH-]=[CH+]1)[C-]1=[CH+][CH-]=[C+](C)[CH-]=[CH+]1)[C-]1=[CH+][CH-]=[C+]([CH-]=[CH+]1)C)[C-]1=[CH+][CH-]=[C+](C)[CH-]=[CH+]1	',
             '2H-CAR-tBu':'c12c3c(-c4nc(-c5c6c(c7c(c(-c8cccc(-c9c(c1cc(c9)C(C)(C)C)[nH]3)n8)cc(C(C)(C)C)c7)[nH]6)cc(C(C)(C)C)c5)ccc4)cc(c2)C(C)(C)C',
             '2H-CAR-iPr':'c12c3c(-c4nc(-c5c6c(c7c(c(-c8cccc(-c9c(c1cc(c9)C(C)C)[nH]3)n8)cc(C(C)C)c7)[nH]6)cc(C(C)C)c5)ccc4)cc(c2)C(C)C	',
             }
print(f"{n_carp} assembler + {n_water} water in {box_size}nm box")

type='2H-CAR-Ph'
smiles=smiles_dict[type]

carp_mol=Molecule.from_pdb_and_smiles(f'{type}.pdb', smiles,allow_undefined_stereo=True)
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
integrator.setFriction(0.5/unit.picosecond)
platform=mm.Platform.getPlatformByName('CPU')
properties={'Threads': '14'} 

simulation=app.Simulation(matched_topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)
size_nm=box_size * unit.nanometer
simulation.context.setPeriodicBoxVectors([size_nm, 0, 0], [0, size_nm, 0], [0, 0, size_nm])

print("minimising energy")
simulation.minimizeEnergy(tolerance=1.0*unit.kilojoule_per_mole / unit.nanometer, maxIterations=2000)

with open(f'topology_{type}.pdb', 'w') as f:
    state=simulation.context.getState(getPositions=True)
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

runs=1000
run_ns=0.01
steps=int((run_ns * 1000000) / timestep)

direc=f"checkpoints_{type}"
if not os.path.exists(direc):
    os.makedirs(direc)

print("nvt")
nvt_steps = 100000  
simulation.step(nvt_steps)
nvt_path = os.path.join(direc, 'nvt.chk')
simulation.saveCheckpoint(nvt_path)

print("npt")
from openmm import MonteCarloBarostat
barostat = MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin, 25)
system.addForce(barostat)
simulation.context.reinitialize(preserveState=True)
npt_steps = 200000  
simulation.step(npt_steps)
npt_path = os.path.join(direc, 'npt.chk')
simulation.saveCheckpoint(npt_path)
print("equi over")

simulation.reporters.append(app.DCDReporter(f'assembly_{type}.dcd', 25000))
simulation.reporters.append(app.StateDataReporter(f'assembly_{type}.log', 25000, step=True, potentialEnergy=True, temperature=True, density=True))

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
print(f"time : {int(hours)} hrs, {mins:.1f} mins")

