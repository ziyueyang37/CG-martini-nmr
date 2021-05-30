import MDAnalysis as mda

u_itp = mda.Universe('../itpfiles/1A2P.itp', topology_format='ITP')
u_pdb = mda.Universe('../mtn/1A2P_CG.pdb')
print(u_itp.atoms.types)
print(len(u_itp.atoms.types))
print(u_pdb.atoms.types)
print(len(u_pdb.atoms.types))