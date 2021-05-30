import MDAnalysis as md
import numpy as np
import tensorflow as tf
#from nmrdata.nmrdata.parse import parse_universe
#from nmrdata.nmrdata import loading
from nmrdata import loading
import pickle
from simtk.openmm import app
import math
import os



'''
def pdb_to_records(path, embedding_dicts, max_atoms=1024):
	#load pdb
    pdb = app.PDBFile(path)
    topo = pdb.topology
    NN = 16
    atoms = np.zeros((max_atoms), dtype=np.int64) # atom type
	#mask = np.zeros((max_atoms), dtype=np.float)
    pos = pdb.getPositions(True) # positions
    print(len(pos))
    position = np.zeros((max_atoms-len(pos), 3), dtype=np.float)
	#peaks = np.zeros((max_atoms), dtype=np.float)
    names = np.zeros((max_atoms), dtype=np.int64)
    nlist = np.zeros((max_atoms, NN, 3), dtype=np.float)
    index = 0
    A = np.zeros((max_atoms, max_atoms), dtype=np.float)
    for residue in topo.residues():
        resname = residue.name
        for bond in residue.bonds():
            A[bond.atom1.index,bond.atom2.index] = 1
            A[bond.atom2.index, bond.atom1.index] = 1
    for atom in topo.atoms():
        atom_name = atom.residue.name + '-' + atom.name
        if atom_name not in embedding_dicts['name']:
            embedding_dicts['name'][atom_name] = len(embedding_dicts['name'])
        names[atom.index] = embedding_dicts['name'][atom_name]
		#mask[atom.index]=float(atom.element.symbol=='H')
        atoms[atom.index] = embedding_dicts['atom'][atom.element.symbol]
    positions =  np.concatenate((pos, position), axis=0)
    frame_nlist = nlist_tf_model(positions, NN)
    for i in range(frame_nlist.shape[0]):
        for j in range(NN):
            nlist[i,j,1] = frame_nlist[i,j,1]
            nlist[i,j,0] = frame_nlist[i,j,0]
            N_index = int(frame_nlist[i,j,1].numpy())
            if A[i,N_index] == 0:
                nlist[i,j,2] = embedding_dicts['nlist']['nonbonded']
            else:
                nlist[i,j,2] = embedding_dicts['nlist'][1]
    return nlist
    '''

    #n_list =n_list.numpy()
    #nlist = np.concatenate((n_list, nlist), axis=2)

    #return make_tfrecord(atoms, mask, nlist, peaks, embedding_dicts['class'][resname], names)

# MAX_ATOM_NUMBER=256
# NEIGHBOR_NUMBER=8
# result = []
# DATA_DIR = '../data'
# pep_dict = [
#   	'ala',
#         'arg',
#         'asn',
#         'asp',
#         'cys',
#         'gln',
#         'glu',
#         'gly',
#         'his',
#         'ile',
#         'leu',
#         'lys',
#         'met',
#         'phe',
#         'pro',
#         'ser',
#         'thr',
#         'trp',
#         'tyr',
#         'val']
#embedding_dicts = load_embeddings(os.path.join(DATA_DIR, 'embeddings.pb'))

# def parse_universe(u, neighbor_number, cutoff=None, pbc=False):
#     '''Converts universe into atoms, edges, nlist
#     '''
#     N = u.atoms.positions.shape[0]
#     dimensions = u.dimensions
#     if cutoff is None:        
#         cutoff = min(dimensions) / 2.01
#         if cutoff == 0:
#             # no box defined
#             bbox = u.atoms.bbox()
#             dimensions = bbox[1] - bbox[0]            
#             cutoff = min(dimensions) / 2.01
#             # make it into proper dimensions
#             dimensions = np.array(list(dimensions) + [90, 90, 90])
#             #u.atoms.wrap(box=dimensions)
#             warnings.warn('Guessing the system dimensions are' + str(dimensions))
#     gridsearch = md.lib.nsgrid.FastNS(cutoff, u.atoms.positions, dimensions, max_gridsize=N**2 // 2, pbc=pbc)
#     results = gridsearch.self_search()
#     ragged_nlist = results.get_indices()
#     ragged_edges = results.get_distances()
#     nlist = np.zeros((N, neighbor_number), dtype=np.int32)
#     edges = np.zeros((N, neighbor_number), dtype=np.float32)
#     #atoms = np.zeros(N, dtype=np.int32) 
#     # check for elements
#     #try:
#         #elements = u.atoms.elements
#     #except md.exceptions.NoDataError as e:
#         #warnings.warn('Trying to guess elements from names')
#         #elements = []
#         #for i in range(N):
#             # find first non-digit character
#             #elements.append([n for n in u.atoms[i].name if not n.isdigit()][0])
#     for i in range(N):
#         nl = ragged_nlist[i][:neighbor_number]        
#         nlist[i, :len(nl)] = nl
#         edges[i, :len(nl)] = ragged_edges[i][:neighbor_number]
#         #try:
#             #atoms[i] = embeddings['atom'][elements[i]]        
#         #except IndexError as e:
#             #warnings.warn('Unparameterized element' + u.atoms[i] + 'will replace with null')
#             #atoms[i] = 0
#     # note we convert atoms to be one hot
#     return edges, nlist

CGpdb = '../mtn/1A2P_CG.pdb'
ATOMpdb = '../PDB-without-ANISOU/1A2P_ATOM.pdb'
neighbor_num = 16
embedding_path = './embeddings_new.pb'
PDBFile = '../shiftx2-trainset-June2011/PDB-training/R001_1QRXA.pdb'
u = md.Universe(CGpdb)
A = loading.parse_universe(u, neighbor_num, embeddings=loading.load_embeddings(), pbc=True)
#A = parse_universe(u, neighbor_num)
#A = pdb_to_records(CGpdb, loading.load_embeddings(embedding_path))
#print(u.atoms)
print(A)