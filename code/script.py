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




CGpdb = '../CGpdb/1A2P_CG.pdb'
neighbor_num = 16
u = md.Universe(CGpdb)
A = loading.parse_universe(u, neighbor_num, embeddings=loading.load_embeddings(), pbc=True)
#A = parse_universe(u, neighbor_num)
#A = pdb_to_records(CGpdb, loading.load_embeddings(embedding_path))
#print(u.atoms)
print(A)