import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import Bio.PDB
import math
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Vector import calc_dihedral
from Bio.PDB.Vector import Vector
import Bio.PDB as PDB


# def get_phi_psi(chain,AA):
#     """
#     return all the phi angle and psi angle for a chain based for given alphabet. The function does not work as desired.
#     """
#     amino_acid = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
#                "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
#
#     # store all the residue for a chain into a list
#     residues = []
#     for res in chain:
#         if res.get_resname() in amino_acid:
#             residues.append(res)
#
#     # get N, CA, C atom's coordinate and store it into correspond list
#     phi_psi_list = []
#     length = len(residues)
#     for i in range(length):
#         res = residues[i]
#         if res.get_resname() in AA:
#             n = res['N'].get_vector()
#             ca = res['CA'].get_vector()
#             c = res['C'].get_vector()
#             # get phi angle and psi angle
#             if i > 0:
#                 try:
#                     prev_res = residues[i - 1]
#                     prev_c = prev_res["C"].get_vector()
#                     phi = PDB.calc_dihedral(prev_c, n, ca, c)
#                 except Exception:
#                     phi = None
#                     continue
#             # the first residue does not have phi angle
#             else:
#                 phi = None
#             if i < (length - 1):
#                 try:
#                     next_res = residues[i + 1]
#                     next_n = next_res["N"].get_vector()
#                     psi = PDB.calc_dihedral(n, ca, c, next_n)
#                 except Exception:
#                     psi = None
#             # the last residue does not have psi angle
#             else:
#                 psi = None
#             phi_psi_list.append((phi, psi))
#
#     new_phi_psi_angles = phi_psi_list[1:-1]
#     return new_phi_psi_angles



def get_phi_psi(structure):
    """
    Calculate phi,psi dihedral angles and return lists.
    Uses the polypeptide class."""

    # Create a list of  polypeptide objects
    ppb = PDB.PPBuilder()
    pp_list = ppb.build_peptides(structure)

    # Get phi and psi angles
    phi_angles_list = []
    psi_angles_list = []

    # Iterate over polypeptide molecules
    for pp in pp_list:

        # Calculate phi and psi angles and unpack list and tuple
        Agg_phi = []
        Agg_psi = []

        for phi,psi in pp.get_phi_psi_list():

            # put them in the lists
            Agg_phi.append(phi)
            Agg_psi.append(psi)

        phi_angles_list.append(Agg_phi)
        psi_angles_list.append(Agg_psi)

    return phi_angles_list, psi_angles_list


# phi_angles_list, psi_angles_list = get_phi_psi()
# print "Phi :"
# print phi_angles_list
# print "Psi :"
# print psi_angles-list

# These are different list that store the values or angles for plotting of the graph.

phi=[]
psi=[]

phi_i=[]
phi_i_plus_one=[]

phi_i_1=[]
psi_i_plus_one=[]

psi_i=[]
phi1_i_plus_one=[]

psi1_i=[]
psi1_i_plus_one=[]


#Helps tp parse the files
parser = Bio.PDB.PDBParser()

#Path of the folder
database = '/home/jammy/PycharmProjects/MyProjects/top100H' #("Please change the path to run the code")

#Helps to open the folder to parse the files.
files = os.listdir(database)

#Implementing the parser

for f in files:
    try:
        s = parser.get_structure(f, database + '/' + f)   #goes through every file in the folder
        phi_angles_list, psi_angles_list = get_phi_psi(s) #calls the function to calculate the angles.

        for i in range(len(phi_angles_list)):
            #Standard Ramachandran data
            if (None not in phi_angles_list[i][1:-1] and None not in psi_angles_list[i][1:-1]): #Indexing for the None values which occur in
                phi.extend(phi_angles_list[i][1:-1])  #Takes the first value                    #specific points in Phi and Psi angles.
                psi.extend(psi_angles_list[i][1:-1])  #Takes the last value

            #Phi "i" & Phi "i+1"
            if (None not in phi_angles_list[i][1:-1] and None not in phi_angles_list[i][2:]):
                phi_i.extend(phi_angles_list[i][1:-1])          #Takes the second values from the phi1 list leaving the None. and then the one before last.
                phi_i_plus_one.extend(phi_angles_list[i][2:])   #Takes the third value till the last ...

            #Phi "i" & Psi "i+1"
            if (None not in phi_angles_list[i][1:-2] and None not in psi_angles_list[i][2:-1]):
                phi_i_1.extend(phi_angles_list[i][1:-2])        #Take the first value  and second from last
                psi_i_plus_one.extend(psi_angles_list[i][2:-1]) #Takes the second value and last value

            #Psi "i" & Phi "i+1"
            if (None not in psi_angles_list[i][:-1] and None not in phi_angles_list[i][1:]):
                psi_i.extend(psi_angles_list[i][:-1])           #Takes from zero and only the last value
                phi1_i_plus_one.extend(phi_angles_list[i][1:])  #Takes the first value and from the last

            #Psi "i" & Psi "i+1"
            if (None not in psi_angles_list[i][:-2] and None not in psi_angles_list[i][1:-1]):
                psi1_i.extend(psi_angles_list[i][:-2])           #Takes from the first value and second from last
                psi1_i_plus_one.extend(psi_angles_list[i][1:-1]) #Takes the first value and one from the last

     #Try to catch the faulty file and puts in the seprate folder
    except Exception as inst:
        with open('/home/jammy/PycharmProjects/MyProjects/error/error' + f + '.txt', 'w') as outfile:
            outfile.write(f)

# Following lines are used for plotting the data.

plt.figure(1)
plt.hist2d(phi, psi, bins=100, norm=LogNorm())
plt.xlabel('phi')
plt.ylabel('psi')
plt.colorbar()
plt.show()


plt.figure(2)
plt.hist2d(phi_i, phi_i_plus_one, bins=100, norm=LogNorm())
plt.xlabel('phi')
plt.ylabel('phi+1')
plt.colorbar()
plt.show()


plt.figure(3)
plt.hist2d(phi_i_1, psi_i_plus_one, bins=100, norm=LogNorm())
plt.xlabel('phi')
plt.ylabel('psi+1')
plt.colorbar()
plt.show()


plt.figure(4)
plt.hist2d(psi_i, phi1_i_plus_one, bins=100, norm=LogNorm())
plt.xlabel('psi')
plt.ylabel('phi+1')
plt.colorbar()
plt.show()


plt.figure(5)
plt.hist2d(psi1_i, psi1_i_plus_one, bins=100, norm=LogNorm())
plt.xlabel('psi')
plt.ylabel('psi+1')
plt.colorbar()
plt.show()



#-----------------------------------------------------------------------------------------------------------------------
#                                                      "RNA CODE"
#=======================================================================================================================


#!/usr/bin/python
# Copyright 2012 Stefan E Seemann <seemann@rth.dk>
import math as m
import sys
import numpy as np
import scipy.stats as st



def pi():
    '''This function takes raw input a sequence. It calculates the nussinov calculations for the M matrix.'''

    # Takes input from the user.
    seq = raw_input("Enter sequence -> ")
    seql=[i for i in seq.upper()]
    l=len(seq)

    #loop size
    h=3

    #basepair scores
    scores={'AU':2, 'UA':2, 'GU':1, 'UG':1, 'GC':3, 'CG':3}

    #initialize array
    m = [[0 for i in range(l)]
                    for j in range(l)]

    #fill scoring matrix
    for j0 in range(h+1,l):
        for i in range(0,l-j0):
            j=i+j0
            #rule 1) i,j paired
            if seql[i]+seql[j] in scores:
                m[i][j]=m[i+1][j-1]+scores[seql[i]+seql[j]]


            #rule 2) i unpaired
            if m[i+1][j] > m[i][j]:
                m[i][j]=m[i+1][j]

            #rule 3) j unpaired
            if m[i][j-1] > m[i][j]:
                m[i][j]=m[i][j-1]

            #rule 4) bifurcation k
            for k in range(i+1+h,j-1-h):
                if m[i][k]+m[k+1][j] > m[i][j]:
                    m[i][j]=m[i][k]+m[k+1][j]

    c = [[0 for i in range(l)]for j in range(l)] # Initiating the c-matrix array.

    #backtracking
    str=['.' for i in range(l)]
    stack=[]
    stack.append([0,l-1])
    while len(stack)>0:
        top=stack.pop(),
        i=top[0][0]
        j=top[0][1]
        if i>=j:
            continue
        elif m[i+1][j] == m[i][j]:
            stack.append([i+1,j])
        elif m[i][j-1] == m[i][j]:
            stack.append([i,j-1])
        elif seql[i]+seql[j] in scores and m[i+1][j-1]+scores[seql[i]+seql[j]] == m[i][j]:
            if scores[seql[i]+seql[j]] == 3:
                c[i][j] = 3
                c[j][i] = 3
            elif scores[seql[i]+seql[j]] == 2:
                c[i][j] = 2
                c[j][i] = 2
            elif scores[seql[i]+seql[j]] == 1:
                c[i][j] = 1
                c[j][i] = 1
            str[i]="("
            str[j]=")"
            stack.append([i+1,j-1])
        else:
            for k in range(i+1+h,j-1-h):
                if m[i][k]+m[k+1][j] == m[i][j]:
                    stack.append([k+1,j])
                    stack.append([i,k])
                    break

    #output
    print seq,l,"\n",''.join(str),"\n","Score:",m[0][l-1], "\n"
    # print "Matrix M \n"
    # for k in range(0,l):  #For printing the matrix form of c matrix
    #     print c[k]


    # Gathering the Pi values for pearson coeffcient.
    pi_list = []
    for i in range(l):
        pi = 0 # Initiates the value of Pi from zero.
        for j in range(l):
            pi += c[i][j] # This line feeds the value from c matrix to pi lists and adds the values.
        pi_list.append(pi)

    # print pi_list
    return pi_list

# Calling the function Pi() for wild and all mutant types.
wild = pi()
mut1 = pi()
mut2 = pi()
mut3 = pi()
mut4 = pi()

def pearson_correlation(wild, mut):
    '''This function calculates Pearson Correlation Coefficient'''

    pearson_correlation = st.pearsonr(wild,mut)
    print pearson_correlation

print pearson_correlation(wild,mut1)
print pearson_correlation(wild,mut2)
print pearson_correlation(wild,mut3)
print pearson_correlation(wild,mut4)







