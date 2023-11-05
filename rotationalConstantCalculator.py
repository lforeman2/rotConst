#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import chemlib
from scipy.constants import h, c, pi, nano, physical_constants
from chemlib import Element


# Get input geometry
def get_input(input_file='geom.xyz'):
    df = pd.read_table(input_file, skiprows=2, delim_whitespace=True, names=['atom', 'x', 'y', 'z'])
    atoms=np.array(df.atom)
    coords=np.array(df[['x', 'y', 'z']])
    return atoms, coords


# Get masses of each atom
def get_masses(atoms):
    masses=[]
    for i in range(0,len(atoms)):
        elementSymbol=atoms[i]
        element = Element(elementSymbol)
        masses.append(element.AtomicMass)
    return masses

# Calculate center of mass
def center_of_mass(atoms, coords, masses):
    RxCM=np.sum(masses*coords[:,0])/sum(masses)
    RyCM=np.sum(masses*coords[:,1])/sum(masses)
    RzCM=np.sum(masses*coords[:,2])/sum(masses)
    xyzCOM=coords-np.array([RxCM*np.ones(len(atoms)), RyCM*np.ones(len(atoms)), RzCM*np.ones(len(atoms))]).T
    return xyzCOM
    
# Calculate moment of inertia tensor
def moment_of_inertia(atoms, xyzCOM):
    I=np.zeros((3,3))
    
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,len(atoms)): # k = 0:O, 1:H, 2:H
                if i==j:
                    I[i,j]=I[i,j]+masses[k]*(sum(xyzCOM[k,:]**2)-xyzCOM[k,i]**2)
                else:
                    I[i,j]=I[i,j]-masses[k]*(xyzCOM[k,i]*xyzCOM[k,j])

    Idiag, evs =np.linalg.eig(I) #amu*A^2

    return Idiag


# Compute rotational constant
def rot_const(units='wavenumber'):
    # Convert moments of inertia units
    m=physical_constants['atomic mass constant'][0]
    Ikgm=Idiag*m*(1E-10)**2 # amu*A^2 --> kg*m^2
    if units=='wavenumber':
        f=0.01
    if units=='gigahertz':
        f=c*nano
    B=f*(h/(8*pi**2*c))*(1/Ikgm)
    return B


def calculator(input_file='geom.xyz', units='wavenumber'):
    atoms, coords=get_input(input_file)
    masses=get_masses(atoms)
    xyzCOM=center_of_mass(atoms, coords, masses)
    Idiag=moment_of_inertia(atoms, xyzCOM)
    B=rot_const(units)
    return B
