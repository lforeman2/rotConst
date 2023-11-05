# rotConst
Calculate rotational constant from .xyz coordinates using python

Input molecular geometry input coordinates ('geom.xyz') should be in units of Angstroms and follow the following format where k is the number of atoms:

  $$ 
  coords=
    \begin{bmatrix}
    x_0 & y_0 & z_0\\
    x_1 & y_1 & z_1\\
    x_2 & y_2 & z_2\\
    ... & ... & ...\\
    x_k & y_k & z_k\\
    \end{bmatrix}
  $$

The rotational constants are given by:  
$$B=\frac{h}{8\pi^2cI}$$  

For a polyatomic molecule, I is the moment of inertia tensor:

  $$ 
  I=
    \begin{bmatrix}
    I_{xx} & I_{xy} & I_{xz}\\
    I_{yx} & I_{yy} & I_{yz}\\
    I_{zx} & I_{zy} & I_{zz}\\
    \end{bmatrix}
  $$

Each element of the moment of inertia tensor can be calculated by:  

$$
    masses=
    \begin{bmatrix}
    m_{0} & m_{1} & m_{2} & ... & m_{k}
    \end{bmatrix}
$$  

$$I_{xx} = \sum_{k=1}^n m_{k} ( y_{k}^2 + z_{k}^2 )$$   

$$I_{yy} = \sum_{k=1}^n m_{k} ( x_{k}^2 + z_{k}^2 )$$  

$$I_{zz} =  \sum_{k=1}^n m_{k} ( x_{k}^2 + y_{k}^2 )$$  

$$I_{ij} = -\sum_{k=1}^n m_{k} ( i_{k}^2 + j_{k}^2 )$$  

## To run the program:
```bCalc(input_file='geom.xyz', units='wavenumber')```  
The units options are 'wavenumber' or 'gigahertz'.  
The input file should have the following format (using water as an example):  
```
  3
Angstrom
  O  0.00000   0.00000   0.11730
  H  0.00000   0.75720  -0.4692
  H  0.00000  -0.75720  -0.4692
```
