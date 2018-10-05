# 1D-enthalpy-method

Enthalpy Method Code  as described in 
A geomorphic enthalpy method: Description and application to the evolution 
of fluvial-deltas under sea-level cycles
William Anderson, Jorge Lorenzo-Trueba, Vaughan Voller

1. Intro
The included code provides three scripts. The first is ‘Enthalpy_Code_Trajectories’ which plots 
trajectories of the shoreline (SH) and alluvial-bedrock transition (ABT) for sea-level change of 
the form Z=A*t^B where A and B are constants. Additionally, the function ‘NonlineSys’ 
provides an analytical solution under sea-level change proportional to the square root of time and 
the script ‘Enthalpy_Code_Movie’ creates a movie of the deltaic system as it grows.

2. Code Description
The script ‘Enthalpy_Code_Trajectories’ takes 3 paramater inputs: Rab, A, and B. Rab is the 
dimensionless parameter qin/(beta*nu), where qin=incoming sediment flux, beta = basement 
slope, and nu = diffusivity. Rab is physically bounded by 0<Rab<1. A and B determine the 
magnitude and rate of sea-level change.
 
If B=0.5, comparison of the numerical prediction of SH/ABT trajectories will be plotted with the 
analytical trajectories. under sea-level change proportional to the square root of time the 
analytical trajectories are determined by the ‘NonlineSys’ function which solves two nonlinear
equations described in Lorenzo-Trueba and Voller 2013.
 
If B=1 and a>0, the numerical trajectories will be plotted with the analytical position of the 
boundaries once the fluvial surface maintains a fixed a geometry. In this scenario the numerical 
trajectories should eventually match the analytical trajectories.
 
Choosing any other values of B will plot only numerical trajectories of the SH and ABT since no 
analytical solution exists.
 
Additionally, setting the parameter m=1 stores several fluvial profiles throughout a run so that a 
movie can be created by the script “Enthalpy_Code_Movie”.



