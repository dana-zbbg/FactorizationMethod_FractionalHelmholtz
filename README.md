# FactorizationMethod_FractionalHelmholtz
This selfcontained code enables to generate far field data for the fractional inhomogeneous Helmholtz equation and contains an implementation of the factorization method to recover the support of the inhomogeneity from far field data. 


The file far_field_main.m computes the far field operator 
with the indicated settings (wave number k, order s, Ninc number of incident fields...).

The far field operator matrix is saved in FF_data.mat file, also containing the inhomogeneity. 

The file Lippmann_SchwingerMatrix_assemble.m is called in far_field_main.m to compute the 
fundamental solution and to assemble the matrix for the discretized Lippmann Schwinger equation.

The files cheval.m, StruveH0.m, StruveH0Y0.m are used to compute the 
fundamental solution in the file LippmannSchwingerMatrix_assemble.m. They can be found at 
https://www.mathworks.com/matlabcentral/fileexchange/37302-struve-functions. 
The license is included the third_party_license.txt. The file hankel_transform.m is also used to
computed the fundamental solution. 

Once the data for the far field operator is obtained, the file FactorizationMethod.m computes the
indicator function W from the far field data, and plots the reconstruction results. The indicator
function W is saved in a 'FM.mat' file.

Under the far_field_data folder, one can find pre-computed far field data for the 'boomerang' 
shape inhomogeneity and the 'cercle+rectangle' inhomogeneity. 

Under the figure folder, one can find pre-computed W indicator functions for the cercle+rectangle
and boomerang inhomogeneities, and files that enables to display figures of the reconstructions. 
