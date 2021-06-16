# BEMA
Code for the paper "Estimation of the Number of Spiked Eigenvalues in a Covariance Matrix by Bulk Eigenvalue Matching Analysis".

BEMA0.R and BEMA.R contain the R code for the two algorithms. By default, the value of beta is 0.1 and M is 500. The value of alpha requires the user's input.  

When applying the BEMA algorithm to large datasets, one needs to adjust the data storage and optimizer used in BEMA according to data size and their computation budget. As an example, 1000Genome_code.ipynb contains the Python code for the BEMA algorithm customized for the 1000 Genome dataset.  
