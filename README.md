stochastic-or-model
===================

So far this repo contains two 
Matlab mex files for stochastic simulation based on the 
Gillespie algorhithm

gillespie_ssim_sparse.cpp is a generic Gillespie solver 
which is optimized for sparse (or small) reaction matricies.

gillespie_ssim_ORexpMod.cpp is a solver for a specific olfactory receptor 
expression model found in arxiv.org/abs/1201.2933

both compile in Matlab with: mex <filename>
