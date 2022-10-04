# R package: gplSVCM

This package is developed for the implementation of the estimation and prediction for generalized partially linear models with spatially varying coefficients models(GPL-SVCM). The model is fitted through a parallel algorithm based on domain decomposition.

## Main functions and package structure

1. fit.gplsvcm.R: fit the GPL-SVCM by the BPST-DDC.

	1.1. fit.gplsvcm.worker.R: fit a local GPL-SVCM at each subregion.

		1.1.1. subdata.R: identify neighbor for each triangle and define a subregion as the triangle and its neighborhood, identify points within subregion and the corresponding triangulation.
		
			dat.dc.R: identify points falling within a triangulation
			ring.R: identify neighborhood triangle for each triangle

		1.1.2. mfit.gplsvcm.R: fit a local GPL-SVCM at each subregion.
		1.1.3. mfit.gsvcm.R: fit a local GSVCM at each subregion.

 
	1.2. ZZ.func.R: aggregate information from each subdomain and generate the covariance matrix of the estimated eta

2. predict.gplsvcm.R: predict the values of model terms for a given new dataset based on the estimated model from *fit.gplsvcm.ddc.R*.  

3. plot.gplsvcm.R: generate raster plots for the estimated coefficient functions. 

4. dat.generator.R: generate a dataset based on simulation studies in our manuscript.

5. beta.func.R: generate the values of true coefficient functions at given points. 