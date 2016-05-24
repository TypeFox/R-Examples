### R code from vignette source 'bbefkr.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: bbefkr.Rnw:64-66
###################################################
library(bbefkr)
options(prompt = "R> ", bbefkr.messages = FALSE, digits = 3)


###################################################
### code chunk number 2: bbefkr.Rnw:198-236 (eval = FALSE)
###################################################
## ## install and load the R package
## install.packages("bbefkr")
## require(bbefkr)
## 
## ## set random seed
## set.seed(123456)
## 
## ## error density is approximated by a kernel density of residuals with a global 
## ## bandwidth, q=2 represents the semi-metric based on the second derivative, 
## ## since functions can be estimated by a basis representation, a B-spline 
## ## basis representation is used with 20 knots
## np_global <- bayMCMC_np_global(data_x=simcurve_smooth_normerr, 
##     data_y=simresp_np_normerr, data_xnew=simcurve_smooth_normerr, 
##     warm=1000, M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## estimated bandwidth parameters (h,b)
## np_global$xpfinalres
## 
## ## estimated values for the regression function
## np_global$mhat
## 
## ## simulation inefficiency factor to evaluate the convergence of MCMC
## ## it measures how many iterations are required to have the iid draws 
## ## from the posterior
## np_global$sif_value
## 
## ## log marginal likelihood computed by Chib (1995)'s method
## np_global$mlikeres
## 
## ## acceptance rates of the random-walk Metropolis sampling algorithm with target of
## ## 0.44 for the bandwidths in the regression function and kernal-form error density
## c(np_global$acceptnwMCMC, np_global$accepterroMCMC)
## 
## ## estimated error density (probability density function)
## np_global$fore.den.mkr
## 
## ## estimated error density (cumulative probability density function)
## np_global$fore.cdf.mkr


###################################################
### code chunk number 3: bbefkr.Rnw:241-267 (eval = FALSE)
###################################################
## ## error density is approximated by a kernel density of residuals with localised 
## ## bandwidths
## np_local <- bayMCMC_np_local(data_x=simcurve_smooth_normerr, 
##     data_y=simresp_np_normerr, data_xnew=simcurve_smooth_normerr, warm=1000, 
##     M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## estimated bandwidth parameters (h, b, badj)
## np_local$xpfinalres
## 
## ## estimated values for the regression function
## np_local$mhat
## 
## ## simulation inefficiency factor to evaluate the convergence of MCMC
## np_local$sif_value
## 
## ## log marginal likelihood computed by Chib (1995)'s method
## np_local$mlikeres
## 
## ## acceptance rates of the random-walk Metropolis algorithm for the bandwidths
## c(np_local$acceptnwMCMC, np_local$accepterroMCMC, np_local$acceptepsilonMCMC)
## 
## ## estimated error density (probability density function)
## np_local$fore.den.mkr
## 
## ## estimated error density (cumulative probability density function)
## np_local$fore.cdf.mkr


###################################################
### code chunk number 4: bbefkr.Rnw:279-303 (eval = FALSE)
###################################################
## ## error density is approximated by a kernel density of residuals with a global 
## ## bandwidth, Xvar is a n by 2 matrix, where each column variable is simulated 
## ## from U(0,1)
## semi_global <- bayMCMC_semi_global(data_x=simcurve_smooth_normerr,
##   data_y=simresp_semi_normerr, data_xnew=simcurve_smooth_normerr, 
##   Xvar=Xvar, warm=1000, M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## estimated regression coefficients
## semi_global$betahat
## 
## ## log marginal likelihood computed by Chib (1995)'s method
## semi_global$mlikeres
## 
## ## error density is approximated by a kernel density of residuals with localised 
## ## bandwidths
## semi_local <- bayMCMC_semi_local(data_x=simcurve_smooth_normerr,
##   data_y=simresp_semi_normerr, data_xnew=simcurve_smooth_normerr, 
##   Xvar=Xvar, warm=1000, M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## estimated regression coefficients
## semi_local$betahat
## 
## ## log marginal likelihood computed by Chib (1995)'s method
## semi_local$mlikeres


###################################################
### code chunk number 5: bbefkr.Rnw:317-331 (eval = FALSE)
###################################################
## ## error density is approximated by a kernel density of residuals with a global 
## ## bandwidth using the semi-metric based on the second derivative
## rough_np_global_deriv <- bayMCMC_np_global(data_x=simcurve_rough_normerr,
##   data_y=simresp_np_normerr, data_xnew=simcurve_rough_normerr, warm=1000, 
##   M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## a global bandwidth using the semi-metric based on three retained 
## ## principal components
## rough_np_global_pca <- bayMCMC_np_global(data_x=simcurve_rough_normerr,
## 	data_y=simresp_np_normerr, data_xnew=simcurve_rough_normerr, warm=1000, 
##   M=1000, semimetric="pca", q=3)
## 
## ## comparing two semi-metrics based on their log marginal likelihoods
## c(rough_np_global_deriv$mlikeres, rough_np_global_pca$mlikeres)


###################################################
### code chunk number 6: bbefkr.Rnw:336-350 (eval = FALSE)
###################################################
## ## error density is approximated by a kernel density of residuals with 
## ## localised bandwidths using the semi-metric based on the second derivative
## rough_np_local_deriv <- bayMCMC_np_local(data_x=simcurve_rough_normerr,
##   data_y=simresp_np_normerr, data_xnew=simcurve_rough_normerr, warm=1000, M=1000,
## 	range.grid=c(0,pi), q=2, nknot=20)
## 
## ## localised bandwidths using the semi-metric based on three retained 
## ## principal components
## rough_np_local_pca <- bayMCMC_np_local(data_x=simcurve_rough_normerr,
## 	data_y=simresp_np_normerr, data_xnew=simcurve_rough_normerr, warm=1000, 
##   M=1000, semimetric="pca", q=3)
## 
## ## comparing two semi-metrics based on their log marginal likelihoods
## c(rough_np_local_deriv$mlikeres, rough_np_local_pca$mlikeres)


###################################################
### code chunk number 7: bbefkr.Rnw:355-383 (eval = FALSE)
###################################################
## ## error density is approximated by a kernel density of residuals with a global 
## ## bandwidth using the semi-metric based on the second derivative
## smooth_np_global_deriv <- bayMCMC_np_global(data_x=simcurve_smooth_normerr,
##   data_y=simresp_np_normerr, data_xnew=simcurve_smooth_normerr, warm=1000, 
##   M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## a global bandwidth using the semi-metric based on functional 
## ## principal components
## smooth_np_global_pca <- bayMCMC_np_global(data_x=simcurve_smooth_normerr,
## 	data_y=simresp_np_normerr, data_xnew=simcurve_smooth_normerr, warm=1000, 
##   M=1000, semimetric="pca", q=3)
## 
## ## comparing two semi-metrics based on their log marginal likelihoods
## c(smooth_np_global_deriv$mlikeres, smooth_np_global_pca$mlikeres)
## 
## ## localised bandwidths using the semi-metric based on the second derivative
## smooth_np_local_deriv <- bayMCMC_np_local(data_x=simcurve_smooth_normerr,
## 	data_y=simresp_np_normerr, data_xnew=simcurve_smooth_normerr, warm=1000, 
##   M=1000, range.grid=c(0,pi), q=2, nknot=20)
## 
## ## localised bandwidths using the semi-metric based on functional 
## ## principal components
## smooth_np_local_pca <- bayMCMC_np_local(data_x=simcurve_smooth_normerr,
## 	data_y=simresp_np_normerr, data_xnew=simcurve_smooth_normerr, warm=1000, 
##   M=1000, semimetric="pca", q=3)
## 
## ## comparing two semi-metrics based on their log marginal likelihoods
## c(smooth_np_local_deriv$mlikeres, smooth_np_local_pca$mlikeres)


###################################################
### code chunk number 8: bbefkr.Rnw:399-444 (eval = FALSE)
###################################################
## ## We use the first 44 pairs of data to estimate the relationship. Based on a 
## ## new curve, we can then predict its response using a functional 
## ## nonparametric regression
## fat_np_global <- bayMCMC_np_global(data_x = specurves[1:44,], data_y = fat[1:44], 
##   		data_xnew = specurves[45,], range.grid=c(0,pi), q=2, nknot=20)
## 
## ## Point forecast
## fat_np_global$pointforecast
## 
## ## 95% prediction interval
## fat_np_global$PI
## 
## ## Using a functional nonparametric regression with localised bandwidths
## fat_np_local <- bayMCMC_np_local(data_x = specurves[1:44,], data_y = fat[1:44], 
##   data_xnew = specurves[45,], range.grid=c(0,pi), q=2, nknot=20)
## 
## ## Point forecast
## fat_np_local$pointforecast
## 
## ## 95% prediction interval
## fat_np_local$PI
## 
## ## Using a semi-functional partial linear regression with a global bandwidth
## fat_semi_global <- bayMCMC_semi_global(data_x = specurves[1:44,], data_y = fat[1:44], 
## 	data_xnew = specurves[45,], Xvar = cbind(protein[1:44], moisture[1:44]), 
## 	Xvarpred = matrix(c(protein[45], moisture[45]), nrow=1), 
## 	range.grid=c(0,pi), q=2, nknot=20)
## 
## ## Point forecast
## fat_semi_global$pointforecast
## 
## ## 95% prediction interval
## fat_semi_global$PI
## 
## ## Using a semi-functional partial linear regression with localised bandwidths
## fat_semi_local <- bayMCMC_semi_local(data_x = specurves[1:44,], data_y = fat[1:44], 
## 	data_xnew = specurves[45,], Xvar = cbind(protein[1:44], moisture[1:44]),
## 	Xvarpred = matrix(c(protein[45], moisture[45]), nrow=1), 
## 	range.grid=c(0,pi), q=2, nknot=20)
## 
## ## Point forecast
## fat_semi_local$pointforecast
## 
## ## 95% prediction interval
## fat_semi_local$PI


###################################################
### code chunk number 9: bbefkr.Rnw:449-459 (eval = FALSE)
###################################################
## fat_np_global_d1 <- bayMCMC_np_global(data_x = specurves[1:44,], data_y = fat[1:44], 
##   data_xnew = specurves[45,], range.grid=c(0,pi), q=1, nknot=20)
##                                    
## fat_np_global_d2 <- bayMCMC_np_global(data_x = specurves[1:44,], data_y = fat[1:44], 
## 	data_xnew = specurves[45,], range.grid=c(0,pi), q=2, nknot=20)                                   
##                                    
## fat_np_global_pca <- bayMCMC_np_global(data_x = specurves[1:44,], data_y = fat[1:44], 
## 	data_xnew = specurves[45,], semimetric="pca", q=3)
## 
## c(fat_np_global_d1$mlikeres, fat_np_global_d2$mlikeres, fat_np_global_pca$mlikeres)


