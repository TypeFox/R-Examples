### This demo is the example from our arXiv preprint on kriging in the SN2011fe dataset
### To do the full size version of this example, you need to set doSmallExample to FALSE
### To do the full likelihood optimization, set doOptim to TRUE

if(require(fields)) {

doSmallExample = TRUE
doOptim = FALSE

if(doSmallExample){
  SN2011fe <- SN2011fe_subset
  SN2011fe_newdata <- SN2011fe_newdata_subset
  SN2011fe_mle <- SN2011fe_mle_subset
  nProc = 3
} else {
  nProc = 210  # set to an appropriate value for your system and for any queue submission request
}

n <- nrow(SN2011fe)
m <- nrow(SN2011fe_newdata)
r <- 100

nu <- 2

inputs <- c(as.list(SN2011fe), as.list(SN2011fe_newdata), nu = nu)

h_n = NULL 
h_m = NULL 

prob <- krigeProblem$new("prob", h_n = h_n, numProcesses = nProc, n = n, m = m, h_m = h_m, predMeanFunction = SN2011fe_predmeanfunc, crossCovFunction = SN2011fe_crosscovfunc,  predCovFunction = SN2011fe_predcovfunc, meanFunction = SN2011fe_meanfunc, covFunction = SN2011fe_covfunc,  inputs = inputs, params = SN2011fe_mle$par, data = SN2011fe$flux, packages = 'fields')

prob$calcLogDens()

if(doOptim) {
  prob$setParams(SN2011fe_initialParams)
  prob$optimizeLogDens(method = "L-BFGS-B", verbose = TRUE, lower = rep(.Machine$double.eps, length(SN2011fe_initialParams)), control = list(parscale = SN2011fe_initialParams))
}

pred <- prob$predict(ret = TRUE, se.fit = TRUE, verbose = TRUE)
realiz <- prob$simulateRealizations(r = r, post = TRUE, verbose = TRUE)

logwavelengths <- unique(SN2011fe_newdata$plogwavelength)
phases <- unique(SN2011fe_newdata$pphase)

par(mfrow = c(2, 2), mai = c(0.5, 0.5, 0.3, 0.1), mgp = c(1.8, 0.7, 0.0))

ord <- order(SN2011fe_newdata$pphase, SN2011fe_newdata$pwavelength)
image(logwavelengths, phases, matrix(pred$fit[ord], length(logwavelengths), length(phases)), xlab = "log wavelength", ylab = "phase", main = 'flux surface')

plot(SN2011fe_newdata$plogwavelength, pred$fit, xlab = "log wavelength", ylab = "flux", type = 'n', main = 'flux stratified by phase')
for(k in seq_along(phases)) {
  lines(SN2011fe_newdata$plogwavelength[SN2011fe_newdata$pphase == phases[k]], pred$fit[SN2011fe_newdata$pphase == phases[k]], col = k)
  text(8.7, pred$fit[SN2011fe_newdata$pphase == phases[k]][1], labels = as.character(phases[k]))
}
  
plot(logwavelengths, pred$fit[SN2011fe_newdata$pphase == 2], type = 'l', xlab = "log wavelength", ylab = "flux", main = 'flux for phase=2 with uncertainty')
lines(logwavelengths, pred$fit[SN2011fe_newdata$pphase == 2] - 2 * pred$se.fit[SN2011fe_newdata$pphase == 2], lty = 2)
lines(logwavelengths, pred$fit[SN2011fe_newdata$pphase == 2] + 2 * pred$se.fit[SN2011fe_newdata$pphase == 2], lty = 2)

plot(logwavelengths, pred$fit[SN2011fe_newdata$pphase == 2], xlab = "log wavelength", ylab = "flux", type = 'n', main = 'realizations of flux for phase=2')
for(rr in 1:r)
  lines(logwavelengths, realiz[SN2011fe_newdata$pphase == 2 , rr], col = 'gray')
lines(logwavelengths, pred$fit[SN2011fe_newdata$pphase == 2])

}
