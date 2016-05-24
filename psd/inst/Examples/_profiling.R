#
# Profile psd core functionality
#
library(psd)
do.prof <- function(PSD){
  nd <- 1e3 
  X.d <- arima.sim(list(order = c(1,1,0), ar = 0.9),n=nd)
  nt <- 8
  Rprof()
  for (i in 1:1000) mypsd <- PSD(X.d, ntaper=nt, plot=FALSE, refresh=TRUE)
  Rprof(NULL)
  prof <- summaryRprof()
  head(prof$by.self)
  head(prof$by.total)
  return(prof)
}
psdprof <- do.prof(PSD=psdcore)
print(str(psdprof))

# a better way:
library(profr)
p <- profr(pspectrum(rnorm(1e4), verbose=FALSE))
plot(p)
