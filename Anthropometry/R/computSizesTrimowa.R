computSizesTrimowa <- function(dataTrim, bust, bustMeasur, nsizes, w, numClust, alpha, niter, algSteps, ah, verbose = FALSE){
 res_trimowa <- list() ; class(res_trimowa) <- "trimowa"
 for (i in 1 : nsizes){
  data = dataTrim[(bust >= bustMeasur[i]) &
                  (bust < bustMeasur[i + 1]), ]
  res_trimowa[[i]] <- trimowa(data, w, numClust, alpha, 
                              niter, algSteps, ah, verbose = FALSE)
 }

  return(res_trimowa)
} 
