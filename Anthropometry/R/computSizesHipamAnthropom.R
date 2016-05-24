computSizesHipamAnthropom <- function(dataHip, bust, bustMeasur, nsizes, maxsplit, orness, type, ah, verbose = FALSE){
 res_hipam <- list() ; class(res_hipam) <- "hipamAnthropom"
 for(i in 1 : nsizes){
  data = dataHip[(bust >= bustMeasur[i]) &
                  (bust < bustMeasur[i + 1]), ]
  dataMat <- as.matrix(data)
  res_hipam[[i]] <- hipamAnthropom(dataMat, maxsplit = maxsplit, orness = orness, 
                                   type = type, ah = ah, verbose = FALSE)
  }
 
 return(res_hipam)
} 