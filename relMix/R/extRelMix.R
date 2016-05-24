extRelMix <-
function(E, AF, marker, R, theta, ped, silent = 0){
  datamatrix <- generate( E, K = NULL, 2)
  AF <- rep( AF, dim(datamatrix)[2]/2)
  datamatrix <- rbind(datamatrix,AF)
  datamatrix <- as.data.frame(datamatrix)
  rownames(datamatrix)[c(1,2)] = c( "CH", "MO")
  alleles <- sort(unique(c(E,AF)))
  dd <- db[db$Marker==marker,]
  freqs <- dd[,3][dd$Allel%in%alleles]
  if(silent == 0) {
     ff <- c(freqs,1-sum(freqs))
     al <- c(alleles,99) #rest=99
     }
  if(silent > 0){
     rest <- 1-sum(freqs)-silent
     ff <- c(freqs,rest,silent)
     al <- c(alleles,99,999) #rest = 99,silent = 999
     }
  locus <- FamiliasLocus( frequencies = ff, name = marker,  
                          allelenames = al, femaleMutationRate = R, maleMutationRate = R)
  relMix( ped, locus, datamatrix,kinship = theta)
}
