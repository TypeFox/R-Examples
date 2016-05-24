#1 - gen.findMRCA
#2 - gen.climbPAR
#3 - gen.getAncestorsPAR
#4 - gen.findFounders
#5 - gen.getFoundersPAR
#6 - gen.findDistance
#7 - gen.find.Min.Distance.MRCA


#####
# gen.findMRCA finds the Most Recent Common Ancestors (MRCA) of the specified individuals 
#  using the given genealogy. It then calculates the distance (number of meiosis)
#  between the individuals, passing by each of the MRCAs previously found.
# This function spawns slave processes that run gen.climbPAR and gen.getAncestorsPAR.
#  
#  - gen -> the genealogy
#  - individuals -> the individuals to consider
#  - NbProcess -> Number of process to use when running this function. default=parallel::detectCores()-1
gen.findMRCA <- function(gen, individuals, NbProcess=parallel::detectCores()-1)
{
 retour = tryCatch({
 # create the processes
 NbProcess = min(length(individuals),NbProcess)
#  cluster = makeCluster(NbProcess, type="SOCK")
  cluster = makePSOCKcluster(NbProcess)
#  registerDoSNOW(cluster)
#  setDefaultCluster(cluster)
  registerDoParallel(cluster)
 
 clusterExport(cluster, c("gen.getAncestorsPAR"))
 
 n = length(individuals)
 x = foreach(i = 1:n) %dopar% { gen.getAncestorsPAR(gen,individuals[i]) }
 
 # intersection of the common ancestors
 indInter <- c()
 if(length(x) > 1) {
  indInter <- intersect(x[[1]], x[[2]])
  if(length(x)>2) z <- lapply(c(3:length(x)), function(i){ indInter <<- intersect(indInter, x[[i]]) })
 }
 else
  indInter <- x[[1]]
 

 # drop the ancestors having children in the intersection
 b_1 <- gen.genout(gen.branching(gen, individuals[1] ))
 
 inter <- subset(b_1, b_1$ind %in% indInter, T)
 inter <- subset(inter, !(inter$ind %in% inter$father | inter$ind %in% inter$mother), T)
 print(paste(dim(inter)[1],"MRCA"))
 
 ### SECOND: Calculate the distance matrix
 clusterExport(cluster, c("gen.climbPAR"))
 y <- foreach(i = 1:length(inter$ind)) %dopar% gen.climbPAR(gen=gen,individuals=individuals,founder=inter$ind[i])

 x <- matrix(nrow=length(individuals), ncol=length(inter$ind))
 for(i in 1:length(inter$ind)){ 
   message <- y[[i]]
   f <- message$founder
   d <- message$distance
   x[,i] <- d
 }
 colnames(x) <- inter$ind
 rownames(x) <- individuals
 x
 }, warning = function(w){ print("warning"); message(w); return(NULL)
 }, error   = function(e){ print("error");   message(e); return(NULL)
 }, finally = { stopCluster(cluster); })
 return(retour)
}

#####
gen.climbPAR <- function(gen, individuals, founder) {
# require(GENLIB)
 dist <- as.numeric(unlist(lapply( individuals, function(deb){ gen.min(gen.branching(gen, deb), founder) })))
 list( founder=founder, distance=dist )
}

#####
gen.getAncestorsPAR <- function(gen, pro) {
# require(GENLIB)
 genTmp <- gen.branching(gen, pro)
 genOut <- gen.genout(genTmp)
 genOut$ind
}


#####
gen.findFounders <- function(gen, individuals, NbProcess=parallel::detectCores()-1) {
 retour = tryCatch({
 
 NbProcess = min(length(individuals),NbProcess)
#  cluster = makeCluster(NbProcess, type="SOCK")
  cluster = makePSOCKcluster(NbProcess)
#  registerDoSNOW(cluster)
#  setDefaultCluster(cluster)
  registerDoParallel(cluster)
 
 clusterExport(cluster, c("gen.getFoundersPAR"))
 
 n = length(individuals)
 x = foreach(i = 1:n) %dopar% { gen.getFoundersPAR(gen,individuals[i]) }
 i=0
 
 res <- c()
 if(length(x) > 1) {
  res <- intersect(x[[1]], x[[2]])
  if(length(x) > 2) z <- lapply(c(3:length(x)), function(i){ res <<- intersect(res, x[[i]]) })
 }
 else if(length(x) == 0) res <- c()
 else				res <- x[[1]]
 
 res
 }, warning = function(w){ print("warning"); message(w); return(NULL)
 }, error   = function(e){ print("error");   message(e); return(NULL)
 }, finally = { stopCluster(cluster); })
 return(retour)
}

#####
# gen.getFoundersPAR is not to be run by itself. 
# It is to be launched through gen.findFounders.
gen.getFoundersPAR <- function(gen, pro) {
# require(GENLIB)
 genTmp <- gen.branching(gen, pro)
 gen.founder(genTmp)
 #found <- gen.founder(genTmp)
 #list(founders=found)
}
 
#####
# gen.findDistance finds the minimum distance (number of meiosis) between the specified
#  individuals and the given founder, using the given genealogy. 
#  - gen = a given genealogy
#  - individuals = individuals in the genealogy, descendants of the given founder.
#  - ancestor = a common ancestor to individuals
gen.findDistance <- function(gen, individuals, ancestor) {
 if(!is.numeric(individuals)) stop("individuals must be numeric");
 if(!is.numeric(ancestor))	stop("ancestor must be numeric");
 if(length(individuals)!=2)	stop("There must be two and only two individuals.")
 if(length(ancestor)>1) 	warning(paste("Only the first element of 'ancestor' will be considered:",ancestor[1]))
 sum(as.numeric(unlist(lapply(individuals, function(deb){ gen.min(gen.branching(gen, deb), ancestor[1]) }))))
}

#####
# gen.find.Min.Distance.MRCA finds the smallest distances (number of meiosis) between
#  pairs of probands given the matrix of MRCAs output by the gen.findMRCA function.
#  - genMatrix = output of gen.findMRCA
#  - individuals = probands to consider in the matrix.
#  - ancestors = MRCAs to consider in the matrix.
gen.find.Min.Distance.MRCA <- function(genMatrix, individuals="ALL", ancestors="ALL") {
 if(ancestors[1] == "ALL") ancestors <- as.numeric(colnames(genMatrix))
 if(     individuals[1] == "ALL") individuals      <- as.numeric(rownames(genMatrix))
 
 if(!is.matrix(genMatrix)){ stop("genMatrix must be a matrix output from gen.findMRCA"); }
 if(sum(unlist(lapply(individuals, function(p){p%in%rownames(genMatrix)}))) != length(individuals)){
  stop("Some individuals are not in the given matrix\n");
 }
 if(sum(unlist(lapply(ancestors, function(f){f%in%colnames(genMatrix)}))) != length(ancestors)){
  stop("Some ancestors are not in the given matrix\n");
 }
 
 res <- matrix(unlist(lapply( ancestors               , function(z) {
			 tmp <- matrix(unlist(lapply( c(1:(length(individuals)-1))  , function(i) {
			  			 lapply(c((i+1):(length(individuals))), function(j) {
						   c( individuals[i], individuals[j], (genMatrix[i,as.character(z)] + genMatrix[j,as.character(z)]) )
						 })
						})), ncol=3, byrow=T)
			 posMin <- grep(min(as.numeric(tmp[,3])),as.numeric(tmp[,3]))
			 if(length(posMin) > 1) t(cbind( rep(z, length(posMin)), tmp[posMin,] ))
			 else 		 c( z, tmp[posMin,] )
			})), ncol=4, byrow=T)
 colnames(res) <- c("founder","proband1","proband2","distance")
 res
}

