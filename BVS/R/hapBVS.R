hapBVS = function(G,min.Hap.freq=.02){

  library(haplo.stats)

  ##Convert G to allele matrix for haplo.stats
  codeAlleles <- function(g) {
	a1 <- ifelse(is.na(g), NA, ifelse(g<=1, 0, 1))
	a2 <- ifelse(is.na(g), NA, ifelse(g>=1, 1, 0))
	alleles <- cbind(a1,a2)
	alleles}

  codeG.AsAlleles <- function(G) {
	v <- {}
	for(m in 1:ncol(G)) { v <- cbind(v, codeAlleles(G[,m])) }
	v}
	
  G.alleles <- as.data.frame(codeG.AsAlleles(G))  
  names(G.alleles) <- rep(names(G), 1, each=2)

  ##Estimate haplotypes
  hap.obj <- haplo.em(G.alleles, locus.label=names(G), miss.val=NA)
  
  ##Create Posterior Haplotype Matrix and Return Estimated Haplotypes  
  createPostHapMatrix <- function(hap.obj, min.Hap.freq=.025) {
	hap.prob <- hap.obj$hap.prob
	haplotypes <- hap.obj$haplotype
 
	hap1 <- hap.obj$hap1code
	hap2 <- hap.obj$hap2code
	indx <- hap.obj$indx.subj
	post <- hap.obj$post
	nreps <- as.vector(hap.obj$nreps)
	uhap  <- sort(unique(c(hap1,hap2)))

	which.haplo <- hap.obj$hap.prob >= min.Hap.freq
	uhap <- uhap[which.haplo]
	x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")

	n.subj <- length(unique(hap.obj$indx.subj))
	n.x <- ncol(x)
	x.post<-matrix(rep(NA, n.subj * n.x), ncol=n.x)

	for(j in 1:n.x){
		x.post[,j] <- tapply(x[,j]*post, indx, sum)
	}
	x.post}
	
	x.post <- createPostHapMatrix(hap.obj, min.Hap.freq=min.Hap.freq)
    hap.index <- (1:length(hap.obj$hap.prob))[hap.obj$hap.prob >= min.Hap.freq]
	haplotypes <- hap.obj$haplotype[hap.index,]
	X.haps <- as.data.frame(x.post[,2:ncol(x.post)])
	names(X.haps) <- paste("Hap", hap.index[2:length(hap.index)], sep="")
	
	return(X.haps)}


