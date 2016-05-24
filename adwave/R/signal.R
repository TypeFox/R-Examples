signal <-
function(table, who=colnames(table), populations,popA=NA, popB=NA,  # for working out group centroids
		 normalize=FALSE,n.pca=5,
		 PCAonly=FALSE, verbose=TRUE,tol=0.001,#tolerance for normalization
              n.signal=NULL,window.size=NULL,  #optionally can window signal 
		 genmap= NULL){ # if supplied, formulate signals in terms of genetic distance along chromosome
  
  #table -- nsnps*people matrix
  #who -- individuals to include for
  #populations -- list containing vector of ID's for each popultion in the analysis
  #popA, popB -- names of ancestral populations (used for forming the axes of variation). Must match one of the names of populations.
  #normalize -- whether to normalize the data matrix. Default is FALSE	
  #no.axes -- number of pca axes to compute (only first PC is used for forming the signals but more may be wanted for visualisation)
  #tol -- tolerance for normalization of admixture signals (see paper)
 #PCAonly -- if true, only compute the PCA, do not compute the signals.
	#optional arguments for windowing the signal:
  #n.signal	 -- number of data points in windowed signal
  #window.size -- size of window specified as proportion of total length. 
	#e.g window.size= 0.002 with signal of length 13000 snps generates window of 0.002*13000=26 SNPs. Does not need to be a round number. 
  #genmap = NULL



  # 1. DATA CHECKING
  #consistancy of input data
  if(sum(is.na(populations[[popA]]))>0 | sum(is.na(populations[[popB]]))>0) stop("Check pop1 pop2. Exiting....")

  if(!is.null(n.signal) & is.null(window.size)) stop("Must supply window.size. Exiting....")	
  if(!is.null(genmap) & is.null(window.size)) stop("Must supply window.size. Exiting....")	
  if(!is.null(genmap) & (length(genmap) != dim(table)[1])) stop("Number of SNPs and genetic map do not match. Exiting....")	


  # Prepare input variables
  individuals <- colnames(table);names(individuals) <- individuals
  axis.who <- individuals[c(populations[[popA]], populations[[popB]])]
  who <- individuals[who]
  all <- unique(c(who,axis.who))
  if(length(all)<length(individuals)) table <- table[,all] 
  snps <- rownames(table)
  n.snps   <- length(snps)

  
  # prepare return object
  result <- list()
  class(result) <- "adsig"
  result$call <- sys.call()
  result$date <- date()
  result$individuals <- who
  result$popA <- individuals[c(populations[[popA]])]
  result$popB <- individuals[c(populations[[popB]])]
	
  result$n.snps <- n.snps
  if(is.null(n.signal)){ result$window.size <- 1}else{ result$window.size <- window.size*n.snps}
  result$pa.ind  <- array(dim=c(length(individuals),n.pca),data=NA, dimnames=list(individuals,1:n.pca))
  result$pa.snp   <-array(dim=c(length(snps),n.pca),dimnames=list(snps,1:n.pca), data=NA)
  result$pc.ind   <- array(dim=c(length(individuals),n.pca),dimnames=list(individuals,1:n.pca),data=NA)
  
  
  if(verbose){
    cat("Number of individuals:\t", length(who),"\n")
    cat("Number of individuals for finding axes:\t", length(axis.who),"\n")
    cat("Number of SNPs:\t\t",   n.snps,"\n")
      }
  
  #DATA PREPARATION	
  #Centering the table- subtract row means so that for a given location along the genome the mean off all snps is zero
  table.mod <- t(apply(table,1,function(x) return(x-mean(x[axis.who],na.rm=TRUE))))
  # Normalize 
  if(normalize){  table.mod <- t(apply(table.mod,1,function(x) return(x/sqrt(var(x[axis.who],na.rm=TRUE)))))  }      
  

	# PCA TO FIND AXES
	result$G <- var(table.mod[,axis.who],na.rm=TRUE)
	spec       <- eigen(result$G,symmetric=TRUE)
	result$ev  <- spec$values # eigenvalues. save for all pca so we can return proportion of explained variance
	if(length(spec$values)<n.pca){n.pca <- length(spec$values)}
	result$pa.ind <- spec$vectors[,1:n.pca] #eigenvectors. store just the most informative, up to n.pca
  
  
  # Find eigenvectors in "snp coords"
  # replace na's by 0 (uninformative)
  tbl.a <- table.mod[,axis.who]
  tbl.a[which(is.na(tbl.a))] <- 0
  for(i in 1:n.pca){
    result$pa.snp[,i] <- apply(tbl.a, 1, function(x) return(sum(x*result$pa.ind[,i],na.rm=TRUE)))
    result$pa.snp[,i] <- result$pa.snp[,i]/sqrt(sum(result$pa.snp[,i]^2))
  }
  rm(tbl.a)
  
  # Find principal components in "individuals coordinates"
  for(i in 1:n.pca){
    result$pc.ind[,i] <- apply(table.mod[,who], 2, function(x) return(sum(x*result$pa.snp[,i],na.rm=TRUE)))/sqrt(n.snps)
  }
  



#Estimating proportion of admixture
#Population means
popnames <- names(populations)
popP <- rep(NA, length(popnames));names(popP)<- popnames
cA <-  mean(result$pc.ind[populations[[popA]],1]) 
cB <-  mean(result$pc.ind[populations[[popB]],1]) 
for(p in popnames){
popP[p]<- (cB- mean(result$pc.ind[populations[[p]],1]) )/(cB-cA)
}
result$popP <- popP
#Individual level
result$indP <- (cB - result$pc.ind[,1])/(cB-cA)

	
	#SIGNAL CREATION
  	if(PCAonly == FALSE){ 

	
	tempsig <- table #start with RAW original input table
	tempsig <- t(apply(tempsig,1,function(x) return(x-mean(x[axis.who],na.rm=TRUE)))) #center
	tempsig[is.na(tempsig)] <- 0 #missing data
	result$sig <- sweep(tempsig,MARGIN=1,result$pa.snp[,1],`*`) #RAW DATA * LOADING
	
	#OPTIONAL WINDOWING 
	if(!is.null(n.signal)){
	#Define filter for windowing
    	func_filt <- function(x){ 
	ws1 <- ceiling(result$n.snps*window.size)
	ws2 <- floor(result$n.snps*window.size)
	filt <- rep(1, ws1)
	if(ws1>ws2){filt[ws1] <- (result$n.snps*window.size - ws2) }
	filt <- filt/sum(filt) 
	filter(x, filter = filt, sides=1)} #define filter
	#take average in a window
	tempsigA <- apply(result$sig,2,FUN = func_filt) # average    
	
	# OPTION A: Even (SNP) sampling
	if(is.null(genmap)){
	tempsigA <- tempsigA[(ceiling(result$n.snps/n.signal)-1):result$n.snps,] #crop to get rid of NA's
	#define function for sampling
	func_interp<- function(x){ res <- approx(1:dim(tempsigA)[1], x,  method="constant", n=n.signal)
                               return(res$y) }
	#Final sampled signal
	result$sig <- apply(t(tempsigA), c(1), FUN=func_interp) 
			}

	if(!is.null(genmap)){
	#OPTION B: Interpolate to genetic distance
	tempsigA[1:(floor(n.snps*window.size)-1),]  <-   rep(tempsigA[(floor(n.snps*window.size)),],each=(floor(n.snps*window.size)-1)) # replicate to deal with NA's
	#Then interpolate to genetic distance
   	 func_interp <- function(x){ res <- approx(genmap, x,  method="linear", n=n.signal)
                               return(res$y) }
    
    result$sig <- apply(t(tempsigA), c(1), FUN=func_interp)  
    #physical positions of interpolated waveform points
    result$gendist <- approx(genmap, tempsigA[,1], method="linear", n=n.signal)$x
    result$step  <- mean(diff(approx(genmap, tempsigA[,1], method="linear", n=n.signal)$x))
				}

						} #END OF SIGNAL CREATION



  #2. NORMALISE - each individual SNP
  fun_normalise <- function(x){(2*x-(UPPER+LOWER))/(UPPER-LOWER)} #define function for nomalisation
  UPPER= rowMeans(result$sig[,populations[[popA]]],na.rm=TRUE)
  LOWER= rowMeans(result$sig[,populations[[popB]]],na.rm=TRUE)
  result$signals <- (apply(result$sig, 2, FUN= fun_normalise))
  diff <- abs(UPPER-LOWER)
  result$signals[diff <= tol] <- 0	
  result$n.tol <- sum(diff < tol)
	
		} # End of signal creation

	result$sig <- NULL # remove non-normalised version to save space

 

  if(verbose){ cat("Finished calculations","\n") }
  

  return(result)
}
