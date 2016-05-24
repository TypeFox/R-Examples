
MS <- function(GENO,niter=10,thetaID="user",params=FALSE,detail=FALSE,neutrality=FALSE, linkage = FALSE,F_ST = FALSE, MSMS = FALSE, big.data=FALSE){


## User can define some additional parameters
PAR <- FALSE
if(class(params)=="test.params"){
PAR <- TRUE
}

if(!(neutrality|linkage|F_ST)){
stop("You have to specify a statistic module. For example: (neutrality=TRUE)")
}


if(neutrality){
numTests    <- 11
testNames   <- c("Tajima.D","n.segregating.sites","Rozas.R_2","Fu.Li.F","Fu.Li.D","Fu.F_S","Strobeck.S","Fay.Wu.H","Zeng.E","theta_Tajima","theta_Watterson") 
F_ST        <- FALSE
linkage     <- FALSE

}

if(linkage){
numTests    <- 5
testNames   <- c("Wall.B","Wall.Q","Rozas.ZA","Rozas.ZZ","Kelly.Z_nS") 
neutrality  <- FALSE
F_ST        <- FALSE
}

if(F_ST){
numTests    <- 6
testNames   <- c("hap.diversity.within","Pi","haplotype.F_ST","nucleotide.F_ST","Nei.G_ST","Hudson.Snn")
neutrality  <- FALSE
linkage     <- FALSE
}


# --------------------------------------------
if(class(GENO)!="GENOME"){
   stop("Input is not a class of GENOME")
}
# ---------------------------------------------


if(neutrality){
 if(!GENO@Pop_Neutrality$calculated){
   stop("Neutrality statistics have to be calculated first !")
 }
}

if(linkage){
 if(!GENO@Pop_Linkage$calculated){
   stop("Linkage Disequilibrium statistics have to be calculated first !")
 }
 if(!GENO@Pop_Neutrality$calculated & !PAR){
  stop("Neutrality Module has to be calculated first! missing theta values")
 }
 
 if(!GENO@Pop_Neutrality$calculated & PAR){
  if(length(params@theta)==0){
  stop("Neutrality Module has to be calculated first! missing theta values")
  }
 }
 
}

if(F_ST){
  if(!GENO@Pop_FSTH$calculated){
   stop("F_ST statistics have to be calculated first !")
 }
  if(!GENO@Pop_Neutrality$calculated & !PAR){
  stop("Neutrality Module has to be calculated first! missing theta values")
 }
 if(!GENO@Pop_Neutrality$calculated & PAR){
  if(length(params@theta)==0){
  stop("Neutrality Module has to be calculated first! missing theta values")
  }
 }
}

if(!PAR){

if(linkage){
  if(length(GENO@Pop_Neutrality$Populations)!= length(GENO@Pop_Linkage$Populations) ){
  stop("You have to define the same Populations as in the Neutrality Module to calculate Linkage Disequilibrium, because the corresponding theta values are calculated in this module.")   
  }
}

if(F_ST){
  if(length(GENO@Pop_Neutrality$Populations)!= length(GENO@Pop_FSTH$Populations) ){
  stop("You have to define the same Populations as in the Neutrality Module to calculate F_ST statistics, because the corresponding theta values are calculated in this module.")   
  }
}

}

if(PAR){

if(length(params@theta)==0){

if(linkage){
  if(length(GENO@Pop_Neutrality$Populations)!= length(GENO@Pop_Linkage$Populations) ){
  stop("You have to define the same Populations as in the Neutrality Module to calculate Linkage Disequilibrium, because the corresponding theta values are calculated in this module.")   
  }
}

if(F_ST){
  if(length(GENO@Pop_Neutrality$Populations)!= length(GENO@Pop_FSTH$Populations) ){
  stop("You have to define the same Populations as in the Neutrality Module to calculate F_ST statistics, because the corresponding theta values are calculated in this module.")   
  }
}
}

}


if(thetaID=="user" & !PAR){
stop("You have to specify theta values in an class of test.params or define a thetaID")
}
if(thetaID=="user" &  PAR ){
  if(length(params@theta)==0){
      stop("You have to specify theta values in an class of test.params")
  }
  if(length(params@theta)!=GENO@genelength){
    stop(paste("The number of loci (", GENO@genelength ,") and the number of value of theta for each loci (", length(params@theta) ,") mismatch."))
  }
}


 # if there are no additional parameter
 nloci                         <-  GENO@genelength 
 locusData                     <-  NULL


### Rufe fuer jedes GEN MS einzelnd auf !

## PROGRESS #########################
 progr <- progressBar()
#####################################


if(PAR){
 if(length(params@migration)!=0){
 migg <- params@migration
}else{migg <- 1}

if(length(params@theta)!=0){
 THETA <- params@theta
 }
}

for(xx in 1:nloci){

if(!PAR) {params   <-  new("test.params")}   # <---------------------------------------- create an own test.params object

#if(length(GENO@region.data@biallelic.matrix[[xx]])!=0){

if(length(popGetBial(GENO,xx))!=0){

 if(neutrality){
 npops        <- length(GENO@region.stats@Pop_Neutrality[[xx]]$Populations)
 populations  <- GENO@region.stats@Pop_Neutrality[[xx]]$Populations
 }
 
 if(linkage){
 npops        <- length(GENO@region.stats@Pop_Linkage[[xx]]$Populations)
 populations  <- GENO@region.stats@Pop_Linkage[[xx]]$Populations
 }

 if(F_ST){
 npops        <- length(GENO@region.stats@Pop_FSTH[[xx]]$Populations)
 populations  <- GENO@region.stats@Pop_FSTH[[xx]]$Populations
 }
 


 ### make the migration parameter
   if(npops>1){
     migration <- vector(,npops)
      for(yy in 1:npops){
         migration[yy] <- length(populations[[yy]])
      }
   }else{migration <- numeric()}
 ################################
     

     params@n.pop      <-  npops
     params@n.loci     <-  1
 
 if(npops>1){
  if(!PAR & length(params@migration)==0){

  	if(MSMS[1]==FALSE){ 
   	   params@migration    <- c(migration,10)
	}else{
	   params@migration    <- c(migration,10)
	}

  }else{
    params@migration <- c(migration,migg)
  }
 }





  params@n.iter     <-  niter
 
  if(npops==1){

          
           if(thetaID=="Tajima")               {params@theta      <-  as.vector(GENO@theta_Tajima[xx,])}
           if(thetaID=="Watterson")            {params@theta      <-  as.vector(GENO@theta_Watterson[xx,])}
           if(thetaID=="user")                 {params@theta      <-  THETA[xx]}         
         
          
               #bial                <- GENO@region.data@biallelic.matrix[[xx]]
                bial                <- popGetBial(GENO,xx) 
               
               if(length(bial)>0){ 
                                
                if(neutrality) {WholePopulations     <- unlist(GENO@region.stats@Pop_Neutrality[[xx]]$Populations)}
                if(linkage)    {WholePopulations     <- unlist(GENO@region.stats@Pop_Linkage[[xx]]$Populations)}
                if(F_ST)       {WholePopulations     <- unlist(GENO@region.stats@Pop_FSTH[[xx]]$Populations)}
                
                 nsammy                              <- length(WholePopulations)

               }                                              
  }
  
## Change ( if there are more than one populations ) get the theta value of the whole Populations

  if(npops>1) 

  { 
       # bial <- GENO@region.data@biallelic.matrix[[xx]]
         bial                <- popGetBial(GENO,xx)   

        if(length(bial)>0){ 

           if(neutrality){WholePopulations  <- unlist(GENO@region.stats@Pop_Neutrality[[xx]]$Populations)}
           if(linkage)   {WholePopulations  <- unlist(GENO@region.stats@Pop_Linkage[[xx]]$Populations)}
           if(F_ST)      {WholePopulations  <- unlist(GENO@region.stats@Pop_FSTH[[xx]]$Populations)}

           nsammy            <- length(WholePopulations)
           erg               <- calc_freqstats(bial,list(WholePopulations))

           if(thetaID=="Tajima")               {thetanew          <- erg$THETA["thetaT",]}
           if(thetaID=="Watterson")            {thetanew          <- erg$THETA["thetaS",]}
           if(thetaID=="user")                 {thetanew          <- THETA[xx]}
        }

  params@theta <- thetanew

  }
 
 ### if there is theta == 0
 if(params@theta==0 | is.na(params@theta)){locusData <- c(locusData,NA);next}
  
  params@n.sam           <-  nsammy
  params@n.sites         <-  GENO@n.sites[xx]

  # Notloesung (aendern)
  # params@fixed.seg.sites <-  GENO@n.biallelic.sites[xx]
  #

  if(neutrality){params@obs.val      <-  getMS(GENO,xx,neutrality=TRUE)}
  if(linkage)   {params@obs.val      <-  getMS(GENO,xx,neutrality=FALSE,linkage=TRUE)}
  if(F_ST)      {params@obs.val      <-  getMS(GENO,xx,neutrality=FALSE,linkage=FALSE,F_ST=TRUE)}
  
  ## CALL MS with the parameters
     locusData         <- c(locusData,coalsimC(params,detail=detail,testNames=testNames,numTests=numTests,neutrality=neutrality,linkage=linkage,F_ST=F_ST,MSMS,big.data))   # list of locus DATA

}else{locusData        <- c(locusData,NA)}  # Wenn ueberhaupt Segregating Sites existieren


    if(nloci > 1000){
    cat(xx ," of ,", nloci, "\n" )
    }

# PROGRESS #######################################################
    progr <- progressBar(xx,nloci, progr)
##################################################################

} # End of ueber alle GENE

#####################################################################################################
# Prepare the summary class over all loci	
############################################
# Init	

if(neutrality){npop                  <- length(GENO@Pop_Neutrality$Populations)}
if(linkage)   {npop                  <- length(GENO@Pop_Linkage$Populations)}
if(F_ST)      {npop                  <- length(GENO@Pop_FSTH$Populations)}


popnames           <- paste("pop",1:npop)
init               <- vector("list",npop)
average            <- as.matrix(init)
variance           <- as.matrix(init)
rownames(average)  <- popnames
rownames(variance) <- popnames
colnames(variance) <- "variance"
colnames(average)  <- "average"
#-----------------------------


for(xx in 1:npop){
	
	# calc average and variance for all loci
	average[[xx]]           <- matrix(ncol=numTests, nrow=niter)
	colnames(average[[xx]]) <- testNames
	
	variance[[xx]]           <- matrix(ncol=numTests, nrow=niter)
	colnames(variance[[xx]]) <- testNames
		
	for (i in 1:niter) {

		for (j in 1:numTests) {
			tmp <- matrix(ncol=1)
			
			for (k in 1:nloci) {	
			  if(class(locusData[[k]])=="loc.stats"){
				   tmp <- rbind(tmp, locusData[[k]]@stats[[xx]][i,][j])
				}
			}

			average[[xx]][i,][j]  <- colMeans(tmp, na.rm=TRUE)
			variance[[xx]][i,][j] <- var(tmp, na.rm=TRUE)
		}
	}
	
} # End of for pops	
  # return(as.matrix(average))
	
        # aggregate all calculated stats in an csstats class
	# data of all each single locus is accumulated in a list of class locstats
	lociData <- new("cs.stats", 

		n.loci		= nloci,
		n.iter		= niter,
		n.pop		= npop,
		
		average		= average,
		variance	= variance,
		
		locus		= locusData
	)

	
#	return(lociData)
	
	init      <- as.matrix(vector("list",npop))
	probLess  <- init
	probEqual <- init
	validIter <- init
	rownames(probLess)   <- popnames
	rownames(probEqual)  <- popnames
	colnames(probLess)   <- "prob.less"
	colnames(probEqual)  <- "prob.equal"
	rownames(validIter)  <- popnames
	colnames(validIter)  <- "valid.iter"

	
	# calc probabilities of all loci if observed value has been given to test against
	
	
  for(xx in 1:npop){
	
		probLess[[xx]]           <- matrix(nrow=nloci, ncol=numTests)
		colnames(probLess[[xx]]) <- testNames
		
		probEqual[[xx]] 	  <- matrix(nrow=nloci, ncol=numTests)
		colnames(probEqual[[xx]]) <- testNames
		
		validIter[[xx]] 	  <- matrix(nrow=nloci, ncol=numTests)
		colnames(validIter[[xx]]) <- testNames
		
		
		# iterate through all locStats objects and create a matrix with all probabilities
		# this will provide a single matrix with average stats for all loci and is made accessible
		# through the csStats class
		
		# return(locusData)
		for (elem in 1:length(locusData)) {
		
		  if(class(locusData[[elem]])=="loc.stats"){
			probLess[[xx]][elem,]  = locusData[[elem]]@loc.prob.less[[xx]]
			probEqual[[xx]][elem,] = rbind(locusData[[elem]]@loc.prob.equal[[xx]])
			validIter[[xx]][elem,] = rbind(locusData[[elem]]@loc.valid.iter[[xx]])
		  }
		}
		
		#return(probLess)
		
  	# calculate average and variance values for probabilties
		averageL        <- colMeans(probLess[[xx]], na.rm=TRUE)
		varianceL       <- apply(probLess[[xx]], 2, var, na.rm=TRUE)
		probLess[[xx]]  <- rbind(probLess[[xx]], averageL, varianceL)
		
		# calculate average and variance values for probabilties
                averageE        <- colMeans(probEqual[[xx]], na.rm=TRUE)
		varianceE       <- apply(probEqual[[xx]], 2, var, na.rm=TRUE)
		probEqual[[xx]] <- rbind(probEqual[[xx]], averageE, varianceE)
		
		# calculate average and variance values for the number of valid iterations
                averageV        <- colMeans(validIter[[xx]], na.rm=TRUE)
		varianceV       <- apply(validIter[[xx]], 2, var, na.rm=TRUE)
		validIter[[xx]] <- rbind(validIter[[xx]], averageV, varianceV)
		
		# assign names to rows in format 'locusX', where X is 1,2,3....
		# for each considered loci
                rnames <- paste("locus", 1:nloci)
                rnames <- c(rnames,"average","variance")
                rownames(probEqual[[xx]]) <- rnames
                rownames(probLess[[xx]])  <- rnames
                rownames(validIter[[xx]]) <- rnames
		
	}# of for npops	
	  	 
		# assign all calculated matrices to the created csStats objects
		lociData@prob.less 	= probLess	
		lociData@prob.equal	= probEqual
		lociData@valid.iter	= validIter
 
		
	# add obsVal to csStats class if present
	#if ( !is.na(obsVal[1][1] )) {
		#colnames(obsVal) <- testNames
		if(neutrality){lociData@obs.val <- getMS(GENO)}
                 if(linkage){lociData@obs.val   <- getMS(GENO,linkage=TRUE,neutrality=FALSE)}
                 if(F_ST){lociData@obs.val      <- getMS(GENO,linkage=FALSE,neutrality=FALSE,F_ST=TRUE)}
                colnames(lociData@obs.val)      <- "obs.val"
	#}
		
	
	
	return(lociData)

}




