#################################################################
#
# stsubpop.R
#
#######################
# stepp subpopulation #
#######################
setClass("stsubpop",
	   representation(win      = "stwin",	# stepp window type
				colvar   = "numeric",   # vector of covariate of interest (V) 
				nsubpop  = "numeric",	# number of subpopulation generated
				subpop   = "ANY",   	# matrix of subpopulations
				npatsub  = "numeric",	# count of each subpopulation
				medianz  = "numeric",	# median of V for each subpopulation
				minz     = "numeric",	# minimum of V for each subpopulation
				maxz	   = "numeric",	# maximum of V for each subpopulation
				init	   = "logical"	# initialized
				)
	   )

setMethod("initialize", "stsubpop",
	function(.Object){
		.Object@init <- FALSE
		return(.Object)
	}
)


setValidity("stsubpop", 
	function(object){
	  if (!is.numeric(object@colvar) || length(object@colvar) < object@win@r2){
		print("Invalid argument")
		return (FALSE)
	  }
	  return(TRUE)
     }
)

setGeneric("generate", function(.Object, win, covariate)
		standardGeneric("generate"))	

setMethod("generate", 
	    signature="stsubpop",
	    definition=function(.Object, win, covariate){
	  	r1 <- win@r1
	  	r2 <- win@r2

	  	# get the overlapping subgroups
	  	#
	  	zvals <- rep(0, length(covariate))
	  	absfreq <- rep(0, length(covariate))
	  	sortedz <- sort(covariate)
	  	zvals[1] <- sortedz[1]
	  	j <- 1
	  	absfreq[1] <- 1
	  	for (i in 2:length(covariate)) {
          	  if (sortedz[i] != zvals[j]) {
                j <- j + 1
                zvals[j] <- sortedz[i]
                absfreq[j] <- 1
          	  }
          	else {
              absfreq[j] <- absfreq[j] + 1
          	}
	    }
	    zvals <- zvals[1:j]
	    absfreq <- absfreq[1:j]
	    cumfreq <- absfreq
	    for (i in 2:length(cumfreq)) cumfreq[i] <- cumfreq[i] + cumfreq[i - 1]
	      I0 <- rep(0, 1000)	# max size of I0 and I1 is 1000; these limit the 
	      I1 <- rep(0, 1000)	# max no. of stepp subpopulation to 1000 also.
	      I0[1] <- 1
	      I1[1] <- sum(cumfreq < r2) + 1
	      stopflag <- 0
	      nsubpop <- 2
	      while (stopflag == 0) {
              indinf <- I0[nsubpop - 1] + 1
              while ((cumfreq[I1[nsubpop - 1]] - cumfreq[indinf - 1]) > r1) {
                indinf <- indinf + 1
              }
  	        I0[nsubpop] <- indinf
  	        indsup <- I1[nsubpop - 1]
  	        while (((cumfreq[indsup] - cumfreq[I0[nsubpop]] + absfreq[I0[nsubpop]]) < r2) && (stopflag == 0)) {
                indsup <- indsup + 1
                stopflag <- 1 * (indsup == length(zvals))
              }
            I1[nsubpop] <- indsup
            nsubpop <- nsubpop + 1
	      }
	    nsubpop <- nsubpop - 1
	    npatsub <- rep(0, nsubpop)
	    minz <- rep(NA, nsubpop)
	    maxz <- rep(NA, nsubpop)
	    npatsub[1] <- cumfreq[I1[1]]
	    for (i in 2:nsubpop) npatsub[i] <- cumfreq[I1[i]] - cumfreq[I0[i] - 1]
	    I0 <- I0[1:nsubpop]
	    I1 <- I1[1:nsubpop]
	    npats <- length(covariate)
	    subpop <- matrix(rep(0, (npats * nsubpop)), ncol = nsubpop)
	    medians <- rep(0, nsubpop)
	    for (i in 1:nsubpop) {
            subpop[, i] <- (covariate >= zvals[I0[i]]) * (covariate <= zvals[I1[i]])
            medians[i] <- round((median(covariate[subpop[, i] == 1])), digits = 2)
            minz[i] <- round(zvals[I0[i]],digits=4)
            maxz[i] <- round(zvals[I1[i]],digits=4)
	    }

	  # update the object
	  .Object@win	<- win
	  .Object@colvar  <- covariate
	  .Object@nsubpop <- nsubpop
	  .Object@subpop  <- subpop
	  .Object@npatsub <- npatsub
	  .Object@medianz <- medians
	  .Object@minz    <- minz
	  .Object@maxz    <- maxz
	  .Object@init    <- TRUE

	  return(.Object)
	  }
)


setMethod("summary", 
	    signature="stsubpop",
	    definition=function(object){
		summary(object@win)
		if (object@init){
		  write(paste("      Number Of Subpopulations Created :", object@nsubpop),file="")
    		  cat("\n")
    		  write("Subpopulation Summary Information",file="")
    		  nper <- apply(object@subpop,2,sum)
    		  temp <- matrix(c(1:object@nsubpop,object@medianz,object@minz,object@maxz,nper),ncol=5)
    		  write("                                  Covariate Summary                  Sample",file="")
    		  write("     Subpopulation        Median       Minimum       Maximum          Size",file="")
    		  for (i in 1:object@nsubpop) {
       	    write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=2),
          	    format(temp[i,3],width=13,nsmall=4),format(temp[i,4],width=13,nsmall=4),
          	    format(temp[i,5],width=13)),file="")
    		  }
		}
		else
		  write("      Subpopulations have not been generated", file="")
	    }	
)

# constructor function for stepp subpopulation
stepp.subpop <- function(swin, cov){
	subp    <- new("stsubpop")		
	subp    <- generate(subp, win=swin, covariate=cov)
	return(subp)
}

