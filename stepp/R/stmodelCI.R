#################################################################
#
# stmodelCI.R
#
#######################
# stepp model: CI     #
#######################
#
cuminc.HR <- function (ftime, fstatus, group, strata, rho = 0, cencode = 0, 
    subset, na.action = na.omit) 
{
#
    d <- data.frame(time = ftime, cause = fstatus, 
				group = as.factor(
						if (missing(group)) rep(1, length(ftime)) else group),
				strata = as.factor(if (missing(strata)) rep(1, length(ftime)) else strata))
    if (!missing(subset)) d <- d[subset, ]
    tmp <- nrow(d)
    d <- na.action(d)
    if (nrow(d) != tmp) 
        cat(format(tmp - nrow(d)), "cases omitted due to missing values\n")
    no <- nrow(d)
    cg <- "  "
    nst <- length(levels(d$strata))
    d <- d[order(d$time), ]
    ugg <- table(d$group)
    d$group <- factor(d$group, names(ugg)[ugg > 0])
    ugg <- levels(d$group)
    censind <- ifelse(d$cause == cencode, 0, 1)
    uc <- table(d$cause[censind == 1])
    if (is.factor(d$cause)) 
        uclab <- names(uc)[uc > 0]
    else uclab <- as.numeric(names(uc)[uc > 0])
    nc <- length(uclab)
    ng <- length(ugg)
    if (ng > 1) {
        ng1 <- ng - 1
        ng2 <- ng * ng1/2
        v <- matrix(0, nrow = ng1, ncol = ng1)
        storage.mode(v) <- "double"
        vt <- double(ng2)
        s <- double(ng1)
    }
    pf <- vector("list", ng * nc)
    stat <- double(nc)
    l <- 0
    for (ii in 1:nc) {
        causeind <- ifelse(d$cause == uclab[ii], 1, 0)
        for (jj in 1:length(ugg)) {
            cg <- c(cg, paste(ugg[jj], uclab[ii]))
            l <- l + 1
            cgind <- d$group == ugg[jj]
            ncg <- length(cgind[cgind])
            n2 <- length(unique(d$time[cgind & causeind == 1]))
            n2 <- 2 * n2 + 2
            tmp <- double(n2)
            z <- .Fortran("cinc", as.double(d$time[cgind]), as.integer(censind[cgind]), 
                as.integer(causeind[cgind]), as.integer(ncg), 
                x = tmp, f = tmp, v = tmp)
#                x = tmp, f = tmp, v = tmp, PACKAGE = "cmprsk")
            pf[[l]] <- list(time = z$x, est = z$f, var = z$v)
        }
        if (ng > 1) {
            causeind <- 2 * censind - causeind
            z2 <- .Fortran("crstm", as.double(d$time), as.integer(causeind), 
                as.integer(d$group), as.integer(d$strata), as.integer(no), 
                as.double(rho), as.integer(nst), as.integer(ng), 
                s, v, as.double(d$time), as.integer(causeind), 
                as.integer(d$group), vt, s, vt, double((4 + 3 * 
                  ng) * ng), integer(4 * ng))
#                  ng) * ng), integer(4 * ng), PACKAGE = "cmprsk")
            stat[ii] <- -1
            a <- qr(z2[[10]])
            if (a$rank == ncol(a$qr)) {
                b <- diag(dim(a$qr)[1])
                stat[ii] <- z2[[9]] %*% qr.coef(a, b) %*% z2[[9]]
            }
#
            if (ii == 1) {
               ome <- (-1)*z2[[9]]
               omevar <- z2[[10]]
            }
        }
    }
    names(pf) <- cg[2:length(cg)]
    if (ng > 1) {
        names(stat) <- uclab
        stat <- list(Tests = cbind(stat = stat, pv = 1 - pchisq(stat, 
            ng - 1), df = rep(ng - 1, length(stat))))
        pf <- c(pf, stat)
#
        omeres <- list(ome=ome,omevar=omevar)
        pf <- c(pf,omeres)
    }
    pf
}


setClass("stmodelCI",
	   representation(coltrt    = "numeric",	# treatment
				coltime   = "numeric",  # time to event
				coltype   = "numeric",	# competing risk type
				trts      = "numeric",	# trt encoding
				timePoint = "numeric"   # evaluated time
				)
	   )


setMethod("estimate",	
	    signature="stmodelCI",
	    definition=function(.Object, sp, ...){
		nsubpop   <- sp@nsubpop
		subpop    <- sp@subpop
		coltrt    <- .Object@coltrt
		survTime  <- .Object@coltime
		trts      <- .Object@trts
		timePoint <- .Object@timePoint
		type      <- .Object@coltype
 
		txassign <- rep(NA, length(coltrt))
		txassign[which(coltrt == trts[1])] <- 1
		txassign[which(coltrt == trts[2])] <- 0

    		ObsCI1          <- rep(NA,nsubpop)
    		ObsCI2          <- rep(NA,nsubpop)
    		ObsCISE1        <- rep(NA,nsubpop)
    		ObsCISE2        <- rep(NA,nsubpop)
    		overallObsCI1   <- NA
    		overallObsCISE1 <- NA
    		overallObsCI2   <- NA
    		overallObsCISE2 <- NA
    		logHR           <- rep(0,nsubpop)
    		logHRSE         <- rep(0,nsubpop)

	    	for (i in 1:nsubpop) {
              result <- cuminc.HR(survTime[subpop[,i]==1],type[subpop[,i]==1],txassign[subpop[,i]==1])
              if (max(result$"1 1"$time) >= timePoint) {
                index       <- sum(result$"1 1"$time <= timePoint)
                ObsCI1[i]   <- result$"1 1"$est[index]
                ObsCISE1[i] <- sqrt(result$"1 1"$var[index])
              }
              if (max(result$"0 1"$time) >= timePoint) {
                index       <- sum(result$"0 1"$time <= timePoint)
                ObsCI2[i]   <- result$"0 1"$est[index]
                ObsCISE2[i] <- sqrt(result$"0 1"$var[index])
              }
              logHR[i]      <- result$ome/result$omevar
              logHRSE[i]    <- sqrt(1/result$omevar)
    		}
		logHRw <- sum(logHR/logHRSE)
		skmw   <- sum((ObsCI1-ObsCI2)/sqrt(ObsCISE1^2 + ObsCISE2^2))

    		result <- cuminc.HR(survTime,type,txassign)
    		if (max(result$"1 1"$time) >= timePoint) {
              index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
              overallObsCI1 <- result$"1 1"$est[index]
              overallObsCISE1 <- sqrt(result$"1 1"$var[index])
    		}
    		if (max(result$"0 1"$time) >= timePoint) {
              index <- sum(result$"0 1"$time <= timePoint)
              overallObsCI2 <- result$"0 1"$est[index]
              overallObsCISE2 <- sqrt(result$"0 1"$var[index])
      	}
    		overallLogHR <- result$ome/result$omevar
    		overallLogHRSE <- sqrt(1/result$omevar)
    		if (sum(is.na(ObsCI1)) != 0 | sum(is.na(ObsCI2)) != 0 | 
                is.na(overallObsCI1) != FALSE | is.na(overallObsCI2) != FALSE) {
              cat("\n")
              print(paste("Unable to estimate survival time at ", 
			  timePoint, " time-unit(s) because there are too few events within one or more subpopulation(s)."))
        	  print(paste("The problem may be avoided by constructing larger subpopulations and/or by selecting a different timepoint for estimation."))
        	  stop()
    		}
		#
		#   create the estimates (Cumulative Incidence) object - to be saved
		#
    		estimate <- list( model		= "CIe",
					sObs1       = ObsCI1,
        		            sSE1        = ObsCISE1,
                  	      oObs1     	= overallObsCI1,
                        	oSE1        = overallObsCISE1,
                        	sObs2       = ObsCI2,
                        	sSE2        = ObsCISE2,
                        	oObs2     	= overallObsCI2,
                        	oSE2        = overallObsCISE2,
					skmw		= skmw,

                        	logHR       = logHR,
                        	logHRSE     = logHRSE,
                        	ologHR      = overallLogHR,
                        	ologHRSE    = overallLogHRSE,
					logHRw	= logHRw
					)

		return(estimate)
	    }
)


setMethod("test",
	    signature="stmodelCI",
	    definition=function(.Object, nperm, sp, effect, showstatus=TRUE, ...){
		test <- NULL
		if (nperm > 0){
		  nsubpop	<- sp@nsubpop
		  subpop	<- sp@subpop
		  coltrt	<- .Object@coltrt
		  survTime  <- .Object@coltime
		  trts	<- .Object@trts
		  timePoint <- .Object@timePoint
		  type      <- .Object@coltype

		  txassign <- rep(NA, length(coltrt))
		  txassign[which(coltrt == trts[1])] <- 1
		  txassign[which(coltrt == trts[2])] <- 0
		  #
		  #   do the permutations
		  #

		  # set up the intial values
    		  pvalue      <- NA
      	  differences <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
		  diffha      <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      	  logHRs      <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      	  logHRha     <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      	  tPerm       <- rep(0, nperm)
      	  no          <- 0
      	  p           <- 0
      	  terminate   <- 0
      	  Ntemp       <- nrow(subpop)
      	  IndexSet1   <- (1:Ntemp)[txassign == 1]
      	  IndexSet2   <- (1:Ntemp)[txassign == 0]
	
	  	  if (showstatus){
		    title <- paste("\nComputing the pvalue with ", nperm)
		    title <- paste(title, "number of permutations\n")
		    cat(title)
	  	    pb <- txtProgressBar(min=0, max=nperm-1, style=3)
		  }

        	  while (no < nperm) {
	    	    ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
	    	    if (showstatus) setTxtProgressBar(pb, no)

          	    Subpop <- as.matrix(subpop)
          	    permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
          	    permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
          	    permuteSubpop[txassign == 0] <- Subpop[sample(IndexSet2),]
          	    subpop      <- permuteSubpop
          	    sCI1        <- rep(NA,nsubpop)
          	    sCI2        <- rep(NA,nsubpop)
          	    sCISE1      <- rep(NA,nsubpop)
          	    sCISE2      <- rep(NA,nsubpop)
          	    overallsCI1 <- NA
          	    overallsCI2 <- NA
          	    slogHRs     <- rep(0, nsubpop)
	    	    slogHRSE    <- rep(0, nsubpop)

          	    for (i in 1:nsubpop) {
                  result <- cuminc.HR(survTime[subpop[,i]==1],type[subpop[,i]==1],txassign[subpop[,i]==1])
                  if (max(result$"1 1"$time) >= timePoint) {
                    index     <- sum(result$"1 1"$time <= timePoint)
                    sCI1[i]   <- result$"1 1"$est[index]
                    sCISE1[i] <- sqrt(result$"1 1"$var[index])
                  }
                  if (max(result$"0 1"$time) >= timePoint) {
                    index     <- sum(result$"0 1"$time <= timePoint)
                    sCI2[i]   <- result$"0 1"$est[index]
                    sCISE2[i] <- sqrt(result$"0 1"$var[index])
                  }
                  slogHRs[i] <- result$ome/result$omevar
                  slogHRSE[i] <- sqrt(1/result$omevar)
          	    }
     	    	    slogHRw <- sum(slogHRs/slogHRSE)
	    	    skmwha  <- sum((sCI1-sCI2)/sqrt(sCISE1^2 + sCISE2^2))

          	    result <- cuminc.HR(survTime,type,txassign)
          	    if (max(result$"1 1"$time) >= timePoint) {
                  index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
                  overallsCI1 <- result$"1 1"$est[index]
          	    }
          	    if (max(result$"0 1"$time) >= timePoint) {
                  index <- sum(result$"0 1"$time <= timePoint)
                  overallsCI2 <- result$"0 1"$est[index]
          	    }
          	    overallSlogHR <- result$ome/result$omevar

          	    if (sum(is.na(sCI1)) == 0 & sum(is.na(sCI2)) == 0 & is.na(overallsCI1) == FALSE & is.na(overallsCI2) == FALSE) {
                  no <- no + 1
                  p <- p + 1
                  for (s in 1:nsubpop) {
                    differences[p,s] <- (sCI1[s]-sCI2[s])-(overallsCI1-overallsCI2)
			  diffha[p,s]      <- (sCI1[s]-sCI2[s])-skmwha
                    logHRs[p,s]      <- slogHRs[s]-overallSlogHR
			  logHRha[p,s]     <- slogHRs[s] - slogHRw
                  }
          	    }
          	    terminate <- terminate + 1
          	    if (terminate >= nperm + 10000) {
                  print(paste("After permuting ", nperm, "plus 10000, or ", 
                 	      nperm + 10000, " times, the program is unable to generate the permutation distribution based on ", 
                	      nperm, "permutations of the data"))
                  print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation"))
                  stop()
                }
              }
	  	  if (showstatus) close(pb)

	  	  # generating the sigmas and pvalues
	  	  sigmas        <- ssigma(differences)
	  	  sigma1        <- sigmas$sigma
	  	  chi2pvalue    <- ppv (differences, sigmas$sigmainv, effect$sObs1-effect$sObs2, effect$oObs1-effect$oObs2, nperm)
	  	  pvalue        <- ppv2(differences,                  effect$sObs1-effect$sObs2, effect$oObs1-effect$oObs2, nperm)

	  	  sigmas        <- ssigma(logHRs)
	  	  HRsigma       <- sigmas$sigma
	  	  #HRchi2pvalue <- ppv (logHRs,      sigmas$sigmainv, effect$logHR, effect$ologHR, nperm)
	  	  HRpvalue      <- ppv2(logHRs,                       effect$logHR, effect$ologHR, nperm)	

	  	  sigmas        <- ssigma(diffha)
	  	  hasigma       <- sigmas$sigma
  	  	  hapvalue      <- ppv (diffha,      sigmas$sigmainv, effect$sObs1-effect$sObs2, effect$skmw, nperm)

	  	  sigmas        <- ssigma(logHRha)
	  	  haHRsigma     <- sigmas$sigma
	  	  haHRpvalue    <- ppv (logHRha,     sigmas$sigmainv, effect$logHR, effect$logHRw, nperm)

	        test = list(model	 = "CIt",
		    	   	  sigma	 = sigma1,
			   	  hasigma	 = hasigma,
			   	  HRsigma    = HRsigma,
			   	  haHRsigma  = haHRsigma,
			   	  pvalue	 = pvalue,
				  chi2pvalue = chi2pvalue,
			   	  hapvalue   = hapvalue,
			   	  HRpvalue   = HRpvalue,
			   	  haHRpvalue = haHRpvalue
			    )
		}

	      return(test)
	    }
)

#
# printing support functions for cumulative incidence model
#
print.estimate.CI <- function(x, timePoint){
	cat("\n")
      #
      #   print out the cumulative incidence results 
      #
      write(paste("Cumulative incidence estimates for treatment group",x@model@trts[1],"at time point",timePoint),file="")
      temp <- matrix(c(1:x@subpop@nsubpop,round(x@effect$sObs1,digits=4),round(x@effect$sSE1,digits=4)),ncol=3)
      write("                        Cumulative",file="")
      write("     Subpopulation      Incidence        Std. Err.",file="")
      for (i in 1:x@subpop@nsubpop) {
        write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
              format(temp[i,3],width=15,nsmall=4)),file="")
      }
      write(paste("        Overall",format(round(x@effect$oObs1,digits=4),nsmall=4,width=16),
            format(round(x@effect$oSE1,digits=4),nsmall=4,width=15)),file="")
      cat("\n")
      write(paste("Cumulative incidence estimates for treatment group",x@model@trts[2],"at time point",timePoint),file="")
      temp <- matrix(c(1:x@subpop@nsubpop,round(x@effect$sObs2,digits=4),round(x@effect$sSE2,digits=4)),ncol=3)
      write("                        Cumulative",file="")
      write("     Subpopulation      Incidence        Std. Err.",file="")
      for (i in 1:x@subpop@nsubpop) {
        write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
              format(temp[i,3],width=15,nsmall=4)),file="")
      }
      write(paste("        Overall",format(round(x@effect$oObs2,digits=4),nsmall=4,width=16),
            format(round(x@effect$oSE2,digits=4),nsmall=4,width=15)),file="")
      cat("\n")
      write(paste("Cumulative incidence differences at time point",timePoint),file="")
      temp <- matrix(c(1:x@subpop@nsubpop,round(x@effect$sObs1-x@effect$sObs2,digits=4),round(sqrt(x@effect$sSE1^2+x@effect$sSE2^2),digits=4)),ncol=3)
      write("                        Cumulative",file="")
      write("                        Incidence",file="")
      write("     Subpopulation      Difference       Std. Err.",file="")
      for (i in 1:x@subpop@nsubpop) {
        write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
              format(temp[i,3],width=15,nsmall=4)),file="")
      }
      write(paste("        Overall",
            format(round(x@effect$oObs1-x@effect$oObs2,digits=4),nsmall=4,width=16),
            format(round(sqrt(x@effect$oSE1^2+x@effect$oSE2^2),digits=4),nsmall=4,width=15)),file="")
      cat("\n")
      write("Hazard ratio estimates",file="")
      temp <- matrix(c(1:x@subpop@nsubpop,round(x@effect$logHR,digits=6),round(x@effect$logHRSE,digits=6),round(exp(x@effect$logHR),digits=2)),ncol=4)
      write("     Subpopulation        Log HR       Std. Err.       Hazard Ratio",file="")
      for (i in 1:x@subpop@nsubpop) {
        write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=6),
              format(temp[i,3],width=14,nsmall=6),format(temp[i,4],width=15,nsmall=2)),file="")
      }
      write(paste("        Overall",format(round(x@effect$ologHR,digits=6),nsmall=6,width=16),
            format(round(x@effect$ologHRSE,digits=6),nsmall=6,width=14),
            format(round(exp(x@effect$ologHR),digits=2),nsmall=2,width=15)),file="")
      cat("\n")

}

print.cov.CI <- function(stobj, timePoint){
      cat("\n")
      write(paste("The covariance matrix of the cumulative incidence differences at", 
            timePoint, "time units for the", stobj@subpop@nsubpop, "subpopulations is:"), 
            file = "")
      print(stobj@result$sigma)

	cat("\n")
	write(paste("The covariance matrix of the log hazard ratios for the", 
          	stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$HRsigma)
      cat("\n")

      write(paste("The covariance matrix (based on the homogeneous association) of the cumulative incidence differences at", 
            timePoint, "time units for the", stobj@subpop@nsubpop, "subpopulations is:"), 
            file = "")
      print(stobj@result$hasigma)

	cat("\n")
	write(paste("The covariance matrix (based on the homogeneous association) of the log hazard ratios for the", 
          	stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$haHRsigma)
      cat("\n")

}

print.stat.CI <- function(pvalue, HRpvalue, chi2pvalue, hapvalue, haHRpvalue){
	cat("\n")
      write(paste("Supremum test results"), file = "")
      write(paste("Interaction P-value based on cumulative incidence estimates :", pvalue), file = "")
      write(paste("Interaction P-value based on hazard ratio estimates :", HRpvalue), file = "")

	cat("\n")
      write(paste("Chi-square test results"), file = "")
      write(paste("Interaction P-value based on cumulative incidence estimates :", 
            chi2pvalue), file = "")

	cat("\n")
      write(paste("Homogeneous association test results"), file = "")
      write(paste("Interaction P-value based on cumulative incidence estimates :", 
            hapvalue), file = "")
      write(paste("Interaction P-value based on hazard ratio estimates :", 
            haHRpvalue), file = "")

	cat("\n")
}

setMethod("print",
	    signature="stmodelCI",
	    definition=function(x, stobj, estimate=TRUE, cov=TRUE, test=TRUE, ...){

		#
		#  1. estimates
		#
	      if (estimate){
		  print.estimate.CI(stobj, x@timePoint)
      	}

  		#
		#   2. covariance matrices
		#
    		if (cov){
		  print.cov.CI(stobj, x@timePoint)
 		}
  
		#
		#   3. Supremum test and Chi-square test results
		#
		if (test){
 	        t <- stobj@result
		  print.stat.CI(t$pvalue, t$HRpvalue, t$chi2pvalue, t$hapvalue, t$haHRpvalue)
		}
 	 }
)

# constructor function for stmodelKM
stepp.CI <- function(coltrt, coltime, coltype, trts, timePoint){
	model <- new("stmodelCI", coltrt=coltrt, coltime=coltime, coltype=coltype,
			trts=trts, timePoint=timePoint)
	return(model)
}


