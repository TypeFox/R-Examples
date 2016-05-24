#################################################################
#
# stmodelKM.R
#
#######################
# stepp model: KM     #
#######################
#
kmest1<-function(y,m,n,ndf,t,s,v,ntpt,tpt,nrr,ndd)
{
  # initialization 
  f<-1.0
  kr<-n
  nrr[1]<-n
  ndd[1]<-0

  for (i in 2:ntpt){
    nrr[i]<-0
    ndd[i]<-0
  }

  ltp<-1
  var<-0.0
  l<-1
  t[1]<-0
  s[1]<-1
  v[1]<-0
  i<-1

  # main loop 
  while (i<= n) {
    k<-i+1
    k2<-0

    while (k2<=0) {
      if (k > n) k2<-1
      else {
	if (y[k] != y[i]) k2<-1 else k<-k+1
      }
	
    } # end while 

    k<-k-1
    nd<-0

    for (j in i:k) nd<-nd+m[j]

    while (ltp<=ntpt && y[i]>tpt[ltp+1]) {
      ltp<-ltp+1
      nrr[ltp]<-kr
    }  # end while 

    ndd[ltp]<-ndd[ltp]+nd
    if (nd>0) {
      t1<- nd / kr  
      f<-f*(1-t1)
      if (nd<kr) var<-var+t1/(kr-nd)
      t[l+1]<-y[i]
      s[l+1]<-s[l]
      v[l+1]<-v[l]
      l<-l+2
      t[l]<-y[i]
      s[l]<-f
      v[l]<-var*f*f
    } # end if 

    i<-k+1
    kr<-n-k     # kr<-n-k-1
    k<-i
  } # end while - main loop 

  l<-l+1
  t[l]<-y[n]
  s[l]<-s[l-1]
  v[l]<-v[l-1]

  # return a list of all the arguments
  kmtest1<-list(y=y,m=m,n=n,ndf=ndf,t=t,s=s,v=v,ntpt=ntpt,tpt=tpt,nrr=as.integer(nrr),ndd=as.integer(ndd))
  kmtest1
  
}


kmest <-
function (time,status,group,tpt,pv=TRUE,pv.strat,pv.sub,rho=0,subset,na.action=na.omit) 
{
  d <- data.frame(time=time,status=status,
    group=as.factor(if (missing(group)) rep(1,length(time)) else group),
    pv.strat=as.factor(if (missing(pv.strat)) rep(1,length(time)) else pv.strat),
    pv.sub = as.factor(if (missing(pv.sub)) rep(TRUE, length(time)) else pv.sub))
  if (!missing(subset)) d <- d[subset,]
  tmp <- nrow(d)
  d <- na.action(d)
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')
  no <- nrow(d)
  if (any(d$status != 0 & d$status != 1)) stop("invalid status values")
  T1 <- table(d$group)
  subgl <- T1>0
  T8 <- d$group[d$status==1]
  if (length(T8) > 0) {
    T8 <- cbind(T1,table(d$group[d$status==1]))
  } else {
    T8 <- cbind(T1,rep(0,length(T8)))
  }
  T8 <- T8[subgl,,drop=FALSE]
  T1 <- names(T1)[subgl]
  time <- ifelse(d$time<=0,.00001,d$time)
  if (missing(tpt)) {
    tpt <- pretty(time)}
  else {
    ym <- round(max(time),2)
    tpt <- c(0,tpt,ym)
  }
  ntpt <- length(tpt)-1
  lev <- vector("character",ntpt)
  for (i in 1:ntpt) lev[i] <- paste(format(tpt[i]),format(tpt[i+1]),sep="-")
  nrd <- rep(0, ntpt)
  o <- order(time)
  Tl <- length(T1)
  z <- as.list(1:Tl)
  for(i in 1:Tl) {
    Ty <- (time[o])[d$group[o]==T1[i]]
    Tm <- (d$status[o])[d$group[o]==T1[i]]
    ndf <- length(unique(Ty[Tm==1]))
    t <- double(2*ndf+2)

    a <- kmest1(Ty, Tm, length(Ty), ndf, t,t,t,ntpt, tpt, nrd, nrd)

    tt <- paste(a[[11]],a[[10]],sep="/")
    names(tt) <- lev
    z[[i]] <- list(time=a[[5]],est=a[[6]],var=a[[7]],tint=tt,nnd=T8[i,])
    }
  names(z) <- T1
  if (pv & Tl > 1) {
    stop("internal error")
    # pv <- logrank(time=d$time,status=d$status,group=d$group,strata=d$pv.strat,rho=rho,subset=(d$pv.sub=='TRUE'))$pv
    # attr(z,"pv") <- pv
  }

  class(z) <- "kmest"
  z
}

#
tpest1<-function(x,n,ind,tp5,ntp)
{
  l <- ntp
  for (i in rev(1:ntp)){
    if (x[n] >= tp5[i]) break;
    ind[l]<-0
    l <- l - 1
  }

  if (l <=0) return (ind);

  if (x[n] == tp5[l]){
    ind[l] <- n
    l <- l - 1
  }

  # assuming unique values sorted in ascending order
  k <- n-1
  loop <- TRUE
  while(loop){
    if (l <= 0) return (ind);

   loop <- FALSE
    for (i in rev(1:k)) {
      if (x[k] <= tp5[l]){
        ind[l] <- k+1
        l <- l - 1
        loop <- TRUE
        break #out of the for loop
      } else 
      k <- k-1
    } # end for loop
  } # end while loop

  # error in the following loop corrected 9-28-04
  for (i in 1:l) ind[i] <- 0
  return (ind);	
}


tpest <-
function(w,times) 
{
  if (!is.null(w$Tests)) w <- w[names(w) != 'Tests']
  ng <- length(w)
  times <- sort(unique(times))
  nt <- length(times)
  storage.mode(times) <- "double"
  storage.mode(nt) <- "integer"
  ind <- matrix(0,ncol=nt,nrow=ng)
  oute <- matrix(NA,ncol=nt,nrow=ng)
  outv <- oute
  storage.mode(ind) <- "integer"
  slct <- rep(TRUE,ng)
  for (i in 1:ng) {
    if (is.null((w[[i]])$est)) { slct[i] <- FALSE} else { 
       z1 <- as.integer(tpest1(w[[i]][[i]], length(w[[i]][[i]]), ind[i,], times, nt))
       ind[i,] <- z1
       oute[i, ind[i,]>0] <- w[[i]][[2]][z1]
       if (length(w[[i]])>2) outv[i,ind[i,]>0] <- w[[i]][[3]][z1]
	}
    }
  
  dimnames(oute) <- list(names(w)[1:ng],as.character(times))
  dimnames(outv) <- dimnames(oute)
  list(est=oute[slct,,drop=FALSE],var=outv[slct,,drop=FALSE])
}

setClass("stmodelKM",
	   representation(coltrt   = "numeric",	# treatment
				survTime = "numeric",   # time to event
				censor   = "numeric",	# time to censor
				trts     = "numeric",	# trt encoding
				timePoint= "numeric"	# evaluated time
				)
	   )

setMethod("estimate",	
	    signature="stmodelKM",
	    definition=function(.Object, sp, ...){

		nsubpop   <- sp@nsubpop
		subpop    <- sp@subpop
		coltrt    <- .Object@coltrt
		trts	    <- .Object@trts
		timePoint <- .Object@timePoint
		survTime  <- .Object@survTime
		censor    <- .Object@censor
 
		txassign <- rep(NA, length(coltrt))
		txassign[which(coltrt == trts[1])] <- 1
		txassign[which(coltrt == trts[2])] <- 0

		skmObs1   <- rep(0, nsubpop)
    		skmObs2   <- rep(0, nsubpop)
    		skmSE1    <- rep(0, nsubpop)
    		skmSE2    <- rep(0, nsubpop)
    		logHR     <- rep(0, nsubpop)
    		logHRSE   <- rep(0, nsubpop)

    		for (i in 1:nsubpop) {
		  trteff1    <- tpest(kmest(survTime[txassign==1 & subpop[,i]==1],
              		          censor[txassign==1 & subpop[,i]==1]),timePoint)
		  trteff2    <- tpest(kmest(survTime[txassign==0 & subpop[,i]==1],
              		          censor[txassign==0 & subpop[,i]==1]),timePoint)
          	  skmObs1[i] <- max(trteff1$est,0)
          	  skmObs2[i] <- max(trteff2$est,0)
          	  skmSE1[i]  <- sqrt(trteff1$var)
          	  skmSE2[i]  <- sqrt(trteff2$var)

          	  LogRank    <- survdiff(Surv(survTime[subpop[,i]==1],censor[subpop[,i]==1])~txassign[subpop[,i]==1])
          	  logHR[i]   <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
          	  logHRSE[i] <- sqrt(1/LogRank$var[1,1])
		}
		overalltrteff1 <- tpest(kmest(survTime[txassign==1],censor[txassign==1]),timePoint)
		overalltrteff2 <- tpest(kmest(survTime[txassign==0],censor[txassign==0]),timePoint)
    		overallSkmObs1 <- max(overalltrteff1$est, 0)
    		overallSkmObs2 <- max(overalltrteff2$est, 0)
    		overallSkmSE1  <- sqrt(overalltrteff1$var)
    		overallSkmSE2  <- sqrt(overalltrteff2$var)
    		LogRank        <- survdiff(Surv(survTime,censor) ~ txassign)
    		overallLogHR   <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
    		overallLogHRSE <- sqrt(1/LogRank$var[1,1])

		# check to make sure that subpopulation estimates are OK.
    		if (sum(is.na(skmObs1)) != 0 | sum(is.na(skmObs2)) != 0 | 
              is.na(overallSkmObs1) != FALSE | is.na(overallSkmObs2) != FALSE) {
          	  cat("\n")
        	  print(paste("Unable to estimate survival time at ", 
                          timePoint, " time-unit(s) because there are too few events within one or more subpopulation(s)."))
        	  print(paste("The problem may be avoided by constructing larger subpopulations and/or by selecting a different timepoint for estimation."))
        	  stop()
    		}

		skmw   <- sum((skmObs1-skmObs2)/sqrt(skmSE1^2 + skmSE2^2))
		logHRw <- sum(logHR/logHRSE)

    		estimate <- list( model	   = "KMe",
					sObs1	   = skmObs1,
        	            	sSE1     = skmSE1,
                        	oObs1    = overallSkmObs1,
                        	oSE1	   = overallSkmSE1,
                        	sObs2	   = skmObs2,
                        	sSE2	   = skmSE2,
                        	oObs2	   = overallSkmObs2,
                        	oSE2	   = overallSkmSE2,
					skmw	   = skmw,

					logHR    = logHR,
					logHRSE  = logHRSE,
					ologHR   = overallLogHR,
					ologHRSE = overallLogHRSE,
					logHRw   = logHRw
				     )
		return(estimate)
	    }
)


setMethod("test",
	    signature="stmodelKM",
	    definition=function(.Object, nperm=2500, sp, effect, showstatus=TRUE, Cox=FALSE, MM=NULL, ...){
		test <- NULL
		if (nperm > 0){
		  nsubpop	<- sp@nsubpop
		  subpop	<- sp@subpop
		  coltrt	<- .Object@coltrt
		  trts	<- .Object@trts
		  survTime  <- .Object@survTime
		  timePoint <- .Object@timePoint
		  censor    <- .Object@censor

		  txassign <- rep(NA, length(coltrt))
		  txassign[which(coltrt == trts[1])] <- 1
		  txassign[which(coltrt == trts[2])] <- 0
		  #
		  #   do the permutations
		  #
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

		  if (Cox) {
		    if (is.null(MM)){
		      fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign"))
		    } else {
		      xnam <- colnames(MM)[2:dim(MM)[2]]
		      fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign + ", paste(xnam, collapse= "+")))
		    }
		  }

		  if (showstatus){
		    title <- paste("\nComputing the pvalue with ", nperm)
		    title <- paste(title, "number of permutations\n")
		    cat(title)
	  	    pb <- txtProgressBar(min=0, max=nperm-1, style=3)
		  }

        	  while (no < nperm) {
		    ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
      	    if (showstatus) setTxtProgressBar(pb, no)

      	    Subpop        <- as.matrix(subpop)
      	    permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
      	    permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
      	    permuteSubpop[txassign == 0] <- Subpop[sample(IndexSet2),]
      	    subpop        <- permuteSubpop
      	    skm1          <- rep(0, nsubpop)
      	    skm2          <- rep(0, nsubpop)
		    skmSE1	      <- rep(0, nsubpop)
		    skmSE2	      <- rep(0, nsubpop)
      	    slogHRs       <- rep(0, nsubpop)
		    slogHRSEs     <- rep(0, nsubpop)

      	    for (i in 1:nsubpop) {
			seff1		<- tpest(kmest(survTime[txassign == 1 & subpop[, i] == 1], censor[txassign == 1 & subpop[, i] == 1]), timePoint)
      		seff2		<- tpest(kmest(survTime[txassign == 0 & subpop[, i] == 1], censor[txassign == 0 & subpop[, i] == 1]), timePoint) 
			skm1[i]     <- max(seff1$est, 0)
        		skm2[i]     <- max(seff2$est, 0)
          		skmSE1[i]   <- sqrt(seff1$var)
          		skmSE2[i]   <- sqrt(seff2$var)
			if (Cox) {
	        	  CoxM       	<- coxph(fmla, subset=(subpop[,i]==1))
		  	  if (is.null(MM)){
	                slogHRs[i]    <- CoxM$coefficient
	                slogHRSEs[i]  <- sqrt(CoxM$var)
		  	  } else {
	                slogHRs[i]    <- CoxM$coefficient[1]
	                slogHRSEs[i]  <- sqrt(CoxM$var[1,1])  
			  }
			} else {
        		  LogRank         <- survdiff(Surv(survTime[subpop[,i] == 1],censor[subpop[,i] == 1]) ~ txassign[subpop[,i] == 1])
        		  slogHRs[i]      <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
			  slogHRSEs[i]    <- sqrt(1/LogRank$var[1,1])
			}

      	    }
		    skmwha               <- sum((skm1-skm2)/sqrt(skmSE1^2 + skmSE2^2))
		    slogHRw		       <- sum(slogHRs/slogHRSEs)

      	    overallSkm1          <- max(tpest(kmest(survTime[txassign == 1], censor[txassign == 1]), timePoint)$est, 0)
      	    overallSkm2          <- max(tpest(kmest(survTime[txassign == 0], censor[txassign == 0]), timePoint)$est, 0)
		    if (Cox) {
			CoxM		       <- coxph(fmla)
		  	if (is.null(MM)){
		    	  overallSLogHR    <- CoxM$coefficient
		        overallLogHRSE   <- sqrt(CoxM$var)
		  	} else {
		    	  overallSLogHR    <- CoxM$coefficient[1]
		    	  overallSLogHRSE  <- sqrt(CoxM$var[1,1])
			}
		    } else {
      	      LogRank            <- survdiff(Surv(survTime,censor) ~ txassign)
      	      overallSLogHR      <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
    		      overallSLogHRSE    <- sqrt(1/LogRank$var[1,1])
		    }
      	    if (sum(is.na(skm1)) == 0 & sum(is.na(skm2)) == 0 & is.na(overallSkm1) == FALSE & is.na(overallSkm2) == FALSE) {
        		no <- no + 1
        		p <- p + 1
        		for (s in 1:nsubpop) {
          	  	  differences[p, s] <- (skm1[s] - skm2[s]) - (overallSkm1 - overallSkm2)
			  diffha[p,s]	  <- (skm1[s] - skm2[s]) - skmwha
          	  	  logHRs[p,s]       <- slogHRs[s] - overallSLogHR
			  logHRha[p,s]      <- slogHRs[s] - slogHRw

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
	        sigmas          <- ssigma(differences)
	        sigma1		<- sigmas$sigma
	        chi2pvalue   	<- ppv (differences, sigmas$sigmainv, effect$sObs1-effect$sObs2, effect$oObs1-effect$oObs2, nperm)
	        pvalue     	<- ppv2(differences,                  effect$sObs1-effect$sObs2, effect$oObs1-effect$oObs2, nperm)

	        sigmas          <- ssigma(diffha)
	        hasigma      	<- sigmas$sigma
  	        hapvalue     	<- ppv (diffha, sigmas$sigmainv, effect$sObs1-effect$sObs2, effect$skmw, nperm)

	  	  sigmas          <- ssigma(logHRs)
	  	  HRsigma    	<- sigmas$sigma
	  	  HRpvalue   	<- ppv2(logHRs, effect$logHR, effect$ologHR, nperm)	

	  	  sigmas          <- ssigma(logHRha)
	  	  haHRsigma    	<- sigmas$sigma
	  	  haHRpvalue   	<- ppv (logHRha, sigmas$sigmainv, effect$logHR, effect$logHRw, nperm)

	        test = list( model	  = "KMt",
				   sigma	  = sigma1,
				   hasigma	  = hasigma,
				   HRsigma    = HRsigma,
				   haHRsigma  = haHRsigma,
				   pvalue	  = pvalue,
				   chi2pvalue = chi2pvalue,
				   hapvalue	  = hapvalue,
				   HRpvalue   = HRpvalue,
				   haHRpvalue = haHRpvalue
				  )
		}

		return(test)
	    }
)

#
# printing support functions for KM model
#
print.estimate.KM <- function(x, timePoint){
	cat("\n")
      write(paste("Survival estimates for treatment group", x@model@trts[1], 
            	"at time point", timePoint), file = "")
        	temp <- matrix(c(1:x@subpop@nsubpop, round(x@effect$sObs1, digits = 4), round(x@effect$sSE1, 
            	digits = 4)), ncol = 3)
        	write("                         Survival", file = "")
        	write("     Subpopulation     Probability      Std. Err.", 
            	file = "")
        	for (i in 1:x@subpop@nsubpop) {
              write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 4), format(temp[i, 3], width = 15, 
                nsmall = 4)), file = "")
        	}
        	write(paste("        Overall", format(round(x@effect$oObs1, 
            	digits = 4), nsmall = 4, width = 16), format(round(x@effect$oSE1, 
            	digits = 4), nsmall = 4, width = 15)), file = "")
        	cat("\n")
        	write(paste("Survival estimates for treatment group", x@model@trts[2], 
            	"at time point", timePoint), file = "")
        	temp <- matrix(c(1:x@subpop@nsubpop, round(x@effect$sObs2, digits = 4), round(x@effect$sSE2, 
            	digits = 4)), ncol = 3)
        	write("                         Survival", file = "")
        	write("     Subpopulation     Probability      Std. Err.", 
            	file = "")
        	for (i in 1:x@subpop@nsubpop) {
              write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 4), format(temp[i, 3], width = 15, 
                nsmall = 4)), file = "")
        	}
        	write(paste("        Overall", format(round(x@effect$oObs2, 
            	digits = 4), nsmall = 4, width = 16), format(round(x@effect$oSE2, 
            	digits = 4), nsmall = 4, width = 15)), file = "")
        	cat("\n")
        	write(paste("Survival differences at time point", timePoint), 
            	file = "")
        	temp <- matrix(c(1:x@subpop@nsubpop, round(x@effect$sObs1 - x@effect$sObs2, digits = 4), 
          	round(sqrt(x@effect$sSE1^2 + x@effect$sSE2^2), digits = 4)), ncol = 3)
        	write("                         Survival", file = "")
        	write("     Subpopulation      Difference      Std. Err.", 
            	file = "")
        	for (i in 1:x@subpop@nsubpop) {
              write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 4), format(temp[i, 3], width = 15, 
                nsmall = 4)), file = "")
        	}
        	write(paste("        Overall", format(round(x@effect$oObs1 - 
            	x@effect$oObs2, digits = 4), nsmall = 4, width = 16), 
            	format(round(sqrt(x@effect$oSE1^2 + x@effect$oSE2^2), 
                	digits = 4), nsmall = 4, width = 15)), file = "")
        	cat("\n")
        	write("Hazard ratio estimates", file = "")
        	temp <- matrix(c(1:x@subpop@nsubpop, round(x@effect$logHR, digits = 6), round(x@effect$logHRSE, 
            	digits = 6), round(exp(x@effect$logHR), digits = 2)), ncol = 4)
        	write("     Subpopulation        Log HR       Std. Err.       Hazard Ratio", 
            	file = "")
        	for (i in 1:x@subpop@nsubpop) {
              write(paste(format(temp[i, 1], width = 12),
				  format(temp[i, 2], width = 19, nsmall = 6),
				  format(temp[i, 3], width = 14, nsmall = 6),
				  format(temp[i, 4], width = 15, nsmall = 2)), 
                		  file = "")
        	}
        	write(paste("        Overall",
			format(round(x@effect$ologHR,      digits = 6), nsmall = 6, width = 16),
			format(round(x@effect$ologHRSE,    digits = 6), nsmall = 6, width = 14),
			format(round(exp(x@effect$ologHR), digits = 2), nsmall = 2, width = 15)),
			file = "")
        	cat("\n")

}

print.cov.KM <- function(stobj, timePoint){
      cat("\n")
	write(paste("The covariance matrix of the Kaplan-Meier differences at", 
           	timePoint, "time units for the", stobj@subpop@nsubpop, "subpopulations is:"), 
           	file = "")
	print(stobj@result$sigma)
      
	cat("\n")
	write(paste("The covariance matrix of the log hazard ratios for the", 
          	stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$HRsigma)
      cat("\n")

	write(paste("The covariance matrix (based on homogeneous association) of the Kaplan-Meier differences at", 
           	timePoint, "time units for the", stobj@subpop@nsubpop, "subpopulations is:"), 
           	file = "")
	print(stobj@result$hasigma)
      
	cat("\n")
	write(paste("The covariance matrix (based on homogeneous association) of the log hazard ratios for the", 
          	stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$haHRsigma)
      cat("\n")

}

print.stat.KM <- function(pvalue, HRpvalue, chi2pvalue, hapvalue, haHRpvalue){
	cat("\n")
      write(paste("Supremum test results"), file = "")
     	write(paste("Interaction P-value based on Kaplan-Meier estimates :", pvalue), file = "")
	write(paste("Interaction P-value based on hazard ratio estimates :", HRpvalue), file = "")

	cat("\n")
      write(paste("Chi-square test results"), file = "")
      write(paste("Interaction P-value based on Kaplan-Meier estimates :", 
     	      chi2pvalue), file = "")

	cat("\n")
      write(paste("Homogeneous association test results"), file = "")
      write(paste("Interaction P-value based on Kaplan-Meier estimates :", 
      	hapvalue), file = "")
      write(paste("Interaction P-value based on hazard ratio estimates :", 
      	haHRpvalue), file = "")

	cat("\n")
}

setMethod("print",
	    signature="stmodelKM",
	    definition=function(x, stobj, estimate=TRUE, cov=TRUE, test=TRUE, ...){

		#
		#  1. estimates
		#
	      if (estimate){
		  print.estimate.KM(stobj, x@timePoint)
      	}

  		#
		#   2. covariance matrices
		#
    		if (cov){
		  print.cov.KM(stobj, x@timePoint)
 		}
  
		#
		#   3. Supremum test and Chi-square test results
		#
		if (test){
 	        t <- stobj@result
		  print.stat.KM(t$pvalue, t$HRpvalue, t$chi2pvalue, t$hapvalue, t$haHRpvalue)
		}
 	 }
)

# constructor function for stmodelKM
stepp.KM <- function(coltrt, survTime, censor, trts, timePoint){
	model <- new("stmodelKM", coltrt=coltrt, survTime=survTime, censor=censor,
			trts=trts, timePoint=timePoint)
	return(model)
}

