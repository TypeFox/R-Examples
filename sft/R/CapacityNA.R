capacityGroup <- function(inData, acc.cutoff=.9, ratio=TRUE, OR=NULL, stopping.rule=c("OR", "AND", "STST"), plotCt=TRUE, ...) {
  subjects <- sort(unique(inData$Subject))
  subjects <- factor(subjects)
  nsubjects <- length(subjects)

  conditions <- sort(unique(inData$Condition))
  conditions <- factor(conditions)
  nconditions <- length(conditions)

  channels <- grep("Channel", names(inData), value=T)
  nchannels <- length(channels)
  if(nchannels < 2) {
    stop("Not enough channels for capacity analysis.")
  }

  # For backward compatibility
  if (!is.null(OR)) {
    if (OR) {
      capacity <- capacity.or
    } else {
      capacity <- capacity.and 
    }
  } else {
    rule <- match.arg(stopping.rule, c("OR","AND","STST"))
    if(rule == "OR") {
      capacity <- capacity.or
    } else if (rule == "AND") {
      capacity <- capacity.and 
    } else if (rule == "STST"){
      capacity <- capacity.stst
    } else  {
      stop("Please choose a valid stopping rule for fPCAcapacity.")
    }
  }

  times <- seq(quantile(inData$RT,.001), quantile(inData$RT,.999), 
              length.out=1000)

  caplist <- vector("list")
  capmodel <- character()


  subj.out <- character()
  cond.out <- character()
  subj.out.g <- character()
  cond.out.g <- character()
  capMat <- numeric()
  if(!ratio) {
    varMat <- numeric()
  }

  devmat <- matrix(NA, nconditions, ceiling(nsubjects/9))

  RTlist <- vector("list", nchannels)
  CRlist <- vector("list", nchannels)

  for ( cn in 1:nconditions ) {
    Zscore <- numeric()
    m <- 0
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    #cond <- conditions[cn]
    condsubjects <- factor(with(inData, sort(unique(Subject[Condition==cond]))))
    ncondsubjects <- length(condsubjects)
    for ( sn in 1:ncondsubjects ) {
      if (is.factor(condsubjects)) {subj <- levels(condsubjects)[sn]} else {subj <- condsubjects[sn] }
      #subj <- condsubjects[sn]
      subj.out <- c(subj.out, subj)
      cond.out <- c(cond.out, cond)

      subj.out.g <- c(subj.out.g, subj)
      cond.out.g <- c(cond.out.g, cond)

      if ( sn %% 9 == 1 ) {
        m <- m+1
        if(plotCt) {
          dev.new()
          par(mfrow=c(3,3))
        }
      }

      ds <- inData$Subject==subj & inData$Condition==cond
      #if ( sum(ds) ==0 ) { next };
      good1 <- TRUE


      if(rule =="STST") {
        usechannel <- ds &  (apply(inData[,channels]>0, 1, sum)==1) & (apply(inData[,channels]<0, 1, sum)>0)
        RTlist[[1]] <- inData$RT[usechannel]
        CRlist[[1]] <- inData$Correct[usechannel]
        if(mean(CRlist[[1]]) < acc.cutoff | sum(CRlist[[1]]) < 10) {
          good1 <- FALSE 
        }

        usechannel <- ds & apply(inData[,channels]>=0, 1, all) & (apply(inData[,channels]!=0, 1, sum)==1)
        RTlist[[2]] <- inData$RT[usechannel]
        CRlist[[2]] <- inData$Correct[usechannel]
        if(mean(CRlist[[2]]) < acc.cutoff | sum(CRlist[[2]]) < 10) {
          good1 <- FALSE 
        }

      } else {
        usechannel <- ds & apply(inData[,channels]>0, 1, all)
        RTlist[[1]] <- inData$RT[usechannel]
        CRlist[[1]] <- inData$Correct[usechannel]
        if(mean(CRlist[[1]]) < acc.cutoff | sum(CRlist[[1]]) < 10) {
          good1 <- FALSE 
        }

        for ( ch in 1:nchannels ) {
          usechannel <- ds & inData[,channels[ch]]>0 & 
                        apply(as.matrix(inData[,channels[-ch]]==0), 1, all)
          RTlist[[ch+1]] <- inData$RT[usechannel]
          CRlist[[ch+1]] <- inData$Correct[usechannel]
          if(mean(CRlist[[ch+1]]) < acc.cutoff | sum(CRlist[[ch+1]]) < 10) {
              good1 <- FALSE
          }
        }
      }




      if( good1 ) {
        n <- length(subj.out)
        caplist[[n]] <- capacity(RT=RTlist, CR=CRlist, ratio=ratio)
        Zscore <- c(Zscore, caplist[[n]]$Ctest$statistic)
        if(caplist[[n]]$Ctest$p.value < .05) {
          if(caplist[[n]]$Ctest$statistic < 0) {
            capmodel <- c(capmodel, "Limited")
          } else {
            capmodel <- c(capmodel, "Super")
          }
        } else {
          capmodel <- c(capmodel, "Nonsignificant")
        }

        capMat <- rbind(capMat, caplist[[n]]$Ct(times))
        if(!ratio) {
          varMat <- rbind(varMat, caplist[[n]]$Var(times))
        }

        if(plotCt) {
          plot(times, tail(capMat,1), type='l',
              xlab="Time", ylab="C(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...)
          if(ratio) {
            abline(1,0, lwd=2)
          } else {
            lines(times, tail(capMat,1)+1.96*sqrt(tail(varMat,1)), lty=2)
            lines(times, tail(capMat,1)-1.96*sqrt(tail(varMat,1)), lty=2)
            abline(0,0, lwd=2)
          } 
          
        }

      } else {

        capmodel <- c(capmodel, NA)
        capMat <- rbind(capMat, rep(NA, length(times)) )
        if(!ratio){
          varMat <- rbind(varMat, rep(NA, length(times)) )
        }

        if(plotCt) {
          plot(c(min(times),max(times)), c(1,1), type='l', lwd=2,
              xlab="Time", ylab="C(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...) 
          text(mean(c(max(times),min(times))),1.2,"Not enough correct.",col='red')

        }
      }
    }

    if(plotCt) {
      dev.new()
      if(sum(cond.out==cond) > 1) {
        matplot(times, t(capMat[cond.out==cond,]), type='l', lty=1,
          main=paste(cond, paste(rule, "Capacity"),sep="\n"), xlab="Time",ylab="C(t)",...)
        if(ratio) {abline(1,0, lwd=2)} else{abline(0,0, lwd=2)} 

      } else {
        plot(times, capMat[cond.out==cond,], type='l', lty=1,
          main=paste(cond, paste(rule, "Capacity"),sep="\n"), xlab="Time",ylab="C(t)",...)
        if(ratio) {abline(1,0, lwd=2)} else{abline(0,0, lwd=2)} 
      }

    }

    subj.out.g <- c(subj.out.g, "Group")
    cond.out.g <- c(cond.out.g, cond)
    mZscore <- mean(Zscore, na.rm=TRUE)
    nZscore <- sum(!is.na(Zscore))
    pZscore <- t.test(Zscore)$p.value
    if( (pZscore < .025) | (pZscore > .975) )  {
      if(mZscore < 0) {
        capmodel <- c(capmodel, "Limited")
      } else {
        capmodel <- c(capmodel, "Super")
      }
    } else {
      capmodel <- c(capmodel, "Nonsignificant")
    }

  }

  overview <- as.data.frame(list(Subject=subj.out.g, Condition=cond.out.g,
      Capacity=capmodel))

  if(ratio){
    #return(list(statistic=Z, Ct.or=capORMat, Ct.and=capANDMat, times=times))
    return(list(overview=overview, Ct.fn=capMat, Ct.fn=capMat, capacity=caplist, times=times))
  } else {
    return(list(overview=overview, Ct.fn=capMat, Ct.var=varMat, capacity=caplist, times=times))
    #return(list(statistic=Z, Ct.or=capORMat, Var.or=varORMat, Ct.and=capANDMat, Var.and=varANDMat, times=times))
  }

}


capacity.or <- function(RT, CR=NULL, ratio=TRUE) {
    if ( is.null(CR) | (length(CR) != length(RT)) ) {
      CR <- vector("list", length(RT))
      for( i in 1:length(RT) ) {
        CR[[i]] <- rep(1, length(RT[[i]]))
      }
    } 
    times <- sort(unique(c(RT, recursive=TRUE))) 
    ncond <- length(RT) - 1 

    # Find Nelson-Aalen Cumulative Hazard Estimates
    numer <- estimateNAH(RT=RT[[1]], CR=CR[[1]])
    denom <- estimateUCIPor(RT=RT[1+(1:ncond)], CR=CR[1+(1:ncond)])

    rmtest <- ucip.test(RT, CR, OR=TRUE)

    if (ratio) {
      C.or <- numer$H(times) / denom$H(times)

      C.or[is.nan(C.or)] <- NA
      C.or[is.infinite(C.or)] <- NA
      C.or <- approxfun(times, C.or)
      return( list(Ct=C.or, Ctest=rmtest) )
    } else {
      C.or <- numer$H(times) - denom$H(times)
      C.or <- approxfun(c(0,times), c(0,C.or))
      Var.or <- numer$Var(times) + denom$Var(times)
      Var.or <- approxfun(c(0,times), c(0,Var.or))
      return( list(Ct=C.or, Var=Var.or, Ctest=rmtest, p.val=rmtest$p.val) )
    }
}


capacity.and <- function(RT, CR=NULL, ratio=TRUE) {
    if ( is.null(CR) | (length(CR) != length(RT)) ) {
      CR <- vector("list", length(RT))
      for( i in 1:length(RT) ) {
        CR[[i]] <- rep(1, length(RT[[i]]))
      }
    } 
    times <- sort(unique(c(RT, recursive=TRUE))) 

    ncond <- length(RT) - 1 

    rmtest <- ucip.test(RT, CR, OR=FALSE)

    # Find Nelson-Aalen Reverse Cumulative Hazard Estimates
    denom <- estimateNAK(RT[[1]], CR[[1]])
    numer <- estimateUCIPand(RT=RT[1+(1:ncond)], CR=CR[1+(1:ncond)])

    # Calculate the and capacity coefficient
    if (ratio) {
      C.and <- numer$K(times) / denom$K(times)
      C.and[is.nan(C.and)] <- NA
      C.and[is.infinite(C.and)] <- NA
      C.and <- approxfun(times, C.and)
      return( list(Ct=C.and, Ctest=rmtest) )
    } else {
      C.and <- denom$K(times) - numer$K(times) 
      C.and <- approxfun(c(times,Inf), c(C.and,0))
      Var.and <- numer$Var(times) + denom$Var(times)
      Var.and <- approxfun(c(times,Inf), c(Var.and,0))
      return( list(Ct=C.and, Var=Var.and, Ctest=rmtest) )
    }
}
