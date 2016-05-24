assessmentGroup <- function(inData, stopping.rule=c("OR", "AND"), correct=c(TRUE, FALSE), fast=c(TRUE,FALSE), detection=TRUE, plotAt=TRUE, ...) {
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

  rule <- match.arg(stopping.rule, c("OR","AND"))
  if(rule == "OR") {
    capacity <- capacity.or
  } else if (rule == "AND") {
    capacity <- capacity.and 
  } else  {
    stop("Please choose a valid stopping rule for assessmentGroup.")
  }
  

  times <- seq(quantile(inData$RT,.001), quantile(inData$RT,.999), 
              length.out=1000)

  atlist <- vector("list")


  subj.out <- character()
  cond.out <- character()
  subj.out.g <- character()
  cond.out.g <- character()
  atMat <- numeric()

  devmat <- matrix(NA, nconditions, ceiling(nsubjects/9))

  RTlist <- vector("list", nchannels)
  CRlist <- vector("list", nchannels)

  for ( cn in 1:nconditions ) {
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
        if(plotAt) {
          dev.new()
          par(mfrow=c(3,3))
        }
      }

      ds <- inData$Subject==subj & inData$Condition==cond
      #if ( sum(ds) ==0 ) { next };
      good1 <- TRUE


      usechannel <- ds & apply(inData[,channels]>0, 1, all)
      RTlist[[1]] <- inData$RT[usechannel]
      CRlist[[1]] <- inData$Correct[usechannel]
      if(sum(CRlist[[1]]) < 10) {
        good1 <- FALSE 
      }

      for ( ch in 1:nchannels ) {
        usechannel <- ds & inData[,channels[ch]]>0 & 
                      apply(as.matrix(inData[,channels[-ch]]==0), 1, all)
        RTlist[[ch+1]] <- inData$RT[usechannel]
        CRlist[[ch+1]] <- inData$Correct[usechannel]
        if(sum(CRlist[[ch+1]]) < 10) {
            good1 <- FALSE
        }
      }


      if( good1 ) {
        n <- length(subj.out)
        atlist[[n]] <- assessment(RT=RTlist, CR=CRlist, stopping.rule=rule, correct=correct, fast=fast, detection=detection )
        atMat <- rbind(atMat, atlist[[n]](times))

        if(plotAt) {
          plot(times, tail(atMat,1), type='l',
              xlab="Time", ylab="A(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...)
          abline(1,0, lwd=2)
        }

      } else {

        atMat <- rbind(atMat, rep(NA, length(times)) )

        if(plotAt) {
          plot(c(min(times),max(times)), c(1,1), type='l', lwd=2,
              xlab="Time", ylab="A(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...) 
          text(mean(c(max(times),min(times))),1.2,"Not enough correct.",col='red')

        }
      }
    }

    if(plotAt) {
      dev.new()
      if(sum(cond.out==cond) > 1) {
        matplot(times, t(atMat[cond.out==cond,]), type='l', lty=1,
          main=paste(cond, paste(rule, "Capacity"),sep="\n"), xlab="Time",ylab="C(t)",...)
        abline(1,0, lwd=2)

      } else {
        plot(times, atMat[cond.out==cond,], type='l', lty=1,
          main=paste(cond, paste(rule, "Capacity"),sep="\n"), xlab="Time",ylab="C(t)",...)
        abline(1,0, lwd=2)
      }

    }


  }


  return(list(At.fn=atMat, assessment=atlist, times=times))

}

assessment <- function(RT, CR, OR=NULL, stopping.rule=c("OR","AND"), correct=c(TRUE, FALSE), fast=c(TRUE, FALSE), detection=TRUE) {
  if (is.null(OR)) {
    rule <- match.arg(stopping.rule, c("OR","AND"))
  } else if (OR ==TRUE) {
    rule <- "OR"
  } else {
    rule <- "AND"
  }

  slow <- !fast
  incorrect <- !correct

  times <- sort(unique(unlist(RT)))
  p.correct <- unlist(lapply(CR, mean))
  
  if (rule=="OR") {
    if (detection) {
      if (correct & fast) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])

        channel1C  <- ecdf.channel1(times) * p.correct[2] * (1-p.correct[3])
        channel2C  <- ecdf.channel2(times) * p.correct[3] * (1-p.correct[2])
        channel12C <- ecdf.channel1(times) * p.correct[2] * (1-ecdf.channel2(times)) * p.correct[3]
        channel21C <- (1-ecdf.channel1(times)) * p.correct[2] * ecdf.channel2(times) * p.correct[3]
        channelCC <- ecdf.channel1(times) * p.correct[2] * ecdf.channel2(times) * p.correct[3]

        numer <- log(channel1C + channel2C + channel12C + channel21C + channelCC)
        denom <- log(ecdf.redundant(times) * p.correct[1])

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Correct and Fast"
      }
      if (correct & slow) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])
        
        channel1C  <- (1-ecdf.channel1(times)) * p.correct[2] * (1-p.correct[3])
        channel2C  <- (1-ecdf.channel2(times)) * p.correct[3] * (1-p.correct[2])
        channelCC <-  (1-ecdf.channel1(times)) * p.correct[2] * (1-ecdf.channel2(times)) * p.correct[3]
        
        numer <- log(channel1C + channel2C + channelCC)
        denom <- log( (1-ecdf.redundant(times)) * p.correct[1] )

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Correct and Slow"
      }
      if (incorrect & fast ) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

        numer <- log( ecdf.channel1(times) * (1-p.correct[2])) + log( ecdf.channel2(times) * (1-p.correct[3]))
        denom <- log( ecdf.redundant(times) * (1-p.correct[1]) )

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Incorrect and Fast"
      }
      if (incorrect & slow) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

        channel1I  <- (1-ecdf.channel1(times)) * (1-p.correct[2]) * (1-p.correct[3])
        channel2I  <- (1-ecdf.channel2(times)) * (1-p.correct[2]) * (1-p.correct[3])
        channelII <-  (1-ecdf.channel1(times)) * (1-ecdf.channel2(times)) * (1-p.correct[2]) * (1-p.correct[3])

        numer <- log( channel1I + channel2I - channelII )
        denom <- log( (1-ecdf.redundant(times)) * (1-p.correct[1]) )

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Incorrect and Slow"
      }
    } else {
      if (correct & fast) {
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==1])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==1] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- cumsum(g)
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.correct[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log(ecdf.redundant(times)*p.correct[1])
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Correct and Fast"
      }
  
      if (correct & slow) {
      
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==1])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==1] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- rev(cumsum(rev(g)))
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.correct[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log( (1-ecdf.redundant(times))*p.correct[1] )
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Correct and Slow"
      }
      
      if (incorrect & fast ) {
        p.incorrect <- 1-unlist(lapply(CR, mean))
      
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==0])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==0] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(1-CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- cumsum(g)
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.incorrect[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log( ecdf.redundant(times)*p.incorrect[1] )
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Incorrect and Fast"
      }
      
      if (incorrect & slow) {
        p.incorrect <- 1-unlist(lapply(CR, mean))
      
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==0])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==0] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(1-CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- rev(cumsum(rev(g)))
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.incorrect[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log( (1-ecdf.redundant(times))*p.incorrect[1] )
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Incorrect and Slow"
      }
    }
  } else if (rule=="AND") {
    if (correct & fast) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])

      numer <- log( ecdf.channel1(times) * p.correct[2] ) + log( ecdf.channel2(times) * p.correct[3] ) 
      denom <- log( ecdf.redundant(times) * p.correct[1] ) 

      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Correct and Fast"
    }
    if (correct & slow) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])

      channel1C  <- (1-ecdf.channel1(times)) * p.correct[2] * p.correct[3]
      channel2C  <- (1-ecdf.channel2(times)) * p.correct[3] * p.correct[2]
      channelCC <-  (1-ecdf.channel1(times)) * p.correct[2] * (1-ecdf.channel2(times)) * p.correct[3]

      numer <- log ( channel1C + channel2C - channelCC ) 
      denom <- log ( (1-ecdf.redundant(times)) * p.correct[1] ) 

      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Correct and Slow"
    }
    if (incorrect & fast ) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

      channel1I  <- ecdf.channel1(times) * (1-p.correct[2]) * p.correct[3]
      channel2I  <- ecdf.channel2(times) * (1-p.correct[3]) * p.correct[2]
      channel12I <- ecdf.channel1(times) * (1-p.correct[2]) * (1-ecdf.channel2(times)) * (1-p.correct[3])
      channel21I <- (1-ecdf.channel1(times)) * (1-p.correct[2]) * ecdf.channel2(times) * (1-p.correct[3])
      channelII <- ecdf.channel1(times) * (1-p.correct[2]) * ecdf.channel2(times) * (1-p.correct[3])

      numer <- log( channel1I + channel2I + channel12I + channel21I + channelII)
      denom <- log( ecdf.redundant(times) * (1-p.correct[1]) )

      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Incorrect and Fast"

    }
    if (incorrect & slow) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

      channel1I  <- (1-ecdf.channel1(times)) * (1-p.correct[2]) * p.correct[3]
      channel2I  <- (1-ecdf.channel2(times)) * (1-p.correct[3]) * p.correct[2]
      channelII <-  (1-ecdf.channel1(times)) * (1-p.correct[2]) * (1-ecdf.channel2(times)) * (1-p.correct[3])

      numer <- log( channel1I + channel2I + channelII )
      denom <- log( (1-ecdf.redundant(times)) * (1-p.correct[1]) )
      
      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Incorrect and Slow"

    }
  }
  return(A)
}
