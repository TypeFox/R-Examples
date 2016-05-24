sicGroup <- function(inData, sictest="ks", domtest="ks", plotSIC=TRUE, ...) {

  subjects <- sort(unique(inData$Subject))
  subjects <- factor(subjects)
  nsubjects <- length(subjects)

  conditions <- sort(unique(inData$Condition))
  conditions <- factor(conditions)
  nconditions <- length(conditions)

  SICnames <- c("ParallelOR", "ParallelAND", "SerialOR", "SerialAND", "Coactive")

  times <- sort(unique(round(inData$RT)))

  sicList <- vector("list")

  subj.out <- character()
  cond.out <- character()
  siresult <- character()
  positive <- character()
  negative <- character()
  mic <- character()
  modelresult <- character()
  predicted <- character()
  rejected <- character()

  allSICfn <- numeric()

  n <- 0

  if( sictest!="ks") {
    cat("Only KS-SIC test is currently implemented.\n")
    return(NA)
  }
  
  for ( cn in 1:nconditions ) {
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    for ( sn in 1:nsubjects ) {
      if (is.factor(subjects)) {subj <- levels(subjects)[sn]} else {subj <- subjects[sn] }
      if (plotSIC & ( sn %% 9 == 1) ) {
        dev.new()
        par(mfrow=c(3,3))
      }

      HH <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==2 & Channel2==2] )
      HL <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==2 & Channel2==1] )
      LH <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==1 & Channel2==2] )
      LL <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==1 & Channel2==1] )

      if ( min( length(HH), length(HL), length(LH), length(LL)) > 10 ) {
        cond.out <- c(cond.out, cond)
        subj.out <- c(subj.out, subj)
        n <- n+1
  
        sicList[[n]] <- sic(HH=HH, HL=HL, LH=LH, LL=LL)
        allSICfn <- rbind(allSICfn, sicList[[n]]$SIC(times))


        if( all(sicList[[n]]$Dominance$p.value[1:4] < .05) & !any(sicList[[n]]$Dominance$p.value[5:8] < .05) ) {
          siresult <- c(siresult, "Pass")
        } else if ( any(sicList[[n]]$Dominance$p.value[5:8] < .05) ) {
          siresult <- c(siresult, "Fail")
        } else {
          siresult <- c(siresult, "Ambiguous")
        }

        rejected.models  <- rep(FALSE,5)

        if (sicList[[n]]$SICtest$positive$p.value < .05) {
          rejected.models[c(2,3)] <- TRUE  # ParallelAND, SerialOR
          positive <- c(positive, "Significant")
        } else {positive <- c(positive, "Nonsignificant")}

        if (sicList[[n]]$SICtest$negative$p.value < .05) {
          rejected.models[c(1,3)] <- TRUE # ParallelOR, SerialOR
          negative <- c(negative, "Significant")
        } else {negative <- c(negative, "Nonsignificant")}

        if (sicList[[n]]$MIC$p.value < .05) {
          rejected.models[c(3,4)] <- TRUE # SerialOR, SerialAND
          if (sicList[[n]]$MIC$statistic > 0) {
            mic <- c(mic, "Positive")
          } else { mic <- c(mic, "Negative") }
        } else {mic <- c(mic,"Nonsignificant")}

        if( sum(rejected.models)==0) {
          predicted <- c(predicted, SICnames[3]) # Serial OR
        } else if( sum(rejected.models)==1) {
          predicted <- c(predicted, SICnames[5]) # Coactive
        } else if( all(rejected.models[1:3] == c(0,1,1) ) ) {
          predicted <- c(predicted, SICnames[1]) # Parallel-OR
        } else if( all(rejected.models[1:3] == c(1,0,1) ) ){
          predicted <- c(predicted, SICnames[2]) # Parallel-AND
        } else if( all(rejected.models[1:4] == c(1,1,1,0) ) ) {
          predicted <- c(predicted, SICnames[4]) # Serial-AND
        } else {
          predicted <- c(predicted, "NA")
        }

        rejected <- c(rejected, paste(SICnames[rejected.models],collapse=",") )
        if(plotSIC) {
          plot(times, sicList[[n]]$SIC(times), type='l',
            main=paste(cond, " Condition\nParticipant ", subj, sep=""), 
            xlab="Time",ylab="SIC(t)",...)
        }
      }
    }

    if(plotSIC) {
      dev.new()
      matplot(times, t(allSICfn[cond.out==cond,]),type='l',lty=1,
        main=paste(cond, " Condition", sep=""), 
        xlab="Time",ylab="SIC(t)",...)
    }
  }

  overview <- as.data.frame(list(Subject=subj.out, Condition=cond.out,
      Selective.Influence=siresult, Positive.SIC=positive, Negative.SIC=negative, 
      MIC=mic, Predicted_by=predicted, Rejected.Models=rejected))
  return(list(overview=overview, sic.fn=allSICfn, sic=sicList, times=times))
}


sic <- function(HH, HL, LH, LL, domtest="ks", sictest="ks", mictest=c("art", "anova")) {
    mictest <- match.arg(mictest)
    RTall <- sort(unique(c(HH, HL, LH, LL)))
    HH.ecdf <- ecdf(HH)
    HL.ecdf <- ecdf(HL)
    LH.ecdf <- ecdf(LH)
    LL.ecdf <- ecdf(LL)

    dominance <- siDominance(HH, HL, LH, LL, method=domtest)

    N<-1/length(HH)+1/length(HL)+1/length(LH)+1/length(LL)
    N <- 1/N

    sicall <- LH.ecdf(RTall) + HL.ecdf(RTall) - HH.ecdf(RTall) - LL.ecdf(RTall)
    SIC <- stepfun(RTall, c(0,sicall))

    sicstat <- sic.test(HH, HL, LH, LL, method=sictest)
    micstat <- mic.test(HH, HL, LH, LL, method=mictest)

    return(list(SIC=SIC, Dominance=dominance, SICtest=sicstat, MICtest=micstat, N=N))
}

sic.test <- function(HH, HL, LH, LL, method="ks") {
  DNAME <- paste("\nHH:", deparse(substitute(HH)), "\tHL:", deparse(substitute(HL)), 
                 "\nLH:", deparse(substitute(LH)), "\tLL:", deparse(substitute(LL)) )
  METHOD <- "Houpt-Townsend KS-SIC test"
  RTall <- sort(unique(c(HH, HL, LH, LL)))
  HH.ecdf <- ecdf(HH)
  HL.ecdf <- ecdf(HL)
  LH.ecdf <- ecdf(LH)
  LL.ecdf <- ecdf(LL)


  N<-1/length(HH)+1/length(HL)+1/length(LH)+1/length(LL)
  N <- 1/N

  sicall <- LH.ecdf(RTall) + HL.ecdf(RTall) - HH.ecdf(RTall) - LL.ecdf(RTall)

  if (method=="ks") {
    Dplus  <- max(0,sicall)
    p.Dplus  <- exp(-2 * N * Dplus ^2 )
    names(Dplus) <- "D^+"
    rp <- list(statistic=Dplus, p.value=p.Dplus, alternative="the SIC is above 0 at some time",
      method=METHOD, data.name=DNAME)
    class(rp) <- "htest"

    Dminus <- abs(min(0,sicall))
    p.Dminus <- exp(-2 * N * Dminus ^2 )
    names(Dminus) <- "D^-"
    rn <- list(statistic=Dminus, p.value=p.Dminus, alternative="the SIC is below 0 at some time",
      method=METHOD, data.name=DNAME)
    class(rn) <- "htest"
  }

  return(list(positive=rp, negative=rn))
}

siDominance <- function(HH, HL, LH, LL, method="ks") {
    #if (method == "ks") {
      
    hh.hl <- suppressWarnings( ks.test(HH,HL,alternative="greater",exact=FALSE) )
    hl.hh <- suppressWarnings( ks.test(HH,HL,alternative="less",exact=FALSE) )
    hh.lh <- suppressWarnings( ks.test(HH,LH,alternative="greater",exact=FALSE) )
    lh.hh <- suppressWarnings( ks.test(HH,LH,alternative="less",exact=FALSE) )
    hl.ll <- suppressWarnings( ks.test(HL,LL,alternative="greater",exact=FALSE) )
    ll.hl <- suppressWarnings( ks.test(HL,LL,alternative="less",exact=FALSE) )
    lh.ll <- suppressWarnings( ks.test(LH,LL,alternative="greater",exact=FALSE) )
    ll.lh <- suppressWarnings( ks.test(LH,LL,alternative="less",exact=FALSE) )
    

    dominance <- as.data.frame(
      list(Test=c("S.hh > S.hl", "S.hh > S.lh", "S.hl > S.ll", "S.lh > S.ll", 
             "S.hh < S.hl", "S.hh < S.lh", "S.hl < S.ll", "S.lh < S.ll"),
      statistic=c(hh.hl$statistic, hh.lh$statistic, hl.ll$statistic, lh.ll$statistic,
                  hl.hh$statistic, lh.hh$statistic, ll.hl$statistic, ll.lh$statistic),
      p.value=c(hh.hl$p.value, hh.lh$p.value, hl.ll$p.value, lh.ll$p.value,
                hl.hh$p.value, lh.hh$p.value, ll.hl$p.value, ll.lh$p.value)) )
      

    #} else if (method=="dp") {
    #  dominance <- c( DPdom(HH,HL)$test, DPdom(HH,LH)$test, 
    #                  DPdom(HL,LL)$test, DPdom(LH,LL)$test)
    #}
    return(dominance)
}


mic.test <- function(HH, HL, LH, LL, method=c("art","anova")) {
  method <- match.arg(method, c("art", "anova"))
  DNAME <- paste("\nHH:", deparse(substitute(HH)), "\tHL:", deparse(substitute(HL)), 
                 "\nLH:", deparse(substitute(LH)), "\tLL:", deparse(substitute(LL)) )

  HH <- HH[!is.na(HH)]
  HL <- HL[!is.na(HL)]
  LH <- LH[!is.na(LH)]
  LL <- LL[!is.na(LL)]

  STATISTIC <- (mean(LL) - mean(LH))  - (mean(HL) - mean(HH))
  names(STATISTIC) <- "MIC"

  allrt <- c(HH, HL, LH, LL)

  if (method=="art") {
    METHOD <- "Adjusted Rank Transform test of the MIC"
    n1 <- length(HH)
    n2 <- length(HL)
    n3 <- length(LH)
    n4 <- length(LL)

    h1 <- c(rep(1, n1+n2), rep(0,n3+n4))
    h2 <- c(rep(1, n1), rep(0, n2), rep(1,n3), rep(0, n4))

    mA0 <- sum( allrt * (1-h1) ) / sum(1-h1)
    mA1 <- sum( allrt * h1 ) / sum(h1)
    mB0 <- sum( allrt * (1-h2) ) / sum(1-h2)
    mB1 <- sum( allrt * h2 ) / sum(h2)

    allrt.m <- allrt - (1-h1)*mA0- h1*mA1 - (1-h2)*mB0 - h2*mB1 
    allrt.m <- round(allrt.m, 15)
    ranks <- rank(allrt.m, ties.method="average")

    stat <- anova(lm(ranks ~ h1*h2))
    PVAL <- stat$"Pr(>F)"[3]
    #return (list(MIC=mic, Df=c(1,length(allrt)-4), "F value"=rval$"F value"[3], p.val=rval$"Pr(>F)"[3]))
  }

  else if(method=="anova"){
    METHOD <- "ANOVA test of the MIC"
    op <- options(contrasts=c("contr.helmert", "contr.poly"))

    h1vec <- c(rep(1, length(HH)), rep(1, length(HL)), 
               rep(0, length(LH)), rep(0, length(LL)))
    h2vec <- c(rep(1, length(HH)), rep(0, length(HL)), 
               rep(1, length(LH)), rep(0, length(LL)))
    stat <- anova(lm(allrt ~h1vec * h2vec))
    options(op)
    PVAL <- stat$"Pr(>F)"[3]
    #return(list(MIC=mic, Df=c(1,length(allrt)-4), "F value"=rval$"F value"[3], p.val=rval$"Pr(>F)"[3]))
  }

  RVAL <- list(statistic=STATISTIC, p.value=PVAL, alternative="the MIC is not zero",
      method=METHOD, data.name=DNAME)
  class(RVAL) <- "htest"
  return(RVAL)

}
