#' @title Proportions Procedures

#' @description Employs the normal approximation to perform test for one or two proportions.
#' 
#' @rdname proptestGC
#' @usage proptestGC(x,data=parent.frame(),n=numeric(),p=NULL,
#'                          alternative=c("two.sided","less","greater"),
#'                          success="yes",first=NULL,conf.level=0.95,
#'                          correct=TRUE,graph=FALSE,verbose=TRUE)
#' @param x Either a formula or a numeric vector.  If formula, it must be of the form ~x
#' indicating the single variable under study, or of the form ~x+y, in which case x is the explanatory grouping variable
#' (categorical with two values) and y is the response categorical variable with two values.
#' When summary data are provided, x is a numeric vector of success counts.
#' @param n When not empty, this is a numeric vector giving the size of each sample.
#' @param p Specifies Null Hypothesis value for population proportion.  If not set, no test is performed.
#' @param data Data frame that supplies the variables x and y.  If any are not in data, then they will be
#' searched for in the parent environment.
#' @param alternative "two.sided" requests computation of a two-sided P-value;  other possible values are "less" and "greater".
#' @param success  When x is a formula, this argument indicates which value of variable x (in case of ~x) or y (in case of ~x+y)
#' is being counted as a success.  When working
#' with formula-data input the value of this parameter MUST be set, even when the variable has only
#' two values.
#' @param first When performing 2-sample procedures, this argument specifies which value of
#' the explanatory variable constitutes the first group.
#' @param conf.level Number between 0 and 1 indicating the confidence-level of the interval supplied.
#' @param correct Applies continuity correction for one-proportion procedures.  It is ignored when
#' when 2-proportions are performed.
#' @param graph If TRUE, plot graph of P-value.
#' @param verbose Indicates how much output goes to the console
#' @return A list, either of class "gcp1test" (one-proportion) or "gcp2test" (two proportions).  
#' Components of this list that may be usefully queried include:  "statistic", "p.value", and "interval".
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' data(m111survey)
#' #2-proportions, formula-data input, 95%-confidence interval only:
#' proptestGC(~sex+seat,data=m111survey,success="2_middle")
#' 
#' #For other confidence levels, use argument conf.level.  For 90%-interval for one proportion p:
#' proptestGC(~sex,data=m111survey,success="male",conf.level=0.90)
#' 
#' #one proportion, formula-data input, confidence interval and two-sided test with H_0:  p = 0.33:
#' proptestGC(~seat,data=m111survey,success="1_front",p=0.33)
#' 
#' #Summary data:
#' #In first sample, 23 successes out of 100 trials.  In second sample, 33 out of 110.
#' proptestGC(x=c(23,33),n=c(100,110))
#' 
#' #Summary data:
#' #In one sample, 40 successes in 100 trials.  Testing whether p = 0.45.
#' proptestGC(x=40,n=100,p=0.45,correct=TRUE)
#' 
#' #Want less output?  Set argument verbose to FALSE:
#' proptestGC(~sex+seat,data=m111survey,success="2_middle",p=0.33,verbose=FALSE)
proptestGC <-
  function(x,data=parent.frame(),n=numeric(),
           p=NULL,
           alternative=c("two.sided","less","greater"),
           success="yes",first=NULL,
           conf.level=0.95,
           correct=TRUE,graph=FALSE,verbose=TRUE)  {
    
    alternative <- match.arg(alternative)
    statistic <- FALSE
    p.value <- FALSE #these will get numerical values if a test is to be performed
    
    #Small Utility Function
    GetP <- function(stat,alternative) {
      switch(alternative,
             less=pnorm(stat),
             greater=pnorm(stat,lower.tail=FALSE),
             two.sided=2*pnorm(abs(stat),lower.tail=FALSE))
    }
    
    #small utility function:conf int for one prop
    GetCI1 <- function(est,se,conf.level,alternative) {
      switch(alternative,
             less=c(lower=0,upper=min(1,est+qnorm(conf.level)*se)),
             two.sided=c(lower=max(0,est+qnorm((1-conf.level)/2)*se),
                         upper=min(1,est+qnorm((1-conf.level)/2,lower.tail=FALSE)*se)),
             greater=c(max(0,lower=est+qnorm(1-conf.level)*se),upper=1)
      )   
    }
    
    #utility function:  conf int for 2 diff of 2 props
    GetCI2 <- function(est,se,conf.level,alternative) {
      switch(alternative,
             less=c(lower=-1,upper=min(1,est+qnorm(conf.level)*se)),
             two.sided=c(lower=max(-1,est+qnorm((1-conf.level)/2)*se),
                         upper=min(1,est+qnorm((1-conf.level)/2,lower.tail=FALSE)*se)),
             greater=c(max(-1,lower=est+qnorm(1-conf.level)*se),upper=1)
      )   
    }
    
    #This function handles formula-data input:
    proptestGC.form <- function(form,data,success="yes",
                                 alternative="two.sided",
                                 conf.level=0.95,
                                 correct=TRUE )  {
      

      
      prsd <- ParseFormula(form)
      
      pullout <- as.character(prsd$rhs)
      
      if (length(pullout)==3) {#we have a bona fide formula
        expname <- as.character(prsd$rhs)[2]
        respname <- as.character(prsd$rhs)[3]
        
        explanatory <- simpleFind(varName=expname,data=data)
        
        expEntries <- unique(explanatory)
        nonTrivial <- length(expEntries[!is.na(expEntries)])
        
        #Check to see that explanatory variable has exactly two values:
        if (nonTrivial != 2) stop(paste(expname,"must have exactly two values."))
        
        
        response <- simpleFind(varName=respname,data=data)
        
        if (!(success %in% unique(response))) {
          stop("No sucesses found.  Did you fail to specify success correctly?")
        }
        
        TwoWayTable <- xtabs(~explanatory+response)
        
        success.tab <- TwoWayTable[,success]
        
        
        #Count all values that are not success as failure:
        failure.names <- colnames(TwoWayTable)[colnames(TwoWayTable)!=success]
        failure.tab <- as.matrix(TwoWayTable[,failure.names])
        comb.failure.tab <- rowSums(failure.tab)
        
        twobytwotab <- cbind(success.tab,comb.failure.tab)
        rownames(twobytwotab) <- NULL
        colnames(twobytwotab) <- NULL
        
        flip <- FALSE
        if (!is.null(first)) {
          if (!(first %in% unique(explanatory))) {
            stop(paste(first,"is not a value of the variable",expname))}
          if (sort(unique(explanatory))[1]!=first) flip <- TRUE
        }
        
        tttp <- prop.table(twobytwotab,margin=1)
        
        if (flip) {
          tttp <- tttp[c(2,1),]
          twobytwotab <- twobytwotab[c(2,1),]
        }
        
        p1.hat <- tttp[1,1]
        p2.hat <- tttp[2,1]
        n1 <- rowSums(twobytwotab)[1]
        n2 <- rowSums(twobytwotab)[2]
        Estimator <- p1.hat-p2.hat
        SE <- sqrt(p1.hat*(1-p1.hat)/n1+p2.hat*(1-p2.hat)/n2)
        
        if (!is.null(p)) {
        statistic <- Estimator/SE
        p.value <- GetP(statistic,alternative)
        }
        
        interval <- GetCI2(Estimator,SE,conf.level,alternative)
        
        successes <- twobytwotab[,1]
        sampsizes <- c(n1,n2)
        SummTab <- cbind(successes,
                         sampsizes,
                         props=successes/sampsizes)
        colnames(SummTab) <- c(success,"n","estimated.prop")
        rownames(SummTab) <- rownames(TwoWayTable)
        if (flip) rownames(SummTab) <- rev(rownames(TwoWayTable))
        
        results <- list(SummTab=SummTab,statistic=statistic,p.value=p.value,interval=interval,graph=graph,
                        explanatory=expname,response=respname,se=SE,conf.level=conf.level,
                        alternative=alternative,verbose=verbose)
        class(results)  <- "gcp2test"
        return(results)
      }
      
      if(length(pullout)==1)  {
        varname <- pullout[1]
        variable <- simpleFind(varName=varname,data=data)
        TallyTable <- xtabs(~variable)
        
        if (!(success %in% unique(variable))) {
          stop("No sucesses found.  Did you fail to specify success correctly?")
        }
        successes <- TallyTable[success]
        
        
        trials <- sum(TallyTable)
        
        Estimator <- successes/trials
        SE <- sqrt(Estimator*(1-Estimator)/trials)
        
        if(!is.null(p)) {
          
        if (correct==TRUE) {
          if (alternative=="less") statistic <- (Estimator+0.5/trials-p)/SE
          if (alternative=="greater") statistic <- (Estimator-0.5/trials-p)/SE
          if (alternative=="two.sided") statistic <- (Estimator-sign(Estimator)*0.5/trials-p)/SE
        } else statistic <- (Estimator-p)/SE
       
        p.value <- GetP(statistic,alternative)
        
          }
        
        interval <- GetCI1(Estimator,SE,conf.level,alternative)
        
        SummTab <- cbind(successes,
                         trials,
                         Estimator)
        colnames(SummTab) <- c(success,"n","estimated.prop")
        rownames(SummTab) <- ""
        
        results <- list(SummTab=SummTab,statistic=statistic,p.value=p.value,interval=interval,graph=graph,
                        variable=varname,se=SE,conf.level=conf.level,
                        alternative=alternative,correct=correct,p=p,verbose=verbose)
        class(results)  <- "gcp1test"
        return(results)
      }
    } #end form function
    
    
    if (is(x,"formula"))  {
      return(proptestGC.form(x,data,alternative=alternative,
                              conf.level=conf.level,
                              success=success,
                              correct=correct))
    }
    
    if (!is(x,"formula")) {  #we have summary data
      if (length(n)!=length(x))  {
        stop("Vector of counts and vector of trials must have same length")
      } 
      
      if (length(x)==1)  {  # one-sample
        successes <- x
        trials <- n
        
        Estimator <- successes/trials
        SE <- sqrt(Estimator*(1-Estimator)/trials)
        
        if(!is.null(p)) {
          
        if (correct==TRUE) {
          if (alternative=="less") statistic <- (Estimator+0.5/trials-p)/SE
          if (alternative=="greater") statistic <- (Estimator-0.5/trials-p)/SE
          if (alternative=="two.sided") statistic <- (Estimator-sign(Estimator)*0.5/trials-p)/SE
        } else statistic <- (Estimator-p)/SE
        
        p.value <- GetP(statistic,alternative)
        
        }
        
        interval <- GetCI1(Estimator,SE,conf.level,alternative)
        
        SummTab <- cbind(successes,
                         trials,
                         Estimator)
        colnames(SummTab) <- c("successes","n","estimated.prop")
        rownames(SummTab) <- ""
        
        results <- list(SummTab=SummTab,statistic=statistic,p.value=p.value,interval=interval,graph=graph,
                        variable=NA,se=SE,conf.level=conf.level,
                        alternative=alternative,correct=correct,p=p,verbose=verbose)
        class(results)  <- "gcp1test"
        return(results)
      }
      
      if (length(x)==2) {
        n1 <- n[1]
        n2 <- n[2]
        p1.hat <- x[1]/n1
        p2.hat <- x[2]/n2
        
        Estimator <- p1.hat-p2.hat
        SE <- sqrt(p1.hat*(1-p1.hat)/n1+p2.hat*(1-p2.hat)/n2)
        
        if(!is.null(p)) {
        statistic <- Estimator/SE  
        p.value <- GetP(statistic,alternative)
        }
        
        
        interval <- GetCI2(Estimator,SE,conf.level,alternative)
        
        SummTab <- cbind(x,
                         n,
                         x/n)
        colnames(SummTab) <- c("successes","n","estimated.prop")
        rownames(SummTab) <- c("Group 1","Group 2")
        
        results <- list(SummTab=SummTab,statistic=statistic,p.value=p.value,interval=interval,graph=graph,
                        explanatory=NA,response=NA,se=SE,conf.level=conf.level,
                        alternative=alternative,verbose=verbose)
        class(results)  <- "gcp2test"
        return(results)
      }
      
    }#end summary data treatment
    
  }#end proptestGC


#' @title Print Function for GC Proportion Test (One-Sample)

#' @description Utility print function
#' @keywords internal
#' 
#' @rdname print.gcp1test
#' @method print gcp1test
#' @usage 
#' \S3method{print}{gcp1test}(x,...)
#' @param x An object of class gcp1test.
#' @param \ldots ignored
#' @return Output to the console.
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @export 
print.gcp1test <- function(x,...)  {
  gcp1test <- x
  verbose <- gcp1test$verbose
  
  odigits <- getOption("digits")
  options(digits=4)
  
  cat("\n\nInferential Procedures for a Single Proportion p:\n")
  if (!is.na(gcp1test$variable)) {
    cat("\tVariable under study is",gcp1test$variable,"\n")
  } else cat("\tResults based on Summary Data\n")
  if (gcp1test$correct==TRUE) {
    cat("\tContinuity Correction Applied to Test Statistic\n")
  }
  
  tab <- gcp1test$SummTab
  
  if (verbose) {
    cat("\n\n")
    cat("Descriptive Results:\n\n")
    
    print(tab)
    
    cat("\n")
    cat("\n")
  }
  checker <- min(tab[1,1],tab[1,2]-tab[1,1])
  if (checker < 10) {
    cat("WARNING:  Either the number of successes or \nthe number of failures is below 10.\nThe normal approximation for confidence intervals\nand P-value may be unreliable\n\n\n")
  }
  
  
  if (verbose) {  
    cat("Inferential Results:\n\n")
    cat("Estimate of p:\t",tab[1,3],"\n")
    cat("SE(p.hat):\t",gcp1test$se,"\n\n")
  }
  cat(gcp1test$conf.level*100,"% Confidence Interval for p:\n\n",sep="")
  int <- gcp1test$interval
  cat(sprintf("%-10s%-20s%-20s","","lower.bound","upper.bound"),"\n")
  cat(sprintf("%-10s%-20f%-20f","",int[1],int[2]),"\n\n")
  
  if (gcp1test$p.value) {
    if (verbose) {
      cat("Test of Significance:\n\n")
      symbol <- switch(gcp1test$alternative,
                       less="<",
                       greater=">",
                       two.sided="!=")
      cat("\tH_0:  p =",gcp1test$p,"\n")
      cat("\tH_a:  p",symbol,gcp1test$p,"\n\n")
    }
    cat("\tTest Statistic:\t\tz =",gcp1test$statistic,"\n")
    cat("\tP-value:\t\tP =",gcp1test$p.value,"\n")
    
    Grapher <- function(stat,alt) {
      rstat <- round(stat,2)
      switch(alt,
             less=invisible(pnormGC(rstat,region="below",graph=T)),
             greater=invisible(pnormGC(rstat,region="above",graph=T)),
             two.sided=invisible(pnormGC(c(-abs(rstat),abs(rstat)),region="outside",graph=T))
      )
    }
    
    if (gcp1test$graph) {
      Grapher(stat=gcp1test$statistic,alt=gcp1test$alternative)
    }
    
  }
  
  options(digits=odigits)
}


#' @title Print Function for GC Proportions Test (Two-Sample)

#' @description Utility print function
#' @keywords internal
#' 
#' @rdname print.gcp2test
#' @method print gcp2test
#' @usage 
#' \S3method{print}{gcp2test}(x,...)
#' @param x An object of class gcp2test.
#' @param \ldots ignored
#' @return Output to the console.
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @export
print.gcp2test <- function(x,...)  {
  gcp2test <- x
  verbose <- gcp2test$verbose
  
  odigits <- getOption("digits")
  options(digits=4)
  
  cat("\n\nInferential Procedures for the Difference of Two Proportions p1-p2:\n")
  if (!is.na(gcp2test$explanatory)) {
    cat("\t",gcp2test$response,"grouped by",gcp2test$explanatory,"\n")
  } else cat("\tResults taken from summary data.\n")
  
  tab <- gcp2test$SummTab
  
  if (verbose) {
    cat("\n\n")
    cat("Descriptive Results:\n\n")
    
    print(tab)
    
    cat("\n")
    cat("\n")
  }
  
  checker <- min(tab[1,1],tab[2,1],tab[1,2]-tab[1,1],tab[2,2]-tab[2,1])
  if (checker < 10) {
    cat("WARNING:  In at least one of the two groups\n",
        "number of successes or number of failures is below 10.\n",
        "The normal approximation for confidence intervals\n",
        "and P-value may be unreliable.\n\n",sep="")
  }
  
  if (verbose) {
    cat("Inferential Results:\n\n")
    cat("Estimate of p1-p2:\t",tab[1,3]-tab[2,3],"\n")
    cat("SE(p1.hat - p2.hat):\t",gcp2test$se,"\n\n")
  }
  cat(gcp2test$conf.level*100,"% Confidence Interval for p1-p2:\n\n",sep="")
  int <- gcp2test$interval
  cat(sprintf("%-10s%-20s%-20s","","lower.bound","upper.bound"),"\n")
  cat(sprintf("%-10s%-20f%-20f","",int[1],int[2]),"\n\n")
  
  if(gcp2test$p.value) {
    if (verbose) { 
      cat("Test of Significance:\n\n")
      symbol <- switch(gcp2test$alternative,
                       less="<",
                       greater=">",
                       two.sided="!=")
      cat("\tH_0:  p1-p2 = 0\n")
      cat("\tH_a:  p1-p2",symbol,"0\n\n")
    }
    cat("\tTest Statistic:\t\tz =",gcp2test$statistic,"\n")
    cat("\tP-value:\t\tP =",gcp2test$p.value,"\n")
    
    Grapher <- function(stat,alt) {
      rstat <- round(stat,2)
      switch(alt,
             less=invisible(pnormGC(rstat,region="below",graph=T)),
             greater=invisible(pnormGC(rstat,region="above",graph=T)),
             two.sided=invisible(pnormGC(c(-abs(rstat),abs(rstat)),region="outside",graph=T))
      )
    }
    
    if (gcp2test$graph) {
      Grapher(stat=gcp2test$statistic,alt=gcp2test$alternative)
    }
    
  }
  options(digits=odigits)
}
