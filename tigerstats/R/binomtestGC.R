#' @title Exact Procedures for a Single Proportion

#' @description Wrapper for binom.test in package \code{stats}.  Employs the binomial distribution 
#' in inferential procedures for a single proportion.
#' 
#' @rdname binomtestGC
#' @usage binomtestGC(x,data=parent.frame(),n=numeric(),p=NULL,
#'                          alternative=c("two.sided","less","greater"),
#'                          success="yes",conf.level=0.95,graph=FALSE,verbose=TRUE)
#' @param x Either a formula or a numeric vector.  If formula, it must be of the form ~x
#' indicating the single variable under study.  When summary data are provided, x is a numeric vector of 
#' success counts.
#' @param data Data frame that supplies the variable x. If not found in data, the variable is searched
#' for in the parent environment.
#' @param n When not empty, this is a numeric vector giving the size of the sample.
#' @param p Specifies Null Hypothesis value for population proportion.  If not set, no test is performed.
#' @param alternative "two.sided" requests computation of a two-sided P-value;  
#' other possible values are "less" and "greater".
#' @param success  When x is a formula, this argument indicates which value of variable x is being counted as a success.  
#' When working with formula-data input the value of this parameter MUST be set, even when the variable has only
#' two values.
#' @param conf.level Number between 0 and 1 indicating the confidence-level of the interval supplied.
#' @param graph If TRUE, plot graph of P-value.  Ignored if no test is performed.
#' @param verbose Determines whether to return lots of information or only the basics
#' @return an object of class GCbinomtest.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #Confidence interval only:
#' binomtestGC(~sex,data=m111survey,success="female")
#' 
#' #Confidence interval and two-sided test, Null Hypothesis p = 0.5:
#' binomtestGC(~sex,data=m111survey,success="female",p=0.5)
#' 
#' #For confidence level other than 95%, use conf.level argument.
#' #For 90% interval:
#' binomtestGC(~sex,data=m111survey,success="female",conf.level=0.90)
#' 
#' #For one-sided test, set alternative argument as desired:
#' binomtestGC(~sex,data=m111survey,p=0.50,
#'     success="female",alternative="greater")
#' 
#' #Summary data:
#' #In one sample, 40 successes in 100 trials.  Testing whether p = 0.45.
#' binomtestGC(x=40,n=100,p=0.45)
binomtestGC <-
  function(x,data=parent.frame(),n=numeric(),
           p=NULL,
           alternative=c("two.sided","less","greater"),
           success="yes",
           conf.level=0.95,
           graph=FALSE,verbose=TRUE)  { 
    
    alternative <- match.arg(alternative)
    
    if (is(x,"formula"))  {
      
      prsd <- ParseFormula(x)
      pullout <- as.character(prsd$rhs)
      
      if (length(pullout) > 1) stop("Incorrect formula")
      
      varname <- pullout[1]
      variable <- simpleFind(varName=varname,data=data)
      TallyTable <- xtabs(~variable)
      
      if (!(success %in% unique(variable))) {
        stop("No sucesses found.  Did you fail to specify success correctly?")
      }
      successes <- TallyTable[success]
      names(successes) <- NULL  #keep name off of p-value in the returned results
      trials <- sum(TallyTable)
      names(trials) <- NULL
      
    }
      
    
    if (!is(x,"formula")) {  #we have summary data
      
        successes <- x
        trials <- n
        
    }
    
    if (is.null(p))  { #no test
      
      res <- stats::binom.test(successes,trials,p=0.5,conf.level=conf.level)
    } else res <- stats::binom.test(successes,trials,p=p,alternative=alternative,conf.level=conf.level)
    
    
    results <- structure(list(
        successes=successes,
        trials=trials,
        p.value=res$p.value,
        conf.int=res$conf.int,
        conf.level=conf.level,
        p=p,
        alternative=alternative,
        input=x,
        varname = ifelse(is(x,"formula"),varname,""),
        success=success,
        verbose=verbose,
        graph=graph),
        class="GCbinomtest")
    
    return(results)
   

  }#end binomtestGC

#' @title Print Function for binomestGC

#' @description Utility print function
#' @keywords internal
#' 
#' @rdname print.GCbinomtest
#' @method print GCbinomtest
#' @usage 
#' \S3method{print}{GCbinomtest}(x,...)
#' @param x An object of class GCbinomtest.
#' @param \ldots ignored
#' @return Output to the console.
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @export
print.GCbinomtest <- function(x,...) {
  
  alternative <- x$alternative
  graph <- x$graph
  p <- x$p
  conf.level <- x$conf.level
  verbose <- x$verbose
  input <- x$input
  varname <- x$varname
  successes <- x$successes
  trials <- x$trials
  n <- trials #in case I fall into the habit of calling trials n
  conf.int <- x$conf.int
  p.value <- x$p.value
  
  cat("Exact Binomial Procedures for a Single Proportion p:\n")
  if (is(input,"formula")) {
    cat("\tVariable under study is",varname,"\n")
  } else cat("\tResults based on Summary Data\n")
  if (verbose==TRUE) {   
    cat("\n")
    cat("Descriptive Results: ",successes,"successes in",trials,"trials\n\n")
  }
  
  p.hat <- successes/trials
  se.phat <- sqrt(p.hat*(1-p.hat)/trials)
  
  if (verbose==TRUE) {
    cat("Inferential Results:\n\n")
    cat("Estimate of p:\t",round(p.hat,4),"\n")
    cat("SE(p.hat):\t",round(se.phat,4),"\n\n")
  }
  cat(conf.level*100,"% Confidence Interval for p:\n\n",sep="")
  int <- conf.int
  cat(sprintf("%-10s%-20s%-20s","","lower.bound","upper.bound"),"\n")
  cat(sprintf("%-10s%-20f%-20f","",int[1],int[2]),"\n\n")
  
  if (!is.null(p)) {
    if (verbose==TRUE) {
      cat("Test of Significance:\n\n")
      symbol <- switch(alternative,
                       less="<",
                       greater=">",
                       two.sided="!=")
      cat("\tH_0:  p =",p,"\n")
      cat("\tH_a:  p",symbol,p,"\n\n")
    }
    cat("\tP-value:\t\tP =",round(p.value,4),"\n")
    
    #deal with graphing in two-sided tests:
    twoSide <- function(count,n,p) {
      x <- 0:n
      p.x <- dbinom(x,size=n,prob=p)
      ourProb <- dbinom(count,size=n,prob=p)
      
      smallOnes <- x[p.x <= ourProb]
      
      if (length(smallOnes) == (n+1)) {
        invisible(pbinomGC(c(0,n),size=n,prob=p,region="between",graph=T))
      } else {
        bigger <- x[p.x > ourProb]
        if (length(bigger) ==1) {
          invisible(pbinomGC(c(bigger,bigger),size=n,prob=p,region="outside",graph=T))
        }
        
        if (length(bigger) > 1) {
          lower <- bigger[1]
          upper <- bigger[length(bigger)]
          invisible(pbinomGC(c(lower,upper),size=n,prob=p,region="outside",graph=T))
        }
      } # end else
    }  #end twoSide
    
    if (graph)   {  
      switch(alternative,
             less=invisible(pbinomGC(successes,size=trials,prob=p,region="below",graph=T)),
             greater=invisible(pbinomGC(successes-1,size=trials,prob=p,region="above",graph=T)),
             two.sided=twoSide(count=successes,n=trials,p=p)
      )
    }
    
  }
}
