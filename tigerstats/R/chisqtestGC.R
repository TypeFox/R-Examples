#' @title Chi-Square Test (GC version)

#' @description Perform chi-square test, either goodness of fit or test for association.  Enter either
#' formula-data input or a summary table.  Simulation is optional.
#' 
#' @rdname chisqtestGC
#' @usage chisqtestGC(x, data = parent.frame(), p = NULL, graph = FALSE,
#'              simulate.p.value = FALSE, B = 2000, verbose = TRUE)
#' @param x Could be a formula.  If so, either ~var (for goodness of fit) or ~var1+var2 (for test for association).
#' Otherwise either a table, matrix or vector of summary data.
#' @param data dataframe supplying variables for formula x.  If variables in x ar enot found in the data,
#' then they will be searched for in the parent environment.
#' @param p For goodness of fit, a vector of probabilities.  This will be automatically scaled so as to sum
#' to 1.  Negative elements result in an error message.
#' @param graph produce relevant graph for P-value (chi-square curve or histogram of simulation results).
#' @param simulate.p.value If FALSE, use a chi-square distribution to estimate the P-value.  Other possible
#' values are "random" and "fixed" and TRUE.  Random effects are suitable for resampling when the data are a random
#' sample from a poulation.  Fixed effects assume that the values of the explanatory variable (row variable for table,
#' var1 in formula ~var1+var2) remain fixed in resampling, and values of response variable are random with null
#' distribution estimated from the data.  When set to TRUE, we implement an equivalent to R's routine.  
#' In our view procedure is
#' most suitable when the data come from a randomized experiment in which the treatment groups are
#' the values of the explanatory variable.
#' @param B number of resamples to take.
#' @param verbose If TRUE, include lots of information in the output.
#' @return an object of class GCchisqtest
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #Goodness of fit test for one factor variable:
#' chisqtestGC(~seat,data=m111survey,p=c(1/3,1/3,1/3))
#' 
#' #Test for relationship between two factor variables:
#' chisqtestGC(~sex+seat,data=m111survey)
#' 
#' #You can input a two-way table directly into chisqtestGC():
#' SexSeat <- xtabs(~sex+seat,data=m111survey)
#' chisqtestGC(SexSeat)
#' 
#' #Several types of simulation are possible, e.g.:
#' chisqtestGC(~weather+crowd.behavior,data=ledgejump,simulate.p.value="fixed",B=2500)
#' 
#' #For less ouptut, set argument verbose to FALSE:
#' chisqtestGC(~sex+seat,data=m111survey,verbose=FALSE)
chisqtestGC <- 
  function (x,data=parent.frame(),p=NULL,graph=FALSE,simulate.p.value=FALSE,B=2000,verbose=TRUE) 
  {
    
###########################################################
#  begin with utiltiy functions
##########################################################
    
    exp.counts <- function(x) (rowSums(x) %*% t(colSums(x)))/sum(x)
    
    chisq.calc <- function(x) {
      expected <- exp.counts(x)
      contributions <- (x - expected)^2/expected
      return(sum(contributions[!is.nan(contributions)]))
    }
    
##########################################################################  
#    simulation for simulate.p.value==TRUE, when test is for association
##########################################################################

    DoubleFixedResampler <- function(x,n) {
      expected <- exp.counts(x)
      csq <- function(x) {
        sum((x-expected)^2/expected)
      }
      statistic <- csq(x)
      nullDist <- numeric(n)
      
      r <- rowSums(x)
      c <- colSums(x)
      
      countOver <- 0
      simsSoFar <- 0
      ourLimit <- 10000 # reps to handle at once
      
      while(simsSoFar <n) {
        reps <- min(n-simsSoFar,ourLimit)
        rtabs <- r2dtable(reps,r=r,c=c)
        sims <- sapply(rtabs,FUN=csq,USE.NAMES=FALSE)
        nullDist[(simsSoFar+1):(simsSoFar+reps)] <- sims
        simsSoFar <- simsSoFar + reps
        countOver <- countOver + length(sims[sims >= statistic])
      }
      return(list(nullDist=nullDist,countOver=countOver))
    }


#################################################   
#   simulation for "fixed" and "random"
#################################################

    RandFixedResampler <- function (x, n, effects = "random") 
    {
      #x is a two-way table, n is number of resamples
      TableResampler <- function(x, n = 1000, effects) {
        rowsampler <- function(x, p) {
          rmultinom(1, size = sum(x), prob = p)
        }
        
      table.samp <- function(x) {
        nullprobs <- colSums(x)/sum(x)
        resamp <- t(apply(x, 1, rowsampler, p = nullprobs))
        rownames(resamp) <- rownames(x)
        colnames(resamp) <- colnames(x)
        as.table(resamp)
      }
      
      rtabsamp <- function(x, n) {
          expected <- exp.counts(x)
          probs <- expected/sum(x)
          resamp.tab <- rmultinom(1, size = n, prob = probs)
          resamp.tab <- matrix(resamp.tab, nrow = nrow(x))
          rownames(resamp.tab) <- rownames(x)
          colnames(resamp.tab) <- colnames(x)
          return(resamp.tab)
      }
      
      resampled.tabs <- array(0, dim = c(nrow(x), ncol(x),n))
      
        if (effects == "fixed") {
          for (i in 1:n) {
            resampled.tabs[, , i] <- table.samp(x)
          }
          return(resampled.tabs)
        }
        if (effects == "random") {
          for (i in 1:n) {
            resampled.tabs[, , i] <- rtabsamp(x, sum(x))
          }
          return(resampled.tabs)
        }
      }
     
      nullDist <- apply(TableResampler(x, n, effects = effects), 
                        3, chisq.calc)
      
      return(nullDist)
    }#end of ChisqResampler
    

###################################################################
#   Simulation in goodness of fit
###################################################################

    GoodnessResampler <- function(x,n,p) {
      
      rowsampler <- function(x, p) {
        rmultinom(1, size = sum(x), prob = p)
      }
      
      chisq.calc.1 <- function(x,p) {
        expected <- sum(x)*p
        return(sum((x - expected)^2/expected))
      }
      
      nullDist <- numeric(n)
      
      for (i in 1:n) {
        resamp.tab <- rowsampler(x,p)
        nullDist[i] <- chisq.calc.1(resamp.tab,p)
      }
      
      return(nullDist)
    }#end of GoodnessResampler
    

 
 ##################################################################
 ####
 ####
 #### Process Input
 ####
 ####
 #################################################################
    
    type <- NULL #will be set to the type of test (association or goodness)
    
    #first see if we have formula-data input, or summary data
    if (is(x,"formula")) #we have formula-data input
    {
    
    prsd <- ParseFormula(x)
    pullout <- as.character(prsd$rhs)
    
    if (length(pullout)==3) #Test for association
      
      {
      type <- "association"
      
      expname <- as.character(prsd$rhs)[2]
      respname <- as.character(prsd$rhs)[3]
      
      explanatory <- simpleFind(varName=expname,data=data)
      response <- simpleFind(varName=respname,data=data)
      data2 <- data.frame(explanatory,response)
      names(data2) <- c(expname,respname)
      tab <- table(data2)
      res <- suppressWarnings(chisq.test(tab))
      
    } #end processing of association test
    
    if(length(pullout)==1)  #goodness of fit
      {
      type <- "goodness"
      varname <- pullout[1]
      

      variable <- simpleFind(varName=varname,data=data)
      tab <- table(variable)
      x <- tab  #provides proper input later if simulation is desired
      res <- suppressWarnings(chisq.test(tab,p=p,rescale.p=TRUE))
      #Note:  if variable has inherited levels (from a previous life) that user no longer
      #expects to see, then length of table will exceed length of p and there will be
      #problems.
      
    } #end processing of goodness of fit test
    
    }#end formula processing
    
  if (!is(x,"formula"))  #we have summary data
  {
    if (length(dim(x))>2) #array more than two dimensions
    {
      stop("For tables with more than two dimensions, use chisq.test()")
    }
    
    tab <- as.table(x)
    if (length(dim(tab))==1) {
      type <- "goodness"
      res <- suppressWarnings(chisq.test(x,p=p,rescale.p=TRUE))
    }#end of goodness of fit processing
    
    if (length(dim(tab))==2) {
      type <- "association"
      res <- suppressWarnings(chisq.test(tab))
    }#end of association processing
    
  }#end processing for summary data


###################################################################################
# preliminary collection of results, these may be altered if simulaion is desired
###################################################################################

  statistic <- res$statistic
  p.value <- res$p.value
  parameter <- res$parameter
  method <- res$method
  observed <- tab
  expected <- res$expected
  residuals <- res$residuals
  if (simulate.p.value==F) {
    sims <- ""
  }
    
################################################
# Handle Simulation, as needed
################################################
  
    if (simulate.p.value %in% c("random","fixed") && type=="association") {
      statistic <- sum((observed-expected)^2/expected) #don't want Yates
      nullDist <- RandFixedResampler(observed,n=B,effects=simulate.p.value)
      p.value <- (length(nullDist[nullDist >= statistic])+1)/(B+1)
      sims <- nullDist
    }
    
    if(simulate.p.value==TRUE  && type=="association") { 
      #user requested R's routines, so give something very close to it
      statistic <- sum((res$observed-res$expected)^2/res$expected) #don't want Yates
      simResults <- DoubleFixedResampler(x=tab,n=B)
      countOver <- simResults$countOver
      p.value <- (countOver+1)/B
      sims <- simResults$nullDist
    }
  
  if (simulate.p.value %in% c("random","fixed") && type=="goodness") {
    stop("For simulation in a goodness of fit test, set simulate.p.value to TRUE")
  }
  
  if(simulate.p.value==TRUE  && type=="goodness") { #we pick up the simulated values ourselves
    expected <- sum(tab)*p
    statistic <- sum((tab-expected)^2/expected) #don't want Yates
    sims <- GoodnessResampler(x,n=B,p=p)
    p.value <- (length(sims[sims >= statistic])+1)/(B+1)
  }

  if (simulate.p.value=="fixed") {
      method <- paste("Pearson's chi-squared test with simulated p-value, fixed row sums\n\t (based on",
                      B,"resamples)")
  }

  if (simulate.p.value=="random") {
      method <- paste("Pearson's chi-squared test with simulated p-value, marginal sums not fixed\n\t (based on",
                  B,"resamples)")
  }

  if (simulate.p.value==TRUE) {
    method <- paste("Pearson's chi-squared test with simulated p-value \n\t (based on",B,"resamples)")
  }


################################################
# create results object
###############################################

    results <- structure(list(
                            statistic=statistic,
                            parameter=parameter,
                            p.value=p.value,
                            method=method,
                            type=type,
                            observed=observed,
                            expected=expected,
                            residuals=residuals,
                            simulate.p.value=simulate.p.value,
                            sims=sims,
                            B=B,
                            verbose=verbose,
                            graph=graph),
                         class="GCchisqtest")

    return(results)
    
    
  }#end chisqtestGC


#' ##########################################################
#' #########################################################


#' @title Print Function for chisqtestGC

#' @description Utility print function
#' @keywords internal
#' 
#' @rdname print.GCchisqtest
#' @method print GCchisqtest
#' @usage 
#' \S3method{print}{GCchisqtest}(x,...)
#' @param x An object of class GCchisqtest.
#' @param \ldots ignored
#' @return Output to the console and to the plot window
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @export
print.GCchisqtest <- function(x,...)  {

  simulate.p.value <- x$simulate.p.value
  B <- x$B
  sims <- x$sims
  verbose <- x$verbose
  type <- x$type
  observed <- x$observed
  expected <- x$expected
  residuals <- x$residuals
  p.value <- x$p.value
  method <- x$method
  statistic <- x$statistic
  graph <- x$graph
  method <- x$method
  
  
  cat(method,"\n\n")
  
  if (verbose) {#print some tables
    
    if (type=="goodness"){
      nice <- cbind(observed,round(expected,2),round(residuals^2,2))
      colnames(nice) <- c("Observed counts","Expected by Null","Contr to chisq stat")
      print(nice)
    }
    
    if (type=="association"){
      cat("Observed Counts:\n")
      print(observed)
      cat("\n")
      cat("Counts Expected by Null:\n")
      print(round(expected,2))
      cat("\n")
      cat("Contributions to the chi-square statistic:\n")
      print(round(residuals^2,2))
    }
    cat("\n\n")
  }#end of verbose table printing
  
  #next, statistic, degrees of freedom and P-value
  cat("Chi-Square Statistic =",round(statistic,4),"\n")
  tab <- observed
  if (length(dim(tab))==1)  {df <- nrow(tab)-1} else {df <- (nrow(tab)-1)*(ncol(tab)-1)}
  cat("Degrees of Freedom of the table =",df,"\n")
  cat("P-Value =",round(p.value,4),"\n\n")
  
  #warn if no simulation and any expected cell counts are below 5
  if (min(expected) <5 && simulate.p.value==FALSE){
    cat(paste0("Some expected cell counts are low:",
            "\n\tthe approximation of the P-value may be unreliable.",
            "\n\tConsider using simulation."))
  }
  
  #finally, the graph
  if (graph && !(simulate.p.value==FALSE)) {
    hist(sims,xlab="Chi-Square Resamples",freq=TRUE,
         main=paste("Distribution of",B,"Resamples"),col="lightblue")
    abline(v=statistic,col="red",lwd=2)
  }
  
  if (graph==TRUE && simulate.p.value==FALSE) {
    invisible(pchisqGC(statistic,region="above",df=df,graph=T))
  }
  
}
