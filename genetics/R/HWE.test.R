# $Id: HWE.test.R 1225 2007-05-29 19:29:10Z warnes $

### Hardy-Weinberg Equilibrium Disequlibrium Estimates, Confidence
### Intervals, and P-values
###



HWE.test <- function(x, ...)
{
	UseMethod("HWE.test")
}


HWE.test.genotype <- function(x, exact=nallele(x)==2,
                              simulate.p.value=!exact, B=10000,
                              conf=0.95, ci.B=1000, ... )
  # future options "bootstrap","exact"
{

  retval <- list()

  # compute disequlibrium
  retval$diseq <- diseq(x)
  
  # compute confidence intervals
  retval$ci <- diseq.ci(x, R=ci.B, conf=conf)

  # do chisq test
  if(exact)
    retval$test  <- HWE.exact(x)
  else
    {
      tab <- retval$diseq$observed.no
      tab  <- 0.5 * (tab + t(tab))   # make symmetric for chisq.test
      retval$test  <- HWE.chisq(x, simulate.p.value=simulate.p.value,B=B,...)
    }

  
  retval$simulate.p.value <- simulate.p.value
  retval$B <- B
  retval$conf <- conf
  retval$ci.B <- ci.B
  retval$test$data.name  <- deparse(substitute(x))
  retval$call  <- match.call()
  class(retval)  <-  c("HWE.test")
  return(retval)
}



print.HWE.test  <-  function(x, show=c("D","D'","r","table"), ...)
  {

    cat("\n")
    cat("\t-----------------------------------\n")
    cat("\tTest for Hardy-Weinberg-Equilibrium\n")
    cat("\t-----------------------------------\n")
    cat("\n")
    if(!is.null(x$locus))
      {
        cat("\n")
        print( x$locus )
      }
    cat("Call: \n")
    print(x$call)
    cat("\n")
    if("D" %in% show)
      {
        cat("Raw Disequlibrium for each allele pair (D)\n")
        cat("\n") 
        print(x$diseq$D)
        cat("\n")
      }
    if("D'" %in% show)
      {
        cat("Scaled Disequlibrium for each allele pair (D')\n")
        cat("\n") 
        print(x$diseq$Dprime)
        cat("\n")
      }
    if("r" %in% show)
      {
        cat("Correlation coefficient for each allele pair (r)\n")
        cat("\n") 
        print(x$diseq$r)
        cat("\n")
      }
    if("table" %in% show)
      {
        cat("Observed vs Expected Allele Frequencies \n")
        cat("\n")
        print(x$diseq$table)
        cat("\n")
      }
    
    if( ncol(x$diseq$r) <= 2 )
      cat("Overall Values\n")
    else
      cat("Overall Values (mean absolute-value weighted by expected allele frequency)\n")
    cat("\n")

    show.tab <- NULL
    
    if("D" %in% show)
      show.tab <- rbind(show.tab, "  D"=x$diseq$D.overall)
    if("D'" %in% show)
      show.tab <- rbind(show.tab, "  D'"=x$diseq$Dprime.overall)
    if("r" %in% show)
      show.tab <- rbind(show.tab, "  r"=x$diseq$r.overall)

    colnames(show.tab) <- "Value"

    print(show.tab)
    
    cat("\n") 

    whichvec <- c("D","D'","r") %in% show

    cat("Confidence intervals computed via bootstrap using", x$ci.B, "samples\n")
    cat("\n")

    if(!is.null(x$ci$warning.text))
      cat(strwrap(paste("WARNING:", x$ci$warning.text), prefix="    * "),"\n",
          sep="\n")
    
    show.tab <- matrix(ncol=4, nrow=4)
    tmp <- format(x$ci$ci[,1:3], digits=getOption("digits"))
    show.tab[,1] <- tmp[,1]  # Observed
    show.tab[,2] <- paste("(", tmp[,2], ", ", tmp[,3], ")", sep="" )
    show.tab[,3] <- x$ci$ci[,4]
    show.tab[,4] <- ifelse(x$ci$ci[,5],"YES","*NO*")

    colnames(show.tab) <- c("Observed", "95% CI", "NA's", "Contains Zero?")
    rownames(show.tab) <- paste("  ", rownames(tmp), sep="")
    
    print(show.tab[whichvec,], quote=FALSE)

    cat("\n")
    cat("Significance Test:\n")
    print(x$test)
    cat("\n")
    cat("\n")
  }

