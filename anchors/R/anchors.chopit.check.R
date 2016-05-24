#######################################################################
##
## Function: anchors.chopit.check()
## Author  : Jonathan Wand <wand(at)stanford.edu>
##           http://wand.stanford.edu
## Created :  2008-04-20
##
## Extracted and refined from former chopit()  
#######################################################################
anchors.chopit.invocation <- function (x) {
    cat("\n\nInvocation options used:\n")
    cat("  normalize   = ",x$normalize ,"\n")
    cat("  random      = ",x$random     ,"\n")
    cat("  linear      = ",x$linear,"\n")
    cat("  method      = ",x$method,"\n")
    cat("  analytic    = ",x$analytic    ,"\n")
    cat("\n")
}

anchors.chopit.check <- function(count, options) {

  if (count$n.vign.set > count$n.self & count$n.self > 0) {
    cat("\nERROR! Number of sets of vignettes exceeds number of self-Q.\n")
    cat(  "ERROR! This does not make any sense. \n")
    cat(  "ERROR! Exiting chopit() function! \n\n")
    stop("")
  }
  ## jw: currently obsolete: count$n.vign.set == 1
  if (count$n.vign.set > 1 & count$n.self == 0 & options$vign.cut == "hetero") {
    cat("\nERROR! CHOPIT() only estimate unconstrained vignette taus one vignette set at a time.\n")
    cat(  "ERROR! Exiting chopit() function! \n\n")
    stop("")
  }
  
  ## Q. should we estimate random effect?
  ## A. If asked for, check to make sure that more than one self-q
  if (options$random & count$n.self < 2) {
    cat("\nWARNING! Cannot estimate random effect with only one self-question\n")
    cat(  "WARNING! DISABLING random effect estimation for this run.\n\n")
    options$random <- FALSE
  }

  ## Q. Can we estimate RE and use linear specification?
  ## A. Currently, no.
  if (options$random & options$linear) {
    cat("\nWARNING! Cannot use linear taus with random effect at this time.\n")
    cat(  "WARNING! Keeping random effect estimation,\n")
    cat(  "WARNING!   BUT SWITCHING to NON-LINEAR taus for this run.\n\n")
    options$linear <- FALSE
  }
  
  ## Q. How many categories?
  ## A. Look at highest index/category
  ##    used in ALL self-questions and 1st set of vignettes
  ##    EXCLUDING ZERO CATEGORY -- missing data

  if (is.na(count$n.cat) || !is.finite(count$n.cat)) {
    cat("ERROR! There is an invalid value in a response variable.\n")
    cat("Here is the list of unique values:",
        sort(count$obs.cat,na.last=TRUE),"\n")
    stop("")
  }
  ## LOOK FOR HOLES IN THE RESPONSE VECTOR
  ## -- this would lead to unidentified threshold/taus
  if ( length(count$obs.cat[count$obs.cat!=0]) != count$n.cat ) {
    cat("WARNING! One or more categories is not used,\n")
    cat("(possibly following listwise deletion of missing values)\n")
    cat("This results in an unidentified threshold---STOPPING\n\n")
    cat("Here is the list of unique values:",
        sort(count$obs.cat,na.last=TRUE),"\n\n")
    stop("")
  }
  if (options$verbose) {
    cat("\nCOUNTING OBSERVATIONS:\n",
        "\nNumber of cases in self likelihood (nobs.self) = ",count$nobs.self,
        "\nNumber of cases in vign likelihood (nobs.vign) = ",count$nobs.vign,
        "\n")
        
    cat("\nCOUNTING TYPES:\n",
        "\nNumber of self-Q (n.self) = ",count$n.self,
#        "\nNumber of sets of vignettes (n.vign.set) = ",count$n.vign.set,
        "\nNumber of categories (n.cat) = ",count$n.cat,
        "\nNumber of variables in gamma (nvars.gamma) = ",count$nvars.gamma,
        "\nNumber of variables in gamma1 (nvars.gamma1) = ",count$nvars.gamma1,
        "\nNumber of vignettes  (n.vign)",count$n.vign,"\n")
#    cat("Mapping between vignettes to self-Q taus (vign.map):\n")
#    cat("  ",vign.map,"\n")
  }


  ## analytical gradient limitations
  if (options$analytic) {
                                        #    if (normalize != "self") {
                                        #      cat("\nAnalytical gradients currently available only for 'normalize=\"self\"'\n")
                                        #      cat("Setting 'analytic=FALSE'\n")
                                        #      analytic <- FALSE
                                        #    }
    if (options$linear != TRUE) {
      cat("\nAnalytical gradients currently available only for 'linear=TRUE'\n")
      cat("Setting 'analytic=FALSE'\n")
      options$analytic <- FALSE
    }
    if (options$random == TRUE) {
      cat("\nAnalytical gradients currently available only for 'random=FALSE' \n")
      cat("Setting 'analytic=FALSE'\n")
      options$analytic <- FALSE
    }
  }
  return(options)
}
