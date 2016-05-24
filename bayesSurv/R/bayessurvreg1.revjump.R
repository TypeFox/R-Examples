####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       bayessurvreg1.revjump.R        ####
####                                            ####
#### FUNCTIONS:  bayessurvreg1.revjump          ####
####################################################

### ======================================
### bayessurvreg1.revjump
### ======================================
## Subfunction for bayessurvreg1.R
##  -> just to make it more readable
##
## Manipulation with the specification of proposal jumps algorithms
##
bayessurvreg1.revjump <- function(prop.revjump)
{
  if (is.null(prop.revjump$algorithm)) prop.revjump$algorithm <- "basic"
  if (is.null(prop.revjump$moody.ring)) prop.revjump$moody.ring <- c(0.5, 0.5)
  if (is.null(prop.revjump$transform.split.combine)) prop.revjump$transform.split.combine <- "richardson.green"
  if (is.null(prop.revjump$transform.split.combine.parms)) prop.revjump$transform.split.combine.parms <-  c(2, 2, 2, 2, 1, 1)
  if (is.null(prop.revjump$transform.birth.death)) prop.revjump$transform.birth.death <- "richardson.green"


  ## Algorithm
  algorithm <- pmatch(prop.revjump$algorithm, table = c("basic", "independent.av", "correlated.av"), nomatch = 0, duplicates.ok = FALSE)[1]
  if (!algorithm) stop ("Unknown algorithm to generate canonical variables for reversible jumps.")
  algorithm <- algorithm - 1    ## NOW: 0 = basic, 1 = independent.av, 2 = correlated.av


  ## Moody ring parameters
  if (algorithm == 0){
    mrparm <- c(0.5, 0.5)
  }    
  else
    if (algorithm == 1){
      if(length(prop.revjump$moody.ring) < 1) stop("Incorrect prop.revjump$moody.ring parameter.")
      mrparm <- c(prop.revjump$moody.ring[1], 0.5)
    }
    else
      if (algorithm == 2){
        if(length(prop.revjump$moody.ring) < 2) stop("Incorrect prop.revjump$moody.ring parameter.")
        mrparm <- c(prop.revjump$moody.ring[1], 0.5)      
      }      
  names(mrparm) <- c("mr.epsilon", "mr.delta")
  if (sum(is.na(mrparm))) stop("Missing moody ring parameters.")
  if (sum(mrparm < 0 | mrparm > 0.5)) stop("Moody ring parameters must lie between 0 and 0.5.")


  ## Transformation for split-combine move
  transsc <- pmatch(prop.revjump$transform.split.combine,
                    table = c("richardson.green", "brooks", "identity"), nomatch = 0, duplicates.ok = FALSE)[1]
  if (!transsc) stop ("Unknown transformation for split-combine move.")
  transsc <- transsc - 1    ## NOW: 0 = richardson.green, 1 = brooks, 2 = identity 


  ## Split-combine transformation parameters
  if (transsc == 2){
    parmssc <- rep(1.0, 6)
  }
  else
    if (transsc == 0 || transsc == 1){
      if(length(prop.revjump$transform.split.combine.parms) < 6) stop("Incorrect prop.revjump$transform.split.combine.parms.")
      parmssc <- prop.revjump$transform.split.combine.parms[1:6]
    }      
  names(parmssc) <- paste("sc", 1:6, sep = "")
  if (sum(is.na(parmssc))) stop("Missing parameters for split-combine transformation.")
  if (sum(parmssc <= 0)) stop("Parameters for split-combine transformation must be positive.")


  ## Transformation for birth-death move
  transbd <- pmatch(prop.revjump$transform.birth.death,
                    table = c("richardson.green"), nomatch = 0, duplicates.ok = FALSE)[1]
  if (!transbd) stop ("Unknown transformation for birth.death move.")
  transbd <- transbd - 1    ## NOW: 0 = richardson.green

  ## Birth-death transformation parameters
  ## in this version: only allocate space of length 6
  parmsbd <- numeric(6)
  names(parmsbd) <- paste("bd", 1:6, sep = "")  

  revjumpi <- c(algorithm, transsc, transbd)
  names(revjumpi) <- c("algorithm", "transformation.sc", "transformation.bd")
  revjumpd <- c(mrparm, parmssc, parmsbd)

  rjdi <- list(integer = revjumpi, double = revjumpd)
  attr(rjdi, "prop.revjump") <- prop.revjump

  return(rjdi)    
}  
