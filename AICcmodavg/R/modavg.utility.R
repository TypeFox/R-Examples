##utility function to format candidate list of models to new class
##extract class from models in list and create new class
formatCands <- function(cand.set) {

  ##extract model class
  if(!is.list(cand.set)) stop("\n\'cand.set\' needs to be a list of candidate models\n")
  n.mods <- length(cand.set)
  all.mods <- lapply(cand.set, class)
  check.class <- unique(all.mods)
  out.class <- NULL
  ##for "coxph", c("coxph.null", "coxph"), c("clogit", "coxph")
  if(all(regexpr("coxph", check.class) != -1)) {
    out.class <- "coxph"
  }

  ##if NULL
  if(is.null(out.class)) {
    if(length(check.class) > 1) stop("\nFunctions do not support mixture of model classes\n")
    out.class <- unlist(check.class)
  }

  ##rename class
  mod.class.new <- c(paste("AIC", paste(out.class, collapse = "."), sep =""))
  
  ##add to list
  new.cand.set <- cand.set

  ##new S3 class
  class(new.cand.set) <- mod.class.new
  return(new.cand.set)
}





##utility functions used with modavg( ) to accomodate different specifications of interaction terms (e.g., A:B, B:A, A*B, B*A)
##in models of same set

####################################################
##function to reverse terms in interaction
reverse.parm <- function(parm) {

  ##check if ":" appears in term
  val <- grep(pattern = ":", x = parm)

  ##set value to NULL
  parm.alt <- NULL
  
  ##if ":" appears, then reverse interaction term
  if(length(val) > 0) {

    ##additional check if interaction involves more than 2 terms
    check.terms <- unlist(strsplit(x = parm, split = ":"))
  
    ##number of terms in interaction
    n.check.terms <- length(check.terms)

    ##issue warning if more than 2 terms are involved
    if(n.check.terms > 2) warning("\nThis function only supports two terms in an interaction:\n",
                                  "for more complex interactions, either create terms manually before analysis \n",
                                  "or double-check that models have been correctly included in model-averaging table\n")

    ##reverse order of interaction
    parm.alt.tmp <- rep(NA, n.check.terms)
    
    for (b in 1:n.check.terms) {
      parm.alt.tmp[b] <- check.terms[n.check.terms - b + 1]
    }

    ##paste terms together
    parm.alt <- paste(parm.alt.tmp, collapse = ":")

    return(parm.alt)
  }
}

#reverse.parm(parm = "BARE:AGE")




####################################################
##function to reverse order of exclude terms with colon or asterisk
reverse.exclude <- function(exclude) {

  ##remove all leading and trailing white space and within parm
  exclude <- lapply(exclude, FUN = function(i) gsub('[[:space:]]+', "", i))
    
  ##determine which terms are interactions with colons
  which.inter <- grep(pattern = ":", x = exclude)
  n.inter <- length(which.inter)

  ##list to hold reverse terms
  excl.list.alt <- list( )
  excl.list.alt2 <- list( )
  inter.star <- list( )
  
  ##if there are interaction terms with colons
  if (n.inter > 0) {
  
    ##create list for interaction
    rev.inter <- exclude[which.inter]

    ##create list to hold results
    excl.rev.list <- list( )
    for (b in 1:length(rev.inter)) {
      excl.rev.list[b] <- strsplit(x = rev.inter[b][[1]], split = ":")
    }

    ##add interaction with asterisk
    inter.star <- lapply(excl.rev.list, FUN = function(i) paste(i, collapse = " * "))
    
    ##additional check if interaction involves more than 2 terms
    n.check.terms <- unlist(lapply(excl.rev.list, length))
    
    ##issue warning if more than 2 terms are involved
    if(any(n.check.terms > 2)) warning("\nThis function only supports two terms in an interaction:\n",
                                       "for more complex interactions, either create terms manually before analysis \n",
                                       "or double-check that models have been correctly excluded in model-averaging table\n")

    ##iterate over each item in excl.rev.list
    for(k in 1:n.inter) {
      
      inter.id <- excl.rev.list[k][[1]]
      n.elements <- length(inter.id)
    
      ##reverse order of interaction
      parm.alt.tmp <- rep(NA, n.elements)
    
      for (b in 1:n.elements) {
        parm.alt.tmp[b] <- inter.id[n.elements - b + 1]
      }

      ##paste terms together
      excl.list.alt[k] <- paste(parm.alt.tmp, collapse = ":")
      excl.list.alt2[k] <- paste(parm.alt.tmp, collapse = " * ")
      
    }
  }


  ##determine which terms are interactions with asterisk
  which.inter.star <- grep(pattern = "\\*", x = exclude)
  n.inter.star <- length(which.inter.star)

  ##set lists to hold values
  inter.space <- list( )
  inter.nospace <- list( )
  ##list to hold reverse terms
  excl.list.alt.star <- list( )
  excl.list.alt.star2 <- list( )

  
  ##if there are interaction terms with asterisks
  if (n.inter.star > 0) {
  
    ##create list for interaction
    rev.inter <- exclude[which.inter.star]

    ##create vector to hold results
    excl.rev.list <- list( )
    for (b in 1:length(rev.inter)) {
      excl.rev.list[b] <- strsplit(x = rev.inter[b][[1]], split = "\\*")
    }

    ##paste interaction term with space
    inter.space <- lapply(excl.rev.list, FUN = function(i) paste(i, collapse = " * "))
    inter.nospace <- lapply(excl.rev.list, FUN = function(i) paste(i, collapse = ":"))
    
    ##additional check if interaction involves more than 2 terms
    n.check.terms <- unlist(lapply(excl.rev.list, length))
  
    ##issue warning if more than 2 terms are involved
    if(any(n.check.terms > 2)) warning("\nThis function only supports two terms in an interaction:\n",
                                       "for more complex interactions, either create terms manually before analysis \n",
                                       "or double-check that models have been correctly excluded in model-averaging table\n")
    
 
    ##iterate over each item in excl.rev.list
    for(k in 1:n.inter.star) {
      
      inter.id <- excl.rev.list[k][[1]]
      n.elements <- length(inter.id)
      
      ##reverse order of interaction
      parm.alt.tmp <- rep(NA, n.elements)
    
      for (b in 1:n.elements) {
        parm.alt.tmp[b] <- inter.id[n.elements - b + 1]
      }

      ##paste terms together
      excl.list.alt.star[k] <- paste(parm.alt.tmp, collapse = " * ")
      excl.list.alt.star2[k] <- paste(parm.alt.tmp, collapse = ":")
 
    }
  }


  ##add step to replicate each term with colon and asterisk
    
  ##combine into exclude
  exclude.out <- unique(c(exclude, excl.list.alt, excl.list.alt2, inter.space, inter.nospace, excl.list.alt.star, excl.list.alt.star2, inter.star))
  if(length(exclude.out) == 0) {exclude.out <- NULL}
  return(exclude.out)
}



##unexported functions
##create function for fixef to avoid importing nlme and lme4
#fixef <- function (mod){
  ##if from lme4
#  if(isS4(mod)) {
#    lme4::fixef(mod)
#  }
  
  ##if from coxme
#  if(identical(class(mod), "coxme")) {
#      mod$coefficients
#    }
  
  ##if from lmekin
#  if(identical(class(mod), "lmekin")) {
#    mod$coefficients$fixef
#  }

#  if(identical(class(mod), "lme")) {
#    nlme::fixef(mod)
#  ##if from nlme, coxme, lmekin
#  }
#}


##create function to identify REML models from lme4
isREML <- function(mod){
  as.logical(mod@devcomp$dims[["REML"]])
}


##extract ranef function
#ranef <- function (mod){
#  ##if from lme4
#  if(isS4(mod)) {
#    lme4::ranef(mod)
#  } else  {
#    mod$coefficients$random
#  } ##if from nlme, coxme, lmekin
#}
