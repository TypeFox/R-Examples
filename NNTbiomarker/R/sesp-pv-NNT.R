# The N's are Number Needed to Treat.  We have a range, such that  if NNT < NNTpos, all patients will be treated, if NNT > NNTneg then all patients will not be treated.
# Suppose N=NNTpos is the number of patients such that if N pts are positive,one will be a true positive.
#The "eq" means that we choose NNTpos so that treating all or not treating all would be equivalent.
#E(loss | treat) = (NNTpos-1) * L[A,H]  = E(loss | wait) = 1 * L[W,D]
#Actually we choose N SMALLER so that TREATing is definitely, comfortably the right thing.
#E(loss | treat) = (NNTpos-1) * L[A,H]  <<  E(loss | wait) = 1 * L[W,D]
# Suppose N=NNTneg is the number of patients such that if N pts are negative, one will be a false negative.
#The "eq" means that we choose NNTneg so that treating all or not treating any would be equivalent.
#E(loss | treat) = (NNTneg-1) * L[A,H]  = E(loss | wait) = 1 * L[W,D]
#Actually we choose N LARGER so that WAITing is definitely, comfortably the right thing.
#E(loss | treat) = (NNTpos-1) * L[A,H]  >> E(loss | wait) = 1 * L[W,D]

#' sesp.from.pv.feasible
#'
#' Computes sensitivity and specificity from  predictive values.
#'
#' @param ppv Positive predictive value
#' @param npv Negative predictive value
#' @param prev Prevalence (prior probability)
#' @param feasible Only return results in [0,1]. Default=TRUE
#' @return c(ppv=ppv, npv=npv, sp=sp, se=se)
#' @details  NNT stands for Number Needed to Treat.  We have a range, such that
#' if NNT < NNTpos, all patients will be treated, if NNT > NNTneg then all
#' patients will not be treated. Suppose N=NNTpos is the number of patients such
#' that if N pts are positive,one will be a true positive. The "eq" means that we
#' choose NNTpos so that treating all or not treating all would be equivalent.
#' E(loss | treat) = (NNTpos-1) * L[A,H]  = E(loss | wait) = 1 * L[W,D] Actually
#' we choose N SMALLER so that TREATing is definitely, comfortably the right
#' thing. E(loss | treat) = (NNTpos-1) * L[A,H]  <<  E(loss | wait) = 1 * L[W,D]
#' Suppose N=NNTneg is the number of patients such that if N pts are negative, one
#' will be a false negative. The "eq" means that we choose NNTneg so that treating
#' all or not treating any would be equivalent. E(loss | treat) = (NNTneg-1) *
#' L[A,H]  = E(loss | wait) = 1 * L[W,D] Actually we choose N LARGER so that
#' WAITing is definitely, comfortably the right thing. E(loss | treat) =
#' (NNTpos-1) * L[A,H]  >> E(loss | wait) = 1 * L[W,D]

sesp.from.pv.feasible = function(ppv, npv, prev, feasible=TRUE) {
  if(feasible) {   ### make sure that feassible results are returned.
    if(ppv < prev) {
      cat("ppv must be increased to prev\n")
      ppv = prev
    }
    if(npv < 1 - prev) {
      cat("npv must be increased to 1 - prev\n")
      npv = 1 - prev
    }
  }
  odds = prev/(1-prev)
  npv.odds = npv/(1-npv)
  ppv.odds = ppv/(1-ppv)
  saved.options = options(digits=4)
  # This is the contra-Bayes theorem.
  sp = npv.odds*(odds-ppv.odds)/(1-npv.odds*ppv.odds)
  se = ppv.odds/odds*(1-sp)
  return(c(ppv=ppv, npv=npv, sp=sp, se=se))
}

#' sesp.from.pv
#'
#' Computes sensitivity and specificity from  predictive values.
#'
#' @param ppv Positive predictive value
#' @param npv Negative predictive value
#' @param pv Alternative input of ppv and npv, in matrix or named vector.
#' @param prev Prevalence (prior probability)
#' @return c(se=se,sp=sp)
#' @aliases pv.to.sesp
#'
sesp.from.pv = pv.to.sesp = function(ppv=0.1, npv=0.7, pv, prev=0.2){
  if(!missing(pv)) {
    if(is.matrix(pv)) {
      if(identical(rownames(pv), c("ppv","npv"))) {
        ppv = pv["ppv", ];  npv = pv["npv", ]
      }
      else if(identical(colnames(pv), c("ppv","npv"))) {
        ppv = pv[ , "ppv"];  npv = pv[ , "npv"]
      }
      else stop("sesp.from.pv: matrix arg pv must have correct rownames or colnames.")
    }
    else if(is.vector(pv)) {
      ppv = pv["ppv"];   npv = pv["npv"]
    }
    else stop("sesp.from.pv: pv must be matrix or vector")
  }
  else {
    if(length(ppv) > 1 & missing(npv))
      warning("length(ppv) > 1 & missing(npv). check usage.")
    if(length(npv) > 1 & missing(ppv))
    warning("length(npv) > 1 & missing(ppv). check usage.")
  }
  ppv.odds= ppv/(1-ppv)
  npv.odds= npv/(1-npv)
  odds = prev/(1-prev)
  sp = as.vector(npv.odds*(odds-ppv.odds)/(1-npv.odds*ppv.odds))
  se = as.vector(ppv.odds/odds*(1-sp))
  if(length(sp) > 1) return(cbind(se=se,sp=sp))
  return(c(se=se,sp=sp))
}

#' pv.from.sesp
#'
#' Computes predictive values from sensitivity and specificity.
#'
#' @param se Positive predictive value
#' @param sp Negative predictive value
#' @param sesp Alternative input for se and sp, as matrix or named vector.
#' @param prev Prevalence (prior probability). Default =  0.001
#' @return c(ppv=ppv, npv=npv)
#' @aliases sesp.to.pv
#'
pv.from.sesp = sesp.to.pv = function(se=0.8, sp=0.8, sesp, prev=0.001) {
  if(!missing(sesp)) {
    if(is.matrix(sesp)) {
      if(identical(rownames(sesp), c("se","sp"))) {
        se = sesp["se", ];  sp = sesp["sp", ]
      }
      else if(identical(colnames(sesp), c("se","sp"))) {
        se = sesp[ , "se"];  sp = sesp[ , "sp"]
      }
      else stop("pv.from.sesp: matrix arg sesp must have correct rownames or colnames.")
    }
    else if(is.vector(sesp)) {
      se = sesp["se"];   sp = sesp["sp"]
    }
    else stop("pv.from.sesp: sesp must be matrix or vector")
  }
  else {
    if(length(se) > 1 & missing(sp))
      warning("length(se) > 1 & missing(sp). check usage.")
    if(length(sp) > 1 & missing(se))
      warning("length(sp) > 1 & missing(se). check usage.")
  }
  ppv = as.vector(prev*se/(prev*se + (1-prev)*(1-sp)))
  npv = as.vector((1-prev)*sp/(prev*(1-se) + (1-prev)*sp))
  if(length(ppv) > 1) return(cbind(ppv=ppv,npv=npv))
  return(c(ppv=ppv, npv=npv))
}


#' NNT.from.pv
#'
#' Compute NNT values from predictive values.
#'
#' @param ppv Positive predictive value
#' @param npv Negative predictive value
#' @param pv Alternative input of ppv and npv, in matrix or named vector.
#' @aliases pv.to.NNT
#' @return c(NNTpos=NNTpos, NNTneg=NNTneg)
#'
NNT.from.pv = pv.to.NNT = function(ppv, npv, pv) {
  if(!missing(pv)) {
    if(is.matrix(pv)) {
      if(identical(rownames(pv), c("ppv","npv"))) {
        ppv = pv["ppv", ];  npv = pv["npv", ]
      }
      else if(identical(colnames(pv), c("ppv","npv"))) {
        ppv = pv[ , "ppv"];  npv = pv[ , "npv"]
      }
      else stop("NNT.from.pv: matrix arg pv must have correct rownames or colnames.")
    }
    else if(is.vector(pv)) {
      ppv = pv["ppv"];   npv = pv["npv"]
    }
    else stop("NNT.from.pv: pv must be matrix or vector")
  }
  else {
    if(length(ppv) > 1 & missing(npv))
      warning("length(ppv) > 1 & missing(npv). check usage.")
    if(length(npv) > 1 & missing(ppv))
      warning("length(npv) > 1 & missing(ppv). check usage.")
  }
  NNTpos = 1/as.vector(ppv)
  NNTneg = 1/(1-as.vector(npv))
  if(length(ppv) > 1 | length(npv) > 1)
    return(list(NNTpos=NNTpos, NNTneg=NNTneg))
  return(c(NNTpos=NNTpos, NNTneg=NNTneg))
}

#' NNT.to.pv
#'
#' Convert between (NNTpos, NNTneg)  and (PPV, NPV).
#'
#' @param NNTpos NNT for a patient with positive test result
#' @param NNTneg NNT for a patient with negative test result
#' @param NNT A matrix or vector of (NNTpos, NNTneg) values.
#' @param prev Prevalence of the "BestToTreat" group before testing.
#' @param calculate.se.sp (default=FALSE) If TRUE, also calculate the sensitivity and specificity using the contra-Bayes theorem.
#' @aliases pv.from.NNT
#' @return For matrix input, cbind(ppv=ppv, npv=npv). For vector input, c(ppv=ppv, npv=npv).
#'
NNT.to.pv = pv.from.NNT = function(NNTpos, NNTneg, NNT, prev, calculate.se.sp=F) {
  if(!missing(NNT)) {
    if(is.matrix(NNT)) {
      if(identical(rownames(NNT), c("NNTpos","NNTneg"))) {
        NNTpos = NNT["NNTpos", ];  NNTneg = NNT["NNTneg", ]
      }
      else if(identical(colnames(NNT), c("NNTpos","NNTneg"))) {
        NNTpos = NNT[ , "NNTpos"];  NNTneg = NNT[ , "NNTneg"]
      }
      else stop("pv.from.NNT: matrix arg NNT must have correct rownames or colnames.")
    }
    else if(is.vector(NNT)) {
      NNTpos = NNT["NNTpos"];   NNTneg = NNT["NNTneg"]
    }
    else stop("pv.from.NNT: NNT must be matrix or vector")
  }
  else {
    if(length(NNTpos) > 1 & missing(NNTneg))
      warning("length(NNTpos) > 1 & missing(NNTneg). check usage.")
    if(length(NNTneg) > 1 & missing(NNTpos))
      warning("length(NNTneg) > 1 & missing(NNTpos). check usage.")
  }
  ppv.odds = 1/(NNTpos-1)
  ppv = as.vector(ppv.odds/(1+ppv.odds))
  npv.odds = (NNTneg-1)/1
  npv = as.vector(npv.odds/(1+npv.odds))
  if(length(ppv) > 1)
    return(cbind(ppv=ppv, npv=npv))
  return(c(ppv=ppv, npv=npv))
}

#' NNT.to.sesp
#'
#' Compute sensitivity aand specificity from NNT values.
#'
#' @param NNTpos NNT for a positive test result
#' @param NNTneg NNT for a negative test result
#' @param NNT Alternative way in input NNT values (matrix or vector)
#' @param prev Prevalence (prior probability)
#' @return c(se=se, sp=sp)
#' @aliases sesp.from.NNT
#'
NNT.to.sesp = sesp.from.NNT = function(NNTpos, NNTneg, NNT, prev) {
  if(!missing(NNT))
    sesp.from.pv(pv=pv.from.NNT(NNT=NNT), prev=prev)
  else
    sesp.from.pv(pv=pv.from.NNT(NNTpos, NNTneg), prev=prev)
}

#' NNT.from.sesp
#'
#' Compute NNT values from sensitivity aand specificity.
#'
#' @param se Positive predictive value
#' @param sp Negative predictive value
#' @param sesp Alternative input for se and sp, as matrix or named vector.
#' @param prev Prevalence (prior probability)
#' @return c(NNTpos=NNTpos, NNTneg=NNTneg)
#' @aliases sesp.to.NNT
#'

NNT.from.sesp = sesp.to.NNT = function(se, sp, sesp, prev) {
  if(!missing(sesp))
    pv.to.NNT(pv=sesp.to.pv(sesp=sesp, prev=prev) )
  else
    pv.to.NNT(pv=sesp.to.pv(se=se, sp=sp, prev=prev) )
}
