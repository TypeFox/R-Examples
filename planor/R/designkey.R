#---------------------------------------------------------------------------
# CLASS designkey and its METHODS
#---------------------------------------------------------------------------

# An S4 class to represent a design key solution in the planor package
# SLOTS
# -  .Data: a list with one key matrix for each prime. Each one is an objet of class keymatrix
# -  factors: the 'designfactors' object that defines the factors
# -  model: the list of components of type c(model,estimate)
# -        containing the model and estimate specifications
# -  nunits: the number of units of the design
# -  recursive: logical, TRUE if the design has been constructed recursively
# Methods of "designkey" : planor.design, summary, show, alias
#---------------------------------------------------------------------------

setClass("designkey",
         contains=c("list"),
         representation(factors="designfactors",
                        model="list",
                        nunits="numeric",
                        recursive="logical"))
##---------------------------------------------------------------------------
## "planor.design.designkey" 
## ----------------------------------------------------
#  Build the design from a design key matrix
# ARGUMENTS
# - key: an object of class  designkey
# - randomize: an optional formula. When set, the final design is randomized according to it.
# RETURN
# An object of class planordesign,
# which  contains the design build from the input.
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#   nunits=3^3, base=~A+B+C, max.sol=2, verbose=TRUE)
# P0 <- planor.design(K0[1]) ##  (planor.design is a method of the class designkey)
# P0.R <- planor.design(K0[1], randomize=~A+B+C+D) ## randomize the final design
# -------------------------------------------------------
planor.design.designkey <- function(key, randomize=NULL, ...){
 key@factors <- key@factors[key@factors@fact.info$model]
    fact.info <- key@factors@fact.info
    Ntf <- nrow(fact.info)
    LIBtf <- rownames(fact.info)
    NIVtf <- fact.info$nlev
    ## Pseudofactors
    pseudo.info <- key@factors@pseudo.info
    FACTtpf <- pseudo.info$parent
    NIVtpf <- pseudo.info$nlev
    LIBtpf <- rownames(pseudo.info)
    Ntpf <- nrow(pseudo.info)

    ## A. Construction of each sub-design associated with each Sylow subgroup

    
    PVuqf <- unique(factorize(key@nunits))

    PVuqf <- PVuqf[order(PVuqf)]
    Nuqf <- length(PVuqf)
    b.pseudodesign.k <- vector("list", length=Nuqf)
    for(k in seq_len(Nuqf)){
        p.k <- PVuqf[k]
        r.k <- nrow(key[[k]])
      aux <-  crossing(rep(p.k,r.k),start=0)
      b.pseudodesign.k[[k]] <- (aux %*% key[[k]]) %%p.k
    }

    ## B. Crossing of the subdesigns
    b.fullpseudodesign <- cross.designs(b.pseudodesign.k)

    ## C. Reordering of the columns by treatment factor
    pseudosorder <- order(NIVtpf, FACTtpf, seq_along(NIVtpf))
      b.fullpseudodesign[,pseudosorder] <- b.fullpseudodesign

    ## D. Calculation of the design with the original treatment factors
       b.back <- matrix(0, nrow= Ntpf, ncol= Ntf)

    for(i in seq_len(Ntpf)){
        select <- (FACTtpf == FACTtpf[i]) & (seq_len(Ntpf) > i)
        b.back[i,FACTtpf[i]] <- prod( NIVtpf[select] )
    }
       b.finaldesign <- as.data.frame(b.fullpseudodesign %*% b.back)

    names(b.finaldesign) <- LIBtf ## columns names
    for(i in seq_len(Ntf)){
        zz <- factor(b.finaldesign[,i])
        levels(zz) <- (key@factors@levels)[[i]]
        b.finaldesign[i] <- zz ## convert column i into factor
    }

    ## E. Randomization, if required
    if (!is.null(randomize)) {
      b.finaldesign <- planor.randomize(randomize, b.finaldesign, ...)
      if( "InitialUNITS" %in% names(b.finaldesign) ){
        key@factors <- bind(planor.factors( factors="InitialUNITS",
                                           nlevels=nrow(b.finaldesign),
                                           dummy=FALSE),
                            key@factors)
      }
    }

    OUT <- new("planordesign",
               design=b.finaldesign,
               factors=key@factors,
               model=key@model,
               designkey=key@.Data,
               nunits= key@nunits,
               recursive=key@recursive )


    return(OUT)
} ## end planor.design.designkey
##------------------------------------------------------------------------
setMethod("planor.design", signature(key="designkey"),
          definition=planor.design.designkey)
##--------------------------------------------------------------------------


# "summary.designkey" 
# Summarises the design properties of a designkey object, by
# printing the summary of each of its key matrices (design key matrix, confounding
# and aliasing relationships)
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# resum <- summary(K0[1])
# ---------------------------------------------

summary.designkey <- function(object,show="dtbw", save="k", ...){
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "summary" in R;

  ## Is some display required?
  isshow <-  (length(show) >0 && show != "" &&
    grepl("[d,t,b,w]", show, ignore.case=TRUE))
  ## Is some output required?
  issave <-  (length(save) >0 && save != "" &&
    grepl("[k,w]", save, ignore.case=TRUE))

  ## Treatment factors
  object@factors <- object@factors[object@factors@fact.info$model]
  fact.info <- object@factors@fact.info
  Ntf <- nrow(fact.info)
  LIBtf <- rownames(fact.info)
  NIVtf <- fact.info$nlev
  BLOCKtf <- fact.info$block
  ## Pseudofactors
  pseudo.info <- object@factors@pseudo.info
  FACTtpf <- pseudo.info$parent
  NIVtpf <- pseudo.info$nlev
  BLOCKtpf <- pseudo.info$block
  ## units factors

  
  PVuqf <- unique(factorize(object@nunits))

  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

  ## Storage preparation of the returned results
  Hgen <- vector(mode="list",length=Nuqf)
  ## Loop on the distinct prime numbers
  for(k in seq_len(Nuqf)){
    p.k <- PVuqf[k]
    if (isshow)
      printpower(p.k)
    retour <- summary.keymatrix(object=object[[k]],
                                  fact=FACTtpf[ NIVtpf == p.k ],
                                  block=BLOCKtpf[ NIVtpf == p.k ],
                                   show, save)

    if (issave) {
      Hgen[[k]] <- retour
    }
  } ## end k

    if (issave) {
      names(Hgen)  <- paste("Prime",  PVuqf[1:Nuqf])
      return(invisible(Hgen))
    } else  return(invisible())
} ## end summary.designkey

##-------------------------------------------------------
setMethod("summary", signature(object="designkey"),
          definition=summary.designkey)
##-----------------------------------------------------------
# "show.designkey" 
# Print the design key matrices of an object of class designkey
#
# RETURN
# show returns an invisible NULL.
# NOTE
# - The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
#
# -  Objects of class designkey are displayed automatically by invocation of ‘show’
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D ,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# ## The method will now be used for automatic printing of a component of K0
# K0[1]
# show(K0[1]) ## idem
# print(K0[1]) ## idem
# ---------------------------------------------
show.designkey <- function(object){
  ## NOTE: the formal argument list "(object)" is
  ## required to be compatible with the generic function  in R;
  ## -----------------------------------------------
  cat("An object of class designkey\n")
  keys <- unclass(object)
  primes <- as.integer( names(keys) )
  

  ## A. Design key matrices
  for(k in seq(length(primes))) {
    printpower(primes[k])
    printgmat(keys[[k]])
  }
  invisible()
} ## end show.designkey

##--------------------------------------------------------------------------
setMethod("show", signature(object="designkey"),
          definition=show.designkey)

##--------------------------------------------------------------------------


# "alias.designkey" 
# Summarise the design properties from a design key matrix.
# Display the design keys matrices and the factorial effects confounded with the mean.
# ARGUMENTS
#    model: an optional model formula (by default the first model in object)
#    ...: ignored
# RETURN
#     invisible
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# alias(K0[1])
# ---------------------------------------------

alias.designkey <- function(object, model, ...){
  if(missing(model)) model <- object@model[[1]][[1]]
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "alias" in R;

  ## Treatment factors
  object@factors <- object@factors[object@factors@fact.info$model]
  fact.info <- object@factors@fact.info
  Ntf <- nrow(fact.info)              #FACT.N
  LIBtf <- rownames(fact.info)
  NIVtf <- fact.info$nlev             #FACT.nlev
  BLOCKtf <- fact.info$block
  ## nbrs of pseudofactors per factor
  NPStf <- fact.info$npseudo          #FACT.npseudo
  ## Pseudofactors
  pseudo.info <- object@factors@pseudo.info
  FACTtpf <- pseudo.info$parent       #PSEUDO.parent
  NIVtpf <- pseudo.info$nlev          #PSEUDO.nlev
  BLOCKtpf <- pseudo.info$block
  ## units factors
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]        #PRIMES
  Nuqf <- length(PVuqf)               #UQUASI.N

  ## Model  (cf. planor.modelterms)
  Mterms <- attributes(terms(model))
  ModelTerms <- Mterms$factors
  if(Mterms$response > 0){ ModelTerms <- ModelTerms[-1,] }
  if(max(ModelTerms)>1){
    stop("Sorry, nesting in the model formulae not possible") }
  ModelTerms <- cbind(MU=0, ModelTerms)
  mLIBtf <- rownames(ModelTerms)
  ## making of coherent rows between model and estimate
  Nmterms <- ncol(ModelTerms)
  ModelFine <- matrix(0, Ntf, Nmterms,
                      dimnames=list(LIBtf, colnames(ModelTerms)))
  ModelFine[mLIBtf,] <- ModelTerms[mLIBtf,]
  ## take off any all-zero column
  ModelFine <- ModelFine[,apply(ModelFine,2,function(x){sum(x)>0})]
   b.modterms <- ModelFine



  ## Decomposition into pseudofactors
  ## the function "planor.ineligibleset" can do that, assuming that
  ## we are in an "independent search" case
  b.modset <- planor.ineligibleset(object@factors, b.modterms)
  ## Loop on the distinct primes
  for(k in seq(Nuqf)){
    p <- PVuqf[k]
    printpower(p)
    ## selection of the adequate rows and columns in the model and estimate matrices
    rows.k <- pseudo.info$nlev == p
    modset <- matrix(b.modset[rows.k, ], ncol=ncol(b.modset))
    model.cols.k <-
      0 < apply(modset, 2, sum)

    model.k <- modset[, model.cols.k, drop=FALSE]
    ## alias calculations for prime p
    alias.keymatrix(object=object[[k]], model=model.k,
                    fact=FACTtpf[rows.k],
                    block=BLOCKtpf[rows.k])
  } ## end k

  return(invisible())
} ## end alias.designkey

##--------------------------------------------------------------
setMethod("alias", signature(object="designkey"),
          definition=alias.designkey)
##--------------------------------------------------
