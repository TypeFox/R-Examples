#---------------------------------------------------------------------------
# CLASS "listofkeyrings" and its METHODS
# An S4 class to store design key solutions when there is only one prime involved or when the solutions are independent between primes
# SLOTS
# - .Data: a list of keyring objects associated with the different primes involved in the design under construction
# - factors: a designfactors object
# - model: a list of components of type c(model,estimate) containing the model and estimate specifications
# - nunits: the number of units of the design.
# METHODS of "listofkeyrings": "[" (or pick), planor.design, summary, show, alias
#---------------------------------------------------------------------------
setClass("listofkeyrings",
         contains=c("list"),
         representation(factors="designfactors",
                        model="list",
                        nunits="numeric"))
#---------------------------------------------------------------------------

# "pick.listofkeyrings" 
#   Extract a single designkey  object (with one key matrix per prime)
#  from an object of class listofkeyrings
# ARGUMENTS
# - keys: an object of class listofkeyrings
# - selection: index vector to select the key matrix for each prime
# RETURN
#   An object of class designkey, which contains  the selected design
# NOTE
#  K <- pick.listofkeyrings(K0,1) can be also be written
# K <- pick(K0,1) or more simply K <- K0[1]
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D, nunits=3^3, base=~A+B+C, max.sol=2)
# K0.1 <- pick(K0,1)##  (pick is a method of the class listofkeyrings)
# K0.1 <- K0[1] ## Another way of extracting ([ is synonym of pick)
# ------------------------------------------------
pick.listofkeyrings <- function(keys,selection){
  if(getOption("verbose")){
    cat( "Extraction of a design key from an object of class listofkeyrings\n")
  }
  Nuqf <- length(keys)
  if ((lgsel <- length(selection)) != Nuqf) {
    stop(
         paste("The length of the selection vector,",
               lgsel,
               ", is different from the number of primes,", Nuqf))
  }

  pickdesign <- vector("list",length=Nuqf)
  names(pickdesign) <- names(keys)

  for(k in seq_len(Nuqf)){
    pickdesign[[k]] <- keys[[k]]@.Data[[selection[k]]]
  }
  pickdesign <- new("designkey", .Data=pickdesign,
                    factors=keys@factors,
                    nunits=keys@nunits,
                    model=keys@model,
                    recursive=FALSE)
  return(pickdesign)
}
# ------------------------------------------------
# "[" Method to return the designkey object
#  of  index i for the  of prime, and index j for the second value, etc, ...,
#   from a listofkeyrings object.
# --------------------------------------
setMethod("[", "listofkeyrings",
          definition=function(x,i,j,...,drop){
            if (missing(j))
              x <- pick.listofkeyrings(x, c(i,...))
            else
              x <- pick.listofkeyrings(x, c(i,j,...))
            x
          })
# --------------------------------------
# "pick" method  for "listofkeyrings
# --------------------------------------
setMethod("pick", signature(keys="listofkeyrings"),
          definition=pick.listofkeyrings)

##------------------------------------------------------------------------
## "planor.design.listofkeyrings" 
# Build one design from an object of class listofkeyrings
# ARGUMENTS
# - key: an object of class listofkeyrings
# - randomize: an optional formula. When set, the final designs are randomized according to it.
# - selection: index vector to select the key matrix for each prime
# RETURN
#   An object of class  planordesign,
# which contains the design built from the selected key matrices
# NOTE
# Restricted to giving a single design
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#              model=~block+(A+B+C+D)^2, estimate=~A+B+C+D, nunits=3^3,
#             base=~A+B+C, max.sol=2, verbose=TRUE)
# P0 <- planor.design(key=K0, select=2) ##  (planor.design is a method of the class listofkeyrings)
# P0.R <- planor.design(key=K0, select=2, randomize=~A+B+C+D) ## Randomize the final designs
# -----------------------------------------------
planor.design.listofkeyrings <- function(key, randomize=NULL, selection, ...){
if(missing(selection)){ selection <- rep(1,length(key)) }
  selected <- pick.listofkeyrings(key,selection)
  OUT <- planor.design.designkey(key=selected, randomize=randomize, ...)
  return(OUT)
}

# --------------------------------------
# "planor.design" method for "listofkeyrings"
# --------------------------------------
setMethod("planor.design", signature(key="listofkeyrings"),
          definition=planor.design.listofkeyrings)
 ##------------------------------------------------------------------------

# "summary.listofkeyrings" 
# Summarizes the design properties from a listofkeyrings object, by
# printing the summary of each key matrix in each keyring
# ARGUMENTS
#  -  object: an object of class listofkeyrings
#  -  show: optional string to identify the type of information to display.
#  -  save: optional string to identify the type of information to return.
#  -  ...: ignored
# NOTE
# The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# print(summary(K0))
# F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=c(6,6,4,2,6) )
# M2 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
# K2 <- planor.designkey(factors=F2, model=M2, nunits=144,
#                        base=~A+B+D, max.sol=2)
# print(summary(K2))
# ---------------------------------------------

summary.listofkeyrings <- function(object, show="tbw", save="kw", ...){
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
  ## Loop on the distinct prime numbers
    sortie <- list()
  for(k in seq_len(Nuqf)){
    p.k <- PVuqf[k]
    if (isshow)
      printpower(p.k)
    ## Loop on the key matrices
    nsol <- length(object[[k]])
    sortie[[k]] <- list()
     for(l in seq_len(nsol)){
      if (isshow)
        cat(paste("--- Solution ", l, " for prime ", p.k, " ---\n\n"))
      retour <- summary.keymatrix(object=object[[k]][[l]],
                        fact=FACTtpf[ NIVtpf == p.k ],
                        block=BLOCKtpf[ NIVtpf == p.k ],
                        show, save)
      if (issave) {
        sortie[[k]][[l]]  <- retour
        names(sortie[[k]])[l]  <- paste("Solution", l, "for prime", p.k)
      }
    } ## end l
  } ## end k
  if ( issave) {
    names(sortie) <- paste("Solution", seq_len(Nuqf))
     return(invisible(sortie))
  }  else    return(invisible())
} ## end summary.listofkeyrings

# --------------------------------------
# "summary" method for "listofkeyrings
# --------------------------------------
setMethod("summary", signature(object="listofkeyrings"),
          definition=summary.listofkeyrings)


##--------------------------------------------------------------------------
# "show.listofkeyrings"
# Print the design key matrices of an object of class listofkeyrings
# ARGUMENT
#   object: an object of class listofkeyrings
# RETURN
# NOTE
# The number of rows and columns of the matrices that are printed
# are limited by the option \code{planor.max.print}
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#   nunits=3^3, base=~A+B+C, max.sol=2)
# ## The method will now be used for automatic printing of a component of K0
# K0
# show(K0) ## idem
# print(K0) ## idem
# ---------------------------------------------


show.listofkeyrings <- function(object){
  ## NOTE: the formal argument list "(object)" is
  ## required to be compatible with the generic function
  ##  in R;
  ## units factors
  cat("An object of class listofkeyrings\n")
  pseudo.info <- object@factors@pseudo.info

  
  PVuqf <- unique(factorize(object@nunits))

  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

  ## A. Design key matrices
  for(k in seq_len(Nuqf)){
    p.k <- PVuqf[k]
    printpower(p.k)
    nsol <- length(object[[k]])
    for(l in seq_len(nsol)){
      cat(paste("--- Solution ", l, " for prime ", p.k, " ---\n\n"))
      printgmat(object[[k]][[l]])
    }
  }

  invisible()
} ## end show.listofkeyrings
# --------------------------------------
# "show" method  for "listofkeyrings"
# --------------------------------------
setMethod("show", signature(object="listofkeyrings"),
          definition=show.listofkeyrings)


##--------------------------------------------------------------------------
# "alias.listofkeyrings" 
# Summarize the design properties from a listofkeyrings object.
# Return the factors, the model and the number of solutions for each prime.
# ARGUMENTS
#  -  object: an object of class listofkeyrings
#  -  model: an optional model formula (by default the first model in object)
#  - ...: ignored
# RETURN
#    A list indexed by the primes p of the object. Each element is a 3-column
#     matrix with one row per solution for prime p. The columns
#     give (i) the number of unaliased treatment effecs; (ii) the number of
#     mutually aliased treatment effects;  (iii) the number of treatment effects
#     aliased with block effects
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# alias(K0)
# F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=c(6,6,4,2,6) )
# M2 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
# K2 <- planor.designkey(factors=F2, model=M2, nunits=144,
#                        base=~A+B+D, max.sol=2)
# alias(K2)
# ---------------------------------------------

alias.listofkeyrings <- function(object, model, ...){
  if(missing(model)) model <- object@model[[1]][[1]]
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "alias" in R;

  ## Treatment factors
  object@factors <- object@factors[object@factors@fact.info$model]
  fact.info <- object@factors@fact.info
  Ntf <- nrow(fact.info)
  LIBtf <- rownames(fact.info)
  NIVtf <- fact.info$nlev
  BLOCKtf <- fact.info$block
  ## nbrs of pseudofactors per factor
  NPStf <- fact.info$npseudo          #FACT.npseudo
  ## Pseudofactors
  pseudo.info <- object@factors@pseudo.info
  FACTtpf <- pseudo.info$parent
  NIVtpf <- pseudo.info$nlev
  BLOCKtpf <- pseudo.info$block
  ## units factors
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

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

  ## output
  alias.stats <- vector(length=Nuqf,mode="list")
  names(alias.stats) <- PVuqf

  ## Loop on the distinct prime numbers
  for(k in seq_len(Nuqf)){
    p <- PVuqf[k]
    printpower(p)
    ## selection of the adequate rows and columns in the model and estimate matrices
    rows.k <- pseudo.info$nlev == p
    model.cols.k <-
      0 < apply(b.modset[rows.k, , drop=FALSE], 2, sum)
    model.k <- b.modset[rows.k, model.cols.k, drop=FALSE]
    ## alias calculations for prime p
    ## Loop on the key matrices
    nsol <- length(object[[k]])
    ## output
    alias.stats.p <- matrix(NA, nrow=nsol, ncol=3)
    colnames(alias.stats.p) <- c("unaliased","trt.aliased","blc.aliased")
    for(l in seq_len(nsol)){
      cat(paste("--- Solution ", l, " for prime ", p, " ---\n\n"))
      alias.stats.p[l,] <- alias.keymatrix(object=object[[k]][[l]],
                                         model=model.k,
                                         fact=FACTtpf[rows.k],
                                         block=BLOCKtpf[rows.k])
    } ## end
    cat("--- Synthesis on the aliased treatment effects for prime ",
        p, " ---\n\n")
    print(alias.stats.p)
    alias.stats[[k]] <- alias.stats.p
  } ## end k

  return(invisible(alias.stats))
} ## end alias.listofkeyrings
# -------------------------------------------------
# Method alias for "listofkeyrings"
setMethod("alias", signature(object="listofkeyrings"),
          definition=alias.listofkeyrings)

