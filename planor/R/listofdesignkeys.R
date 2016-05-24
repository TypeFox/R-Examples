#---------------------------------------------------------------------------
# CLASS "listofdesignkeys" and its METHODS
#  S4 class, typically an output from planor.designkey when the research is recursive
# SLOTS
#  - .Data: a list of design-key solutions; each component
#  of main is a whole solution list across the different primes. It is an object of class designkey
#  - factors: the 'designfactors' object that defines the factors
#  - model: the list of components of type c(model,estimate)
#          containing the model and estimate specifications
#  - nunits: the number of units in the design.
#  METHODS of "listofdesignkeys": "[" (or pick), planor.design, summary, show, alias (not yet)
#---------------------------------------------------------------------------
setClass("listofdesignkeys",
         contains=c("list"),
         representation(factors="designfactors",
                        model="list",
                        nunits="numeric"))
#---------------------------------------------------------------------------
# "pick.listofdesignkeys" 
#   Extract a single designkey object (with one key matrix per prime)
#  from an object of class listofdesignkeys
# ARGUMENTS
# - keys: an object of class listofdesignkeys
# - selection: an integer equal to the position of the required solution
# RETURN
#  An object of class designkey, which contains  the selected design
# NOTE
#  K <- pick.listofdesignkeys(K0,1) can be also be written
# K <- pick(K0,1) or more simply K <- K0[1]
# EXAMPLES
# F2 <- planor.factors( factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2) )
# M2 <- planor.model( model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2 )
# K2 <- planor.designkey(factors=F2, model=M2, nunits=12,
#                       base=~R+C+U, max.sol=2)
# K2.1 <- pick(K2,1)
# K2.1 <- K2[1] ## Another way of extracting ([ is synonym of pick)
# ------------------------------------------------
pick.listofdesignkeys <- function(keys,selection){
  if(getOption("verbose")){
    cat( "Extraction of a design key from an object of class listofdesignkeys\n" )
  }
  if (selection > length(keys)) {
    stop( paste("The selection argument must be smaller than ", length(keys)))
  }



  pickdesign <- keys[[selection]]
  return(pickdesign)
}

# --------------------------------------
# "pick" method for "listofdesignkeys"
# --------------------------------------
setMethod("pick", signature(keys="listofdesignkeys"),
          definition=pick.listofdesignkeys)

#------------------------------------------------------------------------
# "[" method to return the designkey object of  the index solution
#   from a listofdesignkeys object.
# --------------------------------------
setMethod("[", "listofdesignkeys",
          definition=function(x,i,j,...,drop){
            if (missing(j))
              x <- pick.listofdesignkeys(x, c(i,...))
            else
              x <- pick.listofdesignkeys(x, c(i,j,...))
            x
          })


##------------------------------------------------------------------------
## "planor.design.listofdesignkeys"
## ---------------------------------------------------------------
# Build one design from an object of class listofdesignkeys
# ARGUMENTS
# - key: an object of class listofdesignkeys
# - randomize: an optional formula. When set, the final designs are randomized according to it.
# - selection: integer to select the solution
# RETURN
#  An object of class  planordesign,
# which contains the design built from the selected key matrices
# NOTE
# Restricted to giving a single design
# EXAMPLES
# K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2, nunits=12, base=~R+C+U, max.sol=2)
# P0 <- planor.design(key=K0, select=1)
# -----------------------------------------------
planor.design.listofdesignkeys <- function(key, randomize=NULL, selection=1, ...){
    selected <- pick.listofdesignkeys(key,selection)
    OUT <- planor.design.designkey(selected, randomize, ...)
    return(OUT)
}
# --------------------------------------
# "planor.design" method for "listofdesignkeys"
# --------------------------------------
setMethod("planor.design", signature(key="listofdesignkeys"),
          definition=planor.design.listofdesignkeys)


##--------------------------------------------------------------------------
# "summary.listofdesignkeys" 
# Summarizes the design properties of a listofdesignkeys object, by
# printing the summary of each key matrix in each design key (design key matrix, confounding and aliasing relationships)
# ARGUMENTS
#  - object: an object of class listofdesignkeys
#  -  show: optional string to identify the type of information to display.
#  -  save: optional string to identify the type of information to return.
#  -  ...: ignored
# NOTE
# The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
# EXAMPLES
# K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2, nunits=12, base=~R+C+U, max.sol=2)
# print(summary(K0))
# ---------------------------------------------

summary.listofdesignkeys <- function(object, show= "tbw", save="kw", ...){
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
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

## Loop on the solutions
  nsol <- length(object)
  sortie <- list()
  for(l in seq_len(nsol)){
    if (isshow)
    cat("\n********** Solution", l, "**********\n")
    ## Loop on the distinct prime numbers
    sortie[[l]] <- list()
    for(k in seq_len(Nuqf)){
      p.k <- PVuqf[k]
      if (isshow)
        cat(paste("--- Solution", l, "for prime", p.k, " ---\n\n"))
      retour <- summary.keymatrix(object=object[[l]][[k]],
                        fact=FACTtpf[ NIVtpf == p.k ],
                        block=BLOCKtpf[ NIVtpf == p.k ],
                        show, save)
      if (issave) {
        sortie[[l]][[k]] <- retour
        names(sortie[[l]])[k] <- paste("Solution", l, "for prime", p.k)
      }
    } ## end k
  } ## end l
  if ( issave) {
    names(sortie) <- paste("Solution", seq_len(nsol))
    return(invisible(sortie))
  }  else    return(invisible())

} ## end summary.listofdesignkeys


# --------------------------------------
# "summary" method  for "listofdesignkeys"
setMethod("summary", signature(object="listofdesignkeys"),
          definition=summary.listofdesignkeys)

##--------------------------------------------------------------------------
# "show.listofdesignkeys" 
# Print the design key matrices of an object of class listofdesignkeys
# ARGUMENT
#    object: an object of class listofdesignkeys
# RETURN
#    an invisible \sQuote{NULL}.
# NOTE
# The number of rows and columns of the matrices that are printed
# are limited by the option \code{planor.max.print}
# EXAMPLES
# K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2, nunits=12, base=~R+C+U, max.sol=2)
# K0
# show(K0) ## idem
# print(K0) ## idem
# ---------------------------------------------


show.listofdesignkeys <- function(object){
  ## NOTE: the formal argument list "(object)" is
  ## required to be compatible with the generic function
  ##  in R;
  ## units factors
  cat("An object of class listofdesignkeys\n")
  pseudo.info <- object@factors@pseudo.info
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

## Loop on the solutions
  nsol <- length(object)
  for(l in seq_len(nsol)){
    cat("\n********** Solution", l, "**********\n")
    ## Loop on the distinct prime numbers
    for(k in seq_len(Nuqf)){
      p.k <- PVuqf[k]
      cat(paste("--- Solution", l, "for prime", p.k, " ---\n\n"))
      printgmat(object[[l]][[k]])
    }
  }

  invisible()
} ## end show.listofdesignkeys
# --------------------------------------
# "show" method  for "listofdesignkeys
# --------------------------------------
setMethod("show", signature(object="listofdesignkeys"),
          definition=show.listofdesignkeys)

##--------------------------------------------------------------------------
# "alias.listofdesignkeys" 
# Summarize the design properties from a listofdesignkeys object.
# ARGUMENTS
#  - object: an object of class listofdesignkeys
#  - model: an optional model formula (by default the first model in object)
#  - ...: ignored
# RETURN
#    To see FUNCTION NOT YET IMPLEMENTED
# ---------------------------------------------
alias.listofdesignkeys  <- function(object, model, ...){
  stop("NOT YET IMPLEMENTED\n")
  
}
# --------------------------------------
# "alias" method for "listofdesignkeys"
# --------------------------------------
setMethod("alias", signature(object="listofdesignkeys"),
          definition=alias.listofdesignkeys)
