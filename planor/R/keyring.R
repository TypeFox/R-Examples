#---------------------------------------------------------------------------
# CLASS "keyring" and its METHODS
#---------------------------------------------------------------------------
# An S4 class to represent a collection of key matrices associated with the same prime, with each key matrix a possible solution to the same model and estimate specifications. An object of class listofkeyrings is a list of keyring objects associated with the different primes involved in a given factorial design problem
# SLOTS
# -  .Data: a list of keymatrix objects associated with the same prime p and with the same factors
# -  p: a prime number
# - LIB: common row and column names of the key matrices
# - pseudo.info: a dataframe with one row per pseudofactor and with columns containing information on the factors (or pseudofactors) associated with the columns of the key matrices (e.g. 'parent', 'nlev', 'block' 'ordered', 'model', 'basic'; see the description of the class designfactors
# METHODS of "keyring": show, summary
#---------------------------------------------------------------------------
setClass("keyring",
         contains=c("list"),
         representation(p="numeric",
                        pseudo.info="data.frame",
                        LIB="list"))

##--------------------------------------------------------------------------
# "show.keyring" 
# Print an object of class keyring
# ARGUMENT
#  object: an object of class keyring
# RETURN
#  an invisible NULL.
# NOTE
# - The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
# - Non visible slot: pseudo.info
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# ## The method will now be used for automatic printing of a component of K0
# K0[[1]]
# show(K0[[1]]) ## idem
# print(K0[[1]]) ## idem
# ---------------------------------------------
show.keyring <- function(object){
  cat("An object of class keyring\n")
  cat(" Number of solutions:", length(object),
      "for prime", object@p, "\n\n")

  lapply(object, printgmat)

  invisible()
} ## end show.keyring


# --------------------------------------
# "show" method for "keyring"
# --------------------------------------
setMethod("show", signature(object="keyring"),
          definition=show.keyring)


##--------------------------------------------------------------------------
# "summary.keyring" 
# Summarizes the design properties from a keyring object, by
# printing the summary of each of its key matrices
# ARGUMENTS
#   - object: an object of class keyring
#   - show: optional string to identify the type of information to display.
#   - save: optional string to identify the type of information to return.
#   - ...: ignored
# NOTE
# The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#    nunits=3^3, base=~A+B+C, max.sol=2)
# summary(K0[[1]])
# ---------------------------------------------

summary.keyring <- function(object, show="tbw", save ="kw", ...){
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "summary" in R;

  ## Is some display required?
  isshow <-  (length(show) >0 && show != "" &&
    grepl("[d,t,b,w]", show, ignore.case=TRUE))

  ## Pseudofactors
  pseudo.info <- object@pseudo.info
  FACTtpf <- pseudo.info$parent
  NIVtpf <- pseudo.info$nlev
  LIBtpf <- rownames(pseudo.info)
  BLOCKtpf <- pseudo.info$block
  ## units factors
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

  p.k <- object@p
  if (isshow)
    printpower(p.k)
    ## Loop on the key matrices
    nsol <- length(object)
    sortie <- list()


      for(l in seq_len(nsol)){
        if (isshow)
          cat(paste("--- Solution ", l, " ---\n\n"))
       sortie[[l]] <- summary.keymatrix(object=object[[l]],
                        lib=colnames(object[[l]]),
                        fact=FACTtpf[ NIVtpf == p.k ],
                        block=BLOCKtpf[ NIVtpf == p.k ],
                          show, save)
    } ## end l
  if ( length(sortie) >0 ) {
    ## some outut are required
      names(sortie) <- paste("Solution", seq_len(nsol))
    return(invisible(sortie))
    }    else {      return(invisible())}

} ## end summary.keyring

# --------------------------------------
# "summary" method for "keyring"
# --------------------------------------
setMethod("summary", signature(object="keyring"),
          definition=summary.keyring)

