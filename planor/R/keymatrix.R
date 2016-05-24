
#---------------------------------------------------------------------------
# CLASS "keymatrix" and its METHODS
# An S4 class to represent an elementary key matrix in the planor package
# SLOTS
#  - .Data: a matrix of integers modulo \note{p}
#  -  p: a prime number
# METHODS
#    summary, alias, show
# NOTE
# An object of class designkey or keyring is a collection of keymatrix objects
#---------------------------------------------------------------------------
setClass("keymatrix",
         contains=c("matrix"),
         representation(p="numeric"))
##--------------------------------------------------------------------------


# "summary.keymatrix"
# Summarises the main properties of a single key matrix, by displaying the key matrix, the factorial effects confounded with the mean, and the weight profiles of the effects confounded with the mean
#
# ARGUMENTS
#  - object: a key matrix denoted by 'key'
#  - fact: a character or numeric vector of parent factor names for the columns of 'key'
#  - block: a logical vector to identify the columns of 'key' associated with a block factor
#  - ...: ignored
# NOTE
# - The rows of key are associated with units factors (or pseudofactors)
# while its columns are associated with treatment or block factors (or pseudofactors).
# The vectors in the arguments fact and block give information on the treatment and block factors,
# so that their length is expected to be equal to the number of columns of key.
# If missing, fact attributes a distinct parent factor to each column of key and block is set to TRUE for all columns
#
# - The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
# Return
#     A matrix whose columns are generators of the kernel of 'key'
# ---------------------------------------------

summary.keymatrix <- function(object, fact, block,
                              show="dtbw", save="k", ...){

  ## missing arguments
  if(missing(fact)){ fact <- seq(ncol(object)) }
  if(missing(block)){ block <- rep(FALSE, ncol(object)) }

  ## missing row or column names
  if(is.null(rownames(object))){
    rownames(object) <- paste("U",seq(nrow(object)),sep="") }
  if(is.null(colnames(object))){
    colnames(object) <- LETTERS[seq(ncol(object))] }
  ## maximal number of rows and columns that are printed
  maxprint <- getOption("planor.max.print", default=20)

  ## Check some arguments
  if (length(show) >0 && show != "") {
    if (!grepl("[d,t,b,w]", show, ignore.case=TRUE)) {
    ## no character recognized
    warning("No character recognized in argument show. Valid characters are 'd', 't', 'b', 'w' or NULL")
    show <- ""
  }
  } else show <- ""
  if (length(save) >0 && save != "") {
    if ( !grepl("[k,w]",save, ignore.case=TRUE)) {
    ## no character recognized
    warning("No character recognized in argument save. Valid characters are 'k', 'w', or NULL")
    save <- ""
  }
  } else save <- ""

  sortie <- list()
  ## A. Design key matrices
  if (grepl('d', show, ignore.case=TRUE)) {
    cat("DESIGN KEY MATRIX\n")
    printgmat(object )
  }
  p <- object@p
  lib <- colnames(object)
  ## B. Design key kernels
  ## kernel calculation
  b.Hgen <- kernelmatrix.basep(object,p)
  ## b.Hgen is  a matrix with 0 column
  if (ncol(b.Hgen) == 0) {
    warning("Regular matrix: no confounding\n")
    
    save <- sub('w', '', save) 
  } else {
    b.H <- subgroup.basep(b.Hgen,p)
    ## b.H is a matrix
    b.H <- weightorder.basep(b.H, p, fact, block)
    ## printing
    ## no printing of b.H for the moment
    
    ## selection of the columns free from block effects
    selectCol <- apply( b.H[block, , drop=FALSE], 2, sum ) == 0
    ## PLANORlibsk ne modifie rien, ne fait qu'écrire
    if (grepl('t', show, ignore.case=TRUE)) {
      cat("TREATMENT EFFECTS CONFOUNDED WITH THE MEAN\n")
      if(sum(selectCol) > 0){
        H.show <- matrix(b.H[,selectCol], nrow=nrow(b.H))
         .Call("PLANORlibsk", as.integer(nrow(H.show)),
            as.integer(ncol(H.show)),
            H.show,
            as.character(lib),
            as.integer(maxprint))
        cat("\n")
      }
      else cat("nil\n\n")
    } # end grepl
    ## remaining info if there are block factors
    if (grepl('b', show, ignore.case=TRUE)) {
    cat("BLOCK-and-TREATMENT EFFECTS CONFOUNDED WITH THE MEAN\n")
    if(sum(block) > 0){
      H.show <- matrix(b.H[,!selectCol], nrow=nrow(b.H))
       .Call("PLANORlibsk", as.integer(nrow(H.show)),
            as.integer(ncol(H.show)),
            H.show,
            as.character(lib),
            as.integer(maxprint))
      cat("\n")
    }
    else cat("nil\n\n")
  } # end grepl
    ## C. Design key kernels
    if (grepl('w', show, ignore.case=TRUE) ||
        grepl('w', save, ignore.case=TRUE)) {


    Wprint <- rbind(attributes(b.H)$trt.weight,
                    attributes(b.H)$trt.pseudoweight,
                    attributes(b.H)$blc.weight,
                    attributes(b.H)$blc.pseudoweight)

    sortiew<-vector("character")


    if (any(selectCol)) {
      sortiew[1] <- paste(wprofile(Wprint[1,selectCol]), collapse=" ")
    } else {
      sortiew[1] <- "none"
    }
    if (grepl('w', show, ignore.case=TRUE)) {
    cat("WEIGHT PROFILES\n")
    cat("Treatment effects confounded with the mean: ")
    cat(sortiew[1],"\n")
  }

    if (any(!selectCol)) {
      sortiew[2] <- paste(wprofile(Wprint[1,!selectCol]), collapse=" ")
      } else {
         sortiew[2] <- "none"
     }
      if (grepl('w', show, ignore.case=TRUE)) {
    cat("Treatment effects confounded with block effects: ")
    cat(sortiew[2],"\n")

    cat("Treatment pseudo-effects confounded with the mean: ")
    if (any(selectCol)) {
        cat(wprofile(Wprint[2,selectCol]),"\n")} else {
      cat("none\n")
    }

    cat("Treatment pseudo-effects confounded with block effects: ")
    if (any(!selectCol)) {
      cat(wprofile(Wprint[2,!selectCol]),"\n")
    } else {
      cat("none\n")
      }
} # end (grepl('w', show, ignore.case=TRUE)
  } # end grepl on show and save
    cat("\n")
  } ## end  else Hgen with 0 column

    if (grepl('k', save, ignore.case=TRUE)) {
      sortie$k <- b.Hgen
    }
    if (grepl('w', save, ignore.case=TRUE)) {
      sortie$w <- sortiew
      names(sortie$w) <- c("Treatment effects confounded with the mean",
                           "Treatment effects confounded with block effects")
    }

  if (save != "")
    return(invisible(sortie))
  else
    return(invisible())

} ## end summary.keymatrix

# --------------------------------------
# "summary" method for "keymatrix"
# --------------------------------------
 setMethod("summary", signature(object="keymatrix"),
          definition=summary.keymatrix)



##-----------------------------------------------------


# "alias.keymatrix"
# Summarises the aliasing properties of a single key matrix
# ARGUMENTS
# - object: a key matrix denoted by 'key'
# - model: a matrix representing factorial model terms
# - fact: a character or numeric vector of parent factor names for the columns of 'key'
# - block: a logical vector to identify the columns of 'key' associated with a block factor
# - ...: ignored
# NOTE
# The function prints the unaliased treatment effects,
# then the groups of aliased treatment effects, then the treatments effects confounded with block effects
# and finally the unaliased block effects, when considering all the factorial terms that are represented in the \code{model} argument,
# which is set if missing to the identity matrix (main effects only)
# RETURN
#    A vector with
#     (i) the number of unaliased treatment effects; (ii) the number of
#     mutually aliased treatment effects; (iii) the number of treatment effects
#     aliased with block effects
# ---------------------------------------------

alias.keymatrix <- function(object, model,  fact, block, ...){
  
  
  ## missing arguments
  if(missing(fact)){ fact <- seq(ncol(object)) }
  if(missing(block)){ block <- rep(FALSE, ncol(object)) }
  if(missing(model)){ model <- diag(ncol(object)) }
  ## missing row or column names
  if(is.null(rownames(object))){
    rownames(object) <- paste("U",seq(nrow(object)),sep="") }
  if(is.null(colnames(object))){
    colnames(object) <- LETTERS[seq(ncol(object))] }

  p <- object@p
  lib <- colnames(object)
  ## expansion of the model and estimate matrices
  b.model.full <- representative.basep( model, p )
  Nfull <- ncol(b.model.full)

  ## images of the model terms by the key matrix
    b.images.mat <- (object %*% b.model.full) %% p

  ## info on the columns of b.model.full and on the columns of b.images.mat
  trt.log <- rep(NA,Nfull)    # effect with at least one trt pseudofactor
  blc.log <- rep(NA,Nfull)    # effect with at least one blc pseudofactor
  first.ind <- rep(NA,Nfull)  # row index of first pseudofactor involved
  first.val <- rep(1,Nfull)  # coeff of first pseudofactor involved
  first.inv <- rep(1,Nfull)  # inverse of coeff of 1st pseudof involved
  ## Loop on the model effects
  for(j in seq(Nfull)){
    ## identification of the effects including treatment and block pseudofactors
    nonzero.model.j <- b.model.full[,j] != 0
    trt.log[j] <- any( (!block)[nonzero.model.j] )
    blc.log[j] <- any( block[nonzero.model.j] )
    ## normalization of the first non-zero element of the image vectors, if needed
    nonzero.j <- b.images.mat[,j] != 0

    if( any(nonzero.j) & (p!=2) ){

      first.ind[j] <- seq(nrow(b.images.mat))[nonzero.j][1]
      first.val[j] <- b.images.mat[first.ind[j], j]
      ## Raw calculation of the inverses modulo p




      inv <- 0
      inv <- .C("PLANORinv", as.integer(p),
                as.integer(first.val[j]),
                inv= as.integer(inv))$inv
      first.inv[j] <- inv

        b.images.mat[, j] <- (first.inv[j] * b.images.mat[, j]) %%p

        b.model.full[, j] <- (first.inv[j] * b.model.full[, j]) %%p
    } ## end of the normalization
  } ## end of the loop on j
  ## alias.mat construction : matrix with one row per set of confounded effects
  codes <- convertfrom.basep(t(b.images.mat[,]), p)
  codes.u <- unique(codes)
  alias.mat <- matrix(0, length(codes.u), Nfull)
  for(i in seq(length(codes.u))){
    alias.mat[i, codes==codes.u[i]] <- first.inv[codes==codes.u[i]]
  }

  ## characterization of the alias set (rows of alias.mat)
  nb.blc <- c((alias.mat > 0) %*% blc.log)
  nb.trt <- c((alias.mat > 0) %*% trt.log)
  nb.terms <- apply((alias.mat > 0), 1, sum)

  ## initialisation of the output
  alias.stats <- vector(length=3)
  names(alias.stats) <- c("unaliased","trt.aliased","blc.aliased")

  cat("UNALIASED TREATMENT EFFECTS\n")
  select.alias <- (nb.blc == 0) & (nb.terms == 1)
  alias.terms.selected <- alias.mat[select.alias, , drop=FALSE]
  select.terms <- apply(alias.terms.selected, 2, sum) > 0
  alias.stats[1] <- sum(select.terms)
  ## Proceed if there is at least one such term
  if(any(select.terms)){
    model.selected <- b.model.full[,select.terms]
    nb.print <- ncol(model.selected)
    ## Loop on the selected terms
    for(j in seq(nb.print)){
      if(j > 1) cat (" ; ")
      psfact.coeffs <- model.selected[,j]
      psfact.select <- psfact.coeffs > 0
      toprint <- lib[psfact.select]
      psfact.coeffs.sel <- psfact.coeffs[psfact.select]
      if(max(psfact.coeffs > 1)){
        psfact.select2 <- psfact.coeffs.sel > 1
        psfact.pwr <- paste("^",psfact.coeffs.sel[psfact.select2],sep="")
        toprint[psfact.select2] <- paste(toprint[psfact.select2],
                                         psfact.pwr,sep="")
      }
      cat(toprint,sep=":")
    }
    cat("\n\n")
  }
  else cat("nil\n\n")

  cat("ALIASED TREATMENT EFFECTS\n")
  select.alias <- (nb.blc == 0) & (nb.terms > 1)
  ## Proceed if there is at least one alias
  if(any(select.alias)){
    alias.selected <- alias.mat[select.alias,,drop=FALSE]
    nb.alias.selected <- nrow(alias.selected)
    alias.stats[2] <- sum(alias.selected > 0)
    ## Loop on the aliases
    for(a in seq(nb.alias.selected)){
      select.terms <- as.logical(alias.selected[a,])
      model.selected <- b.model.full[,select.terms]
      nb.print <- ncol(model.selected)
      ## Loop on the selected terms
      for(j in seq(nb.print)){
        if(j > 1) cat (" = ")
        psfact.coeffs <- model.selected[,j]
        psfact.select <- psfact.coeffs > 0
        toprint <- lib[psfact.select]
        psfact.coeffs.sel <- psfact.coeffs[psfact.select]
        if(max(psfact.coeffs > 1)){
          psfact.select2 <- psfact.coeffs.sel > 1
          psfact.pwr <- paste("^",psfact.coeffs.sel[psfact.select2],sep="")
          toprint[psfact.select2] <- paste(toprint[psfact.select2],
                                           psfact.pwr,sep="")
        }
        cat(toprint,sep=":")
      }
      cat("\n")
    }
    cat("\n")
  }
  else cat("nil\n\n")

  cat("TREATMENT AND BLOCK EFFECTS CONFOUNDED WITH BLOCK EFFECTS\n")
  select.alias <- (nb.blc > 0) & (nb.terms > 1)
  alias.stats[3] <- 0

  ## Proceed if there is at least one alias
  if(any(select.alias)){
    alias.selected <- alias.mat[select.alias, , drop=FALSE]
    nb.alias.selected <- nrow(alias.selected)
    ## Loop on the aliases
    for(a in seq(nb.alias.selected)){
      select.terms <- as.logical(alias.selected[a,])
      alias.stats[3] <- alias.stats[3] + sum(select.terms[blc.log==0])
      model.selected <- b.model.full[,select.terms]
      nb.print <- ncol(model.selected)
      ## Loop on the selected terms
      for(j in seq(nb.print)){
        if(j > 1) cat (" = ")
        psfact.coeffs <- model.selected[,j]
        psfact.select <- psfact.coeffs > 0
        toprint <- lib[psfact.select]
        psfact.coeffs.sel <- psfact.coeffs[psfact.select]
        if(max(psfact.coeffs > 1)){
          psfact.select2 <- psfact.coeffs.sel > 1
          psfact.pwr <- paste("^",psfact.coeffs.sel[psfact.select2],sep="")
          toprint[psfact.select2] <- paste(toprint[psfact.select2],
                                           psfact.pwr,sep="")
        }
        cat(toprint,sep=":")
      }
      cat("\n")
    }
    cat("\n")
  }
  else cat("nil\n\n")

  cat("UNALIASED BLOCK EFFECTS\n")
  select.alias <- (nb.trt == 0) & (nb.terms == 1)
  alias.terms.selected <- alias.mat[select.alias, , drop=FALSE]
  select.terms <- apply(alias.terms.selected, 2, sum) > 0
  ## Proceed if there is at least one such term
  if(any(select.terms)){
    model.selected <- b.model.full[,select.terms]
    nb.print <- ncol(model.selected)
    ## Loop on the selected terms
    for(j in seq(nb.print)){
      if(j > 1) cat (" ; ")
      psfact.coeffs <- model.selected[,j]
      psfact.select <- psfact.coeffs > 0
      toprint <- lib[psfact.select]
      psfact.coeffs.sel <- psfact.coeffs[psfact.select]
      if(max(psfact.coeffs > 1)){
        psfact.select2 <- psfact.coeffs.sel > 1
        psfact.pwr <- paste("^",psfact.coeffs.sel[psfact.select2],sep="")
        toprint[psfact.select2] <- paste(toprint[psfact.select2],
                                         psfact.pwr,sep="")
      }
      cat(toprint,sep=":")
    }
    cat("\n\n")
  }
  else cat("nil\n\n")

  cat("\n")
  return(invisible(alias.stats))
} ## end alias.keymatrix

##---------------------------------------------------------------------------
# Method alias for "keymatrix"
 setMethod("alias", signature(object="keymatrix"),
          definition=alias.keymatrix)


##--------------------------------------------------------------------------
# "show.keymatrix" 
# Print an object of class keymatrix
# ARGUMENT
#   object: an object of class  keymatrix
# RETURN
#   an invisible NULL.
# NOTE
# The number of rows and columns of the matrices that are printed
# are limited by the option planor.max.print
# EXAMPLES
# K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#   nunits=3^3, base=~A+B+C, max.sol=2)
# show(K0[[1]][[1]])
# ---------------------------------------------

show.keymatrix <- function(object){
  cat("An object of class keymatrix\n")
    p.k <- object@p
    printpower(p.k)
    printgmat(object)
  invisible()
} ## end show.keymatrix
# --------------------------------------
# "show" method for "keymatrix"
# --------------------------------------
setMethod("show", signature(object="keymatrix"),
          definition=show.keymatrix)

