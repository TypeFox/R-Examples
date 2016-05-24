### compact letter displays

cld <- function(object, ...)
    UseMethod("cld")

cld.glht <- function(object, level = 0.05, decreasing = FALSE, ...)
    cld(summary(object), level = level, decreasing = decreasing)

extr <- function(object) {

    stopifnot(object$type == "Tukey")
    

    mf <- model.frame(object$model)
    if (!is.null(attr(mf, "terms"))) {
        tm <- attr(mf, "terms")
    } else {
        tm <- try(terms(object$model))
        if (inherits(tm, "try-error")) stop("no terms component found")
    }
    

    ### <FIXME> not very nice    
    if(any(class(object$model) == "lme")){
      mf <- get_all_vars(tm, mf)
    }
    ### </FIXME>
    
    covar <- (length(attr(tm, "term.labels")) > 1)
    y <- mf[[1L]]
    yname <- colnames(mf)[[1L]]

    stopifnot(length(object$focus) == 1)
    x <- mf[[object$focus]]
    xname <- object$focus

    K <- contrMat(table(x), type = "Tukey")
    comps <- cbind(apply(K, 1, function(k) levels(x)[k == 1]),
                   apply(K, 1, function(k) levels(x)[k == -1]))

    f <- if (inherits(object$model, "coxph")) predict else fitted
    lp <- f(object$model)

    ret <- list(y = y, yname = yname,  
                x = x, xname = xname, 
                weights = model.weights(mf), 
                lp = lp, covar = covar, comps = comps)
    return(ret)
}

cld.summary.glht <- function(object, level = 0.05, decreasing = FALSE, ...) {

    ret <- extr(object)
    signif <- (object$test$pvalues < level)
    # Order the levels according to its mean
    # Tidy up: ret$y[1:length(ret$x)]], cox models concatenates a vector of live/dead
    # I think this way is easier than to deal with gsub later and it's more general
    lvl_order <- levels(ret$x)[order(tapply(as.numeric(ret$y)[1:length(ret$x)], ret$x, mean))]
    # names(signif) <- gsub("\\s", "", rownames(object$linfct))
    ret$signif <- signif
    ret$mcletters <- insert_absorb(signif, decreasing = decreasing, 
                                   comps = ret$comps, lvl_order = lvl_order, ...)
                           
   # start edit
    
    ret$mcletters$Letters <- ret$mcletters$Letters[levels(ret$x)]
    ret$mcletters$monospacedLetters <- ret$mcletters$monospacedLetters[levels(ret$x)]
    ret$mcletters$LetterMatrix <- ret$mcletters$LetterMatrix[levels(ret$x),]
    
   # end edit
                           
                          
    class(ret) <- "cld"
    ret
}

cld.confint.glht <- function(object, decreasing = FALSE, ...) {

    ret <- extr(object)
    ### significant, if confidence interval does not contains 0
    signif <- !(object$confint[, "lwr"] < 0 & object$confint[, "upr"] > 0)
    # Tidy up: ret$y[1:length(ret$x)]], cox models concatenates a vector of live/dead
    # I think this way is easier than to deal with gsub later and it's more general
    lvl_order <- levels(ret$x)[order(tapply(as.numeric(ret$y)[1:length(ret$x)], ret$x, mean))]
    # names(signif) <- gsub("\\s", "", rownames(object$linfct))
    ret$signif <- signif
    ret$mcletters <- insert_absorb(signif, decreasing = decreasing, 
                                   comps = ret$comps, lvl_order = lvl_order, ...)
                           
    # start edit
                           
    ret$mcletters$Letters <- ret$mcletters$Letters[levels(ret$x)]
    ret$mcletters$monospacedLetters <- ret$mcletters$monospacedLetters[levels(ret$x)]
    ret$mcletters$LetterMatrix <- ret$mcletters$LetterMatrix[levels(ret$x),]
                                   
    # end edit                      
                           
    class(ret) <- "cld" 
    ret
}

print.cld <- function(x, ...)
    print(x$mcletters$Letters)

plot.cld <- function(x, type = c("response", "lp"), ...) {
   
    mcletters <- x$mcletters
    ### ms = mono-spaced
    msletters <- mcletters$monospacedLetters
    ### v = vertical
    vletters <- sapply(msletters,
                       function(x) paste(strsplit(x, "")[[1]], "\n", collapse = ""))
    vletters <- vletters[gsub(" ", "", levels(x$x))]
    msletters <- msletters[gsub(" ", "", levels(x$x))]
    type <- match.arg(type)
    dat <- x[c("x", "y", "lp")]
    if (is.null(x$weights)) {
        dat$weights <- rep(1, NROW(x$y))
    } else {
        dat$weights <- x$weights
    }
    dat <- as.data.frame(dat)
    xn <- x$xname
    yn <- x$yname
    if (!is.null(list(...)$xlab)) xn <- list(...)$xlab
    if (!is.null(list(...)$ylab)) yn <- list(...)$ylab

    if (x$covar || type == "lp") {
        ### boxplot to make use of "..." argument
        boxplot(lp ~ x, data = dat, xlab = xn, ylab = "linear predictor", ...)
        axis(3, at = 1:nlevels(dat$x), labels = vletters)
    } else {
        if (is.integer(dat$y)) dat$y <- as.numeric(dat$y)
        switch(class(dat$y), 
            "numeric" = {
                ### boxplot to make use of "..." argument
                boxplot(y ~ x, data = dat, xlab = xn, ylab = yn, ...)
                axis(3, at = 1:nlevels(dat$x), labels = vletters)
            },
            "factor" = {
                at <- xtabs(weights ~ x, data = dat) / sum(dat$weights)
                at <- cumsum(at) - at / 2
                mosaicplot(xtabs(weights ~ x + y, data = dat), main = NULL,
                           xlab = xn, ylab = yn, ...)
                axis(3, at = at, labels = vletters, tick = FALSE)
            },
            "Surv" = {
                plot(survfit(y ~ x, data = dat), lty = 1:nlevels(dat$x), ...)
                nc <- nchar(levels(dat$x))                              
                spaces <- unlist(lapply( max(nc)-nc, function(x) return(paste( rep(" ",x) ,collapse="")))) 
#                old.par <- par(family="mono") 
                legend("topright", lty = 1:nlevels(dat$x), 
                       legend = paste(levels(dat$x), spaces, ": ", msletters, sep=""), 
                       ...)
#                par(old.par)
            })
    }
}

# Function implements the insert-absorb (sweep) heuristic of Piepho 2004:
# "An Algorithm for a Letter-Based Representation of All-Pairwise Comparisons"
#
# x         ... vector of logicals indicating significant comparisons with hyphenated
#               names e.g. A-B, treatmentA-treatmentB, ...
# Letters   ... a set of user defined letters { default is Letters=c(letters, LETTERS) }
# separator ... a separating character used to produce a sufficiently large set of
#               characters for a compact letter display (default is separator=".") in case
#               the number of letters required exceeds the number of letters available
# Decreasing ... Inverse the order of the letters 

insert_absorb <- function( x, Letters=c(letters, LETTERS), separator=".", decreasing = FALSE, 
                           comps = NULL, lvl_order){

  obj_x <- deparse(substitute(x))
  if (is.null(comps)) {
      namx <- names(x)
      namx <- gsub(" ", "", names(x))
      if(length(namx) != length(x))
          stop("Names required for ", obj_x)
      split_names <- strsplit(namx, "-")
      stopifnot( sapply(split_names, length) == 2 )
      comps <- t(as.matrix(as.data.frame(split_names)))
  } 
  rownames(comps) <- names(x)
  lvls <- lvl_order
  n <- length(lvls)
  lmat <- array(TRUE, dim=c(n,1), dimnames=list(lvls, NULL) )

  if( sum(x) == 0 ){                                                        # no differences
    ltrs <- rep(get_letters(1, Letters=Letters, separator=separator), length(lvls) )
    names(ltrs) <- lvls
    colnames(lmat) <- ltrs[1]
    msl <- ltrs
    ret <- list(Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat)
    class(ret) <- "multcompLetters"
    return(ret)
  }
  else{
    signifs <- comps[x,,drop=FALSE]
    
    absorb <- function(m){
      for(j in 1:(ncol(m)-1)){
        for(k in (j+1):ncol(m)){
          if( all(m[which(m[,k]),k] & m[which(m[,k]),j]) ){                 # column k fully contained in column j
            m <- m[,-k, drop=FALSE]
            return(absorb(m))
          }
          else if( all(m[which(m[,j]),k] & m[which(m[,j]),j]) ){            # column j fully contained in column k
            m <- m[,-j, drop=FALSE]
            return(absorb(m))
          }
        }
      }
      return(m)
    }
    for( i in 1:nrow(signifs) ){                                            # insert
      tmpcomp <- signifs[i,]
      wassert <- which(lmat[tmpcomp[1],] & lmat[tmpcomp[2],])               # which columns wrongly assert nonsignificance
      if(any(wassert)){
        tmpcols <- lmat[,wassert,drop=FALSE]
        tmpcols[tmpcomp[2],] <- FALSE
        lmat[tmpcomp[1],wassert] <- FALSE
        lmat <- cbind(lmat, tmpcols)
        colnames(lmat) <- get_letters( ncol(lmat), Letters=Letters,
                                       separator=separator)
        if(ncol(lmat) > 1){                                                 # absorb columns if possible
          lmat <- absorb(lmat)
          colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                         separator=separator )
        }
      }
    }
  }
  lmat <- lmat[,order(apply(lmat, 2, sum))]
  lmat <- sweepLetters(lmat)                                                                  # 1st
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  lmat <- lmat[,order(apply(lmat, 2, sum))]                                                   # 2nd sweep
  lmat <- sweepLetters(lmat)
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x)))), 
                           decreasing = decreasing))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  ltrs <- apply(lmat,1,function(x) return(paste(names(x)[which(x)], sep="", collapse="") ) )
  msl <- matrix(ncol=ncol(lmat), nrow=nrow(lmat))                                             # prepare monospaced letters
  for( i in 1:nrow(lmat) ){
    msl[i,which(lmat[i,])] <- colnames(lmat)[which(lmat[i,])]
    absent <- which(!lmat[i,])
    if( length(absent) < 2 ){
      if( length(absent) == 0 )
        next
      else{
        msl[i,absent] <- paste( rep(" ", nchar(colnames(lmat)[absent])), collapse="" )
      }
    }
    else{
      msl[i,absent] <- unlist( lapply( sapply( nchar(colnames(lmat)[absent]),
                                               function(x) return(rep( " ",x)) ),
                                       paste, collapse="") )
    }
  }
  msl <- apply(msl, 1, paste, collapse="")
  names(msl) <- rownames(lmat)
  ret <- list( Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat, 
               aLetters = Letters, aseparator = separator )
  class(ret) <- "multcompLetters"
  return(ret)
}


# All redundant letters are swept out without altering the information within a LetterMatrix.
#
# mat         ... a LetterMatrix as produced by function insert_absorb()
# start.col   ... either a single integer specifying the column to start with or a vector
#                 of max. length equal to ncol(mat) specifying the column order to be used.
# Letters     ... a set of user defined letters { default is Letters=c(letters, LETTERS) }
# separator   ... a separating character used to produce a sufficiently large set of
#                 characters for a compact letter display (default is separator=".") in case
#                 the number of letters required exceeds the number of letters available 

sweepLetters <- function(mat, start.col=1, Letters=c(letters, LETTERS), separator="."){

  stopifnot( all(start.col %in% 1:ncol(mat)) )
  locked <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))          # 1 indicates that another letter dependes on this entry
  cols <- 1:ncol(mat)
  cols <- cols[c( start.col, cols[-start.col] )]
  if( any(is.na(cols) ) )
    cols <- cols[-which(is.na(cols))]

  for( i in cols){
    tmp <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))
    tmp[which(mat[,i]),] <- mat[which(mat[,i]),]                        # get items of those rows which are TRUE in col "i"
    one <- which(tmp[,i]==1)

    if( all(apply(tmp[,-i,drop=FALSE], 1, function(x) return( any(x==1) ))) ){     # there is at least one row "l" where mat[l,i] is the only item which is TRUE i.e. no item can be removed in this column
      next
    }
    for( j in one ){                                                    # over all 1's
      if( locked[j,i] == 1 ){                                           # item is locked
        next
      }
      chck <- 0
      lck <- list()
      for( k in one ){
        if( j==k ){
          next
        }
        else{                                                           # pair j-k
          rows <- tmp[c(j,k),]
          dbl <- rows[1,] & rows[2,]
          hit <- which(dbl)
          hit <- hit[-which(hit==i)]
          dbl <- rows[1,-i,drop=FALSE] & rows[2,-i,drop=FALSE]
          if( any(dbl) ){
            chck <- chck + 1
            lck[[chck]] <- list(c(j,hit[length(hit)]), c(k,hit[length(hit)]))      # record items which have to be locked, use last column if multiple hits
          }
        }
      }
      if( (chck == (length(one)-1)) && chck != 0 ){                     # item is redundant
        for( k in 1:length(lck) ){                                      # lock items
          locked[ lck[[k]][[1]][1], lck[[k]][[1]][2] ] <- 1
          locked[ lck[[k]][[2]][1], lck[[k]][[2]][2] ] <- 1
        }
        mat[j,i] <- FALSE                                               # delete redundant entry
      }
    }
    if(all(mat[,i]==FALSE)){                                           # delete column where each entry is FALSE and restart
      mat <- mat[,-i,drop=FALSE]
      colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
      return(sweepLetters(mat, Letters=Letters, separator=separator))
    }
  }
  onlyF <- apply(mat, 2, function(x) return(all(!x)))
  if( any(onlyF) ){                                                     # There are columns with just FALSE entries
    mat <- mat[,-which(onlyF),drop=FALSE]
    colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
  }
  return( mat )
}

# Create a set of letters for a letter display. If "n" exceeds the number of letters
# specified in "Letters", they are recycled with one or more separating character(s)
# preceding each recycled letter.
# e.g. get_letters(10, Letters=letters[1:4]) produces:  "a"   "b"   "c"   "d"   ".a"  ".b"  ".c"  ".d"  "..a" "..b"
#
# n             ... number of letters
# Letters       ... the set of characters to be used
# separator     ... a character to be used as separator e.g.
#                   n=5, Letters=c("a","b") => "a", "b", ".a", ".b", "..a"

get_letters <- function( n, Letters=c(letters, LETTERS), separator="." ){

  n.complete <- floor(n / length(Letters))        # number of complete sets of Letters
  n.partial <- n %% length(Letters)               # number of additional Letters
  lett <- character()
  separ=""
  if( n.complete > 0 ){
    for( i in 1:n.complete ){
      lett <- c(lett, paste(separ, Letters, sep="") )
      separ <- paste( separ, separator, sep="" )
    }
  }
  if(n.partial > 0 )
    lett <- c(lett, paste(separ, Letters[1:n.partial], sep="") )
  return(lett)
}
