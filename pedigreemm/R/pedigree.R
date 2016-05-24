#### "pedigree" class methods

#' Constructor for pedigree objects
#'
#' A simple constructor for a pedigree object.  The main point for the
#' constructor is to use coercions to make the calls easier.
#' @param sire integer vector or factor representation of the sires
#' @param dam integer vector or factor representation of the dams
#' @param label character vector of labels
#' @return an pedigree object of class \linkS4class{pedigree}

#' @note \code{sire}, \code{dam} and \code{label} must all have the
#' same length and all labels in \code{sire} and \code{dam} must occur
#' in \code{label}
#' @export
pedigree <- function(sire, dam, label) {
    n <- length(sire)
    labelex <- c(label, NA, 0)
    stopifnot(n == length(dam),
              n == length(label),
              all(sire %in% labelex),
              all(dam %in% labelex))
    sire <- as.integer(factor(sire, levels = label))
    dam <- as.integer(factor(dam, levels = label))
    sire[sire < 1 | sire > n] <- NA
    dam[dam < 1 | dam > n] <- NA
    new("pedigree", sire = sire, dam = dam,
        label = as.character(label))
}

#' Coerce a pedigree to a sparse triangular matrix
#'
#' Create a sparse, unit lower triangular matrix from a pedigree.  The
#' matrix the L factor in the LDL' Cholesky factorization of the
#' inverse of the relationship matrix.
#'
#' @export
setAs("pedigree", "sparseMatrix", # representation as T^{-1}
      function(from) {
	  sire <- from@sire
	  n <- length(sire)
	  animal <- seq_along(sire)
	  j <- c(sire, from@dam)
	  ind <- !is.na(j)
	  as(new("dtTMatrix", i = rep.int(animal, 2)[ind] - 1L,
		 j = j[ind] - 1L, x = rep.int(-0.5, sum(ind)),
		 Dim = c(n,n), Dimnames = list(from@label, NULL),
		 uplo = "L", diag = "U"), "dtCMatrix")
      })

## these data frames are now storage efficient but print less nicely
setAs("pedigree", "data.frame",
      function(from)
      data.frame(sire = from@sire, dam = from@dam,
		 row.names = from@label))

#' Convert a pedigree to a data frame
#'
#' Express a pedigree as a data frame with \code{sire} and
#' \code{dam} stored as factors.  If the pedigree is an object of
#' class \linkS4class{pedinbred} then the inbreeding coefficients are
#' appended as the variable \code{F}
#'
#' @param x a pedigree object of class \linkS4class{pedigree}
#' @return a data frame
#'
ped2DF <- function(x) {
    stopifnot(is(x, "pedigree"))
    lab <- x@label
    lev <- seq_along(lab)
    ans <- data.frame(sire = factor(x@sire, levels = lev, labels = lab),
                      dam  = factor(x@dam,  levels = lev, labels = lab),
                      row.names = lab)
    if (is(x, "pedinbred")) ans <- cbind(ans, F = x@F)
    ans
}

setMethod("show", signature(object = "pedigree"),
	  function(object) print(ped2DF(object)))

setMethod("head", "pedigree", function(x, ...)
	  do.call("head", list(x = ped2DF(x), ...)))

setMethod("tail", "pedigree", function(x, ...)
	  do.call("tail", list(x = ped2DF(x), ...)))

setMethod("chol", "pedigree",
          function(x, pivot, LINPACK) {
              ttrans <- solve(t(as(x, "dtCMatrix")))
              .Call(pedigree_chol, x,
                    as(.Call("Csparse_diagU2N", t(ttrans), PACKAGE = "Matrix"),
                       "dtCMatrix"))
          })

#' Inbreeding coefficients from a pedigree
#'
#' Create the inbreeding coefficients according to the algorithm given
#' in "Comparison of four direct algorithms for computing inbreeding
#' coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
#' Science Journal (2005) 76, 401--406.
#'
#' @param ped an object that inherits from class \linkS4class{pedigree}
#' @return the inbreeding coefficients as a numeric vector
#' @export
inbreeding <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    .Call(pedigree_inbreeding, ped)
}

#' Diagonal of D in the A = TDT' factorization.
#'
#' Determine the diagonal factor in the decomposition of the
#' relationship matrix A as TDT' where T is unit lower triangular.
#'
#' @param ped an object that inherits from class \linkS4class{pedigree}
#' @return a numeric vector
#' @export
Dmat <- function(ped)
{
    F <- inbreeding(ped)
    sire <- ped@sire
    dam <- ped@dam
    Fsire <- ifelse(is.na(sire), -1, F[sire])
    Fdam <-  ifelse(is.na(dam), -1, F[dam])
    ans <- 1 - 0.25 * (2 + Fsire + Fdam)
    names(ans) <- ped@label
    ans
}

#' Relationship factor from a pedigree
#'
#' Determine the right Cholesky factor of the relationship matrix for
#' the pedigree \code{ped}, possibly restricted to the specific labels
#' that occur in \code{labs}. 
#' 
#' @param ped a pedigree that includes the individuals who occur in svec
#' @param labs a character vector or a factor giving the labels to
#'    which to restrict the relationship matrix. If \code{labs} is a
#'    factor then the levels of the factor are used as the labels.
#'    Default is the complete set of labels in the pedigree.
#' @return an object that inherits from \linkS4class{CHMfactor}
#' @export
relfactor <- function(ped, labs)
{
    stopifnot(is(ped, "pedigree"))
    if (missing(labs))                  # square case
        return(Diagonal(x = sqrt(Dmat(ped))) %*% 
               solve(t(as(ped, "sparseMatrix"))))
    labs <- factor(labs) # drop unused levels from a factor
    stopifnot(all(labs %in% ped@label))
#    rect <- Diagonal(x = sqrt(Dmat(ped))) %*% 
#        solve(t(as(ped, "sparseMatrix")), # rectangular factor
#              as(factor(labs, levels = ped@label),"sparseMatrix"))

    rect <- Diagonal(x = sqrt(Dmat(ped))) %*%
        solve(t(as(ped, "sparseMatrix")), # rectangular factor
              as(factor(ped@label, levels = ped@label),"sparseMatrix"))
    tmpA<-crossprod(rect)
    tmp<- ped@label %in% labs 
    tmpA<-tmpA[tmp,tmp]

    orlab <- order(as.numeric(factor(labped<-ped@label[tmp], levels=labs, ordered=T)))
    labped<- as.character(labped[orlab])
    tmpA  <- tmpA[orlab, orlab]
    stopifnot(all.equal(as.character(labped), as.character(labs)))
    relf<-chol(tmpA)
    dimnames(relf)[[1]]<- dimnames(relf)[[2]]<-labs
    relf
}



#' Inverse of the Relationship Matrix
#'
#' @param ped a pedigree that includes the individuals who occur in svec
#'    which to restrict the relationship matrix. If \code{labs} is a
#'    factor then the levels of the factor are used as the labels.
#'    Default is the complete set of labels in the pedigree.
#' @return an object that inherits from \linkS4class{CHMfactor}
#' @export
getAInv <- function(ped)
{
    stopifnot(is(ped, "pedigree"))    
    T_Inv <- as(ped, "sparseMatrix")
    D_Inv <- diag(1/Dmat(ped))
    aiMx<-t(T_Inv) %*% D_Inv %*% T_Inv
    dimnames(aiMx)[[1]]<-dimnames(aiMx)[[2]] <-ped@label
    aiMx
}

#' Additive Relationship Matrix
#'
#' @param ped a pedigree that includes the individuals who occur in svec
#'    which to restrict the relationship matrix. If \code{labs} is a
#'    factor then the levels of the factor are used as the labels.
#'    Default is the complete set of labels in the pedigree.
#' @return an object that inherits from \linkS4class{CHMfactor}
#' @export
getA <- function(ped)
{
    stopifnot(is(ped, "pedigree"))
    aMx<-crossprod(relfactor(ped))
    dimnames(aMx)[[1]]<-dimnames(aMx)[[2]] <-ped@label
    aMx
}

#' Counts number of generations of ancestors for one subject. Use recursion. 
#' 
#' @param pede data frame with a pedigree and a column for the number 
#' of generations of each subject.
#' @param id subject for which we want the number of generations.
#' @return a data frame object with the pedigree and generation of
#'    ancestors for subject id.
#' 
getGenAncestors <- function(pede, id){
    j <- which(pede$id==id)
    parents <- c(pede$sire[j], pede$dam[j])
    parents <- parents[!is.na(parents)]
    np <- length(parents)
    if(np==0)
    {
        pede$gene[j] <-0
        return(pede)
    }
    ## get the number of generations in parent number one
    tmpgenP1 <- pede$gene[pede$id==parents[1]]
    if( is.na(tmpgenP1))
    {
        pede <- getGenAncestors(pede, parents[1])
        genP1  <- 1 + pede$gene[pede$id==parents[1]]
    }  else {
        genP1 <- 1 + tmpgenP1
    }
    ## find out if there is a parent number two
    if (np==2){
        tmpgenP2 <- pede$gene[pede$id==parents[2]]
        if( is.na(tmpgenP2))
        {
            pede <- getGenAncestors(pede, parents[2])
            genP2  <- 1 + pede$gene[pede$id==parents[2]]
        }  else {
            genP2 <- 1 + tmpgenP2
        }
        genP1 <- max(genP1, genP2)
    }
    pede$gene[j] <- genP1
    ## print(paste('id:', id, ', gen:', genP1, ', row:', j))
    pede
}

#' Edits a disordered or incomplete pedigree, 
#'    1_ add labels for the sires and dams not listed as labels before.
#'    2_ order pedigree based on recursive calls to getGenAncestors.
#' @param sire integer vector or factor representation of the sires
#' @param dam integer vector or factor representation of the dams
#' @param label character vector of labels
#' @param verbose logical to print the row of the pedigree that the 
#'    function is ordering. Default is FALSE.
#' @return a data frame with the pedigree ordered.
#' @export
editPed <- function(sire, dam, label, verbose = FALSE)
{
    nped <- length(sire)
    if (nped != length(dam))  stop("sire and dam have to be of the same length")
    if (nped != length(label)) stop("label has to be of the same length than sire and dam")
    tmp <- unique(sort(c(as.character(sire), as.character(dam))))
    
    missingP <-NULL
    if(any(completeId <- ! tmp %in% as.character(label))) missingP <- tmp[completeId]
    labelOl <- c(as.character(missingP),as.character(label))
    sireOl <- c(rep(NA, times=length(missingP)),as.character(sire))
    damOl  <- c(rep(NA, times=length(missingP)),as.character(dam))
    sire <- as.integer(factor(sireOl, levels = labelOl))
    dam  <- as.integer(factor(damOl, levels = labelOl))
    nped <-length(labelOl)
    label <-1:nped
    sire[!is.na(sire) & (sire<1 | sire>nped)] <- NA
    dam[!is.na(dam) & (dam < 1 | dam > nped)] <- NA
    pede <- data.frame(id= label, sire= sire, dam= dam, gene=rep(NA, times=nped))
    noParents <- (is.na(pede$sire) & is.na(pede$dam))
    pede$gene[noParents] <- 0
    for(i in 1:nped){
        if(verbose) print(i)
        if(is.na(pede$gene[i])){
            id <-pede$id[i]
            pede <-getGenAncestors(pede, id)
        }}
    ord<- order(pede$gene)
    ans<-data.frame(label=labelOl, sire=sireOl, dam=damOl, gene=pede$gene,
                    stringsAsFactors =F)
    ans[ord,]
}

pedigreemm <-
    function(formula, data, family = NULL, REML = TRUE, pedigree = list(),
             control = list(), start = NULL, verbose = FALSE, 
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    gaus <- FALSE
    if (is.null(family)) {
        gaus <- TRUE
    } else {
        ## copied from glm()
        if (is.character(family)) 
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family)) 
            family <- family()
        if (!inherits(family, "family")) stop("unknown family type")
        gaus <- family$family == "gaussian" && family$link == "identity"
    }
    mc <- match.call()
    lmerc <- mc                         # create a call to lmer
    lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
    lmerc$pedigree <- NULL
    if (!gaus) lmerc$REML <- NULL

    if (!length(pedigree))              # call [g]lmer instead
        return(eval.parent(lmerc))

    stopifnot(is.list(pedigree),        # check the pedigree argument
              length(names(pedigree)) == length(pedigree),
              all(sapply(pedigree, is, class2 = "pedigree")))
    
    lmf <- eval(lmerc, parent.frame())

    
    relfac <- pedigree          # copy the pedigree list for relfactor
    pnms <- names(pedigree)
    pp <- lmf@pp
    resp <- lmf@resp
    fl <- lmf@flist
    stopifnot(all(pnms %in% names(fl)))
    asgn <- attr(fl, "assign")
    Zt <- pp$Zt
    for (i in seq_along(pedigree)) {
        tn <- which(match(pnms[i], names(fl)) == asgn)
        if (length(tn) > 1)
            stop("a pedigree factor must be associated with only one r.e. term")
        ind <- (lmf@Gp)[tn:(tn+1L)]
        rowsi <- (ind[1]+1L):ind[2]
        relfac[[i]] <- relfactor(pedigree[[i]], rownames(Zt)[rowsi])
        Zt[rowsi,] <- relfac[[i]] %*% Zt[rowsi,]
    }
    reTrms <- list(Zt=Zt,theta=lmf@theta,Lambdat=pp$Lambdat,Lind=pp$Lind,
                   lower=lmf@lower,flist=lmf@flist,cnms=lmf@cnms, Gp=lmf@Gp)
    dfl <- list(fr=lmf@frame, X=pp$X, reTrms=reTrms, start=lmf@theta)
    if (gaus) {
        dfl$REML = resp$REML > 0L
        devfun <- do.call(mkLmerDevfun,dfl)
        opt <- optimizeLmer(devfun, optimizer="Nelder_Mead",...)
    } else {
        dfl$family <- family
        devfun <- do.call(mkGlmerDevfun,dfl)
        opt <- optimizeGlmer(devfun, optimizer="Nelder_Mead",...)
    }
    mm <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)
    cls <- if (gaus) "lmerpedigreemm" else "glmerpedigreemm"
    ans <- do.call(new, list(Class=cls, relfac=relfac,
                             frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                             theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                             devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
    ans@call <- evalq(mc)
    ans
}

setMethod("ranef", signature(object = "pedigreemm"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), pedigree = TRUE, ...)
      {
          if ((postVar <- as.logical(postVar)) && (pedigree <- as.logical(pedigree)))
              stop("code for applying pedigree and posterior variances not yet written")
          ans <- ranef(as(object, "merMod"), postVar, drop = FALSE)
          ans <- ans[whichel]
          if (pedigree) {
              if (postVar)
                  stop("postVar and pedigree cannot both be true")
              rf <- object@relfac
              for (nm in names(rf)) {
                  dm <- data.matrix(ans[[nm]])
                  cn <- colnames(dm)
                  rn <- rownames(dm)
                  dm <- as.matrix(rf[[nm]] %*% dm)
                  colnames(dm) <- cn
                  rownames(dm) <- rn
                  ans[[nm]] <- data.frame(dm, check.names = FALSE)
              }
          }
          if (drop)
              ans <- lapply(ans, function(el)
                        {
                            if (ncol(el) > 1) return(el)
                            pv <- drop(attr(el, "postVar"))
                            el <- drop(as.matrix(el))
                            if (!is.null(pv))
                                attr(el, "postVar") <- pv
                            el
                        })
          ans
      })


setMethod("fitted", signature(object = "pedigreemm"),
          function(object, ...) {
              stop("fitted() not applicable to pedigreemm objects")
          })


setMethod("residuals", signature(object = "pedigreemm"),
          function(object, ...) {
              stop("residuals() not applicable to pedigreemm objects")
          })

