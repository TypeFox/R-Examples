`numPerms` <- function(object, control = how()) {
  ## constant holding types where something is permuted
  TYPES <- c("free","grid","series","none")

  ## expand object if a numeric or integer vector of length 1
  if((is.numeric(object) || is.integer(object)) &&
     (length(object) == 1))
    object <- seq_len(object)
  ## number of observations in data
  n <- nobs(object)

  ## get the permutation levels from control
  WI <- getWithin(control)
  PL <- getPlots(control)
  BL <- getBlocks(control)

  ## any strata to permute within / blocking?
  BLOCKS <- getStrata(control, which = "blocks")
  PSTRATA <- getStrata(control, which = "plots")
  typeP <- getType(control, which = "plots")
  typeW <- getType(control, which = "within")

  ## mirroring?
  mirrorP <- getMirror(control, which = "plots")
  mirrorW <- getMirror(control, which = "within")

  ## constant - i.e. same perm within each plot?
  constantW <- getConstant(control)

  ## grid dimensions
  colW <- getCol(control, which = "within")
  colP <- getRow(control, which = "plots")

  ## Some checks; i) Plot strata must be of same size when permuting strata
  ##                 or having the same constant permutation within strata
  ##             ii) In grid designs, grids must be of the same size for all
  ##                 strata
  ##
  ## FIXME - this probably should be in check()!
  if(!is.null(PSTRATA)) {
    tab <- table(PSTRATA)
    same.n <- length(unique(tab))
    if((typeP != "none" || isTRUE(constantW)) && same.n > 1) {
      stop("All levels of strata must have same number of samples for chosen scheme")
    }
    if(typeP == "grid" && same.n > 1) {
      stop("Unbalanced grid designs are not supported")
    }
  }

  ## the various designs allowed imply multipliers to number of samples
  ## for the restricted permutations

  mult.p <- mult.wi <- 1

  ## within types
  if(typeW %in% c("series","grid")) {
    mult.wi <- 2
    if(isTRUE(all.equal(typeW, "grid")) && !is.null(colW) && colW > 2) {
      mult.wi <- 4
    } else {
      if(isTRUE(all.equal(n, 2)))
        mult.wi <- 1
    }
  }
  ## plot-level types
  if(typeP %in% c("series","grid")) {
    mult.p <- 2
    if(isTRUE(all.equal(typeP, "grid")) && !is.null(colP) && colP > 2) {
      mult.p <- 4
    } else {
      if(isTRUE(all.equal(length(tab), 2))) # was all.equal(n, 2)
        mult.p <- 1
    }
  }

  ## within
  ## another check - shouldn't this be moved? FIXME
  if(!typeW %in% TYPES) {
    stop("Ambiguous permutation type in 'control$within$type'")
  }

  ## calculate the number of possible permutations

  ## Compute number of permutations for each block
  if(is.null(BLOCKS))
      BLOCKS <- factor(rep(1, n))

  ## split an index vector
  indv <- seq_len(n)
  spl <- split(indv, BLOCKS)

  ## loop over the components of spl & apply doNumPerms
  np <- sapply(spl, doNumPerms, mult.p, mult.wi, typeP, typeW, PSTRATA,
               mirrorP, mirrorW, constantW)

  ## multiply up n perms per block
  prod(np)
}

`doNumPerms` <- function(obs, mult.p, mult.wi, typeP, typeW, PSTRATA,
                         mirrorP, mirrorW, constantW) {
    n <- nobs(obs) ## obs is index vector for object, split by blocks

    if(!is.null(PSTRATA)) {
        ## take only the PSTRATA needed for this block, drop unused levels
        PSTRATA <- droplevels(PSTRATA[obs])

        ## need only those strata for the current block. As obs is the index
        ## vector, split by block, this now gives nobs per plot strata
        tab <- table(PSTRATA)
        same.n <- length(unitab <- unique(tab))
    }

    ## plots
    num.p <- if(isTRUE(all.equal(typeP, "free"))) {
        exp(lfactorial(length(levels(PSTRATA))))
    } else if(typeP %in% c("series", "grid")) {
        if(isTRUE(mirrorP)) {
            mult.p * length(tab)
        } else {
            length(tab)
        }
    } else {
        1
    }

    num.wi <- if(isTRUE(all.equal(typeW, "none"))) {
        ## no within permutations. note we multiply num.p by this
        ## values so it is 1 not 0!!
        1
    } else if(isTRUE(all.equal(typeW, "free"))) {
        if(!is.null(PSTRATA)) {
            if(constantW) {
                factorial(tab[1])
            } else {
                prod(factorial(tab))
            }
        } else {
            exp(lfactorial(n))
        }
    } else {
        if(!is.null(PSTRATA)) {
            if(same.n > 1) {
                multi <- rep(2, length = length(tab))
                multi[which(tab == 2)] <- 1
                if(mirrorW) {
                    prod(multi * tab)
                } else {
                    prod(tab)
                }
            } else {
                if(mirrorW) {
                    if(constantW) {
                        mult.wi * unitab[1]
                    } else {
                        prod(mult.wi * tab)
                    }
                } else {
                    if(constantW) {
                        unitab[1] ## FIXME: unitab[1]?? (unique(tab)[1])
                    } else {
                        prod(tab)
                    }
                }
            }
        } else {
            if(mirrorW)
                mult.wi * n
            else
                n
        }
    }

    ## return
    num.p * num.wi
}
