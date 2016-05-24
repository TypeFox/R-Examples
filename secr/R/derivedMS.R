############################################################################################
## package 'secr'
## derived.R
## 2010-10-22 allow object to be secrlist
## 2011-03-27 adjustments for zero capthist
## 2012-11-03 CLdensity and CLgradient moved from functions.R
## 2013-06-24 fixed bug in esa dummy grp (0 should be 1) that caused intermittent fault in derived
## 2014-04-05 fixed bug mapply SIMPLIFY = FALSE
## 2014-10-28 adjust for linear mask
############################################################################################


CLdensity <- function (beta, object, individuals, sessnum)
## object is a fitted secr object (CL=T)
## individuals is vector indexing the subset of a to be used
# Return the density for given g0, sigma, z in beta
# Only 1 session
{
    sum(1 / esa (object, sessnum, beta)[individuals])
}
############################################################################################

CLgradient <- function (object, individuals, sessnum, eps=0.001, clust)
## object is a fitted secr object (CL=T)
## individuals is vector indexing the subset of a to be used
{
  beta <- object$fit$par
  if (object$detectfn %in% c(0,2,9)) {    ## halfnormal, exponential and binary SS have just 2 parameters
      est <- beta[1:2]
      g   <- double(2)
  } else {                       ## other detectfn have 3 parameters
      est <- beta[1:3]
      g   <- double(3)
  }

  est <- beta
  g   <- beta

  ## consider replacing this with packaged, optimized gradient function (nlme fnHess?)

  grad <- function(i, est, eps) {
          temp     <- est[i]
          if (temp != 0.0) delta <- eps * abs(temp)
          else             delta <- eps
          est[i]  <- temp - delta
          fminus  <- CLdensity (est, object, individuals, sessnum)
          est[i]  <- temp + delta
          fplus   <- CLdensity (est, object, individuals, sessnum)
          (fplus - fminus) / (2.0 * delta)
  }
  if (is.null(clust)) {
      sapply(1:length(est), grad, est = est, eps = eps)
  }
  else {
      clusterExport(clust, c("esa", "CLdensity","object","individuals","sessnum"), environment())
      parSapply(clust, 1:length(est), grad, est = est, eps = eps)
  }
}
############################################################################################

CLmeanesa <- function (beta, object, individuals, sessnum, noccasions = NULL)
## object is a fitted secr object (CL=T)
## individuals is vector indexing the subset of a to be used
# Return the weighted mean esa for given g0, sigma, z in beta
# Only 1 session
{
## mean (esa (object, sessnum, beta)[individuals])
## modified 2010-11-30 after suggestion of DLB

##  noccasions = NULL added 2011-04-04

    a <- esa (object, sessnum, beta, noccasions=noccasions)[individuals]
    length(a) / sum (1/a)
}
############################################################################################

esagradient <- function (object, individuals, sessnum, noccasions = NULL, eps=0.001, clust)
##  noccasions = NULL added 2011-04-04
{
  beta <- object$fit$par
  if (object$detectfn %in% c(0,2,9)) {    ## halfnormal, exponential and binary SS
                                          ## have just 2 parameters
      est <- beta[1:2]
      g   <- double(2)
  } else {                       ## other detectfn have 3 parameters
      est <- beta[1:3]
      g   <- double(3)
  }

  est <- beta
  g   <- beta

  ## consider replacing this with fdHess from package nlme

  grad <- function (i, est, eps) {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- CLmeanesa (est, object, individuals, sessnum, noccasions)
      est[i]  <- temp + delta
      fplus   <- CLmeanesa (est, object, individuals, sessnum, noccasions)
      (fplus - fminus) / (2.0 * delta)
  }
  if (is.null(clust))
      sapply(1:length(est), grad, est=est, eps=eps)
  else {
      clusterExport(clust, c("esa", "CLmeanesa","object","individuals","sessnum","noccasions"), environment())
      parSapply(clust, 1:length(est), grad, est = est, eps = eps)
  }
}

############################################################################################

derived <- function (object, sessnum = NULL, groups=NULL, alpha=0.05, se.esa = FALSE,
    se.D = TRUE, loginterval = TRUE, distribution = NULL, ncores = 1) {

## Generate table of derived parameters from fitted secr object

## modified 2009 07 21 for multiple sessions
## modified 2009 08 27 to report variance components

## multi-session -> separate each session, for the present
## is.null(sessnum) implies all possible sessions
## is.null(groups) implies all possible individuals

## multi-session -> pooled sessions not implemented

## groups within sessions
## groups to be found in covariates(object$capthist) cf

    if (inherits(object, 'secrlist')) {
        lapply(object, derived, sessnum, groups, alpha, se.esa, se.D,
               loginterval, distribution)
    }
    else {
        if (!is.null(distribution)) {
            if (tolower(distribution) %in% c('poisson','binomial'))
                object$details$distribution <- distribution
            else stop ("distribution not recognised")
        }
        if (!inherits(object, 'secr')) {
            warning ("requires fitted secr model")
            return (NULL)
        }
        else

        if (is.list(object$capthist) & is.null(sessnum)) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            jj <- match (sessnames, session(object$capthist))
            if (ncores > 1) {
                clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
                clusterEvalQ(clust, requireNamespace('secr'))
                output <- parLapply(clust, jj, derived, object = object, groups = groups,
                    alpha = alpha, se.esa = se.esa, se.D = se.D, loginterval = loginterval,
                    distribution = distribution, ncores = 1)
                stopCluster(clust)
            }
            else {
                output <- lapply(jj, derived, object = object, groups = groups,
                                 alpha = alpha, se.esa = se.esa, se.D = se.D,
                                 loginterval = loginterval, distribution = distribution, ncores = 1)
            }
            names(output) <- sessnames
            output
        }
        else {

            se.deriveD <- function (selection, selected.a, asess) {
                A <-  if (inherits(mask, 'linearmask'))
                    masklength(mask)
                else
                    ## maskarea(object$mask, asess)
                    maskarea(mask)
                s2 <- switch (tolower(object$details$distribution),
                   poisson  = sum (1/selected.a^2),
                   binomial = sum (( 1 - selected.a / A) / selected.a^2))
                CLg  <- CLgradient (object, selection, asess, clust=clust)
                varDn <- CLg %*% object$beta.vcv %*% CLg
                list(SE=sqrt(s2 + varDn), s2=s2, varDn=varDn)
            }
            se.deriveesa <- function (selection, asess) {
                CLesa  <- esagradient (object, selection, asess, clust=clust)
                sqrt(CLesa %*% object$beta.vcv %*% CLesa)
            }
            weighted.mean <- function (a) {
                ## allows for varying a 2010-11-30
                length(a) / sum(1/a)
            }

            getderived <- function (selection, NT = 0) {
                if (is.null(sessnum)) sessnum <- 1
                selected.a <- esa(object, sessnum)[selection]
                derivedmean <- c(weighted.mean(selected.a), sum(1/selected.a) )
                derivedSE <- c(NA, NA)
                varcomp1 <- c(NA, NA)
                varcomp2 <- c(NA, NA)
                if (se.esa) derivedSE[1] <- se.deriveesa(selection, sessnum)
                if (se.D) {
                    varDlist <- se.deriveD(selection, selected.a, sessnum)
                    derivedSE[2] <- varDlist$SE
                    varcomp1[2] <- varDlist$s2
                    varcomp2[2] <- varDlist$varDn
                }
                A <-  if (inherits(mask, 'linearmask'))
                    masklength(mask)
                else
                    maskarea(mask)
                temp <- data.frame (
                                    row.names = c('esa','D'),
                                    estimate = derivedmean + c(0,NT/A),  ## NT in 'telemetry' below
                                    SE.estimate = derivedSE)

                temp <- add.cl(temp, alpha, loginterval)

                temp$CVn <- varcomp1^0.5 / temp$estimate
                temp$CVa <- varcomp2^0.5 / temp$estimate
                temp$CVD <- temp$SE.estimate / temp$estimate
                temp$CVD[1] <- NA   ## not for esa

                nmash <- attr(capthist, 'n.mash')
                ## no need to allow for unmash as Density not a parameter
                if (!is.null(nmash)) {
                        temp[2,1:4] <- temp[2,1:4] / length(nmash)
                        ## message ("D was adjusted for ", length(nmash), " mashed clusters\n")
                    }
                temp
            }


            ## mainline

            if (ncores > 1) {
                clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
                clusterEvalQ(clust, requireNamespace('secr'))
            }
            else {
                clust <- NULL
            }

            if (is.null(sessnum)) {
                capthist <- object$capthist
                mask <- object$mask
            }
            else {
                capthist <- object$capthist[[sessnum]]
                mask <- object$mask[[sessnum]]
            }

            grp <- group.factor(capthist, groups)
            ind <- 1:nrow(capthist)

            if (nrow(capthist)>0)
                individuals <- split (ind, grp)
            else
                individuals <-  split (numeric(0), grp) ## list of empty grp levels
            n <- length(individuals)   ## number of groups

            xyl <- telemetryxy(capthist)   ## attr(capthist,'xylist')
            if ( is.null(xyl) | is.null(object$details$telemetrytype) )
                object$details$telemetrytype <- 'none'
            if (object$details$telemetrytype == 'concurrent') {
                telem <- telemetered(capthist)
                nottelem <- split(!telem, grp)
                NT <- tapply(telem, grp, sum)
                ## drop telemetry animals
                individuals <- mapply ('[', individuals, nottelem, SIMPLIFY = FALSE)
            }
            else {
                NT <- numeric(n)  ## zeros
            }
            if ( n > 1)
                out <- mapply (getderived, individuals, NT, SIMPLIFY = FALSE)
            else {
                if (n == 1)
                    out <- getderived(individuals[[1]], NT)
                else
                    out <- getderived(numeric(0), NT)
            }
            if (ncores>1) {
                stopCluster(clust)
            }
            out
        }
    }
}
############################################################################################

