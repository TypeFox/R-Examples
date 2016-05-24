#' Calculate the Standardized total gain
#'
#' Computes the standardized total gain for the risk difference. Optionally produces bootstrap standard errors, and permutation test statistic. The standardized total gain is the area between the risk difference curve and the horizontal line at the marginal risk difference.
#'
#' @param psdesign A psdesign object. It must contain a risk model, an
#'   integration model, and estimated parameters. Bootstrapped parameters are
#'   optional
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param sig.level Significance level for bootstrap confidence intervals
#' @param n.samps The number of samples to take over the range of S.1 at which
#'   the VE is calculated
#' @param bootstraps If true, and bootstrapped estimates are present, will
#'   calculate bootstrap standard errors and confidence bands.
#' @param permute If true, will do permutation test for whether the STG is different from 0
#' @param permute.times Numeric, number of permutations to run
#' @param progress.bar Logical, if true will display a progress bar in the console
#'
#' @export
#'
calc_STG <- function(psdesign, t, sig.level = .05, n.samps = 5000, bootstraps = TRUE, permute = TRUE, permute.times = 2000, progress.bar = TRUE){

  stopifnot("estimates" %in% names(psdesign))

  impped <- psdesign$integration.models$S.1$icdf_sbarw(runif(n.samps))
  randrows <- sample(1:nrow(impped), ncol(impped), replace = TRUE)
  integrated <- impped[cbind(randrows, 1:n.samps)]
  obss <- psdesign$augdata$S.1

  trueobs <- sample(obss[!is.na(obss)],
                    floor(n.samps * mean(!is.na(obss))),
                    prob = psdesign$augdata$cdfweights[!is.na(obss)], replace = TRUE)


  if(is.factor(obss)){
    integrated <- factor(integrated, levels = levels(obss))
  }
  Splot <- sort(unlist(list(integrated, trueobs)))

  dat1 <- data.frame(S.1 = Splot, Z = 1)
  dat0 <- data.frame(S.1 = Splot, Z = 0)


  if(is.null(psdesign$risk.model$args$model)){

    others <- NULL

  } else {

    mod <- eval(psdesign$risk.model$args$model)
    mterms0 <- rownames(attr(terms(mod), "factors"))[-1]
    mterms <- unlist(lapply(colnames(psdesign$augdata), function(x){
      if(length(grep(x, mterms0, fixed = TRUE) > 0)){
        return(x)
      } else return(NULL)
    }))

    others <- mterms[!mterms %in% c("S.1", "Z")]
    for(j in others){

      tempval <- median(psdesign$augdata[, j], na.rm = TRUE)
      dat1[, j] <- tempval
      dat0[, j] <- tempval

    }

  }

  risks <- riskcalc(psdesign$risk.function, psdesign$augdata$Y, psdesign$estimates$par, t, dat0, dat1)

  obsdelta <- data.frame(S.1 = Splot, delta = risks$R1 - risks$R0)
  obstheta <- mean(obsdelta$delta)

  obsSTG <- stg(risks$R1, risks$R0, TRUE)

  retSTG <- list(obsSTG = obsSTG, bootstraps = NULL, permutation = NULL)

  if(bootstraps && "bootstraps" %in% names(psdesign)){

    bsests <- psdesign$bootstraps
    bootSTGs <- matrix(NA, nrow = nrow(bsests), ncol = 2)

    for(i in 1:nrow(bsests)){

      thispar <- as.numeric(bsests[i, -ncol(bsests)])

      thisrisk <- riskcalc(psdesign$risk.function, psdesign$augdata$Y,
                           thispar, t, dat0, dat1)

      bootSTGs[i, ] <- c(stg(thisrisk$R1, thisrisk$R0, TRUE), bsests[i, "convergence"])

    }

    bootSTGs <- as.data.frame(bootSTGs)
    colnames(bootSTGs)[ncol(bootSTGs)] <- "convergence"
    A1 <- as.data.frame(summarize_bs(bootSTGs, obsSTG, sig.level = sig.level, CI.type = "pointwise")$table)

    colnames(A1) <- gsub("%", "", paste("STG", colnames(A1), sep = "."), fixed = TRUE)

    retSTG$bootstraps <- A1

  }

  if(permute){

    if(progress.bar){
      cat(paste("Permuting", permute.times, "replicates:\n"))
      pb <- txtProgressBar(min = 1, max = permute.times)
    }

    perm.STG <- rep(NA, permute.times)
    for(i in 1:length(perm.STG)){

#       env.copy <- new.env(parent = parent.env(environment(psdesign$likelihood)))
#       ps.copy <- psdesign
#
#       objs <- ls(environment(psdesign$likelihood))
#       for(j in objs){
#         assign(j, get(j, environment(psdesign$likelihood)), envir = env.copy)
#       }
#
#       assign("Y.trt", mixup(env.copy$Y.trt, env.copy$trtmat[, "Z"]), envir = env.copy)
#       assign("Y.untrt", mixup(env.copy$Y.untrt, env.copy$untrt.expand[, "Z"]), envir = env.copy)
#
#       environment(ps.copy$likelihood) <- env.copy
#       perm.est <- ps.copy + eval(ps.copy$estimate.call)

      psdesign.0 <- psdesign
      psdesign.0$augdata$Y <- mixup(psdesign.0$augdata$Y, psdesign.0$augdata$Z)

      if(is.factor(psdesign.0$augdata$S.1)){

        if(length(unique(as.numeric(psdesign.0$augdata$S.1))) != length(unique(as.numeric(psdesign$augdata$S.1)))){
          bootpar[[i]] <- c(rep(NA, psdesign$nparam), convergence = 11)
          next
        }

      }
      ## re-call integration models

      psdesign2 <- psdesign.0
      for(intj in psdesign$integration.models){
        psdesign2 <- psdesign2 + do.call(as.character(intj$model$args[[1]]), intj$model$args[-1])
      }

      ## re-call risk model

      psdesign3 <- psdesign2 + do.call(as.character(psdesign$risk.model$args[[1]]), psdesign$risk.model$args[-1])

      # estimate

      perm.est <- psdesign3 + eval(psdesign$estimate.call)

      risks.perm <- riskcalc(perm.est$risk.function, perm.est$augdata$Y, perm.est$estimates$par, t, dat0, dat1)

      perm.STG[i] <- stg(risks.perm$R1, risks.perm$R0, TRUE)

      if(progress.bar){
        setTxtProgressBar(pb, value = i)
        flush.console()
        }

    }
    if(progress.bar){
      close(pb)
  }
    perm.p <- mean(perm.STG > abs(obsSTG))
    retSTG$permutation <- structure(list(p.value = perm.p, permuted.stats = perm.STG), class = "permutation")


  }

  retSTG

}

#' @export

print.permutation <- function(x, ...){

  print(paste0("permutation p = ", x$p.value))

}

mixup <- function(x, z){

  z0 <- x[z == 0]
  z1 <- x[z == 1]

  x[z == 0] <- z0[sample(1:length(z0), length(z0), replace = FALSE)]
  x[z == 1] <- z1[sample(1:length(z1), length(z1), replace = FALSE)]
  x

}

#' Compute the standardized total gain
#' @param R1 Risk in the treatment group
#' @param R0 Risk in the control group
#' @param stand Standardize?

stg <- function(R1, R0, stand = TRUE){

  delt <- R0 - R1
  if(stand){
    den <- mean(delt) * (1 - mean(delt)) * 2
  } else {
    den <- 1
  }
  mean(abs(delt - mean(delt))) / den

}