#' Summary for posterior of a \code{pogit} object
#' 
#' Returns basic information about the model and the priors, MCMC details and 
#' (model averaged) posterior means with 95\%-HPD intervals for the regression 
#' effects and estimated posterior inclusion probabilities. 
#' 
#' To assess mixing and efficiency of MCMC sampling, the effective sample size 
#' (ESS) and the integrated autocorrelation time (IAT) are computed. ESS 
#' estimates the equivalent number of independent draws corresponding to the 
#' dependent MCMC draws and is defined as ESS = \eqn{M}/\eqn{\tau}, where \eqn{\tau} 
#' is the IAT and \eqn{M} is the number of MCMC iterations after the burn-in phase. 
#' IAT is computed as \eqn{\tau = 1 + 2 \sum_{k=1}^K \rho(k)}
#' using the initial monotone sequence estimator (Geyer, 1992) for K and 
#' \eqn{\rho(k)} is the empirical autocorrelation at lag \eqn{k}.
#' 
#' @param object an object of class \code{pogit}
#' @param IAT if \code{TRUE}, integrated autocorrelation times (IAT) and 
#'   effective samples sizes (ESS) of the MCMC samples are computed (see 
#'   details); defaults to \code{FALSE}. 
#' @param printRes if \code{TRUE}, model averaged posterior means for the 
#'  reporting probabilities and risks are computed for the Pogit model; 
#'  defaults to \code{FALSE}.
#' @param ... further arguments passed to or from other methods (not used)
#'   
#' @author Michaela Dvorzak <m.dvorzak@@gmx.at>
#' @return an object of class \code{summary.pogit}
#' 
#' @references Geyer, C. J. (1992). Practical Markov Chain Monte Carlo. 
#'  \emph{Statistical Science}, \strong{7}, 473-483.
#'   
#' @aliases print.summary.pogit
#' @export
 
summary.pogit <- function(object, IAT = FALSE, printRes = FALSE, ...){
  stopifnot(class(object) == "pogit")
  
  if (object$family %in% c("logit", "pogit")){
    samples <- lapply(object$samplesL, thinMCMC, start = object$mcmc$burnin + 1, 
                      object$mcmc)
    samples$absThetaAlpha <- NULL
    if(!is.null(samples$thetaAlpha)) samples$absThetaAlpha <- abs(samples$thetaAlpha)
    postMeansL <- lapply(samples, function(x){
      if (!is.null(x[[1]])){
        return(colMeans(x))
      }
    })  
    hpdL <- lapply(samples[c("alpha", "absThetaAlpha")], function(x){
      if (!is.null(x[[1]])){
        return(t(apply(x, MARGIN = 2, hpdMCMC)))
      }
    })
    
    if (IAT){
      iatL <- lapply(samples[c("alpha", "pdeltaAlpha", "pgammaAlpha")], function(x){
        if (!is.null(x[[1]])) return(t(apply(x, MARGIN = 2, iatMCMC)))
      })
    }
    
    if (printRes && object$family == "pogit"){
      muL <- object$data$W%*%postMeansL$alpha
      if (object$model.logit$ri==1){
        linp <- muL + postMeansL$ai
      } else linp <- muL
      p.est <- exp(linp)/(1 + exp(linp))
    }
  }
  
  if (object$family %in% c("pogit", "poisson")){
    samples <- lapply(object$samplesP, thinMCMC, start = object$mcmc$burnin + 1, object$mcmc)
    samples$absThetaBeta <- NULL
    if(!is.null(samples$thetaBeta)) samples$absThetaBeta <- abs(samples$thetaBeta)
    postMeansP <- lapply(samples, function(x){
      if (!is.null(x[[1]])){
        return(colMeans(x))
      }
    })
    hpdP <- lapply(samples[c("beta", "absThetaBeta")], function(x){
      if (!is.null(x[[1]])){
        return(t(apply(x, MARGIN = 2, hpdMCMC)))
      }
    })
    
    if (IAT){
      iatP <- lapply(samples[c("beta", "pdeltaBeta", "pgammaBeta")], function(x){
        if (!is.null(x[[1]])) return(t(apply(x, MARGIN = 2, iatMCMC)))
      })
    }
    
    if (printRes && object$family == "pogit"){
      muP <- object$data$X%*%postMeansP$beta
      if (object$model.pois$ri == 1){
        linp <- muP + postMeansP$bi
      } else linp <- muP
      lambda.est <- exp(linp)
    }   
  }
  
  if (object$family == "negbin"){
    samples <- lapply(object$samplesNB, thinMCMC, start = object$mcmc$burnin + 1, 
                      object$mcmc)
    postMeansNB <- lapply(samples, function(x){
      if (!is.null(x[[1]])){
        return(colMeans(x))
      }
    })
    hpdNB <- lapply(samples[c("beta", "rho")], function(x){
      if (!is.null(x[[1]])){
        return(t(apply(x, MARGIN = 2, hpdMCMC)))
      }
    })
    
    if (IAT){
      iatNB <- lapply(samples[c("beta", "pdeltaBeta", "rho")], function(x){
        if (!is.null(x[[1]])) return(t(apply(x, MARGIN = 2, iatMCMC)))
      })
    }
    
  }
  
  if (object$family %in% c("logit", "pogit")){
    tabrow <- object$model.logit$d + object$model.logit$ri + 1
    tabcol <- 3 + as.numeric(object$BVS)
    resL <- c(postMeansL$alpha, postMeansL$absThetaAlpha)
    resHPDL <- rbind(hpdL$alpha, hpdL$absThetaAlpha)
    
    resmod <- data.frame(matrix(NA, nrow = tabrow, ncol = tabcol))
    rownames(resmod) <- names(resL)
    rownames(resmod)[1] <- "(Intercept)"
    colnres <- switch(as.character(object$BVS),
                      "TRUE" = c("Estimate","P(.=1)","95%-HPD[l]","95%-HPD[u]"),
                      "FALSE" = c("Estimate","95%-HPD[l]","95%-HPD[u]"))
    colnames(resmod) <- colnres
    resmod$Estimate <- round(resL, 3)
    if (tabcol > 3){
      resmod[-1,2] <- round(c(postMeansL$pdeltaAlpha, postMeansL$pgammaAlpha), 3)
    }
    resmod[,(tabcol-1):tabcol] <- round(resHPDL, 3)
    tabLogit <- as.matrix(resmod)
    
    if (IAT){
      if (!object$BVS) iatP$pdeltaAlpha <- NULL
      resIAT <- round(data.frame(rbind(iatL$alpha, iatL$pdeltaAlpha, iatL$pgammaAlpha)), 2) 
      IATLogit <- as.matrix(resIAT)
    }
    
    if (printRes && object$family == "pogit"){
      probs <- data.frame(prob = round(p.est, 3))
      rownames(probs) <- as.character(seq_len(nrow(object$data$W)))
    }
  }
  
  if (object$family %in% c("pogit", "poisson")){
    tabrow <- object$model.pois$d + object$model.pois$ri + 1
    tabcol <- 3 + as.numeric(object$BVS)    
    resP <- c(postMeansP$beta, postMeansP$absThetaBeta)
    resHPDP <- rbind(hpdP$beta, hpdP$absThetaBeta)
    
    resmod <- data.frame(matrix(NA, nrow = tabrow, ncol = tabcol))
    rownames(resmod) <- names(resP)
    rownames(resmod)[1] <- "(Intercept)"
    colnres <- switch(as.character(object$BVS),
                      "TRUE" = c("Estimate","P(.=1)","95%-HPD[l]","95%-HPD[u]"),
                      "FALSE" = c("Estimate","95%-HPD[l]","95%-HPD[u]"))
    colnames(resmod) <- colnres
    resmod$Estimate <- round(resP, 3)
    if (tabcol > 3){
      resmod[-1,2] <- round(c(postMeansP$pdeltaBeta, postMeansP$pgammaBeta), 3)
    }
    resmod[,(tabcol-1):tabcol] <- round(resHPDP, 3)
    tabPois <- as.matrix(resmod)
    
    if (IAT){
      if (!object$BVS) iatP$pdeltaBeta <- NULL
      resIAT <- round(data.frame(rbind(iatP$beta, iatP$pdeltaBeta, iatP$pgammaBeta)), 2) 
      IATPois <- as.matrix(resIAT)
    }
    
    if (printRes && object$family == "pogit"){
      risks <- data.frame(risk = round(lambda.est, 3))
      rownames(risks) <- as.character(seq_len(nrow(object$data$X)))
      if (all(object$data$E == 1)) colnames(risks) <- "intensity"
    }
  }
  
  if (object$family == "negbin"){
    tabrow <- object$model.nb$d + 2
    tabcol <- 3 + as.numeric(object$BVS)    
    resNB   <- c(postMeansNB$beta, postMeansNB$rho)
    resHPDNB <- rbind(hpdNB$beta, hpdNB$rho)
    
    resmod <- data.frame(matrix(NA, nrow = tabrow, ncol = tabcol))
    rownames(resmod) <- names(resNB)
    rownames(resmod)[1] <- "(Intercept)"
    colnres <- switch(as.character(object$BVS),
                      "TRUE" = c("Estimate","P(.=1)","95%-HPD[l]","95%-HPD[u]"),
                      "FALSE" = c("Estimate","95%-HPD[l]","95%-HPD[u]"))
    colnames(resmod) <- colnres
    resmod$Estimate <- round(resNB, 3)
    if (tabcol > 3){
      resmod[-c(1,tabrow),2] <- round(postMeansNB$pdeltaBeta, 3)
    }
    resmod[,(tabcol-1):tabcol] <- round(resHPDNB, 3)
    tabNB <- as.matrix(resmod)
    
    if (IAT){
      if (!object$BVS) iatNB$pdeltaBeta <- NULL
      resIAT <- round(data.frame(rbind(iatP$beta, iatP$pdeltaBeta, iatP$rho)), 2) 
      IATNB <- as.matrix(resIAT)
    }
    
  }
  
  if (object$family == "pogit"){
    rownames(tabLogit)[1] <- "(Intercept) Logit"
    rownames(tabPois)[1]  <- "(Intercept) Poisson"
    tabPogit <- rbind(tabPois, tabLogit)
    if (IAT) IATPogit <- rbind(IATPois, IATLogit)
  }
  
  modTable <- switch(as.character(object$family),
                     "logit"   = tabLogit,
                     "poisson" = tabPois,
                     "pogit"   = tabPogit, 
                     "negbin"  = tabNB, "\n")
  
  if (IAT){
    iatTable <- switch(as.character(object$family),
                       "logit"   = IATLogit,
                       "poisson" = IATPois,
                       "pogit"   = IATPogit, 
                       "negbin"  = IATNB, "\n")
  } else iatTable <- NULL
  
  if (printRes && object$family == "pogit"){
    resTable <- list(probs = probs, risks = t(risks))
  } else resTable <- NULL
  
  return(structure(
    c(list(modTable = modTable, iatTable = iatTable, resTable = resTable), 
      IAT = IAT, printRes = printRes, object), class = "summary.pogit"))
}  




#' @rdname summary.pogit 
#' @param x a \code{summary.pogit} object produced by \code{summary.pogit()} 
#' @export

print.summary.pogit <- function(x, ...){
  if (x$family=="logit"){
    bin <- "binomial "
    if (all(x$data$N==1)) bin <- ""
  }
  cat(paste(
    bvs <- switch(as.character(x$BVS), 
                  "TRUE"  = "Bayesian variable selection",
                  "FALSE" = "MCMC"), 
    "for the", switch(as.character(x$family),
                  "logit"   = paste(bin, "logit", sep=""),
                  "pogit"   = switch(as.character(x$fun), 
                                     "select_poisson" = "Pogit",
                                     "select_poissonOD" = "overdispersed Pogit"),
                  "poisson" = switch(as.character(x$fun), 
                                     "select_poisson" = "Poisson",
                                     "select_poissonOD" = "overdispersed Poisson"),
                  "negbin"  = "negative binomial"), "model:\n"))
  
  if (typeof(x$IAT) != "logical") stop("invalid 'IAT' argument") 
  
  cat("\nCall:\n")
  print(x$call)
  
  cat("\n\nMCMC:")
  cat("\nM =", x$mcmc$M, "draws after a burn-in of", x$mcmc$burnin)
  
  if (x$BVS) cat("\nBVS started after", x$mcmc$startsel, "iterations")
  cat("\nThinning parameter:", x$mcmc$thin)
  
  
  if(x$family == "negbin"){
    cat(paste0("\n\nAcceptance rate for rho:\n", x$acc.rho, "%"))
  }
  
  cat("\n\nPrior:")
  if (x$family %in% c("poisson", "pogit")){
    pr.txt <- switch(as.character(x$BVS),
                     "TRUE" = paste("spike-and-slab prior with",
                                    switch(x$prior.pois$slab,
                                           "Normal" = paste0(x$prior.pois$slab, " slab"),
                                           "Student" = paste0(x$prior.pois$slab, "-t slab")
                                    )
                     ),
                     "FALSE" = paste(switch(x$prior.pois$slab, 
                                            "Normal" = paste0(x$prior.pois$slab, " prior"),
                                            "Student" = paste0(x$prior.pois$slab, "-t prior")
                     )
                     )
    )
    
    if (x$family=="pogit"){
      cat("\n")
      cat(paste0("- Poisson: ", pr.txt, " [V=", x$prior.pois$V, "]\n"))
    } else {
      cat(paste0(" ", pr.txt, " [V=", x$prior.pois$V, "]\n"))
    }
    
    if (x$model.pois$d > 0){  
      pM <- as.character(paste0("b0[", 0:(x$model.pois$d + x$model.pois$ri), "]"))
      priorSet2 <- round(with(x$prior.pois, c(priorMean = c(unname(m0), unname(aj0)))), 3)
      names(priorSet2) <- pM
      cat("\n")
      print(priorSet2)
    }
    
    priorSet <- with(x$prior.pois, 
                     c("w[a]" = unname(w[1]), "w[b]" = unname(w[2]),
                       "pi[a]" = unname(pi[1]), "pi[b]" = unname(pi[2]))
    )
    if (x$BVS) print(priorSet)
  }
  
  if (x$family %in% c("logit", "pogit")){
    bvsl <- as.character(x$BVS)
    if (x$family=="pogit" && x$method=="infprior") bvsl <- "FALSE"
    pr.txt <- switch(bvsl,
                     "TRUE" = paste(" spike-and-slab prior with", 
                                    switch(x$prior.logit$slab,
                                           "Normal" = paste0(x$prior.logit$slab, " slab"),
                                           "Student" = paste0(x$prior.logit$slab, "-t slab")
                                    )
                     ),
                     "FALSE" = paste(switch(x$prior.logit$slab,
                                            "Normal" = paste0(x$prior.logit$slab, " prior"),
                                            "Student" = paste0(x$prior.logit$slab, "-t prior")
                                    )
                     )
    )
    
    if (x$family == "pogit"){
      cat("\n")
      cat(paste0("- Logit: ", pr.txt, " [V=", x$prior.logit$V, "]\n"))
    } else {
      cat(paste0(" ", pr.txt, " [V=", x$prior.logit$V, "]\n"))
    }
    
    priorSet <- with(x$prior.logit, 
                     c("w[a]" = unname(w[1]), "w[b]" = unname(w[2]),
                       "pi[a]" = unname(pi[1]), "pi[b]" = unname(pi[2]))
    )
    
    if (x$model.logit$d > 0){               
      pM <- as.character(paste0("a0[", 0:(x$model.logit$d + x$model.logit$ri), "]"))
      priorSet2 <- round(with(x$prior.logit, c(priorMean = c(unname(m0), unname(aj0)))), 3)
      names(priorSet2) <- pM
      cat("\n")
      print(priorSet2)
    }
    
    if (x$BVS) print(priorSet)  
  }
  
  if (x$family == "negbin"){
    pr.txt <- switch(as.character(x$BVS),
                     "TRUE" = paste("spike-and-slab prior with",
                                    switch(x$prior.nb$slab,
                                           "Normal" = paste0(x$prior.nb$slab, " slab"),
                                           "Student" = paste0(x$prior.nb$slab, "-t slab")
                                    )
                     ),
                     "FALSE" = paste(switch(x$prior.nb$slab, 
                                            "Normal" = paste0(x$prior.nb$slab, " prior"),
                                            "Student" = paste0(x$prior.nb$slab, "-t prior")
                     )
                     )
    )
    
    cat(paste0(" ", pr.txt, " [V=", x$prior.nb$V, "]\n"))
    
    if (x$model.nb$d > 0){  
      pM <- as.character(paste("b0[", 0:(x$model.nb$d), "]",sep=""))
      priorSet2 <- round(with(x$prior.nb, c(priorMean = c(unname(m0), unname(aj0)))), 3)
      names(priorSet2) <- pM
      cat("\n")
      print(priorSet2)
    }
    
    if (x$BVS){
      priorSet <- with(x$prior.nb, 
                       c("w[a]" = unname(w[1]), "w[b]" = unname(w[2]))
      )
    } else priorSet <- NULL
    
    priorSet0 <- with(x$prior.nb, 
                      c("c0" = unname(c0), "C0" = unname(C0))
    )
    
    print(c(priorSet, priorSet0))
  }
  
  
  cat(paste(strwrap(paste(
    switch(as.character(x$BVS),
           "TRUE"  = "\n\nModel averaged posterior means, estimated posterior 
           inclusion probabilities",
           "FALSE" = "\n\nPosterior means"),
    "and 95%-HPD intervals:"), exdent = 1), collapse = "\n"))
  cat("\n\n")
  print(x$modTable)
  
  
  # print integrated autocorrelation times of samples 
  if (x$IAT){
    cat(paste(
      strwrap(paste("\n\nIntegrated autocorrelation times (IAT) and effective
                    sample size (ESS):"), 
              exdent = 1), collapse = "\n"))
    cat("\n\n")
    print(x$iatTable)
  }
  
  if (x$printRes && x$family=="pogit"){
    cat("\n\nEstimated reporting probabilities:\n")
    print(x$resTable$probs)
    cat("\n\nEstimated risks:\n")
    print(x$resTable$risks)
  }
  
  cat("\n")
  invisible(x)  
}