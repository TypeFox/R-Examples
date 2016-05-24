`gsci.ratio` <-
function(est, vcmat, Num.Contrast, Den.Contrast,
                      degfree=NULL, conf.level=0.95, alternative="two.sided", adjusted=TRUE){
    
    
    if (length(est) != ncol(Num.Contrast) | length(est) != ncol(Den.Contrast) | length(est) != ncol(vcmat)) stop("The length of est and the column number of \n\tvcmat, Num.Contrast, and Den.Contrast\n\thave to be the same!")
    if (nrow(Den.Contrast) != nrow(Num.Contrast)) stop("Num.Contrast and Den.Contrast need the same row number!")

    n.comp <- nrow(Num.Contrast)
    gammaC.vec <- (Num.Contrast %*% est)/Den.Contrast %*% est
    CorrMat.plug <- matrix(rep(NA, n.comp * n.comp), nrow = n.comp)
    for (i in 1:n.comp) {
      for (j in 1:n.comp) {
        CorrMat.plug[i, j] <- (gammaC.vec[i] * Den.Contrast[i, ] -
                Num.Contrast[i, ]) %*% vcmat %*% (gammaC.vec[j] * Den.Contrast[j,
                ] - Num.Contrast[j, ])/(sqrt((gammaC.vec[i] * Den.Contrast[i,
                ] - Num.Contrast[i, ]) %*% vcmat %*% (gammaC.vec[i] * Den.Contrast[i,
                ] - Num.Contrast[i, ])) * sqrt((gammaC.vec[j] * Den.Contrast[j,
                ] - Num.Contrast[j, ]) %*% vcmat %*% (gammaC.vec[j] * Den.Contrast[j,
                ] - Num.Contrast[j, ])))
      }
    }
    Quad.root <- function(Aj, Bj, Cj) {
      Discrimi <- Bj^2 - 4 * Aj * Cj
      if ((Aj > 0) & (Discrimi >= 0)){
          Limit.s <- (-Bj + plus.minus * sqrt(Discrimi))/(2 * Aj)
      } else Limit.s <- NA
      return(Limit.s)
    }
    if (alternative == "two.sided") {
      side <- 2
      plus.minus <- c(-1, 1)
      if (adjusted){
        if (is.null(degfree)){
          Cplug <- qmvnorm(conf.level, interval = c(0, 10),
            corr = CorrMat.plug, mean = rep(0, n.comp), tail = "both.tails")$quantile
        } else {
          Cplug <- qmvt(conf.level, interval = c(0, 10), df = as.integer(degfree),
            corr = CorrMat.plug, delta = rep(0, n.comp), tail = "both.tails")$quantile
        }
      } else {
        if (is.null(degfree)){
          Cplug <- qnorm(1-((1-conf.level)/2))
        } else {
          Cplug <- qt(1-((1-conf.level)/2), df = degfree)
        }
      }
    }
    if ((alternative == "less") | (alternative == "greater")) {
      side <- 1
      if (alternative == "less") plus.minus <- 1 else plus.minus <- -1
      if (adjusted){
        if (is.null(degfree)){
          Cplug <- qmvnorm(conf.level, interval = c(0, 10),
            corr = CorrMat.plug, mean = rep(0, n.comp), tail = "lower.tail")$quantile
        } else {
          Cplug <- qmvt(conf.level, interval = c(0, 10), df = as.integer(degfree),
            corr = CorrMat.plug, delta = rep(0, n.comp), tail = "lower.tail")$quantile
        }
      } else {
        if (is.null(degfree)){
          Cplug <- qnorm(conf.level)
        } else {
          Cplug <- qt(conf.level, df = degfree)
        }
      }
    }
    quant <- Cplug
    PlugCL <- matrix(rep(NA, side * n.comp), nrow = n.comp)
    for (j in 1:n.comp) {
      AjPlug <- (Den.Contrast[j, ] %*% est)^2 - (Cplug^2) * Den.Contrast[j, ] %*% vcmat %*% Den.Contrast[j, ]
      BjPlug <- -2 * ((Num.Contrast[j, ] %*% est) * (Den.Contrast[j,] %*% est) -
                (Cplug^2) * Num.Contrast[j, ] %*% vcmat %*% Den.Contrast[j, ])
      CjPlug <- (Num.Contrast[j, ] %*% est)^2 - (Cplug^2) * Num.Contrast[j, ] %*% vcmat %*% Num.Contrast[j, ]
      PlugCL[j,] <- Quad.root(AjPlug[1,1], BjPlug[1,1], CjPlug[1,1])
    }
    sci.table <- data.frame(PlugCL)
    SCIs <- cbind(gammaC.vec, PlugCL)
    if (alternative == "two.sided") colnames(SCIs) <- c("estimate", "lower", "upper")
    if (alternative == "less") colnames(SCIs) <- c("estimate", "upper")
    if (alternative == "greater") colnames(SCIs) <- c("estimate", "lower")
    CI <- matrix(SCIs[,-1], nrow=nrow(Num.Contrast))
    colnames(CI) <- colnames(SCIs)[-1]
    estimate <- cbind(estimate=SCIs[,1])
    if (all(!is.null(rownames(Num.Contrast)))){
      rownames(CI) <- rownames(Num.Contrast)
      rownames(estimate) <- rownames(Num.Contrast)
    }

    out <- list()
    out$estimate <- estimate
    out$CorrMat.est <- CorrMat.plug
    out$Num.Contrast <- Num.Contrast
    out$Den.Contrast <- Den.Contrast
    out$conf.int <- CI
    out$NSD <- any(is.na(SCIs))
    out$alternative <- alternative
    out$conf.level <- conf.level
    out$compnames <- rownames(SCIs)
    out$methodname <- if (adjusted) paste("Simultaneous ",signif(conf.level*100,4),"% confidence intervals", sep="") else paste("Unadjusted ",signif(conf.level*100,4),"% confidence intervals", sep="")
    out$type <-  "User defined"
    out$method <- if (adjusted) "Plug" else "Unadj"
    out$df <-  degfree
    out$quantile <- quant
    class(out) <- c("sci.ratio","gsci.ratio")
    out
}

