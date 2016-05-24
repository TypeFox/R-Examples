smc <- function(object, ncomps = object$ncomp, corrected = F) {
  Coefs <- as.matrix(object$coefficients)[, ncomps]
  new.t.smc <- as.matrix(object$Xdata) %*% ((Coefs) / sqrt(crossprod((Coefs))))
  X.hat.smc <- new.t.smc %*% t(as.vector(((Coefs) / sqrt(crossprod((Coefs))))))
  X.error.smc <- object$Xdata - X.hat.smc
  if(corrected == T) {
    My.Cols <- 1:ncol(X.error.smc)
    Corrections.sMC <- llply(My.Cols, function(this.col) {
      acf.ci <- function(COLUMN, ci) {
        clim0 <- qnorm((1 + ci) / 2)/sqrt(length(COLUMN))
        clim <- clim0 * sqrt(cumsum(c(1, 2 * acf(COLUMN, plot = F)$acf[-1, 1, 1]^2)))
        ACF.LB <- max(-clim)
        ACF.UB <- min(clim)
        return(c(ACF.LB, ACF.UB))
      }
      ACF <- Significant <- NULL
      Sig.ACFs <- data.frame(acf(X.error.smc[, this.col], plot = F)$acf[-1, 1, 1], 
                      ifelse(acf(X.error.smc[, this.col], plot = F)$acf[-1, 1, 1] > acf.ci(X.error.smc[, this.col], .95)[1] &
                      acf(X.error.smc[, this.col], plot = F)$acf[-1, 1, 1] < acf.ci(X.error.smc[, this.col], .95)[2], "Not.Sig", "Sig"))
      names(Sig.ACFs)[c(1:2)] <- c("ACF", "Significant")
      Sig.ACFs.F <- subset(Sig.ACFs, Significant == "Sig" & ACF > 0)
      Sig.ACFs.F2 <- Sig.ACFs.F[diff(c(0, as.numeric(rownames(Sig.ACFs.F)))) == 1, ][, 1]
      if(length(Sig.ACFs.F2) == 0) {
        K <- 0
      } else {
        K <- 1:length(Sig.ACFs.F2)
      }
      Num.adj <- 1
      ACF.MEAN <- (1 / sum((1 - (K / nrow(X.error.smc))))) * sum((1 - (K / nrow(X.error.smc))) * (Sig.ACFs.F2))
      if(length(Sig.ACFs.F2) == 0) {
        Dem.adj <- 1
      } else {
        Dem.adj <- ((1 + ACF.MEAN) / (1 - ACF.MEAN))
      }
      c(Num.adj, Dem.adj)
    }); Corrections.sMC <- do.call("rbind", Corrections.sMC)
    names(Corrections.sMC) <- c("Num.adj", "Dem.adj")
    smc.df <- data.frame(smc = ((apply(X.hat.smc, 2, function(x) sum((x)^2)) / 1)) / 
      (((apply(X.error.smc, 2, function(x) sum((x)^2))) / (nrow(model.frame(object)) - 2)) * Corrections.sMC[, 2]))
  } else {
    smc.df <- data.frame(smc = ((apply(X.hat.smc, 2, function(x) sum((x)^2)) / 1)) / 
        (((apply(X.error.smc, 2, function(x) sum((x)^2))) / (nrow(model.frame(object)) - 2))))
  }
  smc.df$p.value <- 1-pf(smc.df$smc, 1, (nrow(X.error.smc) - 2))
  smc.df$f.value <- qf(.95,1, (nrow(X.error.smc) - 2))
  smc.df$Significant <- ifelse(smc.df$p.value < 0.05, "Yes", "No")
  
Results <- list(smc = smc.df, modeled = X.hat.smc, error = X.error.smc)
class(Results) <- "smc"
Results
}

print.smc <- function(x, ...) {
  cat("\nsmc summary:\n")
  print(x$smc[order(x$smc$smc, decreasing = TRUE), ])
}

smc.modeled <- function(x, ...) {
  cat("\nsmc model space:\n")
  print(x$modeled)
}

smc.error <- function(x, ...) {
  cat("\nsmc model error:\n")
  print(x$error)
}

plot.smc <- function(x, variables = "all", ...) {
  if(variables == "all") {
    df <- x$smc[order(x$smc$smc, decreasing = TRUE), ]
  } else {
    df <- x$smc[order(x$smc$smc, decreasing = TRUE), ][1:variables, ]
  }
  df$Variables <- factor(row.names(df), levels = row.names(df)[order(df$smc)])
  print(with(df, ggplot(df, aes(smc, Variables)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point(aes(col = Significant), size = 5) + 
    ggtitle("sMC Results") + 
    xlab("F-statistics") + 
    theme(plot.title = element_text(size = 20)) + 
    theme(axis.title.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20, angle = 90)) + 
    theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
    theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
}



