defacto <-
function (model, plot = TRUE, axes = c(1, 2), select.covar = NULL,
   select.annot = NULL, lim.b = 0.01, lab = TRUE, cex = 1) {
   palette(c("black", "red", "green3", "blue", "cyan", "magenta",
       "darkgray", "darkgoldenrod", "darkgreen", "violet", "turquoise",
       "orange", "lightpink", "lavender", "yellow", "lightgreen",
       "lightgrey", "lightblue", "darkkhaki", "darkmagenta",
       "darkolivegreen", "lightcyan", "darkorange", "darkorchid",
       "darkred", "darksalmon", "darkseagreen", "darkslateblue",
       "darkslategray", "darkslategrey", "darkturquoise", "darkviolet",
       "lightgray", "lightsalmon", "lightyellow", "maroon"))
   if (class(model)[1] != "FAMTmodel")
       stop("Class of model should be FAMTmodel")
   m = nrow(model$adjdata$expression)
   annot = data.frame(model$adjdata$annotations)
   covar = data.frame(model$adjdata$covariates)
   if (model$nbf <= 1)
       plot = FALSE
   if (length(axes) != 2)
       stop(paste("axes should be 2 values in 1:", model$nbf,
           sep = ""))
   if (!any(is.element(axes, 1:model$nbf)))
       stop(paste("axes should be in 1:", model$nbf, sep = ""))
   if (model$nbf == 0)
       stop("Number of factors should be at least 1.")
   if (model$nbf == 1)
       plot = FALSE
   if ((lim.b <= 0) | (lim.b > 1))
       stop("lim.b should be in ]0,1]")
   if (is.null(select.covar))
       select.covar = (1:ncol(model$adjdata$covariates))[-model$x]
   if (!is.null(select.annot)) {
       if (!any(is.element(1:ncol(annot), select.annot)))
           stop(paste("select.annot should be NULL or in 1:",
               ncol(annot), sep = ""))
   }
   if (is.null(select.annot))
       select.annot = (1:ncol(model$adjdata$annotations))[which(names(model$adjdata$annotations) !=
           "ID")]
   if (lab == TRUE)
       idlab = (1:ncol(model$adjdata$annotations))[which(names(model$adjdata$annotations) ==
           "ID")]
   if (!any(is.element(1:ncol(covar), select.covar)))
       stop(paste("select.covar should be NULL or in 1:", ncol(covar),
           sep = ""))
   select.covar = select.covar[!is.element(select.covar, model$idcovar)]
   isfactor.annot = unlist(lapply(select.annot,function(k,data) is.factor(data[,k]),data=annot))
   select.annot = select.annot[isfactor.annot]
   isfactor.covar = unlist(lapply(select.covar,function(k,data) is.factor(data[,k]),data=covar))
   if (!is.null(isfactor.covar)) quanti.sup = select.covar[!isfactor.covar]
   if (is.null(isfactor.covar)) quanti.sup = integer(0)
   quali.sup = select.covar[isfactor.covar]
   if (plot == TRUE) {
       b2 = apply(model$FA$B[, axes]^2, 1, sum)
       limb2 = quantile(b2, 1 - lim.b)
       select = b2 >= limb2
       coord = model$FA$B[select, axes]
       plot(0, 0, xlab = paste("Factor", axes[1]), ylab = paste("Factor",
           axes[2]), xlim = c(-1, 1), ylim = c(-1, 1), col = "white",
           asp = 1, cex = cex, main = paste("Loadings", axes[1],
               "and", axes[2]))
       x.cercle <- seq(-1, 1, by = 0.01)
       y.cercle <- sqrt(1 - x.cercle^2)
       lines(x.cercle, y = y.cercle)
       lines(x.cercle, y = -y.cercle)
       abline(v = 0, lty = 2, cex = cex)
       abline(h = 0, lty = 2, cex = cex)
       points(coord, pch = 16, col = "blue", cex = 0.75)
       pos = rep(1, sum(select))
       pos[(abs(coord[, 1]) > abs(coord[, 2])) & (coord[, 1] >=
           0)] = 4
       pos[(abs(coord[, 1]) > abs(coord[, 2])) & (coord[, 1] <
           0)] = 2
       pos[(abs(coord[, 1]) <= abs(coord[, 2])) & (coord[, 2] >=
           0)] = 3
       text(coord, labels = annot[select, idlab], pos = pos,
           cex = cex)
       if (length(quali.sup)>0) {
          for (k in 1:length(quali.sup)) {
          dev.new()
          coord = model$FA$Factors[, axes]
          plot(coord, xlab = paste("Factor", axes[1]), ylab = paste("Factor",
               axes[2]), col = "white", asp = 1, cex = 1, main = paste("factors",
               axes[1], "and", axes[2]), bty = "l")
          pos = rep(1, ncol(model$adjdata$expression))
          pos[(abs(coord[, 1]) > abs(coord[, 2])) & (coord[,
               1] >= 0)] = 4
          pos[(abs(coord[, 1]) > abs(coord[, 2])) & (coord[,
               1] < 0)] = 2
          pos[(abs(coord[, 1]) <= abs(coord[, 2])) & (coord[,
               2] >= 0)] = 3
          if (lab == TRUE)
               text(coord, label = names(model$adjdata$expression),
                 cex = cex, pos = pos)
          points(coord, col = as.numeric(covar[, quali.sup[k]]),
               pch = 16)
          legend("bottomright", pch = rep(16, length(levels(covar[,
               quali.sup[k]]))), col = as.numeric(as.factor(levels(covar[,
               quali.sup[k]]))), legend = paste(names(covar)[quali.sup[k]],
               levels(covar[, quali.sup[k]])))
          }
       }       

if (length(quali.sup)==0) {
          dev.new()
          coord = model$FA$Factors[, axes]
          plot(coord, xlab = paste("Factor", axes[1]), ylab = paste("Factor",
               axes[2]), col = "white", asp = 1, cex = 1, main = paste("factors",
               axes[1], "and", axes[2]), bty = "l")
          pos = rep(1, ncol(model$adjdata$expression))
          pos[(abs(coord[, 1]) > abs(coord[, 2])) & (coord[,
               1] >= 0)] = 4
          pos[(abs(coord[, 1]) > abs(coord[, 2])) & (coord[,
               1] < 0)] = 2
          pos[(abs(coord[, 1]) <= abs(coord[, 2])) & (coord[,
               2] >= 0)] = 3
          if (lab == TRUE)
               text(coord, label = names(model$adjdata$expression),
                 cex = cex, pos = pos)
          points(coord, pch = 16)
       }
   }
   resloadings = vector(length = model$nbf, "list")
   for (j in 1:model$nbf) {
       limit = quantile(abs(model$FA$B[, j]), 1 - lim.b)
       select = (1:m)[abs(model$FA$B[, j]) > limit]
       ord = order(model$FA$B[select, j])
       resloadings[[j]] = data.frame(Row = select, B = model$FA$B[select,
           j], ID = annot$ID[select])[ord, ]
   }
   if (length(select.annot)==0) resannot = NULL
   if (length(select.annot)>0) {
      resannot = matrix(0, nrow = model$nbf, ncol = length(select.annot))
      colnames(resannot) = names(annot)[select.annot]
      rownames(resannot) = paste("Loadings", 1:model$nbf)
      for (j in 1:model$nbf) {
         for (k in 1:length(select.annot)) {
            blm = lm(model$FA$B[, j] ~ annot[, select.annot[k]])
            resannot[j, k] = anova(blm)[1, 5]
         }
       }
    }
  if (length(select.covar) == 0) rescovar = NULL
  if (length(select.covar) > 0) {
     rescovar = matrix(0, nrow = model$nbf, ncol = length(select.covar))
     colnames(rescovar) = names(covar)[select.covar]
     rownames(rescovar) = paste("Factor", 1:model$nbf)
     for (j in 1:model$nbf) {
        for (k in 1:length(select.covar)) {
           flm = lm(model$FA$Factors[, j] ~ covar[, select.covar[k]])
           rescovar[j, k] = anova(flm)[1, 5]
        }
     }
   }
   res = list(loadings = resloadings, covariates = rescovar,annotations = resannot)
   return(res)
}
