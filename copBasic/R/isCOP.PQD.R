"isCOP.PQD" <- function(cop=NULL, para=NULL, uv=NULL, empirical=FALSE, verbose=TRUE, ...) {
      tau <- rho <- gin <- blm <- NULL
      blm <- blomCOP(cop=cop, para=para, ...)
      if(verbose) message("  Blomqvist's Beta = ",round(blm, digits=4))
      if(empirical) { # sample versions
         U <- para[,1]; V <- para[,2]; n <- length(U)
         gin <- ( sum(abs(rank(U)+rank(V)-n-1)) -
                  sum(abs(rank(U)-rank(V)))       ) / as.integer(n^2/2)
         if(verbose) message("      Gini's Gamma = ",round(gin, digits=4))

         rho <- cor(U, V, method="spearman")
         if(verbose) message("    Spearman's Rho = ",round(rho, digits=4))

         tau <- cor(U, V, method="kendall")
         if(verbose) message("     Kendall's Tau = ",round(tau, digits=4))
      } else {
         try(gin <- giniCOP(cop=cop, para=para, ...), silent=TRUE)
         if(is.null(gin)) gin <- giniCOP(cop=cop, para=para, by.concordance=TRUE, brute=TRUE, ...)
         if(verbose) message("      Gini's Gamma = ",round(gin, digits=4))

         try(rho <- rhoCOP(cop=cop, para=para, ...), silent=TRUE)
         if(is.null(rho)) rho <- rhoCOP(cop=cop, para=para, brute=TRUE, ...)
         if(verbose) message("    Spearman's Rho = ",round(rho, digits=4))

         try(tau <- tauCOP(cop=cop, para=para, ...), silent=TRUE)
         if(is.null(tau)) tau <- tauCOP(cop=cop, para=para, brute=TRUE, ...)
         if(verbose) message("     Kendall's Tau = ",round(tau, digits=4))
     }

    global.PQD <- ifelse(3*tau >= rho & rho >= 0 & gin >= 0 & blm >= 0, TRUE, FALSE)
    if(is.null(uv)) {
       return(global.PQD)
    } else {
       U <- uv[,1]; V <- uv[,2]
       uvs  <- U*V
       cops <- cop(U, V, para=para, ...)
       deltaCops <- (cops - uvs)
       local.PQD <- deltaCops >= 0
       local.NQD <- deltaCops <= 0
       zeroDelta.Copula.to.Pi <- local.PQD == local.NQD
       zz <- list(global.PQD=global.PQD, Delta.Copula.to.Pi=deltaCops,
                  zeroDelta.Copula.to.Pi=zeroDelta.Copula.to.Pi,
                  local.PQD=local.PQD, local.NQD=local.NQD,
                  beta=blm, gamma=gin, rho=rho, tau=tau)
       return(zz)
    }
}
