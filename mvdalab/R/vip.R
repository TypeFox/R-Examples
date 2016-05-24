vip <- function(object, ncomp = object$ncomp, conf = .95) {
  if(object$val.method == "none" | object$val.method == "loo") {
    A <- c(object$y.loadings2)[1:ncomp]^2
    B <- apply(as.matrix(object$weights[, 1:ncomp]), 2, function(x) sum(x^2))
    C <- sweep(as.matrix(object$weights[, 1:ncomp])^2, 2, A / B, "*")
    if(ncomp > 1) {
    VIPa <- data.frame(t(sqrt(nrow(C) * apply(C, 1, cumsum) / cumsum(A))))
    } else {
      VIPa <- data.frame((sqrt(nrow(C) * apply(C, 1, cumsum) / cumsum(A))))
    }
    names(VIPa) <- 1:ncomp
    VIP <- reshape(VIPa, idvar = "variables", ids = row.names(VIPa),
            times = names(VIPa[1:ncomp]), timevar = "ncomp",
            varying = list(1:ncomp), direction = "long")
    names(VIP)[2] <- "VIP"
    row.names(VIP) <- NULL
    Results <- list(VIP = VIP, VIP.s = NULL, vip.boots = NULL, 
                    val.method = object$val.method, ncomp = ncomp)
    class(Results) <- "vip"
    Results
  } else {
    A.all <- c(object$y.loadings2)[1:ncomp]^2
    B.all <- apply(as.matrix(object$weights[, 1:ncomp]), 2, function(x) sum(x^2))
    C.all <- sweep(as.matrix(object$weights[, 1:ncomp])^2, 2, A.all / B.all, "*")
    if(ncomp > 1) {
      VIP.all <- data.frame(t(sqrt(nrow(C.all) * apply(C.all, 1, cumsum) / cumsum(A.all))))
    } else {
      VIP.all <- data.frame((sqrt(nrow(C.all) * apply(C.all, 1, cumsum) / cumsum(A.all))))
    }
    names(VIP.all) <- 1:ncomp
    VIP.All.Data <- reshape(VIP.all, idvar = "variables", ids = row.names(VIP.all),
                   times = names(VIP.all[1:ncomp]), timevar = "ncomp",
                   varying = list(1:ncomp), direction = "long")
    names(VIP.All.Data)[2] <- "VIP"
    row.names(VIP.All.Data) <- NULL
    vip.boots <- plyr::llply(1:object$validation$bootstraps, function(x) {
      A <- object$validation$y.loadings2[x, 1:ncomp]^2
      B <- apply(as.matrix(object$validation$weights[[x]][, 1:ncomp]), 2, function(x) sum(x^2))
      C <- sweep(as.matrix(object$validation$weights[[x]][, 1:ncomp])^2, 2, A / B, "*")
      sqrt(nrow(C) * apply(C, 1, cumsum) / cumsum(A))
    })
    A <- plyr::llply(1:ncomp, function(x) {
      B <- plyr::llply(1:object$validation$bootstraps, function(y) {
        if(ncomp > 1) {
              vip.boots[[y]][x, ]
        } else vip.boots[[y]]
          })
    B <- do.call("rbind", B)
    })
    Upper <- 1 - (((1 - conf)/2))
    Lower <- 1 - Upper
    C <- plyr::llply(1:ncomp, function(x) {
        C.a <- data.frame(t(apply(A[[x]], 2, function(x) quantile(x, c(Lower, Upper), na.rm = T))))
        C.a$ncomp <- x
        C.a$variables <- row.names(C.a)
        C.a$boot.mean <- apply(A[[x]], 2, function(x) mean(x, na.rm = T))
        C.a$'Bootstrap Error' <- apply(A[[x]], 2, function(x) sd(x, na.rm = T))
        C.a$Skewness <- apply(A[[x]], 2, function(x) skewness(x, na.rm = T))
        names(C.a)[1:2] <- paste(c(Lower * 100, Upper * 100), "%", sep = "")
        row.names(C.a) <- NULL
        C.a
    })
    vips.quants <- do.call("rbind", C)
    row.names(vips.quants) <- NULL
    vips.quants$Actual <- VIP.All.Data[, 2]
    vips.quants$Bias <- vips.quants$boot.mean - vips.quants$Actual
    vips.quants$'t value (VIP = 1)' <- (vips.quants$Actual - 1) / vips.quants$'Bootstrap Error'
    vips.quants$'bias corrected t value (VIP = 1)' <- ((vips.quants$Actual - (vips.quants$Bias)) - 1) / vips.quants$'Bootstrap Error'
    vips.quants$'bias t value' <- vips.quants$Bias / vips.quants$'Bootstrap Error'
    VIP.Out <- vips.quants[, c(3, 4, 8, 1:2, 5, 7, 9, 6, 10:11)]
    VIP <- data.frame(stack(vips.quants[, 1:2]), ncomp = vips.quants[, 3], variables = vips.quants[, 4])
    Results <- list(VIP = vips.quants, VIP.s = VIP, vip.boots = A, VIP.All.Data = VIP.All.Data, 
                    val.method = object$val.method, ncomp = ncomp, VIP.Out = VIP.Out)
    class(Results) <- "vip"
    Results
  }
}

print.vip <- function(x, ncomp = x$VIP$ncomp, ...) {
  if(x$val.method == "none" | x$val.method == "loo") {
    V <- x$VIP[x$VIP$ncomp %in% ncomp, ]
    print(V)
  } else {
    V <- x$VIP.Out[x$VIP$ncomp %in% ncomp, ]
    print(V)
  }
}

vip.boots <- function(x, ...) {
  if (is.null(x$vip.boots)) {
    stop("No bootstrapping was done for this model")
  }
  cat("\nVIP Boots:\n")
  print(x$vip.boots)
}





