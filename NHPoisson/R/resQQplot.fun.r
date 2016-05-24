resQQplot.fun<-
function (nsim, objres, covariates, clevel = 0.95, cores = 1, 
     tit = "", fixed.seed=NULL, histWgraph=TRUE) 
{
    typeI <- objres$typeI
    typeRes <- objres$ScaRes$typeRes
    fittedlambda <- objres$fittedlambda
    emplambda <- objres$emplambda
    lambda <- objres$mlePP@lambdafit
    beta <- as.list(objres$mlePP@coef)
    tind <- objres$mlePP@tind
    if (is.null(tit) | (tit == "")) 
        tit <- objres$mlePP@tit
    t <- objres$mlePP@t
    h <- objres$res
    lint <- objres$lint
    if (is.null(t)) 
        t <- c(1:length(lambda))
    if (typeRes == "Raw") 
        res <- objres$RawRes
    else res <- objres$ScaRes$ScaRes
    cl <- makeCluster(cores)
  if (!is.null(fixed.seed)) clusterSetRNGStream(cl=cl,iseed=fixed.seed)
    clusterExport(cl, c("resSim.fun", "simNHP.fun", "buscar", 
        "fitPP.fun", "CalcResD.fun", "CalcRes.fun"))
    writeLines("Calculating...", sep = " ")
    Mres <- parSapply(cl, c(1:nsim), FUN = resSim.fun, lambda = lambda, 
        covariates = covariates, beta = beta, lint = lint, t = t, 
        tind = tind, typeI = typeI, typeRes = typeRes, h = h)
    stopCluster(cl)

    resinf <- apply(Mres, FUN = quantile, MARGIN = 1, p = 1 - 
        clevel, na.rm = TRUE)
    ressup <- apply(Mres, FUN = quantile, MARGIN = 1, p = clevel, 
        na.rm = TRUE)
    resmed <- apply(Mres, FUN = mean, MARGIN = 1, na.rm = TRUE)
    aux <- sort.list(resmed)
    scay <- c(min(resmed, ressup, resinf, na.rm = T), max(resmed, 
        ressup, resinf, na.rm = T))
    writeLines("Done")
    marca <- as.numeric((res <= ressup) & (res >= resinf))

    if ((.Platform$OS.type=="windows")&(histWgraph==TRUE))
	dev.new(record=TRUE)
    par(mfrow=c(1,1))
    plot(resmed, res, ylab = "empirical quantile of the residuals", 
        ylim = scay, xlab = "mean quantile of simulated residuals", 
        type = "n")
    points(resmed[marca == 1], res[marca == 1], pch = 16, cex = 0.6)
    points(resmed[marca == 0], res[marca == 0], pch = 16, cex = 0.6, 
        col = "red")
    lines(resmed[aux], resmed[aux])
    lines(resmed[aux], ressup[aux], col = "blue")
    lines(resmed[aux], resinf[aux], col = "blue")
    return(list(resmed = resmed, ressup = ressup, resinf = resinf, 
        nsim = nsim, objres = objres,fixed.seed=fixed.seed))
}
