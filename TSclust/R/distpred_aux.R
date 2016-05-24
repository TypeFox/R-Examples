
###########################################################################
# Nonparametric approach to generate bootstrap predictions: 
# Option 1 - Autoregressive bootstrap (AB), if metAB = 1
# Option 2 - Conditional bootstrap (CB), if metCB = 1
###########################################################################

pred.AB.CB <- function( serie, tipo.b.NW, l.grid, ancho.CV, cota, n.rep.boot, l.sobrantes, factor.segunda.b, l.horizon, metAB, metCB)
    
    # Auxiliar functions required:
    # h_cv, nadaraya.watson.mod, bootstrap.res in AUXILIARES1 file.
    # dpill, dpik in KernSmooth package.
    
{
    prediction.boot <- matrix(0, l.horizon, n.rep.boot)
    prediction.boot.cond <- matrix(0, l.horizon, n.rep.boot)
    
    if ( (metAB == 1) |  (metCB == 1) ) 
    {
        long <- length(serie)
        nube <- cbind(serie[1:(long-1)],serie[2:long])
        
        # STEP 1. 
        #########
        
        # Computing the bandwidth in accordance with the value of tipo.b.NW
        
        if (tipo.b.NW == 1) 		# cross-validation
        {
            LG <- seq(abs(max(nube[,1])-min(nube[,1]))/1000, abs(max(nube[,1])-min(nube[,1]))/4, length.out=l.grid)
            g_m <- h_cv(nube,ancho.CV,LG)					
        }
        
        
        recursive.trycatch.bw <- function( trimval) {
            if (trimval > 0.3) stop(paste("Coud not find appropiate bandwith for a series,", trimval))
            g_m <- NULL
            tryCatch( {  #catch common dpill problem and solve it be increasing trim
                g_m <- dpill(nube[,1],nube[,2], trim=trimval) 
            }, error = function(e) {
                gm <- recursive.trycatch.bw(trimval+0.01) 
            })
            g_m
        }
        
        if (tipo.b.NW == 2) {
            g_m <- -1 
            gm <- recursive.trycatch.bw(0.0)
            #tryCatch( {  #catch common dpill problem an solve it be increasing trim
            #    g_m <- dpill(nube[,1],nube[,2]) 
            #}, error = function(e) {
            #    g_m <- dpill(nube[,1],nube[,2], trim=0.02) 
            #})            
        }
        
        h <- factor.segunda.b * g_m 
        
        # Obtaining the truncated Nadaraya-Watson estimator for the autoregression
        
        n.w <- nadaraya.watson.mod(nube,nube[,1],g_m,kind.of.kernel = "dnorm",cota)
        
        # STEP 2.  Computing nonparametric residuals and estimating their density
        #########
        
        eps <- nube[,2] - n.w
        eps.centred <- eps - mean(eps)
        g_f <- dpik(eps.centred, scalest="minim", level=2, kernel="normal", gridsize=401)										
        
        # STEP 3. 	Drawing bootstrap resamples from the nonparametric estimator 
        ######### 	of the residuals density with bandwidth g_f
        
        eps.boot <- bootstrap.res(eps.centred, g_f, kind.of.kernel="dnorm", n.rep.boot, length(eps.centred)+l.sobrantes)				
    }
    
    # Hereafter, steps 4 and 5 depend on wheteher AB or CB are used
    
    if (metAB == 1) 
    {
        
        # STEP 4 (AB)	Constructing the bootstrap series based on the 				#			autoregression estimated in STEP 1 and the bootstrap 
        #############	residuals resamples obtained in STEP 3  
        
        resample.boot <- matrix(0, nrow(eps.boot), ncol(eps.boot))
        resample.boot[1,] <- eps.boot[1,]
        for (i in 2:nrow(eps.boot) )  
            resample.boot[i,] <- nadaraya.watson.mod(nube, resample.boot[i-1,] ,g_m, kind.of.kernel = "dnorm", cota) +eps.boot[i,]			
        resample.boot <- resample.boot[l.sobrantes:nrow(eps.boot),]
        
        # STEP 5 (AB) 	Drawing new bootstrap resamples from the nonparametric 
        #			estimator of the residuals density (with g_f) to obtain 		#############	the bootstrap prediction-paths.
        
        eps.boot.pred <- bootstrap.res(eps.centred, g_f, kind.of.kernel="dnorm", n.rep.boot, l.horizon)								
        
        # Computing the AB-bootstrap prediction-paths
        
        resample.boot.datos <- array(0, dim=c(long-1,2,n.rep.boot))
        for (i in 1:n.rep.boot) 
            resample.boot.datos[,,i] <- cbind(resample.boot[1:(long-1),i],resample.boot[2:long,i])
        
        for (j in 1:n.rep.boot)		
            prediction.boot[1,j] <- nadaraya.watson.mod(resample.boot.datos[,,j], serie[long], h, kind.of.kernel = "dnorm", cota) + eps.boot.pred[1,j] 	
        #
        if (l.horizon > 1)
            for (i in 2:l.horizon) 
                for (j in 1:n.rep.boot)		
                    prediction.boot[i,j] <- nadaraya.watson.mod(resample.boot.datos[,,j], prediction.boot[(i-1),j], h, kind.of.kernel = "dnorm",cota) + eps.boot.pred[i,j]						
        
    }   
    
    if (metCB == 1) 
    {
        
        # STEP 4 (CB) 	Drawing new bootstrap resamples from the nonparametric 
        #			estimator of the residuals density (with g_f) to obtain 		#############	the bootstrap prediction-paths.
        
        eps.boot.pred.cond <- bootstrap.res(eps.centred, g_f, kind.of.kernel = "dnorm", n.rep.boot, l.horizon)						
        
        start <- nadaraya.watson.mod(nube,serie[long], h, kind.of.kernel = "dnorm", cota)
        
        # Computing the CB-bootstrap prediction-paths
        
        prediction.boot.cond[1,] <- start + eps.boot.pred.cond[1,]
        if (l.horizon > 1) 
            for (i in 2:l.horizon ) 
                prediction.boot.cond[i,] <- nadaraya.watson.mod(nube, prediction.boot.cond[(i-1),], h ,kind.of.kernel = "dnorm",cota ) + eps.boot.pred.cond[i,]								
        
    }   
    
    list( AB = prediction.boot, CB = prediction.boot.cond)
}


###########################################################################
# Parametric approach to generate bootstrap predictions: Sieve bootstrap 
###########################################################################

pred.SB <- function( serie, n.rep.boot, l.sobrantes, l.horizon, metSB)
    
    # Auxiliar functions required:
    # bootstrap.res in AUXILIARES1 file.
    # dpik in KernSmooth package.
    
{
    prediction.boot.ar <- matrix(0, l.horizon, n.rep.boot)
    if (metSB == 1) 
    {
        long <- length(serie)
        # An AR(1) model is assumed so that AIC is not actually used
        
        # STEP 1.   AR(1) fit, estimated AR coeff. and estimated residuals var.
        #########
        
        yt.sim.ar <- ar.ols(serie, aic=FALSE, demean=FALSE, order.max=1)	
        phi.ar <- yt.sim.ar$ar[1,1,1]		# estimated autoregressive coeff.
        eps.ar <- yt.sim.ar$resid[2:long]	# residuals vector (length: long-1)
        var.innov.ar <- yt.sim.ar$var.pred[1,1]	# estimated innovations variance
        
        # STEP 2. 	Nonparametric estimator of residuals density
        #########
        
        eps.centred.ar <- eps.ar - mean(eps.ar)
        g_f.ar <- dpik(eps.centred.ar, scalest="minim", level=2, kernel="normal", gridsize=401)
        
        # STEP 3. 	Drawing bootstrap resamples from the nonparametric estimator 
        ######### 	of the residuals density with bandwidth g_f.ar
        
        eps.boot.ar <- bootstrap.res(eps.centred.ar, g_f.ar, kind.of.kernel = "dnorm", n.rep.boot, long+l.sobrantes+1)
        
        # STEP 4. 	Define recursively the bootstrap series and obtain the new 
        #########	estimated AR(1) coeff. 
        
        resample.boot.ar <- matrix(0, nrow(eps.boot.ar), ncol(eps.boot.ar))
        resample.boot.ar[1,] <- eps.boot.ar[1,]
        for ( i in 2:(long+l.sobrantes+1) )  
            resample.boot.ar[i,] <- phi.ar * resample.boot.ar[(i-1),] + eps.boot.ar[i,]
        resample.boot.ar <- resample.boot.ar[(l.sobrantes+2):nrow(eps.boot.ar),]
        phi.boot.ar <- array(0,dim=c(1,1,ncol(eps.boot.ar))) 
        for ( i in 1:ncol(eps.boot.ar) )  
            phi.boot.ar[,,i] <- ar.ols(resample.boot.ar[,i], aic=FALSE,order.max=1,demean=FALSE)$ar[1,1,1]	# an AR(1) is fitted 
        
        # STEP 5. 	The bootstrap prediction-paths are constructed by using new 
        #		bootstrap resamples of residuals (from the density estimated 	#########	with g_f.ar) and the AR(1) models based on phi.boot.ar coeffs.
        
        eps.boot.pred.ar <- bootstrap.res(eps.centred.ar, g_f.ar, kind.of.kernel = "dnorm", n.rep.boot, l.horizon)
        
        prediction.boot.ar[1,] <- phi.boot.ar * serie[long] + eps.boot.pred.ar[1,]
        if (l.horizon > 1)	
            for ( i in 2:l.horizon )  
                prediction.boot.ar[i,] <- phi.boot.ar * prediction.boot.ar[i-1,] + eps.boot.pred.ar[i,]	
    }
    
    SB <- prediction.boot.ar
    return(SB) 
}





backtransf.TRAMO <- function (t.Predicciones, M.L.Series, Logaritmos, Diferencias, 
                              Medias) 
{
    lts <- dim(M.L.Series)[1]
    nhorizon <- dim(t.Predicciones)[1]
    nboot <- dim(t.Predicciones)[2]
    nmetodos <- 1
    ns <- dim(t.Predicciones)[3]
    Pred <- t.Predicciones
    dimnames(Pred) <- dimnames(t.Predicciones)
    vector.xi <- vector("list", ns)
    for (j4 in 1:ns) for (j3 in 1:nmetodos) for (j2 in 1:nboot) {
        if (Diferencias[j4] > 0) {
            vector.xi[[j4]] <- M.L.Series[(lts - Diferencias[j4] + 
                                             1):lts, j4]
            aux <- diffinv(t.Predicciones[, j2,  j4], differences = Diferencias[j4], 
                           xi = vector.xi[[j4]])
            Pred[, j2, j4] <- aux[-(1:Diferencias[j4])]
        }
        Pred[, j2, j4] <- Pred[, j2,  j4] + Medias[j4]
        if (Logaritmos[j4] == 0) 
            Pred[, j2, j4] <- exp(Pred[, j2, j4])
    }
    #    Centres <- array(0, dim = c(nhorizon, ns))
    #   dimnames(Centres)[[2]] <- dimnames(t.Predicciones)[[3]]
    #   for (j3 in 1:nmetodos) for (j4 in 1:ns) {
    #       aux <- Pred[, , j4]
    #       Centres[, j4 ] <- apply(aux, 1, mean)
    #   }
    Pred
    #   list(Predicciones = Pred, CentrosPred = Centres)
}


transf.TRAMO <- function (Series, Logaritmos, Diferencias) 
{
    Series <- as.matrix(Series)
    ns <- dim(Series)[2]
    ls <- dim(Series)[1]
    l.Series <- Series
    for (i in 1:ns) if (Logaritmos[i] == 0)
        l.Series[, i] <- log(Series[, i])
    m.l.Series <- l.Series
    medias <- numeric(length = ns)
    for (i in 1:ns) {
        medias[i] <- mean(l.Series[, i])
        m.l.Series[, i] <- l.Series[, i] - medias[i]
    }
    MD <- max(Diferencias)
    t.Series <- array(0, dim = c(ls - MD, ns))
    for (i in 1:ns) if (Diferencias[i] > 0) {
        aux <- diff(m.l.Series[, i], dif = Diferencias[i])
        linf <- length(aux) - (ls - MD) + 1
        t.Series[, i] <- aux[linf:length(aux)]
    }
    else {
        t.Series[, i] <- m.l.Series[(1 + MD):ls, i]
    }
    colnames(t.Series) <- colnames(Series)
    list(T.Ser = t.Series, M.L.Series=m.l.Series, Medias.log.series = medias)
}



