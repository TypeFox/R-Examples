.packageName <- "CCP"



"WilksLambda" <-
function(rho,p,q) 
{
    minpq <- min(p, q)
    WilksLambda <- numeric(minpq)
    for (rhostart in 1:minpq) {
        WilksLambda[rhostart] <- prod((1 - rho^2)[rhostart:minpq])
    }
    invisible(WilksLambda)
}







"HotellingLawleyTrace" <-
function (rho, p, q) 
{
    minpq <- min(p, q)   
    HotellingLawleyTrace <- numeric(minpq)
    for (rhostart in 1:minpq) {
        rhosq = rho^2
        HotellingLawleyTrace[rhostart] = sum((rhosq/(1-rhosq))[rhostart:minpq])        
    }
    invisible(HotellingLawleyTrace)
}




"PillaiBartlettTrace" <-
function (rho, p, q) 
{
    minpq <- min(p, q)   
    PillaiBartlettTrace <- numeric(minpq)
    for (rhostart in 1:minpq) {
        PillaiBartlettTrace[rhostart] = sum((rho^2)[rhostart:minpq])        
    }
    invisible(PillaiBartlettTrace)
}





"RaoF.stat" <-
function (rho, N, p, q) 
{
    minpq <- min(p, q)   
    WilksLambda <- WilksLambda(rho,p,q)
    
    RaoF <- numeric(minpq)
    df1  <- numeric(minpq)
    df2  <- numeric(minpq)
    de = N - 1.5 - (p + q)/2   
    for (rhostart in 1:minpq) {
        k = rhostart - 1  # 0,1,2
        df1[rhostart] = (p-k)*(q-k)
        if ((p-k == 1)||(q-k == 1))
            nu = 1
        else         
            nu = sqrt((df1[rhostart]^2-4)/((p-k)^2+(q-k)^2-5)) 
        df2[rhostart] = de*nu - df1[rhostart]/2 + 1
        nu1 = 1/nu
        w = WilksLambda[rhostart]^nu1
        RaoF[rhostart] = df2[rhostart]/df1[rhostart]*(1-w)/w        
    }
    invisible(list(stat = WilksLambda, approx = RaoF, df1 = df1, df2 = df2))
}



"Hotelling.stat" <-
function (rho, N, p, q) 
{
    minpq <- min(p, q)
    HotellingLawleyTrace <- HotellingLawleyTrace(rho, p, q)
   
    Hotelling <- numeric(minpq)
    df1       <- numeric(minpq)
    df2       <- numeric(minpq)
    for (rhostart in 1:minpq) {
        k = rhostart - 1
        df1[rhostart]  = (p - k)*(q - k)
        df2[rhostart]  = minpq*(N - 2 - p - q + 2*k) + 2
        Hotelling[rhostart] = HotellingLawleyTrace[rhostart]/minpq/df1[rhostart] * df2[rhostart]
    }
    invisible(list(stat = HotellingLawleyTrace, approx = Hotelling, df1 = df1, df2 = df2))
}



"Pillai.stat" <-
function (rho, N, p, q) 
{
    minpq <- min(p, q)
    PillaiBartlettTrace <- PillaiBartlettTrace(rho, p, q)
   
    Pillai <- numeric(minpq)
    df1    <- numeric(minpq)
    df2    <- numeric(minpq)
    for (rhostart in 1:minpq) {
        k = rhostart - 1
        df1[rhostart]  = (p - k)*(q - k)
        df2[rhostart]  = minpq*(N - 1 + minpq - p - q + 2*k) 
        Pillai[rhostart] = PillaiBartlettTrace[rhostart] / df1[rhostart] * df2[rhostart] / ( minpq - PillaiBartlettTrace[rhostart] )
    }
    invisible(list(stat = PillaiBartlettTrace, approx = Pillai, df1 = df1, df2 = df2))
}





"p.Roy" <-
function (rho, N, p, q) 
{
    stat  <- rho[1]^2  # Roy
    df1  <- q
    df2  <- N - 1 - q
    approx   <- stat / df1 * df2 / ( 1 - stat )

    p.value  <- 1 - pf(approx, df1, df2)

    invisible(list(id="Roy", stat = stat, approx = approx, df1 = df1,  df2 = df2, p.value = p.value))
}





"p.asym" <-
function (rho, N, p, q, tstat = "Wilks") 
{
    minpq  <- min(p, q)

    if (length(rho) != minpq) 
        stop(" Function p.asym: Improper length of vector containing the canonical correlations.\n")

    if (tstat == "Wilks")
        {
        out  <- RaoF.stat(rho, N, p, q)
        head <- paste("Wilks' Lambda, using F-approximation (Rao's F):\n") 
	  }           
    else if ( tstat == "Hotelling" )
        {
        out <- Hotelling.stat(rho, N, p, q)
        head <- paste(" Hotelling-Lawley Trace, using F-approximation:\n")
        }        
    else if ( tstat == "Pillai" )
        {
        out <- Pillai.stat(rho, N, p, q)
        head <- paste(" Pillai-Bartlett Trace, using F-approximation:\n")
        } 
    else if ( tstat == "Roy" )
        {
        minpq = 1
        out <- p.Roy(rho, N, p, q)
        head <- paste(" Roy's Largest Root, using F-approximation:\n")
        }         
    else
        stop(" Function p.asym: test statistic must be of type 'Wilks', 'Hotelling', 'Pillai', or 'Roy'.\n")

    stat   <- out$stat
    approx <- out$approx
    df1    <- out$df1
    df2    <- out$df2

    p.value <- numeric(minpq)
    for (rhostart in 1:minpq) {
        p.value[rhostart] <- 1 - pf(approx[rhostart], df1[rhostart], df2[rhostart])
    }

    tab <- cbind(stat, approx, df1, df2, p.value)
    rn = character(length = minpq)
    for (k in 1:minpq) rn[k] = paste(k," to ",minpq,": ",sep = "")
    rownames(tab) = rn
    cat(head)
    print(tab)

    if ( tstat == "Roy" )
        cat("\n F statistic for Roy's Greatest Root is an upper bound.\n")

    invisible(list(id = tstat, stat = stat, approx = approx, df1 = df1, df2 = df2, p.value = p.value))
}






"p.perm" <-
function(X, Y, nboot = 999, rhostart = 1, type = "Wilks")
{
    if( nrow(Y) != nrow(X)) 
        stop(" Function p.perm: Input arrays must not have unequal number of rows.\n") 
    if( (type == "Roy") && (rhostart > 1) ) 
        stop(" Function p.perm: When using Roy's Largest Root, parameter rhostart can only be 1.\n")
 
    N <- nrow(Y)               	
    ind <- 1:N
    p = dim(X)[2]
    q = dim(Y)[2]
    minpq = min(p,q)

   if( rhostart > minpq ) 
        stop(" Function p.perm: Parameter rhostart too big.\n")

    rho0  <- cancor(X,Y[ind,])$cor	   	    
    
    if (type == "Wilks")
        stat0 <- WilksLambda(rho0,p,q)[rhostart]
    else if ( type == "Hotelling" )
        stat0 <- HotellingLawleyTrace(rho0,p,q)[rhostart]
    else if ( type == "Pillai" )
        stat0 <- PillaiBartlettTrace(rho0,p,q)[rhostart]
    else if ( type == "Roy" )
        stat0 <- p.Roy(rho0,N,p,q)$stat
    else
        stop(" Function p.perm: Illegal type of test statistic.\n")

    stat <- numeric(nboot)
    for (i in 1:nboot ) {
        ind.mixed <- sample(ind, size = N, replace = FALSE)  	      
        rho    <- cancor(X,Y[ind.mixed,])$cor                         
        
        if (type == "Wilks")
            stat[i] <- WilksLambda(rho,p,q)[rhostart]
        else if ( type == "Hotelling" )
            stat[i] <- HotellingLawleyTrace(rho,p,q)[rhostart]
        else if ( type == "Pillai" )
            stat[i] <- PillaiBartlettTrace(rho,p,q)[rhostart]
        else if ( type == "Roy" )
            stat[i] <- p.Roy(rho,N,p,q)$stat
    }

    if (type == "Wilks") {
        nexcess <- sum(stat <= stat0)
        p <- mean(stat <= stat0)
        }
    else {    
        nexcess <- sum(stat >= stat0) 
        p <- mean(stat >= stat0)
        }
    
    mstat = mean(stat)

    tab <- cbind(stat0,mstat,nboot,nexcess,p)
    rownames(tab) = ""
    cat(" Permutation resampling using",type,"'s statistic:\n")
    print(tab)    

    invisible(list(id="Permutation", type = type, stat0 = stat0, stat = stat, nexcess = nexcess, p.value = p))
}


 


"plt.asym" <-
function(p.asym.out, rhostart = 1)
{
    if ("id" %in% names(p.asym.out))
        type <- p.asym.out$id
    else
        stop(" Function plt.asym: Use output of 'p.asym' as input for 'plt.asym'.\n")

    if ( rhostart > length(p.asym.out$stat) ) 
        stop(" Function plt.asym: Parameter 'rhostart' too big for this dataset.\n")      
 
    df1 <- p.asym.out$df1[rhostart]
    df2 <- p.asym.out$df2[rhostart]
    stat <- p.asym.out$approx[rhostart]
    pval <- p.asym.out$p.value[rhostart]

    minpq = length(p.asym.out$df1)
    
    if (type == "Wilks")
        head <- paste("F-approximation for Wilks Lambda, rho =", rhostart, "to", minpq)        
    else if ( type == "Hotelling" )
        head <- paste("F-approximation for Hotelling-Lawley Trace, rho=", rhostart, "to", minpq)             
    else if ( type == "Pillai" )
        head <- paste("F-approximation for Pillai-Bartlett Trace, rho=", rhostart, "to", minpq)           
    else if ( type == "Roy" )
        head <- paste("F-approximation for Roys Largest Root, rho=", rhostart, "to", minpq)            

    pv = format(pval, scientific = FALSE, digits = 3)
    dv2 = format(df2, scientific = FALSE, digits = 3)
    statf = format(stat, scientific = FALSE, digits = 3)
    subtext = paste("F=", statf, ", df1=", df1, ", df2=", dv2,", p=", pv)
    x = seq(qf(0.005,df1,df2), qf(0.995,df1,df2), length.out= 200)
    plot(x, df(x,df1,df2), main = head, ylab = "", xlab = "", lwd = 2, type = "l", col = "grey", sub = subtext)
    abline(v = stat, lwd = 3, col = "red", lty=3)
    abline(h = 0, lwd = 2, col = "grey", lty = 3)
}




"plt.perm" <-
function(p.perm.out)
{
    if ("id" %in% names(p.perm.out))
        id = p.perm.out$id
    else 
        stop(" Function plot.perm: Use output of p.perm as input for plot.perm(1).\n")
    if (id != "Permutation") 
        stop(" Function plot.perm: Use output of p.perm as input for plot.perm(2).\n")
    
    stat0 = p.perm.out$stat0
    stat  = p.perm.out$stat
    pval  = p.perm.out$p.value
    nexcess  = p.perm.out$nexcess
    type = p.perm.out$type

    head <- "Permutation distribution"  		
    hist(stat, breaks = 30, freq = FALSE, main = head)
    sform = format(stat0, scientific = FALSE, digits = 3)
    pv = format(pval, scientific = FALSE, digits = 3)
    mtext(paste("test=",type, ", original test statistic=", sform, ", p=",pv), side = 3)  			
    abline(v=mean(stat),lwd=2,col="grey",lty=3)       		 
    abline(v=stat0,lwd=3,col="red",lty=3)  				

}








