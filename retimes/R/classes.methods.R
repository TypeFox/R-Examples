setClass("timefit",
    representation(
        x = "vector",
        samples = "matrix",
        par = "vector",
        bootPar = "matrix",
        bootValue = "vector",
        sigmaValid = "vector",
        start = "vector",
        logLik = "numeric",
        AIC = "numeric",
        BIC = "numeric"
    )
)

setMethod("show","timefit",
    function(object)
    {
        cat("    mu:",sprintf("%.4f",object@par[1]),"\n")
        cat(" sigma:",sprintf("%.4f",object@par[2]),"\n")
        cat("   tau:",sprintf("%.4f",object@par[3]),"\n")
        cat("---\n")
        cat(" logLik:",sprintf("%.4f",object@logLik)," ")
        cat("AIC:",sprintf("%.4f",object@AIC),"\n")
    }
)

setMethod("plot","timefit",
    function(object,x="timefit",y="timefit")
    {
        D <- density(object@x)
        lim <- list(x=range(D$x),y=range(D$y))
        lim$y[2] <- lim$y[2]+lim$y[2]/2.5
        plot(D,xlim=lim$x,ylim=lim$y,main="",xlab="",ylab="",lty=1)
        par(new=TRUE)
        curve(dexgauss(x,mu=object@par[1],sigma=object@par[2],tau=object@par[3]),
            xlim=lim$x,ylim=lim$y, main="Distribution",xlab="Observed data",ylab="Density",lty=2)
    }
)

setMethod("logLik","timefit",
    function(object)
        return(object@logLik)
)

setMethod("AIC","timefit",
    function(object)
        return(object@AIC)
)

setMethod("BIC","timefit",
    function(object)
        return(object@BIC)
)
