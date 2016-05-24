lavielle <- function(x, ...)
{
    UseMethod("lavielle")
}

lavielle.default <- function(x, Lmin, Kmax, ld = 1,
                             type=c("mean","var","meanvar"), ...)
{
    if (any(is.na(x)))
        stop("No missing values allowed in the series")
    series2 <- scale(x)
    if ((Kmax)*Lmin>length(x))
        stop(paste("Please be careful:\n With a maximum of Kmax =", Kmax,
                   "segments in the series\n and a minimum segment length",
                   "Lmin = ", Lmin, "observations,\n the series should",
                   "contain at least Kmax * Lmin = ", Kmax*Lmin, "observations.\n",
                   "It actually contains", length(x), "observations.\n",
                   "Either decrease Kmax or Lmin"))
    type <- match.arg(type)
    Kmax <- Kmax+1
    type <- which(type==c("mean","var","meanvar"))
    res <- .Call("contrastM", series2, Lmin, type, ld, PACKAGE="adehabitatLT")
    res2 <- .Call("dynprog", res, Kmax, PACKAGE="adehabitatLT")
    K <- Kmax-1
    outp <- list(contmat=res,
                 sumcont=res2[[1]], matpath=res2[[2]], Kmax=Kmax-1, Lmin = Lmin,
                 ld=ld,
                 series=x)
    class(outp) <- "lavielle"
    attr(outp, "typeseg") <- "default"
    return(invisible(outp))
}

print.lavielle <- function(x, ...)
{
    cat("Segmentation of a series using the method of Lavielle\n",
        "$contmat: contrast matrix\n",
        "$sumcont: optimal contrast\n",
        "$matpath: matrix of the path\n",
        "$Kmax:    maximum number of segments\n",
        "$Lmin:    Minimum number of obs. to build a segment\n",
        "$ld:      The size of the subsampling grid\n",
        "$series:  the series\n\n")
}

chooseseg <- function(lav, S=0.75, output=c("full","opt"),
                      draw=TRUE)
{
    if (!inherits(lav,"lavielle"))
        stop("lav should be of class lavielle")
    output <- match.arg(output)
    res2 <- lav$sumcont
    Kmax <- lav$Kmax
    Jk <- res2[1:Kmax,ncol(res2)]
    Jkt <- ((Jk[Kmax]-Jk)/(Jk[Kmax]-Jk[1]))*(Kmax-1)+1
    D <- c(Inf, sapply(2:(length(Jkt)-1), function(K) Jkt[K-1] - 2*Jkt[K] + Jkt[K+1]))
    df <- data.frame(K=1:(Kmax-1), D=D, Jk = Jk[-Kmax])
    cons <- c(as.numeric(Jk[-1]<Jk[-length(Jk)]))
    if (draw) {
        plot(1:(Kmax-1), df$Jk, ty="b",
             xlab="K", ylab="Jk")
    }
    if (length(df$K[df$D>S&cons==1])>=1) {
        Kopt <- max(df$K[df$D>S&cons==1])
    } else {
        Kopt <- 1
    }
    if (output=="full")
        return(df)
    return(Kopt)
}

findpath <- function(lav, K, plotit=TRUE)
{
    if (!inherits(lav,"lavielle"))
        stop("lav should be of class lavielle")
    ## Find path
    if (attr(lav, "typeseg") == "default") {
        series <- lav$series
        pat <- .Call("findpath", lav$matpath, K, as.integer(lav$Kmax), PACKAGE="adehabitatLT")
        ## conversion according to the subsampling grid ld:
        pat <- (pat*lav$ld) - lav$ld + 1
        if (plotit) {
            plot(series, ty="l")
            abline(v=pat, col="red")
        }
        pat <- rev(pat)
        pat <- lapply(2:(K+1), function(i) {
            if (i != (K+1)) {
                return(c(pat[i-1], pat[i]-1))
            } else {
                return(c(pat[i-1], length(lav$series)))
            }
        })
        return(pat)
    } else {
        nna.places <- attr(lav, "nna.places")
        ltraj <- attr(lav, "ltraj.ori")
        x <- ld(ltraj)
        attr(lav, "typeseg") <- "default"
        lipla <- findpath(lav, K, plotit=plotit)
        lipla <- lapply(1:length(lipla), function(i) {
            c(nna.places[lipla[[i]][1]], nna.places[lipla[[i]][2]])
        })
        if (lipla[[1]][1]!=1)
            lipla[[1]][1] <- 1
        if (lipla[[length(lipla)]][2]!=nrow(x))
            lipla[[length(lipla)]][2] <- nrow(x)
        for (i in 2:length(lipla)) {
            if (lipla[[i]][1]!=(lipla[[i-1]][2]+1)){
                lipla[[i-1]][2] <- lipla[[i]][1]-1
            }
        }
        res <- do.call("c.ltraj", lapply(1:length(lipla), function(i) {
            y <- lipla[[i]]
            ff <- dl(x[y[1]:y[2],])
            burst(ff) <- paste("Segment", i, sep=".")
            return(ff)
        }))

        return(res)
    }
}


lavielle.ltraj <- function(x, Lmin, Kmax, ld = 1, which="dist",
                           type=c("mean","var","meanvar"), ...)
{
    if (!inherits(x, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (length(x)>1)
        stop("only implemented for one trajectory")
    y <- x[[1]]
    if (!is.null(attr(y, "infolocs")))
        y <- cbind(y, attr(y, "infolocs"))
    ex <- parse(text = which)
    coin <- eval(ex, envir = y)
    cons <- !is.na(coin)
    nna.places <- c(1:length(coin))[cons]
    coin <- coin[nna.places]
    lav <- lavielle(x=coin, type=type, Lmin=Lmin, Kmax=Kmax, ld = ld)
    attr(lav, "nna.places") <- nna.places
    attr(lav, "ltraj.ori") <- x
    attr(lav, "typeseg") <- "ltraj"
    return(lav)
}

