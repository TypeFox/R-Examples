modpartltraj <- function(tr, limod)
{
    if (!inherits(tr,"ltraj"))
        stop("tr should be of class ltraj")
    if (length(tr) > 1)
        stop("only implemented for only one trajectory")
    a <- tr[[1]]

    re <- do.call("cbind", lapply(limod, function(mod) {
        ex <- parse(text = mod)
        return(eval(ex, envir=a))
    }))
    re[re<=.Machine$double.xmin] <- .Machine$double.xmin
    cons <- apply(re,1,function(x) all(!is.na(x)))
    re <- re[cons,]
    re <- t(apply(re,1,function(x) x/sum(x)))
    colnames(re) <- paste("mod", 1:length(limod), sep=".")
    indiceNA <- c(1:nrow(a))[cons]
    attr(re, "nna.places") <- indiceNA
    class(re) <- "modpartltraj"
    return(re)
}





print.modpartltraj <- function(x, ...)
{
    if (!inherits(x, "modpartltraj"))
        stop("non convenient x")
    cat("***************************************\n")
    cat("* Probabilities computed for a trajectory \n")
    cat("* with the following models: \n\n")
    print(colnames(x))
}



bestpartmod <- function(mods, Km=30, plotit=TRUE,
                        correction = TRUE, nrep=100)
{
    if (!inherits(mods,"modpartltraj"))
        stop("mods should be of class modpartltraj")
    if (Km > nrow(mods))
        stop("too large number of segments required")

    indiceNA <- attr(mods, "nna.places")


    toto <- .C("optcutr", as.double(t(as.matrix(mods))),
               double(Km), as.integer(Km),
               as.integer(nrow(mods)), as.integer(ncol(mods)),
               double(Km-1), PACKAGE = "adehabitatLT")
    mk <- toto[[2]]
    names(mk) <- as.character(1:length(mk))
    if (correction) {
        yyb <- list()
        for (r in 1:nrep) {
            modb <- mods[sample(1:nrow(mods)),]
            toto <- .C("optcutr", as.double(t(as.matrix(modb))),
                       double(Km), as.integer(Km),
                       as.integer(nrow(mods)), as.integer(ncol(mods)),
                       double(Km-1), PACKAGE = "adehabitatLT")
            yyb[[r]] <- toto[[2]]
        }
        hh <- lapply(yyb, function(o) mk-o)
        yyc2 <- do.call("rbind",hh)
        colnames(yyc2) <- as.character(1:length(mk))
        med <- apply(yyc2, 2, median, na.rm=TRUE)

        if (plotit) {
            fac <- gl(ncol(yyc2), nrow(yyc2))
            plot(c(yyc2)~fac, range=0)
            lines(1:Km,med, col="black", lwd=2)
            abline(v=which.max(med), col="grey")
        }
        cat("Maximum likelihood for K = ", which.max(med), "\n")
        invisible(list(mk = mk, correction = yyc2))
    } else {
        names(mk) <- as.character(1:length(mk))
        if (plotit) {
            plot(1:length(mk), mk,
                 xlab="Number of partitions", ylab="log(Likelihood)", ty="l")
            abline(v=which.max(mk), col="grey")
        }
        cat("Maximum likelihood for K = ", which.max(mk), "\n")
        invisible(list(mk = mk, correction = "none"))
    }
}



partmod.ltraj <- function(tr, npart, mods, na.manage=c("prop.move","locf"))
{
    if (!inherits(tr, "ltraj"))
        stop("tr should be of class \"ltraj\"")
    if (length(tr)>1)
        stop("only one trajectory can be passed")
    if (!inherits(mods,"modpartltraj"))
        stop("mods should be of class modpartltraj")
    na.manage <- match.arg(na.manage)


    cor <- tr[[1]]
    indiceNA <- attr(mods, "nna.places")

    if (npart > nrow(mods))
        stop("too large number of segments required")

    toto <- .C("partrajr", as.double(t(as.matrix(mods))),
               double(npart), integer(npart), integer(npart+1),
               as.integer(nrow(mods)), as.integer(ncol(mods)),
               as.integer(npart), PACKAGE="adehabitatLT")

    curloc <- rev(toto[[4]])
    curloc[2:length(curloc)] <- curloc[2:length(curloc)]+1 ## Where the
                                        # segments begin
    curmod <- rev(toto[[3]])
    curma <- rev(toto[[2]])

    ## Cuts the burst:
    filo <- curloc[-length(curloc)]
    lalo <- curloc[-1]
    lalo[length(lalo)] <- nrow(cor)
    resltr <- lapply(1:length(lalo), function(i) {
        if (i ==1) {
            xyt <- cor[1:indiceNA[lalo[i]],c("x","y","date")]
        } else {
            if (i == length(lalo)) {
                xyt <- cor[indiceNA[filo[i]]:nrow(cor),c("x","y","date")]
            } else {
                xyt <- cor[indiceNA[filo[i]]:indiceNA[lalo[i]],
                           c("x","y","date")]
            }
        }
        return(as.ltraj(xyt[,c("x","y")], xyt[,c("date")],
                        id = id(tr), burst = i))
    })


    ## function cseq to split a segment into runs
    cseq <- function(x)
    {
        id <- diff(c(1,c(1:length(x))[abs(c(0,diff(x)))>0], length(x)+1))
        split(x,unlist(sapply(1:length(id), function(i) rep(i, id[i]))))
    }


    ## Attributes randomly the missing values between two segments
    ## to these segments according to the proportion of missing values
    ## observed in each model
    if (na.manage=="prop.move") {
        ## First count the "internal" and "external" missing values
        ## (within a segment and at a border)
        nadf <- do.call("rbind", lapply(1:length(resltr), function(i) {
            nas <- is.na(resltr[[i]][[1]]$dist[-nrow(resltr[[i]][[1]])])
            vec <- cseq(nas)
            beg <- sum(vec[[length(vec)]])
            intern <- sum(nas) - beg
            return(c(beg,intern))
        }))
        nadf <- as.data.frame(nadf)

        ## proportion of internal missing values for each model
        typmod <- tapply(nadf[,2], factor(curmod), sum)
        typmod <- typmod/sum(typmod)


        ## attribution of the relocations
        for (i in 2:length(resltr)) {
            gg <- resltr[[i-1]][[1]]
            gg2 <- resltr[[i]][[1]]
            gg <- gg[-nrow(gg),]
            ff <- cseq(is.na(gg$dist))
            nna <- sum(ff[[length(ff)]])
            if (nna>1) {
                prot <- sum(typmod[names(typmod)%in%c(curmod[i - 1], curmod[i])])
                if (prot > 1e-08) {
                    nna1 <- floor(nna * typmod[names(typmod)==curmod[i - 1]]/prot)
                } else {
                    nna1 <- floor(nna/2)
                }
                nna2 <- nna-nna1
                gg2 <- rbind(gg[(nrow(gg)-nna2):nrow(gg),], gg2)
                gg <- gg[1:(nrow(gg)-nna2),]

                resltr[[i-1]] <- as.ltraj(gg[,c("x","y")], gg[,c("date")],
                                          id = id(tr),
                                          burst = i-1)
                resltr[[i]] <- as.ltraj(gg2[,c("x","y")], gg2[,c("date")],
                                        id = id(tr),
                                        burst = i)

            }
        }
    }

    resltr <- do.call("c.ltraj", resltr)

    resu <- list(ltraj=resltr,
                 stats=list(locs=curloc, Mk = curma, mod=curmod,
                 which.mod=colnames(mods)[curmod]))
    attr(resu, "nna.places") <- indiceNA
    class(resu) <- "partltraj"
    return(resu)
}


print.partltraj <- function(x, ...)
{
    if (!inherits(x,"partltraj"))
        stop("non convenient x")
    cat(paste("Number of partitions:", length(x$ltraj)))
    cat("\nPartition structure:\n")

    indiceNA <- attr(x, "nna.places")
    loc <- indiceNA[x$stat$loc]
    loc[1] <- 1
    loc[length(loc)] <- sum(unlist(lapply(x$ltraj, nrow))) -
        length(x$ltraj) +1

    rel <- c(c(unlist(t(cbind(loc,
                              rep("|", length(loc)))))))
    rel <- rel[-length(rel)]
    modn <- c(c(unlist(t(cbind(rep("---", length(x$stat$mod)),
                               x$stat$mod
                               )))), "---")
    modn2 <- c(c(unlist(t(cbind(rep("-------", length(x$stat$mod)),
                                x$stat$which.mod
                                )))), "-------")

    print(data.frame(relocation=rel,Num=modn,Model=modn2))
    cat("\nThe segments are contained in the component $ltraj of the list\n\n")
}




plot.partltraj <- function(x, col, addpoints=TRUE, lwd=2,...)
{
    if (!inherits(x, "partltraj"))
        stop("non convenient class of object")

    mod <- x$stat$mod
    if (missing(col)) {
        rg <- max(mod)-1
        col <- grey((1:max(mod))/(max(mod)+rg/5))
    }

    xy <- do.call("rbind", lapply(x$ltraj,
                                  function(i) data.frame(i$x,i$y)))
    plot(xy, asp=1, xlab="x",ylab="y", ty="n")
    colt <- lapply(1:length(mod), function(i) {
        rep(col[mod[i]], length(x$ltraj[[i]]$x))
    })
    lapply(1:length(x$ltraj), function(g) {
        if (addpoints)
            points(x$ltraj[[g]][,c("x","y")], pch=16, col=colt[[g]])
        lines(x$ltraj[[g]][,c("x","y")], pch=16, col=colt[[g]], lwd=lwd)
    })
    points(xy[1,], pch=16, cex=2, col="blue")
}
