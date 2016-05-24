goodmankruskal_gamma <- function(P,...) {
    nr <- nrow(P); nc <- ncol(P)
    Pconc <- 0
    for (i in seq_len(nr-1)) {
        h <- seq(i+1,nr)
        for (j in seq_len(nc-1)) {
                k <- seq(j+1,nc)
                Pconc <- Pconc+2*P[i,j]*sum(P[h,k])
            }
    }
    Pdisc <- 0
    for (i in seq_len(nr-1)) {
        h <- seq(i+1,nr)
        for (j in (seq_len(nc-1)+1)) {
            k <- seq(1,j-1)
            Pdisc <- Pdisc+2*P[i,j]*sum(P[h,k])
        }
    }
    list(C=Pconc,D=Pdisc,gamma=(Pconc-Pdisc)/(Pconc+Pdisc))
}


##' @export
gkgamma <- function(x,data=parent.frame(),strata=NULL,all=FALSE,iid=TRUE,...) {
    if (inherits(x,"formula")) {
        xf <- getoutcome(x,sep="|")
        xx <- attr(xf,"x")
        if (length(xx)==0) stop("Not a valid formula")
        yx <- update(as.formula(paste0(xf,"~.")),xx[[1]])
        if (length(xx)>1) {
            strata <- interaction(model.frame(xx[[2]],data=data))
            x <- yx
        } else {
            x <- model.frame(yx,data=data)
        }
    }
    if (!is.null(strata)) {
        dd <- split(data,strata)
        gam <- lapply(dd,function(d,...) gkgamma(x,data=d,...),
                      ...,
                      iid=TRUE,
                      keep=1:2)
        mgam <- Reduce(function(x,y,...) merge(x,y,...),gam)
        res <- estimate(mgam,function(p,...) {
            k <- length(p)/2
            cd <- lapply(seq(k),function(x) p[(1:2)+2*(x-1)])
            dif <- unlist(lapply(cd,function(x) x[1]-x[2]))
            tot <- unlist(lapply(cd,function(x) x[1]+x[2]))
            gam <- dif/tot
            ##w <- tot/sum(tot)
            ##pgam <- sum(w*gam)            
            c(gam,pgamma=sum(dif)/sum(tot))
        },labels=c(paste0("\u03b3:",names(dd)),"pgamma"),
        iid=iid)
        k <- length(dd)
        if (!iid) {
            for (i in seq_along(gam))
                gam[[i]][c("iid","id")] <- NULL
        }
        homtest <- estimate(res,lava::contr(seq(k),k+1),iid=FALSE)
        return(structure(list(cl=match.call(),k=k,n=unlist(lapply(dd,nrow)),
                              strata=gam,gamma=res,homtest=homtest),
                         class="gkgamma"))
    }
    if (is.table(x) || is.data.frame(x) || is.matrix(x)) {
        x <- multinomial(x)
    }
    if (!inherits(x,"multinomial")) stop("Expected table, data.frame or multinomial object")    
    estimate(x,function(p) {
        P <- x$position; P[] <- p[x$position]
        goodmankruskal_gamma(P)
    },iid=iid,...)
}

##' @export
print.gkgamma <- function(x,...) {
    for (i in seq_along(x$strata)) {
        cat(names(x$strata)[i]," (n=",x$n[i],"):\n",sep="")
        e <- x$strata[[i]]
        print(e)
    }
    printline(50)
    cat("\nGamma coefficient:\n\n")
    print(x$gamma)
    printline(50)
    cat("\nHomogeneity test:\n")
    with(x$homtest$compare,
         cat("chisq = ",statistic,
             ", df = ",parameter,
             ", p-value = ",p.value,"\n"))
}




