indval <- function(x, ...)
{
    UseMethod("indval")
}

indval.default <- function(x,clustering,numitr=1000, ...)
{
    if (!is.data.frame(x)) x <- data.frame(x)

    clustering <- clustify(clustering)
    levels <- levels(clustering)
    clustering <- as.integer(clustering)

    if (any(apply(x>0,2,sum)==0)) stop('All species must occur in at least one plot')
    numplt <- nrow(x)
    numspc <- ncol(x)
    numcls <- as.integer(length(table(clustering)))
    maxcls <- rep(0,numspc)
    relfrq <- matrix(0,nrow=numspc,ncol=numcls)
    relabu <- matrix(0,nrow=numspc,ncol=numcls)
    indval <- matrix(0,nrow=numspc,ncol=numcls)
    indcls <- rep(0,numspc)
    pval <- rep(0,numspc)
    tmpfrq <- rep(0.0,numcls)
    tmpabu <- rep(0.0,numcls)
    pclass <- rep(0,numplt)
    tclass <- rep(0,numplt)
    errcod <- 0
    tmp <- .Fortran("duleg",
        as.double(as.matrix(x)),
        as.integer(numplt),
        as.integer(numspc),
        as.integer(factor(clustering)),
        as.integer(table(clustering)),
        as.integer(numcls),
        as.integer(numitr),
        relfrq = relfrq,
        relabu = relabu,
        indval = indval,
        pval = pval,
        indcls = indcls,
        maxcls = as.integer(maxcls),
        as.double(tmpfrq),
        as.double(tmpabu),
        as.integer(pclass),
        as.integer(tclass),
        errcod = as.integer(errcod),
        PACKAGE='labdsv')
    out <- list(relfrq=data.frame(tmp$relfrq),relabu=data.frame(tmp$relabu),
              indval=data.frame(tmp$indval),maxcls=tmp$maxcls,indcls=tmp$indcls,
              pval=tmp$pval,error=tmp$errcod)
    row.names(out$relfrq) <- names(x)
    row.names(out$relabu) <- names(x)
    row.names(out$indval) <- names(x)
    names(out$maxcls) <- names(x)
    names(out$indcls) <- names(x)
    names(out$pval) <- names(x)
    names(out$relfrq) <- levels
    names(out$relabu) <- levels
    names(out$indval) <- levels
    class(out) <- 'indval'
    if (out$error == 1) cat('WARNING: one or more species not assigned to any cluster\n')
    out
}

indval.stride <- function(x,taxa,numitr=1,...)
{
    res <- rep(NA,ncol(x$clustering))
    for (i in 1:ncol(x$clustering)) {
        res[i] <- mean(indval(taxa,x$clustering[,i],numitr=numitr)$indcls)
    }
    clusters <- x$seq
    indval <- res
    out <- data.frame(clusters,indval)
    out
}

summary.indval <- function (object, p = 0.05, type='short', digits=2, show=p, sort=FALSE, 
                            too.many = 100, ...) 
{
    if (object$error == 1) cat('WARNING: one or more species not assigned to any cluster\n')
    if (type == 'short') {
        tmp <- data.frame(object$maxcls[object$pval <= p], round(object$indcls[object$pval <= 
            p], 4), object$pval[object$pval <= p])
        names(tmp) <- c("cluster", "indicator_value", "probability")
        if (nrow(tmp) <= too.many) print(tmp[order(tmp$cluster, -tmp$indicator_value), ])
        cat(paste("\nSum of probabilities                 = ", sum(object$pval),
            "\n"))
        cat(paste("\nSum of Indicator Values              = ", round(sum(object$indcls),digits=2),
            "\n"))
        cat(paste("\nSum of Significant Indicator Values  = ", round(sum(tmp$indicator_value),digits=2),
            "\n"))
        cat(paste("\nNumber of Significant Indicators     = ", nrow(tmp),"\n"))
        cat(paste("\nSignificant Indicator Distribution\n"))
        print(table(tmp$cluster))

    } else {
        tmp <- format(round(object$indval,digits=digits))
        keep <- apply(object$indval,1,function(x){max(x)>show})
        tmp <- tmp[keep,] 
        tmp[tmp < show] <- substring(" .  ",1,nchar(tmp[1,1]))
        print(tmp)
    }
    if (sort) {
        repeat {
            plots <- readline(' enter the plots    : ')
            if (plots == "") {
                break
            } else {
                pnt <- readline(' in front of        : ')
            }
            for (i in strsplit(plots,",")[[1]]){
                ord <- 1:nrow(tmp)
                x <- match(i,row.names(tmp))
                print(paste(i,x))
                if (!is.na(x)) {
                    ord <- ord[-x]
                    y <- match(pnt,row.names(tmp[ord,]))
                    print(y)
                    if (!is.na(y)) {
                            if (y==1) {
                                ord <- c(x,ord)
                            } else {
                                first <- ord[1:(y-1)]
                                last <- ord[y:length(ord)]
                                ord <- c(first,x,last)
                            }
                            tmp <- tmp[ord,]
                            print(tmp)
                        } else {
                            print(paste('species',pnt,'does not exist'))
                        }
                    } else {
                        print(paste('species',i,'does not exist'))
                    }
                }
            }
        invisible(tmp)
    }
}
