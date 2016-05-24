`.alignment` <-
function(targets, actual, align.method="spline"){
    if(align.method == "spline"){
        require(splines)
    }
    ret <- list()
    for(x in 1:dim(actual)[1]){
        acts <- as.numeric(actual[x,!is.na(actual[x,])])
        targs <- as.numeric(targets[!is.na(actual[x,])])
        pts <- length(acts)
        fr <- switch(as.character(pts),
            "0" = data.frame(Y=1:2, X=1:2),
            "1" = data.frame(Y=1:2+targs, X=1:2+acts),
            data.frame(Y=targs, X=acts)
        )
        if(pts <= 2 || align.method == "affine"){
            ret[[x]] <- lm(Y~X, dat=fr)
        } else if (pts == 3 || align.method == "PL") {
            form <- paste("Y ~ X + ", paste("I(pmax(X-", acts[2:(pts-1)], ", 0))",
                collapse=" + ", sep=""), sep="")
            ret[[x]] <- lm(as.formula(form), dat=fr)
        } else {
            ret[[x]] <- interpSpline(acts, targs)        
        }
    }
    names(ret) <- rownames(actual)
    ret
}

`.benj.hoch` <-
function(pvals, FDR=.1){
    all(diff(pvals)>=0) || stop("p-values must be sorted in ascending order")
    max(c(0,which(pvals<=FDR*(1:length(pvals))/length(pvals))))
}

`.biweight.FTICRMS` <-
function(x, K=9, max.iter=20, na.rm=TRUE){
    if(!na.rm && sum(is.na(x))>0){
        return(NA)
    } else {
        x <- x[!is.na(x)]
    }
    iter <- 0
    ind.old <- 0
    ind.new <- 1
    center.bw <- median(x)
    scale.bw <- median(abs(x-center.bw))
    while(!identical(ind.old,ind.new) && iter<max.iter){
        iter <- iter + 1
        ind.old <- ind.new
        u <- (x-center.bw)/scale.bw
        ind.new <- (abs(u) <= K)
        center.bw <- sum(ind.new * x * (1-(u/K)^2)^2)/sum(ind.new * (1-(u/K)^2)^2)
        scale.bw <- median(abs(x-center.bw))
    }
    if(iter == max.iter && !identical(ind.old,ind.new)){
        warning(paste("Iteration limit of", max.iter, "reached without convergence."))
    }
    list(center = center.bw, scale = scale.bw, iter = iter)
}

`.break.clusters` <-
function(clust){
    if(max(table(clust$File))==1){
        list(clust)
    } else {
        cuts <- c()
        tmp <- lapply(levels(clust$File), function(x){which(clust$File==x)})
        tmp <- do.call(rbind,lapply(tmp, function(x){if(length(x)<=1){c()} else
            if(length(x)==2){x} else {cbind(x[-length(x)], x[-1])}}))
        if(is.null(dim(tmp))){
            tmp <- matrix(tmp,nrow=1)
        }
        tmp <- tmp[order(tmp[,1]), , drop=FALSE]
        while(any(diff(tmp[,2])<0)){
            tmp <- tmp[c(diff(tmp[,2]),1)>0, , drop=FALSE]
        }
        while(min(dim(tmp))){
            cuts <- c(cuts, mean(clust$Center_hat[c(tmp[1,2],max(tmp[tmp[,1]<tmp[1,2],1]))]))
            tmp <- tmp[tmp[,1]>=tmp[1,2], , drop=FALSE]
        }
        by(clust,findInterval(clust$Center_hat,cuts),function(x){x})
    }
}

`.break.dups` <-
function(clust.list){
    ret <- unlist(lapply(clust.list,.break.clusters),recursive=FALSE)
    names(ret) <- NULL
    ret
}

`.centered.coeffs` <-
function(M){
    Y <- M[,2]
    X <- M[,1]-mean(M[,1])
    Xb <- cbind(1,X,X^2)
    as.vector(solve(t(Xb) %*% Xb, t(Xb) %*% Y))
}

`.cluster.matrix` <-
function(clust.peaks, filenames=levels(clust.peaks[[1]]$File)){
    tmp <- data.frame(matrix(NA,ncol=length(filenames),nrow=length(clust.peaks)))
    colnames(tmp) <- as.character(filenames)
    rownames(tmp) <- lapply(clust.peaks,function(x){paste(round(range(x$Center_hat),4),
        collapse=",")})

    for(i in 1:length(clust.peaks)){
        tmp[i,as.character(clust.peaks[[i]]$File)] <- clust.peaks[[i]]$Max_hat
    }

    tmp
}

`.cluster.peaks` <-
function(peaks, clust.method="ppm", clust.constant=10){
    tmp <- peaks
    tmp <- tmp[order(tmp$Center_hat),]
    rownames(tmp) <- 1:dim(tmp)[1]
    if(clust.method=="constant"){
        diffs <- clust.constant
    } else if(clust.method=="ppm"){
        diffs <- clust.constant*tmp$Center_hat[-1]/10^6
    } else if(clust.method=="usewidth"){
        mod <- lm(log(Width_hat) ~ log(Center_hat), tmp)
        diffs <- exp(predict(mod)[-dim(tmp)[1]])*clust.constant/exp(predict(mod, 
            data.frame(Center_hat = median(tmp$Center_hat))))
        diffs <- pmax(diffs, clust.constant)
    }
    clust <- cumsum(c(TRUE,diff(tmp$Center_hat) > diffs))
    by(tmp,clust,function(x){x})
}

`.get.sp.masses` <-
function(all.masses, sp.masses, iso.dist, inds=TRUE){
    isos <- sort(as.numeric(outer((0:iso.dist)*1.003, sp.masses, "+")))
    ends <- as.numeric(outer(c(-.02,.02), isos, "+"))
    wh <- which(diff(ends)<0)
    if(length(wh)){
        ends <- ends[-c(wh,wh+1)]
    }
    rowmass <- sapply(all.masses, function(x)
        mean(as.numeric(strsplit(x, ",")[[1]])))
    ret <- which(findInterval(rowmass, ends) %% 2 == 1)
    if(!inds){
        massends <- do.call(rbind, strsplit(all.masses[ret],","))
        massends <- as.numeric(apply(massends,1,as.numeric))
        ends <- matrix(ends,ncol=2,byrow=TRUE)
        ends <- ends[!apply(apply(ends, 2, findInterval, vec=massends) %% 2 == 1, 
            1, any),]
        ret <- apply(ends,1,paste,collapse=",")
    }
    ret
}

`.loc.maxes` <-
function(vec){
    1+which(vec[-c(1,length(vec))]>vec[-c(1,2)] & vec[-c(1,length(vec))] > 
        vec[-c(length(vec),length(vec)-1)])
}

`.loc.mins` <-
function(vec){
    .loc.maxes(-vec)
}

`.peak.parab` <-
function(pts, num.pts, R2.thresh){
    Rsq <- sapply(1:(dim(pts)[1]-num.pts+1), function(x){.Rsquared(pts[x+(0:(num.pts-1)),])})
    if(max(Rsq)>=R2.thresh){
        locs <- which.max(Rsq)+(0:(num.pts-1))
        coeffs <- .centered.coeffs(pts[locs,])
        return(c(mean(pts[locs,1])-coeffs[2]/2/coeffs[3], coeffs[1]-coeffs[2]^2/4/coeffs[3],
            -1/coeffs[3]))
    } else {
        return(rep(0,3))
    }
}

`.Rsquared` <-
function(M){
    Y <- M[,2]
    X <- M[,1]-mean(M[,1])
    Q <- poly(X,2)
    (t(Y) %*% Q %*% t(Q) %*% Y)/(t(Y) %*% (diag(length(Y))-1/length(Y)) %*% Y)
}

