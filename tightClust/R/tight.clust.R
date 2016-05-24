tight.clust <-
function(x, target, k.min, alpha=.1, beta=.6, top.can=7, seq.num=2, resamp.num=10, samp.p=.7, nstart=1, remain.p=.1, k.stop=5, standardize.gene=TRUE, random.seed=NULL) {

    find.candidates<-function(x, k, alpha=.1, top.can=7, resamp.num=10, samp.p=.7, nstart=1) {

        kmeans.classify<-function(x, k, p=.7, ...) {
            size<-dim(x)[1]
            partial.size<-round(p*size)
            partial<-sample(rownames(x), partial.size)
            centers<-kmeans(x[partial,], k, ...)$centers
            return(suppressWarnings(kmeans(x,centers,iter.max=1,algorithm="Lloyd"))$cluster)
        }

        find.candidates.one<-function(x) {
            tmp<-apply(x==1,1,sum)
            return(which(x[,which(tmp==max(tmp))[1]]==1))
        }

        extend.candidate<-function(D, can, alpha=.1) {
            can.ex<-which(apply(as.matrix(D[,can]>=1-alpha),1,all))
            D.temp<-D[can.ex,can.ex]
            if(!is.matrix(D.temp)) {
                D.temp<-as.matrix(D.temp)
                colnames(D.temp)<-names(can.ex)
            }
            D.bad<-apply(as.matrix(D.temp<1-alpha),1,sum)
            while(sum(D.bad)>0) {
                index<-which(D.bad==max(D.bad))[1]
                D.temp<-D.temp[-index,-index]
                D.bad<-apply(as.matrix(D.temp<1-alpha),1,sum)
            }
            return(can.ex[colnames(D.temp)])
        }

        N<-dim(x)[1]
        Dbar<-matrix(0,N,N)
        for(i in 1:resamp.num) {
            cl<-kmeans.classify(x, k, p=samp.p, nstart=nstart, iter.max=100)
            D<-outer(cl,cl,function(a,b) a==b)
            Dbar=Dbar+D
        }
        Dbar=Dbar/resamp.num
        colnames(Dbar)<-1:N
        rownames(Dbar)<-1:N
        i=1
        D.temp<-Dbar
        res<-list()
        while(i<=top.can*2 && dim(D.temp)[1]>0) {
            candidate.one<-find.candidates.one(D.temp)
            candidate<-extend.candidate(D.temp, candidate.one, alpha=alpha)
            D.temp<-D.temp[-candidate,-candidate]
            res[[i]]<-names(candidate)
            mode(res[[i]])<-"numeric"
            i=i+1
        }
        res<-res[order(unlist(lapply(res,length)),decreasing=TRUE)][1:top.can]
        return(res)
    }

    if(!is.null(random.seed)) set.seed(random.seed)
    original.data<-x
    if(standardize.gene) x<-t(scale(t(x)))
    k.max<-k.min+10
    id<-rownames(x)
    N<-dim(x)[1]
    write(paste("Number of points:",N,"\tDimension:",dim(x)[2],"\n"),"")
    rownames(x)<-1:N
    index.m<-as.matrix(expand.grid(lapply(1:seq.num, function(x) 1:top.can)))
    remain<-N
    nfound<-0
    found<-TRUE
    k0<-k.min
    k<-k0
    candidates<-list()
    tclust<-list()
    while(nfound<target && remain/N>=remain.p && (found || k<=k.max)) {
        if(found) {
            write(paste("Looking for tight cluster",nfound+1,"..."),"")
            k<-k0
            for(i in 1:seq.num) {
                write(paste("k =",k+i-1),"")
                candidates[[i]]<-find.candidates(x,k+i-1,alpha=alpha,top.can=top.can,resamp.num=resamp.num,samp.p=samp.p,nstart=nstart)
            }
        } else {
            candidates<-candidates[-1]
            candidates[[seq.num]]<-find.candidates(x,k+seq.num-1,alpha=alpha,top.can=top.can,resamp.num=resamp.num,samp.p=samp.p,nstart=nstart)
        }

        calc.beta<-function(y) {
            temp<-lapply(1:seq.num, function(z) candidates[[z]][[y[z]]])
            i.temp<-temp[[1]]
            u.temp<-temp[[i]]
            for(j in 2:seq.num) {
                i.temp<-intersect(i.temp,temp[[j]])
                u.temp<-union(u.temp,temp[[j]])
            }
            return(length(i.temp)/length(u.temp))
        }

        beta.temp<-unlist(apply(index.m, 1, calc.beta))
        if(any(beta.temp>=beta)) {
            found=TRUE
            nfound=nfound+1
            write(paste(nfound, "tight cluster(s) found!"),"")
            if(k0>k.stop) k0=k0-1
            found.temp<-candidates[[seq.num]][[index.m[which(beta.temp>=beta)[1],seq.num]]]
            tclust[[nfound]]<-rownames(x)[found.temp]
            mode(tclust[[nfound]])<-"numeric"
            x<-x[-found.temp,]
            remain<-remain-length(tclust[[nfound]])
            write(paste("Cluster size:",length(tclust[[nfound]]),"\tRemaining number of points:",remain,"\n"),"")
        } else {
            found=FALSE
            k=k+1
        }
    }
    clust.id<-rep(-1,N)
    size<-unlist(lapply(tclust,length))
    for(i in 1:length(tclust)) clust.id[tclust[[i]]]<-i
    res<-list(data=original.data,cluster=clust.id,size=size)
    class(res)<-"tight.clust"
    return(res)
}
