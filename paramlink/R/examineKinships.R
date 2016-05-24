examineKinships = function(x, who='all', interfam=c('founders', 'none', 'all'), makeplot=T, ...) {
    
    list_of_linkdats = !inherits(x, 'linkdat') && is.list(x) && all(sapply(x, inherits, 'linkdat'))
    if(list_of_linkdats && any(duplicated(sapply(x, function(xx) xx$famid))))
        stop("When input is a list of linkdat objects, they must have different family ID (famid)")
    
    rel.groups = c('parents', 'siblings', 'grandparents', 'cousins', 'distant', 'unrelated')
    if(identical(who, 'all')) {
        groups = rel.groups
        if(list_of_linkdats) groups = groups[groups != 'distant']
    }    
    else if(identical(who, 'close'))
        groups = rel.groups[1:4]
    else 
        groups = sapply(who, match.arg, rel.groups[-5])

    
    if(makeplot)  {
        #oldpar = par(no.readonly = TRUE)
        if(!list_of_linkdats) {
            layout(rbind(1:2), widths=c(.4,.6))
            plot(x, av=T, cex=.8, title="", margins=c(1,1,3.1,1))
        }
        .plot_IBDtriangle(0.6, mar=c(4.1,4.1,3.1,.1))
        leg_txt = c('Parent-offspring','Siblings', 'Halfsib/Uncle/Grand', '1st cousins', 'Distant', 'Unrelated')[rel.groups %in% groups]
        leg_col = c(2,4,3,6,"cyan","darkgray")[rel.groups %in% groups]
        legend("topright", title=" According to pedigree:", title.adj=0, legend=leg_txt, col=leg_col, pch=16)
    }
    
    dat = .ibdPrep(x)
    p = s = g = co = d = u = NULL
    
    if ('parents' %in% groups) {
        P = related.pairs(x, 'parents', available=T)
        p = IBDestimate(P, dat=dat, plot.action=makeplot, pointcol=2, ...)
    }
    if ('cousins' %in% groups) {
        C = rbind(related.pairs(x, 'cousins', half=F, available=T), 
                  related.pairs(x, 'grandparents', degree=3, available=T), 
                  related.pairs(x, 'nephews_nieces', half=T, available=T), 
                  related.pairs(x, 'nephews_nieces', removal=2, half=F, available=T))
        co = IBDestimate(C, dat=dat, plot.action=makeplot, pointcol=6, ...)
    }
    if ('grandparents' %in% groups) {
        G = rbind(related.pairs(x, 'siblings', half=T, available=T), 
              related.pairs(x, 'grandparents', available=T), 
              related.pairs(x, 'nephews_nieces', half=F, available=T))
        g = IBDestimate(G, dat=dat, plot.action=makeplot, pointcol=3, ...)
    }          
    if ('siblings' %in% groups) {
        S = related.pairs(x, 'siblings', half=F, available=T)
        s = IBDestimate(S, dat=dat, plot.action=makeplot, pointcol=4, ...)
    }
    if ('unrelated' %in% groups) {
        interfam = match.arg(interfam)
        U = related.pairs(x, 'unrelated', available=T, interfam=interfam)
        u = IBDestimate(U, dat=dat, plot.action=makeplot, pointcol="darkgray", ...)
    }
    if ('distant' %in% groups) {
        rest = t(combn(x$available, 2))
        taken = rbind(P,U,S,G,C)
        D = rest[!(rest[,1]*1000 + rest[,2]) %in% (taken[,1]*1000 + taken[,2]),, drop=F]
        d = IBDestimate(D, dat=dat, plot.action=makeplot, pointcol="cyan", ...)
    }
    #if(makeplot) par(oldpar)
    res = list(parents=p, siblings=s, grandparents_uncles_halfsibs=g, cousins=co, distant=d, unrelated=u)
    res = res[!sapply(res, is.null)]
    invisible(res)
}

.IBDsuspects = function(x, radius) { #x output from examineKinships, radius=tolerated distance
    
    extract_distant = function(df, x0, y0, rad) {
        if(is.null(df)) return(NULL)
        subs = subset(df, sqrt((df$k0-x0)^2+(df$k2-y0)^2) > rad)
        if(nrow(subs)>0) subs else NULL
    }
    
    susp = list('suspects-parents'                       = extract_distant(x$parents, 0, 0, radius), 
                'suspects-siblings'                      = extract_distant(x$siblings, 0.25, 0.25, radius),
                'suspects-halfsibs/grandparents/uncles'  = extract_distant(x$grandparents, 0.5, 0, radius), 
                'suspects-cousins/halfuncles/granduncles'= extract_distant(x$cousins, 0.75, 0, radius),
                'suspects-distant'                       = extract_distant(x$distant, 0.9, 0, radius),
                'suspects-unrelated'                     = extract_distant(x$unrelated, 1, 0, radius))
    susp = susp[!sapply(susp, is.null)]
    if(length(susp)==0) susp = NULL
}


IBDestimate = function(g1, g2=NULL, dat, f=NULL, error=NULL, plot.action=2, pointcol=2, cex=1, ...) {
    if(inherits(dat, 'linkdat') || (is.list(dat) && all(sapply(dat, inherits, 'linkdat')))) dat = .ibdPrep(dat)
    if(is.null(f)) f = attr(dat, 'frqs')
    nams = names(dat)
    if(is.matrix(g1)) {
        pairs = match(g1, nams)
        dim(pairs) = dim(g1)
    }
    else {
        g1 = unlist(lapply(g1, grep, nams))
        g2 = if(is.null(g2)) g1 else unlist(lapply(g2, grep, nams))
        pairs = cbind(rep(g1, each=length(g2)), rep(g2, length(g1)))
        pairs = pairs[pairs[,1] != pairs[,2], , drop=F]
        pairs[pairs[,1] > pairs[,2],] <- pairs[pairs[,1] > pairs[,2], 2:1] #sort
        pairs = unique(pairs)
    }
    if(nrow(pairs)==0) return()
    
    arrays = .make_arrays(f, extended=!is.null(error))
    if(plot.action==2) .plot_IBDtriangle()
    
    res = apply(pairs, 1, function(p) {
        k = round(.ibd_estim(p[1], p[2], dat, 'ml', arrays=arrays, error=error, ...), 4)
        if(plot.action > 0) points(k[1], k[2], pch=16, col=pointcol, cex=cex)
        k
    })
    data.frame(ID1 = nams[pairs[,1]], ID2 = nams[pairs[,2]], N = res[3,], k0 = res[1,], k1 = 1-res[1,]-res[2,], k2 = res[2,], row.names=NULL, stringsAsFactors=F)
}



.ibd_estim = function(x, y, dat, f, method=c('ml', 'discrete'), error=NULL, errormodel=1, arrays=NULL, start=c(0.99,0.001), constr=T, tol=1e-8, verbose=F) {
    stopifnot(is.numeric(x) && length(x)==1, is.numeric(y) && length(y)==1) 
    G1 = dat[[x]]; G2 = dat[[y]]
    nonmissing = G1>=0 & G2>=0
    if(!all(nonmissing)) {
        G1 = G1[nonmissing]; G2=G2[nonmissing]
        dat = dat[nonmissing,]
        f = f[nonmissing]
        if(!is.null(arrays))
            arrays = lapply(arrays, function(a) 
                switch(length(dim(a))+1, a[nonmissing], NULL, a[, nonmissing, drop=F], a[, , nonmissing, drop=F]))
    }
    
    if(is.null(arrays))  arrays = .make_arrays(f, extended=!is.null(error))
    
    if(verbose) {
        cat(nrow(dat), "nonmissing genotypes.\nJoint genotype distribution:\n")
        print(tt<-table(G1,G2,dnn=list(x,y)))
        cat("Mean frequencies in each cell:\n")
        for(a in 0:2) for(b in 0:2) tt[a+1, b+1] = mean(arrays$FRQ[G1==a & G2==b])
        print(round(tt,3))
    }
        
    len = length(G1)
    if(is.numeric(error)) {
        err_array = .errorModel(error, errormodel, verbose=verbose)
        
        ind = cbind(1:3, rep(G1+1,each=9), rep(1:3,each=3), rep(G2+1,each=9))
        E = err_array[ind]
        loglik = function(k) #k0=k[1]; k2=k[2]; k1 = 1-k0-k2
            sum(log(.colSums(E*(k[1]*arrays$UNREL + (1-k[1]-k[2])*arrays$PARKID + k[2]*arrays$MZ), 9, len)))
    }
    else {
        A = c(c(G1,G2)+1, seq_len(len))
        dim(A) = c(len,3)
        pG1 = arrays$GT_FREQ[A[, c(1,3)]]
        pG2 = arrays$GT_FREQ[A[, c(2,3)]]
        kidpar = arrays$TRANSMIT[A]
        mz = G1==G2
        loglik = function(k) sum(log(pG1*(k[1]*pG2 + (1-k[1]-k[2])*kidpar + k[2]*mz))) # k0=k[1]; k2=k[2]; k1 = 1-k0-k2
    }    
    if(constr) constraints=list(ineqA=matrix(c(1,0,-1,0,1,-1),3,2), ineqB=c(0,0,1)) else constraints = NULL
    k = switch(match.arg(method),
        discrete = {
                REL_LIST = list('U'=c(1,0),'Parent'=c(0,0),'MZ'=c(0,1),'HS/Gr/Unc'=c(.5,0),'Sib'=c(.25,.25),'FC'=c(.75,0),'SC'=c(15/16,0))
                REL_LIST[[which.max(sapply(REL_LIST, loglik))]]},
        ml = {
            ML = maxLik(loglik, start=start, constraints=constraints, tol=tol)
            ML$estimate}
        )  
    k = c(k, len); names(k) = c('k0', 'k2', 'N')
    if(verbose) cat("Estimate: k = (", paste(round(k[1:2],3),collapse=", "), ')\n', sep="")
    k    
}


.plot_IBDtriangle = function(cex_text=1, mar) {
    # plots the empty triangle with some precomputed points.
    REL_LIST = list('U'=c(1,0),'Parent'=c(0,0),'MZ'=c(0,1),'HS/Gr/Unc'=c(.5,0),'Sib'=c(.25,.25),'FC'=c(.75,0),'SC'=c(15/16,0))
    if(!missing(mar)) par(mar=mar)
    plot(NULL, xlim=c(-.1,1.1), ylim=c(-.1,1.1), xlab="k0", ylab="k2", )        
    fixed = do.call(rbind,REL_LIST)
    points(fixed, pch=16, lwd=2)
    text(fixed, labels=names(REL_LIST), pos=c(1,1,3,1,3,1,1), cex=cex_text)
    abline(v=0);abline(h=0);abline(a=1,b=-1)
    kk0 = seq(0,1,length=100); kk2 = 1+kk0-2*sqrt(kk0)
    points(kk0,kk2,type='l', lty=2)
}  
 

.make_arrays = function(f, extended) {
    # precomputes arrays that are used repeatedly in the likelihood computation
    gf = rbind((1-f)^2, 2*f*(1-f), f^2)
    
    transmit = array(0, dim=c(3,3,length(f)))
    transmit[1,1,] = transmit[3,2,] = 1-f
    transmit[1,2,] = transmit[3,3,] = f
    transmit[2,2,] = 0.5; 
    transmit[2,1,] = 0.5*(1-f); transmit[2,3,] = 0.5*f
    
    if(!extended) 
        return(list(FRQ = f, GT_FREQ = gf, TRANSMIT=transmit))
    
    unrel_arr = parkid_arr = mz_arr = array(0, dim=c(3,3,length(f)))
    for(i in 1:3) {
        for(j in 1:3) {
            unrel_arr[i,j,] = gf[i,] * gf[j,]
            parkid_arr[i,j,] = transmit[i,j, ]*gf[i,]
        }
        mz_arr[i,i,] = gf[i,]
    }        
    list(FRQ=f, UNREL=unrel_arr, PARKID=parkid_arr, MZ=mz_arr)
}


.errorModel = function(e, errormodel=1, verbose=T) {
    ERR = switch(errormodel,
            c(1-e-e^2, e, e^2, e, 1-2*e, e, e^2, e, 1-e-e^2),
            c(1-2*e, e, e, e, 1-2*e, e, e, e, 1-2*e),
            c(1-e, e, 0, e, 1-2*e, e, 0, e, 1-e))
    dim(ERR) = c(3,3)
    if(verbose) {dimnames(ERR) = rep(list(c('AA','AB','BB')),2); cat('Error model:\n'); print(ERR)}
    return(ERR %o% ERR)
}


# SKETCH
# .rel_resamp = function(G1, G2, f, samp=NULL, method=c('ml', 'discrete', 'both'), start=c(0.99,0.001), plot=T) {
    # method = match.arg(method)
    # if(method=='both') {
        # par(mfrow=c(1,2))
        # rel_resamp(G1, G2, f, samp=samp, method='ml', start=start, plot=plot)
        # rel_resamp(G1, G2, f, samp=samp, method='dis', start=start, plot=plot)
        # return(NULL)
    # }
    # if(plot) {
        # .plot_IBDtriangle()
        # fixed = do.call(rbind,REL_LIST2)
    # }
    # res = sapply(1:samp[1], function(i) {
        # sampl = sample.int(length(f), samp[2], replace=F)
        # ibd_estim(G1[sampl],G2[sampl],f=f[sampl], method=method, start=start)})
    
    # switch(method,
    # discrete = {
        # pr = sapply(1:length(REL_LIST), function(i) mean(res==i)*100)
        # av = REL_LIST2[[which.max(pr)]]
        # conf = max(pr)
        # if(plot)  text(fixed[pr>0,,drop=F], sprintf("%.0f%%", pr[pr>0]), pos=4-c(1,1,3,1,3,1,1)[pr>0], cex=.7)
    # }, ml = {
        # av = rowMeans(res)
        # conf = apply(res,1,quantile,c(.05,.95))
        # if(plot) points(res[1,], res[2,], cex=.9, col="gray")
    # })
    
    # if(plot) {
        # points(av[1],av[2],col=2,pch=4,cex=2,lwd=3)
        # points(fixed, pch=16, lwd=2)
    # }
    # list(av=c(k0=av[1], k1=1-sum(av), k2=av[2]), conf=conf)
# }


.ibdPrep = function(x, removebad=T) {
    ### if x is a list of linkdats
    if(!inherits(x, 'linkdat') && is.list(x) && all(sapply(x,inherits, 'linkdat'))) {
        alldata = lapply(x, .ibdPrep, removebad = F)
        alldata = lapply(seq_along(x), function(i) {dd = alldata[[i]]; names(dd) = paste(x[[i]]$famid, names(dd), sep="-"); dd})
        dat = do.call(cbind, alldata)
        frqs = do.call(cbind, lapply(alldata, attr, 'frqs'))
        if(removebad) {
            chrom = vapply(x[[1]]$markerdata, function(m) attr(m, 'chrom'), numeric(1))
            badmarkers = (!is.na(chrom) & chrom > 22) | apply(frqs, 1, function(fr) any(fr==0 | fr==1))  
            # ad hoc way of removing markers that are (mistakenly) annotated as uniallelic in *some* patients
            dat = dat[!badmarkers,]; frqs = frqs[!badmarkers, ]
        }
        attr(dat, 'frqs') = rowMeans(frqs)
        return(dat)
    }
    
    ### a single linkdat object
    frqs = 1 - unlist(lapply(x$markerdata, function(m) attr(m, 'afreq')[1]))  #NB: Samsvarer med eksom-format, der frekv oppgis for ALT.
    if(removebad) {
        annot = vapply(x$markerdata, function(m) c(attr(m, 'nalleles'), attr(m, 'chrom')), numeric(2))
        badmarkers = annot[1,]!=2 | (!is.na(annot[2,]) & annot[2,] > 22) | frqs==0 | frqs==1 
        x = removeMarkers(x, which(badmarkers))
        frqs = frqs[!badmarkers]    
    }
     
    M = as.matrix(x, inc=F)[, -(1:6), drop=F]
    even = 2*(seq_len(x$nMark))
    MM = paste(M[, even-1], M[, even], sep="") #TODO: SLOW!
    dat = matrix(-1, nrow=x$nInd, ncol=x$nMark)
    #dat[MM=='00' | MM=='01' | MM=='10'] = -1 
    dat[MM=='11'] = 0; dat[MM=='12' | MM=='21'] = 1; dat[MM=='22'] = 2
    dat = as.data.frame(t.default(dat))
    names(dat) = x$orig.ids
    attr(dat, 'frqs') = frqs
    dat
}



.setSNPfreqs = function(x, newfreqs) {
    stopifnot(all(vapply(x$markerdata, function(m) attr(m, 'nalleles'), numeric(1))==2))
    newfreqs = rep(newfreqs, length=x$nMark)
    for(i in seq_len(x$nMark)) attr(x$markerdata[[i]], 'afreq') = c(newfreqs[i], 1-newfreqs[i])
    x
}
