# .packageName<-'MixGHD'


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                        Coaleased                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################




###################################### Main ##########################################




MainMCGHD=function(data=NULL, gpar0=NULL, G=2, max.iter=100, eps=1e-2,  label=NULL, method="km",nr=NULL){
    pcol=ncol(data)
    
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgparC(data=data, g=G, w=z,l=lc)
        
    }
    else{
        if (is.null(gpar0)) gpar  = rmgpar(g=G,p=ncol(data),data=data, method=method,nr=nr)
        else{ gpar = gpar0
            for(i in 1:G){
                gpar[[i]]$gam   = eigen( gpar0[[i]]$sigma)$vectors
                gpar[[i]]$phi   = eigen( gpar0[[i]]$sigma)$values}
        }
        
    }
    
    
    
    loglik = numeric(max.iter)
    for (i in 1:3) {
        gpar = EMgrstep(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llik(data, gpar)
    }
    
    while ( ( getall(loglik[1:i]) > eps) & (i < (max.iter) ) )  {
        i = i+1
        gpar = EMgrstep(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llik(data, gpar)
        
    }
    if(i<max.iter){loglik=loglik[-(i+1:max.iter)]}
    #BIC=2*loglik[max.iter]-log(nrow(data))*(2*(G-1)+G*(2*pcol+0.5*pcol*(pcol-1))+G*2*pcol+G*2)=loglik[i]
    BIC=2*loglik[i]-log(nrow(data))*((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2))
    AIC=2*loglik[i]-2*((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2))
    AIC3=2*loglik[i]-3*((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2))
    z=weights(data=data, gpar= gpar)
    ICL=BIC+sum(log(apply(z,1,max)))
    par=partrue(gpar,G)
    val = list(loglik= loglik[1:i], gpar=gpar,par=par, z=z, map=MAP(data=data, gpar= gpar, label=label),BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3 )
    return(val)
    
}



EMgrstep <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL, vg=NULL,it=NULL) {
    
    G= length(gpar$pi);
    
    if (is.null(w)) w = weights(data=data, gpar=gpar,v=v)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    #print("****")
    for (k in 1:G ) {
        vg = iweights(data=data, gpar=gpar[[k]],v=v)
        gpar[[k]] = updatemaScplp(y=data, par=gpar[[k]], weights=w[,k], iweights=vg, alpha.known=NULL, v=v, it=it)
        meanvg = apply( vg, 2, weighted.mean, w=w[,k] )
        if (meanvg[2] < 0.01 ) meanvg = c(.99, .01)
        gpar[[k]]$wg = meanvg
        #		print(gpar[[k]]$wg)
    }
    gpar$pi = apply(w,2,mean)
    #	print(gpar$pi)
    return(gpar)
}

#####Stopping criteria

getall <- function(loglik) {
    if (length(loglik) <3) stop("must have at least 3 likelihood values")
    n = length(loglik)
    lm1 = loglik[n]
    lm  = loglik[(n-1)]
    lm_1  = loglik[(n-2)]
    am = (lm1 - lm)/(lm - lm_1)
    lm1.Inf = lm + (lm1 - lm)/(1-am)
    val = lm1.Inf - lm
    if (is.nan(val)) val=0
    if (val < 0) val= 1
    return( val )
}

################ rotate the parameters


partrue<-function(gpar,G=2){
    gparT=gpar
    grpraT=list()
    for(i in 1:G){
        par=gpar[[i]]
        sort=sort.int(par$phi,decreasing=T,index.return=T)
        par$phi=sort$x
        par$gam=par$gam[,sort$ix]
        
        gam=gpar[[i]]$gam
        
        par$mu=par$mu%*%t(gam)
        par$alpha=par$alpha%*%t(gam)
        #par$Sigma=gam%*%diag(par$sigma)%*%t(gam)
        
        gparT[[i]]=par
    }
    return(gparT)
    
}




##################################density functions###################################

#####Mixture of GH and MSGH


ddmsghyp <- function(data=NULL, par=NULL) {
    if (par$wg[2] != 0) {
        val1 = ddghyp(y=data,par=par,log=TRUE)
        val2 = dmsghyp(y=data,par=par,log=TRUE)
        val =  par$wg[1]*exp(val1) + par$wg[2]*exp(val2)
    } else val = ddghyp(y=data,par=par,log=FALSE)
    
    return( val )
}



##Multivariate GH rotated

ddghyp <- function(y=NULL, par=NULL, log=FALSE, invS=NULL) {
    # x  is a n * p matrix
    x     = y %*% (par$gam)
    mu = par$mu
    
    alpha = par$alpha
    
    if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
    d = length(mu)
    omega =  par$cpl0[1];
    lambda = par$cpl0[2];
    
    pa = omega + sum(alpha^2/par$phi )
    xmu = sweep(x,2,mu)
    mx = omega + apply( xmu^2, 1, weighted.sum, wt=1/par$phi)
    kx = sqrt(mx*pa)
    
    lvx = matrix(0, nrow=nrow(x), 4)
    lvx[,1] = (lambda - d/2)*log(kx)
    lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
    lvx[,3] = apply(xmu, 1, weighted.sum,  wt=alpha/par$phi )
    
    lv = numeric(6)
    lv[1] = -1/2*sum( log( par$phi) ) -d/2*(log(2)+log(pi))
    lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
    lv[3] = -lambda/2*( log(1) )
    lv[4] = lambda*log(omega) * 0
    lv[5] = (d/2 - lambda)*log( pa )
    
    val = apply(lvx,1,sum) + sum(lv)
    if (!log) val = exp( val )
    
    return(val)
}


### MS GH



dmsghyp <- function(y, par, log=FALSE) {
    # x is a n x p matrix
    x     = y %*% (par$gam)
    mu    = par$mu; phi = par$phi;
    alpha = par$alpha;
    d = length(mu); chi = par$cpl[,1]; psi = par$cpl[,1];
    lambda = par$cpl[,2];
    
    xmu = sweep(x,2,mu,"-")
    
    pa = psi + alpha^2/phi   # p numeric
    mx = sweep(sweep(xmu^2,2,1/phi, FUN="*"), 2,chi, "+") # n x p matrix
    kx = sqrt(sweep(mx, 2, pa, "*")) # nxp matrix
    
    lx1 = sweep( sweep(log(mx),2,log(pa),"-"), 2, (lambda - 1/2)/2, "*")
    lx2 = t(apply(kx, 1, function(z,lam=NULL) { log(besselK( z, nu=lambda-1/2, expon.scaled =TRUE)) - z }, lam=lambda ))
    lx3 = sweep(xmu, 2, alpha/phi, FUN="*")
    
    lv  = matrix(0, nrow=d, ncol=3)
    lv[,1] = - 1/2*log( phi ) - 1/2*(log(2)+log(pi)) # d=1
    lv[,2] = - (log(besselK( sqrt(chi*psi), nu=lambda, expon.scaled =TRUE)) - sqrt(chi*psi) )
    lv[,3] = lambda/2*(log(psi)-log(chi) )
    
    if(ncol(y)==1){lx2=t(lx2)}
    val = apply(lx1 + lx2 + lx3, 1, sum) + sum(lv)
    
    if (!log) val = exp( val )
    
    return(val)
}





##################################### Log likelihood ##################################
llik <- function(data=NULL, gpar=NULL, label=NULL) {
    fx = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:length(gpar$pi) ) fx[,k] = ddmsghyp(data=data, par=gpar[[k]])
    
    lablik = matrix(1, nrow=nrow(data), ncol=length(gpar$pi))
    if (!is.null(label)) lablik = combinewk(lablik, label= label)
    
    val = apply(fx*lablik, 1, weighted.sum, wt=gpar$pi)
    val = sum(log(val))
    
    return(val)
}




####################################labeling##########################################

MAP <- function(data, gpar, label=NULL) {
    w = weights(data=data, gpar=gpar, v=1)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
    z = as.numeric(z)
    return( z)
}


######################################### E step #######################################

##############zig  model weights

combinewk <- function(weights=NULL, label=NULL)	{
    # known is a numeric with
    # 0 if unknown group membership
    # 1,2,3,.. for label of known group
    if (is.null(label)) stop('label is null')
    kw     = label !=0
    for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
    return(weights)
}



weights <- function(data=NULL, gpar=NULL, v=1) {
    G = length(gpar$pi)
    if (G > 1) {
        fx = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
        for (k in 1:G ) fx[,k] = ddmsghyp(data=data, par=gpar[[k]] )
        
        w = t(apply(fx, 1, function(z,wt,v) {
            x= (z*wt)^v;
            
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
    } else w = matrix(1,nrow=nrow(data), ncol=G)
    return(w)
}


##############uig  mixtures weights


iweights <- function(data=NULL, gpar=NULL, v=1) {
    fx = matrix(0, nrow=nrow(data), ncol=2)
    
    fx[,1] =   ddghyp(y=data, par=gpar, log=FALSE)
    fx[,2] =  dmsghyp(y=data, par=gpar, log=FALSE)
    
    if (gpar$wg[2] != 0 ) {
        iw = t(apply(fx, 1, function(z=NULL,wt=NULL,v=NULL) {
            x= (z*wt)^v;
            
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$wg,v=v ))
    } else iw =cbind( rep(1, nrow(data)), rep(0, nrow(data)) )
    
    return(iw)
}




############Expected values GIG


###Bessel function
logbesselKv <- function(x, y) { log(besselK(x=y, nu=x, expon.scaled=TRUE)) - log(y)}


logbesselKvFA <- function(x, y) {
    val = log(besselK(x=y, nu=x, expon.scaled=FALSE))
    sun = is.infinite(val)
    val[sun] = besselK.nuAsym(x=y[sun], nu=abs(x[sun]), k.max=4, log=TRUE, expon.scaled=FALSE)
    return(val)
}


besselKv    <- function(x, y) { besselK(x=y, nu=x)}


######GIG univariate
gig <- function(a=NULL,b=NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3
    sab =  sqrt(a*b)
    kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
    kv  = besselK( sab, nu=v, expon.scaled =TRUE)
    kv12 = kv1/kv
    
    sb.a = sqrt(b/a)
    w    = kv12*sb.a
    invw = kv12*1/sb.a - 2*v/b
    logw = log(sb.a) + grad( logbesselKvFA, x=rep(v,length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    
    
    val =cbind(w,invw,logw)
    return(val)
}

gig2 <- function(x=NULL, par=NULL, invS=NULL) {
    # returns the same as gig
    d = length(par$mu)
    
    omega = par$cpl0[1]
    a1 = omega + sum(par$alpha^2/par$phi )
    xmu= 	sweep(x,2,par$mu, "-")
    b1 = omega + apply( 	xmu^2, 1, weighted.sum, wt=1/par$phi)
    v1 = par$cpl0[2]-d/2
    
    val = gig(b=b1,a=a1, v=v1)
    return(val)
}

gig20 <- function(sr2=NULL, par=NULL) {
    # returns the same as gig
    d = length(par$mu)
    
    omega = par$cpl0[1]
    a1 = omega + sum(par$alpha^2/par$phi )
    b1 = omega + sr2
    v1 = par$cpl0[2]-d/2
    
    val = gig(b=b1,a=a1, v=v1)
    return(val)
}


#####GIG MS
gigp <- function(a=NULL,B =NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3
    SaB =  sqrt(sweep(B, 2, a, "*" ))
    
    Kv12 = t(apply(SaB, 1, function(x=NULL, v=NULL) {
        kv1 = besselK( x, nu=v+1, expon.scaled =TRUE)
        kv0 = besselK( x, nu=v, expon.scaled =TRUE)
        kv12 = kv1/kv0
        return(kv12)
    }, v=v ) )
    
    SB.a =  sqrt(sweep(B, 2, 1/a, "*" ))
    if(nrow(Kv12)==1) {Kv12=t(Kv12) }
    W    = Kv12*SB.a
    invW = Kv12/SB.a - 2*sweep(1/B, 2, v, "*" )
    #print(cbind(rep(v,each=nrow(SaB)),as.numeric(SaB)) )
    
    logW = log(SB.a) + matrix(grad( logbesselKvFA, x=rep(v,each=nrow(SaB)), y=as.numeric(SaB), method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)), nrow=nrow(SaB), ncol=ncol(SaB) )
    
    val = list(W=W,invW=invW,logW=logW)
    return(val)
}


gig2p <- function(sr2=NULL, par=NULL) {
    # returns the same as gig
    
    omega = par$cpl[,1]
    a1 = omega + par$alpha^2/par$phi   # vector
    B1 = sweep(sr2, 2, omega, "+") # matrix
    v1 = par$cpl[,2]-1/2 # vector
    
    val = gigp(a=a1, B=B1, v=v1)
    return(val)
}



weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )

########################### M step Parameter estimation ##################################

################Initialization
rmgpar <- function(g=NULL, p=NULL,data=NULL, method="kmeans",n=10, nr=10) {
    if(g==1){
        # lk=kmeans(data,g)
        lc=as.matrix(t(apply(data,2,mean)))
        l=as.vector(rep(1,nrow(data)))#lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  try({ gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)}, TRUE)
        return(gpar)
    }
    else{    if(method=="modelBased"){
        lk=gpcm(data,  G=g, mnames=c("VVV"))
        l=lk$map
        lc=lk$gpar[[1]]$mu
        for(i in 2:g){
            lc=rbind(lc,lk$gpar[[i]]$mu)
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)
        # l=BAR(data,l)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        
        lc=apply(data[l==1,],2,mean)
        for(i in 2:g){
            lc=rbind(lc,apply(data[l==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n)  gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)
        
        #l=BAR(data,l)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        llkO=-Inf
        for(i in 1:nr){
            l=round(runif(nrow(data))*(g-1)+1)
            
            lc=apply(data[l==1,],2,mean)
            for(i in 2:g){
                lc=rbind(lc,apply(data[l==i,],2,mean))
            }
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gparO  = rgparC(data=data, g=g, w=z,l=lc)
            
            loglik = numeric(100)
            for (i in 1:3) {
                gparO = EMgrstep(data=data, gpar=gparO, v=1, label = l,it=i)
                loglik[i] = llik(data, gparO)
            }
            
            while ( ( getall(loglik[1:i]) > 1) & (i < (100) ) )  {
                i = i+1
                gparO = EMgrstep(data=data, gpar=gparO, v=1, label = l,it=i)
                loglik[i] = llik(data, gparO)
                
            }
            llk=llik(data,gparO)
            if(llk>llkO){
                llkO=llk
                gpar=gparO
                
            }
        }
        return(gpar)
    }
    
    else if(method=="kmedoids"){
        lk=pam(data,g)
        lc=lk$medoids
        l=lk$clustering#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n) gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)}
    
    
    else{ lk=kmeans(data,g)
        lc=lk$centers
        l=lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n) gpar = EMgrstep(data=data, gpar=gpar, v=1, label= NULL,  w=z, it=j)
        return(gpar)}}
    
    
}



rgparC<- function(data,g=NULL, w=NULL,l=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    # l=rmgparMS(g, p,data,method=method)
    for (k in 1:g) val[[k]] = rmparC(data=data,wt=w[,k],k,lc=l)
    
    val$pi = rep(1/g,g)
    return(val)
}

rmparC <- function(data,wt,k,lc=NULL) {
    par = list()
    p=ncol(data)
    par$mu =lc[k,]#rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )#l[k,] #rnorm(p,0,.01);
    #par$phi  =  rnorm(p, 0, sqrt(1/nrow(data)) )#rep(1,p);
    par$alpha = rep(0,p)# rnorm(p, 0, sqrt(1/nrow(data)) )#rnorm(p,0,.01);
    
    sigma = ( ( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(sigma)=abs(diag(sigma))
    if(p==1){sigma=var(data)}
    
    for(i in 1:p){if(sigma[i,i]<0.1){sigma[i,i]=0.1}}
    if (any(eigen(sigma)$values <= 0 ) ) par$sigma =  diag(apply(data,2,var))
    par$gam   = eigen( sigma)$vectors
    par$phi  = eigen( sigma)$values
    par$sigma=sigma
    par$cpl0 = c(1,-1/2)
    par$cpl = cbind( rep(1,p), rep(-1/2,p))
    wg <- 0.5#runif(1)
    par$wg = c(wg,1-wg)
    return(par)
}
########################## MAIN



updatemaScplp <- function(y=NULL, par=NULL, weights=NULL, iweights=NULL, alpha.known=NULL, v=1,it=NULL) {
    p=ncol(y)
    n=nrow(y)
    if (is.null(weights)) weights=rep(1,n)
    if (is.null(iweights)) iweights=matrix(c(0.5),nrow=n,ncol=2)
    
    x = y %*% (par$gam)
    sr2 = sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    abc  = gig20(sr2=apply(sr2,1,sum), par=par)
    
    sumwi = apply(iweights,2,mean)
    sumw = sum(weights)
    
    if (sumwi[2] != 0 ) {
        abcM = gig2p(sr2=sr2, par=par)
        S1 = abc[,1]*iweights[,1] +    abcM$W*iweights[,2]
        S2 = abc[,2]*iweights[,1] + abcM$invW*iweights[,2]
    } else {
        S1 = matrix(abc[,1], nrow=n, ncol=p)
        S2 = matrix(abc[,2], nrow=n, ncol=p)
    }
    
    meanS1 = apply(S1, 2, weighted.mean, w= weights)
    meanS2 = apply(S2, 2, weighted.mean, w= weights)
    S      = (sweep(S2,2, meanS1, "*") -1)*weights
    sums   = apply(S, 2, sum)
    
    mu.new    = as.vector(apply(x*S,2,sum)/sums)
    alpha.new = as.vector(apply(x*sweep(-S2, 2, meanS2,"+")*weights, 2, sum )/sums)
    
    #PHI
    x.mu     = sweep(x,2,mu.new, "-")
    xmu.mean = apply(x.mu, 2, weighted.mean, w= weights)
    xmu2     = apply(x.mu^2*S2, 2, weighted.mean, w= weights)
    phi.new  = xmu2 -2*xmu.mean*alpha.new + meanS1*alpha.new^2
    
    ################CPLMs
    if( sumwi[2] == 0.01) {
        A = apply(abcM$W,    2, weighted.mean, w=weights*iweights[,2]  )
        B = apply(abcM$invW, 2, weighted.mean, w=weights*iweights[,2] )
        C = apply(abcM$logW, 2, weighted.mean, w=weights*iweights[,2] )
        cpl.new = t(apply(cbind(par$cpl, A, B, C), 1, function(z) {
            temp = updateol(ol= z[1:2], ABC=z[3:5], n=2)
            return( temp )
        }))
    } else cpl.new=par$cpl
    
    
    ##############CPL uni
    ABC = apply(abc, 2, weighted.mean, w=weights*iweights[,1])
    cpl.newU =  updateol(ol=par$cpl0, ABC=ABC, n=2)
    
    if (it %% 2 ==0) new.gam = updategam1MS(gam0=par$gam, y=y, sigma= phi.new, alpha= alpha.new, mu= mu.new, wt=weights, invW= S2)
    else new.gam = updategam2MS(gam0=par$gam, y=y, sigma= phi.new, alpha= alpha.new, mu= mu.new, wt=weights, invW= S2)
    
    new.par = list(mu=mu.new, phi=phi.new, alpha=alpha.new, cpl=cpl.new, gam= new.gam,  cpl0= cpl.newU )
    return(new.par)
}




######################## MS cpl
updateol <- function(ol=NULL, ABC=NULL, n=1) {
    #print( c(88888888, ol,ABC) )
    ol0 = ol
    for (i in 1:n) {
        if (ABC[3] == 0) {
            ol[2] = 0
        } else {
            #print( c(1, ol[2], ol[1] ) )
            #if(ol[1]<=0){ol[1]=eps}
            #      print(ol)
            bv = grad( logbesselKvFA, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
            ol[2] = ABC[3]*(ol[2]/bv)
        }
        
        Rp = Rlam(ol[1],lam=+ol[2])
        Rn = Rlam(ol[1],lam=-ol[2])
        f1 = Rp + Rn - (ABC[1]+ABC[2])
        f2 = ( Rp^2 - (2*ol[2]+1)/ol[1] *Rp -1 ) + ( Rn^2 - (2*(-1*ol[2])+1)/ol[1] *Rn -1 )
        
        if ( ol[1] - f1/f2 > 0) ol[1] = ol[1] - f1/f2
        
    }
    return(ol)
}



Rlam <- function(x, lam=NULL) {
    v1 = besselK(x, nu= lam+1, expon.scaled=FALSE)
    v0 = besselK(x, nu= lam, expon.scaled=FALSE)
    val = v1/v0
    
    set= is.infinite(v1)|is.infinite(v0)|v0==0|v1==0
    if(any(set)){
        lv1=besselK.nuAsym(x=x[set],nu=abs(lam[set]+1),k.max=4, expon.scaled=FALSE,log=TRUE)
        lv0=besselK.nuAsym(x=x[set],nu=abs(lam[set]),k.max=4, expon.scaled=FALSE,log=TRUE)
        val[set]=exp(lv1-lv0)
    }
    return(val)
}





################################## gamma


updategam1 <- function(gam0=NULL, y=NULL, phi=NULL, alpha=NULL, mu=NULL, wt=NULL, ut=NULL, invW=NULL, abc=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/phi, "*" )
    #    Bs2 = matrix(abc[,2], nrow(x),ncol(x) )/phi
    Bs2 = sweep(matrix(abc[,2], nrow(x),ncol(x) ), 2, 1/phi, "*" )
    
    # u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt*ut[,2], "*")
    wty2 = sweep(y,1, wt*ut[,1],"*")
    
    F0t = t( x * Bs ) %*% wty + t(x*Bs2) %*% wty2
    e2 = apply(Bs,1,max)*wt*ut[,2]
    e1 = max(Bs2)*wt*ut[,1]
    
    if( sum(e2)==0 ) {    F1t =  ((cov.wt(y, center=rep(0,ncol(y)), wt=e1, method="ML" )$cov*sum(e1)) %*% (gam0))
        
    } else {
        if(sum(e1)==0){F1t = ((cov.wt(y, center=rep(0,ncol(y)), wt=e2, method="ML" )$cov*sum(e2)) %*% (gam0))
            
        } else {
            F1t = ((cov.wt(y, center=rep(0,ncol(y)), wt=e2, method="ML" )$cov*sum(e2)) %*% (gam0)) + ((cov.wt(y, center=rep(0,ncol(y)), wt=e1, method="ML" )$cov*sum(e1)) %*% (gam0))
            
        }
    }
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
    v1 = (Bs2*mu + alpha/phi)
    A0 = t(v)  %*% wty
    A1 = t(v1) %*% wty2
    
    print("ryan1")
    F2 = F0t - t(F1t) - A0 - A1
    z  = svd(-F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    #gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
    
    #   if (obj.gam(gam0=gam,y=y, phi=phi, alpha=alpha, mu=mu,  wt=wt, ut=ut, invW=invW,abc=cpl0) <  obj.gam(gam0=gam0,y=y, phi=phi, alpha=alpha, mu=mu, wt=wt, ut=ut, invW=invW,abc=cpl0) ) {
    #		print('not good enough')
    #       gam = gam0
    #   }
    
    return(gam)
}


updategam2 <- function(gam0=NULL, y=NULL,  phi=NULL, alpha=NULL, mu=NULL, wt=NULL, ut=NULL, invW=NULL, abc=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/phi, "*" )
    Bs2 =  matrix(abc[,2],nrow(x),ncol(x))/phi
    # u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt*ut[,2], "*")
    wty2 = sweep(y,1,wt*ut[,1],"*")
    F0t = t( x * Bs ) %*% wty + t(x*Bs2) %*% wty2
    
    e1 = apply(y^2, 1, sum)*wt*ut[,1]
    e2 = apply(y^2, 1, sum)*wt*ut[,2]
    
    if(ncol(y)==1){F1t = (apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0) + weighted.sum(Bs2,wt=e1)*t(gam0)}
    else{
        F1t = diag(apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0) + weighted.sum(Bs2,wt=e1)*t(gam0)}
    
    
    #v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+") + (Bs2*mu + alpha/phi)
    #A0 = t(v)  %*% wty
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
    v1 = (Bs2*mu + alpha/phi)
    A0 = t(v)  %*% wty
    A1 = t(v1) %*% wty2
    
    
    F2 = F0t - F1t - A0 -A1
    z  = svd(F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    # gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
    
    #   if (obj.gam(gam0=gam,y=y, phi=phi, alpha=alpha, mu=mu,  wt=wt, ut=ut, invW=invW,abc=cpl0) <  obj.gam(gam0=gam0,y=y, phi=phi, alpha=alpha, mu=mu, wt=wt, ut=ut, invW=invW,abc=cpl0) ) {
    #		print('not good enough')
    #       gam = gam0
    #   }
    
    return(gam)
}


updategam2MS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/sigma, "*" )
    #	u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt, "*")
    F0t = t( x * Bs ) %*% wty
    
    e2  = apply(y^2, 1, sum)*wt
    
    if( ncol(y)==1 ) F1t = (apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0)
    else F1t = diag(apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0)
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
    A0 = t(v)  %*% wty
    
    F2 = F0t - F1t - A0
    z  = svd( -F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    
    if (objgamMS(gam0=gam,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) <  objgamMS(gam0=gam0,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) ) {
        #		print('not good enough')
        gam = gam0
    }
    
    return(gam)
}



tr <- function(A=NULL) { sum(diag(A)) }



obj.gam <- function(gam0=NULL, y=NULL, phi=NULL, alpha=NULL, mu=NULL, wt=NULL, ut=NULL, invW=NULL,abc=NULL) {
    
    invUmat=matrix(abc[,2],nrow(y),ncol(y))
    ###MS part
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs  = sweep(invW, 2, 1/phi, "*" )
    wtx = sweep(x, 1, (wt*ut[,2]), "*")
    F0t = t( x * Bs ) %*% wtx
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
    A0 = t(v)  %*% wtx
    
    valMS = -1*(tr(F0t ) - 2*tr( A0  ))/sumw
    
    
    ###Multivariate part
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs  = sweep(invUmat, 2, 1/phi, "*" )
    wtx = sweep(x, 1, (wt*ut[,1]), "*")
    F0t = t( x * Bs ) %*% wtx
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
    A0 = t(v)  %*% wtx
    
    valMUL = -1*(tr(F0t ) - 2*tr( A0  ))/sumw
    
    val=valMS+valMUL
    
    return(val)
}



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                       Factor analyzer                                                     ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

mainMGHFA<-function(data=NULL, gpar0, G, n, label  , eps, method ,q,nr=nr ) {
    pcol=ncol(data)
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgpar(data=data, g=G, w=z,l=lc)
        
    }
    else{
        if (is.null(gpar0)) gpar = igparM(data=data, g=G,q=q,method=method,nr=nr)
        else gpar  = gpar0}
    loglik = numeric(n)
    for (i in 1:3) {
        gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)	###parameter estimation
        loglik[i] = llikFA(data, gpar) ##likelyhood
    }
    while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
        i = i+1
        gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)	###parameter estimation
        loglik[i] = llikFA(data, gpar) ##likelyhood
    }
    if(i<n){loglik=loglik[-(i+1:n)]}
    BIC=2*loglik[i]-log(nrow(data))*((G-1)+G*(3*pcol+2+pcol*q-q*(q-1)/2))
    val = list(loglik= loglik, gpar=gpar, z=weightsFA(data=data, gpar= gpar), map=MAPFA(data=data, gpar= gpar, label=label) , BIC=BIC)
    return(val)
}
####################### Univ CPL



RRlamz <- function(x,lam=NULL,z=1) {
    val =RlamU(x, lam=-lam)+ RlamU(x, lam= lam)
    zval = val - z
    return(zval)
}

RlamU <- function(x, lam=NULL) {
    
    v1 = besselK(x, nu= lam+1, expon.scaled=TRUE)
    v0 = besselK(x, nu= lam, expon.scaled=TRUE)
    val = v1/v0
    return(val)
}

updateolU <- function(ol=NULL, ABC=NULL, n=1) {
    
    for (i in 1:n) {
        if (ABC[3] == 0) {
            ol[2] = 0
        } else {
            # if(ol[1]<=0){ol[1]=eps}
            bv = grad( logbesselKv, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
            ol[2] = ABC[3]*(ol[2]/bv)
        }
        
        
        lam0 = ol[2]
        omg0 = ol[1]
        Rp = RlamU(omg0,lam=+lam0)
        Rn = RlamU(omg0,lam=-lam0)
        f1 = Rp + Rn - (ABC[1]+ABC[2])
        f2 = ( Rp^2 - (2*lam0+1)/omg0 *Rp -1 ) + ( Rn^2 - (2*(-1*lam0)+1)/omg0 *Rn -1 )
        # note, it is suppose to be f1/2 and f2/2
        if ( ol[1] - f1/f2 > 0 ) ol[1] = ol[1] - f1/f2
    }
    
    return(ol)
}





weightsFA <- function(data=NULL, gpar=NULL, v=1) {
    G = length(gpar$pi)
    if (G > 1) {
        zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
        for (k in 1:G ) zlog[,k] =  ddghypFA(x=data, par=gpar[[k]], log=TRUE)
        w = t(apply(zlog, 1, function(z,wt,v) {
            x=exp(v*(z + log(wt)) );
            sun=is.infinite(x)
            x[sun]=1
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
    } else w = matrix(1,nrow=nrow(data), ncol=G)
    return(w)
}



###Function1




igparM <- function(data=NULL, g=NULL,q=2,method="kmeans",nr=NULL) {
    ##initialization
    gpar = igpar3M(data=data, g=g, n=10,q=q,method=method,nr=nr)
    return(gpar)
}


####Function for igparM
igpar3M <- function(data=NULL, g=NULL, n=10,label=NULL,q=2,method="kmeans",nr=10) {
    if(g==1){
        # lk=kmeans(data,g)
        lc=as.matrix(t(apply(data,2,mean)))
        l=as.vector(rep(1,nrow(data)))#lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgpar(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  try({ gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label= l,  w=z)}, TRUE)
        return(gpar)
    }
    else{if(method=="modelBased"){
        lk=gpcm(data,  G=g, mnames=c("VVV"))
        l=lk$map
        lc=lk$gpar[[1]]$mu
        for(i in 2:g){
            lc=rbind(lc,lk$gpar[[i]]$mu)
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgpar(data=data, g=g, w=z,l=lc,q=q)
        
        for (j in 1:n)  gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label= l,  w=z)
        return(gpar)
        # l=BAR(data,l)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        lc=apply(data[l==1,],2,mean)
        for(i in 2:g){
            lc=rbind(lc,apply(data[l==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgpar(data=data, g=g, w=z,l=lc,q=q)
        
        for (j in 1:n)  gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label= l,  w=z)
        return(gpar)
        
        #l=BAR(data,l)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        llkO=-Inf
        for(i in 1:nr){
            l=round(runif(nrow(data))*(g-1)+1)
            
            lc=apply(data[l==1,],2,mean)
            for(i in 2:g){
                lc=rbind(lc,apply(data[l==i,],2,mean))
            }
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gparO  = rgpar(data=data, g=g, w=z,l=lc,q=q)
            
            
            
            loglik = numeric(100)
            for (i in 1:3) {
                gparO = EMgrstepFA(data=data, gpar=gparO, v=1, label = l,w=z)
                loglik[i] = llikFA(data, gparO)
            }
            
            while ( ( getall(loglik[1:i]) > 1) & (i < (100) ) )  {
                i = i+1
                gparO = EMgrstepFA(data=data, gpar=gparO, v=1, label = l,w=z)
                loglik[i] = llikFA(data, gparO)
                
            }
            llk=llikFA(data,gparO)
            if(llk>llkO){
                llkO=llk
                gpar=gparO
                
            }
        }
        return(gpar)
    }
    
    else if(method=="kmedoids"){
        lk=pam(data,g)
        lc=lk$medoids
        l=lk$clustering#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgpar(data=data, g=g, w=z,l=lc,q=q)
        
        for (j in 1:n) gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label= l,  w=z)
        return(gpar)}
    
    
    
    else{ lk=kmeans(data,g)
        lc=lk$centers
        l=lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgpar(data=data, g=g, w=z,l=lc,q=q)
        
        for (j in 1:n)  gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label= l,  w=z)
        return(gpar)}}
    
}
rgpar <- function(data, g=2, w=NULL,l=NULL, q=2) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    for (k in 1:g) val[[k]] = iparM(data=data, wt=w[,k],q=q,lc=l[k,])
    val$pi = rep(1/g,g)
    return(val)
}


iparM <- function(data, wt=NULL,q=2,lc=NULL) {
    if (is.null(wt)) wt = rep(1,nrow(data))
    p = ncol(data)
    val = list()
    val$mu    = lc#rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )
    val$alpha = rep(0,p)#rnorm(p, 0, sqrt(1/nrow(data)) )
    e=eigen(cov.wt(data, wt = wt, method="ML")$cov)
    dia=diag(e$values[1:q])
    if(q==1){val$Lambda=as.matrix(e$vectors[,1],p,q)
        dia=e$vectors[,1:q]
        val$sigma =val$Lambda%*%t(val$Lambda)
    }
    else{
        val$Lambda=e$vectors[,1:q]
        val$sigma =val$Lambda%*%dia%*%t(val$Lambda)######Cris###################################################################################
    }
    val$err=diag(cov.wt(data, wt = wt, method="ML")$cov-val$sigma)*diag(p)########Cris########################################################################################################################
    #  val$sigma = diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    if (any(eigen(val$sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    val$cpl   = c(1,-1/2)
    return(val)
    
}
#####end function for igpar

EMgrstepFA <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL) {
    ##parameter estimtion
    if (is.null(w)) w = weightsFA(data=data, gpar=gpar,v=v)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    
    G= length(gpar$pi);##number clusters
    d= length(gpar[[1]]$mu);##number variables
    #Sk = array(0, c(d,d,G) )
    for (k in 1:G ) {
        
        gpar[[k]] = updatemaScplM1(x=data, par=gpar[[k]], weights=w[,k], invS=NULL, alpha.known=NULL, v=v)
        w = weightsFA(data=data, gpar=gpar,v=v)
        gpar[[k]] = updatemaScplM2(x=data, par=gpar[[k]], weights=w[,k], invS=NULL, alpha.known=NULL, v=v)
        #Sk[,,k]   = gpar[[k]]$sigma
        
    }
    gpar$pi = apply(w,2,mean)
    return(gpar)
}








##################################################NEW###############################################################################

##################################################NEW###############################################################################






updatemaScplM1 <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2FA(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        A = ABC[1]
        B = ABC[2]
        u = (B - abc[,2])*weights
        t = (A*abc[,2]-1)*weights
        T = sum(t)
        
        mu.new    = apply(x, 2, weighted.sum, w=t)/T
        alpha.new = apply(x, 2, weighted.sum, w=u)/T
    } else {
        alpha.new = alpha.known
        mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) - alpha.new/ABC[2]
    }
    alpha.new = alpha.new*v
    
    
    
    
    par.old = c(  log(par$cpl[1]) , par$cpl[2])
    #a = optim( par.old, loglikgig4, ABC=ABC)
    par.ol = c(exp(par.old[1]), par.old[2])
    a = updateolU(ol=par.ol, ABC=ABC, n=2)
    #if ( loglikgig4(par.old, ABC) < loglikgig4(par.old, ABC) ) cpl.new = par$cpl
    #  else cpl.new =  c( rep(exp(a$par[1]), 2), a$par[2])
    #if(a[1]<=0){a[1]=eps}
    cpl.new =  c( a[1], a[2])
    
    
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=par$sigma, cpl=cpl.new, Lambda=par$Lambda, err=par$err )
    
    return(new.par)
}












updatemaScplM2 <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2FA(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    mu.new   = par$mu
    
    alpha.new = par$alpha
    
    
    cpl.new =  par$cpl
    
    
    
    
    #####second estep
    d=ncol(x)
    q=ncol(par$Lambda)
    var=par$Lambda
    fi=par$err
    fi1=ginv(fi)
    fi1=diag(d)*as.vector(fi1)
    dia=diag(q)
    if( length(q)==0){dia=q}
    inv=fi1-fi1%*%var%*%ginv(dia+t(var)%*%fi1%*%var)%*%t(var)%*%fi1
    beta=t(var)%*%inv
    term1=sweep(x,2, FUN ="-", STATS=mu.new)-matrix(alpha.new,nrow(x),d,1)*matrix(abc[,1],nrow(x),d)
    uhat=beta%*%t(term1)
    term1=sweep(x,2, FUN ="-", STATS=mu.new)*matrix(abc[,2],nrow(x),d)-matrix(alpha.new,nrow(x),d,1)
    
    uhatb=beta%*%t(sweep(x,2, FUN ="-", STATS=mu.new)*matrix(abc[,2],nrow(x),d)-matrix(alpha.new,nrow(x),d,1))
    
    
    A = cov.wt(x, wt=abs(abc[,2]*weights), center=mu.new, method="ML")$cov*ABC[2] #returns a list containing weighted covariance matrix and the mean of the data
    r = apply(x, 2, weighted.sum, wt=weights)/sumw
    R = A - (outer( r - mu.new, alpha.new) + outer(alpha.new, r - mu.new)) + outer(alpha.new,alpha.new)*ABC[1]
    
    
    euu=sum(abc[,2]*weights)*diag(q)-sum(abc[,2]*weights)*beta%*%var+beta%*%R%*%t(beta)*sumw
    
    #euu=diag(q)-beta%*%var+beta%*%t(term1)%*%term1%*%t(beta)
    va.new=(t(sweep(x,2, FUN ="-", STATS=mu.new)*(weights))%*%t(uhatb)-t(matrix(alpha.new,nrow(x),d,1)*(weights))%*%t(uhat))%*%ginv(euu)#/sum(abc[,2]*weights)
    
    # va.new=(t(sweep(x,2, FUN ="-", STATS=mu.new)*weights)%*%t(uhatb)-(alpha.new)%*%t(apply(weights*uhat,1,"sum")))%*%ginv(euu)
    
    
    #euu=diag(q)-beta%*%var+beta%*%t(term1*abc[,2]*weights)%*%term1%*%t(beta)
    fip1=t(sweep(x,2, FUN ="-", STATS=mu.new)*weights)%*%t(uhatb)%*%t(var)
    fip2=(alpha.new)%*%t(apply(weights*uhat,1,"sum"))%*%t(var)+(var%*%euu%*%t(var))
    
    fi.new=diag(R+(2/sumw)*(-fip1+fip2))
    # fi.new=diag(R+(t(-2*t(sweep(x,2, FUN ="-", STATS=mu.new)*weights))%*%t(uhatb)%*%t(var)+2*(alpha.new)%*%t(apply(weights*uhat,1,"sum"))%*%t(var)+(var%*%euu%*%t(var)))/sumw)
    # euu=diag(q)-beta%*%var+beta%*%t(term1*abc[,2]*weights)%*%term1%*%t(beta)
    #fi.new=diag(R-(t(2*abc[,2]*term1*weights)%*%t(uhat)%*%t(var)-sum(abc[,2]*weights)*(var%*%euu%*%t(var)))/sumw)
    #fi.new=diag((t(abc[,2]*term1*weights)%*%(term1)-t(2*abc[,2]*term1*weights)%*%t(uhat)%*%t(var)+sum(abc[,2]*weights)*(var%*%euu%*%t(var)))/sumw)
    # fi.new=1/fi.new
    #fi.new=par$err
    
    #va.new=R%*%t(inv)%*%var
    # inv=fi1-fi1%*%va.new%*%ginv(dia+t(va.new)%*%fi1%*%va.new)%*%t(va.new)%*%fi1
    #fi.new=(R%*%t(inv))*diag(d)%*%fi
    
    var=va.new%*%t(va.new)+diag(d)*fi.new
    
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=var, cpl=cpl.new, Lambda=va.new, err=fi.new )
    
    return(new.par)
}









MAPFA <- function(data, gpar, label=NULL) {
    w = weightsFA(data=data, gpar=gpar, v=1)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
    z = as.numeric(z)
    return( z)
}




gig2FA <- function(x=NULL, par=NULL, invS=NULL) {
    #  returns a matrix with dim length(a) x 3
    d = length(par$mu)
    
    #if (is.null(invS)) invS = ginv(par$sigma)
    alpha = par$alpha
    omega = exp(  log(par$cpl[1])  )
    
    a1 = omega + as.numeric( alpha %*% invS %*% alpha )
    b1 = omega + as.numeric(mahalanobis(x, center=par$mu, cov=invS, inverted=TRUE))
    v1 = par$cpl[2]-d/2
    
    val = gigFA(b=b1,a=a1, v=v1)
    return(val)
}

gigFA <- function(a=NULL,b=NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3 stima le yi
    sab =  sqrt(abs(a*b))
    kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
    kv  = besselK( sab, nu=v, expon.scaled =TRUE)
    kv12 = kv1/kv
    #sun=is.nan(kv12)
    #kv12[sun]=1
    
    sb.a = sqrt(abs(b/a))
    w    = kv12*sb.a
    invw = kv12*1/sb.a - 2*v/b
    
    #sqr1w=sqrt(sb.a)*besselK( sab, nu=v-0.5, expon.scaled =TRUE)/kv
    logw = log(sb.a) + grad( logbesselKv, x=rep(abs(v),length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    
    val =cbind(w,invw,logw)#,sqr1w,w2)
    return(val)
}






loglikgig4 <- function(ol=NULL, ABC=ABC) {
    val = numeric(3)
    omega  =  exp(ol[1])
    val[1] = -log(besselK( x=omega, nu=ol[2], expon.scaled=TRUE)) + omega
    val[2] = (ol[2]-1)*ABC[3]
    val[3] = -1/2* omega*(ABC[2] + ABC[1])
    val = -1*sum(val)
    return(val)
}



###Function2

##for funtion 2





ddghypFA <- function(x, par, log=FALSE, invS=NULL) {
    # x  is a n * p matrix
    # generalized hyperbolic distribution
    
    mu    = par$mu
    sigma = par$sigma
    alpha = par$alpha
    
    #check for input error
    if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
    d = length(mu)
    omega = exp(log(par$cpl[1]))
    lambda = par$cpl[2];
    
    
    #distances between points and mu
    if (is.null(invS)) invS = ginv(sigma)
    pa = omega + alpha %*% invS %*% alpha
    mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
    
    kx = sqrt(abs(mx*pa))####################################################################################################################
    xmu = sweep(x,2,mu)
    
    lvx = matrix(0, nrow=nrow(x), 4)
    lvx[,1] = (lambda - d/2)*log(kx)
    lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
    lvx[,3] = xmu %*% invS %*% alpha
    
    
    lv = numeric(6)
    
    lv[1] = -1/2*log(abs(det( sigma ))) -d/2*(log(2)+log(pi))##################################################################################
    lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
    lv[3] = -lambda/2*( log(1) )
    lv[4] = lambda*log(omega) *0
    lv[5] = (d/2 - lambda)*log(abs( pa) )############################################################################################
    
    val = apply(lvx,1,sum) + sum(lv)
    if (!log) val = exp( val )
    
    return(val)
}



#### combinewk see Function for igpar




##end for function 2


###Function3


llikFA <- function(data,gpar, delta=0) {
    ##log likelyhood estimation
    logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:length(gpar$pi) ) logz[,k] = ddghypFA(x=data, par=gpar[[k]], log=TRUE)
    val = sum(log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
    return(val)
}


weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                              GHD                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

mainMGHD<-function(data=NULL, gpar0, G, n, label  , eps, method  ,nr=NULL) {
    
    pcol=ncol(data)
    
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgparGH(data=data, g=G, w=z,l=lc)
        
    }
    else{
        if (is.null(gpar0)) gpar = try(igpar(data=data, g=G, method=method,nr=nr))
        else gpar  = gpar0}
    
    
    loglik = numeric(n)
    for (i in 1:3) {
        gpar = try(EMgrstepGH(data=data, gpar=gpar, v=1, label = label))	###parameter estimation
        loglik[i] = llikGH(data, gpar)}
    
    while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
        i = i+1
        gpar = try(EMgrstepGH(data=data, gpar=gpar, v=1, label = label))	###parameter estimation
        
        loglik[i] = llikGH(data, gpar) ##likelyhood
    }
    if(i<n){loglik=loglik[-(i+1:n)]}
    BIC=2*loglik[i]-log(nrow(data))*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
    z=weightsGH(data=data, gpar= gpar)
    ICL=BIC+sum(log(apply(z,1,max)))
    AIC=2*loglik[i]-2*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
    AIC3=2*loglik[i]-3*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
    val = list(loglik= loglik, gpar=gpar, z=z, map=MAPGH(data=data, gpar= gpar, label=label),BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3 )
    return(val)
    
}

###Function1

igpar <- function(data=NULL, g=NULL, method="kmeans",nr=NULL) {
    ##initialization
    gpar = igpar3(data=data,g=g, n=10, method=method,nr=nr)
    return(gpar)
}


####Function for igpar
igpar3 <- function(data=NULL, g=NULL, n=10,label=NULL, method="kmeans",nr=10) {
    if(g==1){
        # lk=kmeans(data,g)
        lc=as.matrix(t(apply(data,2,mean)))
        l=as.vector(rep(1,nrow(data)))#lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparGH(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= l,  w=z)}, TRUE)
        return(gpar)
    }
    else{
        if(method=="modelBased"){
            lk=gpcm(data,  G=g, mnames=c("VVV"))
            l=lk$map
            lc=lk$gpar[[1]]$mu
            for(i in 2:g){
                lc=rbind(lc,lk$gpar[[i]]$mu)
            }
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gpar  = rgparGH(data=data, g=g, w=z,lc)
            
            for (j in 1:n)  try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= l,  w=z)}, TRUE)
            return(gpar)
            # l=BAR(data,l)
        }
        else if(method=="hierarchical"){
            l=(cutree(hclust(dist(data),"ward.D"), k=g))
            #lk=gpcm(data,  G=g, mnames=c("VVV"),label=l)
            lc=apply(data[l==1,],2,mean)
            for(i in 2:g){
                lc=rbind(lc,apply(data[l==i,],2,mean))
            }
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gpar  = rgparGH(data=data, g=g, w=z,lc)
            
            for (j in 1:n)  try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= l,  w=z)}, TRUE)
            return(gpar)
            
            #l=BAR(data,l)
        }
        else if(method=="random"){
            #l=kmeans(data,g)
            llkO=-Inf
            for(i in 1:nr){
                l=round(runif(nrow(data))*(g-1)+1)
                lc=apply(data[l==1,],2,mean)
                for(i in 2:g){
                    lc=rbind(lc,apply(data[l==i,],2,mean))
                }
                z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
                gparO  = rgparGH(data=data, g=g, w=z,l=lc)
                
                
                loglik = numeric(100)
                for (i in 1:3) {
                    gparO = try({EMgrstepGH(data=data, gpar=gparO, v=1, label = l,  w=z)}, TRUE)
                    loglik[i] = llikGH(data, gparO)
                }
                
                while ( ( getall(loglik[1:i]) > 1) & (i < (100) ) )  {
                    i = i+1
                    gparO = try({EMgrstepGH(data=data, gpar=gparO, v=1, label = l,  w=z)}, TRUE)
                    loglik[i] = llikGH(data, gparO)
                    
                }
                
                
                llk=llikGH(data,gparO)
                if(llk>llkO){
                    llkO=llk
                    gpar=gparO
                    
                }
            }
            return(gpar)
        }
        
        else if(method=="kmedoids"){
            lk=pam(data,g)
            lc=lk$medoids
            l=lk$clustering#$centers}
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gpar  = rgparGH(data=data, g=g, w=z,l=lc)
            
            for (j in 1:n) try({gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= l,  w=z)}, TRUE)
            return(gpar)}
        
        else{ lk=kmeans(data,g)
            lc=lk$centers
            l=lk$cluster#$centers}
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gpar  = rgparGH(data=data, g=g, w=z,lc)
            
            for (j in 1:n)  try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= l,  w=z)}, TRUE)
            return(gpar)}
    }
    
}



rgparGH <- function(data, g=2, w=NULL,l=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    for (k in 1:g) val[[k]] = ipar(data=data, wt=w[,k],lc=l[k,])
    val$pi = rep(1/g,g)
    return(val)
}


ipar <- function(data, wt=NULL,lc=NULL) {
    if (is.null(wt)) wt = rep(1,nrow(data))
    p = ncol(data)
    val = list()
    val$mu    = lc#rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )
    val$alpha = rep(0,p)#rnorm(p, 0, sqrt(1/nrow(data)) )
    val$sigma = diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(val$sigma)=abs(diag(val$sigma))
    if(p==1){val$sigma=var(data)}
    
    for(i in 1:p){if(val$sigma[i,i]<0.1){val$sigma[i,i]=0.1}}
    if (any(eigen(val$sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    val$cpl   = c(1,-1/2)
    return(val)
}
#####end function for igpar

EMgrstepGH <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL) {
    ##parameter estimtion
    if (is.null(w)) w = try(weightsGH(data=data, gpar=gpar,v=v))
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    
    G= length(gpar$pi);##number clusters
    d= length(gpar[[1]]$mu);##number variables
    Sk = array(0, c(d,d,G) )
    for (k in 1:G ) {
        #	print(gpar[[k]]$cpl)
        gpar[[k]] = updatemaScpl(x=data, par=gpar[[k]], weights=w[,k], invS=NULL, alpha.known=NULL, v=v)
        Sk[,,k]   = gpar[[k]]$sigma
        
        #		print(gpar[[k]]$cpl)
        #		print(c("******"))
    }
    gpar$pi = apply(w,2,mean)
    return(gpar)
}







updatemaScpl <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2GH(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        A = ABC[1]
        B = ABC[2]
        u = (B - abc[,2])*weights
        t = (A*abc[,2]-1)*weights
        T = sum(t)
        
        mu.new    = apply(x, 2, weighted.sum, w=t)/T
        alpha.new = apply(x, 2, weighted.sum, w=u)/T
    } else {
        alpha.new = alpha.known
        mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) - alpha.new/ABC[2]
    }
    alpha.new = alpha.new*v
    
    A = cov.wt(x, wt=abc[,2]*weights, center=mu.new, method="ML")$cov*ABC[2] #returns a list containing weighted covariance matrix and the mean of the data
    #	R = A - outer(alpha.new,alpha.new)*ABC[1]
    #	r = sweep(x, 2, mu.new)
    r = apply(x, 2, weighted.sum, wt=weights)/sumw
    R = A - (outer( r - mu.new, alpha.new) + outer(alpha.new, r - mu.new)) + outer(alpha.new,alpha.new)*ABC[1]
    for(i in 1:ncol(R)){if(R[i,i]<0.00001){R[i,i]=0.00001}}
    par.old = c(  log(par$cpl[1]) , par$cpl[2])
    #a = optim( par.old, loglikgig4, ABC=ABC)
    par.ol = c(exp(par.old[1]), par.old[2])
    a = updateol(ol=par.ol, ABC=ABC, n=2)
    #if ( loglikgig4(par.old, ABC) < loglikgig4(par.old, ABC) ) cpl.new = par$cpl
    #  else cpl.new =  c( rep(exp(a$par[1]), 2), a$par[2])
    
    cpl.new =  c( a[1], a[2])
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=R, cpl=cpl.new )
    return(new.par)
}






gig2GH <- function(x=NULL, par=NULL, invS=NULL) {
    #  returns a matrix with dim length(a) x 3
    d = length(par$mu)
    
    #if (is.null(invS)) invS = ginv(par$sigma)
    alpha = par$alpha
    omega = exp(  log(par$cpl[1])  )
    
    a1 = omega + as.numeric( alpha %*% invS %*% alpha )
    b1 = omega + as.numeric(mahalanobis(x, center=par$mu, cov=invS, inverted=TRUE))
    v1 = par$cpl[2]-d/2
    
    val = gigGH(b=b1,a=a1, v=v1)
    return(val)
}

gigGH <- function(a=NULL,b=NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3 stima le yi
    sab =  sqrt(a*b)
    kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
    kv  = besselK( sab, nu=v, expon.scaled =TRUE)
    kv12 = kv1/kv
    
    sb.a = sqrt(b/a)
    w    = kv12*sb.a
    invw = kv12*1/sb.a - 2*v/b
    logw = log(sb.a) + grad( logbesselKv, x=rep(v,length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    
    
    val =cbind(w,invw,logw)
    return(val)
}





###Function2


##for funtion 2


weightsGH <- function(data=NULL, gpar=NULL, v=1) {
    G = length(gpar$pi)
    if (G > 1) {
        zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
        for (k in 1:G )zlog[,k] =  ddghypGH(x=data, par=gpar[[k]], log=TRUE)
        w = t(apply(zlog, 1, function(z,wt,v) {
            fstar=v*(z + log(wt))-max(v*(z + log(wt)))
            x=exp(fstar);
            
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
    } else w = matrix(1,nrow=nrow(data), ncol=G)
    return(w)
}





ddghypGH <- function(x, par, log=FALSE, invS=NULL) {
    # x  is a n * p matrix
    # generalized hyperbolic distribution
    
    mu    = par$mu
    sigma = par$sigma
    alpha = par$alpha
    
    #check for input error
    if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
    d = length(mu)
    omega = exp(log(par$cpl[1]))
    lambda = par$cpl[2];
    
    
    #distances between points and mu
    if (is.null(invS)) invS = ginv(sigma)
    pa = omega + alpha %*% invS %*% alpha
    mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
    
    kx = sqrt(mx*pa)
    xmu = sweep(x,2,mu)
    
    lvx = matrix(0, nrow=nrow(x), 4)
    lvx[,1] = (lambda - d/2)*log(kx)
    lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
    lvx[,3] = xmu %*% invS %*% alpha
    
    
    lv = numeric(6)
    if(is.nan(log(det( sigma )))){sigma =diag(ncol(mu))}
    
    lv[1] = -1/2*log(det( sigma )) -d/2*(log(2)+log(pi))
    lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
    lv[3] = -lambda/2*( log(1) )
    lv[4] = lambda*log(omega) *0
    lv[5] = (d/2 - lambda)*log( pa )
    
    val = apply(lvx,1,sum) + sum(lv)
    if (!log) val = exp( val )
    
    return(val)
}


##end for function 2


###Function3


llikGH <- function(data,gpar, delta=0) {
    ##log likelyhood estimation
    logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:length(gpar$pi) ) logz[,k] = ddghypGH(x=data, par=gpar[[k]], log=TRUE)
    val = sum(log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
    return(val)
}


weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )


###function 4

MAPGH <- function(data, gpar, label=NULL) {
    w = weightsGH(data=data, gpar=gpar, v=1)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
    z = as.numeric(z)
    return( z)
}

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                            MSGHD                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

mainMSGHD<-function(data=NULL, gpar0, G, n, label  , eps, method ,nr=NULL ) {
    pcol=ncol(data)
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgparMS(data=data, g=G, w=z,l=lc)
        
    }
    else{
        if (is.null(gpar0)) gpar = rmgparMS(g=G,p=ncol(data),data,method=method,nr=nr)
        else gpar  = gpar0}
    
    loglik = numeric(n)
    for (i in 1:3) {
        gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llikMS(data, gpar)
    }
    while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
        i = i+1
        gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llikMS(data, gpar)
    }
    if(i<n){loglik=loglik[-(i+1:n)]}
    BIC=2*loglik[i]-log(nrow(data))*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    z=weightsMS(data=data, gpar= gpar)
    ICL=BIC+sum(log(apply(z,1,max)))
    AIC=2*loglik[i]-2*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    AIC3=2*loglik[i]-3*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    val = list(loglik= loglik, gpar=gpar, z=z, map=MAPMS(data=data, gpar= gpar, label=label), BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3)
    return(val)
}


rmgparMS <- function(g=NULL, p=NULL,data=NULL, method="kmeans",n=10,nr=10) {
    if(g==1){
        # lk=kmeans(data,g)
        lc=as.matrix(t(apply(data,2,mean)))
        l=as.vector(rep(1,nrow(data)))#lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  try({ gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)}, TRUE)
        return(gpar)
    }
    
    else{if(method=="modelBased"){
        lk=gpcm(data,  G=g, mnames=c("VVV"))
        l=lk$map
        lc=lk$gpar[[1]]$mu
        for(i in 2:g){
            lc=rbind(lc,lk$gpar[[i]]$mu)
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)
        # l=BAR(data,l)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        lc=apply(data[l==1,],2,mean)
        for(i in 2:g){
            lc=rbind(lc,apply(data[l==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n)  gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)
        
        #l=BAR(data,l)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        llkO=-Inf
        for(i in 1:nr){
            l=round(runif(nrow(data))*(g-1)+1)
            lc=apply(data[l==1,],2,mean)
            for(i in 2:g){
                lc=rbind(lc,apply(data[l==i,],2,mean))
            }
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gparO  = rgparMS(data=data, g=g, w=z,l=lc)
            
            
            
            loglik = numeric(100)
            for (i in 1:3) {
                gparO = EMgrstepMS(data=data, gpar=gparO, v=1, label = l,it=i)
                loglik[i] = llikMS(data, gparO)
            }
            
            while ( ( getall(loglik[1:i]) > 1) & (i < (100) ) )  {
                i = i+1
                gparO = EMgrstepMS(data=data, gpar=gparO, v=1, label = l,it=i)
                loglik[i] = llikMS(data, gparO)
                
            }
            llk=llikMS(data,gparO)
            if(llk>llkO){
                llkO=llk
                gpar=gparO
                
            }
        }
        return(gpar)
    }
    
    else if(method=="kmedoids"){
        lk=pam(data,g)
        lc=lk$medoids
        l=lk$clustering#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n) gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)}
    
    
    else{ lk=kmeans(data,g)
        lc=lk$centers
        l=lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n)  gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)}}
    
}


##########Inizialization
rgparMS <- function(g=NULL,data, w=NULL,l=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    # l=rmgparMS(g, p,data,method=method)
    for (k in 1:g) val[[k]] = rparMS(data=data,wt=w[,k],k,lc=l)
    val$pi = rep(1/g,g)
    return(val)
}

rparMS <- function(data,wt,k,lc=NULL) {
    par = list()
    p=ncol(data)
    par$mu =lc[k,]
    #par$phi  =  rnorm(p, 0, sqrt(1/nrow(data)) )#rep(1,p);
    par$alpha =rep(0,p)
    
    sigma = cov.wt(data, wt = wt, method="ML")$cov
    diag(sigma)=abs(diag(sigma))
    if(p==1){sigma=var(data)}
    
    for(i in 1:p){if(sigma[i,i]<0.1){sigma[i,i]=0.1}}
    if (any(eigen(sigma)$values <= 0 ) ) par$sigma =  diag(apply(data,2,var))
    par$gam   = eigen( sigma)$vectors
    par$phi  = eigen( sigma)$values
    par$cpl = cbind( rep(1,p), rep(-1/2,p))
    #	par$gam   = diag( rep(1,p) )
    #par$gam   = eigen( cov(matrix(rnorm((p+1)*p), p+1, p)))$vectors
    
    return(par)
}


#########Main

EMgrstepMS <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL,it=NULL) {
    if (is.null(w)) w = weightsMS(data=data, gpar=gpar,v=v)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    
    G= length(gpar$pi);
    for (k in 1:G ) {
        w = weightsMS(data=data, gpar=gpar,v=v)
        if (!is.null(label)) w = combinewk(weights=w, label= label)
        gpar[[k]] = updatemaScplpMS(y=data, par=gpar[[k]], weights=w[,k], alpha.known=NULL, v=v,it=it)
    }
    gpar$pi = apply(w,2,mean)
    return(gpar)
}



weightsMS <- function(data=NULL, gpar=NULL, v=1) {
    G = length(gpar$pi)
    if (G > 1) {
        zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
        for (k in 1:G ) zlog[,k] =  dmsghypMS(y=data, par=gpar[[k]], log=TRUE)
        w = t(apply(zlog, 1, function(z,wt,v) {
            x=exp(v*(z + log(wt)) );
            
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
    } else w = matrix(1,nrow=nrow(data), ncol=G)
    return(w)
}




dmsghypMS <- function(y, par, log=FALSE) {
    # x is a n x p matrix
    x     = y %*% (par$gam)
    mu    = par$mu;
    sigma = par$phi;
    alpha = par$alpha;
    d = length(mu); chi = par$cpl[,1]; psi = par$cpl[,1];
    lambda = par$cpl[,2];
    
    xmu = sweep(x,2,mu,"-")
    
    pa = psi + alpha^2/sigma   # p numeric
    mx = sweep(sweep(xmu^2,2,1/sigma, FUN="*"), 2,chi, "+") # n x p matrix
    kx = sqrt(sweep(mx, 2, pa, "*")) # nxp matrix
    
    lx1 = sweep( sweep(log(mx),2,log(pa),"-"), 2, (lambda - 1/2)/2, "*")
    lx2 = t(apply(kx, 1, function(z,lam=NULL) { log(besselK( z, nu=lambda-1/2, expon.scaled =TRUE)) - z }, lam=lambda ))
    lx3 = sweep(xmu, 2, alpha/sigma, FUN="*")
    
    lv = numeric(3)
    lv1 = -1/2*log( sigma ) -1/2*(log(2)+log(pi)) # d=1
    lv2 = - (log(besselK( sqrt(chi*psi), nu=lambda, expon.scaled =TRUE)) - sqrt(chi*psi) )
    lv3 = lambda/2*(log(psi)-log(chi) )
    if(ncol(y)==1){lx2=t(lx2)}
    lx=lx1 +lx2 + lx3
    #val = apply(lx1 + lx2 + lx3, 1, sum) + sum(lv1 + lv2 + lv3)
    val = apply(lx, 1, sum) + sum(lv1 + lv2 + lv3)
    if (!log) val = exp( val )
    
    return(val)
}




gig2pMS <- function(sr2=NULL, par=NULL) {
    # returns the same as gig
    
    omega = exp(log(par$cpl[,1]) )
    a1 = omega + par$alpha^2/par$phi   # vector
    B1 = sweep(sr2, 2, omega, "+") # matrix
    v1 = par$cpl[,2]-1/2 # vector
    
    val = gigp(a=a1, B=B1, v=v1)
    return(val)
}



updatemaScplpMS <- function(y=NULL, par=NULL, weights=NULL, alpha.known=NULL, v=1,it=NULL) {
    if (is.null(weights)) weights=rep(1,nrow(x))
    
    x = y %*% (par$gam)
    sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    abc = gig2pMS(sr2=sr2, par=par)
    
    #	new.gam= par$gam
     if (it %% 2 ==0) new.gam = updategam2MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    else new.gam = updategam1MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    #new.gam = update.gam1(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    x = y %*% (new.gam)
    #	sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    #print(par$cpl)
    #	abc = gig2p(sr2=sr2, par=par)
    
    
    sumw = sum(weights)
    A = apply(abc$W,   2,weighted.sum, wt=weights)/sumw
    B = apply(abc$invW,2,weighted.sum, wt=weights)/sumw
    C = apply(abc$logW,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        u = sweep(sweep(-abc$invW, 2, B, "+"), 1, weights, "*")
        t = sweep(sweep(abc$invW, 2, A, "*")-1, 1, weights, "*")
        T = apply(t,2,sum)
        
        mu.new    = apply(x*t,2,sum)/T
        alpha.new = apply(x*u,2,sum)/T
        
    } else {
        alpha.new = alpha.known
        mu.new    = apply(x*abc$invW, 2, weighted.sum, wt=weights)/sumw - alpha.new/B
    }
    alpha.new = alpha.new*v
    
    # update for sigma
    Ax = apply(sweep(x, 2, mu.new, "-")^2*abc$invW, 2, weighted.sum, wt=weights)/sumw
    ax = apply(x, 2, weighted.sum, wt=weights)/sumw
    sigma.new = Ax - 2*(ax - mu.new)*alpha.new + alpha.new^2*A
    
    omega=  exp(log(par$cpl[,1]) )
    test  = cbind(omega, lambda=par$cpl[,2], A, B, C)
    cpl.new = t(apply(test, 1, function(z) {
        temp = updateol(ol= z[1:2], ABC=z[3:5], n=2)
        return( c( ( temp[1]), temp[2]) )
    }))
    
    #	new.gam = update.gam1(gam0=par$gam, x=x, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=par$gam, y=y, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=new.gam, x=y %*% (new.gam), sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    
    new.par = list(mu=mu.new, alpha=alpha.new,phi=sigma.new, cpl=cpl.new, gam= new.gam )
    return(new.par)
}


updategam2MS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/sigma, "*" )
    #	u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt, "*")
    F0t = t( x * Bs ) %*% wty
    
    e2  = apply(y^2, 1, sum)*wt
    
    if(ncol(y)==1){F1t = (apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0)}
    else{
        F1t = diag(apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0)}
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
    A0 = t(v)  %*% wty
    
    F2 = F0t - F1t - A0
    z  = svd( F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
    
    if (objgamMS(gam0=gam,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) <  objgamMS(gam0=gam0,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) ) {
        #		print('not good enough')
        gam = gam0
    }
    
    return(gam)
}






updategam1MS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/sigma, "*" )
    u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt, "*")
    F0t = t( x * Bs ) %*% wty
    
    e2  = apply(Bs, 1, max)*wt
    F1t = (cov.wt(y, center=rep(0,ncol(y)), wt=e2, method="ML" )$cov*sum(e2)) %*% (gam0)
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
    A0 = t(v)  %*% wty
    
    F2 = F0t - t(F1t) - A0
    z  = svd( F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
    
    if (objgamMS(gam0=gam,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) <  objgamMS(gam0=gam0,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) ) {
        #		print('not good enough')
        gam = gam0
    }
    
    return(gam)
}




objgamMS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs  = sweep(invW, 2, 1/sigma, "*" )
    wtx = sweep(x, 1, wt, "*")
    F0t = t( x * Bs ) %*% wtx
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
    A0 = t(v)  %*% wtx
    
    val = -1*(tr(F0t ) - 2*tr( A0  ))/sumw
    
    return(val)
}


weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )


################likelihood

llikMS <- function(data, gpar) {
    logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:length(gpar$pi) ) logz[,k] = dmsghypMS(y=data, par=gpar[[k]], log=TRUE)
    val = sum(log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
    return(val)
}


#############MAP

MAPMS <- function(data, gpar, label=NULL) {
    w = weightsMS(data=data, gpar=gpar, v=1)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
    z = as.numeric(z)
    return( z)
}

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                            cMSGHD                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################




maincMSGHD<-function(data=NULL, gpar0=NULL, G, n, label  ,eps, method,nr=NULL ) {
    pcol=ncol(data)
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgparMSr(data=data, g=G, w=z,l=lc)
        
    }
    else{
        if (is.null(gpar0)) gpar = rmgparMSr(g=G,data,method=method,nr=nr)
        else gpar  = gpar0}
    
    loglik = numeric(n)
    for (i in 1:3) {
        gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llikMS(data, gpar)
    }
    while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
        i = i+1
        gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llikMS(data, gpar)
    }
    if(i<n){loglik=loglik[-(i+1:n)]}
    BIC=2*loglik[i]-log(nrow(data))*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    z=weightsMS(data=data, gpar= gpar)
    map=MAPMS(data=data, gpar= gpar, label=label)
    ICL=BIC+sum(log(apply(z,1,max)))
    AIC=2*loglik[i]-2*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    AIC3=2*loglik[i]-3*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    
    val = list(loglik= loglik, gpar=gpar, z=z, map=map, BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3)
    return(val)
    
    
    
    
}




##########Inizialization

rmgparMSr <- function(g=NULL, data=NULL, method="kmeans",n=10,nr=10) {
    if(g==1){
        # lk=kmeans(data,g)
        lc=as.matrix(t(apply(data,2,mean)))
        l=as.vector(rep(1,nrow(data)))#lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  try({ gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)}, TRUE)
        return(gpar)
    }
    else{if(method=="modelBased"){
        lk=gpcm(data,  G=g, mnames=c("VVV"))
        l=lk$map
        lc=lk$gpar[[1]]$mu
        for(i in 2:g){
            lc=rbind(lc,lk$gpar[[i]]$mu)
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,lc)
        
        for (j in 1:n)  gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)
        # l=BAR(data,l)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        lc=apply(data[l==1,],2,mean)
        for(i in 2:g){
            lc=rbind(lc,apply(data[l==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n)  gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)
        
        #l=BAR(data,l)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        llkO=-Inf
        for(i in 1:nr){
            l=round(runif(nrow(data))*(g-1)+1)
            lc=apply(data[l==1,],2,mean)
            for(i in 2:g){
                lc=rbind(lc,apply(data[l==i,],2,mean))
            }
            z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
            gparO  = rgparMSr(data=data, g=g, w=z,l=lc)
            
            
            
            
            
            loglik = numeric(100)
            for (i in 1:3) {
                gparO = EMgrstepMSr(data=data, gpar=gparO, v=1, label = l,  w=z,it=i)
                loglik[i] = llikMS(data, gparO)
            }
            
            while ( ( getall(loglik[1:i]) > 1) & (i < (100) ) )  {
                i = i+1
                gparO = EMgrstepMSr(data=data, gpar=gparO, v=1, label = l,  w=z,it=i)
                loglik[i] = llikMS(data, gparO)
                
            }
            
            llk=llikMS(data,gparO)
            if(llk>llkO){
                llkO=llk
                gpar=gparO
                
            }
        }
        return(gpar)
    }
    
    else if(method=="kmedoids"){
        lk=pam(data,g)
        lc=lk$medoids
        l=lk$clustering#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n) gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)}
    
    
    
    else{ lk=kmeans(data,g)
        lc=lk$centers
        l=lk$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=lc)
        
        for (j in 1:n)  gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z,it=j)
        return(gpar)}}
    
}



####################################change
rgparMSr<- function(g=NULL,data, w=NULL,l=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    # l=rmgparMS(g, p,data,method=method)
    for (k in 1:g) val[[k]] = rparMSr(data=data,wt=w[,k],k,lc=l)
    val$pi = rep(1/g,g)
    return(val)
}

####################################change
rparMSr <- function(data,wt,k,lc=NULL) {
    par = list()
    p=ncol(data)
    par$mu =lc[k,]#rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )#l[k,] #rnorm(p,0,.01);
    #par$phi  =  rnorm(p, 0, sqrt(1/nrow(data)) )#rep(1,p);
    par$alpha =rep(0,p)#   rnorm(p, 0, sqrt(1/nrow(data)) )#rnorm(p,0,.01);
    
    sigma = cov.wt(data, wt = wt, method="ML")$cov#diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(sigma)=abs(diag(sigma))
    if(p==1){sigma=var(data)}
    
    for(i in 1:p){if(sigma[i,i]<0.1){sigma[i,i]=0.1}}
    if (any(eigen(sigma)$values <= 0 ) ) par$sigma =  diag(apply(data,2,var))
    par$gam   = eigen( sigma)$vectors
    par$phi  = eigen( sigma)$values
    par$cpl = cbind( rep(1,p), rep(1,p))
    return(par)
}


#########Main

EMgrstepMSr <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL,it=NULL) {
    if (is.null(w)) w = weightsMS(data=data, gpar=gpar,v=v)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    
    G= length(gpar$pi);
    for (k in 1:G ) {
        w = weightsMS(data=data, gpar=gpar,v=v)
        if (!is.null(label)) w = combinewk(weights=w, label= label)
        
        gpar[[k]] = updatemaScplpMSr(y=data, par=gpar[[k]], weights=w[,k], alpha.known=NULL, v=v,it=it)
    }
    gpar$pi = apply(w,2,mean)
    return(gpar)
}


####################################change
updatemaScplpMSr <- function(y=NULL, par=NULL, weights=NULL, alpha.known=NULL, v=1,it=NULL) {
    if (is.null(weights)) weights=rep(1,nrow(x))
    
    x = y %*% (par$gam)
    sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    abc = gig2pMS(sr2=sr2, par=par)
    
    #	new.gam= par$gam
     if (it %% 2 ==0) new.gam = updategam2MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    else new.gam = updategam1MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    #new.gam = update.gam1(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    x = y %*% (new.gam)
    #	sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    #print(par$cpl)
    #	abc = gig2p(sr2=sr2, par=par)
    
    
    sumw = sum(weights)
    A = apply(abc$W,   2,weighted.sum, wt=weights)/sumw
    B = apply(abc$invW,2,weighted.sum, wt=weights)/sumw
    C = apply(abc$logW,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        u = sweep(sweep(-abc$invW, 2, B, "+"), 1, weights, "*")
        t = sweep(sweep(abc$invW, 2, A, "*")-1, 1, weights, "*")
        T = apply(t,2,sum)
        
        mu.new    = apply(x*t,2,sum)/T
        alpha.new = apply(x*u,2,sum)/T
        
    } else {
        alpha.new = alpha.known
        mu.new    = apply(x*abc$invW, 2, weighted.sum, wt=weights)/sumw - alpha.new/B
    }
    alpha.new = alpha.new*v
    
    # update for sigma
    Ax = apply(sweep(x, 2, mu.new, "-")^2*abc$invW, 2, weighted.sum, wt=weights)/sumw
    ax = apply(x, 2, weighted.sum, wt=weights)/sumw
    sigma.new = Ax - 2*(ax - mu.new)*alpha.new + alpha.new^2*A
    
    omega=  exp(log(par$cpl[,1]) )
    test  = cbind(omega, lambda=par$cpl[,2], A, B, C)
    cpl.new = t(apply(test, 1, function(z) {
        temp = updateol(ol= z[1:2], ABC=z[3:5], n=2)
        return( c( ( temp[1]), temp[2]) )
    }))
    
    #	new.gam = update.gam1(gam0=par$gam, x=x, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=par$gam, y=y, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=new.gam, x=y %*% (new.gam), sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    for(i in 1:nrow(cpl.new)){
        if( cpl.new[i,2]<1){cpl.new[i,2]=1}}
    new.par = list(mu=mu.new, alpha=alpha.new, phi=sigma.new, cpl=cpl.new, gam= new.gam )
    return(new.par)
}

