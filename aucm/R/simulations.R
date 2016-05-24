# simulate data: two covariates.
# should have been named sim.pepe
sim.dat.1=function(n, seed, add.outliers=FALSE,std.dev = 0.2){    
    
    eta=function(x1,x2) 4*x1 - 3*x2 - (x1-x2)^3
    
    # simulate predictors
    sd.= std.dev # most important parameter
    set.seed(seed)
    x=mvtnorm::rmvnorm(n, sigma = matrix(c(1,rep(.9,2),1)*sd.,2,2) )
    
    # add outliers
    if(add.outliers)    {
        mix.p=.05
        n2=ceiling(n*mix.p)
        while(TRUE){
            x.out=mvtnorm::rmvnorm(n*mix.p*10, sigma = matrix(c(2,0,0,2),2,2) )
            x.out=x.out[x.out[,1]-x.out[,2]>0,,drop=FALSE]
            if (nrow(x.out)>=n2) break;
        }        
        x=rbind(x[1:(n-n2),],x.out[1:n2,])    
    }    
    
    x1s=x[,1]; x2s=x[,2]
    dat=data.frame(y=rbern(nrow(x), expit(eta(x1s,x2s))), x1=x1s, x2=x2s, eta=eta(x1s,x2s))
    dat  
}

# this function actually simulates 4 covariates
sim.pepe.2=function(n, seed, add.outliers=FALSE){    
    
    std.dev = 0.2
    eta=function(x1,x2,x3,x4) (4*x1 - 3*x2 - (x1-x2)^3) + x3*x4
    
    # simulate predictors
    sd.= std.dev # most important parameter
    set.seed(seed)
    x=mvtnorm::rmvnorm(n, sigma = matrix(c(1,rep(.9,2),1)*sd.,2,2) )
    
    # add outliers
    if(add.outliers)    {
        mix.p=.05
        n2=ceiling(n*mix.p)
        while(TRUE){
            x.out=mvtnorm::rmvnorm(n*mix.p*10, sigma = matrix(c(2,0,0,2),2,2) )
            x.out=x.out[x.out[,1]-x.out[,2]>0,,drop=FALSE]
            if (nrow(x.out)>=n2) break;
        }        
        x=rbind(x[1:(n-n2),],x.out[1:n2,])
    }    
    
    x1s=x[,1]; x2s=x[,2]
    x3s=rnorm(n); x4s=rnorm(n)
    linear.pred=eta(x1s,x2s,x3s,x4s)
    dat=data.frame(y=rbern(n, expit(linear.pred)), x1=x1s, x2=x2s, x3=x3s, x4=x4s, eta=linear.pred)
    
    dat  
}





sim.easy=function(n, seed, add.outliers=FALSE){    
    
    set.seed(seed)
    n=n/2
    x=mvtnorm::rmvnorm(n, mean=c(0,0))
    z=mvtnorm::rmvnorm(n, mean=c(1,1))
#    plot(rbind(x,z))
#    points(x, col=2)
    x=rbind(z,x)
    x1s=x[,1]; x2s=x[,2]
    dat=data.frame(y=rep(1:0,each=n), x1=x1s, x2=x2s)
    
    dat
    
}


# simulates skin of the orange data
skin.orange = function(n,seed,noise=FALSE, add.outliers=FALSE){
    if (noise & add.outliers) stop ("we don't support both noise and add.outliers at this time")
    
    n=n/2
    d=4
    set.seed(seed)
    x1=mvtnorm::rmvnorm(n, sigma = diag(rep(1,d)) )
    x2=mvtnorm::rmvnorm(n*50, sigma = diag(rep(1,d)) )
    l=apply(x2, 1, function(x) sum(x**2))
    cond=l>9&l<16
    if (sum(cond)<n) stop("not enough x2")
    x2=x2[cond,][1:n,]
    ## graph it, not particulary helpful in 2D, since the skin is in 4D
    #plot(x1[,1], x1[,2], xlab="dim 1", ylab="dim 2")
    #points(x2[,1], x2[,2], col=2)
    dat = data.frame(y=rep(0:1, each=n), x=rbind(x1,x2))
    names(dat)[-1]="x"%+%1:d
    
    if(noise){
        dat=data.frame(dat, x=mvtnorm::rmvnorm(n*2, sigma = diag(rep(1,6)) ))
        names(dat)[-1]="x"%+%1:(ncol(dat)-1)    
    }
    
    if (add.outliers){
    
        mix.p = 0.2
        n2 = ceiling(n * mix.p)
        dat$y[sample(n + 1:n, n2)] = 0
    
        mix.p=.2
        n2=ceiling(n*mix.p)
        dat$y[sample(1:n,n2)]=1
    }
    
    # scramble row order so that it is not all noncases followed by all cases
    dat=dat[sample(1:nrow(dat)),]
    
    dat
}


# simulates skin of the orange data
sim.ring = function(n,seed,add.outliers=FALSE, outliers.perc=NULL){
    
    set.seed(seed)
    
    # simulate data on the unit disc uniformly
    X=NULL
    n1=n/2
    while(TRUE){
        X1=matrix(runif(n1*2, -1, 1), n1, 2)
        indisc=apply(X1,1,function(x) sum(x**2))<1
        X=rbind(X, X1[indisc,])
        if (nrow(X)>n1) break;
    }
    X=X[1:n1,]
    
    # simulate data on a ring uniformly
    X2=NULL
    n1=n/2
    while(TRUE){
        X1=matrix(runif(n1*2, -1.5, 1.5), n1, 2)
        d=apply(X1,1,function(x) sum(x**2))
        inring= d<1.5^2 & d>1
        X2=rbind(X2, X1[inring,])
        if (nrow(X2)>n1) break;
    }
    X2=X2[1:n1,]
    
    dat=rbind(cbind(y=0,X), cbind(y=1,X2))
    
    if (add.outliers) {
            
        X3=NULL
        n1=n/2 * outliers.perc
        while(TRUE){
            X1=matrix(runif(n1*2, -1.5, 1.5), n1, 2)
            d=apply(X1,1,function(x) sum(x**2))
            inring= d<1.3^2 & d>1.1^2
            X3=rbind(X3, X1[inring,])
            if (nrow(X3)>n1) break;
        }
        X3=X3[1:n1,]
    
        dat=rbind(dat, cbind(y=0, X3))
    }
    
    dat=data.frame(dat)
    names(dat)[2]="x1"
    names(dat)[3]="x2"
    
    dat
    
}



# simulate data from Wu and Liu
sim.disc=function(n, seed, add.outliers=FALSE, flip.rate=0.1) {
    
    set.seed(seed)
    
    # simulate data on the unit disc uniformly
    X=NULL
    n1=n*2
    while(TRUE){
        X1=matrix(runif(n1*2, -1, 1), n1, 2)
        indisc=apply(X1,1,function(x) sum(x**2))<1
        X=rbind(X, X1[indisc,])
        if (nrow(X)>n) break;
    }
    X=X[1:n,]
    #plot(X[,1], X[,2])
    
    # let the first and third quadron be case
    y=as.numeric((X[,1]>0 & X[,2]>0) | (X[,1]<0 & X[,2]<0))
    
    if (add.outliers){
        flip=sample(1:n, n*flip.rate)
        y[flip]=1-y[flip]
    }
    
    res=data.frame(y,X)
    names(res)[2:3]=c("x1","x2")
    res
    
}
## test
#dat=sim.disc(n=100, seed=1, add.outliers=TRUE)
#plot(dat$x1, dat$x2, col=ifelse(dat$y==1,"red","black"), pch=19)


# simulate data from Wu and Liu
sim.d20c=function(n, seed, add.outliers=FALSE) {
    sim.disc (n, seed, add.outliers=add.outliers, flip.rate=0.2)
}


sim.MH1 = function(n, seed) {
    set.seed(seed)
    beta=c(1,-1,1,-1)
    x = cbind(x1=rnorm(n, 0,1), x2=rnorm(n, 0,1), x3=rbern(n, 0.5), x4=rbern(n, 0.5))
    y    = rbern(n, expit( x %*% beta) ) 
    data.frame(y=y,x)
} 
#sim.MH1(n,seed)

sim.MH2 = function(n, seed) {
    set.seed(seed)
    beta=c(1,-1,1,-1)
    x = cbind(x1=rnorm(n, 0,1), x2=rnorm(n, 0,1), x3=rnorm(n, 0,1), x4=rnorm(n, 0,1))
    y    = rbern(n, expit( x %*% beta) ) 
    data.frame(y=y,x)
} 
#sim.MH1(n,seed)

# only first two covariates of MH1
sim.MH11 = function(n, seed, beta=c(1,-1), alpha=0) {
    
    # this generates identical datasets to Shuxin's code
    set.seed(seed)
    x = cbind(x1=rnorm(n, 0,1), x2=rnorm(n, 0,1))
    y    = rbern(n, expit( alpha + x %*% beta ) ) 
    data.frame(y=y,x)
    
} 
#sim.MH1(n,seed)



# Y Huang
sim.YH1 = function(n, seed) {
        
    set.seed(seed)
    C = rbern(n, .5)
    n1=sum(C); n0=n-n1
    x = rbind ( cbind(x1=rnorm(n1, 1/2, 1), x2=rnorm(n1, -1/2, 1)), 
            cbind(x1=rnorm(n0, -1/2, 1), x2=rnorm(n0, 1/2, 1)) )
    beta = c(1,-1)
    y    = rbern(n, expit( x %*% beta ) ) 
    data.frame(y=y,x)
    
} 

#dat=sim.YH1(n,seed)
#plot(x2~x1, dat[dat$y==1,])
#plot(x2~x1, dat[dat$y==0,])
#plot(x2~x1, dat)
#calAUC(compute.roc (dat$x1-dat$x2, dat$y)) # compare to pnorm(1)
#calAUC(compute.roc (dat$x1-2*dat$x2, dat$y)) # compare to pnorm(3/sqrt(10))


# case control sampling
sim.YH2 = function(n, seed, mu1=c(0,0), mu0=c(1,1)) {
        
    set.seed(seed)
    n1=n/2; n0=n/2
    x = rbind ( cbind(x1=rnorm(n1, mu1[1], 1), x2=rnorm(n1, mu1[2], 1)), 
            cbind(x1=rnorm(n0, mu0[1], 1), x2=rnorm(n0, mu0[2], 1)) )
    y    = rep(c(1,0), each=n/2) 
    data.frame(y=y,x)
    
} 
#dat=sim.YH2(n=1000,seed=1)
#calAUC(compute.roc (dat$x1+dat$x2, dat$y)) # compare to pnorm(1)
#dat=sim.YH2(n=1000, seed=1, mu1=c(1/2,-1/2), mu0=c(-1/2,1/2))







sim.p1=function(n, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    x0=runif(n, 20,80)
    x0=scale(x0)
    
    x1=rnorm(n)
    x1[I(x0<0)] = x1[I(x0<0)]/2 + x0[I(x0<0)]/2 
    x2=rnorm(n)
    x2[I(x0>0)] = x2[I(x0>0)]/2 + x0[I(x0>0)]/2 
    
    X=cbind(1, x0, x1, x2)
    colnames(X)=c("intercept", "x0", "x1", "x2")
    beta=c(1, 3, 0, 0)
    eta=X %*% beta
    y=rbern(n, expit(eta))
    dat=data.frame(y, X[,-1,drop=FALSE])
    
    dat
}


# not of much use
sim.pauc.old1=function(n, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    x0=runif(n, 20,80)
    x0=scale(x0)
    
    x1=x0+rnorm(n)
    # mean(x1[x0>quantile(x0, 0.8)]) # roughly 1.5
    x1=x1-1.5
    x1.x0 = I(x0>quantile(x0, 0.8)) * x1
    
    X=cbind(1, x0, x1, x1.x0)
    colnames(X)=c("intercept", "x0", "x1", "x1.x0")
    beta=c(1, 1/2, 0, 2)
    eta=X %*% beta
    y=rbern(n, expit(eta))
    dat=data.frame(y, X[,-1,drop=FALSE])
    
    dat
}


sim.mix=function(n, type, seed=NULL, add.outliers=FALSE, plot=FALSE) {
    
    if (type==1) {
        # generate means of the gaussian components
        set.seed(4)
        Sigma=diag(4); Sigma[1,2]<-Sigma[2,1]<-0
        n.centers=8
        centers = mvtnorm::rmvnorm(n.centers, rep(0,4), Sigma)
        
        if (!is.null(seed)) set.seed(seed)    
        dat=NULL
        for (i in 1:n.centers) {
            dat.tmp = mvtnorm::rmvnorm(n/n.centers, centers[i,], 0.7*diag(4))
            y = 1 - i %% 2 
            dat=rbind(dat, data.frame(y=y, dat.tmp))
        }    
        
    } else if (type==2) {
        # generate means of the gaussian components
        set.seed(4)
        Sigma=diag(4); Sigma[1,2]<-Sigma[2,1]<-0
        n.centers=8
        centers = mvtnorm::rmvnorm(n.centers, rep(0,4), Sigma)
        
        if (!is.null(seed)) set.seed(seed)    
        dat=NULL
        for (i in 1:n.centers) {
            dat.tmp = mvtnorm::rmvnorm(n/n.centers, centers[i,], 0.6*diag(4))
            y = 1 - i %% 2 
            dat=rbind(dat, data.frame(y=y, dat.tmp))
        }    
    
    } else stop("type not supported")
    
    
    names(dat)=tolower(names(dat))
    
    if (plot) {
        plot(centers[,1], centers[,2], pch=2, col=rep(1:2,each=n.centers/2), type="p", xlim=c(-3,3), ylim=c(-3,3))
        points(x2~x1, dat, col=dat$y+1)
    }
        
    dat$case = as.factor(ifelse(dat$y==1,"yes","no"))
    dat
    
} 
#dat=sim.mix(200,1, plot=T)



sim.NL=function(n, type, seed=NULL, add.outliers=FALSE) {
    
    if (!is.null(seed)) set.seed(seed)
    
    df=4
    if (type==51 | type==111) {
        # correlated student's t
        sigma=diag(4)
        sigma[1,2]<-sigma[2,1]<-sigma[3,4]<-sigma[4,3]<-.5
        x = mvtnorm::rmvt(n, mean=rep(0,4), sigma=sigma, df = df)
        dimnames(x)=list(NULL, c("x1","x2","x3","x4"))
    } else if (type==52) {
        # correlated lognormal
        sigma=.2*diag(4)
        sigma[1,2]<-sigma[2,1]<-sigma[3,4]<-sigma[4,3]<-.5*sigma[1,1]  
        x = exp(mvtnorm::rmvnorm(n, mean=rep(0,4), sigma=sigma))
        dimnames(x)=list(NULL, c("x1","x2","x3","x4"))
    } else {
        # independent student's t
        x = cbind(x1=rt(n, df), x2=rt(n, df), x3=rt(n, df), x4=rt(n, df))
    }
    #x = cbind(x1=rnorm(n, 0,1), x2=rnorm(n, 0,1), x3=rnorm(n, 0,1), x4=rnorm(n, 0,1))
    
    # 1-4 are for quadratic kernel simulation study
    if (type==1) {
        k=10
        if (!add.outliers) {
            eta = k*(x[,1]*x[,2])
        } else {
            eta = ifelse (x[,1]^2+x[,2]^2<4, .3*(x[,1]*x[,2]) + x[,3] + x[,3]*x[,4], k*(-x[,1]*x[,2]) + x[,3] + x[,3]*x[,4])            
        }
        eta = eta + x[,3] + x[,3]*x[,4]  #x[,3] + x[,4]^2 # the alternative here makes gam better but does not highlight the diff bt raucq and svmq
    } else if (type==2) {
        k=10
        if (!add.outliers) {
            eta = k*(x[,1]*x[,2])
        } else {
            eta = ifelse (x[,1]^2+x[,2]^2<4, .3*(x[,1]^2 + x[,2]^2) + x[,3] + x[,3]*x[,4], k*(x[,1]^2 - x[,2]^2) + x[,3] + x[,3]*x[,4])
        }
        eta = eta + x[,3] + x[,3]*x[,4]  #x[,3] + x[,4]^2 # the alternative here makes gam better but does not highlight the diff bt raucq and svmq
    } else if (type==3) {
        k=10
        if (!add.outliers) {
            eta = k*(x[,1]*x[,2])
        } else {
            # outliers boundary is a vertical line
            eta = ifelse (x[,1]^2+x[,2]^2<6, k*(x[,1]*x[,2]), k*(x[,1]) ) 
        }
        eta = eta + x[,3] + x[,3]*x[,4]  #x[,3] + x[,4]^2 # the alternative here makes gam better but does not highlight the diff bt raucq and svmq
    } else if (type==4) {
        k=10
        if (!add.outliers) {
            eta = k*(x[,1]*x[,2])
        } else {
            # rotate the class boundary by 45 degree, quadatric kernel seems to be able to handle it
            eta = ifelse (x[,1]^2+x[,2]^2<6, k*(x[,1]*x[,2]), k*(x[,1]^2-x[,2]^2) )  
        }
        eta = eta + x[,3] + x[,3]*x[,4]  #x[,3] + x[,4]^2 # the alternative here makes gam better but does not highlight the diff bt raucq and svmq
    
    } else if (type==5 | type==51 | type==52) {    
        # for rbf kernel
        k=10
        v=.5
        
        # for plotting
#        if (!add.outliers) {
#            eta = k*(sin(x[,1]*v*pi)+sin(x[,2]*v*pi))
#        } else {
#            eta = ifelse (x[,1]^2+x[,2]^2<1, k*(sin(x[,1]*v*pi)+sin(x[,2]*v*pi)), k*(cos(x[,1]*v*pi)+cos(x[,2]*v*pi)) )  
#        }
    
        if (!add.outliers) {
            eta =                                             k*(sin(x[,1]*v*pi)+sin(x[,2]*v*pi)+sin(x[,3]*v*pi)+sin(x[,4]*v*pi))
        } else {
            eta = ifelse (x[,1]^2+x[,2]^2+x[,3]^2+x[,4]^2<10, k*(sin(x[,1]*v*pi)+sin(x[,2]*v*pi)+sin(x[,3]*v*pi)+sin(x[,4]*v*pi)), 
                                                              k*(cos(x[,1]*v*pi)+cos(x[,2]*v*pi)+cos(x[,3]*v*pi)+cos(x[,4]*v*pi)) 
                  )  
        }
        
    } else if (type==6) {    
        k=10
        v=.5
        eta = k* (sin(x[,1]*v*pi)+cos(x[,2]*v*pi)+sin(x[,3]*v*pi)*sin(x[,4]*v*pi))
        
    } else if (type==7) {    
        k=10; v=.5
        eta = k* (sin(x[,1]*v*pi)+cos(x[,2]*v*pi)+x[,3]*x[,4])
        
    } else if (type==8) {    
        k=10; v=.5
        eta = k* (sin(x[,1]*v*pi)+exp(x[,2])+5*x[,3]*x[,4] + x[,3]**2 + x[,4]**2)
        #eta = k* (sin(x[,1]*v*pi)+exp(x[,2])+ x[,3]*x[,4] + 0.5*x[,3]**2 + 0.5*x[,4]**2) # linear is pretty good
        
    } else if (type==9) {    
        k=8; v=.5
        eta = k* ( sin(x[,1]*v*pi)+ sin(x[,2]*v*pi)*cos(x[,1]*v*pi) + x[,3]*x[,4] + x[,3]**2 + x[,4]**2)
        
    } else if (type==10) {    
        eta = 5*x[,1]*cos(x[,2]) + x[,3]**3 + 5*x[,3]*x[,4]
        
    } else if (type==11 | type==111) {    
        k=8; v=.5
        eta = k* ( sin(x[,1]*v*pi)+ cos(x[,1]*x[,2]*pi) + 3*x[,3]*x[,4] + x[,3]**2 + x[,4]**2)
                
    } else stop("type not supported: "%+% type)
    y = rbern(n, expit(eta))
    
    x=scale(x)
    dat=data.frame(y=y,x,eta=eta)
    dat$case = as.factor(ifelse(y==1,"yes","no"))
    dat
    
} 
