is.whole <- function(x, tol = .Machine$double.eps^0.5){all(abs(x - round(x)) < tol)}

regularGrid <- function(height, width) {
  if(! is.whole(height) || ! is.whole(width) || length(height) > 1 || length(width) > 1){
    warning("Height and width must be integers")
  } else{
    locations <- matrix(data = 0, nrow = height * width, ncol = 2)
    locations[, 1] <- rep(1:height, width)
    locations[, 2] <- rep(1:width, each = height)
    return( locations )
  }
}

coordinates <- function(locations, indices){
    return(lapply(c(1:nrow(indices)), function(i) locations[indices[i,] == TRUE,]))
}


Kernel<-function(a,tau){return((tau + 1)*(a^tau))}
tildeL<-function(ranks,k, tau, x, n = nrow(ranks)){
    d <- ncol(ranks)
    return(.C("tildeLGeneral", as.integer(as.vector(ranks)),as.integer(d),as.double(k),as.double(tau),
                  as.double(x),as.double(n),result=double(1), PACKAGE = "tailDepFun")$result)
}
deltaL<-function(ranks,k,a,x,tau){tildeL(ranks,k,tau,a*x)/a - tildeL(ranks,k,tau,x)}
rhoL<-function(ranks,k,x,tau,a,r){
    res<-1 - (1/log(r))*log(abs(deltaL(ranks,k,a,r*x,tau)/deltaL(ranks,k,a,x,tau)))
    return(ifelse(res > -0.1, -1, res))
}

alphaL<-function(ranks,k,x,rhot,n=nrow(ranks)){
    d <- ncol(ranks)
    return(.C("alphaGeneral",as.integer(as.vector(ranks)),as.integer(d),as.double(k),as.double(n),
                  as.double(x),as.double(rhot),result=double(1), PACKAGE = "tailDepFun")$result)
}

############################ Brown-Resnick: M-estimator #########################################
tauDerivMatrix <- function(pars){
    C22 <- (-1/(pars[1]^3))*2*((cos(pars[2])^2) + (pars[3]^2)*(sin(pars[2])^2))
    C32 <- (-1/(pars[1]^3))*2*((sin(pars[2])^2) + (pars[3]^2)*(cos(pars[2])^2))
    C42 <- (-1/(pars[1]^3))*(pars[3]^2 - 1)*sin(2*pars[2])
    C3 <- ((pars[3]^2 - 1)/(pars[1]^2))*c(0, sin(2*pars[2]), -sin(2*pars[2]), cos(2*pars[2]))
    C4 <- (pars[3]/(pars[1]^2))*c(0, 2*(sin(pars[2])^2), 2*(cos(pars[2])^2), sin(2*pars[2]))
    result <- cbind(c(1,0,0,0),c(0, C22, C32, C42),C3,C4)
    colnames(result) <- NULL
    return(result)
}

parsTauToBR<-function(tau){
    wortel<-sqrt(1 + (4*(tau[3]^2))/((tau[2]-tau[1])^2))
    rat<-(tau[2]+tau[1])/(tau[1]-tau[2])
    beta<-0.5*atan((2*tau[3])/(tau[2]-tau[1]))
    if(beta >=0){
        cc<-sqrt((rat-wortel)/(rat+wortel))
        rho<-sqrt((cc^2 - 1)/(wortel*(tau[2]-tau[1])))
    } else{
        beta<-0.5*(atan((2*tau[3])/(tau[2]-tau[1]))+pi)
        cc<-1/sqrt((rat-wortel)/(rat+wortel))
        rho<-sqrt((cc^2 - 1)/(wortel*(tau[1]-tau[2])))
    }
    return(c(rho,beta,cc))
}

parsBRtoTau<-function(pars){
    tau11<-(1/(2*(pars[1]^2)))*(1 + pars[3]^2 - (pars[3]^2 - 1)*cos(2*pars[2]))
    tau22<-(1/(2*(pars[1]^2)))*(1 + pars[3]^2 + (pars[3]^2 - 1)*cos(2*pars[2]))
    tau12<-(((pars[3]^2 - 1)/(2*(pars[1]^2)))*sin(2*pars[2]))
    return(c(tau11,tau22,tau12))
}

tailBRInt <- function(loc, pars){ #pars = (alpha,rho) or pars = (alpha,tau)
    locdiff<-loc[2,]-loc[1,]
    if(length(pars) == 2){
        a <- sqrt(2)*((sqrt(locdiff[1]^2 + locdiff[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        a <- sqrt(2)*((t(locdiff) %*% tau %*% locdiff)^(pars[1]/4))
    }
    logres<-(a^2 + stats::pnorm((-3*a)/2,log.p=TRUE) - log(3))
    return(as.vector(stats::pnorm(a/2) + exp(logres)))
}

MestimatorMinimizeBR<-function(pars,totlist,pairs,w){
    res <- unlist(lapply(pairs,tailBRInt,pars)) - totlist
    return(t(res) %*% w %*% res)
}

ellBR<-function(t,loc,pars){
    if((loc[1]==0) && (loc[2]==0)){
        return(pmax(t[1],t[2]))
    } else{
        locdiff<-loc[3:4]-loc[1:2]
        if(length(pars) == 2){
            a <- sqrt(2)*((sqrt(loc[1]^2 + loc[2]^2)/pars[2])^(pars[1]/2))
        } else{
            tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
            a <- sqrt(2)*((t(loc) %*% tau %*% loc)^(pars[1]/4))
        }
        res<-.C("ellsmith2",as.double(t),as.double(a),result=double(1),
                PACKAGE = "tailDepFun")$result
        return(res)
    }
}

ellBR4d<-function(t,rlist,gamlist,d){
    result<-vector(length=d)
    for(i in 1:d){
        const<-(gamlist[[i]]/2 + log(t[i]/t[-i])/gamlist[[i]])
        result[i]<-t[i]*pmvnorm(lower=rep(-Inf,d-1),upper=const,corr=rlist[[i]],algorithm=TVPACK)[1]
    }
    return(sum(result))
}

gfunc<-function(pars,loc){
    if(length(pars) == 2){
        a2 <- (sqrt(loc[1]^2 + loc[2]^2)/pars[2])^pars[1]
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        a2 <- (t(loc) %*% tau %*% loc)^(pars[1]/2)
    }
    return(a2)
}

corrfunc<-function(pars,loc,d,i){
    result<-matrix(1,nrow=(d-1),ncol=(d-1))
    for(j in 1:(d-2)){
        for(k in (j+1):(d-1)){
            lo<-loc[,-i]
            result[j,k]<-result[k,j]<-((gfunc(pars,loc[,i]-lo[,j]) + gfunc(pars,loc[,i]-lo[,k]) -
                                            gfunc(pars,lo[,j]-lo[,k]))/(2*sqrt(gfunc(pars,loc[,i]-lo[,j])*gfunc(pars,loc[,i]-lo[,k]))))
        }
    }
    return(result)
}

tailBR <- function(t, loc, pars, d){
    if(d==2){
        return(ellBR(t,loc[,2]-loc[,1],pars))
    } else{
        rlist<-gamlist<-vector('list',length=d)
        for(i in 1:d){
            rlist[[i]]<-corrfunc(pars,loc,d,i)
            gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,loc[,i]-loc[,j])))
        }
        return(ellBR4d(t,rlist,gamlist,d))
    }
}

ellBR4dv2<-function(t,rlist,gamlist,d,ind){
    return(2*t[ind]*ellBR4d(t,rlist,gamlist,d))
}
ellBR4dv3<-function(t,coord,pars,ind){
    return(ellBR(c(pmax(t[ind[1]],t[ind[2]]),t[ind[3]]),coord[,2]-coord[,1],pars))
}
ellBR4dv4<-function(t,coord,pars){
    return(4*t[1]*t[2]*ellBR(t,coord[,2]-coord[,1],pars))
}

IBR<-function(a){ #  This function can be implemented in C as well.
    return(stats::pnorm(a/2) + exp(a*a)*(1/3)*stats::pnorm(-(3/2)*a))
}
I2BR<-function(x,a){  #  This function can be implemented in C as well.
    return(0.5*stats::pnorm(a/2 - log(x)/a) + x*stats::pnorm(a/2 + log(x)/a) + 0.5*x*x*exp(a*a)*stats::pnorm(-(3/2)*a - log(x)/a))
}
I3BR<-function(x,a){
    res<-.C("I3smith2",as.double(x),as.double(a),result=double(1), PACKAGE = "tailDepFun")$result
    return(res)
}

I1<-function(u,a,rlist,gamlist,d,ind){ return(I3BR(u[ind],a)*ellBR4d(u,rlist,gamlist,d))}
I1v2<-function(u,a,pars,coordt,ind,ind2){ return(I3BR(u[ind],a)*ellBR4dv3(u,coordt,pars,ind2))}
I1v3<-function(u,a,rlist,gamlist,d){return(I3BR(u[1],a)*ellBR4d(c(u[1],u[3],u[2]),rlist,gamlist,d))}
I1v4<-function(u,a,rlist,gamlist,d){return(I3BR(u[1],a)*ellBR4d(c(u[2],u[1],u[3]),rlist,gamlist,d))}
I3<-function(u,coord,pars,coordt){
    loc1<-coord[,2]-coord[,1]
    loc2<-coord[,4]-coord[,3]
    if(length(pars)==2){
        anr1 <- sqrt(2)*((sqrt(loc1[1]^2 + loc1[2]^2)/pars[2])^(pars[1]/2))
        anr2 <- sqrt(2)*((sqrt(loc2[1]^2 + loc2[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(pars[1]/4))
        anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(pars[1]/4))
    }
    return(I3BR(u[1],anr1)*I3BR(u[2],anr2)*ellBR(u,coordt[,2]-coordt[,1],pars))
}

####### case 1: locations s,t,u,v all different
intBRcase1<-function(coord,pars,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
    loc1<-coord[,2]-coord[,1]
    loc2<-coord[,4]-coord[,3]
    if(length(pars)==2){
        anr1 <- sqrt(2)*((sqrt(loc1[1]^2 + loc1[2]^2)/pars[2])^(pars[1]/2))
        anr2 <- sqrt(2)*((sqrt(loc2[1]^2 + loc2[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(pars[1]/4))
        anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(pars[1]/4))
    }
    ########### T1
    d<-4
    rlist<-gamlist<-vector('list',length=d)
    for(i in 1:d){
        rlist[[i]]<-corrfunc(pars,coord,d,i)
        gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coord[,i]-coord[,j])))
    }
    T1<-(IBR(anr1) + IBR(anr2) - adaptIntegrate(ellBR4d,lowerLimit=c(0,0,0,0),upperLimit=c(1,1,1,1),rlist=rlist,gamlist=gamlist,d=4,tol=Tol)$integral)
    ############# T2
    d<-3
    rlistpt1<-rlistpt2<-gamlistpt1<-gamlistpt2<-vector('list',length=d)
    coordt1<-coord[,-4]
    coordt2<-coord[,-3]
    for(i in 1:d){
        rlistpt1[[i]]<-corrfunc(pars,coordt1,d,i)
        gamlistpt1[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt1[,i]-coordt1[,j])))
        rlistpt2[[i]]<-corrfunc(pars,coordt2,d,i)
        gamlistpt2[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt2[,i]-coordt2[,j])))
    }
    T2<-(2*IBR(anr1)*(I2BR(1,anr2) -0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlistpt1,gamlist=gamlistpt1,d=3,ind=3,tol=Tol)$integral -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlistpt2,gamlist=gamlistpt2,d=3,ind=3,tol=Tol)$integral)
    ############# T3
    d<-3
    rlistpt1<-rlistpt2<-gamlistpt1<-gamlistpt2<-vector('list',length=3)
    coordt1<-coord[,-2]
    coordt2<-coord[,-1]
    for(i in 1:d){
        rlistpt1[[i]]<-corrfunc(pars,coordt1,d,i)
        gamlistpt1[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt1[,i]-coordt1[,j])))
        rlistpt2[[i]]<-corrfunc(pars,coordt2,d,i)
        gamlistpt2[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt2[,i]-coordt2[,j])))
    }
    T3<-(2*IBR(anr2)*(I2BR(1,anr1) -0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1), a=anr1,rlist=rlistpt1,gamlist=gamlistpt1,d=3,ind=1,tol=Tol)$integral -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1), a=anr1,rlist=rlistpt2,gamlist=gamlistpt2,d=3,ind=1,tol=Tol)$integral)
    ################ T4
    T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,3)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,4)])$integral)
    ################# T5
    T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,4)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,3)])$integral)
    return(T1 - T2 - T3 + T4 + T5)
}

## case 2: t = u
intBRcase2<-function(coord,pars,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
    loc1<-coord[,2]-coord[,1]
    loc2<-coord[,4]-coord[,3]
    if(length(pars)==2){
        anr1 <- sqrt(2)*((sqrt(loc1[1]^2 + loc1[2]^2)/pars[2])^(pars[1]/2))
        anr2 <- sqrt(2)*((sqrt(loc2[1]^2 + loc2[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(pars[1]/4))
        anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(pars[1]/4))
    }
    ########### T1
    d<-3
    rlist<-gamlist<-vector('list',length=d)
    coordt<-coord[,-3]
    for(i in 1:d){
        rlist[[i]]<-corrfunc(pars,coordt,d,i)
        gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt[,i]-coordt[,j])))
    }
    T1<-(IBR(anr1) + IBR(anr2) -
             adaptIntegrate(ellBR4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),rlist=rlist,gamlist=gamlist,d=3,ind=2,tol=Tol)$integral)
    ############# T2
    T2<-(2*IBR(anr1)*(I2BR(1,anr2) -0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,pars=pars,coordt=coord[,1:2],ind=3,ind2=c(2,3,1),tol=Tol)$integral -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlist,gamlist=gamlist,d=3,ind=3,tol=Tol)$integral)
    ############# T3
    T3<-(2*IBR(anr2)*(I2BR(1,anr1) -0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1, rlist=rlist,gamlist=gamlist,d=3,ind=1,tol=Tol)$integral -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,pars=pars,coordt=coord[,3:4],ind=1,ind2=c(1,2,3),tol=Tol)$integral)
    ################ T4
    T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,3)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,4)])$integral)
    ################# T5
    T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,4)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,3)])$integral)
    return(T1 - T2 - T3 + T4 + T5)
}

## case 3: s = u, t neq v
intBRcase3<-function(coord,pars,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
    loc1<-coord[,2]-coord[,1]
    loc2<-coord[,4]-coord[,3]
    if(length(pars)==2){
        anr1 <- sqrt(2)*((sqrt(loc1[1]^2 + loc1[2]^2)/pars[2])^(pars[1]/2))
        anr2 <- sqrt(2)*((sqrt(loc2[1]^2 + loc2[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(pars[1]/4))
        anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(pars[1]/4))
    }
    ########### T1
    d<-3
    rlist<-gamlist<-vector('list',length=d)
    coordt<-coord[,-3]
    for(i in 1:d){
        rlist[[i]]<-corrfunc(pars,coordt,d,i)
        gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt[,i]-coordt[,j])))
    }
    T1<-(IBR(anr1)+IBR(anr2)-adaptIntegrate(ellBR4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),rlist=rlist,gamlist=gamlist,d=3,ind=1,tol=Tol)$integral)
    ############# T2
    T2<-(2*IBR(anr1)*(I2BR(1,anr2) - 0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,pars=pars,coordt=coord[,1:2],ind=3,ind2=c(1,3,2),tol=Tol)$integral -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlist,gamlist=gamlist,d=3,ind=3,tol=Tol)$integral)
    ############# T3
    T3<-(2*IBR(anr2)*(I2BR(1,anr1) -0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,pars=pars,coordt=coord[,c(1,4)],ind=1,ind2=c(1,2,3),tol=Tol)$integral -
             adaptIntegrate(I1v4,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,rlist=rlist,gamlist=gamlist,d=3,tol=Tol)$integral)
    ################ T4
    T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,3)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,4)])$integral)
    ################# T5
    T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,4)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,3)])$integral)
    return(T1 - T2 - T3 + T4 + T5)
}

## case 4: s neq u, t = v
intBRcase4<-function(coord,pars,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
    loc1<-coord[,2]-coord[,1]
    loc2<-coord[,4]-coord[,3]
    if(length(pars)==2){
        anr1 <- sqrt(2)*((sqrt(loc1[1]^2 + loc1[2]^2)/pars[2])^(pars[1]/2))
        anr2 <- sqrt(2)*((sqrt(loc2[1]^2 + loc2[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        anr1 <- sqrt(2)*((t(loc1) %*% tau %*% loc1)^(pars[1]/4))
        anr2 <- sqrt(2)*((t(loc2) %*% tau %*% loc2)^(pars[1]/4))
    }
    ########### T1
    d<-3
    rlist<-gamlist<-vector('list',length=d)
    coordt<-coord[,-4]
    for(i in 1:d){
        rlist[[i]]<-corrfunc(pars,coordt,d,i)
        gamlist[[i]]<-sapply(c(1:d)[-i], function(j) sqrt(2*gfunc(pars,coordt[,i]-coordt[,j])))
    }
    T1<-(IBR(anr1)+IBR(anr2)-adaptIntegrate(ellBR4dv2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),rlist=rlist,gamlist=gamlist,d=3,ind=2,tol=Tol)$integral)
    ############# T2
    T2<-(2*IBR(anr1)*(I2BR(1,anr2) - 0.5) + 2*I2BR(1,anr2) - 2*IBR(anr2) -
             adaptIntegrate(I1,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,rlist=rlist,gamlist=gamlist,d=3,ind=3,tol=Tol)$integral -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr2,pars=pars,coordt=coord[,c(1,2)],ind=3,ind2=c(2,3,1),tol=Tol)$integral)
    ############# T3
    T3<-(2*IBR(anr2)*(I2BR(1,anr1) - 0.5) + 2*I2BR(1,anr1) - 2*IBR(anr1) -
             adaptIntegrate(I1v3,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1, rlist=rlist,gamlist=gamlist,d=3,tol=Tol)$integral -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=anr1,pars=pars,coordt=coord[,c(2,3)],ind=1,ind2=c(1,3,2),tol=Tol)$integral)
    ################ T4
    T4<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,3)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,4)])$integral)
    ################# T5
    T5<-((IBR(anr1) - I2BR(1,anr1))*(1 - 2*I2BR(1,anr2)) + (IBR(anr2) -  I2BR(1,anr2))*(1 - 2*I2BR(1,anr1)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,4)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,3)])$integral)
    return(T1 - T2 - T3 + T4 + T5)
}

## case 5: s = u, t = v
intBRcase5<-function(coord,pars,Tol){ #coord=(s,t,u,v) where s,t,u,v, are 2-dim columnvectors
    loc<-coord[,2]-coord[,1]
    if(length(pars)==2){
        a <- sqrt(2)*((sqrt(loc[1]^2 + loc[2]^2)/pars[2])^(pars[1]/2))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        a <- sqrt(2)*((t(loc) %*% tau %*% loc)^(pars[1]/4))
    }
    ########### T1
    coordt<-coord[,1:2]
    T1<-(2*IBR(a)  - adaptIntegrate(ellBR4dv4,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coordt,pars=pars,tol=Tol)$integral)
    ############# T2
    T2<-(2*IBR(a)*(I2BR(1,a) - 0.5) + 2*I2BR(1,a) - 2*IBR(a) -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,pars=pars,
                            coordt=coordt,ind=3,ind2=c(1,3,2),tol=Tol)$integral -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,pars=pars,
                            coordt=coordt,ind=3,ind2=c(2,3,1),tol=Tol)$integral)
    ############# T3
    T3<-(2*IBR(a)*(I2BR(1,a) - 0.5) + 2*I2BR(1,a) - 2*IBR(a) -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,pars=pars,
                            coordt=coordt,ind=1,ind2=c(1,2,3),tol=Tol)$integral -
             adaptIntegrate(I1v2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),a=a,pars=pars,
                            coordt=coordt,ind=1,ind2=c(1,3,2),tol=Tol)$integral)
    ################ T4
    T4<-((IBR(a) - I2BR(1,a))*(1 - 2*I2BR(1,a)) + (IBR(a) -  I2BR(1,a))*(1 - 2*I2BR(1,a)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,3)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,4)])$integral)
    ################# T5
    T5<-((IBR(a) - I2BR(1,a))*(1 - 2*I2BR(1,a)) + (IBR(a) -  I2BR(1,a))*(1 - 2*I2BR(1,a)) -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(1,4)])$integral -
             adaptIntegrate(I3,lowerLimit=c(0,0),upperLimit=c(1,1),coord=coord,pars=pars,coordt=coord[,c(2,3)])$integral)
    return(T1 - T2 - T3 + T4 + T5)
}

MestimatorpsidotBR<-function(pars,pairs){ #nrow(pairs) = number of pairs, ncol(pairs)=4
    npairs<-length(pairs)
    loc<-matrix(unlist(lapply(pairs,function(i) i[2,]-i[1,])),ncol=2,byrow=TRUE)
    psi<-matrix(,nrow=npairs,ncol=4)
    deriv<-vector(length=4)
    tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
    for(i in 1:npairs){
        a <- sqrt(2)*((t(loc[i,]) %*% tau %*% loc[i,])^(pars[1]/4))
        temp<-sqrt(2)*(pars[1]/4)*((t(loc[i,]) %*% tau %*% loc[i,])^(pars[1]/4 - 1))
        deriv[1]<-(a*log((t(loc[i,]) %*% tau %*% loc[i,])^(1/4)))
        deriv[2]<-temp*(loc[i,1]^2)
        deriv[3]<-temp*(loc[i,2]^2)
        deriv[4]<-temp*(2*loc[i,1]*loc[i,2])
        psi[i,]<-deriv*((1/2)*stats::dnorm(a/2) - (1/2)*exp(a^2)*stats::dnorm(-(3/2)*a) + (2/3)*a*exp(a^2)*stats::pnorm(-(3/2)*a))
    }
    return(psi)
}

MestimatorpsidotBRIso<-function(pars,pairs){ #nrow(pairs) = number of pairs, ncol(pairs)=4
    npairs<-length(pairs)
    loc<-matrix(unlist(lapply(pairs,function(i) i[2,]-i[1,])),ncol=2,byrow=TRUE)
    psi<-matrix(,nrow=npairs,ncol=2)
    deriv<-vector(length=2)
    for(i in 1:npairs){
        temp <- sqrt(loc[i,1]^2 + loc[i,2]^2)
        a <- sqrt(2)*((temp/pars[2])^(pars[1]/2))
        deriv[1]<-(a*log(temp/pars[2])/2)
        deriv[2]<-(-pars[1]/sqrt(2))*((temp/pars[2])^(pars[1]/2 - 1))*(temp/(pars[2]^2))
        psi[i,]<-deriv*((1/2)*stats::dnorm(a/2) - (1/2)*exp(a^2)*stats::dnorm(-(3/2)*a) + (2/3)*a*exp(a^2)*stats::pnorm(-(3/2)*a))
    }
    return(psi)
}

MestimatorAsymVarBR<-function(pairs,pars,Tol){
    npairs<-length(pairs)
    BRMat<-matrix(,nrow=npairs,ncol=npairs)
    for(i in 1:npairs){
        for(j in i:npairs){
            temp<-cbind(t(pairs[[i]]),t(pairs[[j]]))
            if(all(temp[,1]==temp[,4]) && (! all(temp[,2]==temp[,3]))){
                BRMat[i,j]<-intBRcase2(cbind(temp[,3:4],temp[,1:2]),pars,Tol=Tol)
            }else if(all(temp[,1]==temp[,3]) && all(temp[,2]==temp[,4])){
                BRMat[i,j]<-intBRcase5(temp,pars,Tol=Tol)
            }else if(all(temp[,2]==temp[,4]) && (! all(temp[,1]==temp[,3]))){
                BRMat[i,j]<-intBRcase4(temp,pars,Tol=Tol)
            }else if(all(temp[,1]==temp[,3]) && (! all(temp[,2]==temp[,4]))){
                BRMat[i,j]<-intBRcase3(temp,pars,Tol=Tol)
            }else if(all(temp[,2]==temp[,3]) && (! all(temp[,1]==temp[,4]))){
                BRMat[i,j]<-intBRcase2(temp,pars,Tol=Tol)
            }else if((! all(temp[,1]==temp[,3])) && (! all(temp[,2]==temp[,4]))){
                BRMat[i,j]<-intBRcase1(temp,pars,Tol=Tol)
            }
            BRMat[j,i]<-BRMat[i,j]
        }
    }
    return(BRMat)
}
############################ Brown-Resnick: WLS estimator #########################################
gammaf<-function(pars,loc){ #pars=(alpha,rho) or (alpha,tau11,tau22,tau12)
    if(length(pars) == 2){
        return((sqrt(loc[1]^2 + loc[2]^2)/pars[2])^(pars[1]))
    } else{
        tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
        return((t(loc) %*% tau %*% loc)^(pars[1]/2))
    }
}
corrfuncEC<-function(pars,loc,d,i){
    result<-matrix(1,nrow=(d-1),ncol=(d-1))
    for(j in 1:(d-2)){
        for(k in (j+1):(d-1)){
            lo<-loc[-i,]
            result[j,k]<-result[k,j]<-((gammaf(pars,loc[i,]-lo[j,]) + gammaf(pars,loc[i,]-lo[k,]) -
             gammaf(pars,lo[j,]-lo[k,]))/(2*sqrt(gammaf(pars,loc[i,]-lo[j,])*gammaf(pars,loc[i,]-lo[k,]))))
        }
    }
    return(result)
}

ecBR2<-function(loc,pars){return(2*stats::pnorm(sqrt(gammaf(pars,(loc[1,]-loc[2,]))/2)))}

fell<-function(loc,pars,d){ # pars = c(alpha,rho)
    if(d == 2){
        result <- stats::pnorm(sqrt(gammaf(pars,(loc[1,]-loc[2,]))/2))
        return(rep(result,2)) #sum of result gives ell_J
    } else{
        rlist <- lapply(c(1:d), function(i) corrfuncEC(pars,loc,d,i))
        gamlist <- lapply(c(1:d), function(i) sapply(c(1:d)[-i],function(j) gammaf(pars,(loc[i,]-loc[j,]))))
        alg <- ifelse(d == 3 | d == 4, "TVPACK", "GenzBretz")
        result <- sapply(c(1:d), function(i) pmvnorm(lower=rep(-Inf,d-1),upper=sqrt(gamlist[[i]]/2),
                                                     corr=rlist[[i]],algorithm = alg)[1])
        return(result)  #sum of result gives ell_J
    }
}


AsymBREC <- function(J1,J2,pars){
    d <- nrow(J1)
    ellderJ1 <- fell(J1,pars,d)
    ellderJ2 <- fell(J2,pars,d)
    ellJ1 <- sum(ellderJ1)
    ellJ2 <- sum(ellderJ2)
    J12<- unique(rbind(J1,J2))
    T1 <- ellJ1 + ellJ2 - sum(fell(J12,pars,nrow(J12)))
    J1pl <- lapply(c(1:d), function(i) unique(rbind(J1,J2[i,])))
    J2pl <- lapply(c(1:d), function(i) unique(rbind(J2,J1[i,])))
    T2 <- sum(sapply(c(1:d), function(i) ellderJ1[i]*(1 + ellJ2 - sum(fell(J2pl[[i]],pars,nrow(J2pl[[i]]))))))
    T3 <- sum(sapply(c(1:d), function(i) ellderJ2[i]*(1 + ellJ1 - sum(fell(J1pl[[i]],pars,nrow(J1pl[[i]]))))))
    Jsingle <- lapply(c(1:d), function(i) lapply(c(1:d), function(j) unique(rbind(J1[i,],J2[j,]))))
    card <- lapply(c(1:d), function(i) lapply(c(1:d), function(j) nrow(Jsingle[[i]][[j]])))
    T4 <- sum(sapply(c(1:d), function(i) sapply(c(1:d), function(j){
        ellj1j2<-ifelse(card[[i]][[j]] > 1, sum(fell(Jsingle[[i]][[j]],pars, card[[i]][[j]])), 1)
        return(ellderJ1[i]*ellderJ2[j]*(2 - ellj1j2))})))
    return(T1 - T2 - T3 + T4)
}

AsymBRDiagEC <- function(J1, pars){
    d <- nrow(J1)
    ellderJ1 <- fell(J1, pars,d)
    T1 <- sum(ellderJ1)
    Jsingle <- lapply(c(1:d), function(i) lapply(c(1:d), function(j) unique(rbind(J1[i,],J1[j,]))))
    card <- lapply(c(1:d), function(i) lapply(c(1:d), function(j) nrow(Jsingle[[i]][[j]])))
    T4 <- sum(sapply(c(1:d), function(i) sapply(c(1:d), function(j){
        ellj1j2<-ifelse(card[[i]][[j]] > 1, sum(fell(Jsingle[[i]][[j]],pars,card[[i]][[j]])), 1)
        return(ellderJ1[i]*ellderJ1[j]*(2 - ellj1j2))})))
    return(T4 - T1)
}

phidotBREC<-function(J,pars){ #nrow(pairs) = number of pairs, ncol(pairs)=4
    npairs<-length(J)
    loc<-matrix(unlist(lapply(J,function(i) i[2,]-i[1,])),ncol=2,byrow=TRUE)
    psi<-matrix(,nrow=npairs,ncol=4)
    deriv<-vector(length=4)
    tau <- rbind(c(pars[2], pars[4]), c(pars[4], pars[3]))
    for(i in 1:npairs){
        a <- sqrt(2)*((t(loc[i,]) %*% tau %*% loc[i,])^(pars[1]/4))
        temp<-sqrt(2)*(pars[1]/4)*((t(loc[i,]) %*% tau %*% loc[i,])^(pars[1]/4 - 1))
        deriv[1]<-(a*log((t(loc[i,]) %*% tau %*% loc[i,])^(1/4)))
        deriv[2]<-temp*(loc[i,1]^2)
        deriv[3]<-temp*(loc[i,2]^2)
        deriv[4]<-temp*(2*loc[i,1]*loc[i,2])
        psi[i,]<-deriv*stats::dnorm(a/2)
    }
    return(psi)
}

phidotBRECIso<-function(J,pars){ #nrow(pairs) = number of pairs, ncol(pairs)=4
    npairs<-length(J)
    loc<-matrix(unlist(lapply(J,function(i) i[2,]-i[1,])),ncol=2,byrow=TRUE)
    psi<-matrix(,nrow=npairs,ncol=2)
    deriv<-vector(length=2)
    for(i in 1:npairs){
        temp <- sqrt(loc[i,1]^2 + loc[i,2]^2)
        a <- sqrt(2)*((temp/pars[2])^(pars[1]/2))
        deriv[1]<-(a*log(temp/pars[2])/2)
        deriv[2]<-(-pars[1]/sqrt(2))*((temp/pars[2])^(pars[1]/2 - 1))*(temp/(pars[2]^2))
        psi[i,]<-deriv*stats::dnorm(a/2)
    }
    return(psi)
}


WLSAsymVarBR <- function(J, pars){
    q <- length(J)
    resmat <- tempdiag <- diag(sapply(c(1:q), function(j) AsymBRDiagEC(J[[j]],pars)))
    resmat[lower.tri(resmat)] <- unlist(sapply(c(1:(q-1)),
                 function(i) sapply(c((i+1):q), function(j) AsymBREC(J[[i]],J[[j]],pars))))
    M<-t(resmat) + resmat - tempdiag
    return(M)
}

WLSminimizeBR<-function(pars,J,totlist,w){
    if(pars[1] <= 0.01 || pars[1] >= 1.99){
        return(10^6)
    } else{
        res<-unlist(lapply(J, function(i) ecBR2(i,pars))) - totlist
        return(t(res) %*% w %*% res)
    }
}

WLSminimizeBRcu<-function(pars,J,totlist){
    if(pars[1] <= 0.01 || pars[1] >= 1.99){
        return(10^6)
    } else{
        res<-unlist(lapply(J, function(i) ecBR2(i,pars))) - totlist
        Sigma <- WLSAsymVarBR(J,pars)
        return(t(res) %*% solve(Sigma) %*% res)
    }
}


##################### Max-linear: Mestimator ################################
factorfunc<-function(Bpars){
    maxCol <- apply(Bpars,2,max)
    minCol <- apply(Bpars,2,min)
    result <- sapply(which(maxCol>0), function(i) (minCol[i]^2)/(6*maxCol[i]) + maxCol[i]/2)
    return(sum(result))
}

MestimatorMinimizeML<-function(pars,totlist,indices,Bmatrix){
    if(any(pars < 0.01) || any(pars > 0.99)){
        return(10^6)
    } else{
        Bpars <- Bmatrix(pars)
        factor <- lapply(c(1:nrow(indices)), function(i) factorfunc(Bpars[indices[i,]==TRUE,]))
        res <- unlist(factor) - totlist
        return((t(res) %*% res))
    }
}


###################### Max-linear: WLS estimator ##############################
fecML<-function(th, cst){return(sum(apply(th, 2, function(l) max(l*cst))))}
ecderML<-function(th, cst, d){
    return(sapply(c(1:d), function(k) sum(apply(th, 2, function(col) col[k]*(col[k]*cst[k]==max(col*cst))))))
}
AsymML <- function(C1,C2,Bpars,d){
    ec1 <- fecML(Bpars, C1)
    ec2 <- fecML(Bpars, C2)
    der1 <- ecderML(Bpars,C1,d)
    der2 <- ecderML(Bpars,C2,d)
    T1 <- ec1 + ec2 - fecML(Bpars, pmax(C1,C2))
    T2 <- sum(sapply(c(1:d), function(k) {
        Ctemp <- rep(0,d)
        Ctemp[k] <- C2[k]
        der2[k]*(ec1 + C2[k] - fecML(Bpars,pmax(C1,Ctemp)))
    }))
    T3 <- sum(sapply(c(1:d), function(k) {
        Ctemp <- rep(0,d)
        Ctemp[k] <- C1[k]
        der1[k]*(ec2 + C1[k] - fecML(Bpars,pmax(Ctemp,C2)))
    }))
    T4 <- sum(sapply(c(1:d), function(k) sapply(c(1:d), function(l){
        Ctemp1 <- Ctemp2 <- rep(0,d)
        Ctemp1[k] <- C1[k]
        Ctemp2[l] <- C2[l]
        der1[k]*der2[l]*(C1[k] + C2[l] - fecML(Bpars,pmax(Ctemp1,Ctemp2)))
    })))
    return(T1 - T2 - T3 + T4)
}

AsymVarMLwls <- function(Bpars,indices,d){
    q <- nrow(indices)
    resmat <- matrix(0,nrow=q,ncol=q)
    resmat[lower.tri(resmat, diag = TRUE)] <- unlist(sapply(c(1:q), function(i) sapply(c(i:q),
                    function(j) AsymML(indices[i,],indices[j,],Bpars,d))))
    final <- t(resmat) + resmat - diag(diag(resmat))
    return(final)
}

phidotBasic <- function(indices, pars){
    q <- nrow(indices)
    Bpars <- cbind(pars,1-pars)
    result <- sapply(c(1:q), function(j){
        maxj <- apply(Bpars, 2, function(l) max(l*indices[j,]))
        sapply(seq(1:length(pars)), function(l){
            return(indices[j,l]*((Bpars[l,1]*indices[j,l] == maxj[1]) - (Bpars[l,2]*indices[j,l] == maxj[2])))
        })
    })
    return(t(result))
}

WLSminimizeML<-function(pars,totlist, indices, Bmatrix, w){
    Bpars <- Bmatrix(pars)
    if(any(rowSums(Bpars) > 1) || any(c(Bpars) < 0)){
        return(10^10)
    } else{
        res <- apply(indices, 1, function(j) sum(apply(Bpars, 2, function(l) max(l*j))))
        diff<- res - totlist
        return(t(diff) %*% w %*% diff)
    }
}
WLSminimizeMLcu<-function(pars,totlist, indices, Bmatrix){
    Bpars <- Bmatrix(pars)
    if(any(rowSums(Bpars) > 1) || any(c(Bpars) < 0)){
        return(10^10)
    } else{
        res <- apply(indices, 1, function(j) sum(apply(Bpars, 2, function(l) max(l*j))))
        diff<- res - totlist
        Sigma <- AsymVarMLwls(Bpars,indices,d = ncol(indices))
        while(rcond(Sigma) < 1e-05){Sigma <- Sigma + (0.001)*diag(length(totlist))}
        return(t(diff) %*% solve(Sigma) %*% diff)
    }
}

################# Logistic: M-estimator
logfunc<-function(x,theta){return((sum(x^(1/theta)))^theta)}

MestimatorMinimizeGumbel<-function(pars,q,totlist){
    res<-adaptIntegrate(logfunc,lowerLimit=c(0,0),upperLimit=c(1,1),theta=pars)$integral
    diff<-rep(res,q) - totlist
    return(t(diff) %*% diff)
}

ellLog <- function(x,theta){return((sum(x^(1/theta)))^theta)}
ellLogv2 <- function(x2,x1,theta){return((x1^(1/theta) + x2^(1/theta))^theta)}
ellLogv3 <- function(x2,x1,theta){return((x1^(1/theta - 1))*((x1^(1/theta) + x2^(1/theta))^(theta-1)))}
ellLogCase2<-function(x,theta){return(x[1]*ellLog(x,theta=theta))} #dim(x)=3
ellLogCase3<-function(x,theta){return(x[1]*x[2]*ellLog(x,theta=theta))} #dim(x)=2

Ilog<-function(theta){return(adaptIntegrate(ellLog,lowerLimit=c(0,0),upperLimit=c(1,1),theta=theta)$integral)}
I2log<-function(x,theta){stats::integrate(ellLogv2,lower=0,upper=1,x1=x,theta=theta,abs.tol=0)$value}
I3log<-function(x,theta){stats::integrate(ellLogv3,lower=0,upper=1,x1=x,theta=theta,abs.tol=0)$value}

I1finlog<-function(u,theta){return(I3log(u[3],theta)*ellLog(u,theta))} #dim(u)=3
I1finlogcase2<-function(u,theta){return(I3log(u[3],theta)*ellLog(c(max(u[1],u[3]),u[2]),theta))}
I3finlog<-function(u,theta){return(I3log(u[1],theta)*I3log(u[2],theta)*ellLog(u,theta))} #dim(u)=2
I3finlogcase2<-function(u,theta){return(I3log(u[1],theta)*I3log(u[2],theta)*max(u[1],u[2]))} #dim(u)=2

intLog1<-function(theta){ # all four different
    C1<-Ilog(theta)
    C2<-I2log(1,theta)
    T1<-(2*C1 - adaptIntegrate(ellLog,lowerLimit=rep(0,4),upperLimit=rep(1,4),theta=theta)$integral)
    T2<-(4*C1*(C2 -0.5) + 4*C2 - 4*C1 -
             4*adaptIntegrate(I1finlog,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),theta=theta)$integral)
    T4<-(4*(C1 - C2)*(1 - 2*C2)  -
             4*adaptIntegrate(I3finlog,lowerLimit=c(0,0),upperLimit=c(1,1),theta=theta)$integral)
    return(T1 - T2 + T4)
}
intLog2<-function(theta){ # s = u or t = v or t = u
    C1<-Ilog(theta)
    C2<-I2log(1,theta)
    T1<-(2*C1 - 2*adaptIntegrate(ellLogCase2,lowerLimit=rep(0,3),upperLimit=rep(1,3),theta=theta)$integral)
    T2<-(4*C1*(C2 -0.5) + 4*C2 - 4*C1 -
             2*adaptIntegrate(I1finlog,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),theta=theta)$integral -
             2*adaptIntegrate(I1finlogcase2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),theta=theta)$integral)
    T4<-(4*(C1 - C2)*(1 - 2*C2)  -
             3*adaptIntegrate(I3finlog,lowerLimit=c(0,0),upperLimit=c(1,1),theta=theta)$integral -
             adaptIntegrate(I3finlogcase2,lowerLimit=c(0,0),upperLimit=c(1,1),theta=theta)$integral)
    return(T1 -T2 + T4)
}

intLog3<-function(theta){ # s=u,t=v
    C1<-Ilog(theta)
    C2<-I2log(1,theta)
    T1<-(2*C1 - 4*adaptIntegrate(ellLogCase3,lowerLimit=rep(0,2),upperLimit=rep(1,2),theta=theta)$integral)
    T2<-(4*C1*(C2 -0.5) + 4*C2 - 4*C1 -
             4*adaptIntegrate(I1finlogcase2,lowerLimit=c(0,0,0),upperLimit=c(1,1,1),theta=theta)$integral)
    T4<-(4*(C1 - C2)*(1 - 2*C2)  -
             2*adaptIntegrate(I3finlog,lowerLimit=c(0,0),upperLimit=c(1,1),theta=theta)$integral -
             2*adaptIntegrate(I3finlogcase2,lowerLimit=c(0,0),upperLimit=c(1,1),theta=theta)$integral)
    return(T1 -T2 + T4)
}

psiLogToInt<-function(x, theta){
    y <- x^(1/theta)
    return((sum(y)^(theta-1))*(log(sum(y))*sum(y) - (1/theta)*(y[1]*log(x[1]) + y[2]*log(x[2]))))
}

psiMestimatorGumbel<-function(theta,npairs){
    val <- adaptIntegrate(psiLogToInt, lower = c(0,0), upper = c(1,1), theta = theta)$integral
    return(rep(val,npairs))
}

AsymVarMestimatorGumbel<-function(indices,theta){
    nind<-nrow(indices)
    logMat<-matrix(,nrow=nind,ncol=nind)
    case1<-intLog1(theta=theta)
    case2<-intLog2(theta=theta)
    case3<-intLog3(theta=theta)
    for(i in 1:nind){
        for(j in i:nind){
            temp<-c(which(indices[i,]==1),which(indices[j,]==1))
            if(length(unique(temp))==2){
                logMat[i,j]<-case3
            }else if(length(unique(temp))==3){
                logMat[i,j]<-case2
            }else if(length(unique(temp))==4){
                logMat[i,j]<-case1
            }
            logMat[j,i]<-logMat[i,j]
        }
    }
    return(logMat)
}

###################### Logistic: WLS estimator ######################################
WLSminimizeGumbel<-function(theta,indices,totlist){
    res <- apply(indices, 1, function(i) sum(i^(1/theta))^theta)
    diff<- rep(res,nrow(indices))-totlist
    return(t(diff) %*% diff)
}

fecGum<-function(theta, cst){return(sum(cst^(1/theta))^theta)}
ecderGum <- function(theta,cst){return((sum(cst^1/theta)^(theta - 1))*(cst^(1/theta - 1)))}

AsymWLSGum <- function(C1,C2,theta,d){
    ec1 <- fecGum(theta, C1)
    ec2 <- fecGum(theta, C2)
    der1 <- ecderGum(theta,C1)
    der2 <- ecderGum(theta,C2)
    T1 <- ec1 + ec2 - fecGum(theta, pmax(C1,C2))
    T2 <- sum(sapply(c(1:d), function(k) {
        Ctemp <- rep(0,d)
        Ctemp[k] <- C2[k]
        der2[k]*(ec1 + C2[k] - fecGum(theta,pmax(C1,Ctemp)))
    }))
    T3 <- sum(sapply(c(1:d), function(k) {
        Ctemp <- rep(0,d)
        Ctemp[k] <- C1[k]
        der1[k]*(ec2 + C1[k] - fecGum(theta,pmax(Ctemp,C2)))
    }))
    T4 <- sum(sapply(c(1:d), function(k) sapply(c(1:d), function(l){
        Ctemp1 <- Ctemp2 <- rep(0,d)
        Ctemp1[k] <- C1[k]
        Ctemp2[l] <- C2[l]
        der1[k]*der2[l]*(C1[k] + C2[l] - fecGum(theta,pmax(Ctemp1,Ctemp2)))
    })))
    return(T1 - T2 - T3 + T4)
}

AsymVarWLSGumbel <- function(theta,indices){
    qd <- dim(indices)
    resmat <- matrix(0,nrow=qd[1],ncol=qd[1])
    resmat[lower.tri(resmat, diag = TRUE)] <- unlist(sapply(c(1:qd[1]), function(i) sapply(c(i:qd[1]),
                               function(j) AsymWLSGum(indices[i,],indices[j,],theta,qd[2]))))
    if(all(dim(resmat) == c(1,1))){
        return(resmat)
    } else{
        return(t(resmat) + resmat - diag(diag(resmat)))
    }
}

psiWLSGumbel<-function(theta, indices){
    result <- apply(indices, 1, function(x){
        y <- x^(1/theta)
        yP <- y[y>0]
        return((sum(yP)^(theta - 1))*(log(sum(yP))*sum(yP) - (1/theta)*sum(yP*log(x[x>0]))))
    })
    return(matrix(result,ncol=1))
}

psiMestimatorGumbel<-function(theta,npairs){
    val <- adaptIntegrate(psiLogToInt, lowerLimit = c(0,0), upperLimit = c(1,1), theta = theta)$integral
    return(rep(val,npairs))
}

############### Data application
BmatrixEURO <- function(theta){
    row6<-c(max(theta[1]*theta[5],theta[3]*theta[6]), (1-theta[1])*theta[5], (1-theta[3])*theta[6])
    row7<-c(max(theta[1]*theta[7],theta[3]*theta[8]), (1-theta[1])*theta[7], (1-theta[3])*theta[8])
    row8<-c(max(theta[2]*theta[9],theta[3]*theta[10]), (1-theta[2])*theta[9], (1-theta[3])*theta[10])
    row9<-c(max(theta[2]*theta[11],theta[4]*theta[12]), (1-theta[2])*theta[11], (1-theta[4])*theta[12])
    row10<-c(max(theta[1]*theta[13],theta[4]*theta[14]), (1-theta[1])*theta[13], (1-theta[4])*theta[14])
    return(rbind(c(1,0,0,0,0,0,0,0,0,0),
                 c(theta[1],1-theta[1],0,0,0,0,0,0,0,0),
                 c(theta[2],0,1-theta[2],0,0,0,0,0,0,0),
                 c(theta[3],0,0,1-theta[3],0,0,0,0,0,0),
                 c(theta[4],0,0,0,1-theta[4],0,0,0,0,0),
                 c(row6[1],row6[2],0,row6[3],0,1-sum(row6),0,0,0,0),
                 c(row7[1],row7[2],0,row7[3],0,0,1-sum(row7),0,0,0),
                 c(row8[1],0,row8[2],row8[3],0,0,0,1-sum(row8),0,0),
                 c(row9[1],0,row9[2],0,row9[3],0,0,0,1-sum(row9),0),
                 c(row10[1],row10[2],0,0,row10[3],0,0,0,0,1-sum(row10))))
}

LdotEURO <- function(pars, grid){
    q <- nrow(grid)
    u12 <- pars[1] ; u13 <- pars[2]; u14 <- pars[3]; u15 <- pars[4]
    u26 <- pars[5]; u46 <- pars[6]; u27 <- pars[7]; u47 <- pars[8]
    u38 <- pars[9]; u48 <- pars[10]; u39 <-pars[11]; u59 <- pars[12]; u210 <- pars[13]; u510 <- pars[14]
    result <- sapply(c(1:q), function(i){
        deriv <- rep(0,14)
        M <- max(grid[i,1], u12*grid[i,2], u13*grid[i,3], u14*grid[i,4], u15*grid[i,5], u12*u26*grid[i,6],
                 u14*u46*grid[i,6], u12*u27*grid[i,7], u14*u47*grid[i,7], u13*u38*grid[i,8], u14*u48*grid[i,8],
                 u13*u39*grid[i,9], u15*u59*grid[i,9], u12*u210*grid[i,10], u15*u510*grid[i,10])
        deriv[1] <- grid[i,2]*(M==u12*grid[i,2]) + u26*grid[i,6]*(M==u12*u26*grid[i,6]) + u27*grid[i,7]*(M==u12*u27*grid[i,7]) +
            grid[i,10]*u210*(M==u12*u210*grid[i,10]) - max(grid[i,2], u26*grid[i,6], u27*grid[i,7],u210*grid[i,10]) +
            u26*grid[i,6]*(u12*u26 < u14*u46) + grid[i,7]*u27*(u12*u27 < u14*u47) + grid[i,10]*u210*(u12*u210 < u15*u510)
        deriv[2] <- grid[i,3]*(M==u13*grid[i,3]) + u38*grid[i,8]*(M==u13*u38*grid[i,8]) + u39*grid[i,9]*(M==u13*u39*grid[i,9]) -
            max(grid[i,3],u38*grid[i,8],u39*grid[i,9]) + u38*grid[i,8]*(u13*u38 < u14*u48) +
            u39*grid[i,9]*(u13*u39 < u15*u59)
        deriv[3] <- grid[i,4]*(M==u14*grid[i,4]) + u46*grid[i,6]*(M==u14*u46*grid[i,6]) + u47*grid[i,7]*(M==u14*u47*grid[i,7]) +
            u48*grid[i,8]*(M==u14*u48*grid[i,8]) - max(grid[i,4],u46*grid[i,6],u47*grid[i,7],u48*grid[i,8]) +
            u46*grid[i,6]*(u14*u46 < u12*u26) + u47*grid[i,7]*(u14*u47 < u12*u27) + u48*grid[i,8]*(u14*u48 < u13*u38)
        deriv[4] <- grid[i,5]*(M==u15*grid[i,5]) + u59*grid[i,9]*(M==u15*u59*grid[i,9]) + u510*grid[i,10]*(M==u15*u510*grid[i,10]) -
            max(grid[i,5],u59*grid[i,9],u510*grid[i,10]) + u59*grid[i,9]*(u15*u59 < u13*u39) +
            u510*grid[i,10]*(u15*u510 < u12*u210)
        deriv[5] <- u12*grid[i,6]*(M==u12*u26*grid[i,6]) - grid[i,6] + u12*grid[i,6]*(u12*u26 < u14*u46) +
            (1-u12)*grid[i,6]*(u26*grid[i,6] > max(grid[i,2],u27*grid[i,7],u210*grid[i,10]))
        deriv[6] <- u14*grid[i,6]*(M==u14*u46*grid[i,6]) - grid[i,6] + u14*grid[i,6]*(u14*u46 < u12*u26) +
            (1-u14)*grid[i,6]*(u46*grid[i,6] > max(grid[i,4],u47*grid[i,7],u48*grid[i,8]))
        deriv[7] <- u12*grid[i,7]*(M==u12*u27*grid[i,7]) - grid[i,7] + u12*grid[i,7]*(u12*u27 < u14*u47) +
            (1-u12)*grid[i,7]*(u27*grid[i,7] > max(grid[i,2],u26*grid[i,6],u210*grid[i,10]))
        deriv[8] <- u14*grid[i,7]*(M==u14*u47*grid[i,7]) - grid[i,7] + u14*grid[i,7]*(u14*u47 < u12*u27) +
            (1-u14)*grid[i,7]*(u47*grid[i,7] > max(grid[i,4],u46*grid[i,6],u48*grid[i,8]))
        deriv[9] <- u13*grid[i,8]*(M==u13*u38*grid[i,8]) - grid[i,8] + u13*grid[i,8]*(u13*u38 < u14*u48) +
            (1-u13)*grid[i,8]*(u38*grid[i,8] > max(grid[i,3],u39*grid[i,9]))
        deriv[10]<- u14*grid[i,8]*(M==u14*u48*grid[i,8]) - grid[i,8] + u14*grid[i,8]*(u14*u48 < u13*u38) +
            (1-u14)*grid[i,8]*(u48*grid[i,8] > max(grid[i,4], u46*grid[i,6],u47*grid[i,7]))
        deriv[11]<- u13*grid[i,9]*(M==u13*u39*grid[i,9]) - grid[i,9] + u13*grid[i,9]*(u13*u39 < u15*u59) +
            (1-u13)*grid[i,9]*(u39*grid[i,9] > max(grid[i,3], u38*grid[i,8]))
        deriv[12]<- u15*grid[i,9]*(M==u15*u59*grid[i,9]) - grid[i,9] + u15*grid[i,9]*(u15*u59 < u13*u39) +
            (1-u15)*grid[i,9]*(u59*grid[i,9] > max(grid[i,5], u510*grid[i,10]))
        deriv[13]<- u12*grid[i,10]*(M==u12*u210*grid[i,10]) - grid[i,10] + u12*grid[i,10]*(u12*u210 < u15*u510) +
            (1-u12)*grid[i,10]*(u210*grid[i,10] > max(grid[i,2], u26*grid[i,6], u27*grid[i,7]))
        deriv[14]<- u15*grid[i,10]*(M==u15*u510*grid[i,10]) - grid[i,10] + u15*grid[i,10]*(u15*u510 < u12*u210) +
            (1-u15)*grid[i,10]*(u510*grid[i,10] > max(grid[i,5],u59*grid[i,9]))
        return(deriv)
    })
    return(t(result))
}

