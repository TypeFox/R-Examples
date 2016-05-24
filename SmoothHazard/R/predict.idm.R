# 0 : health state
# 1 : illness state
# 2 : death state

##### Fonction qui calcule les predictions avec leurs intervalles de confiance entre les temps s et t

#' Predictions for an illness-death model using either a penalized likelihood
#' approach or a Weibull parametrization.
#' 
#' Predict transition probabilities and cumulative probabilities from an object
#' of class \code{idmSplines} with confidence intervals are calculated.
#' 
#' @param object an \code{idm} class objects returned by a call to the
#' \code{\link{idm}} function
#' @param s time at prediction.
#' @param t time for prediction.
#' @param Z01 vector for the values of the covariates on the transition 0 --> 1
#' (in the same order as the covariates within the call. The default values are
#' all 0.
#' @param Z02 vector for the values of the covariates on the transition 0 --> 2
#' (in the same order as the covariates within the call. The default values are
#' all 0.
#' @param Z12 vector for the values of the covariates on the transition 1 --> 2
#' (in the same order as the covariates within the call. The default values are
#' all 0.
#' @param nsim number of simulations for the confidence intervals calculations.
#' The default is 2000.
#' @param CI boolean: with (\code{TRUE}) or without (\code{WRONG}) confidence
#' intervals for the predicted values. The default is \code{TRUE}.
#' @param \dots other parameters.
#' @return a list containing the following predictions with pointwise
#' confidence intervals: \item{p00}{the transition probability \eqn{p_{00}}.}
#' \item{p01}{the transition probability \eqn{p_{01}}.} \item{p11}{the
#' transition probability \eqn{p_{11}}.} \item{p12}{the transition probability
#' \eqn{p_{12}}.} \item{p02_0}{the probability of direct transition from state
#' 0 to state 2.} \item{p02_1}{the probability of transition from state 0 to
#' state 2 via state 1.} \item{p02}{transition probability \eqn{p_{02}}. Note
#' that \code{p02}=\code{p_02_0}+\code{p02_1}.} \item{F01}{the lifetime risk of
#' disease. \code{F01}=\code{p01}+\code{p02_1}.} \item{F0.}{the probability of
#' exit from state 0. \code{F0.}=\code{p02_0}+\code{p01}+\code{p02_1}.}
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{idm}}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' data(Paq1000)
#' library(prodlim)
#' fit <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000) 
#' 
#' pred <- predict(fit,s=70,t=80,Z01=c(1),Z02=c(1),Z12=c(1))
#' pred
#' 
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' 
#' pred <- predict(fit.splines,s=70,t=80,Z01=c(1),Z02=c(1))
#' pred
#' 
#' }
#'
#' @S3method predict idm
predict.idm <- function(object,s,t,Z01,Z02,Z12,nsim=2000,CI=TRUE,...) {
    x <- object
    if (x$method=="Splines"){
        # if covariates: cov=c(cov1,cov2,cov3,...)
        if (inherits(x,"idmSplines")) {
            nvar01 <- x$NC[1]
            nvar02 <- x$NC[2]
            nvar12 <- x$NC[3]
            if(missing(Z01) && nvar01>0) warning("value(s) used for covariate(s) on transition 01: 0 \n")
            if(missing(Z02) && nvar02>0) warning("value(s) used for covariate(s) on transition 02: 0 \n")
            if(missing(Z12) && nvar12>0) warning("value(s) used for covariate(s) on transition 12: 0 \n")
            nz01 <- x$nknots01
            nz02 <- x$nknots02
            nz12 <- x$nknots12
            zi01 <- x$knots01
            zi02 <- x$knots02
            zi12 <- x$knots12
            the01 <- x$theta01
            the02 <- x$theta02
            the12 <- x$theta12
            if (!missing(Z01))  {
                if (length(Z01) != x$NC[1]) {
                    stop("The length of the Z01 arguments must match the number of covariates on the ij transition.")
                }else{
                    beta01 <- x$coef[1:nvar01]
                    bZ01 <- t(beta01)%*%Z01
                }
            }else{
                if (nvar01 != 0) { beta01 <- x$coef[1:nvar01] }
                else { beta01 <- NULL }
                bZ01 <- 0
            }
            if (!missing(Z02))  {
                if (length(Z02) != x$NC[2]) {
                    stop("The length of the Z02 arguments must match the number of covariates on the ij transition.")
                }else{
                    beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)]
                    bZ02 <- t(beta02)%*%Z02
                }
            }else{
                if (nvar02 != 0) { beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)] }
                else { beta02 <- NULL }
                bZ02 <- 0
            }
            if (!missing(Z12))  {
                if (length(Z12) != x$NC[3]) {
                    stop("The length of the Z12 arguments must match the number of covariates on the ij transition.")
                }else{
                    beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)]
                    bZ12 <- t(beta12)%*%Z12
                }
            }else{
                if (nvar12 != 0) { beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)] }
                else { beta12 <- NULL }
                bZ12 <- 0
            }
            res <- Predict0.idmPl(s,t,zi01,nz01,the01^2,zi12,nz12,the12^2,zi02,nz02,the02^2,bZ01,bZ12,bZ02)
            if (CI!=FALSE) {
		### CI prediction by Monte-Carlo
                Vmean <- c(the01,the02,the12,beta01,beta02,beta12) # vector of estimates
                Mvar = x$V # covariance matrix
                Xtheta01 <- as.list(NULL)
                Xtheta02 <- as.list(NULL)
                Xtheta12 <- as.list(NULL)
                XbZ01 <- as.list(NULL)
                XbZ02 <- as.list(NULL)
                XbZ12 <- as.list(NULL)
                set.seed(21)
                X <- mvtnorm::rmvnorm(nsim,Vmean,Mvar) 
                # 1 set of simulated parameters for each element of the list
                for(i in (1:nsim) ) {
                    Xtheta01[[i]]=X[i,1:(nz01+2)]^2
                    Xtheta02[[i]]=X[i,(nz01+3):(nz01+nz02+4)]^2
                    Xtheta12[[i]]=X[i,(nz01+nz02+5):(nz01+nz02+nz12+6)]^2
                }
                for(i in (1:nsim) ) {
                    if (!missing(Z01)) {			
                        XbZ01[[i]]=X[i,(nz01+nz02+nz12+7):(nz01+nz02+nz12+6+nvar01)] %*% Z01
                    } else {
                        XbZ01[[i]]=0
                    }
                    if (!missing(Z02)) {	
                        XbZ02[[i]]=X[i,(nz01+nz02+nz12+6+nvar01+1):(nz01+nz02+nz12+6+nvar01+nvar02)] %*% Z02
                    } else {
                        XbZ02[[i]]=0
                    }
                    if (!missing(Z12)) {	
                        XbZ12[[i]]=X[i,(nz01+nz02+nz12+6+nvar01+nvar02+1):(nz01+nz02+nz12+6+nvar01+nvar02+nvar12)] %*% Z12
                    } else {
                        XbZ12[[i]]=0
                    }
                }
                Xres1 <- mapply(function(Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02) Predict0.idmPl(s,t,zi01,nz01,Xtheta01,zi12,nz12,Xtheta12,zi02,nz02,Xtheta02,XbZ01,XbZ12,XbZ02),Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02)
                Xres2 <- apply(Xres1,2,as.numeric) # transforme le type des elts de la matrice de 'list' en 'numeric'
                Xres3 <- t(apply(Xres2,1,sort)) # classement des valeurs
                iinf=((nsim/100)*2.5) + 1
                isup=(nsim/100)*97.5
                if (iinf != round(iinf)){
                    delta <- (ceiling(iinf)-iinf < iinf-floor(iinf))
                    iinf <- ceiling(iinf)*delta + floor(iinf)*(1-delta)
                }
                if (isup != round(isup)){
                    delta <- (ceiling(isup)-isup < isup-floor(isup))
                    isup <- ceiling(isup)*delta + floor(isup)*(1-delta)
                }
                Xres4 <- cbind(Xres3[,iinf],Xres3[,isup]) # 1ere colonne = bornes inf pour chaque valeur ; 2eme colonne = borne sup pour chaque valeur
                return(list(p00=c(res$p00,Xres4[1,]),p01=c(res$p01,Xres4[2,]),p11=c(res$p11,Xres4[3,]),
                            p12=c(res$p12,Xres4[4,]),p02_0=c(res$p02_0,Xres4[5,]),p02_1=c(res$p02_1,Xres4[6,]),
                            p02=c(res$p02,Xres4[7,]),F01=c(res$F01,Xres4[8,]),F0.=c(res$F0.,Xres4[9,])))	
            }
            else 
                return(res)	
	}
    }else{
        nvar01 <- x$NC[1]
        nvar02 <- x$NC[2]
        nvar12 <- x$NC[3]
        if(missing(Z01) && nvar01>0) warning("value(s) used for covariate(s) on transition 01: 0 \n")
        if(missing(Z02) && nvar02>0) warning("value(s) used for covariate(s) on transition 02: 0 \n")
        if(missing(Z12) && nvar12>0) warning("value(s) used for covariate(s) on transition 12: 0 \n")
        a01 <- x$modelPar[1]
        b01 <- x$modelPar[2]
        a02 <- x$modelPar[3]
        b02 <- x$modelPar[4]
        a12 <- x$modelPar[5]
        b12 <- x$modelPar[6]
        if (!missing(Z01))  {
            if (length(Z01) != x$NC[1]) {
                stop("The length of the Z01 arguments must match the number of covariates on the ij transition.")
            }else{
                beta01 <- x$coef[1:nvar01]
                bZ01 <- t(beta01)%*%Z01
            }
        }else{
            if (nvar01 != 0) { beta01 <- x$coef[1:nvar01] 
                           }else{ beta01 <- NULL 
                              }
            bZ01 <- 0
        }
        if (!missing(Z02))  {
            if (length(Z02) != x$NC[2]) {
                stop("The length of the Z02 arguments must match the number of covariates on the ij transition.")
            }else{
                beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)]
                bZ02 <- t(beta02)%*%Z02
            }
        }else{
            if (nvar02 != 0){ beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)] 
                          }else{ beta02 <- NULL
                             }
            bZ02 <- 0
        }
        if (!missing(Z12))  {
            if (length(Z12) != x$NC[3]) {
                stop("The length of the Z12 arguments must match the number of covariates on the ij transition.")
            }else{
                beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)]
                bZ12 <- t(beta12)%*%Z12
            }
        }else{
            if (nvar12 != 0){ beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)] 
                          }else{ beta12 <- NULL 
                             }
            bZ12 <- 0
        }
        res <- Predict0.idmWeib(s,t,a01,1/b01,a02,1/b02,a12,1/b12,bZ01,bZ02,bZ12)
        if (CI!=FALSE) {
            ### CI prediction by Monte-Carlo
            Vmean <- c(sqrt(a01),sqrt(b01),sqrt(a02),sqrt(b02),sqrt(a12),sqrt(b12),beta01,beta02,beta12) # vector of estimates
            Mvar = x$V # covariance matrix
            Xa01 <- as.list(NULL)
            Xb01 <- as.list(NULL)
            Xa02 <- as.list(NULL)
            Xb02 <- as.list(NULL)
            Xa12 <- as.list(NULL)
            Xb12 <- as.list(NULL)
            XbZ01 <- as.list(NULL)
            XbZ02 <- as.list(NULL)
            XbZ12 <- as.list(NULL)
            set.seed(21)
            X <- mvtnorm::rmvnorm(nsim,Vmean,Mvar) 
            # 1 set of simulated parameters for each element of the list
            for(i in (1:nsim) ) {
                Xa01[[i]]=X[i,1]^2
                Xb01[[i]]=1/(X[i,2]^2)
                Xa02[[i]]=X[i,3]^2
                Xb02[[i]]=1/(X[i,4]^2)
                Xa12[[i]]=X[i,5]^2
                Xb12[[i]]=1/(X[i,6]^2)
            }
            for(i in (1:nsim) ) {
                if (!missing(Z01)) {			
                    XbZ01[[i]]=X[i,(7:(6+nvar01))] %*% Z01
                } else {
                    XbZ01[[i]]=0
                }
                if (!missing(Z02)) {	
                    XbZ02[[i]]=X[i,((6+nvar01+1):(6+nvar01+nvar02))] %*% Z02
                } else {
                    XbZ02[[i]]=0
                }
                if (!missing(Z12)) {	
                    XbZ12[[i]]=X[i,((6+nvar01+nvar02+1):(6+nvar01+nvar02+nvar12))] %*% Z12
                } else {
                    XbZ12[[i]]=0
                }
            }
            Xres1 <- mapply(function(Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12)
                            Predict0.idmWeib(s,t,a01=Xa01,b01=Xb01,a02=Xa02,b02=Xb02,a12=Xa12,b12=Xb12,bZ01=XbZ01,bZ02=XbZ02,bZ12=XbZ12),Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12)
            Xres2 <- apply(Xres1,2,as.numeric) # transforme le type des elts de la matrice de 'list' en 'numeric'
            Xres3 <- t(apply(Xres2,1,sort)) # classement des valeurs
            iinf=((nsim/100)*2.5) + 1
            isup=(nsim/100)*97.5
            if (iinf != round(iinf)){
                delta <- (ceiling(iinf)-iinf < iinf-floor(iinf))
                iinf <- ceiling(iinf)*delta + floor(iinf)*(1-delta)
            }
            if (isup != round(isup)){
                delta <- (ceiling(isup)-isup < isup-floor(isup))
                isup <- ceiling(isup)*delta + floor(isup)*(1-delta)
            }
            Xres4 <- cbind(Xres3[,iinf],Xres3[,isup]) # 1ere colonne = bornes inf pour chaque valeur ; 2eme colonne = borne sup pour chaque valeur
            return(list(p00=c(res$p00,Xres4[1,]),p01=c(res$p01,Xres4[2,]),p11=c(res$p11,Xres4[3,]),
                        p12=c(res$p12,Xres4[4,]),p02_0=c(res$p02_0,Xres4[5,]),p02_1=c(res$p02_1,Xres4[6,]),
                        p02=c(res$p02,Xres4[7,]),F01=c(res$F01,Xres4[8,]),F0.=c(res$F0.,Xres4[9,])))	
        }
        else 
            return(res)	
    }
}


Predict0.idmPl <- function(s,t,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01=0,bZ12=0,bZ02=0) {
    p11 <- rep(0,length(t))
    p12 <- rep(0,length(t))
    p00 <- rep(0,length(t))
    p02_0 <- rep(0,length(t))
    p01 <- rep(0,length(t))
    p02_1 <- rep(0,length(t))
    p02 <- rep(0,length(t))
    if (all(s>t)) {stop("You must respect the condition 's<t' to calculate p(s,t)")} 
    if (s>(min(zi01[nz01+6],zi02[nz02+6],zi12[nz12+6]))) {stop("argument s is out of bornes")}    
    if (any(t>zi12[nz12+6])) {stop("argument t is out of bornes")}
    p11 <- S.pl(s,t,zi12,nz12,the12,bZ12)
    p12 <- 1-p11
    p00 <- S.pl(s,t,zi01,nz01,the01,bZ01)*S.pl(s,t,zi02,nz02,the02,bZ02)
    p02_0 <- sapply(t,function(t) {integrate(f=function(x)
                                             {S.pl(s,x,zi01,nz01,the01,bZ01)*S.pl(s,x,zi02,nz02,the02,bZ02)*susp(x,zi02,nz02,the02,bZ02)$intensity}
                                             ,lower=s,upper=t)$value})
    p01 <- sapply(t,function(t) {integrate(f=function(x)
                                           {S.pl(s,x,zi01,nz01,the01,bZ01)*S.pl(s,x,zi02,nz02,the02,bZ02)*susp(x,zi01,nz01,the01,bZ01)$intensity*S.pl(x,t,zi12,nz12,the12,bZ12)}
                                           ,lower=s,upper=t)$value})
    p02_1 <- 1-p00-p02_0-p01
    p02 <- p02_0+p02_1
    return(list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=p01+p02_1,F0.=p02_0+p01+p02_1))
}




Predict0.idmWeib <- function(s,t,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {
	p11=rep(0,length(t))
	p12=rep(0,length(t))
	p00=rep(0,length(t))
	p02_0=rep(0,length(t))
	p01=rep(0,length(t))
	p02_1=rep(0,length(t))
	p02=rep(0,length(t))
 	if (all(s>t)) {stop("You must respect the condition 's<t' to calculate p(s,t)")
	}else{
		p11 = S.weib(s,t,a12,b12,bZ12)
		p12 = 1-p11
		p00 = S.weib(s,t,a01,b01,bZ01)*S.weib(s,t,a02,b02,bZ02)
		p02_0 = sapply(t,function(t) {integrate(f=function(x)
		{S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a02,b02,bZ02)}
		,lower=s,upper=t)$value })
		p01 = sapply(t,function(t) {integrate(f=function(x)
		{S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)*S.weib(x,t,a12,b12,bZ12)}
		,lower=s,upper=t)$value})
		p02_1 = 1-p00-p02_0-p01
		p02 = p02_0+p02_1
 	}
return(list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=p01+p02_1,F0.=p02_0+p01+p02_1))
}


### fonction d'intensité de transition
# a = shape parameter
# b = scale parameter
# bz = (vecteur des coeffs de régression)^T * (vecteur des variables) 
iweibull <- function(x,a,b,bZ=0) {
  res = (a/b) * (x/b)**(a-1) * exp(bZ)
  return(res) 
}

#############################################


### fonction de survie entre 2 temps s et t
# S(s,t) = S(t)/S(s)
# if S(s)=0, S(s,t)=0
S.weib <- function(s,t,a,b,bZ=0) {	
	res <- 0
	St <- (1-pweibull(t,shape=a,scale=b))^(exp(bZ))
	Ss <- (1-pweibull(s,shape=a,scale=b))^(exp(bZ))
	if (length(s)==1){
		if (Ss==0){res <- 0}
		else{res <- St/Ss}
	}else{
		idx0 <- which(Ss==0)
		idx <- which(Ss!=0)
		res[idx0] <- 0 
		res[idx] <- St/Ss[idx]
	}
	return(res)
}

susp <- function(x,zi,nz,the,bZ=0) {
    gl=rep(0,length(x))   # risque cumule
    lam=rep(0,length(x))  # risque
    su=rep(0,length(x))   # survie
    TF=rep(0,length(x)) # T si z[i-1]<=x[.]<z[i], F sinon
    som=0
    for (i in 5:(nz+3)) {
        TF = ( (zi[i-1]<=x) & (x<zi[i]) )
        if (sum(TF) != 0) { 
            ind = which(TF) 
            mm3=rep(0,length(ind))
            mm2=rep(0,length(ind))
            mm1=rep(0,length(ind))
            mm=rep(0,length(ind))
            im3=rep(0,length(ind))
            im2=rep(0,length(ind))
            im1=rep(0,length(ind))
            im=rep(0,length(ind))
            j = i-1
            if (j>4) { 
                som = sum(the[1:(j-4)])
            }
            ht = x[ind]-zi[j] #
            htm = x[ind]-zi[j-1] #
            h2t = x[ind]-zi[j+2] #
            ht2 = zi[j+1]-x[ind] #
            ht3 = zi[j+3]-x[ind] #
            hht = x[ind]-zi[j-2] #
            h = zi[j+1]-zi[j]
            hh = zi[j+1]-zi[j-1]
            h2 = zi[j+2]-zi[j]
            h3 = zi[j+3]-zi[j]
            h4 = zi[j+4]-zi[j]
            h3m = zi[j+3]-zi[j-1]
            h2n = zi[j+2]-zi[j-1]
            hn= zi[j+1]-zi[j-2]
            hh3 = zi[j+1]-zi[j-3]
            hh2 = zi[j+2]-zi[j-2]
            mm3[ind] = ((4*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2[ind] = ((4*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4*h2t*htm*ht2)/(hh2*h2n*hh*h))+((4*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1[ind] = (4*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4*htm*ht*h2t)/(h3m*h2*h*h2n))+((4*ht3*ht*ht)/(h3m*h3*h2*h))
            mm[ind] = 4*(ht*ht*ht)/(h4*h3*h2*h)
            im3[ind] = (0.25*(x[ind]-zi[j-3])*mm3[ind])+(0.25*hh2*mm2[ind])+(0.25*h3m*mm1[ind])+(0.25*h4*mm[ind])
            im2[ind] = (0.25*hht*mm2[ind])+(h3m*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
            im1[ind] = (htm*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
            im[ind] = ht*mm[ind]*0.25
            gl[ind] = som +(the[j-3]*im3[ind])+(the[j-2]*im2[ind])+(the[j-1]*im1[ind])+(the[j]*im[ind])
            lam[ind] = (the[j-3]*mm3[ind])+(the[j-2]*mm2[ind])+(the[j-1]*mm1[ind])+(the[j]*mm[ind])
        } # fin if (sum(TF) != 0)
    } # fin for
    TF = (x>=zi[nz+3])
    if (sum(TF) != 0) {
        ind = which(TF)
        som = sum(the[1:(nz+2)])
        gl[ind] = som
        lam[ind] = 4*the[nz+2]/(zi[nz+3]-zi[nz+2])
    }
    TF = (x<zi[4])
    if (sum(TF) != 0) {
        ind = which(TF)
        gl[ind] = 0
        lam[ind] = 0
    }
    e = exp(bZ)
    lam=lam*e
    gl=gl*e
    su = exp(-gl)
    return(list(intensity=lam,cumul.intensity=gl,survival=su))
}

A <- function(s,t,zi,nz,the,bZ=0) {
    res=rep(0,length(t))
    TF = (t>=zi[length(zi)])
    ind = which(TF)
    if (sum(TF)!=0) {res[ind]=susp((zi[nz+6]-10^-5),zi,nz,the,bZ)$cumul.intensity-susp(s,zi,nz,the,bZ)$cumul.intensity}
    TF = (t<zi[length(zi)])
    ind = which(TF)
    if (sum(TF)!=0) {res[ind]=susp(t[ind],zi,nz,the,bZ)$cumul.intensity-susp(s,zi,nz,the,bZ)$cumul.intensity}
    return(res)
}

### fonction de survie entre 2 temps s et t
# S(s,t) = S(t)/S(s)
#        = exp(-A(s,t))
S.pl <- function(s,t,zi,nz,the,bZ=0) {
	if (length(t)>=length(s)){
		res=rep(0,length(t))
		TF = (t>zi[length(zi)])
		ind = which(TF)
		if (sum(TF)!=0) {res[ind]=0}
		TF = (t<=zi[length(zi)])
		ind = which(TF)
		if (sum(TF)!=0) {res[ind]=susp(t[ind],zi,nz,the,bZ)$survival/susp(s,zi,nz,the,bZ)$survival}
	}else{		
		res=rep(0,length(s))
		if (t>zi[length(zi)]) {res=0}
		else {res=susp(t,zi,nz,the,bZ)$survival/susp(s,zi,nz,the,bZ)$survival}
	}
	return(res)
}

