#' Predictions of life expectancies from an illness-death.
#' 
#' Predict life expectancies from an object of class \code{idm}. Life
#' expectancies are calculated at time \code{s} for a subject who has the
#' covariates values Z01, Z02, Z12. Confidence intervals are calculated.
#' 
#' 
#' @param x an object as returned by a call to the \code{\link{idm}}
#' function.
#' @param s time at prediction.
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
#' The default is 1000.
#' @param CI boolean: with (\code{TRUE}) or without (\code{WRONG}) confidence
#' intervals for the life expectancies. The default is \code{TRUE}.
#' @param \dots others arguments.
#' @return a list containing: \item{life.in.0.expectancy}{life expectancy in
#' state 0 and confidence interval.} \item{life.expectancy.nondis}{life
#' expectancy of a non-diseased subject and confidence interval.}
#' \item{life.expectancy.dis}{life expectancy of a diseased subject and
#' confidence interval.}
#' @author C. Touraine
#' @seealso \code{\link{idm}}
#' @keywords methods
##' @examples
##' \dontrun{
##' library(lava)
##' library(prodlim)
##' set.seed(17)
##' d <- simulateIDM(100)
##' table(d$seen.ill,d$seen.exit)
##' fitIC <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'              conf.int=FALSE)
##' try(lifexpect(fitIC,s=10),silent=TRUE)
##' 
##' 
##'     data(Paq1000)
##' 
##'     fit <- idm(formula02=prodlim::Hist(time=t,event=death,entry=e)~certif,
##'                formula01=prodlim::Hist(time=list(l,r),event=dementia)~certif,
##'                method="Splines",
##'                data=Paq1000,
##'                conf.int=FALSE) 
##' 
##'     pred <- lifexpect(fit,s=70,t=80,Z01=c(1),Z02=c(1),Z12=c(1))
##' 
##' }
 
#' @export lifexpect
lifexpect <- function(x,s,Z01,Z02,Z12,nsim=1000,CI=TRUE,...) {
    if (x$method=="Weib"){
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
            if (length(Z01) != nvar01) {
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
            if (length(Z02) != nvar02) {
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
            if (length(Z12) != nvar12) {
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
        res <- lifexpect0.idmWeib(s,a01,1/b01,a02,1/b02,a12,1/b12,bZ01,bZ02,bZ12)
        if (CI==TRUE) {
            ### CI prediction by Monte-Carlo
            Vmean <- c(sqrt(a01),sqrt(b01),sqrt(a02),sqrt(b02),sqrt(a12),sqrt(b12),beta01,beta02,beta12) # vector of parameters
            Mcov = x$V
            # une simulation pour chaque element d'une liste
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
            X <- mvtnorm::rmvnorm(nsim,Vmean,Mcov) 
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
            Xres1 <- mapply(function(Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12) lifexpect0.idmWeib(s,a01=Xa01,b01=Xb01,a02=Xa02,b02=Xb02,a12=Xa12,b12=Xb12,bZ01=XbZ01,bZ02=XbZ02,bZ12=XbZ12),Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12)
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
            return(list(life.in.0.expectancy=c(res$life.in.0.expectancy,Xres4[1,]),
                        life.expectancy.nondis=c(res$life.expectancy.nondis,Xres4[2,]),
                        life.expectancy.dis=c(res$life.expectancy.dis,Xres4[3,])))	
        }
        else 
            return(res)
    } else {
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
            if (length(Z01) != nvar01) {
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
            if (length(Z02) != nvar02) {
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
            if (length(Z12) != nvar12) {
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
        res <- lifexpect0(s,zi01,nz01,the01^2,zi12,nz12,the12^2,zi02,nz02,the02^2,bZ01,bZ12,bZ02)
        if (CI==TRUE) {
            ### CI prediction by Monte-Carlo
            Vmean <- c(the01,the02,the12,beta01,beta02,beta12) # vector of parameters
            Mcov = x$V
            # une simulation pour chaque element d'une liste
            Xtheta01 <- as.list(NULL)
            Xtheta02 <- as.list(NULL)
            Xtheta12 <- as.list(NULL)
            XbZ01 <- as.list(NULL)
            XbZ02 <- as.list(NULL)
            XbZ12 <- as.list(NULL)
            set.seed(21)
            X <- mvtnorm::rmvnorm(nsim,Vmean,Mcov) 
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
            Xres1 <- mapply(function(Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02) lifexpect0(s,zi01,nz01,Xtheta01,zi12,nz12,Xtheta12,zi02,nz02,Xtheta02,XbZ01,XbZ12,XbZ02),Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02)
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
            return(list(life.in.0.expectancy=c(res$life.in.0.expectancy,Xres4[1,]),
                        life.expectancy.nondis=c(res$life.expectancy.nondem,Xres4[2,]),
                        life.expectancy.dis=c(res$life.expectancy.dem,Xres4[3,])))	
        }
        else 
            return(res)
  		

            }
}

lifexpect0.idmWeib <- function(s,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {
    ## print("lifexpect0.idmWeib")
    ET12 = integrate(
        f=function(x) {
            S.weib(s,x,a12,b12,bZ12)
        },s,Inf)
    ET0. = integrate(f=function(x) {
        S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)
    },s,Inf)
    ET01 = integrate(f=function(x) {
        sapply(x,function(x) {integrate(f=function(y)
                                        {
                                            S.weib(s,y,a01,b01,bZ01)*S.weib(s,y,a02,b02,bZ02)*iweibull(y,a01,b01,bZ01)*S.weib(y,x,a12,b12,bZ12)
                                        }
                                        ,lower=s,upper=x)$value})
    },s,Inf)
    return(list(life.in.0.expectancy=ET0.$value,life.expectancy.nondis=ET01$value+ET0.$value,life.expectancy.dis=ET12$value))

}

lifexpect0 <- function(s,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01=0,bZ12=0,bZ02=0) {
    ET12 = integrate(f=function(x) {
        Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p11
    },s,zi12[nz12+6])
    ET0. = integrate(f=function(x) {
        Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p00
    },s,zi02[nz02+6])
    ET01 = integrate(f=function(x) {
        Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p01
    },s,zi01[nz01+6])
    return(list(life.in.0.expectancy=ET0.$value,life.expectancy.nondis=ET01$value+ET0.$value,life.expectancy.dis=ET12$value))

}
