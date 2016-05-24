#-- required package --#
# library(survival) ; 

#######################################################
# Function: Est.Cval: return point estimte
#           ver.002  --- 2011.2.15 --  can handle ties 
#           ver.003  --- 2011.9.21 --  FORTRAN (conc)
#           ver.003b --- 2013.2.13 --  add nofit option
#######################################################
Est.Cval<-function(mydata, tau, nofit=FALSE){

    ## =============== ##
    ## Weight          ##
    ## =============== ##
    cens <- kmcens(mydata[,1],mydata[,2], tau)
    GXi<-cens$surv[match(mydata[,1], cens$distinct, nomatch=1)]
    Wi<-1/GXi/GXi*mydata[,2]*as.numeric(mydata[,1]<tau) 

    ## ======================== ##
    ## Fit Cox if nofit==FALSE  ##
    ## ======================== ##

    if(nofit){
       rs = mydata[,3]
    }else{
       fit.cox <- Est.PH(mydata)
       rs=as.vector(fit.cox$rs)
    }
 
    ## =============== ##
    ## C-stat (D hat)  ##
    ## =============== ##
    cstat=conc(mydata[,1],mydata[,2],Wi,rs)
    
    ## =============== ##
    ## OUTPUT          ##
    ## =============== ##
    Z=list()
    Z$Dhat=cstat
    Z$rs=rs
    if(nofit==FALSE){
     Z$beta=fit.cox$beta
     Z$beta.var=fit.cox$var
     Z$rs=fit.cox$rs
     Z$Ui=fit.cox$Ui 
     Z$ft=fit.cox$ft
    }
    Z$cens.surv=cens$surv
    Z$cens.psii=cens$psii 
    Z$distinct=cens$distinct
    Z$wt=Wi
    Z$nofit=nofit

    return(Z)

}



#######################################################
# Function: Inf.Cval: give CI
#     ver.003 --- 2011.9.21 --  FORTRAN (conc, unoU2P)
#     ver.003b -- 2013.2.13 --  seed = null (default)
#######################################################
Inf.Cval<-function(mydata, tau, itr=1000, seed=NULL){

	if(!is.null(seed)) {set.seed(seed)}
	
	emp<-Est.Cval(mydata, tau)
	n<-nrow(mydata); p<-ncol(mydata)-2 ; 

	#-------------
	A<-as.matrix(emp$beta.var)
    distinct<-emp$distinct
    cens.surv<-emp$cens.surv

	#-------------
	time<-mydata[,1] ; status<-mydata[,2] ; covs<-mydata[,-1:-2] ;
    UJ<-emp$Ui
    cens.psii<-emp$cens.psii
    Wi<-emp$wt

	#-------------
	temp<-rep(0,itr)
#    Wa.stock<-Wb.stock<-WG.stock<-c()
#    b.star.stock<-c()
#    G.star.stock<-c()

	for (k in 1: itr){

 	    #----- random number ---#
        xi<-rexp(n)

        b.star = emp$beta + n*A%*%unoU2P(emp$Ui,xi) 
        rs.star<-as.matrix(covs)%*%as.vector(b.star)

        G.wk<-rep(0,length(distinct))
        G.wk=unoU2P(cens.psii,xi)
        G.star <- cens.surv - G.wk*cens.surv
        GXi<-G.star[match(mydata[,1], distinct, nomatch=1)]
        Wi.star<-1/GXi/GXi*mydata[,2]*as.numeric(mydata[,1]<tau)

        CW=unoCW(time, status, Wi, Wi.star, emp$rs, rs.star, xi, emp$Dhat)

        #----- result  ---#
	    temp[k] <-  (CW$Wa + CW$Wb + CW$Wg)*sqrt(n)

        #----- for check ----#
        #     b.star.stock<-rbind(b.star.stock, b.star)
        #     G.star.stock<-rbind(G.star.stock, G.star)
        #     Wa.stock<-c(Wa.stock, CW$Wa)
        #     Wb.stock<-c(Wb.stock, CW$Wb)
        #     WG.stock<-c(WG.stock, CW$Wg)
	}
	se=sqrt(var(temp)/n)
	low95<-emp$Dhat - 1.959*se
	upp95<-emp$Dhat + 1.959*se
	# return(list(Dhat=emp$Dhat, low95=low95, upp95=upp95, out=temp, b.star.stock=b.star.stock, G.star.stock=G.star.stock, Wa=Wa.stock, Wb=Wb.stock, WG=WG.stock, se=se))
	return(list(Dhat=emp$Dhat, low95=low95, upp95=upp95, se=se, ft=emp$ft))
}


#######################################################
# Function: Inf.Cval.Delta: give Delta C and its CI
#     ver.003b -- 2013.2.13 --  seed = null (default)
#######################################################
Inf.Cval.Delta<-function(mydata, covs0, covs1, tau, itr=1000, seed=NULL){
   
   mydata1=data.frame(mydata, covs1)
   mydata0=data.frame(mydata, covs0)


	if(!is.null(seed)) {set.seed(seed)}

	n<-nrow(mydata1); p<-ncol(mydata1)-2 ; 

	emp1<-Est.Cval(mydata1, tau)
	emp0<-Est.Cval(mydata0, tau)
	
	#-------------
    distinct<-emp1$distinct
    cens.surv<-emp1$cens.surv

	A1<-as.matrix(emp1$beta.var)
	A0<-as.matrix(emp0$beta.var)

	#-------------
	time<-mydata1[,1] ; status<-mydata1[,2] ; 
	covs1<-mydata1[,-1:-2] ;
	covs0<-mydata0[,-1:-2] ;

    cens.psii<-emp1$cens.psii
    Wi<-emp1$wt

    UJ1<-emp1$Ui
    UJ0<-emp0$Ui

	#-------------
	temp=temp1=temp0=rep(0,itr)

	for (k in 1: itr){


        xi<-rexp(n)
        
        b1.star = emp1$beta + n*A1%*%unoU2P(emp1$Ui,xi) 
        rs1.star<-as.matrix(covs1)%*%as.vector(b1.star)

        b0.star = emp0$beta + n*A0%*%unoU2P(emp0$Ui,xi) 
        rs0.star<-as.matrix(covs0)%*%as.vector(b0.star)

        G.wk<-rep(0,length(distinct))
        G.wk=unoU2P(cens.psii,xi)
        G.star <- cens.surv - G.wk*cens.surv
        GXi<-G.star[match(mydata1[,1], distinct, nomatch=1)]
        Wi.star<-1/GXi/GXi*mydata1[,2]*as.numeric(mydata1[,1]<tau)

        CW1=unoCW(time, status, Wi, Wi.star, emp1$rs, rs1.star, xi, emp1$Dhat)
        CW0=unoCW(time, status, Wi, Wi.star, emp0$rs, rs0.star, xi, emp0$Dhat)

        #----- result  ---#
	    temp1[k] <-  (CW1$Wa + CW1$Wb + CW1$Wg)*sqrt(n)
	    temp0[k] <-  (CW0$Wa + CW0$Wb + CW0$Wg)*sqrt(n)
	    temp[k] <-  temp1[k]-temp0[k]

 	}
	se=sqrt(var(temp)/n)
	Delta.hat<-emp1$Dhat - emp0$Dhat
	low95<-Delta.hat - 1.959*se
	upp95<-Delta.hat + 1.959*se

	se1=sqrt(var(temp1)/n)
	c1low95<-emp1$Dhat - 1.959*se1
	c1upp95<-emp1$Dhat + 1.959*se1

	se0=sqrt(var(temp0)/n)
	c0low95<-emp0$Dhat - 1.959*se0
	c0upp95<-emp0$Dhat + 1.959*se0

    C0=c(emp0$Dhat,se0, c0low95, c0upp95)
    C1=c(emp1$Dhat,se1, c1low95, c1upp95)
    DeltaC=c(Delta.hat, se, low95, upp95)
    out=as.matrix(rbind(C1,C0,DeltaC))
    rownames(out)=c("Model1","Model0","Delta")
    colnames(out)=c("Est","SE","Lower95","Upper95")
    return(out)   
#	return(list(C0=C0, C1=C1, DeltaC=DeltaC))
}



#######################################################
## Function: Est.PH: Fit Cox and return beta'x
## --fixed in version 2
#######################################################
Est.PH<-function(mydata){

	covs<-as.matrix(mydata[,-1:-2])
	ft<-coxph(Surv(mydata[,1],mydata[,2]) ~ covs)
	rs<-as.matrix(covs)%*%as.vector(ft$coefficient)
	
	#---- score ----#
	time <- mydata[,1]
	n<-nrow(covs); p<-ncol(covs)
	S0<-rep(0, n) ; S1<-matrix(0, n, p)
	for (i in 1:n){
		S0[i]=as.numeric(time>=time[i])%*%exp(rs)/n ;
		S1[i,]=apply(VEC2MAT(as.numeric(time>=time[i]),p)*VEC2MAT(exp(rs),p)*t(as.matrix(covs)), 1, mean) ;
		}
 	 Ui<-(covs - S1/S0)*mydata[,2]
	
	return(list(beta=ft$coefficient, var=ft$var, rs=rs, Ui=Ui, ft=ft))
	}


#============================================================
VEC2MAT<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=TRUE)
    }


#######################################################
## Function: kmcens: gives Keplan-Meier of Censoring ## 
#######################################################
## I: time      (nx1)
##    status    (nx1)  1=fail, 0=censor  
##    tau       (1x1)  truncate point P(C>tau)>0
##=====================================================
## O: surv	(tx1) G(t) <---------Kaplan-Meier
##    nelson	(tx1) Lambda_C(t) <--Nelson-Alan 
##    distinct	(tx1) t <------------Unique time point at observed event
##    pi_0(t)   (tx1) pr[1{X>=t}] 
##    pi_X(t)   (tx1) pi_X = pr[ 1{X>=t} * G(t) ]
##    pi_T(t)   (tx1) pi_T = pr[ 1{X>=t} / G(t) ]
##    Fn(t)     (tx1)   Fn = pr[ 1{X<=t} & stuatus==1]
##    Mi(t)     (nxt) Martingale for censoring
##    psii(t)	(nxt) iid term
#######################################################
kmcens<-function(time, status, tau){

#-- initial --#

 distinct<-unique(sort(time))
 #--- do not do this ---#
 # distinct<-unique(sort(time[status==1 & time < tau]))

 t<-length(distinct)
 n<-length(time)
 surv<-rep(0,t)
 nel.wk<-rep(0,t)
 nelson<-rep(0,t)

#- i = 1 -#
   yi<-sum(as.numeric(time>=distinct[1]))
   di<-sum(as.numeric(time==distinct[1] & status==0))
   surv[1]<- 1 * (1-di/yi)
   nel.wk[1]<-di/yi

#- i = 2 to t -#
 for (i in 2:t){
   yi<-sum(as.numeric(time>=distinct[i]))
   di<-sum(as.numeric(time==distinct[i] & status==0))
   surv[i]<-surv[i-1]*(1-di/yi)
   nel.wk[i]<-di/yi
   }

 #- shift -#
  surv[2:t]<-surv[1:(t-1)]
  surv[1]<-1
  nel.wk[2:t]<-nel.wk[1:(t-1)]
  nel.wk[1]<-0
 #---------#

  nelson<-cumsum(nel.wk)

 #--------------------------------------------------------#
 # pi_0(t) --- (tx1 vector)--- pr[ 1{X>=t}]              -#
 # pi_X(t) --- (tx1 vector)--- pr[ 1{X>=t} * G(t) ]      -#
 # pi_T(t) --- (tx1 vector)--- pr[ 1{X>=t} / G(t) ]      -#
 #   Fn(t) --- (tx1 vector)--- pr[ 1{X<=t} & stuatus==1] -#
 #--------------------------------------------------------#
 pi_0<-rep(0,t)
 pi_X<-rep(0,t)
 pi_T<-rep(0,t)
   Fn<-rep(0,t)
 for (i in 1:t){
  yi<- as.numeric(time >= distinct[i])
  ni<- as.numeric(time <= distinct[i] & status==1)
  pi_0[i]<-mean(yi)
  pi_X[i]<-mean(yi)*surv[i]
  pi_T[i]<-mean(yi)/surv[i]
    Fn[i]<-mean(ni)
 }

 #-------------------------------#
 #  Mi(t) --- (n x t) Martingale #
 #-------------------------------#
 Mi<-matrix(0, n, t)
 wk1<-matrix(0, n, t)
 wk2<-matrix(0, n, t)

 for (j in 1:t){ 
   wk1[,j]<- as.numeric(time <= distinct[j] & status==0)
   wk2[,j]<- as.numeric(time >= distinct[j])
 }
 for (k in 1:n){
   Mi[k,] <- wk1[k,] - wk2[k,]*nelson
 }
 

 #--------------------------------#
 # psii(t)	   (n x t matrix) #
 #--------------------------------#
 psii<-matrix(0, n, t) 
  wk1<-matrix(0, n, t)
  wk2<-matrix(0, n, t)

 for (j in 1:t){ 
   wk1[,j]<- as.numeric(time == distinct[j] & status==0)
   wk2[,j]<- as.numeric(time >= distinct[j])
 }

for (i in 1:n){ 
   psii[i,]<- (cumsum(wk1[i,]/pi_0) - cumsum(nel.wk*wk2[i,]/pi_0))
 }


 #----------------------------------------------------------#
 wk1<-NULL ; wk2<-NULL
 return(list(surv=surv, nelson=nelson, distinct=distinct, pi_0=pi_0, pi_X=pi_X, pi_T=pi_T, Mi=Mi, psii=psii, Fn=Fn))
}


#######################################################
# Function: cvC: return cross-validation estimates 
#######################################################
cvC=function(mydata, tau, cvK=10, Rep=10){

time=mydata[,1]
fail=mydata[,2]
X=as.matrix(mydata[,c(-1,-2)])

Dhat=rep(0,cvK)
Cest=rep(0, Rep)	

for(r in 1:Rep){
    group=sample(1:cvK,nrow(mydata),replace=TRUE)

    for(kk in 1:cvK){
     mod0=coxph(Surv(time,fail)~X,subset=group!=kk)
     XTest=as.matrix(X[group==kk,])
     scores=XTest%*%mod0$coef
     D=data.frame(cbind(time[group==kk],fail[group==kk],scores))
     Dhat[kk]=Est.Cval(D,tau)$Dhat
    }
    Cest[r]=weighted.mean(Dhat,table(group))
 }
  return(mean(Cest))
}


