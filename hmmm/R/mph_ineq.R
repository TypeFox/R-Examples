#
# 
#
######################   Begin mphineq.fit ########################################## 
mphineq.fit  <- function(y, Z, ZF=Z, h.fct=0,derht.fct=0,d.fct=0,derdt.fct=0,
                   L.fct=0,derLt.fct=0,X=NULL,formula=NULL,names=NULL,lev=NULL,E=NULL,
maxiter=100,step=1,
                   norm.diff.conv=1e-5,norm.score.conv=1e-5,
                   y.eps=0,chscore.criterion=2,
                   m.initial=y,mup=1)
{

start.time <- proc.time()[3]  
version <- "mphineq.fit, version 1.0.1, 5/10/06"
Zlist<-cocadise(Z,formula=formula,lev=lev,names=names)
if(is.null(X)){X<-0}
inv <- ginv 
#solve       
y<-as.matrix(y) 
 
lenh<-0
if ((missing(h.fct))&(sum(abs(X)) != 0))
# 
# if the user inputs the linear predictor model design matrix,
# but does not input the constraint function h.fct then do... 

{
    if(is.null(E)){U <- create.U(X)}
else{U<-t(E)}
	
    if (sum(abs(U)) == 0) {h.fct <- 0}
    else {
       h.fct <- function(m) {
          t(U)%*%L.fct(m)
           }
   }

} 
# end "if missing h.fct and X != 0"
else 

{
   U <- "Not created within the program."
}
#
	
if ((is.function(derht.fct)==FALSE)&(sum(abs(X)) != 0)&(is.function(derLt.fct)==TRUE))
# 
# if the user inputs the linear predictor model design matrix,
# but does not input the constraint function h.fct then do... 
{
    U <- create.U(X)
    if (sum(abs(U)) == 0) {derht.fct <- 0}
    else {
       derht.fct <- function(m) {
          derLt.fct(m)%*%U
          

       }
   }
} 




 
	 
#########################################################################
 
if ((is.function(h.fct)==TRUE)||(is.function(d.fct)==TRUE)||(class(formula)
=="formula"))  {

#
   lenm <- length(y);
   m <- as.matrix(c(m.initial)) + y.eps 
   m[m==0] <- 0.01
  
   p<-m*c(1/Z%*%t(Z)%*%m)

		  

    
  # Zlist<-cocadise(Z,formula=formula,lev=lev)
   xi <- Zlist$IMAT%*%log(m)
if ((is.function(d.fct)==FALSE)&(is.function(h.fct)==FALSE)) {
p  <- as.matrix(exp(Zlist$DMAT%*%xi))
         p<-(p/sum(p))
        m<-p*c(Z%*%t(Z)%*%y)
}

   if (is.function(h.fct)==TRUE){
   	h <- hobs <- h.fct(m)
   lenh <- length(h) 

if (is.function(derht.fct)==FALSE) {
     H <- num.deriv.fct(h.fct,m) 
   }
   else {
     H <- derht.fct(m)
   }  
   #HtDHinv <-
 HtDHinvobs <- inv(t(H)%*%(H*c(m))) #HtDHinv <- inv(t(H)%*%Dm%*%H)
  # HHtDHinv <- H%*%HtDHinv
    }
 	if (is.function(d.fct)==TRUE)
	{
#
# if d(.) is not the zero function, i.e. there is at least one constraint, then do...
# 
   
   d <- dhobs <- d.fct(m)
 lend <- length(d)

 
 ##}
  
  if (is.function(derdt.fct)==FALSE)
   {
     DH <- num.deriv.fct(d.fct,m) 
   }
   else {
     DH <- derdt.fct(m)
 

   }  
}	
#p  <- as.matrix(exp(Zlist$DMAT%*%xi))
        # p<-(p/sum(p))
        # m<-p*c(Z%*%t(Z)%*%y)
		  
  
		
  	
	
  	 
   #lagrange multiplers not needed and computed by this program
   lam <- matrix("NA",lenh,1)
   Dm <- diag(c(m+1e-08))-((ZF*c(m))%*%t(ZF*c(p)))
   
   if (!is.matrix(Dm)) {return ("unable to reach convergence")}






   
   norm.score <- 999999   
   
   theta<-xi
  
   iter <- 0
   step.iter <- 0
   norm.diff <-  10
  
  # cat("\n",version,", running...\n")
   while ( ((norm.diff > norm.diff.conv)||(norm.score > norm.score.conv))
     &(iter< maxiter))
   {
     
       
       

	 

	   			   #programazione quadratica programazione quadratica  programazione quadratica programazione quadratica
	   
	qpmatr<-t(Zlist$DMAT)%*%Dm%*%Zlist$DMAT
  
	   if ((is.function(d.fct)==TRUE)&(is.function(h.fct)==TRUE)) {
          Amat<-cbind(t(Zlist$DMAT)%*%Dm%*%H,t(Zlist$DMAT)%*%Dm%*%DH)
          bvec<-rbind(-h,-d ) 
          }
########################################################################
		  else {
               if (is.function(h.fct)==TRUE){ 			 
              Amat<- t(Zlist$DMAT)%*%Dm%*%H
               bvec<- -h         }
               if (is.function(d.fct)==TRUE) {
             Amat<-t(Zlist$DMAT)%*%Dm%*%DH
              bvec<- -d }
if ((is.function(d.fct)==FALSE)&(is.function(h.fct)==FALSE)) {
           Amat<-matrix(0,nrow(qpmatr),1)
           bvec<-0
}
              }
     #           qpmatr<-t(Zlist$DMAT)%*%Dm%*%Zlist$DMAT
		
if (any(is.null(qpmatr))||any(is.na(qpmatr))) {
print("matrix in quadratic programming not positive def.")
return("matrix in quadratic programming not positive def.")}	  
	   	
As<- solve.QP(qpmatr,t(Zlist$DMAT)%*%(y-m), Amat, bvec, meq=lenh, factorized=FALSE)
# -------------------------------------------------------------------------------

if(is.null(As)){
print("Error in solve.QP")
break}
#-------------------------------------------------------------------------------------	
#function used by optimize

 ff.fct<-function(steptemp){

 theta.temp <- theta + steptemp*matrix(As$solution)
        
         p  <- as.matrix(exp(Zlist$DMAT%*%theta.temp))
         #p<-p/sum(p) MODIFICA IN PROVA PER STRATI
         p<-p*c(1/Z%*%t(Z)%*%p)
         m<-p*c(Z%*%t(Z)%*%y) 
        
      if (is.function(h.fct)==FALSE) {
          h<-0}
      else{
         h  <- h.fct(m)
           }
        if (is.function(d.fct)==FALSE) {
         dd<-100 
          }
         else {
        dd <- d.fct(m)
         }  




        
        

        
         norm.score.temp <-
		 as.matrix(2/sum(y)*sum(y[y>0]*(log(y[y>0])-log(m[y>0]))))+
		 mup*sum(abs(h))
            -mup*sum(pmin(dd,dd*0))
         norm.score.temp
		 
}




stepco<-optimize(ff.fct, c(step*0.5^5, step), tol = 0.0001)
step.temp<-stepco$minimum
step.iter<-step.temp



     
         theta.temp <- theta + step.temp*matrix(As$solution)
		 
        

         norm.diff  <- sqrt(sum((theta-theta.temp)*(theta-theta.temp)))
        p  <- as.matrix(exp(Zlist$DMAT%*%theta.temp))
         #p<-(p/sum(p))  MODINPROVAPERSTRATI
          p<-p*c(1/Z%*%t(Z)%*%p)
         m<-p*c(Z%*%t(Z)%*%y)
		  
        if (is.function(h.fct)==FALSE) {
          h<-0}
        else{
         h  <- h.fct(m)
if (is.function(derht.fct)==FALSE) {
            H <- num.deriv.fct(h.fct,m)
         }
         else {
           H <- derht.fct(m)
         }

}

  if (is.function(d.fct)==FALSE) {
         d<-100 
          }
         else {
        d <- d.fct(m)
         }  


         Dm <- diag(c(m+1e-08))-((ZF*c(m))%*%t(ZF*c(p)))
       
         if (!is.matrix(Dm)) {return ("unable to reach convergence")}

        
		
        norm.score <- sum(abs(h))-sum(pmin(d,d*0))
         	 

       

       theta <- theta.temp 


       iter <- iter + 1

       if(chscore.criterion==0){
      # cat("  iter=",iter, "[",step.iter,"]", 
       #    " norm.diff=",norm.diff," norm.score=", norm.score,"\n")
                                 }
   }
}

satflag<-dim(Zlist$DMAT)[1]-dim(Zlist$DMAT)[2]

if ((is.function(h.fct)==TRUE)||(  (class(formula)=="formula")&(satflag >1  ) )){
##################

if ((is.function(h.fct)==FALSE)&(class(formula)=="formula")){

M<-cbind(Zlist$DMAT,matrix(1,nrow(Zlist$DMAT)  ))
H<-create.U(M)
H<-diag(1/c(m))%*%H
hobs<-t(H)%*%log(m)
lenh <- length(hobs)
#print(lenh)
 lam <- matrix("NA",lenh,1)

HtDHinvobs <- inv(t(H)%*%(H*c(m))) 
#H<-t(H)
}


if ((is.function(h.fct)==TRUE)&(class(formula)=="formula")){
#lenh <- length(t(H)%*%log(m))
M<-cbind(Zlist$DMAT,matrix(1,nrow(Zlist$DMAT)  ))
H2<-create.U(M)
H2<-diag(1/c(m))%*%H2
H<-cbind(H,H2)
hobs<-t(H)%*%log(m)
lenh <- length(hobs)

 lam <- matrix("NA",lenh,1)
lenh<-lenh-dim(H2)[2]
#print(lenh)
HtDHinvobs <- inv(t(H)%*%(H*c(m))) 
#H<-t(H)
}


 HtDHinv <- inv(t(H)%*%(H*c(m)))     #HtDHinv <- inv(t(H)%*%Dm%*%H)
         HHtDHinv <- H%*%HtDHinv


                                          #Ninv <- diag(c(1/Z%*%t(Z)%*%y))
   p <- m*c(1/Z%*%t(Z)%*%y)                #p <- Ninv%*%m
   resid <- y-m
   covresid <- (H*c(m))%*%HtDHinv%*%t(H*c(m)) 
                                    #covresid <- Dm%*%H%*%HtDHinv%*%t(H)%*%Dm
   covm.unadj <- covm <- Dm -  covresid
   if (sum(ZF) != 0) {
       covm <- covm.unadj -  ((ZF*c(m))%*%t(ZF*c(m)))*c(1/Z%*%t(Z)%*%y)
      #covm <- covm.unadj - Ninv%*%Dm%*%ZF%*%t(ZF)%*%Dm
   }
    covp <- t(t((covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))*
               c(1/Z%*%t(Z)%*%y))* c(1/Z%*%t(Z)%*%y)) 
   #covp <- Ninv%*%(covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))%*%Ninv
   #covp <- Ninv%*%(covm.unadj-Ninv%*%Dm%*%Z%*%t(Z)%*%Dm)%*%Ninv 
   
 # Compute adjusted residuals...
     dcovresid <- diag(covresid)
     dcovresid[abs(dcovresid)<1e-8] <- 0
     adjresid <- resid
     adjresid[dcovresid > 0] <- resid[dcovresid>0]/sqrt(dcovresid[dcovresid>0])
 # end compute adjusted residuals.

   presid <- resid/sqrt(m) 
#----------------------------
   covlam <- HtDHinv
   Gsq <- as.matrix(2*sum(y[y>0]*(log(y[y>0])-log(m[y>0]))))
   Xsq <- as.matrix(t(y-m)%*%((y-m)*c(1/m)))
   #Xsq <- as.matrix(t(y-m)%*%Dminv%*%(y-m))
###do not compute wsd if inequalities are present
if(is.function(d.fct)==FALSE){
   Wsq <- as.matrix(t(hobs)%*%HtDHinvobs%*%hobs)}
else  {Wsq<-as.matrix("NA")}
   beta <- "NA"
   covbeta <-  "NA"
   covL <-  "NA"
   L <- "NA"
   Lobs <- "NA"
   Lresid <- "NA"
   
   if (sum(abs(X)) != 0) {
       L <- L.fct(m)
       Lobs <- L.fct(y+y.eps)
       if (is.function(derLt.fct)==FALSE) {
          derLt <- num.deriv.fct(L.fct,m)
       }
       else {
          derLt <- derLt.fct(m)
       }
       PX <- inv(t(X)%*%X)%*%t(X) 
       beta <- PX%*%L 
       covL <- t(derLt)%*%covm%*%derLt
       covbeta <- PX%*%covL%*%t(PX)
       Lres <- Lobs - L
       covLres <- t(derLt)%*%covresid%*%derLt
       dcovLres <-   diag(covLres)
       dcovLres[abs(dcovLres)<1e-8] <- 0
       Lresid <- Lres
       Lresid[dcovLres > 0] <- Lres[dcovLres>0]/sqrt(dcovLres[dcovLres>0])  
       # end compute adjusted Link residuals.
       

       lbeta <- ll <- c()
       for (i in 1:length(beta)) { 
          lbeta <- c(lbeta,paste("beta",i,sep=""))
       }  
       for (i in 1:length(L)) {
          ll <- c(ll,paste("link",i,sep=""))
       } 
       dimnames(beta) <- list(lbeta,"BETA")  
     if(!is.null(colnames(X))){dimnames(beta) <- list(colnames(X),"BETA")}
       dimnames(covbeta) <- list(lbeta,lbeta)   
       dimnames(L) <- list(ll,"ML LINK") 
       dimnames(Lobs) <- list(ll,"OBS LINK") 
       dimnames(covL) <- list(ll,ll)
       dimnames(Lresid) <- list(ll,"LINK RESID")

   }
} # end "if h.fct is not the zero function"
else {
#


# else, if h.fct = 0, i.e. there are no constraints, then do...
# 
   lenh <- 0 
   lenm <- length(y)
   if(is.function( d.fct)==FALSE){
   m <- as.matrix(c(m.initial))+y.eps 
   m[m==0] <- 0.01 
   xi <- log(m) 
   Dm <- diag(c(m))
   Dminv <- diag(c(1/m))
   s <- y-m
   norm.score <- sqrt(sum(s*s))
   theta <- xi 
   lentheta <- length(theta)
   iter <- 0
   norm.diff <- 10
  # cat("\n",version,", running...\n") 
   while ( ((norm.diff > norm.diff.conv)||(norm.score > norm.score.conv))
     &(iter< maxiter))
   {     
       
       A <- Dminv      
       thetanew <- theta + step*(s*c(1/m))   #thetanew <- theta + step*A%*%s
       norm.diff <- sqrt(sum((theta-thetanew)*(theta-thetanew)))
       theta <- thetanew
       m <- exp(theta)          
       Dm <- diag(c(m))
       Dminv <- diag(c(1/m))
       s <- y-m 
       norm.score <- sqrt(sum(s*s))   
       iter <- iter + 1
      # cat("  iter=",iter, " norm.diff=",norm.diff," norm.score=", norm.score,"\n")
   } 
}
                                            #Ninv <- diag(c(1/Z%*%t(Z)%*%y)) 
   p <- m*c(1/Z%*%t(Z)%*%y)                  #p <- Ninv%*%m
   resid <- 0*y 
   covm.unadj <- covm <- Dm 
   covresid <- 0*covm 
   if (sum(ZF) != 0) {
       covm <- covm.unadj -  ((ZF*c(m))%*%t(ZF*c(m)))*c(1/Z%*%t(Z)%*%y)
      #covm <- covm.unadj - Ninv%*%Dm%*%ZF%*%t(ZF)%*%Dm
   }
   covp <- t(t((covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))*
               c(1/Z%*%t(Z)%*%y))* c(1/Z%*%t(Z)%*%y)) 
   #covp <- Ninv%*%(covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))%*%Ninv
   #covp <- Ninv%*%(covm.unadj-Ninv%*%Dm%*%Z%*%t(Z)%*%Dm)%*%Ninv 
                                            
   adjresid <-  0*y 
   presid <- 0*y
   covlam <- as.matrix(0);
#-----------------------------------
lam <- as.matrix(0)
   Gsq <- as.matrix(2*sum(y[y>0]*(log(y[y>0])-log(m[y>0]))))
   Xsq <- as.matrix(t(y-m)%*%((y-m)*c(1/m)))   
   #Xsq <- as.matrix(t(y-m)%*%Dminv%*%(y-m))
  #do not compute wsd if inequalities are present
if(is.function(d.fct)==FALSE){
 Wsq <- as.matrix(0) }
else  {Wsq<-as.matrix("NA")}
   beta <- "NA"
   covbeta <-  "NA"
   covL <-   "NA"
   L <-  "NA"
   Lresid <- "NA"
   Lobs <- "NA"
   if (sum(abs(X)) != 0) {
       L <- L.fct(m)
       Lobs <- L.fct(y)
       if (is.function(derLt.fct)==FALSE)  {
         derLt <- num.deriv.fct(L.fct,m)
       }
       else {
         derLt <- derLt.fct(m)
       }
       PX <- inv(t(X)%*%X)%*%t(X) 
       beta <- PX%*%L 
       covL <- t(derLt)%*%covm%*%derLt
       Lresid <- 0*L
       covbeta <- PX%*%covL%*%t(PX) 
       lbeta <- ll <- c() 
       for (i in 1:length(beta)) { 
          lbeta <- c(lbeta,paste("beta",i,sep=""))
       }  
       for (i in 1:length(L)) {
          ll <- c(ll,paste("link",i,sep=""))
       } 
       dimnames(beta) <- list(lbeta,"BETA")  
       if(!is.null(colnames(X))){dimnames(beta) <- list(colnames(X),"BETA")}
      dimnames(covbeta) <- list(lbeta,lbeta)
       dimnames(Lobs) <- list(ll,"OBS LINK")   
       dimnames(L) <- list(ll,"ML LINK")  
       dimnames(covL) <- list(ll,ll) 
       dimnames(Lresid) <- list(ll,"LINK RESID") 
   }
}# end "else if h.fct = 0"
#
# ASSIGN LABELS...
#
  lm <- ly <- lp <- lbeta <- lr <- lar <- lpr <- ll <- llam <- c()
 


if(is.null(rownames(y))){
  for (i in 1:lenm) {
     lm <- c(lm,paste("m",i,sep=""))
     ly <- c(ly,paste("y",i,sep=""))
     lp <- c(lp,paste("p",i,sep=""))
     lr <- c(lr,paste("r",i,sep=""))
     lar <- c(lar,paste("adj.r",i,sep=""))
     lpr <- c(lpr,paste("pearson.r",i,sep=""))
  }
 }

else{ly<-c( paste("y(",rownames(y),")"))
     lm<- c( paste("m(",rownames(y),")"))
     lp <-c( paste("p(",rownames(y),")"))
     lr<- c(paste("r(",rownames(y),")"))
     lar<-c(paste("a.r(",rownames(y),")"))
     lpr<-c(paste("p.r(",rownames(y),")"))}





  for (i in 1:length(lam)) {
     llam <- c(llam,paste("lambda",i,sep=""))
  } 
  dimnames(y) <- list(ly,"OBS")
  dimnames(m) <- list(lm,"FV")
  dimnames(p) <- list(lp,"PROB")
  dimnames(resid) <- list(lr,"RAW RESIDS")
  dimnames(presid) <- list(lpr,"PEARSON RESIDS")
  dimnames(adjresid) <- list(lar, "ADJUSTED RESIDS")  
  dimnames(lam) <- list(llam,"LAGRANGE MULT") 
  dimnames(covm) <- list(lm,lm)
  dimnames(covp) <- list(lp,lp)
  dimnames(covresid) <- list(lr,lr)
  dimnames(covlam) <- list(llam,llam)
  dimnames(Xsq) <- list("","PEARSON SCORE STATISTIC")
  dimnames(Gsq) <- list("","LIKELIHOOD RATIO STATISTIC")
  dimnames(Wsq) <- list("","GENERALIZED wALD STATISTIC")

if (is.function(derht.fct)==FALSE) {derht.fct <- "Numerical derivatives used."}
if (is.function(derLt.fct)==FALSE) {derLt.fct <- "Numerical derivatives used."}
#cat("\n")
#cat(" Time Elapsed:", proc.time()[3]-start.time,"seconds")
#cat("\n")
lenh<-lenh*(is.function(h.fct)==TRUE)
modlist<-list(y=y,m=m,covm=covm,p=p,covp=covp, 
lambda=lam,covlambda=covlam,
resid=resid,presid=presid,adjresid=adjresid,covresid=covresid,
Gsq=Gsq,Xsq=Xsq,Wsq=Wsq,df=lenh,
beta=beta,covbeta=covbeta, Lobs=Lobs, L=L,covL=covL,Lresid=Lresid,
iter=iter, norm.diff=norm.diff,norm.score=norm.score,
h.fct=h.fct,derht.fct=derht.fct,L.fct=L.fct,derLt.fct=derLt.fct,
d.fct=d.fct,derdt.fct=derdt.fct,
X=X,U=U,Z=Z,ZF=ZF,Zlist=Zlist,version=version)
class(modlist)="mphfit"
modlist
}




