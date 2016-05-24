robustRegBS<-function(formula,data,tune=4.685,m=TRUE,max.it=1000,tol=1e-5,anova.table=FALSE){

bi<-FALSE
if(m==FALSE){bi<-TRUE}


modelFrame=model.frame(formula,data)
X=model.matrix(formula,data)
y=model.extract(modelFrame,"response")

model=lm.fit(X,y)
b=model$coefficients
n<-length(y)
p<-length(b)
if(bi){
 tune<-(tune*sqrt(2*p*n))/(n-2*p)
 hii<-lm.influence(model)$hat
 pi<-(1-hii)/sqrt(hii)}

convergence<-FALSE
for(i in 1:max.it){
 b_temp<-b
 r<-y-fit_rcpp(X,b) #replaced r<-y-X%*%b
 s<-mad_rcpp(r) #replaced s<-median(abs(r-median(r)))/.6745

 if(m){rstar<-(r/s)}else{rstar<-r/(s*pi)}
 
 psiBS<-psiBS_rcpp(rstar,tune)
 w<-psiBS/rstar #replaced 
 b<-lm.wfit(x=X,y=y,w=w[,1])$coefficients

 
 if(i>4){
  if(sum(abs(b-b_temp))<tol){
   cat("\nRobust Regression with Bisquare Function\n")
   cat("Convergence achieved after:",i,"iterations\n")
   convergence<-TRUE
   break}}
}

if(convergence==FALSE){b<-NULL;w<-NULL}
if(convergence && anova.table){

derivPsiBiSquare<-function(r,c){
true<-abs(r)<=c
false<-abs(r)>c
psi<-(1-r^2/c^2)*(1-(5*r^2)/c^2)*true+false*0
return(psi)
}
	
r3<-function(x){return(round(x,3))}
r2<-function(x){return(round(x,2))}

 ybarw<-sum(y*w)/sum(w)
 ytild=fit_rcpp(X,b) #replaced ytild<-X%*%b
 ssreg<-sum(w*(ytild-ybarw)^2)
 sserr<-sum(w*(y-ytild)^2)
 sstot<-sum(w*(y-ybarw)^2)
 dfr<-length(b)-1
 dferr<-length(y)-dfr-1
 dftot<-length(y)-1
 msr<-ssreg/dfr
 mse<-sserr/dferr
 sbsq<-(s^2*(n^2/(n-length(b)))*sum(psiBS^2))/sum(derivPsiBiSquare(rstar,tune))^2
 F<-msr/sbsq
 W<-Diagonal(x=w[1:length(w),]) 
 c<-diag(x=solve(a=t(X)%*%W%*%X))
 sec<-sqrt(sbsq*c)
 t<-b/sec
 
 cat("source  ","\t","SS","\t","\t","df","\t","MS","\t","\t","F","\n")
 cat("model   ","\t",r2(ssreg),"\t",dfr,"\t",r2(msr),"\t",r2(F),"\n")
 cat("error   ","\t",r2(sserr),"\t",dferr,"\n")
 cat("tot     ","\t",r2(sstot),"\t",dftot,"\n")
 cat("rsquared","\t",r3(ssreg/sstot),"\n")
 cat("mse     ","\t",mse,"\n")
 estimates<-cbind(b,sec,r2(t))
 #colnames(estimates)<-c("estimate","se","t")
 print(estimates)
 }

 object=list("coefficients"=b,"weights"=w)
 return(invisible(object)) 
}

