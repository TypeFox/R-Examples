ICBayes.default <-
function(
L,
R,
model,
status,
xcov,
x_user=NULL,
order=2,
sig0=10,
coef_range=5,
v0=0.1,
a_eta=1,
b_eta=1,
knots=NULL,
grids=NULL,
conf.int=0.95,
plot_S=TRUE,
chain.save=FALSE,
dd1,
niter=11000,
burnin=1000,
thin=1,
...){

   if (model=='case2ph'){
est<-case2ph(L,R,status,xcov,x_user,order,sig0,coef_range,
a_eta,b_eta,knots,grids,niter)
   } else
   if (model=='case1ph'){
est<-case1ph(L,R,status,xcov,x_user,order,sig0,coef_range,
a_eta,b_eta,knots,grids,niter)
   } else
   if (model=='case1po'){
est<-case1po(L,R,status,xcov,x_user,order,sig0,coef_range,
a_eta,b_eta,knots,grids,niter)
   } else
   if (model=='case2probit'){
est<-case2probit(L,R,status,xcov,x_user,order,
v0,a_eta,b_eta,knots,grids,niter)
   }

   p=ncol(as.matrix(xcov))
   if (chain.save==TRUE){
write(t(est$parbeta),dd1,ncolumns=p)
   }
   wbeta=as.matrix(est$parbeta[seq((burnin+thin),niter,by=thin),],ncol=p) # thinned beta samples
   wparsurv0=as.matrix(est$parsurv0[seq((burnin+thin),niter,by=thin),],ncol=length(grids))
   wparsurv=as.matrix(est$parsurv[seq((burnin+thin),niter,by=thin),],ncol=length(grids))
   coef<-apply(wbeta,2,mean)
   coef_ssd<-apply(wbeta,2,sd)
   credinterv<-array(rep(0,p*2),dim=c(p,2))
   colnames(credinterv)<-c(paste(100*(1-conf.int)/2,"%CI"),
paste(100*(0.5+conf.int/2),"%CI"))
   for (r in 1:p){
credinterv[r,]<-quantile(wbeta[,r],c((1-conf.int)/2,0.5+conf.int/2))
   }
   CPO=1/apply(est$parfinv[seq((burnin+thin),niter,by=thin),],2,mean)
   LPML=sum(log(CPO))      # bigger is better   

   if (plot_S==FALSE){
   output<-list(coef = coef,
    coef_ssd = coef_ssd,
    coef_ci = credinterv,
    LPML = LPML)} else {
   output<-list(coef = coef,
    coef_ssd = coef_ssd,
    coef_ci = credinterv,
    grids = est$grids,
    S0_m = apply(wparsurv0,2,mean),
    S_m = apply(wparsurv,2,mean),
    LPML = LPML)}
   output$call<-match.call()
   class(output)<-"ICBayes"
   output
   }
