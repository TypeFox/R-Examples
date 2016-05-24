pssm.power<-function(parameter=c("N","Delta","Power","Power-Simulation"),
                     endpoint=c("Progression","Survival"),
                     progression.delta=.2,#approx %improvement or treatment
                     nsamp=100,alpha=0.025,pow=0.80,
                     tsurv=3,h0=0,prior.surv=0,
                     alloc=c(1,1),phaz=log(.2),shaz=log(.2),
                     accrual=3,followup=2,inc=0.05,
                     delta=0,intervals=1,nsims=500,n=5000,seed=NULL){
  
  #note assumes progression.delta is negative indicated less hazard of progression
  m=accrual+followup
  realparms=c(shaz*rep(1,m*(m+1)/2),phaz*rep(1,m),-progression.delta,0)
  #browser()
  methods=c("N","Delta","Power","Power-Simulation")
 
  visits=seq(from=inc,to=m-inc,by=inc)
  df1=pssm.generate.data(theta1=-progression.delta,theta2=0,
                         phaz.progression=(phaz)*rep(1,m),
                         phaz.survival=(shaz)*rep(1,m*(m+1)/2),
                         accrual=accrual,followup=followup,
                         m=m,n=n,times=visits,delta=delta,alloc=alloc,
                         seed=seed)
  
  ps=pssm(Surv(tprog0,tprog1)~rx,Surv(tdeath,cdeath)~rx,
          dat=df1,intervals=intervals,start=NULL,rescale=1)
  vst=data.frame(parameter="Progression",endpoint="Proportion",
                 value=mean(is.na(df1$tprog1)))
  
  vst=rbind(vst,data.frame(parameter="Death",endpoint="Proportion",
                 value=mean(df1$cdeath))
                 )
#   fr=alloc[1]/sum(alloc)
#   fr2=fr*(1-fr)
  cc=matrix(c(0,1),2,1)
  mul=matrix(c(-1,1),1,2)
  isi<-function(i,j){parameter[i]==methods[j]}
  for (i in 1:length(parameter)){
    for (j in 1:length(endpoint)){
        prog=(endpoint[j]=="Progression")
       if(prog)  {
        sig1=n*(ps@se.estimates.progression)^2
        delta1=progression.delta
      } else {
          us=pssm.survivalcurv(ps,cov1=cc,cov2=cc)
          vsf<-function(ns,h0) {
         
            vs=us(tsurv,prior.surv*n/ns)
            delta1=mul%*%vs$estimate+h0
            sig1=(n/ns)*(mul%*%attributes(vs)$covariance%*%t(mul)) 
            return(c(delta1/sqrt(sig1),((qnorm(1-alpha)+qnorm(pow))^2)*sig1/(delta1^2)))            
          }
      } 
      
      if (isi(i,1)){
        if (prog) v1=((qnorm(1-alpha)+qnorm(pow))^2)*sig1/(delta1^2)
        else v1=uniroot(function(x) vsf(x,h0)[2]-1,c(1,100000))$root    
       vst=rbind(vst,data.frame(parameter=methods[1],endpoint=endpoint[j],value=v1))
      }
      
      if (isi(i,2)){
  
        if (prog) lab="Delta" else lab="HO"
        if (prog) v1=(qnorm(1-alpha)+qnorm(pow))*sqrt(sig1/(nsamp))
        else v1=uniroot(function(x) vsf(nsamp,x)[2]-1,c(-.05,.05),extendInt="yes")$root
        vst=rbind(vst,data.frame(parameter=lab,endpoint=endpoint[j],value=v1))
        }
          
      if (isi(i,3)){
        if (prog) v1=pnorm(delta1/sqrt(sig1/(nsamp))-qnorm(1-alpha))
        else v1=pnorm(vsf(nsamp,h0)[1]-qnorm(1-alpha))
        vst=rbind(vst,data.frame(parameter=methods[3],endpoint=endpoint[j],value=v1)) 
          }
      if(isi(i,4)&&(!(is.element("Power-Simulation",vst$parameter)))){
        k=qnorm(alpha)
        st1=NULL
        powprog=0
        powsurv=0
        #rs=matrix(0,nsims,4)
        for(isim in 1:nsims){

          df1=pssm.generate.data(theta1=-progression.delta,theta2=0,
                                 phaz.progression=(phaz)*rep(1,m),
                                 phaz.survival=(shaz)*rep(1,m*(m+1)/2),
                                 accrual=accrual,followup=followup,
                                 m=m,n=nsamp,times=visits,delta=delta)
          ps=pssm(Surv(tprog0,tprog1)~rx,Surv(tdeath,cdeath)~rx,
                  dat=df1,intervals=intervals,
                  start=st1,rescale=1)
          st1=ps@estimates
          powprog=powprog+(-ps@estimates.progression/ps@se.estimates.progression>-k)
          us=pssm.survivalcurv(ps,cov1=cc,cov2=cc)
          vs=us(tsurv,prior.surv)
          delta1=mul%*%vs$estimate+h0
          powsurv=powsurv+(delta1/sqrt(mul%*%attributes(vs)$covariance%*%t(mul))>-k)
#           rs[isim,]=c(delta1,
#                       vdelta1=mul%*%attributes(vs)$covariance%*%t(mul),
#                       delta1/sqrt(mul%*%attributes(vs)$covariance%*%t(mul)),
#                        -ps@estimates.progression/ps@se.estimates.progression
#                        )
        
        }
        vst=rbind(vst,data.frame(parameter=methods[4],endpoint='Progression',
                                  value=powprog/nsims))
        vst=rbind(vst,data.frame(parameter=methods[4],endpoint='Survival',
                                  value=powsurv/nsims))
        
#         vst=rbind(vst,data.frame(parameter="mean",endpoint="delta1",
#                                  value=mean(rs[,1])))
#         vst=rbind(vst,data.frame(parameter="SD",endpoint='delta1',
#                                  value=sqrt(var(rs[,1]))))
#         vst=rbind(vst,data.frame(parameter="Exp SD",endpoint='delta1',
#                                  value=mean(sqrt(rs[,2]))))
#         vst=rbind(vst,data.frame(parameter="Exp SD",endpoint='delta1',
#                                  value=mean(sqrt(rs[,3]))))
#         vst=rbind(vst,data.frame(parameter="Exp SD",endpoint='delta1',
#                                  value=sqrt(var(rs[,3]))))
        
      }}}
  #return(list(vst,rs))
   return(vst)
}
#vu=pssm.power(parameter=c("N","Delta","Power"),endpoint=c("Progression","Survival"),progression.delta=.405,nsamp=500,h0=.04)
