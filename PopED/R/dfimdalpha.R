
dfimdalpha <- function(alpha, model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped.db,ha){
  bpop=bpopdescr[,2,drop=F]
  bpop[bpopdescr[,1,drop=F]!=0]=alpha[1:sum(bpopdescr[,1,drop=F]!=0)]
  d=ddescr[,2,drop=F]
  d[ddescr[,1]!=0]=alpha[(sum(bpopdescr[,1,drop=F]!=0)+1):length(alpha)]
  d=getfulld(d,covd)
  fim=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,poped.db)
  fim <- fim$ret
  grad=array(0,dim=c(size(fim,1),size(fim,2),length(alpha)))
  for(i in 1:length(alpha)){
    alpha_plus=alpha
    alpha_plus[i]=alpha_plus[i]+ha
    bpop=bpopdescr[,2,drop=F]
    bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
    d=ddescr[,2,drop=F]
    d[ddescr[,1]!=0]=alpha_plus[(sum(bpopdescr[,1,drop=F]!=0)+1):length(alpha_plus)]
    d=getfulld(d,covd)
    fim_plus=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,poped.db)
    fim_plus <- fim_plus$ret
    if((i>sum(bpopdescr[,1]!=0))){
      grad[,,i]=(fim_plus-fim)/ha
    } else {
      #central differences for fixed effects
      alpha_minus=alpha
      alpha_minus[i]=alpha_minus[i]-ha
      bpop=bpopdescr[,2,drop=F]
      bpop[bpopdescr[,1,drop=F]!=0]=alpha_minus[1:sum(bpopdescr[,1,drop=F]!=0)]
      fim_minus=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,poped.db)
      fim_minus <- fim_minus$ret
      grad[,,i]=(fim_plus-fim_minus)/(2*ha)
    }
  }
  return(list( grad= grad,fim =fim )) 
}
  