
## require(MCMCpack);

splitParameters<-function(paramList){
  param.fixed=list();
  param.free.init =list();
  param.free.weight =list();
  param.free.lower =list();
  param.free.upper =list();
  for(iter in seq(along=paramList)){
    if(class(paramList[[iter]])!="clineParameter")
      stop(paste("Item",paramList[[iter]],"not a cline parameter!"));
    parName<-attr(paramList[[iter]],"param");
    parValue<-paramList[[iter]]$val;
    
    if(attr(paramList[[iter]],"fixed")){
      param.fixed[[parName]]<-parValue;
    }else{
      param.free.init[[parName]]<-parValue;
      param.free.weight[[parName]]<-paramList[[iter]]$w;
      param.free.lower[[parName]]<-attr(paramList[[iter]],"limit.lower");
      param.free.upper[[parName]]<-attr(paramList[[iter]],"limit.upper");
    }
  }
  param<-list(init=param.free.init,tune=param.free.weight,
              lower=param.free.lower,upper=param.free.upper,
              fixed=param.fixed);
  return(param);
}
