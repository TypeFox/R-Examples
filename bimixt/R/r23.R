##################
# Bimixt Package #
##################


### Terminology
### There are two scales that are relevant: the original data scale (OS) and the transformed scale (TS)
### The "lower" component refers to the component with the smaller mean
### The "upper" component refers to the component with the greater mean


library(pROC)

bimixt.model <- function(case,control,type="binorm",start.vals=NULL){
  # binorm: binormal model
  # 4c: four component model
  # 2cu: two component unconstrained model
  # 2cc: two component constrained model
  case=na.omit(case)
  control=na.omit(control)
  
  n.cs=length(case)
  n.ctrl=length(control)
  
  if(n.cs<3 || n.ctrl<3){stop("there are insufficient data to fit the model.")}
  
  types=c("binorm","4c","2cu","2cc")
  if(is.null(case)||is.null(control)){stop("must specify case and control values")}
  if(!(type %in% types)){stop("model type not supported")}
  if(!is.numeric(case) || !is.numeric(control)){stop("case and control values must be numeric")}
  
  if(type=="binorm"){
    if(!is.null(start.vals)){warning("start values not used for binorm model")}
    out <- bc.binorm(case,control)
    class(out) <- "model"
    if((n.cs+n.ctrl)<=5){warning("The number of data points is less than or equal to the number of free parameters.")}
  }
  
  if(type=="4c"){
    if(!is.null(start.vals)){
      if(!is.vector(start.vals$pi.cs)||!is.vector(start.vals$sig.cs)||!is.vector(start.vals$mu.cs)||!is.vector(start.vals$pi.ctrl)||!is.vector(start.vals$sig.ctrl)||!is.vector(start.vals$mu.ctrl)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(!is.numeric(start.vals$pi.cs)||!is.numeric(start.vals$sig.cs)||!is.numeric(start.vals$mu.cs)||!is.numeric(start.vals$pi.ctrl)||!is.numeric(start.vals$sig.ctrl)||!is.numeric(start.vals$mu.ctrl)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(is.null(start.vals$pi.cs)||is.null(start.vals$mu.cs)||is.null(start.vals$sig.cs)||is.null(start.vals$pi.ctrl)||is.null(start.vals$mu.ctrl)||is.null(start.vals$sig.ctrl)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(length(start.vals$pi.cs)!=2||length(start.vals$mu.cs)!=2||length(start.vals$sig.cs)!=2||length(start.vals$pi.ctrl)!=2||length(start.vals$mu.ctrl)!=2||length(start.vals$sig.ctrl)!=2){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(sum(start.vals$pi.cs)!=1||sum(start.vals$pi.ctrl)!=1){stop("Starting proportions must sum to 1.")}
      start.vals.controls=list(mu=start.vals$mu.ctrl,sig=start.vals$sig.ctrl,pi=start.vals$pi.ctrl)
      start.vals.cases=list(mu=start.vals$mu.cs,sig=start.vals$sig.cs,pi=start.vals$pi.cs)
    }
    
    if(is.null(start.vals)){
      start.vals.controls=NULL
      start.vals.cases=NULL
    }
    
    out <- bc.fourcomp(case,control,start.vals.cases=start.vals.cases,start.vals.controls=start.vals.controls)
    
    class(out) <- "model"
    if((n.cs+n.ctrl)<=11){warning("The number of data points is less than or equal to the number of free parameters.")}
  }
  
  if(type=="2cu"){
    if(!is.null(start.vals)){
      if(!is.numeric(start.vals$pi.ctrl)||!is.vector(start.vals$pi.cs)||!is.vector(start.vals$sig)||!is.vector(start.vals$mu)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(!is.numeric(start.vals$pi.ctrl)||!is.numeric(start.vals$pi.cs)||!is.numeric(start.vals$sig)||!is.numeric(start.vals$mu)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(is.null(start.vals$pi.cs)||is.null(start.vals$pi.ctrl)||is.null(start.vals$mu)||is.null(start.vals$sig)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(length(start.vals$pi.cs)!=2||length(start.vals$pi.ctrl)!=2||length(start.vals$mu)!=2||length(start.vals$sig)!=2){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(sum(start.vals$pi.cs)!=1 || sum(start.vals$pi.ctrl)!=1){stop("Starting proportions must sum to 1.")}
    }
    out <- bc.twocomp(case,control,constrained=F,start.vals=start.vals)
    class(out) <- "model"
    if((n.cs+n.ctrl)<=7){warning("The number of data points is less than or equal to the number of free parameters.")}
  }
  
  if(type=="2cc"){
    if(!is.null(start.vals)){
      if(!is.vector(start.vals$pi)||!is.vector(start.vals$sig)||!is.vector(start.vals$mu)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(!is.numeric(start.vals$pi)||!is.numeric(start.vals$sig)||!is.numeric(start.vals$mu)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(is.null(start.vals$pi)||is.null(start.vals$mu)||is.null(start.vals$sig)){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(length(start.vals$pi)!=2||length(start.vals$mu)!=2||length(start.vals$sig)!=2){stop("Invalid start values provided. Start values must be a list of numeric vectors of length 2. List names must match those given in help(bimixt.model).")}
      if(sum(start.vals$pi)!=1){stop("Starting proportions must sum to 1.")}
    }
    out <- bc.twocomp(case,control,constrained=T,start.vals=start.vals)
    class(out) <- "model"
    if((n.cs+n.ctrl)<=6){warning("The number of data points is less than or equal to the number of free parameters.")}
  }
  out
}

mn <- function(model,transformed=F){
  if(!transformed){out=list(untransformed.cases=model$mu.cases.unt,untransformed.controls=model$mu.controls.unt)}
  if(transformed){out=list(transformed.cases=model$mu.cases,transformed.controls=model$mu.controls)}  
out
  }

stdev <- function(model,transformed=F){
  if(!transformed){out=list(untransformed.cases=model$sig.cases.unt,untransformed.controls=model$sig.controls.unt)}
  if(transformed){out=list(transformed.cases=model$sig.cases,transformed.controls=model$sig.controls)}  
out
  }


prop <- function(model) list(cases=model$pi.cases,controls=model$pi.controls)

lambda <- function(model) model$lambda

type <- function(model) model$type

maxll <- function(model) model$max.loglike

#### summary method
summary.model <- function(object,...){
  model=object
  if(model$type=="binorm")
  {
    df=data.frame(matrix(c(model$mu.cases.unt,model$sig.cases.unt,model$mu.cases,model$sig.cases,1,0,
                           model$mu.controls.unt,model$sig.controls.unt,model$mu.controls,model$sig.controls,0,1),nrow=2,ncol=6,byrow=T))
    df=df[order(df$X1),]
    names(df)=c("Untrans. Mean","Untrans. Stand. Dev.","Trans. Mean","Trans. Stand. Dev.","Case Prop.","Control Prop.")
    row.names(df)=c("component 1", "component 2")
  }
  
  if(model$type=="2cc")
  {
    n=which(model$mu.cases==model$mu.controls)
    df=data.frame(matrix(c(model$mu.cases.unt[n],model$sig.cases.unt[n],model$mu.cases[n],model$sig.cases[n],model$pi.cases[n],1,
                           model$mu.cases.unt[-n],model$sig.cases.unt[-n],model$mu.cases[-n],model$sig.cases[-n],model$pi.cases[-n],0),nrow=2,ncol=6,byrow=T))
    df=df[order(df$X1),]
    names(df)=c("Untrans. Mean","Untrans. Stand. Dev.","Trans. Mean","Trans. Stand. Dev.","Case Prop.","Control Prop.")
    row.names(df)=c("component 1", "component 2")
  }
  
  if(model$type=="2cu")
  {
    n=which(model$mu.controls==model$mu.cases[1])
    df=data.frame(matrix(c(model$mu.cases.unt[1],model$sig.cases.unt[1],model$mu.cases[1],model$sig.cases[1],model$pi.cases[1],model$pi.controls[n],
                           model$mu.cases.unt[-1],model$sig.cases.unt[-1],model$mu.cases[-1],model$sig.cases[-1],model$pi.cases[-1],model$pi.controls[-n]),nrow=2,ncol=6,byrow=T))
    df=df[order(df$X1),]
    names(df)=c("Untrans. Mean","Untrans. Stand. Dev.","Trans. Mean","Trans. Stand. Dev.","Case Prop.","Control Prop.")
    row.names(df)=c("component 1", "component 2")
  }
  
  if(model$type=="4c")
  {
    df=data.frame(matrix(c(model$mu.cases.unt[1],model$sig.cases.unt[1],model$mu.cases[1],model$sig.cases[1],model$pi.cases[1],0,
                           model$mu.cases.unt[-1],model$sig.cases.unt[-1],model$mu.cases[-1],model$sig.cases[-1],model$pi.cases[-1],0,
                           model$mu.controls.unt[1],model$sig.controls.unt[1],model$mu.controls[1],model$sig.controls[1],0,model$pi.controls[1],
                           model$mu.controls.unt[-1],model$sig.controls.unt[-1],model$mu.controls[-1],model$sig.controls[-1],0,model$pi.controls[-1]
    ),nrow=4,ncol=6,byrow=T))
    df=df[order(df$X1),]
    names(df)=c("Untrans. Mean","Untrans. Stand. Dev.","Trans. Mean","Trans. Stand. Dev.","Case Prop.","Control Prop.")
    row.names(df)=c("component 1", "component 2", "component 3", "component 4")
  }
  return(df)
}

#### plot method
plot.model <- function(x,histogram=T, breaks = "Sturges",main=model$type, cols=c("#008ED6","#990033"),ylab="Density",xlab="",...){
  model=x
  
  case=model$case
  control=model$control
  prop.cs=length(case)/length(c(case,control))
  prop.ctrl=length(control)/length(c(case,control))
  
  h=hist(c(case,control),breaks=breaks,plot=F,warn.unused=F)
  x.lb=min(h$breaks)
  x.ub=max(h$breaks)
  lam=lambda(model)
  x=seq(x.lb,x.ub,length=1000)
  
  h.cs=hist(case,breaks=h$breaks,plot=F,warn.unused=F)
  h.cs$density=prop.cs*h.cs$density
  h.ctrl=hist(control,breaks=h$breaks,plot=F,warn.unused=F)
  h.ctrl$density=prop.ctrl*h.ctrl$density
  
  if(model$type == "binorm"){
    
    mu1=model$mu.cases
    mu3=model$mu.controls
    sig1=model$sig.cases
    sig3=model$sig.controls
    p1=prop.cs*1
    p2=prop.ctrl*1
    
    cases=p1*boxcox.inv.density(x,lam,mu1,sig1)
    controls=p2*boxcox.inv.density(x,lam,mu3,sig3) 
  }
  
  if(model$type == "2cc"){
    
    mu1=model$mu.controls 
    mu3=model$mu.cases[1]
    mu4=model$mu.cases[2]
    sig1=model$sig.controls
    sig3=model$sig.cases[1]
    sig4=model$sig.cases[2]
    p1=prop.ctrl*model$pi.controls
    p3=prop.cs*model$pi.cases[1]
    p4=prop.cs*model$pi.cases[2]
    
    controls=p1*boxcox.inv.density(x,lam,mu1,sig1)
    cases=p3*boxcox.inv.density(x,lam,mu3,sig3)+p4*boxcox.inv.density(x,lam,mu4,sig4)
  }
  
  if(model$type == "2cu" || model$type == "4c"){
    
    mu1=model$mu.controls[1]
    mu2=model$mu.controls[2] 
    mu3=model$mu.cases[1]
    mu4=model$mu.cases[2]
    sig1=model$sig.controls[1]
    sig2=model$sig.controls[2]
    sig3=model$sig.cases[1]
    sig4=model$sig.cases[2]
    p1=prop.ctrl*model$pi.controls[1]
    p2=prop.ctrl*model$pi.controls[2]
    p3=prop.cs*model$pi.cases[1]
    p4=prop.cs*model$pi.cases[2]
    
    controls=p1*boxcox.inv.density(x,lam,mu1,sig1)+p2*boxcox.inv.density(x,lam,mu2,sig2)
    cases=p3*boxcox.inv.density(x,lam,mu3,sig3)+p4*boxcox.inv.density(x,lam,mu4,sig4) 
  }
  
  d.y.ub=range(c(cases,controls),finite=T)[2] 
  y.lb=0
  
  if(histogram){
    h.y.ub=max(c(h.cs$density,h.ctrl$density))
    y.ub=max(h.y.ub,d.y.ub)
    plot(h.ctrl,freq=F,xlim=c(x.lb,x.ub),ylim=c(y.lb,y.ub),main=paste(main),ylab=paste(ylab),xlab=xlab,border=cols[1])
    plot(h.cs,freq=F,border=cols[2],add=T)
    lines(x,controls,type="l",col=cols[1])
  }
  if(!histogram){
    plot(x,controls,type="l",main=paste(main),ylab=paste(ylab),xlab=xlab,xlim=c(x.lb,x.ub), ylim=c(y.lb,d.y.ub),col=cols[1])
  }
  lines(x,cases,col=cols[2])    
}


#### print method
print.model <- function(x,...){
  model=x
  if(model$type == "binorm"){
    out <- matrix(c("Model Type",format(type(model)),
                    "","Untransformed Case Mean", format(model$mu.cases.unt),
                    "","Transformed Case Means", format(model$mu.cases),
                    "","Untransformed Case Standard Deviation",format(model$sig.cases.unt),
                    "","Transformed Case Standard Deviation",format(model$sig.cases),
                    "","Case Proportion",format(model$pi.cases),
                    "","Untransformed Control Mean",paste0(format(model$mu.controls.unt)),
                    "","Transformed Control Mean",paste0(format(model$mu.controls)),
                    "","Untransformed Control Standard Deviation", paste0(format(model$sig.controls.unt)),
                    "","Transformed Control Standard Deviation", paste0(format(model$sig.controls)),
                    "","Control Proportion", paste0(format(model$pi.controls)),
                    "","Lambda",format(lambda(model)),
                    "","Maximum Loglikelihood",format(model$max.loglike)),nrow=38,ncol=1,dimnames=list(rep("",38),""))
  }
  if(model$type == "2cc"){
    out <- matrix(c("Model Type",format(type(model)),
                    "","Untransformed Case Means", paste0(format(model$mu.cases.unt[1]),", ",format(model$mu.cases.unt[2])),
                    "","Transformed Case Means", paste0(format(model$mu.cases[1]),", ",format(model$mu.cases[2])),
                    "","Untransformed Case Standard Deviations",paste0(format(model$sig.cases.unt[1]),", ",format(model$sig.cases.unt[2])),
                    "","Transformed Case Standard Deviations",paste0(format(model$sig.cases[1]),", ",format(model$sig.cases[2])),
                    "","Case Proportions",paste0(format(model$pi.cases[1]),", ",format(model$pi.cases[2])),
                    "","Untransformed Control Mean",paste0(format(model$mu.controls.unt)),
                    "","Transformed Control Mean",paste0(format(model$mu.controls)),
                    "","Untransformed Control Standard Deviation", paste0(format(model$sig.controls.unt)),
                    "","Transformed Control Standard Deviation", paste0(format(model$sig.controls)),
                    "","Control Proportion", paste0(format(model$pi.controls)),
                    "","Lambda",format(lambda(model)),
                    "","Maximum Loglikelihood",format(model$max.loglike)),nrow=38,ncol=1,dimnames=list(rep("",38),""))
  }
  
  if(model$type == "2cu" || model$type == "4c"){
    out <- matrix(c("Model Type",format(type(model)),
                    "","Untransformed Case Means",paste0(format(model$mu.cases.unt[1]),", ",format(model$mu.cases.unt[2])),
                    "","Transformed Case Means",paste0(format(model$mu.cases[1]),", ",format(model$mu.cases[2])),
                    "","Untransformed Case Standard Deviations", paste0(format(model$sig.cases.unt[1]),", ",format(model$sig.cases.unt[2])), 
                    "","Transformed Case Standard Deviations",paste0(format(model$sig.cases[1]),", ",format(model$sig.cases[2])),
                    "","Case Proportions",paste0(format(model$pi.cases[1]),", ",format(model$pi.cases[2])),
                    "","Untransformed Control Means",paste0(format(model$mu.controls.unt[1]),", ",format(model$mu.controls.unt[2])),
                    "","Transformed Control Means",paste0(format(model$mu.controls[1]),", ",format(model$mu.controls[2])),
                    "","Untransformed Control Standard Deviations", paste0(format(model$sig.controls.unt[1]),", ",format(model$sig.controls.unt[2])), 
                    "","Transformed Control Standard Deviations", paste0(format(model$sig.controls[1]),", ",format(model$sig.controls[2])),
                    "","Control Proportions",paste0(format(model$pi.controls[1]),", ",format(model$pi.controls[2])),
                    "","Lambda",format(lambda(model)),
                    "","Maximum Loglikelihood",format(model$max.loglike)
    ),nrow=38,ncol=1,dimnames=list(rep("",38),""))
  }
  
  print(out,quote=F)
}


#### Likelihood Ratio Test
lr.test <- function(model1,model2){
  
  if(model1$control != model2$control || model1$case != model2$case) stop("The data are not consistent across models. LR test is not applicable.")
  if(model1$type==model2$type) warning("models are equivalent")
  
  my.chi <- function(ll0,lla,df){
    test.stat=2*abs(lla-ll0)
    p.val=1-pchisq(test.stat,df)
    return(p.val)
  }
  
  if(model1$type=="binorm"){
    
    if(model2$type=="binorm"){
      p=my.chi(model1$max.loglike,model2$max.loglike,0)
    }
    
    if(model2$type=="2cc"){
      if(model1$max.loglike>model2$max.loglike){p=NA}
      else{p=my.chi(model1$max.loglike,model2$max.loglike,1)}
    }
    
    if(model2$type=="2cu"){
      if(model1$max.loglike>model2$max.loglike){p=NA}
      else{p=my.chi(model1$max.loglike,model2$max.loglike,2)}
    }
    
    if(model2$type=="4c"){
      if(model1$max.loglike>model2$max.loglike){p=NA}
      else{p=my.chi(model1$max.loglike,model2$max.loglike,6)}
    }
  }
  
  
  if(model1$type=="2cc"){
    if(model2$type=="binorm"){
      if(model1$max.loglike<model2$max.loglike){p=NA}
      else{p=my.chi(model2$max.loglike,model1$max.loglike,1)}
    }
    
    if(model2$type=="2cc"){
      p=my.chi(model1$max.loglike,model2$max.loglike,0)
    }
    
    if(model2$type=="2cu"){
      if(model1$max.loglike>model2$max.loglike){p=NA}
      else{p=my.chi(model1$max.loglike,model2$max.loglike,1)}
    }
    
    if(model2$type=="4c"){
      if(model1$max.loglike>model2$max.loglike){p=NA}
      else{p=my.chi(model1$max.loglike,model2$max.loglike,5)}
    }
    
  }
  
  if(model1$type=="2cu"){
    if(model2$type=="binorm"){
      if(model1$max.loglike<model2$max.loglike){p=NA}
      else{p=my.chi(model2$max.loglike,model1$max.loglike,2)}
    }
    
    if(model2$type=="2cc"){
      if(model1$max.loglike<model2$max.loglike){p=NA}
      else{p=my.chi(model2$max.loglike,model1$max.loglike,1)}
    }
    
    if(model2$type=="2cu"){
      p=my.chi(model2$max.loglike,model1$max.loglike,0)
    }
    
    if(model2$type=="4c"){
      if(model1$max.loglike>model2$max.loglike){p=NA}
      else{p=my.chi(model1$max.loglike,model2$max.loglike,4)}
    }
    
  }
  
  if(model1$type=="4c"){
    if(model2$type=="binorm"){
      if(model1$max.loglike<model2$max.loglike){p=NA}
      p=my.chi(model2$max.loglike,model1$max.loglike,6)
    }
    
    if(model2$type=="2cc"){
      if(model1$max.loglike<model2$max.loglike){p=NA}
      else{p=my.chi(model2$max.loglike,model1$max.loglike,5)}
    }
    
    if(model2$type=="2cu"){
      if(model1$max.loglike<model2$max.loglike){p=NA}
      else{p=my.chi(model2$max.loglike,model1$max.loglike,4)}
    }
    
    if(model2$type=="4c"){
      p=my.chi(model1$max.loglike,model2$max.loglike,0)
    }
  }
  
  if(is.na(p)){warning("The smaller model has a larger maximum likelihood value than the larger model.")}
  names(p)="p value"
  
  return(p)
}


ROCcoords <- function(model, direction = 'auto', x, input){
  directions <- c("auto","<",">")
  inputs <- c("se", "sp", "t", "sens", "spec", "thr", "sensitivity", "specificity", "threshold")
  if(!(direction %in% directions)){stop("direction not supported")}
  if(!(input %in% inputs)){stop("input not supported")}
  if(input %in% c("se","sens","sensitivity")){if(x>1 || x<0){stop("sensitivity must be between 0 and 1")}}
  if(input %in% c("sp","spec","specificity")){if(x>1 || x<0){stop("specificity must be between 0 and 1")}}

  lambda <- model$lambda
  ctrl<-boxcox(model$control,lambda)
  cs<-boxcox(model$case,lambda)
  
  if(direction=="auto")
  {
    if(median(cs)>=median(ctrl)){
      direction="<"
    }
    else{
      direction=">"
    }
  }
  
  if(model$type == "binorm"){
    mn1=model$mu.controls
    ms1=model$mu.cases   
    sn1=model$sig.controls
    ss1=model$sig.cases
    pn1=model$pi.controls
    ps1=model$pi.cases
    
    if(input %in% c("se","sens","sensitivity")){
      if(direction=="<"){
        t <- ms1-ss1*qnorm(x/ps1)
        FPF <- pn1*pnorm((mn1-t)/sn1)
        t <- boxcox.inv(t,lambda)
        out <- c(t, 1-FPF, x)
      }
      if(direction==">"){
        x <- 1-x
        t <- ms1-ss1*qnorm(x/ps1)
        FPF <- pn1*pnorm((mn1-t)/sn1)
        t <- boxcox.inv(t,lambda)
        out <- c(t, FPF, 1-x)     
      }
    }
    
    if(input %in% c("sp","spec","specificity")){
      if(direction=="<"){
        x <- 1-x
        t <- mn1-sn1*qnorm(x/pn1)
        TPF <- ps1*pnorm((ms1-t)/ss1)
        t <- boxcox.inv(t,lambda)
        out <- c(t, 1-x, TPF)        
      }
       if(direction==">"){
         t <- mn1-sn1*qnorm(x/pn1)
         TPF <- ps1*pnorm((ms1-t)/ss1)
         t <- boxcox.inv(t,lambda)
         out <- c(t, x, 1-TPF) 
       }
    }
    
    if(input %in% c("t","thr","threshold")){
      t <- boxcox(x,lambda)
      FPF <- pn1*pnorm((mn1-t)/sn1)
      TPF <- ps1*pnorm((ms1-t)/ss1)
      
      if(direction=="<"){
        out <- c(x,1-FPF,TPF)
      }
      
      if(direction==">"){
        out <- c(x,FPF,1-TPF)
      }
    }
  }
  
  if(model$type == "2cc"){
    mn1=model$mu.controls
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sn1=model$sig.controls
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pn1=model$pi.controls
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    if(input %in% c("t","thr","threshold")){
      t <- boxcox(x,lambda)
      FPF <- pn1*pnorm((mn1-t)/sn1)
      TPF <- ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
      
      if(direction=="<"){
        out <- c(x,1-FPF,TPF)
      }
      
      if(direction==">"){
        out <- c(x,FPF,1-TPF)
      }
    }
    
    if(input %in% c("sp","spec","specificity")){
      if(direction=="<"){
        x <- 1-x
        t <- mn1-sn1*qnorm(x/pn1)
        TPF <- ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
        t <- boxcox.inv(t,lambda)
        out <- c(t, 1-x, TPF)        
      }
      if(direction==">"){
        t <- mn1-sn1*qnorm(x/pn1)
        TPF <- ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
        t <- boxcox.inv(t,lambda)
        out <- c(t, x, 1-TPF) 
      }
    }
    
    if(input %in% c("se","sens","sensitivity")){
      tmin=min(ms1-5*ss1,ms2-5*ss2,mn1-5*sn1)
      tmax=max(ms1+5*ss1,ms2+5*ss2,mn1+5*sn1)

      coord.help1 <- function(t,sens,ps1,ms1,ss1,ps2,ms2,ss2){x=ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
      out<- (sens-x)^2}
      
      if(direction=="<"){
        t.opt <- optimize(f=coord.help1,interval=c(tmin,tmax),sens=x,ps1,ms1,ss1,ps2,ms2,ss2)
        t <- t.opt$minimum
        FPF <- (pn1*pnorm((mn1-t)/sn1))
        t <- boxcox.inv(t,lambda)
        out <- c(t, 1-FPF, x)
      }
      if(direction==">"){
        x <- 1-x
        t.opt <- optimize(f=coord.help1,interval=c(tmin,tmax),sens=x,ps1,ms1,ss1,ps2,ms2,ss2)
        t <- t.opt$minimum
        FPF <- (pn1*pnorm((mn1-t)/sn1))
        t <- boxcox.inv(t,lambda)
        out <- c(t, FPF, 1-x)     
      }
    }
  }
  
  if(model$type == "2cu" || model$type=="4c"){
    mu.n=model$mu.controls
    mn1=mu.n[1]
    mn2=mu.n[2]
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sig.n=model$sig.controls
    sn1=sig.n[1]
    sn2=sig.n[2]
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pi.n=model$pi.controls
    pn1=pi.n[1]
    pn2=pi.n[2]
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    tmin=min(ms1-5*ss1,ms2-5*ss2,mn1-5*sn1,mn2-5*sn2)
    tmax=max(ms1+5*ss1,ms2+5*ss2,mn1+5*sn1,mn2+5*sn2)

    if(input %in% c("se","sens","sensitivity")){
      coord.help2 <- function(t,sens,ps1,ms1,ss1,ps2,ms2,ss2){x=ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
      out<- (sens-x)^2}
      
      if(direction=="<"){
        t.opt <- optimize(f=coord.help2,interval=c(tmin,tmax),sens=x,ps1,ms1,ss1,ps2,ms2,ss2)
        t <- t.opt$minimum
        FPF <- (pn1*pnorm((mn1-t)/sn1)+pn2*pnorm((mn2-t)/sn2))
        t <- boxcox.inv(t,lambda)
        out <- c(t, 1-FPF, x)
      }
      if(direction==">"){
        x <- 1-x
        t.opt <- optimize(f=coord.help2,interval=c(tmin,tmax),sens=x,ps1,ms1,ss1,ps2,ms2,ss2)
        t <- t.opt$minimum
        FPF <- (pn1*pnorm((mn1-t)/sn1)+pn2*pnorm((mn2-t)/sn2))
        t <- boxcox.inv(t,lambda)
        out <- c(t, FPF, 1-x)     
      }
    }
    
    if(input %in% c("sp","spec","specificity")){
      
      coord.help3 <- function(t,spec,pn1,mn1,sn1,pn2,mn2,sn2){x=pn1*pnorm((mn1-t)/sn1)+pn2*pnorm((mn2-t)/sn2)
      out<- (spec-x)^2}
      
      if(direction=="<"){
        x <- 1-x
        t.opt <- optimize(f=coord.help3,interval=c(tmin,tmax),spec=x,pn1,mn1,sn1,pn2,mn2,sn2)
        t <- t.opt$minimum
        TPF <- ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
        t <- boxcox.inv(t,lambda)
        out <- c(t, 1-x, TPF)        
      }
      if(direction==">"){
        t.opt <- optimize(f=coord.help3,interval=c(tmin,tmax),spec=x,pn1,mn1,sn1,pn2,mn2,sn2)
        t <- t.opt$minimum
        TPF <- ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
        t <- boxcox.inv(t,lambda)
        out <- c(t, x, 1-TPF) 
      }
    }
    
    if(input %in% c("t","thr","threshold")){
      t <- boxcox(x,lambda)
      FPF <- pn1*pnorm((mn1-t)/sn1)+pn2*pnorm((mn2-t)/sn2)
      TPF <- ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
      
      if(direction=="<"){
        out <- c(x,1-FPF,TPF)
      }
      
      if(direction==">"){
        out <- c(x,FPF,1-TPF)
      }
    }
  }  
  
  names(out) <- c("threshold","specificity", "sensitivity")
  return(out)
}


#### AUC method
ROCauc <- function(model,direction="auto"){
  
  directions <- c("auto","<",">")
  if(!(direction %in% directions)){stop("direction not supported")}
  ctrl=model$control
  cs=model$case
  emp.AUC=pROC::auc(roc(controls=ctrl,cases=cs,direction=direction))
  if(direction=="auto")
  {
    if(median(cs)>=median(ctrl)){
      direction="<"
    }
    else{
      direction=">"
    }
  }
  
  if(model$type == "binorm"){
    mn1=model$mu.controls
    ms1=model$mu.cases   
    sn1=model$sig.controls
    ss1=model$sig.cases
    pn1=model$pi.controls
    ps1=model$pi.cases
    if(direction==">"){
      smooth.AUC=ps1*pn1*pnorm((mn1-ms1)/sqrt(ss1^2+sn1^2)) 
    }
    if(direction=="<"){
      smooth.AUC=ps1*pn1*pnorm((ms1-mn1)/sqrt(ss1^2+sn1^2)) 
    }    
  }
  
  if(model$type == "2cc"){
    mn1=model$mu.controls
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sn1=model$sig.controls
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pn1=model$pi.controls
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    if(direction==">"){
      smooth.AUC=ps1*pn1*pnorm((mn1-ms1)/sqrt(ss1^2+sn1^2)) + ps2*pn1*pnorm((mn1-ms2)/sqrt(ss2^2+sn1^2)) 
    }
    
    if(direction=="<"){
      smooth.AUC=ps1*pn1*pnorm((ms1-mn1)/sqrt(ss1^2+sn1^2)) + ps2*pn1*pnorm((ms2-mn1)/sqrt(ss2^2+sn1^2)) 
    }
  }
  
  if(model$type == "2cu" || model$type=="4c"){
    mu.n=model$mu.controls
    mn1=mu.n[1]
    mn2=mu.n[2]
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sig.n=model$sig.controls
    sn1=sig.n[1]
    sn2=sig.n[2]
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pi.n=model$pi.controls
    pn1=pi.n[1]
    pn2=pi.n[2]
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    if(direction==">"){
      smooth.AUC=ps1*pn1*pnorm((mn1-ms1)/sqrt(ss1^2+sn1^2)) + ps1*pn2*pnorm((mn2-ms1)/sqrt(ss1^2+sn2^2)) + ps2*pn1*pnorm((mn1-ms2)/sqrt(ss2^2+sn1^2)) + ps2*pn2*pnorm((mn2-ms2)/sqrt(ss2^2+sn2^2))
    }
    
    if(direction=="<"){
      smooth.AUC=ps1*pn1*pnorm((ms1-mn1)/sqrt(ss1^2+sn1^2)) + ps1*pn2*pnorm((ms1-mn2)/sqrt(ss1^2+sn2^2)) + ps2*pn1*pnorm((ms2-mn1)/sqrt(ss2^2+sn1^2)) + ps2*pn2*pnorm((ms2-mn2)/sqrt(ss2^2+sn2^2))
    }
  }
  out=c(emp.AUC,smooth.AUC)
  names(out)=c("empirical AUC","smooth AUC")
  out
}

#### ROC pAUC method
ROCpauc <- function(model,spec.lower=.95,spec.upper=1,direction="auto"){
  
  if(spec.upper>1 || spec.upper < 0 || spec.lower > 1 || spec.lower < 0){stop("specificity values must be between 0 and 1")}
  directions <- c("auto","<",">")
  if(!(direction %in% directions)){stop("direction not supported")}
  ctrl=model$control
  cs=model$case
  emp.pauc=pROC::auc(roc(controls=ctrl,cases=cs,direction=direction),partial.auc=c(spec.lower,spec.upper))
  if(direction=="auto")
  {
    if(median(cs)>=median(ctrl)){
      direction="<"
    }
    else{
      direction=">"
    }
  }
  
  if(model$type == "binorm"){
    mn1=model$mu.controls
    ms1=model$mu.cases   
    sn1=model$sig.controls
    ss1=model$sig.cases
    pn1=model$pi.controls
    ps1=model$pi.cases
    
    tmin=min(ms1-5*ss1,mn1-5*sn1)
    tmax=max(ms1+5*ss1,mn1+5*sn1)
    t=seq(tmin,tmax,length=10000)
    tn=length(t)
    
    FPF=(pn1*pnorm((mn1-t)/sn1))
    TPF=ps1*pnorm((ms1-t)/ss1)
  }
  
  if(model$type == "2cc"){
    mn1=model$mu.controls
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sn1=model$sig.controls
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pn1=model$pi.controls
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    tmin=min(ms1-5*ss1,ms2-5*ss2,mn1-5*sn1)
    tmax=max(ms1+5*ss1,ms2+5*ss2,mn1+5*sn1)
    t=seq(tmin,tmax,length=10000)
    tn=length(t)
    
    FPF=(pn1*pnorm((mn1-t)/sn1))
    TPF=ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
  }
  
  if(model$type == "2cu" || model$type=="4c"){
    mu.n=model$mu.controls
    mn1=mu.n[1]
    mn2=mu.n[2]
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sig.n=model$sig.controls
    sn1=sig.n[1]
    sn2=sig.n[2]
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pi.n=model$pi.controls
    pn1=pi.n[1]
    pn2=pi.n[2]
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    tmin=min(ms1-5*ss1,ms2-5*ss2,mn1-5*sn1,mn2-5*sn2)
    tmax=max(ms1+5*ss1,ms2+5*ss2,mn1+5*sn1,mn2+5*sn2)
    t=seq(tmin,tmax,length=10000)
    tn=length(t)
    
    FPF=(pn1*pnorm((mn1-t)/sn1)+pn2*pnorm((mn2-t)/sn2))
    TPF=ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
  }
  
  FPF[tn]=0
  TPF[tn]=0
  FPF[1]=1
  TPF[1]=1
  
  if(direction==">"){
    FPF = 1-FPF
    TPF = 1-TPF
  }
  
  FPF.upper=1-spec.lower
  FPF.lower=1-spec.upper
  parts=which(FPF>=FPF.lower & FPF<=FPF.upper)
  trap.n=length(parts)+1
  trap.end=min(parts)
  trap.start=max(parts)
  pAUC.trap.calc=rep(NA,trap.n)
  
  TPF.lower=approx(FPF,TPF,xout=FPF.lower,rule=2)$y
  TPF.upper=approx(FPF,TPF,xout=FPF.upper,rule=2)$y
  pAUC.trap.calc[1]=abs(FPF[trap.start]-FPF.lower)*abs(TPF[trap.start]+TPF.lower)/2
  pAUC.trap.calc[trap.n]=abs(FPF.upper-FPF[trap.end])*abs(TPF[trap.end]+TPF.upper)/2
  if(direction==">"){
    pAUC.trap.calc[1]=abs(FPF[trap.start]-FPF.upper)*abs(TPF[trap.start]+TPF.lower)/2
    pAUC.trap.calc[trap.n]=abs(FPF.lower-FPF[trap.end])*abs(TPF[trap.end]+TPF.upper)/2
  }
  
  if(trap.n>2)
  {
    for(i in 2:(trap.n-1))
    {
      pAUC.trap.calc[i]=(abs(FPF[trap.start-i+2]-FPF[trap.start-i+1])*abs(TPF[trap.start-i+2]+TPF[trap.start-i+1])/2)
    }
  }
  
  pAUC.trap.calc=sum(pAUC.trap.calc)
  
  out=c(emp.pauc,pAUC.trap.calc, spec.upper, spec.lower)
  names(out)=c("empirical pAUC","smooth pAUC","upper specificity", "lower specificity")
  out
}

#### ROC plot method
ROCplot <- function(model,direction="auto"){
  directions <- c("auto","<",">")
  if(!(direction %in% directions)){stop("direction not supported")}
  ctrl=boxcox(model$control,model$lambda)
  cs=boxcox(model$case,model$lambda)
  proc=pROC::roc(controls=ctrl,cases=cs,direction=direction)
  plot(proc,legacy.axes=T)
  if(direction=="auto")
  {
    if(median(cs)>=median(ctrl)){
      direction="<"
    }
    else{
      direction=">"
    }
  }
  
  if(model$type == "binorm"){
    mn1=model$mu.controls
    ms1=model$mu.cases   
    sn1=model$sig.controls
    ss1=model$sig.cases
    pn1=model$pi.controls
    ps1=model$pi.cases
    
    tmin=min(ms1-5*ss1,mn1-5*sn1)
    tmax=max(ms1+5*ss1,mn1+5*sn1)
    t=seq(tmin,tmax,length=10000)
    tn=length(t)
    
    FPF<-rep(NA,tn)
    TPF<-rep(NA,tn)
    
    FPF<-pn1*pnorm((mn1-t)/sn1)
    TPF<-ps1*pnorm((ms1-t)/ss1)
    
    if(direction=="<"){
      lines(1-FPF,TPF,type="l",lty=5,lwd=2,col="blue")        
    }
    if(direction==">"){
      lines(FPF,1-TPF,type="l",lty=5,lwd=2,col="blue")         
    }
  }
  
  if(model$type == "2cc"){
    mn1=model$mu.controls
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sn1=model$sig.controls
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pn1=model$pi.controls
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    tmin=min(ms1-5*ss1,ms2-5*ss2,mn1-5*sn1)
    tmax=max(ms1+5*ss1,ms2+5*ss2,mn1+5*sn1)
    t=seq(tmin,tmax,length=10000)
    tn=length(t)
    
    FPF=rep(NA,tn)
    TPF<-rep(NA,tn)
    
    FPF=(pn1*pnorm((mn1-t)/sn1))
    TPF=ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
    
    if(direction=="<"){
      lines(1-FPF,TPF,type="l",lty=5,lwd=2,col="blue")        
    }
    if(direction==">"){
      lines(FPF,1-TPF,type="l",lty=5,lwd=2,col="blue")        
    }
  }
  
  if(model$type == "2cu" || model$type=="4c"){
    mu.n=model$mu.controls
    mn1=mu.n[1]
    mn2=mu.n[2]
    
    mu.s=model$mu.cases
    ms1=mu.s[1]
    ms2=mu.s[2]
    
    sig.n=model$sig.controls
    sn1=sig.n[1]
    sn2=sig.n[2]
    
    sig.s=model$sig.cases
    ss1=sig.s[1]
    ss2=sig.s[2]
    
    pi.n=model$pi.controls
    pn1=pi.n[1]
    pn2=pi.n[2]
    
    pi.s=model$pi.cases
    ps1=pi.s[1]
    ps2=pi.s[2]
    
    tmin=min(ms1-5*ss1,ms2-5*ss2,mn1-5*sn1,mn2-5*sn2)
    tmax=max(ms1+5*ss1,ms2+5*ss2,mn1+5*sn1,mn2+5*sn2)
    t=seq(tmin,tmax,length=10000)
    tn=length(t)
    
    FPF=rep(NA,tn)
    TPF<-rep(NA,tn)
    
    FPF=(pn1*pnorm((mn1-t)/sn1)+pn2*pnorm((mn2-t)/sn2))
    TPF=ps1*pnorm((ms1-t)/ss1)+ps2*pnorm((ms2-t)/ss2)
    
    if(direction=="<"){
      lines(1-FPF,TPF,type="l",lty=5,lwd=2,col="blue")        
    }
    if(direction==">"){
      lines(FPF,1-TPF,type="l",lty=5,lwd=2,col="blue")         
    }
  }
}


rmix <- function(n,mu1,s1,mu2,s2,p1){
  if(p1>1 || p1< 0){stop("proportion should be between 0 and 1")}
  if(n<1){stop("n must be greater than 0")}
  if(s1 < 0 || s2 < 0){stop("standard deviations must be non-negative")}
  z <- rbinom(n,size=1,prob=p1)
  z * rnorm(n,mean=mu1,sd=s1) + (1-z)*rnorm(n,mean=mu2,sd=s2)
}



em.twocomp.m1 <- function(x.all,case.indicator,max.iters=1000,errtol=1e-9,control.comp=1,start.vals=NULL){
  #### em.twocomp.m1
  #### EM algorithm for two-component mixture
  #### in which all controls are in component control.comp and cases are a mixture of the two components
  #### components are ordered by their means
  
  k <- 2   ### two components
  
  n.cases <- sum(case.indicator)
  n.all <- length(case.indicator)
  n.controls <- n.all - n.cases
  tau <- matrix(0,nrow=n.all,ncol=k)
  ep.sig <- 1e-5
  
  ### create starting values for mu, sigma and pi

    comp.sig <- start.vals$sig
    comp.pi <- start.vals$pi
    comp.mu <- start.vals$mu

  
  iter <- 0
  current.error <- 1
  max.loglike <- -1e6
  while (iter < max.iters && current.error > errtol){
    iter <- iter + 1
    last.ll <- max.loglike
    
    ### E step
    for(i in 1:n.all){
      for(j in 1:k){
        if (i <= n.cases){
          tau[i,j] = comp.pi[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
        } else {
          tau[i,j] <- 1 * (j==control.comp)  ## each control has prob 1 of being in component component.map
        }
      }
      sum.tau <- sum(tau[i,])
      if (sum.tau > 0){
        tau[i,] <- tau[i,] / sum.tau
      } else {
        tau[i,] <- comp.pi
      }
    }
    
    ### M step
    for(j in 1:k){
      comp.mu[j] <- sum(tau[,j]*x.all)/sum(tau[,j])
      if (j > 1){
        if (comp.mu[j] < comp.mu[(j-1)]){  #### enforce constraint on means
          comp.mu[j] <- comp.mu[(j-1)] + 0.01*comp.sig[(j-1)]
        }
      }
      comp.sig[j] <- sqrt(sum(tau[,j]*(x.all-comp.mu[j])^2)/sum(tau[,j]))
    }
    comp.sig[comp.sig < ep.sig] <- ep.sig
    comp.pi <- apply(tau[1:n.cases,],2,mean)
    
    if(sum(!is.finite(comp.sig))||sum(!is.finite(comp.mu))||sum(!is.finite(comp.pi))){stop("Failed to converge. Reject starting values.")}
    
    
    tmp.ll <- rep(NA,n.all)
    for(i in 1:n.cases){
      tmp.ll[i] <- 0
      for(j in 1:k){
        tmp.ll[i] = tmp.ll[i] + comp.pi[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
      }
    }
    for(i in (n.cases+1):n.all){
      tmp.ll[i] = dnorm(x.all[i],mean=comp.mu[control.comp],sd=comp.sig[control.comp])
    }
    max.loglike <- sum(log(tmp.ll))
    if (is.na(max.loglike) || is.nan(max.loglike)){
      stop("loglike is na")
    }
    current.error <- mean((last.ll - max.loglike)^2)
  }
  
  list(max.loglike=max.loglike,mu=comp.mu,sig=comp.sig,pi=comp.pi,n.iters=iter,control.comp=control.comp)
}


em.twocomp.m2 <- function(x.all,max.iters=1000,errtol=1e-9,start.vals=NULL){
  #### em.twocomp.m2
  #### EM algorithm for two-component mixture
  #### used twice for four-component model, with unconstrained components
  #### components are ordered by their means
  
  k <- 2   ### two components
  
  n.all <- length(x.all)
  tau <- matrix(0,nrow=n.all,ncol=k)
  ep.sig <- 1e-5
  
  ### create starting values for mu, sigma and pi
  comp.sig <- start.vals$sig
  comp.pi <- start.vals$pi
  comp.mu <- start.vals$mu

  iter <- 0
  current.error <- 1
  max.loglike <- -1e6
  while (iter < max.iters && current.error > errtol){
    iter <- iter + 1
    last.ll <- max.loglike
    
    ### E step
    for(i in 1:n.all){
      for(j in 1:k){
        tau[i,j] = comp.pi[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
      }
      sum.tau <- sum(tau[i,])
      if (sum.tau > 0){
        tau[i,] <- tau[i,] / sum.tau
      } else {
        tau[i,] <- comp.pi
      }
    }
    
    ### M step
    for(j in 1:k){
      comp.mu[j] <- sum(tau[,j]*x.all)/sum(tau[,j])
      if (j > 1){
        if (comp.mu[j] < comp.mu[(j-1)]){  #### enforce constraint on means
          comp.mu[j] <- comp.mu[(j-1)] + 0.01*comp.sig[(j-1)]
        }
      }
      comp.sig[j] <- sqrt(sum(tau[,j]*(x.all-comp.mu[j])^2)/sum(tau[,j]))
    }
    comp.sig[comp.sig < ep.sig] <- ep.sig
    comp.pi <- apply(tau,2,mean)
    
    if(sum(!is.finite(comp.sig))||sum(!is.finite(comp.mu))||sum(!is.finite(comp.pi))){stop("Failed to converge. Reject starting values.")}
    
    
    tmp.ll <- rep(0,n.all)
    for(i in 1:n.all){
      for(j in 1:k){
        tmp.ll[i] = tmp.ll[i] + comp.pi[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
      }
    }
    max.loglike <- sum(log(tmp.ll))
    if (is.na(max.loglike) || is.nan(max.loglike)){
      stop("loglike is na")
    }
    current.error <- mean((last.ll - max.loglike)^2)
  }
  
  list(max.loglike=max.loglike,mu=comp.mu,sig=comp.sig,pi=comp.pi,n.iters=iter)
}



em.twocomp.m3 <- function(x.all,case.indicator,max.iters=1000,errtol=1e-9,control.comp=1,start.vals=NULL){
  #### em.twocomp.m3
  #### EM algorithm for two-component mixture
  #### in which controls and cases are mixtures of two components
  #### components are ordered by their means
  
  k <- 2   ### two components
  
  n.cases <- sum(case.indicator)
  n.all <- length(case.indicator)
  n.controls <- n.all - n.cases
  ep.sig <- 1e-5
  tau <- matrix(0,nrow=n.all,ncol=k)
  
  ### create starting values for mu, sigma and pi

    comp.sig <- start.vals$sig
    comp.pi.cs <- start.vals$pi.cs
    comp.pi.ctrl <- start.vals$pi.ctrl
    comp.mu <- start.vals$mu

  
  iter <- 0
  current.error <- 1
  max.loglike <- -1e6
  while (iter < max.iters && current.error > errtol){
    iter <- iter + 1
    last.ll <- max.loglike
    
    ### E step
    for(i in 1:n.all){
      for(j in 1:k){
        if (i <= n.cases){
          tau[i,j] = comp.pi.cs[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
        } else {
          # edit
          tau[i,j] <- comp.pi.ctrl[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
        }
      }
      sum.tau <- sum(tau[i,])
      if (sum.tau > 0){
        tau[i,] <- tau[i,] / sum.tau
      } else {
        if(i <= n.cases){
          tau[i,] <- comp.pi.cs
        } else{
          tau[i,] <-comp.pi.ctrl
        }
      }
    }
    
    ### M step
    for(j in 1:k){
      comp.mu[j] <- sum(tau[,j]*x.all)/sum(tau[,j])
      if (j > 1){
        if (comp.mu[j] < comp.mu[(j-1)]){  #### enforce constraint on means
          comp.mu[j] <- comp.mu[(j-1)] + 0.01*comp.sig[(j-1)]
        }
      }
      comp.sig[j] <- sqrt(sum(tau[,j]*(x.all-comp.mu[j])^2)/sum(tau[,j]))
    }
    comp.sig[comp.sig < ep.sig] <- ep.sig
    comp.pi.cs <- apply(tau[1:n.cases,],2,mean)
    comp.pi.ctrl <- apply(tau[(1+n.cases):n.all,],2,mean)
    
    if(sum(!is.finite(comp.sig))||sum(!is.finite(comp.mu))||sum(!is.finite(comp.pi.cs))||sum(!is.finite(comp.pi.ctrl))){stop("Failed to converge. Reject starting values.")}
    
    tmp.ll <- rep(NA,n.all)
    for(i in 1:n.cases){
      tmp.ll[i] <- 0
      for(j in 1:k){
        tmp.ll[i] = tmp.ll[i] + comp.pi.cs[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
      }
    }
    for(i in (n.cases+1):n.all){
      tmp.ll[i] <- 0
      for(j in 1:k){
        tmp.ll[i] = tmp.ll[i] + comp.pi.ctrl[j]*dnorm(x.all[i],mean=comp.mu[j],sd=comp.sig[j])
      }    
    }
    max.loglike <- sum(log(tmp.ll))
    if (is.na(max.loglike) || is.nan(max.loglike)){
      stop("loglike is na")
    }
    current.error <- mean((last.ll - max.loglike)^2)
  }
  
  list(max.loglike=max.loglike,mu=comp.mu,sig=comp.sig,pi.cs=comp.pi.cs,pi.ctrl=comp.pi.ctrl,n.iters=iter,control.comp=control.comp)
}


boxcox <- function(x,lambda){
  ### boxcox transformation of x
  ### x is assumed to be positive
  ### if x is non-positive, should return NA
  #### NOTE: This is always a monotonically increasing function
  
  ep.lambda <- 1e-10
  y <- 0
  if (abs(lambda) < ep.lambda){
    y <- log(x)
  } else {
    y <- (x^lambda - 1)/lambda
  }
  y
}


boxcox.inv <- function(y,lambda){
  ### inverse boxcox transformation of y
  
  ep.lambda <- 1e-10
  x <- 0
  if (abs(lambda) < ep.lambda){
    x <- exp(y)
  } else {
    x <- (lambda*y + 1)^(1/lambda)
  }
  x
}


boxcox.deriv <- function(x,lambda){
  ep.lambda <- 1e-10
  if(abs(lambda)<ep.lambda){
    y <- 1/x
  } else{
    y <- x^(lambda-1)
  }
  y
}

boxcox.inv.density <- function(y,lambda,mu,sig){
  ep.lambda <- 1e-10
  x <- 0
  if (abs(lambda) < ep.lambda){
    x <- dnorm(log(y),mu,sig)*abs(1/y)
  }
  else{
    x <-dnorm((y^lambda-1)/lambda,mu,sig)*abs(y^(lambda-1))
  }
  x
}


bc.binorm <- function(case,control,lambda.bounds=c(-5,5))
{
  t0 <- proc.time()[3]
  if (min(c(case,control)) <= 0){
    stop("All input values must be positive.")
  }
  if (sum(is.finite(c(case,control))) != length(is.finite(c(case,control)))){
    stop("All input values must be finite.")
  }
  
  binorm <- function(lambda,case,control)
  {
    transform.cs=boxcox(case,lambda)
    transform.ctrl=boxcox(control,lambda)
    
    n <- length(case)
    mu.cs <- mean(transform.cs)
    mu.ctrl <- mean(transform.ctrl)
    sig.cs <- sd(transform.cs)
    sig.ctrl <- sd(transform.ctrl)
    
    ll.case <- -(n/2)*log(2*pi*sig.cs^2) - sum((transform.cs-mu.cs)^2)/(2*sig.cs^2) + (lambda - 1)*sum(log(case))
    ll.control <- -(n/2)*log(2*pi*sig.ctrl^2) - sum((transform.ctrl-mu.ctrl)^2)/(2*sig.ctrl^2) + (lambda - 1)*sum(log(control))
    ll=ll.case+ll.control
    
    res <- list()
    res$mu.cases <- mu.cs
    res$sig.cases <- sig.cs
    res$mu.controls <- mu.ctrl
    res$sig.controls <- sig.ctrl
    res$max.loglike <- ll
    res
  }
  
  help.bc.binorm <- function(lambda,case,control)
  {
    out <- binorm(lambda,case,control)
    ll <- out$max.loglike
    ll
  }
  
  bc.binorm.optim <- optimize(help.bc.binorm, case=case, control=control, lower=lambda.bounds[1], upper=lambda.bounds[2],maximum=T)
  res <- list()
  
  res$lambda <-bc.binorm.optim$maximum
  res$type <- "binorm"
  model <- binorm(res$lambda, case, control)
  res$mu.cases <- model$mu.cases
  res$sig.cases <- model$sig.cases
  res$pi.cases=1
  res$mu.controls <- model$mu.controls
  res$sig.controls <- model$sig.controls
  res$pi.controls=1
  res$max.loglike=bc.binorm.optim$objective
  res$case=case
  res$control=control
  
  cs.estimates <- rnorm(1e6,res$mu.cases,res$sig.cases)
  ctrl.estimates <- rnorm(1e6,res$mu.controls,res$sig.controls)
  
  if(res$lambda>0){
    cs.estimates[which(cs.estimates<=-1/res$lambda)]<-NA
    ctrl.estimates[which(ctrl.estimates<=-1/res$lambda)]<-NA
  }
  
  if(res$lambda<0){
    cs.estimates[which(cs.estimates>=-1/res$lambda)]<-NA
    ctrl.estimates[which(ctrl.estimates>=-1/res$lambda)]<-NA
  }
  
  monte.carlo.cases <- boxcox.inv(cs.estimates,res$lambda)
  monte.carlo.controls <- boxcox.inv(ctrl.estimates,res$lambda)
  
  res$mu.cases.unt <- mean(monte.carlo.cases,na.rm=T)
  res$mu.controls.unt <- mean(monte.carlo.controls,na.rm=T)
  res$sig.cases.unt <- sd(monte.carlo.cases,na.rm=T)
  res$sig.controls.unt <- sd(monte.carlo.controls,na.rm=T)
  res$time <- proc.time()[3] - t0
  res
}


bc.twocomp <- function(x.cases,x.controls,constrained=T,lambda.bounds=c(-5,5),control.comp=1,start.vals=NULL){
  ### control.comp dictates which component controls should be assigned to
  k <- 2
  t0 <- proc.time()[3]
  
  #### Helper function for boxcox two component model
  #### in which all controls are in one component (control.comp) and cases are a mixture of the two components
  
  x.all <- c(x.cases,x.controls)
  #### important that the cases are place first (assumed in em.twocomp.m1)!!!
  
  if (min(x.all) <= 0){
    stop("All input values must be positive.")
  }
  if (sum(is.finite(c(x.cases,x.controls))) != length(is.finite(c(x.cases,x.controls)))){
    stop("All input values must be finite.")
  }
  
  n.cases <- length(x.cases)
  n.controls <- length(x.controls)
  case.indicator <- c(rep(1,n.cases),rep(0,n.controls))
  n.all <- n.cases + n.controls
  ep.sig <- 1e-5
  
  if(is.null(start.vals)){
    clusters <- kmeans(x.all,k)
    comp.mu <-as.vector(as.numeric(clusters$centers))
    comp.sig <- c(sd(x.all[which(clusters$cluster==1)]),sd(x.all[which(clusters$cluster==2)]))
    comp.sig[comp.sig < ep.sig] <- ep.sig
    
    if(constrained){
      comp.pi <- c(length(x.all[which(clusters$cluster==1 & case.indicator==1)])/n.cases,
                     length(x.all[which(clusters$cluster==2 & case.indicator==1)])/n.cases)
      start.vals <- list(mu=comp.mu,sig=comp.sig,pi=comp.pi)
    } 
    
    if(!constrained){
      comp.pi.cs <-c(length(x.all[which(clusters$cluster==1 & case.indicator==1)])/n.cases,
                     length(x.all[which(clusters$cluster==2 & case.indicator==1)])/n.cases)
      comp.pi.ctrl <-c(length(x.all[which(clusters$cluster==1 & case.indicator==0)])/n.controls,
                       length(x.all[which(clusters$cluster==2 & case.indicator==0)])/n.controls)
      start.vals <- list(mu=comp.mu,sig=comp.sig,pi.cs=comp.pi.cs,pi.ctrl=comp.pi.ctrl)
    }
    
    start.vals$sig[is.na(start.vals$sig)]<-ep.sig
  }
  
  if(constrained){
    if(start.vals$mu[1]>start.vals$mu[2]){
      start.vals$mu <- rev(start.vals$mu)
      start.vals$sig <- rev(start.vals$sig)
      start.vals$pi <- rev(start.vals$pi)
    }
  }
    
    if(!constrained){
      if(start.vals$mu[1]>start.vals$mu[2]){
        start.vals$mu <- rev(start.vals$mu)
        start.vals$sig <- rev(start.vals$sig)
        start.vals$pi.cs <- rev(start.vals$pi.cs)
        start.vals$pi.ctrl <- rev(start.vals$pi.ctrl)
      }
    }
  

  help.bc <- function(lambda,x.all,constrained,case.indicator,control.comp=1,start.vals){
    k <- 2  ### two components
    y.all <- boxcox(x.all,lambda)
    
    var <- (start.vals$sig)^2*(boxcox.deriv(start.vals$mu,lambda))^2
    start.vals$sig <- sqrt(var)
    start.vals$mu <- boxcox(start.vals$mu,lambda)

    if(constrained){
      em.output <- em.twocomp.m1(y.all,case.indicator=case.indicator,control.comp=control.comp,start.vals=start.vals)
    }
    if(!constrained){
      em.output <- em.twocomp.m3(y.all,case.indicator=case.indicator,control.comp=control.comp,start.vals=start.vals)
    }
    sum.log <- sum(log(x.all))
    plus <- (lambda - 1)*sum.log   ### adjustment to log likelihood due to BC transformation
    -1*em.output$max.loglike - plus
  }
  
  help.output <- optimize(help.bc,x.all=x.all,constrained=constrained,case.indicator=case.indicator,
                          lower=lambda.bounds[1],upper=lambda.bounds[2],control.comp=control.comp,start.vals=start.vals)
  
  res <- list()
  res$lambda <- help.output$minimum
  y.all <- boxcox(x.all,res$lambda)
  
  
  var <- (start.vals$sig)^2*(boxcox.deriv(start.vals$mu,res$lambda))^2
  start.vals$sig <- sqrt(var)
  start.vals$mu <- boxcox(start.vals$mu,res$lambda)

  
  if(constrained){
    res$type <- "2cc"
    em.output <- em.twocomp.m1(y.all,case.indicator=case.indicator,control.comp=control.comp,start.vals=start.vals)   
  }
  if(!constrained){
    res$type <- "2cu"
    em.output <- em.twocomp.m3(y.all,case.indicator=case.indicator,control.comp=control.comp,start.vals=start.vals)   
  }
  
  res$mu.cases <- em.output$mu
  res$sig.cases <- em.output$sig
  
  if(constrained){
    res$pi.cases <- em.output$pi
    res$mu.controls<- em.output$mu[control.comp]
    res$sig.controls <- em.output$sig[control.comp]
    res$pi.controls <- 1
    
    
    n=which(res$mu.cases!=res$mu.controls)	
    cs.estimates <- rnorm(1e6,res$mu.cases[n],res$sig.cases[n])
    ctrl.estimates <- rnorm(1e6,res$mu.controls,res$sig.controls)
    
    
    if(res$lambda>0){
      cs.estimates[which(cs.estimates<=-1/res$lambda)]<-NA
      ctrl.estimates[which(ctrl.estimates<=-1/res$lambda)]<-NA
    }
    
    if(res$lambda<0){
      cs.estimates[which(cs.estimates>=-1/res$lambda)]<-NA
      ctrl.estimates[which(ctrl.estimates>=-1/res$lambda)]<-NA
    }
    
    monte.carlo.controls <- boxcox.inv(ctrl.estimates,res$lambda)
    res$mu.controls.unt <- mean(monte.carlo.controls,na.rm=T)
    res$sig.controls.unt <- sd(monte.carlo.controls,na.rm=T)
    
    monte.carlo.cases <- boxcox.inv(cs.estimates,res$lambda)
    res$mu.cases.unt <- c(NA,NA)
    res$sig.cases.unt <- c(NA,NA)
    res$mu.cases.unt[n] <- mean(monte.carlo.cases,na.rm=T)
    res$mu.cases.unt[-n] <- res$mu.controls.unt
    res$sig.cases.unt[n] <- sd(monte.carlo.cases,na.rm=T)
    res$sig.cases.unt[-n] <- res$sig.controls.unt
  }
  
  if(!constrained){
    res$pi.cases <- em.output$pi.cs
    res$mu.controls <- em.output$mu
    res$sig.controls <- em.output$sig
    res$pi.controls <- em.output$pi.ctrl
    
    ctrl1.estimates <- rnorm(1e6,res$mu.controls[1],res$sig.controls[1])
    ctrl2.estimates <- rnorm(1e6,res$mu.controls[2],res$sig.controls[2])
    
    if(res$lambda>0){
      ctrl1.estimates[which(ctrl1.estimates<=-1/res$lambda)]<-NA
      ctrl2.estimates[which(ctrl2.estimates<=-1/res$lambda)]<-NA
    }
    
    if(res$lambda<0){
      ctrl1.estimates[which(ctrl1.estimates>=-1/res$lambda)]<-NA
      ctrl2.estimates[which(ctrl2.estimates>=-1/res$lambda)]<-NA
    }
    
    monte.carlo.controls1 <- boxcox.inv(ctrl1.estimates,res$lambda)
    monte.carlo.controls2 <- boxcox.inv(ctrl2.estimates,res$lambda)
    res$mu.controls.unt <- c(mean(monte.carlo.controls1,na.rm=T),mean(monte.carlo.controls2,na.rm=T))
    res$sig.controls.unt <- c(sd(monte.carlo.controls1,na.rm=T),sd(monte.carlo.controls2,na.rm=T))
    
    n=which(res$mu.cases==res$mu.controls[1])
    res$mu.cases.unt <- c(NA,NA)
    res$sig.cases.unt <- c(NA,NA)
    res$mu.cases.unt[n] <- res$mu.controls.unt[1]
    res$mu.cases.unt[-n] <- res$mu.controls.unt[-1]
    res$sig.cases.unt[n] <- res$sig.controls.unt[1]
    res$sig.cases.unt[-n] <- res$sig.controls.unt[-1]
    
  }
  
  res$max.loglike <- em.output$max.loglike + (res$lambda - 1)*sum(log(x.all))
  res$case=x.cases
  res$control=x.controls
  
  res
  
  res$time <- proc.time()[3] - t0
  res
}


bc.fourcomp <- function(x.cases,x.controls,lambda.bounds=c(-5,5),start.vals.cases=NULL,start.vals.controls=NULL){
  k <- 2
  t0 <- proc.time()[3]
  
  if (min(c(x.cases,x.controls)) <= 0){
    stop("All input values must be positive.")
  }
  if (sum(is.finite(c(x.cases,x.controls))) != length(is.finite(c(x.cases,x.controls)))){
    stop("All input values must be finite.")
  }
  
  n.cases <- length(x.cases)
  n.controls <- length(x.controls)
  ep.sig <- 1e-5
  
  if(is.null(start.vals.cases)){
    clusters.cs <- kmeans(x.cases,k)
    comp.mu.cs <-as.vector(as.numeric(clusters.cs$centers))
    comp.sig.cs <- c(sd(x.cases[which(clusters.cs$cluster==1)]),sd(x.cases[which(clusters.cs$cluster==2)]))
    comp.sig.cs[comp.sig.cs < ep.sig] <- ep.sig
    comp.pi.cs <- c(length(x.cases[which(clusters.cs$cluster==1)])/n.cases,
                 length(x.cases[which(clusters.cs$cluster==2)])/n.cases)
    start.vals.cases <- list(mu=comp.mu.cs,sig=comp.sig.cs,pi=comp.pi.cs)
    
    clusters.ctrl <- kmeans(x.controls,k)
    comp.mu.ctrl <-as.vector(as.numeric(clusters.ctrl$centers))
    comp.sig.ctrl <- c(sd(x.controls[which(clusters.ctrl$cluster==1)]),sd(x.controls[which(clusters.ctrl$cluster==2)]))
    comp.sig.ctrl[comp.sig.ctrl < ep.sig] <- ep.sig
    comp.pi.ctrl <- c(length(x.controls[which(clusters.ctrl$cluster==1)])/n.controls,
                      length(x.controls[which(clusters.ctrl$cluster==2)])/n.controls)
    start.vals.controls <- list(mu=comp.mu.ctrl,sig=comp.sig.ctrl,pi=comp.pi.ctrl)
    
    start.vals.controls$sig[is.na(start.vals.controls$sig)]<-ep.sig
    start.vals.cases$sig[is.na(start.vals.cases$sig)]<-ep.sig
  }
  
  if(start.vals.cases$mu[1]>start.vals.cases$mu[2]){
    start.vals.cases$mu <- rev(start.vals.cases$mu)
    start.vals.cases$sig <- rev(start.vals.cases$sig)
    start.vals.cases$pi <- rev(start.vals.cases$pi)
  }
  
  if(start.vals.controls$mu[1]>start.vals.controls$mu[2]){
    start.vals.controls$mu <- rev(start.vals.controls$mu)
    start.vals.controls$sig <- rev(start.vals.controls$sig)
    start.vals.controls$pi <- rev(start.vals.controls$pi)
  }

  help.bc <- function(lambda,x.cases,x.controls,start.vals.cases,start.vals.controls){
    k <- 2  ### two components
    y.cases <- boxcox(x.cases,lambda)
    y.controls <- boxcox(x.controls,lambda)
    
    var.cs <- (start.vals.cases$sig)^2*(boxcox.deriv(start.vals.cases$mu,lambda))^2
    start.vals.cases$sig <- sqrt(var.cs)
    start.vals.cases$mu <- boxcox(start.vals.cases$mu,lambda)
    var.ctrl <- (start.vals.controls$sig)^2*(boxcox.deriv(start.vals.controls$mu,lambda))^2
    start.vals.controls$sig <- sqrt(var.ctrl)
    start.vals.controls$mu <- boxcox(start.vals.controls$mu,lambda)
    
    em.cases <- em.twocomp.m2(y.cases,start.vals=start.vals.cases)
    em.controls <- em.twocomp.m2(y.controls,start.vals=start.vals.controls)
    sum.log <- sum(log(x.cases)) + sum(log(x.controls))
    plus <- (lambda - 1)*sum.log   ### adjustment to log likelihood due to BC transformation
    -1*em.cases$max.loglike - em.controls$max.loglike - plus
  }
  
  help.output <- optimize(help.bc,x.cases=x.cases,x.controls=x.controls, 
                          start.vals.cases=start.vals.cases,start.vals.controls=start.vals.controls,
                          lower=lambda.bounds[1],upper=lambda.bounds[2])
  
  res <- list()
  res$lambda <- help.output$minimum
  res$type <- "4c"
  y.cases <- boxcox(x.cases,res$lambda)
  y.controls <- boxcox(x.controls,res$lambda)
  
  var.cs <- (start.vals.cases$sig)^2*(boxcox.deriv(start.vals.cases$mu,res$lambda))^2
  start.vals.cases$sig <- sqrt(var.cs)
  start.vals.cases$mu <- boxcox(start.vals.cases$mu,res$lambda)
  var.ctrl <- (start.vals.controls$sig)^2*(boxcox.deriv(start.vals.controls$mu,res$lambda))^2
  start.vals.controls$sig <- sqrt(var.ctrl)
  start.vals.controls$mu <- boxcox(start.vals.controls$mu,res$lambda)

  em.cases <- em.twocomp.m2(y.cases,start.vals=start.vals.cases)
  em.controls <- em.twocomp.m2(y.controls,start.vals=start.vals.controls)
  res$mu.cases <- em.cases$mu
  res$sig.cases <- em.cases$sig
  res$pi.cases <- em.cases$pi
  res$max.loglike.cases <- em.cases$max.loglike + (res$lambda - 1)*sum(log(x.cases))
  res$mu.controls <- em.controls$mu
  res$sig.controls <- em.controls$sig
  res$pi.controls <- em.controls$pi
  res$max.loglike.controls <- em.controls$max.loglike + (res$lambda - 1)*sum(log(x.controls))
  res$max.loglike <- res$max.loglike.controls+res$max.loglike.cases
  
  
  cs1.estimates <- rnorm(1e6,res$mu.cases[1],res$sig.cases[1])
  cs2.estimates <- rnorm(1e6,res$mu.cases[2],res$sig.cases[2])
  ctrl1.estimates <- rnorm(1e6,res$mu.controls[1],res$sig.controls[1])
  ctrl2.estimates <- rnorm(1e6,res$mu.controls[2],res$sig.controls[2])
  
  
  if(res$lambda>0){
    cs1.estimates[which(cs1.estimates<=-1/res$lambda)]<-NA
    cs2.estimates[which(cs2.estimates<=-1/res$lambda)]<-NA
    ctrl1.estimates[which(ctrl1.estimates<=-1/res$lambda)]<-NA
    ctrl2.estimates[which(ctrl2.estimates<=-1/res$lambda)]<-NA
  }
  
  if(res$lambda<0){
    cs1.estimates[which(cs1.estimates>=-1/res$lambda)]<-NA
    cs2.estimates[which(cs2.estimates>=-1/res$lambda)]<-NA
    ctrl1.estimates[which(ctrl1.estimates>=-1/res$lambda)]<-NA
    ctrl2.estimates[which(ctrl2.estimates>=-1/res$lambda)]<-NA
  }
  
  monte.carlo.controls1 <- boxcox.inv(ctrl1.estimates,res$lambda)
  monte.carlo.controls2 <- boxcox.inv(ctrl2.estimates,res$lambda)  
  monte.carlo.cases1 <- boxcox.inv(cs1.estimates,res$lambda)
  monte.carlo.cases2 <- boxcox.inv(cs2.estimates,res$lambda)
  res$mu.cases.unt <- c(mean(monte.carlo.cases1,na.rm=T),mean(monte.carlo.cases2,na.rm=T))
  res$sig.cases.unt <- c(sd(monte.carlo.cases1,na.rm=T),sd(monte.carlo.cases2,na.rm=T))
  res$mu.controls.unt <- c(mean(monte.carlo.controls1,na.rm=T),mean(monte.carlo.controls2,na.rm=T))
  res$sig.controls.unt <- c(sd(monte.carlo.controls1,na.rm=T),sd(monte.carlo.controls2,na.rm=T))
  
  res$case=x.cases
  res$control=x.controls
  res$time <- proc.time()[3] - t0
  
  res
}
