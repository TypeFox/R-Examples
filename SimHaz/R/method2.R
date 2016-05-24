#' @export
tdSim.method2 <-function(N, duration, lambda12, lambda23=NULL, lambda13, HR=NULL, 
                         exp.prop,rateC, min.futime = 0, min.postexp.futime = 0){
  try(if(is.null(lambda23) & is.null(HR)){stop("either lambda23 or HR(Hazard ratio) must be set")})
  if(is.null(lambda23) & !is.null(HR)){
    lambda23 = lambda13 * HR
  }
  expose<-rbinom(n=N,size=1,prob=exp.prop)
  t12 <- rep(NA, N)
  t13 <- rep(NA, N)
  t23 <- rep(NA, N)
  t12[as.logical(expose)] <- rexp(n=sum(expose), rate=lambda12)
  t13[as.logical(expose)] <- Inf
  t23[as.logical(expose)] <- rexp(n=sum(expose), rate=lambda23)
  t12[!as.logical(expose)] <- Inf
  t13[!as.logical(expose)] <- rexp(n=N-sum(expose), rate=lambda13)
  t23[!as.logical(expose)] <- Inf
  #  C <- runif(n=N, min=0, max=1000)
  C <- rexp(n=N, rate=rateC)
  C=pmin(C,rep(duration,length(C)))
  df <- data.frame(matrix(ncol = 5, nrow = N))
  df[,1] = seq(N)
  index <- C < pmin(t12,t13)
  if(sum(index) > 0){
    df[index, c(2,3)] <- data.frame(C[index],C[index])
    df[index, c(4,5)] <- cbind(rep(0,sum(index)),rep(0,sum(index)))}
  index <- t13 <= pmin(t12, C)
  if(sum(index) > 0){
    df[index, c(2,3)] <- data.frame(t13[index], t13[index])
    df[index, c(4,5)] <- cbind(rep(0,sum(index)),rep(1,sum(index)))}
  index <- (t12 <= pmin(t13, C) & (t12+t23) > C)
  if(sum(index) > 0){
    df[index, c(2,3)] <- data.frame(t12[index], C[index])
    df[index, c(4,5)] <- cbind(rep(1,sum(index)),rep(0,sum(index)))}
  index <- (t12 < t13 & (t12+t23) <= C)
  if(sum(index) > 0){
    df[index, c(2,3)] <- data.frame(t12[index], t12[index]+t23[index])
    df[index, c(4,5)] <- cbind(rep(1,sum(index)),rep(1,sum(index)))}
  colnames(df) = c('id',"exp.time",'end','exp','status')
  if(min.futime>0){
    df <- df[df$end>min.futime,]
    df$id <- seq(nrow(df))
  }
  df_exp <- df[df$exp==1, ]
  if(min.postexp.futime>0){
    if(sum(df_exp$end-df_exp$exp.time > min.postexp.futime) == 0){
      print('no exposure left')
    }
    df_exp <- df_exp[df_exp$end-df_exp$exp.time > min.postexp.futime,]
  }
  df1 <- df_exp
  df1$start <- 0
  df1$stop <- df1$exp.time
  df1$status <- 0
  df1$x <- 0
  df2 <- df_exp
  df2$start <- df2$exp.time
  df2$stop <- df2$end
  df2$status <- df2$status
  df2$x <- 1
  df3 <- df[df$exp==0, ]
  df3$start <- 0
  df3$stop <- df3$end
  df3$status <- df3$status
  df3$x <- 0
  merge_exposed <- merge(df1[,c("id","start","stop","status","x")],df2[,c("id","start","stop","status","x")],all.x=TRUE,all.y=TRUE)
  merged_df <- merge(merge_exposed, df3[,c("id","start","stop","status","x")],all.x=TRUE,all.y=TRUE)
  return(merged_df)
}

#' @export
getpower.method2=function(nSim=500, N, duration=24, scenario,lambda12, lambda23=NULL, lambda13, HR=NULL,exp.prop,rateC, 
                          min.futime, min.postexp.futime,output.fn, simu.plot=FALSE) 
{ set.seed(999)
  try(if(is.null(lambda23) & is.null(HR)){stop("either lambda23 or HR(Hazard ratio) must be set")})
  if(is.null(lambda23) & !is.null(HR)){
    lambda23 = lambda13 * HR
  }
  #N=400;duration=24;medTTEC=24;rho=1;medTTC=14;b=0.3;er=0.2;s1=1;s2=6;fA=4;fB=4
  res=matrix(0,nSim,10)
  colnames(res)=c("N.eff","N.effexp.p","betahat","HR","signif","events",
                  "events_c","events_exp","medsurvt_c","medsurvt_exp")
  alpha=.05
  if(simu.plot){
    dat <- tdSim.method2(N, duration, lambda12=lambda12, lambda23=lambda23, lambda13=lambda13, exp.prop=exp.prop,rateC=rateC, min.futime=min.futime, min.postexp.futime=min.postexp.futime)
    plot_simuData(dat)
  }
  for(k in 1:nSim)
  {
    dat <- tdSim.method2(N, duration, lambda12=lambda12, lambda23=lambda23, lambda13=lambda13, exp.prop=exp.prop,rateC=rateC, min.futime=min.futime, min.postexp.futime=min.postexp.futime)    
    fit <- coxph(Surv(start,stop, status) ~ factor(x), data=dat)
    sfit <- survfit(Surv(start,stop, status) ~ factor(x), data=dat)
    res[k,"N.eff"] <- length(unique(dat$id))
    res[k,"N.effexp.p"] <- sum(dat$x)/length(unique(dat$id))
    res[k,"betahat"] <- summary(fit)$coef[,"coef"]
    res[k,"HR"] <- summary(fit)$coef[,"exp(coef)"]
    res[k,"signif"] <- ifelse(summary(fit)$coef[,"Pr(>|z|)"]<alpha,1,0)
    res[k,"events"] <- sum(dat$status)
    res[k,"events_c"] <- summary(sfit)$table[1,'events']
    res[k,"events_exp"] <- summary(sfit)$table[2,'events']
    res[k,"medsurvt_c"] <- summary(sfit)$table[1,'median']
    res[k,"medsurvt_exp"] <- summary(sfit)$table[2,'median']
  }
  df=data.frame(i_scenario=scenario,
                i_N=N,
                i_min.futime=min.futime,
                i_min.postexp.futime=min.postexp.futime,
                i_exp.prop=exp.prop,
                i_lambda12=lambda12,
                i_lambda23=lambda23,
                i_lambda13=lambda13,
                i_rateC=rateC,
                N_eff=mean(res[,"N.eff"]),
                N_effexp_p=mean(res[,"N.effexp.p"]),
                bhat=mean(res[,"betahat"]),
                HR=mean(res[,"HR"]),                     
                d=mean(res[,"events"]),
                d_c=mean(res[,"events_c"]),
                d_exp=mean(res[,"events_exp"]),
                mst_c=mean(na.omit(res[,"medsurvt_c"])),
                mst_exp=mean(na.omit(res[,"medsurvt_exp"])),
                pow=mean(res[,"signif"])
  )
  if(file.exists(output.fn)){
    write.table(df,file=output.fn,row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
  }
  else{
    write.table(df,file=output.fn,row.names=FALSE,col.names=TRUE,sep=",")
  }
  return(df)
}
