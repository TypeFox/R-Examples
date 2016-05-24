mcppb20 <- function(formula, data, tr = 0.2, nboot = 599, crit = NA){
  #
 
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  con=0
  alpha=.05
  grp=NA
  WIN=FALSE
  win=.1
  
  con<-as.matrix(con)
  if(is.matrix(x)){
    xx<-list()
    for(i in 1:ncol(x)){
      xx[[i]]<-x[,i]
    }
    x<-xx
  }
  if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
  if(!is.na(sum(grp))){  # Only analyze specified groups.
    xx<-list()
    for(i in 1:length(grp))xx[[i]]<-x[[grp[1]]]
    x<-xx
  }
  J<-length(x)
  tempn<-0
  for(j in 1:J){
    temp<-x[[j]]
    temp<-temp[!is.na(temp)] # Remove missing values.
    tempn[j]<-length(temp)
    x[[j]]<-temp
  }
  Jm<-J-1
  d<-ifelse(sum(con^2)==0,(J^2-J)/2,ncol(con))
  if(is.na(crit) && tr != .2){
    stop("A critical value must be specified when the amount of trimming differs from .2")
  }
  if(WIN){
    if(tr < .2){
      warning("When Winsorizing, the amount of trimming should be at least .2")
    }
    if(win > tr)stop("Amount of Winsorizing must <= amount of trimming")
    if(min(tempn) < 15){
      warning("Winsorizing with sample sizes less than 15 can result in poor control over the probability of a Type I error.")
    }
    for (j in 1:J){
      x[[j]]<-winval(x[[j]],win)
    }
  }
  if(is.na(crit)){
    if(d==1)crit<-alpha/2
    if(d==2 && alpha==.05 && nboot==1000)crit<-.014
    if(d==2 && alpha==.05 && nboot==2000)crit<-.014
    if(d==3 && alpha==.05 && nboot==1000)crit<-.009
    if(d==3 && alpha==.05 && nboot==2000)crit<-.0085
    if(d==3 && alpha==.025 && nboot==1000)crit<-.004
    if(d==3 && alpha==.025 && nboot==2000)crit<-.004
    if(d==3 && alpha==.01 && nboot==1000)crit<-.001
    if(d==3 && alpha==.01 && nboot==2000)crit<-.001
    if(d==4 && alpha==.05 && nboot==2000)crit<-.007
    if(d==5 && alpha==.05 && nboot==2000)crit<-.006
    if(d==6 && alpha==.05 && nboot==1000)crit<-.004
    if(d==6 && alpha==.05 && nboot==2000)crit<-.0045
    if(d==6 && alpha==.025 && nboot==1000)crit<-.002
    if(d==6 && alpha==.025 && nboot==2000)crit<-.0015
    if(d==6 && alpha==.01 && nboot==2000)crit<-.0005
    if(d==10 && alpha==.05 && nboot<=2000)crit<-.002
    if(d==10 && alpha==.05 && nboot==3000)crit<-.0023
    if(d==10 && alpha==.025 && nboot<=2000)crit<-.0005
    if(d==10 && alpha==.025 && nboot==3000)crit<-.001
    if(d==15 && alpha==.05 && nboot==2000)crit<-.0016
    if(d==15 && alpha==.025 && nboot==2000)crit<-.0005
    if(d==15 && alpha==.05 && nboot==5000)crit<-.0026
    if(d==15 && alpha==.025 && nboot==5000)crit<-.0006
  }
  if(is.na(crit) && alpha==.05)crit<-0.0268660714*(1/d)-0.0003321429
  if(is.na(crit))crit<-alpha/(2*d)
  if(d> 10 && nboot <5000){
    warning("Suggest using nboot = 5000 when the number of contrasts exceeds 10.")
  }
  icl<-round(crit*nboot)+1
  icu<-round((1-crit)*nboot)
  if(sum(con^2)==0){
    con<-matrix(0,J,d)
    id<-0
    for (j in 1:Jm){
      jp<-j+1
      for (k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  psihat<-matrix(0,ncol(con),6)
  dimnames(psihat)<-list(NULL,c("con.num","psihat","se","ci.lower",
                                "ci.upper","p-value"))
  bvec<-matrix(NA,nrow=J,ncol=nboot)
  #set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  for(j in 1:J){
 #   paste("Working on group ",j)
    data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
    bvec[j,]<-apply(data,1,mean,tr) # Bootstrapped trimmed means for jth group
  }
  test<-NA
  for (d in 1:ncol(con)){
    top<-0
    for (i in 1:J){
      top<-top+con[i,d]*bvec[i,]
    }
    test[d]<-(sum(top>0)+.5*sum(top==0))/nboot
    test[d]<-min(test[d],1-test[d])
    top<-sort(top)
    psihat[d,4]<-top[icl]
    psihat[d,5]<-top[icu]
  }
  for (d in 1:ncol(con)){
    psihat[d,1]<-d
    testit<-lincon1(x,con[,d],tr,pr=FALSE)
    psihat[d,6]<-2*test[d]
    psihat[d,2]<-testit$psihat[1,2]
    psihat[d,3]<-testit$test[1,4]
  }
 
 list(psihat=psihat,crit.p.value=2*crit,con=con)
 
 fnames <- as.character(unique(mf[,2]))
 
 groups <- t(apply(con, 2, function(cc) {
   c(which(cc == 1), which(cc == -1))
 }))
 
 psihat1 <- cbind(groups, psihat[, -c(1, 3)])
 colnames(psihat1)[1:2] <- c("Group")
 
 result <- list(comp = psihat1, fnames = fnames, call = cl)
 class(result) <- "mcp1"
 result
 
}
