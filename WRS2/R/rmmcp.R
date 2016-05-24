rmmcp<-function(y, groups, blocks, tr = 0.2){
  #
  # MCP on trimmed means with FWE controlled with Rom's method
  #
  
  cols1 <- deparse(substitute(y))
  cols2 <- deparse(substitute(groups))
  cols3 <- deparse(substitute(blocks))
  dat <- data.frame(y, groups, blocks)
  colnames(dat) <- c(cols1, cols2, cols3)
  cl <- match.call()
  
  x <- reshape(dat, idvar = cols3, timevar = cols2, direction = "wide")[-1]  ## wide format
  grp <- c(1:length(x))
  
  con = 0
  alpha = 0.05
  dif=TRUE
  flagcon=F
  if(!is.matrix(x))x<-matl(x)
  if(!is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
  con<-as.matrix(con)
  J<-ncol(x)
  xbar<-vector("numeric",J)
  x<-elimna(x)  # Remove missing values
  nval<-nrow(x)
  h1<-nrow(x)-2*floor(tr*nrow(x))
  df<-h1-1
  for(j in 1: J)xbar[j]<-mean(x[,j],tr)
  if(sum(con^2!=0))CC<-ncol(con)
  if(sum(con^2)==0)CC<-(J^2-J)/2
  ncon<-CC
  if(alpha==.05){
    dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
    if(ncon > 10){
      avec<-.05/c(11:ncon)
      dvec<-c(dvec,avec)
    }}
  if(alpha==.01){
    dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
    if(ncon > 10){
      avec<-.01/c(11:ncon)
      dvec<-c(dvec,avec)
    }}
  if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
  if(sum(con^2)==0){
    flagcon<-T
    psihat<-matrix(0,CC,5)
    dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
    test<-matrix(NA,CC,6)
    dimnames(test)<-list(NULL,c("Group","Group","test","p.value","p.crit","se"))
    temp1<-0
    jcom<-0
    for (j in 1:J){
      for (k in 1:J){
        if (j < k){
          jcom<-jcom+1
          q1<-(nrow(x)-1)*winvar(x[,j],tr)
          q2<-(nrow(x)-1)*winvar(x[,k],tr)
          q3<-(nrow(x)-1)*wincor(x[,j],x[,k],tr)$cov
          sejk<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
          if(!dif){
            test[jcom,6]<-sejk
            test[jcom,3]<-(xbar[j]-xbar[k])/sejk
            temp1[jcom]<-2 * (1 - pt(abs(test[jcom,3]), df))
            test[jcom,4]<-temp1[jcom]
            psihat[jcom,1]<-j
            psihat[jcom,2]<-k
            test[jcom,1]<-j
            test[jcom,2]<-k
            psihat[jcom,3]<-(xbar[j]-xbar[k])
          }
          if(dif){
            dv<-x[,j]-x[,k]
            test[jcom,6]<-trimse(dv,tr)
            temp<-trimci(dv,alpha=alpha/CC,pr=FALSE,tr=tr)
            test[jcom,3]<-temp$test.stat
            temp1[jcom]<-temp$p.value
            test[jcom,4]<-temp1[jcom]
            psihat[jcom,1]<-j
            psihat[jcom,2]<-k
            test[jcom,1]<-j
            test[jcom,2]<-k
            psihat[jcom,3]<-mean(dv,tr=tr)
            psihat[jcom,4]<-temp$ci[1]
            psihat[jcom,5]<-temp$ci[2]
          }
        }}}
    temp2<-order(0-temp1)
    zvec<-dvec[1:ncon]
    sigvec<-(test[temp2]>=zvec)
    if(sum(sigvec)<ncon){
      dd<-ncon-sum(sigvec) #number that are sig.
      ddd<-sum(sigvec)+1
      zvec[ddd:ncon]<-dvec[ddd]
    }
    test[temp2,5]<-zvec
    if(!dif){
      psihat[,4]<-psihat[,3]-qt(1-alpha/(2*CC),df)*test[,6]
      psihat[,5]<-psihat[,3]+qt(1-alpha/(2*CC),df)*test[,6]
    }}
  if(sum(con^2)>0){
    if(nrow(con)!=ncol(x))warning("The number of groups does not match the number of contrast coefficients.")
    ncon<-ncol(con)
    psihat<-matrix(0,ncol(con),4)
    dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
    test<-matrix(0,ncol(con),5)
    dimnames(test)<-list(NULL,c("con.num","test","p.value","p.crit","se"))
    temp1<-NA
    for (d in 1:ncol(con)){
      psihat[d,1]<-d
      if(!dif){
        psihat[d,2]<-sum(con[,d]*xbar)
        sejk<-0
        for(j in 1:J){
          for(k in 1:J){
            djk<-(nval-1)*wincor(x[,j],x[,k], tr)$cov/(h1*(h1-1))
            sejk<-sejk+con[j,d]*con[k,d]*djk
          }}
        sejk<-sqrt(sejk)
        test[d,1]<-d
        test[d,2]<-sum(con[,d]*xbar)/sejk
        test[d,5]<-sejk
        temp1[d]<-2 * (1 - pt(abs(test[d,2]), df))
      }
      if(dif){
        for(j in 1:J){
          if(j==1)dval<-con[j,d]*x[,j]
          if(j>1)dval<-dval+con[j,d]*x[,j]
        }
        temp1[d]<-trimci(dval,tr=tr,pr=FALSE)$p.value
        test[d,1]<-d
        test[d,2]<-trimci(dval,tr=tr,pr=FALSE)$test.stat
        test[d,5]<-trimse(dval,tr=tr)
        psihat[d,2]<-mean(dval,tr=tr)
      }}
    test[,3]<-temp1
    temp2<-order(0-temp1)
    zvec<-dvec[1:ncon]
    sigvec<-(test[temp2,3]>=zvec)
    if(sum(sigvec)<ncon){
      dd<-ncon-sum(sigvec) #number that are sig.
      ddd<-sum(sigvec)+1
      #zvec[ddd:ncon]<-dvec[ddd]
    }
    test[temp2,4]<-zvec
    psihat[,3]<-psihat[,2]-qt(1-test[,4]/2,df)*test[,5]
    psihat[,4]<-psihat[,2]+qt(1-test[,4]/2,df)*test[,5]
  }
  if(flagcon)num.sig<-sum(test[,4]<=test[,5])
  if(!flagcon)num.sig<-sum(test[,3]<=test[,4])

  fnames <- as.character(unique(groups))
  psihat1 <- cbind(psihat, test[,4:5])
  result <- list(comp = psihat1, fnames = fnames, call = cl)
  class(result) <- "mcp2"
  result
}
