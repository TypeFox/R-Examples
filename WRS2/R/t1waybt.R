t1waybt <-
function(formula, data, tr = 0.2, nboot = 599){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  if (tr==0.5) warning("Comparing medians should not be done with this function!")
  
  grp <- 1:length(x) 
  J <- length(x)

  for(j in 1:J){
    temp<-x[[j]]
    x[[j]]<-temp[!is.na(temp)] # Remove any missing values.
  }
  bvec  <-  array(0,c(J,2,nboot))
  hval <- vector("numeric",J)
  #if(SEED)set.seed(2) # set seed of random number generator so that
  for(j in 1:J){
    hval[j] <- length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
    xcen <- x[[grp[j]]]-mean(x[[grp[j]]],tr)
    data <- matrix(sample(xcen,size=length(x[[grp[j]]])*nboot,replace=TRUE),nrow=nboot)
    bvec[j,,] <- apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
  }
  m1 <- bvec[,1,]  # J by nboot matrix containing the bootstrap trimmed means
  m2 <- bvec[,2,]  # J by nboot matrix containing the bootstrap sq standard errors
  wvec <- 1/m2  # J by nboot matrix of w values
  uval <- apply(wvec,2,sum)  # Vector having length nboot
  blob <- wvec*m1
  xtil <- apply(blob,2,sum)/uval # nboot vector of xtil values
  blob1 <- matrix(0,J,nboot)
  for (j in 1:J)blob1[j,] <- wvec[j,]*(m1[j,]-xtil)^2
  avec <- apply(blob1,2,sum)/(length(x)-1)
  blob2 <- (1-wvec/uval)^2/(hval-1)
  cvec <- apply(blob2,2,sum)
  cvec <- 2*(length(x)-2)*cvec/(length(x)^2-1)
  testb <- avec/(cvec+1)
  ct <- sum(is.na(testb))
  if(ct>0){
    warning("Some bootstrap estimates of the test statistic could not be computed.")
  }
  neff <- sum(!is.na(testb))
  test <- t1wayv2(x,tr=tr,grp=grp)
  pval <- mean(test$TEST<=testb,na.rm=TRUE)
  result <- list(test=test$TEST,p.value=pval,Var.Explained=test$Var.Explained,Effect.Size=test$Effect.Size, nboot.eff = neff, call = cl)
  class(result) <- c("t1waybt")
  result
}
