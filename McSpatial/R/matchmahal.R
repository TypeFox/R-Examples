matchmahal <- function(form,data=NULL,discard="none",distance="logit",m.order="none",nclose=0,ytreat=1) {
  if (distance=="Mahalanobis"|distance=="mahalanobis"|distance=="Mahal") {distance="mahal"}
 
  if (identical(data,NULL)) {data <- model.frame(form)}
  yname <- names(model.frame(form,data=data))[1]
  ynum <- which(names(data)==yname)
  y <- data[,ynum]
  yvect <- levels(factor(y))
  n = nrow(data)
  data$origobs <- seq(1:nrow(data))
  data$matchobs <- array(0,dim=n)

  newform <- as.formula(form,env=smalldata)	
  if (identical(yvect,c("FALSE","TRUE"))) {yvect=c(0,1)}
  if (identical(yvect,c("TRUE","FALSE"))) {yvect=c(1,0)}

  for (j in yvect[yvect!=ytreat]) {
    smalldata <- data[y==ytreat|y==j,]
    smalldata[,ynum] <- ifelse(smalldata[,ynum]==ytreat,1,0)

    if (distance!="none") {
      if (distance!="mahal") {
        fit <- glm(newform,data=smalldata,family=binomial(link=distance))
        p <- fitted(fit)
      }
      if (distance=="mahal") {
        xmat <- as.matrix(model.matrix(newform,data=smalldata)[,-1])
        p <- mahalanobis(xmat,colMeans(xmat),cov(xmat))
      }
      yvar <- smalldata[,ynum]

      minp = max(min(p[yvar==1]), min(p[yvar==0]))
      maxp = min(max(p[yvar==1]), max(p[yvar==0]))
      discard1 <- yvar==1&(p<minp|p>maxp)
      discard0 <- yvar==0&(p<minp|p>maxp)
      sampvar <- array(TRUE,dim=nrow(smalldata)) 
      if (discard=="both") {sampvar <- ifelse(discard1==TRUE|discard0==TRUE, FALSE, sampvar)}
      if (discard=="control") {sampvar <- ifelse(discard0==TRUE, FALSE, sampvar)}
      if (discard=="treat") {sampvar <- ifelse(discard1==TRUE, FALSE, sampvar)}
      smalldata <- smalldata[sampvar==TRUE,]
      p <- p[sampvar==TRUE]

      if (m.order=="largest")  {o <- order(p,decreasing=TRUE)  }
      if (m.order=="smallest") {o <- order(p,decreasing=FALSE) }
    }
      if (m.order=="random")   {o <- order(rnorm(nrow(smalldata))) }
      if (m.order!="none")  {smalldata <- smalldata[o,] }

    xmat0 <- as.matrix(model.matrix(newform,data=smalldata)[smalldata[,ynum]==0, -1])
    xmat1 <- as.matrix(model.matrix(newform,data=smalldata)[smalldata[,ynum]==1, -1])
    n0 = nrow(xmat0)
    n1 = nrow(xmat1)
    origobs0 <- smalldata$origobs[smalldata[,ynum]==0]
    origobs1 <- smalldata$origobs[smalldata[,ynum]==1]

    xname0 <- colnames(xmat0)
    xname1 <- colnames(xmat1)
    xname <- intersect(xname0,xname1)
    if (length(xname)<length(xname0)|length(xname)<length(xname1)) {
      xmat0 <- as.matrix(xmat0[,xname])
      xmat1 <- as.matrix(xmat1[,xname])
      xname <- sort(union(xname0[!xname0%in%xname],xname1[!xname1%in%xname]))
      cat("The following variables were dropped from the distance calculations:","\n")
      cat(xname,"\n")
    }
    cmat <- cov(rbind(xmat0,xmat1))
    obs0 <- array(0,dim=n1)
    obs1 <- array(0,dim=n1)
    maxd = -1
    dist01 <- array(maxd,dim=n1)
    for (i in seq(1:n1)) { 
      mmat <- mahalanobis(xmat0,xmat1[i,],cmat)
      jj = which.min(mmat)
      obs0[i] = origobs0[jj]
      obs1[i] = origobs1[i]
      dist01[i] = mmat[jj]
      if (nrow(xmat0)==1) {break}
      xmat0 <- as.matrix(xmat0[-jj,])
      if (nrow(xmat0)>1&ncol(xmat0)==1&ncol(xmat1)>1){xmat0 <- t(xmat0)}
      origobs0 <- origobs0[-jj]
    }
    if (nclose>0) {
      maxd = max(dist01)
      dist01 <- ifelse(dist01<0,maxd,dist01)
      o <- order(dist01)
      obs0[o>nclose|dist01==maxp] <- 0
      obs1[o>nclose|dist01==maxp] <- 0
    }
    obs0 <- obs0[obs0>0]
    obs1 <- obs1[obs1>0]
    data$matchobs[obs1] <- data$origobs[obs1]
    data$matchobs[obs0] <- obs1
  }

  data <- data[data$matchobs>0,]
  return(data)
}

