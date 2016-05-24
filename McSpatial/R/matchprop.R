matchprop <- function(form,data=NULL,distance="logit",discard="both",reestimate="FALSE",m.order="none",nclose=0,ytreat=1) {
  if (distance=="Mahalanobis"|distance=="mahalanobis"|distance=="Mahal") {distance="mahal"}
 
  if (identical(data,NULL)) {data <- model.frame(form)}
  yname <- names(model.frame(form,data=data))[1]
  ynum <- which(names(data)==yname)
  y <- data[,ynum]
  yvect <- levels(factor(y))
  n = nrow(data)
  data$origobs <- seq(1:nrow(data))
  data$matchobs <- array(0,dim=n)

  if (identical(yvect,c("FALSE","TRUE"))) {yvect=c(0,1)}
  if (identical(yvect,c("TRUE","FALSE"))) {yvect=c(1,0)}

  newform <- as.formula(form,env=smalldata)	

  for (j in yvect[yvect!=ytreat]) {
    smalldata <- data[y==ytreat|y==j,]
    smalldata[,ynum] <- ifelse(smalldata[,ynum]==ytreat,1,0)

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
    yvar <- yvar[sampvar==TRUE]
    if (reestimate==TRUE&distance!="mahal") {
      fit <- glm(newform,data=smalldata,family=binomial(link=distance))
      p <- fitted(fit)
    }
    if (reestimate==TRUE&distance=="mahal") {
      xmat <- as.matrix(model.matrix(newform,data=smalldata)[,-1])
      p <- mahalanobis(xmat,colMeans(xmat),cov(xmat))
    }

    smalldata <- data.frame(smalldata$origobs,p,yvar)
    names(smalldata) <- c("origobs","p","yvar")
    o <- seq(1,nrow(smalldata))
    if (m.order=="decreasing")  {o <- order(smalldata$p,decreasing=TRUE)  }
    if (m.order=="increasing") {o <- order(smalldata$p,decreasing=FALSE) }
    if (m.order=="random")   {o <- order(rnorm(nrow(smalldata))) }
    if (m.order!="none")  {smalldata <- smalldata[o,] }
    data0 <- smalldata[smalldata$yvar==0,]
    data1 <- smalldata[smalldata$yvar==1,]
    n0 = nrow(data0)
    n1 = nrow(data1)

    obs0 <- array(0,dim=n1)
    obs1 <- array(0,dim=n1)
    maxd = (max(data1$p) - min(data0$p)) + 1
    dist01 <- array(maxp,dim=n1)
    for (i in seq(1:n1)) { 
      jj = which.min(abs(data1$p[i]-data0$p))
      obs0[i] = data0$origobs[jj]
      obs1[i] = data1$origobs[i]
      dist01[i] = abs(data1$p[i]-data0$p[jj])
      if (nrow(data0)==1) {break}
      data0 <- data0[-jj,]
    }
    if (nclose>0) {
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

   

