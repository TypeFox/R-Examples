dbkGrad  <- function(obsq, limx, limy, exposures=NULL, transformation=c("none","log","logit","Gompertz"), bwtypex=c("FX","VC","EX"), bwtypey=c("FX","VC","EX"), adaptx=c("a","b","ab"), adapty=c("a","b","ab"), hx=0.002, hy=0.002, sx=0.2,sy=0.2, cvres=c("propres","res"), cvhx=FALSE,cvhy=FALSE, cvsx=FALSE,cvsy=FALSE, alpha=0.05){
  obsq <- as.matrix(obsq)
  if (missing(limx)) limx <-c(1,nrow(obsq)) 
  if (length(limx)==1) limx<- c(limx,limx)
  if (missing(limy)) limy <-c(1,ncol(obsq))
  if (length(limy)==1) limy<- c(limy,limy)
  if(any(is.na(as.integer(rownames(obsq))))| is.null(rownames(obsq))){
    print("Warning: rownames have been overwritten")
    rownames(obsq) <- 0:(nrow(obsq)-1)
  }
  if(any(is.na(as.integer(colnames(obsq))))| is.null(colnames(obsq))){
    print("Warning: colnames have been overwritten")
    colnames(obsq) <- 0:(ncol(obsq)-1)
  }
  obsq <- .extract(obsq,rows=limx[1]:limx[2], columns=limy[1]:limy[2],byage=TRUE)
  if (is.null(dimnames(obsq))) dimnames(obsq) <- list(limx[1]:limx[2],1:(ncol(obsq)))
  if(!is.null(exposures)) exposures <- .extract(as.matrix(exposures),rows=limx[1]:limx[2], columns=limy[1]:limy[2],byage=TRUE)
  transformation <- match.arg(transformation)
  bwtypex <- match.arg(bwtypex)
  bwtypey <- match.arg(bwtypey)
  cvres     <- match.arg(cvres)
  if(any(bwtypex ==c("EX","VC")) | any(bwtypey ==c("EX","VC"))){
    if (is.null(exposures)) stop("no value specified for the exposures parameter") 
    else if(any(dim(exposures)!=dim(obsq))) stop("obsq and exposures have different dimension")
  }
  adaptx <- match.arg(adaptx)
  adapty <- match.arg(adapty)
  maxx      <- nrow(obsq)
  maxy      <- ncol(obsq)
  iandj     <- expand.grid(i=1:maxx,j=1:maxy)
  if (bwtypex=="VC"){
    xweights     <- .vcf(rowSums(exposures),rowSums(obsq*exposures)/rowSums(exposures))
    xweights     <- xweights/sum(xweights)
  }
  else if (bwtypex=="EX"){
    xweights <- rowSums(sum(exposures)/exposures)
    xweights <- xweights/max(xweights)
  }
  else  {
    xweights <- rep(1, maxx)
  }

  if (bwtypey=="VC"){
    yweights <- .vcf(colSums(exposures),colSums(obsq*exposures)/colSums(exposures))
    yweights     <- xweights/sum(xweights)
  }
  else if (bwtypey=="EX"){
    yweights <- colSums(sum(exposures)/exposures)
    yweights <- yweights/max(yweights)
  }
  else if (bwtypey=="FX"){
    yweights <- rep(1, maxy)
  }

  obsq  <- .transformation(obsq,transformation)
                        
  if (any(cvhx,cvhy,cvsx,cvsy)) {
    nls.out <-  nls.lm(par=list(hx=hx,hy=hy,sx=sx,sy=sy),hx=hx,hy=hy,sx=sx,sy=sy,cvhx=cvhx,cvhy=cvhy, cvsx=cvsx,cvsy=cvsy, maxx=maxx,maxy=maxy,cvres=cvres, obsq=obsq,xweights=xweights, yweights=yweights,iandj=iandj,adaptx=adaptx,adapty=adapty, fn=.dkCV, lower=c(1e-200,1e-200,0,0), upper=c(Inf,Inf,1,1), control = nls.lm.control(maxiter=1000,nprint=2))
    hx <- nls.out$par$hx
    hy <- nls.out$par$hy
    sx <- nls.out$par$sx
    sy <- nls.out$par$sy
    cvRSS <- nls.out$rsstrace[nls.out$niter]
    if (is.nan(cvRSS)){
      cat("Warning: cross-validation returned nan; specify different initial values.")}
  }
  else{
    cvRSS <- .dkCV(par=NULL, hx=hx,hy=hy,sx=sx,sy=sy,cvhx=FALSE,cvhy=FALSE, cvsx=FALSE,cvsy=FALSE, maxx=maxx,maxy=maxy,cvres=cvres, obsq=obsq,xweights=xweights, yweights=yweights,iandj=iandj,adaptx=adaptx,adapty=adapty)
    if (any(is.nan(cvRSS))){
       cat("Warning: for the specified h", ifelse(bwtypex!="FX","and s,",""),"no smoothing was applied at ages:")
      cat((0:limx[2])[is.nan(cvRSS)],sep = ",")
    }
   cvRSS <- sum(cvRSS)^2
  }
  q         <-  lapply(1:nrow(iandj),FUN=".ddbeta2",iandj=iandj,maxx=maxx,maxy=maxy,xweights=xweights,yweights=yweights,obsq=obsq,exposures=exposures, cv=FALSE,adaptx=adaptx,adapty=adapty,hx=hx,hy=hy,sx=sx,sy=sy)
  qxest     <-  matrix(unlist(sapply(q,"[",1)),maxx,maxy, dimnames=dimnames(obsq))
  kernels   <-  array(unlist(sapply(q,"[",2)),c(maxx,maxy,maxx*maxy))
  qxest     <- .btransformation(qxest,transformation) 
  obsq      <- .btransformation(obsq,transformation)

  result <- list(fitted.values=qxest, hx=hx, sx=sx, residuals=obsq-qxest, obsq=obsq, limx=limx, limy=limy,cvRSS=cvRSS, exposures=exposures, call=match.call()) #,cvRSS=cvRSS, kernels=kernels
  if(!is.null(exposures) && all(dim(exposures)==dim(obsq))){
    varqx <- qxest*(1-qxest)/exposures
    qxstderr <- matrix(unlist(lapply(1:(maxx*maxy),FUN=".qx.stderr",kernels, varqx)),maxx,maxy)
    
    halfCIlength <- qnorm(1-alpha/2) *qxstderr
    upperbound <- matrix(unlist(lapply(qxest+halfCIlength,function(x)min(x,1))),maxx,maxy)
    lowerbound <- matrix(unlist(lapply(qxest-halfCIlength,function(x)max(x,1e-100))),maxx,maxy)
    halfCIlength <- qnorm(1-alpha/(maxx)/2) *qxstderr
    bonferroniUB <- matrix(unlist(lapply(qxest+halfCIlength,function(x)min(x,1))),maxx,maxy)
    bonferroniLB <- matrix(unlist(lapply(qxest-halfCIlength,function(x)max(x,1e-100))),maxx,maxy)
    halfCIlength <- qnorm((1-alpha)^(1/2)/2) *qxstderr
    sidakUB <- matrix(unlist(lapply(qxest+halfCIlength,function(x)min(x,1))),maxx,maxy)
    sidakLB <- matrix(unlist(lapply(qxest-halfCIlength,function(x)max(x,1e-100))),maxx,maxy)
    result <- c(result,list(upperbound=upperbound,lowerbound=lowerbound,bonferroniUB=bonferroniUB,bonferroniLB=bonferroniLB,sidakUB=sidakUB,sidakLB=sidakLB))
  } 
  class(result) <- "dbkGrad"
  result
}

.trans <- function(xmin=0,xmax){
  startgrid <- xmin:xmax
  transgrid <- (startgrid-xmin+0.5)/(xmax-xmin+1)
  return(transgrid)
}
.transformation <- function(q,transformation)
{
  if (transformation %in%  c("log","logit","Gompertz"))  q <- replace(q, q==0, 1e-6)
  if (transformation %in%  c("logit","Gompertz"))  q <- replace(q, q==1, 1-1e-6)
  out <- switch(transformation,
          "none" = q,
          "log" = log(q),
          "logit" = log(q/(1-q)),
          "Gompertz" = log(-log(1-q)),
          "arcsin" =  asin(sqrt(q))
          )
  out
}  
.btransformation <- function(q,transformation)
{
  switch(transformation,
          "none" = q,
          "log" = exp(q),
          "logit"= exp(q)/(1+exp(q)),
          "Gompertz" =  1-exp(-exp(q)),
          "arcsin" = sin(q)^2
  )
}
.qx.stderr <- function(l, kernels, varqx){
  sum(kernels[,,l]^2*varqx)^(1/2)  
}

.ddbeta2 <- function(l,iandj,maxx,maxy,xweights,yweights,obsq,cv,exposures=NULL,adaptx,adapty,hx,hy,sx,sy){
  # iandj is a dataframe with all combinations of age and years
  # l points to the age and year to be estimated (index of iandj)
  # returns the esitmated value and the kernel matrix used to estimate it
  f  <- function(s,t){
    bx <- switch(adaptx, a=hx*xweights[iandj$i[l]]^sx, b= hx*xweights[iandj$i[s]]^sx, ab=hx*(xweights[iandj$i[s]]/xweights[iandj$i[l]])^sx)
    by <- switch(adapty, a=hy*yweights[iandj$j[l]]^sy, b= hy*yweights[iandj$j[t]]^sy, ab=hx*(yweights[iandj$j[t]]/yweights[iandj$j[l]])^sy)
    dbx <- .ddbeta(x[s],xmin=0,xmax=maxx-1,m=iandj$i[l]-1,h=bx)
    dby <- .ddbeta(y[t],xmin=0,xmax=maxy-1,m=iandj$j[l]-1,h=by)
    db  <- dbx*dby
    db  }

  x <- 0:(maxx-1)
  y <- 0:(maxy-1)
  K <- outer(x+1,y+1,f)
  if (cv) K[iandj$i[l],iandj$j[l]]  <- 0
  K <- K/sum(K)
  qxest <- sum(K*obsq)
  list(qxest=qxest,K=K) 
}

.dkCV <- function(par,hx,hy, sx, sy,cvhx, cvhy, cvsx, cvsy, xweights, yweights, maxx, maxy, obsq, cvres,iandj,adaptx,adapty) {
  if (cvhx) hx <- par$hx
  if (cvhy) hy <- par$hy
  if (cvsx) sx <- par$sx
  if (cvsy) sy <- par$sy

  qxest      <- matrix(0,(maxx),(maxy))             
  qxest      <- lapply(1:nrow(iandj),FUN=".ddbeta2",iandj=iandj,maxx=maxx,maxy=maxy,xweights=xweights,yweights=yweights,obsq=obsq, cv=TRUE,adaptx=adaptx,adapty=adapty,hx=hx,hy=hy,sx=sx,sy=sy)
  qxest      <- matrix(unlist(sapply(qxest,"[",1)),maxx,maxy)
  if (cvres=="propres") { CV = qxest/obsq-1 }
  else {CV = qxest-obsq}
  CV
}

.vcf <- function(n,p){
  sqrt((1-p)/(n*p))
}


.ddbeta <- function(x,xmin=0,xmax,m,h){
  y      <- (x-xmin+0.5)/(xmax-xmin+1)
  a      <- (m-xmin+0.5)/(h*(xmax-xmin+1))+1
  b      <- (xmax-m+0.5)/(h*(xmax-xmin+1))+1
  db     <- dbeta(x=y,shape1=a,shape2=b,ncp = 0,log = FALSE)
  dbdeno <- dbeta(.trans(xmin=xmin,xmax=xmax),shape1=a,shape2=b,ncp = 0,log = FALSE)
  db/sum(dbdeno)
}

as.data.frame.dbkGrad <- function(x, row.names = x$limx[1]:x$limx[2], optional = FALSE, ...)
{
  l <- list(obsq=x$obsq, fitted.values=x$fitted.values)
  if(!is.null(x$exposures)){
    l <- c(l,list(exposures=x$exposures, lowerbound=x$lowerbound, upperbound=x$upperbound, bonferroniUB=x$bonferroniUB,bonferroniLB=x$bonferroniLB,sidakUB=x$sidakUB,sidakLB=x$sidakLB))
  }
  df <- as.data.frame(l,row.names=row.names, optional=optional, ...)
  df
}
print.dbkGrad <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nGraduated Rates:\n")
  print(x$fitted.values)
}

residuals.dbkGrad <- function(object,type=c("working","proportional","response","deviance", "pearson"),...){
  type <- match.arg(type)
  o <- object$obsq
  f <- object$fitted.values
  e <- object$exposures
  
  res <- switch(type,
                working = o - f,
                proportional= (o/f)-1,
                response= (o - f)*e,
                deviance=sign(o - f) * sqrt(2*e*o*log(ifelse(o == 0, 1,o/f)) +2*e*(1-o)* log(ifelse(o==1,0,(1-o)/(1-f)))),
                pearson =(o - f)*e/sqrt(e*f*(1-f))
  )
  as.matrix(res)
}

.extract<- function(mat,rows, columns,byage){
  y <- mat[rows,columns]
  if (byage)  {
    if (is.vector(y))  y <- as.matrix(y) 
    if (length(rows) ==1) y <- t(y)
    rownames(y) <- rownames(mat)[rows]
    colnames(y) <- colnames(mat)[columns]
  }else {
    if (is.vector(y))  y <- as.matrix(y) else y <- t(y)
    if (length(columns) ==1) y <- t(y)
    rownames(y) <- colnames(mat)[columns]
    colnames(y) <- rownames(mat)[rows]
  }
  as.matrix(y)
}

plot.dbkGrad <- function(x, plottype=c("obsfit","fitted", "observed", "exposure","residuals","checksd"),plotstyle=c("mat", "level","persp"),restype=c("working","proportional","response","deviance", "pearson"),
                         byage=TRUE, columns, rows, CI=TRUE, CBBonf=FALSE, CBSidak=FALSE, logscale=TRUE,alphares=0.05,col,...)
{
  def.par <- par(no.readonly = TRUE) 
  if (missing(columns)) columns <- 1:ncol(x$fitted.values)
  if (missing(rows)) rows <- 1:nrow(x$fitted.values)
  plottype <- match.arg(plottype)
  restype <- match.arg(restype)
  if (missing(plotstyle)) plotstyle <- ifelse(length(rows)>1 && length(columns)>1,"level","mat")
  else plotstyle <-  match.arg(plotstyle)
  if (!plottype=="checksd") plottype <- paste0(plottype,".",plotstyle)
  rows <- sort(rows)
  columns <- sort(columns)
  eval(parse(text=paste0(".",plottype,".plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,alphares=alphares,restype=restype,col=col,...)")))
  if(any(plottype==c("obsfit.mat","fitted.mat") & !is.null(x$exposures))){
    if (CI).ci.plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,type="CI",col=col)
    if (CBBonf).ci.plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,type="CBBonf",col=col)
    if (CBSidak).ci.plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,type="CBSidak",col=col)
  }
  par(def.par)
}
.level.plot<- function(grid,formula,key.labels,logscale,xlab="year",ylab="age",col=1,pch,lwd,legend,
                       palette="terrain.colors",at= quantile(grid$z, prob = seq(0,1, length.out=11)),
                       digits=5,sub="",...){
  if (missing(key.labels)){
    if(logscale) key.labels<-round(exp(at),digits=digits)
    else key.labels<-round(at,digits=digits) 
  } 
  col.regions  <- get(palette)(length(at))
  print(levelplot(x=formula, grid, xlab=xlab, ylab=ylab, sub=sub, col.regions=col.regions, colorkey = list(col = col.regions, labels=as.character(key.labels)), at = at,...))
}

.obsfit.level.plot <-   function(x,rows,columns,byage,logscale,alphares,restype,col,...){
  y <- .extract(x$obsq,rows,columns, byage=T)
  xx <- as.numeric(rownames(y))
  yy <- as.numeric(colnames(y))
  grid <- expand.grid(list(x = xx, y = yy, type = c("Actual", "Fitted")))
  grid$z <- log(c(as.vector(y),as.vector(.extract(x$fitted.values,rows,columns, byage=T))))
  .level.plot(grid=grid,formula=z ~ y * x | type,layout = c(2, 1),logscale=TRUE,col=col,...)
}
.observed.level.plot <- function(x,rows,columns,byage,logscale,alphares,restype,col,...){
  y <- .extract(x$obsq,rows,columns, byage=TRUE)
  xx <- as.numeric(rownames(y))
  yy <- as.numeric(colnames(y))
  grid <- expand.grid(list(x = xx, y = yy))
  grid$z <- log(as.vector(y))
  .level.plot(grid,formula=z ~ y * x,logscale=TRUE,col=col,...)
}
.fitted.level.plot <-   function(x,rows,columns,byage,logscale,alphares,restype,col,...){
  y <- .extract(x$fitted.values,rows,columns, byage=TRUE)
  xx <- as.numeric(rownames(y))
  yy <- as.numeric(colnames(y))
  grid <- expand.grid(list(x = xx, y = yy))
  grid$z <- log(as.vector(y))
  .level.plot(grid,formula=z ~ y * x,logscale=TRUE,col=col,...)
}
.exposure.level.plot <- function(x,rows,columns,byage,logscale,alphares,restype,at,col,...){
  y <- .extract(x$exposure,rows,columns, byage=TRUE)
  xx <- as.numeric(rownames(y))
  yy <- as.numeric(colnames(y))
  grid <- expand.grid(list(x = xx, y = yy))
  grid$z <- as.vector(y)
  .level.plot(grid=grid,formula=z ~ y * x,logscale=FALSE,digits=0,col=col,...)
}
.residuals.level.plot <-function(x,rows,columns,byage,logscale,alphares,restype,ylab=paste(restype, "residuals"),sub,col,...)
{
  res <- .extract(mat=residuals(x, type=restype) ,rows, columns,byage=TRUE)  
  rows <- 1:nrow(res)
  columns <- 1:ncol(res)
  xx <- as.numeric(rownames(res))
  yy <- as.numeric(colnames(res))
  grid <- expand.grid(list(x = xx, y = yy))
  grid$z <- as.vector(res)
  if (missing(sub))  sub <-  paste(restype, "residuals")
  def.par <- par(no.readonly = TRUE)
  .level.plot(grid,formula=z ~ y * x,logscale=FALSE,sub=sub,col=col,...)
  readline("Press <Enter> to continue")
  par(def.par) 
  .residuals.density.plot(res,xlab=sub,logscale=FALSE,sub=sub,col=col)
}
.residuals.density.plot <- function(res,rows,columns,byage,logscale,alphares,restype,xlab,main,sub=sub,col,...)
{
  if (missing(col)) col <- 1
  plot(density(res),main="",xlab=xlab,col=col,...)
}

.residuals.mat.plot <- function(x,rows,columns,byage,logscale,alphares,restype,xlab,ylab=paste(restype, "residuals"),col,...)    
{
  def.par <- par(no.readonly = TRUE) 
  res <- .extract(mat=residuals(x, type=restype) ,rows, columns,byage)
  for (i in 1:ncol(res) ){
    .mat.plot(mat= res,logscale=FALSE,ylab=ylab,x=x,byage=TRUE,rows=1:nrow(res), columns=i,type="h",col=col,...)
    abline(qnorm(alphares/2),0, col="gray", lty=2)
    abline(0,0, col="gray", lty=2)
    abline(qnorm(1-alphares/2),0, col="gray", lty=2)
    if (i< ncol(res)) readline("Press <Enter> to continue")
  }
  readline("Press <Enter> to continue")
  par(def.par) 
  .residualsVsFitted.mat.plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,restype=restype,alphares=alphares,col=col,...)
  par(def.par) 
  readline("Press <Enter> to continue")
  .residuals.density.plot(res,xlab=ylab,col=col,...)
}

.residualsVsFitted.mat.plot <- function(x,rows,columns,byage,logscale,alphares,restype,xlab="fitted values",ylab=paste(restype, "residuals"),col,mar=c(5, 4, 4, 1) + 0.1,cex.axis=0.8,cex=0.8,...)
{
  res    <- as.vector(.extract(mat=residuals(x, type=restype),rows, columns,byage))
  fitted <- as.vector(.extract(mat=x$fitted.values,rows, columns,byage))
  atx <- xlabels <-  round(fitted,digits=5)
  if (missing(col)) col <- 1
  if (logscale)  {
      fitted <- log(fitted)
      xlab <- paste(xlab,"(log scale)")
      atx <- log(atx)
  }
  aty <- round(seq(min(res),max(res),length.out =5),digits=5)
  layout(cbind(1,1))
  par(mar=mar,cex.axis=cex.axis,cex=cex)
  matplot(y=res,x=fitted,xlab=xlab,ylab=ylab,type="h",axes = FALSE,col=col,...) #
  axis(1,at = atx,labels = xlabels)
  axis(2, at = aty, labels =aty )
  abline(qnorm(alphares/2),0, col="gray", lty=2)
  abline(0,0, col="gray", lty=2)
  abline(qnorm(1-alphares/2),0, col="gray", lty=2)
}

.checksd.plot <- function(x,rows,columns,byage,logscale,alphares,restype,col,...){
  if (missing(col)) col <- 1
  y <- .extract(mat=residuals(x, type=restype),rows, columns,byage)
  for (i in 1:ncol(y) ){
    acf(y[,i],main="",sub=paste(restype,"residuals",colnames(y)[i]),col=col,...)
    readline("Press <Enter> to continue")
    ADF(y[,i],main="",sub=paste(restype,"residuals",colnames(y)[i]),col=col,...)
    if (i< ncol(y)) readline("Press <Enter> to continue")
  }
}



.ci.plot <-function(x,logscale,rows, columns, byage,col,type,...)
{
  f <- function(d,x0,x1,y0,y1,col,type,lty) {
    segments(x0=x0, y0=y0[,d],x1=x0,y1=y1[,d],col=col[d],lty,...)
    segments(x0=x0-l, y0=y0[,d],x1=x0+l,y1=y0[,d],col=col[d],...)
    segments(x0=x0-l, y0=y1[,d],x1=x0+l,y1=y1[,d],col=col[d],...)
  }
  
  if (type=="CI") {
    y0 <- x$lowerbound
    y1 <- x$upperbound
    lty <- "solid"
  } 
  else if (type=="CBBonf") {
    y0 <- x$bonferroniLB
    y1 <- x$bonferroniUB
    lty <- "dotted"
  } 
  else {
    y0 <- x$sidakLB
    y1 <- x$sidakUB
    lty <- "dotted"
  }
  rownames(y0) <- rownames(x$fitted)
  colnames(y0) <- colnames(x$fitted)
  
  y0 <- .extract(y0,rows,columns,byage=byage)
  y1 <- .extract(y1,rows,columns,byage=byage)
  
  if(logscale) {
    y0 <- log(y0)
    y1 <- log(y1)
  }

  l <-0.05
  if (missing(col)) col <- 1:length(columns)
  x0 <-as.numeric(rownames(y0))
  foo <-lapply(1:ncol(y0), FUN=f,x0=x0,x1=x0, y0=y0,y1=y1,col=col,lty=lty,...)
}

.exposure.mat.plot<- function(x,rows,columns,byage,logscale,alphares,restype,ylab="exposures",col,pch="",lwd=1.5,...)
{
  y <- .extract(mat=x$exposures,rows, columns,byage)
  for (i in 1:ncol(y) ){
    yat <- seq(0,max(y[,i]),length=5) 
    ylabels <- round(yat,digits=0)
    .mat.plot(mat=y,logscale=FALSE,ylab=ylab,col=col,x=x,byage=TRUE,rows=1:nrow(y),columns=i,type="h",pch=pch,lwd=lwd,ylabels=ylabels,yat=yat,ylim=range(yat),...)
    if (i< ncol(y)) readline("Press <Enter> to continue")
  }
}
  
.obsfit.mat.plot<- function(x,rows,columns,byage,logscale,alphares,restype,xlab,ylab,col,pch,lwd=1,legend,type="p",mar=c(5, 4, 4, 1) + 0.1,cex.axis=0.8,cex=0.8,...)
{
  if(missing(ylab)) ylab <- "mortality rates"
  if (!byage) {
    y <-(cbind(t(x$fitted.values),t(x$obsq))[columns,c(rows,nrow(x$fitted.values)+rows)])
    if (is.vector(y)) y <- t(y)
    rownames(y) <-colnames(x$fitted.values)[columns]
    if(missing(xlab)) xlab <- "Years"
  } else 
  {
    y <-(cbind(x$fitted.values,x$obsq)[rows,c(columns,ncol(x$fitted.values)+columns)])
    if (is.vector(y)) y <- t(y)
    if(missing(xlab)) xlab <- "age" 
    rownames(y) <-rownames(x$fitted.values)[rows]
  }
  legend <- paste(colnames(y),c(rep("fitted",(ncol(y)/2)),rep("observed",(ncol(y)/2))))
  columns <-1:(length(columns)*2)
  rows <-1:(length(rows))
  if(missing(pch)) pch <- c(rep(16,(ncol(y)/2)),rep(1,(ncol(y)/2)))
  if(missing(col)) col <- c(1:(ncol(y)/2),1:(ncol(y)/2))
  if(logscale) {
    y <- log(y)
    yat <-ylabels <- seq(min(y),max(y),length.out =5)
    ylabels <- exp(yat)
    ylab <- paste(ylab,"(log scale)")
  }
  else yat <-ylabels <- seq(min(y),max(y),length.out =5)

  if(missing(legend)) legend <- colnames(y)
  
  layout(cbind(2,1), widths=c(6,1))
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("left", legend=legend, col=col, pch=pch, h=FALSE,cex=0.8, bty="n",inset=c(0,0)) 
  par(mar=mar,cex.axis=cex.axis,cex=cex)
  matplot(y=y,x=as.numeric(rownames(y)),pch=pch,xlab=xlab,ylab=ylab, col=col,type=type,axes = FALSE,...) #
  axis(1,at = as.numeric(rownames(y)),labels = round(as.numeric(rownames(y))))
  axis(2,at = yat, labels = round(ylabels,digits=5))
}

.dupl <- function(x) unlist(lapply(x,FUN=rep,2))


.fitted.mat.plot<- function(x,rows,columns,byage,logscale,alphares,restype,ylab="graduated mortality rates",col,...){
  .mat.plot(mat=x$fitted.values,logscale=logscale,ylab=ylab,x=x,byage=byage,rows=rows,columns=columns,col=col, ...)
}
.observed.mat.plot<- function(x,rows,columns,byage,logscale,alphares,restype,ylab="observed mortality rates",col,...)
{
  .mat.plot(mat=x$obsq,logscale=logscale,ylab=ylab,x=x,byage=byage,rows=rows,columns=columns,col=col,...)
}

.mat.plot<- function(mat,x,rows,columns,byage,logscale,alphares,restype,xlab,ylab,
                     col,lwd=1,pch=16,legend,at,by,type,xlabels,ylabels,yat,mar=c(5, 4, 4, 1) + 0.1,cex.axis=0.8,cex=0.8,...)
{
  y <- mat[rows,columns]
  if (length(rows)==1) byage <- FALSE
  if (byage)  {
    if (is.vector(y))  y <- as.matrix(y) 
    if (length(rows) ==1) y <- t(y)
    rownames(y) <- rownames(mat)[rows]
    colnames(y) <- colnames(mat)[columns]
    if(missing(xlab)) xlab <- "Age" 
    if (missing(type)) type <- "b"
  } else {
    if (is.vector(y))  y <- as.matrix(y) else y <- t(y)
    if (length(columns) ==1) y <- t(y)
    rownames(y) <- colnames(mat)[columns]
    colnames(y) <- rownames(mat)[rows]
    if(missing(xlab)) xlab <-  "Years" 
    if (missing(type)) type <- "b"
  }
  if (missing(col)) col <- 1:ncol(y)
  if(logscale) {
    y <- log(y)
    if (missing(yat)) yat <- seq(min(y),max(y),length.out =5)
    if (missing(ylabels)) ylabels <- round(exp(yat),digits=5) 
    ylab <- paste(ylab,"(log scale)")
  }
  else
  {
    if (missing(yat)) yat <- seq(min(y),max(y),length.out =5)
    if (missing(ylabels)) ylabels <- round(yat,digits=5)
  }
  layout(cbind(2,1), widths=c(6,1))# widths=c(1,3)
  par(mar=c(0, 0, 0, 0))
  plot.new()
  if(missing(legend)) legend <- colnames(y)
  if (missing(xlabels)) xlabels <- round(as.numeric(rownames(y)),digits=3)
  legend("left", legend=legend, col=col, pch=pch, h=FALSE,cex=0.8, bty="n",inset=c(0,0)) 
  par(mar=mar,cex.axis=cex.axis,cex=cex)
  matplot(y=y,x=as.numeric(rownames(y)),pch=pch,xlab=xlab,ylab=ylab, col=col,type=type,axes = FALSE,lwd=lwd,...) #
  axis(1,at = as.numeric(rownames(y)),labels = xlabels)
  axis(2,at = yat, labels = ylabels)
}   
.observed.persp.plot<- function(x,rows,columns,byage,logscale,alphares,restype,zlab="observed mortality rates",col,...){
  z <- .extract(x$obs,rows,columns,TRUE)
  .persp.plot(z=z,logscale=logscale,zlab=zlab,col=col,...)
}
.fitted.persp.plot<- function(x,rows,columns,byage,logscale,alphares,restype,zlab="fitted mortality rates",col,...){
  z <- .extract(x$fitted.values,rows,columns,TRUE)
  .persp.plot(z=z,logscale=logscale,zlab=zlab,col=col,...)
}
.exposure.persp.plot<- function(x,rows,columns,byage,logscale,alphares,restype,zlab,col,...){
  z <- .extract(x$exposure,rows,columns,TRUE)
  if(missing(zlab)) zlab <- "esposures"
 # zat <- seq(0,max(z),length=5) 
#  zlabels <- round(zat,digits=0)
  .persp.plot(z=z,logscale=FALSE,rows=rows,columns=columns,zlab=zlab,segments=TRUE,zlim=c(0,max(z)),col=col,...)
}

.residuals.persp.plot <- function(x,rows,columns,byage,logscale,alphares,restype,xlab,zlab=paste(restype, "residuals"),col,...)    
{
  def.par <- par(no.readonly = TRUE) 
  z  <- .extract(mat=residuals(x, type=restype) ,rows, columns,TRUE)  
  .persp.plot(z=z,logscale=FALSE,zlab=zlab,segments=TRUE,col=col,...)
  readline("Press <Enter> to continue")
  par(def.par) 
  .residuals.density.plot(res=z,xlab=zlab,rows=rows,columns=columns,byage=byage,logscale=logscale,alphares=alphares,restype=restype,col=col,...)
  readline("Press <Enter> to continue")
  par(def.par) 
  .residualsVsFitted.persp.plot(x=x,rows=rows,columns=columns,zlab=zlab,restype=restype,byage=byage,logscale=logscale,segments=TRUE,col=col,...)
}


.residualsVsFitted.persp.plot <- function(x,rows,columns,byage,logscale,alphares,restype,
                                          xlab="fitted values",ylab=ifelse(byage,"age","year"),zlab=paste(restype, "residuals")
                                          ,col,...)    
{
  z <- .extract(mat=residuals(x, type=restype),rows, columns,!byage)
  x <- .extract(mat=x$fitted.values,rows, columns,!byage)
  if (logscale) x <- log(x)
  for (i in 1:ncol(x)){
    ord <- order(x[,i])
    z[,i] <- z[ord,i]
    x[,i] <- x[ord,i]
  }
  .persp.plot(z=z,x=x,byage=byage,logscale=FALSE,xlab=xlab,ylab=ylab,zlab=zlab,col=col,...)
}  
.obsfit.persp.plot <- function(x,rows,columns,byage,logscale,alphares,restype,col,...)
{
  .observed.persp.plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,alphares=alphares,restype=restype,col=col,...)
  readline("Press <Enter> to continue")
  .fitted.persp.plot(x=x,rows=rows,columns=columns,byage=byage,logscale=logscale,alphares=alphares,restype=restype,col=col,...)
}
.persp.plot<- function(x=matrix(rep(as.numeric(rownames(z)),ncol(z)),nrow(z),ncol(z)),
                       y= matrix(rep(as.numeric(colnames(z)),nrow(z)),nrow(z),ncol(z),byrow=TRUE),
                       z,rows,columns,byage,logscale,alphares,restype,
                       theta=-30,col,pch=20,cex=0.6,lwd=2,
                       zlab="mortality rates",ylab="years",xlab ="age",legend,axes=TRUE, box=TRUE,border=FALSE,
                       segments=FALSE,xlim=range(x),ylim=range(y),zlim=range(z),ticktype="detailed",mar=c(5, 4, 4, 1) + 0.1,cex.axis=0.8,...){
  if(logscale) {
    z <- log(z)
    zlab <- paste(zlab,"(log scale)")
  }
  if (missing(col)) col <- 1:ncol(z)
  else col <- rep(col, ncol(z))
  par(mar=mar,cex.axis=cex.axis,cex=cex)
 # plot.new()
  persp(x=x[,1], y=y[1,], z=z,
        xlab = xlab,
        ylab = ylab,
        zlab = zlab,
        xlim = xlim,
        ylim=ylim,
        zlim=zlim,
        border = border,
        theta = theta, 
        axes = axes, 
        ticktype = ticktype,
        box = box,
        ...
  ) -> res
  
  for(i in 1:ncol(y)){
    p0 <- trans3d(x = x[,i],y =y[,i],z = z[,i], pmat = res)
    if (!segments) points(p0, col = col[i], lwd=lwd, cex=cex, pch=pch)
    else 
    {
      p1 <- trans3d(x = x[,i],y =y[,i],z = 0, pmat = res)  
      segments(x0=p0$x,y0=p0$y,x1=p1$x,y1=p1$y, col = col[i], lwd=lwd, cex=cex, pch=pch)
      s0 <- trans3d(x = min(x),y =y[,i],z = 0, pmat = res) 
      s1 <- trans3d(x = max(x),y =y[,i],z = 0, pmat = res)  
      segments(x0=s0$x,y0=s0$y,x1=s1$x,y1=s1$y, col = col[i], lwd=1, cex=cex, pch=pch, lty=3)
    }
  }
}  
