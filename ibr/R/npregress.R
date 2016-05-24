npregress <- function(x,y,criterion="rmse",bandwidth=NULL,kernel="g",control.par=list(),cv.options=list()) {
  kern <- c("g","q","e","u")
  kernel <- match.arg(kernel,kern)
  if (any(is.na(x))) stop("NA's in x\n")
  if (any(is.na(y))) stop("NA's in y\n")
  if (!is.numeric(x)&(is.data.frame(x))) {
    x <- x[,1]
    if (!is.numeric(x)) stop("first column of data-frame is not numeric\n")
  }
  if (is.matrix(x)) {
    x <- as.vector(x)
  }
  if (!is.numeric(x)) stop("x must be a numeric vector (or a data-frame with first column of numeric type)\n")
  if (!is.numeric(y)) stop("y must be numeric\n")
  if (! is.vector(y)) {
    y <- as.vector(y)
  }
  if (missing(bandwidth)) bandwidth <- NULL
  if ((!is.null(bandwidth))&(!is.numeric(bandwidth) || (bandwidth<0))) stop("invalid bandwidth\n")
  contr.sp <- list(bandwidth=NULL,degree=0)
  contr.sp[(names(control.par))] <- control.par
  if (is.null(bandwidth)) {
    if ((!is.null(contr.sp$bandwidth))&(!is.numeric(contr.sp$bandwidth) || (contr.sp$bandwidth<0))) stop("invalid bandwidth\n")
    bandwidth <- contr.sp$bandwidth
  }
  if (contr.sp$degree>1) stop("Not implemented. Please consider using KernSmooth or another library for degree greater or equal to 2\n")
  if (contr.sp$degree==0) nom <- "reg" else nom <- "regpol"
  crit <-c("rmse","map")
  criterion <- match.arg(criterion,crit)
  n <- length(x)
  if ((criterion%in%c("rmse","map"))&(is.null(bandwidth))) {
    mini <- 1/n
    kmax=floor(log(n*diff(range(x))/3)/log(1+1/n))
    gridbw <- 1/n*(1+1/n)^(0:kmax)
    cv <- list(gridbw=gridbw,ntest=1,ntrain=NULL,Kfold=TRUE,type="consecutive",seed=NULL,npermut=NULL)
    cv[(names(cv.options))] <- cv.options
    if ((!is.numeric(cv$gridbw))||(!is.vector(cv$gridbw))) stop("invalid gridbw component of cv.options: must be a numeric vector\n")
    if (!all(sapply(cv[c(2,3,6,8)], FUN=function(x) is.numeric(x)||is.null(x)))) stop("invalid cv parameters: must be numeric or NULL\n")
    if (any(names(cv.options)=="ntrain")) cv$ntest <- NULL
    sel <- cvobs(n,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed)
    ordre <- unlist(sel)
    xord <- x[ordre]
    yord <- y[ordre]
    nj <- unlist(lapply(sel,length))
    effold <- c(0,cumsum(nj))
    neffold <- length(sel)
    nom1 <- paste(nom,kernel,"cv",sep="")
    prov <- .C(nom1,as.double(xord),as.integer(length(xord)),as.double(yord),as.double(gridbw),as.integer(length(gridbw)),as.integer(effold),as.integer(neffold),double(length(gridbw)),double(length(gridbw)))
    rmse <- sqrt(prov[[8]]/sum(n-nj))
    map <- prov[[9]]/sum(n-nj)
    choixbw <- list(gridbw=gridbw,rmse=rmse,map=map)
    bandwidth <- gridbw[switch(criterion,rmse=which.min(rmse),map=which.min(map))]
  } else {
    choixbw <- NULL
    criterion <- "user"
  }
  nom2 <- paste(nom,kernel,sep="")
  if (contr.sp$degree==1) {
    prov <- .C(nom2,as.double(x),as.integer(n),as.double(y),as.double(bandwidth),as.double(x),as.integer(n),double(n),double(1),double(n))
   } else {
   prov <- .C(nom2,as.double(x),as.integer(n),as.double(y),as.double(bandwidth),as.double(x),as.integer(n),double(n),double(1),double(n))
  }
  fit <- prov[[7]]
  df <- prov[[8]]
  residuals <- y- fit
  
  res <- list(bandwidth=bandwidth,residuals=residuals,fitted=fit,df=df,call=list(x=x,y=y,criterion=criterion,kernel=kernel,degree=contr.sp$degree),criteria=choixbw)
  class(res) <- c("npregress", "list")
  return(res)
}

     
 
