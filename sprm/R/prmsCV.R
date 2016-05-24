prmsCV <-
function(formula, data, as, nfold=10, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999, center = "median", scale = "qn", usesvd=FALSE, plot=TRUE, numit=100, prec=0.01, alpha=0.15){
  
#  library(cvTools)
#  library(fields)
#  source("prms.r")
#  source("predict.prm.r")
 
  if (missing(as)){
    stop("Specify vector as")
  }
  if(!class(formula)=="formula"){formula <- formula(formula)}  
  if(is.data.frame(data) | is.list(data)){
    mt <- terms(formula, data=data)
    yname <- dimnames(attr(mt,"factors"))[[1]][1]
    if(is.list(data)){
      datnames <- names(data)
    } else {
      datnames <- colnames(data)
    }
    ic <- attr(mt, "intercept")
    if (ic==0){
      data <- tryCatch({data <- cbind(data[[which(datnames==yname)]], model.matrix(mt, data))},
                       error=function(err){
                         error <- TRUE
                         return(error)
                       }) 
    } else{
      data <- tryCatch({data <- cbind(data[[which(datnames==yname)]],model.matrix(mt, data)[,-1])},
                       error=function(err){
                         error <- TRUE
                         return(error)
                       }) 
    }
    if (is.logical(data)){
      stop("Data cannot be matched with formula.")
    } else {
      colnames(data)[1] <- dimnames(attr(mt,"factors"))[[1]][1]
      data <- as.data.frame(data)
    }    
  } else {
    stop("Wrong data fromat.")
  }

  n <- dim(data)[1]
  p <- dim(data)[2]-1
  
  if(any(as>(n-n/nfold))|any(as>p)){
    stop("Maximal number of components is too large.")
  }
  if (any(as<=0)){
    stop("The number of components has to be positive.")
  }
  if(!any(fun == c("Hampel", "Huber", "Fair"))){
    stop("Invalid weighting function. Choose Hampel, Huber or Fair for parameter fun.")
  }
  if(probp1>1|probp1<=0){
    stop("Parameter probp1 is a probability. Choose a value between 0 and 1")
  }
  if(fun=="Hampel"){
    if (!(probp1<hampelp2 & hampelp2<hampelp3 & hampelp3<=1)){
      stop("Wrong choise of parameters for Hampel function. Use 0<probp1<hampelp2<hampelp3<=1")
    }
  }

  folds <- cvFolds(n, K = nfold, R = 1, type = "random")
  spe <- matrix(nrow=n, ncol=length(as))
  for (i in c(1:length(as))){
    a <- as[i]
    for (f in 1:nfold){
      dtrain <- data[folds$which!=f,]
      dtest <- data[folds$which==f,]
      trainmod <- prms(formula,data=dtrain, a=a, fun, probp1,
                       hampelp2, hampelp3,  center=center, scale=scale, usesvd=usesvd, numit, prec)
      spe[folds$which==f,i] <- (as.vector(predict.prm(trainmod, newdata=dtest)) - dtest[,1])^2
    }
  }
  
  mspe <- apply(spe, 2, function(x) mean(sort(x)[1:(length(x)*(1-alpha))]))
  mrspe <- apply(spe,2, function(x) mean(sqrt(sort(x)[1:(length(x)*(1-alpha))])))
  
  optind <- which.min(mspe)[1]  

  a <- as[optind]
  
  prmspFit <- prms(formula=formula,data=data, a=a, fun, probp1,
                    hampelp2, hampelp3,  center, scale, usesvd, numit, prec)
  
  if (plot==TRUE){
    
		q75 <- apply(spe, 2, function(x) quantile(x, 3/4))
		q25 <- apply(spe, 2, function(x) quantile(x, 1/4))
		mpe <- mspe
		ylab <- "squared prediction error"
	  
    plotdat <- data.frame(as=as, mspe=mpe, q25=rep(q25,length(as)), q75=rep(q75,length(as)))
    plotmspe <- ggplot(plotdat, aes(x=as, y=mspe)) + geom_line() + geom_point() + 
    geom_point(data=subset(plotdat,as==a), aes(x=as, y=mspe), color="red") +
    geom_errorbar(aes(x=as, ymin=q25, ymax=q75), width=0.5, linetype=2) + 
    scale_x_continuous(breaks=as)+
	ylab(ylab)
    print(plotmspe)
  }
  return(list(opt.mod=prmspFit, spe=spe))
}
