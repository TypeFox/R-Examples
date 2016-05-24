sprmsCV <-
function(formula, data, as, etas, nfold=10, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999, center = "median", scale = "qn", plot=TRUE, numit=100, prec=0.01, alpha=0.15){
  
#  library(cvTools)
#  library(reshape)
#  library(ggplot2)
#  source("sprms.r")
#  source("predict.sprm.r")
  MSPE <- number <- NULL
  if (missing(as)){
    stop("Specify vector as")
  }
  if (missing(etas)){
    etas <- seq(0,0.9,0.1) 
  }
  
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
    stop("The values in vector as have to be positive.")
  }
  if(any(etas<0)|any(etas>=1)){
    stop("Parameter etas has to be a vector with values coming from the intervall [0,1)")
  }
  if(!any(fun == c("Hampel", "Huber", "Fair"))){
    stop("Invalid psi function. Choose Hampel, Huber or Fair for parameter fun.")
  }
  if(probp1>1|probp1<=0){
    stop("Parameter probp1 is a probability. Choose a value between 0 and 1")
  }
  if(fun=="Hampel"){
    if (!(probp1<hampelp2 & hampelp2<hampelp3 & hampelp3<=1)){
      stop("Wrong choise of parameters for Hampel function. Use 0<probp1<hampelp2<hampelp3<=1")
    }
  }
  if(plot==TRUE & length(as)==1){
    warning("Plot option is only availagle for multiple values of as.")
    plot <- FALSE
  }

  folds <- cvFolds(n, K = nfold, R = 1, type = "random")
  spe <- array(NA, dim = c(n, length(as), length(etas)))
  dimnames(spe) <- list(1:n, 1:length(as), 1:length(etas))
  nzcoef <- array(NA, dim = c(nfold, length(as), length(etas)))
  dimnames(nzcoef) <- list(1:nfold, 1:length(as), 1:length(etas))
  for (f in 1:nfold){
    dtrain <- data[folds$which!=f,]
    dtest <- data[folds$which==f,]
    for (i in 1:length(as)){
      a <- as[i]
      for (j in 1:length(etas)){
        e <- etas[j]
        res <- tryCatch({     
          trainmod <- sprms(formula,data=dtrain, a=a, eta=e, fun, probp1,
                            hampelp2, hampelp3,  center, scale, print=FALSE, numit, prec)
          }, error=function(err){
            error <- TRUE
            return(error)
        })          
        if (is.logical(res)){
          print(paste("CV broke off for a=", a,"and eta=", e))
        } else{
          trainmod <- res
          spe[folds$which==f,i,j] <- (predict.sprm(trainmod, newdata=dtest) - dtest[,1])^2
          nzcoef[f,i,j] <- sum(trainmod$coefficients!=0)
        }
      } 
    }  
  }  
  
  mspe <- apply(spe, c(2,3), function(x) mean(sort(x)[1:(length(x)*(1-alpha))]))
  sdspe <- apply(spe, c(2,3), sd)
  mnzcoef <- apply(nzcoef, c(2,3), mean)

  optind <- which(mspe==min(mspe, na.rm=TRUE), arr.ind=TRUE)[1,]
  
  eta <- etas[optind[2]]
  a <- as[optind[1]]
  
  if (plot==TRUE){
    print(paste("optimal model: a =", a, "eta =", eta))
    mycol1 <- colorRampPalette(c("blue2","red", "orange", "yellow"))(64)
    mycol2 <- colorRampPalette(c("darkgreen", "yellow"))(64)    
    ggmspe <- mspe
    rownames(ggmspe) <- as
    colnames(ggmspe) <- etas
    ggmspe <- melt(ggmspe)
    names(ggmspe) <- c("a", "eta", "MSPE")
    mspeplot <- ggplot(ggmspe, aes(x=as.factor(a), y=as.factor(eta), fill=MSPE)) +
      geom_tile() +
      scale_fill_gradientn(colours=mycol1) +
      ggtitle(paste0("Trimmed Mean SPE (minimum at a= ", a, ", eta= ", eta, ")")) +
      xlab("a") + 
      ylab("eta")
      
    
    ggmnzcoef <- mnzcoef
    rownames(ggmnzcoef) <- as
    colnames(ggmnzcoef) <- etas
    ggmnzcoef <- melt(ggmnzcoef)
    names(ggmnzcoef) <- c("a", "eta", "number")
    mnzcoefplot <- ggplot(ggmnzcoef, aes(x=as.factor(a), y=as.factor(eta), fill=number))+
      geom_tile() +
      scale_fill_gradientn(colours=mycol2)+
      ggtitle("# nonzero coefficients") +
      xlab("a") + 
      ylab("eta")
    
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2,1)))
    print(mspeplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
    print(mnzcoefplot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
  }
  
  prmspFit <- sprms(formula=formula, data=data, a=a, eta=eta, fun, probp1,
                    hampelp2, hampelp3,  center, scale, print=FALSE, numit, prec)
  
  return(list(opt.mod=prmspFit, spe=spe, nzcoef=nzcoef))
}
