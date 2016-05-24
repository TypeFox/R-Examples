prmdaCV <-
  function (formula, data, as, nfold=10, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999, probp4=0.01, yweights=TRUE, class=c("regfit", "lda"), prior=c(0.5,0.5), center = "median", scale = "qn", plot=TRUE, numit=100, prec=0.01) {
    ## 14/10/09 IH
    ## Cross validation for partial robust M-regression for discriminanat analysis
    ##
    ## uses prmsDA, predict.prmda, weig
    ##
    ## formula .... an object of class formula.
    ## data ....... a data frame or list which contains the variables given in formula.
    ## as ......... a vector with positive integers, which are the number of PRMS components to be estimated in the models.
    ## nfold ...... the number of folds used for cross validation, default is nford=10 for 10-fold CV.
    ## fun ........ an internal weighting function for case weights. Choices are "Hampel" (preferred), "Huber" or "Fair".
    ## probp1 ..... the 1-alpha value at which to set the first outlier cutoff for the weighting function.
    ## hampelp2 ... the 1-alpha values for second cutoff. Only applies to fun="Hampel".
    ## hampelp3 ... the 1-alpha values for third cutoff. Only applies to fun="Hampel".
    ## yweights ... logical; if TRUE y weights are calculated.
    ## class ...... type of classification; choices are "regfit" or "lda".
    ## prior ...... vector of length 2 with proir porbabilities of the groups; only used if class="lda".
    ## center ..... type of centering of the data in form of a string that matches an R function, e.g. "mean" or "median".
    ## scale ...... type of scaling for the data in form of a string that matches an R function, e.g. "sd" or "qn" or alternatively "no" for no scaling.
    ## plot ....... logical, default is TRUE. If TRUE a plot is generated with a measure of the prediction accuracy for each model (see Details).
    ## numit ...... the number of maximal iterations for the convergence of the coefficient estimates.
    ## prec ....... a value for the precision of estimation of the coefficients.
    
    # require(vegan)
    # require(penalizedLDA)
    
    if(!class(formula)=="formula"){formula <- formula(formula)}
    if(is.data.frame(data) | is.list(data)){
      mt <- terms(formula, data=data)
      yname <- dimnames(attr(mt,"factors"))[[1]][1]
      ic <- attr(mt, "intercept")
      if (ic==0){
        data <- tryCatch({data <- cbind(data[,which(colnames(data)==yname)], model.matrix(mt, data))},
                         error=function(err){
                           error <- TRUE
                           return(error)
                         }) 
      } else{
        data <- tryCatch({data <- cbind(data[,which(colnames(data)==yname)],model.matrix(mt, data)[,-1])},
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
    
    yorig <- data[,1]
    
    if (sum(unique(data[,1])%in%c(-1,1))<2){
      yfac <- as.factor(data[,1])
      levels(yfac) <- c("-1", "1")
      data[,1] <- as.numeric(as.character(yfac))
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
    if (sum(prior)!=1|any(prior<=0)|length(prior)!=2){
      stop("Invalid prior probabilities. Choose two values between 0 and 1 with sum 1 for parameter prior.")
    }
    if(!any(class == c("regfit", "lda"))){
      stop("Invalid classification method. Choose regfit or lda for parameter class.")
    }
    if (!yweights){
      probp4 <- FALSE
    }
    
    folds <- balancedfolds(data[,1], nfolds=min(min(table(data[,1])), nfold))    
    pcm <- matrix(nrow=n, ncol=length(as))
    wtest <- matrix(nrow=n, ncol=length(as))
    for (i in c(1:length(as))){
      for (f in 1:nfold){
        dtrain <- data[-folds[[f]],]
        dtest <- data[folds[[f]],]
        a <- as[i]
        trainmod <- prmda(formula,data=dtrain, a=a, fun, probp1,
                          hampelp2, hampelp3, probp4, yweights=yweights, class=class, prior=prior,center=center, scale=scale, numit, prec)
        pcm[folds[[f]],i] <- (((predict(trainmod, newdata=dtest)>0)-0.5)*2)
        wtest[folds[[f]],i] <- weig(trainmod, dtest, a, yweights=probp4)
      }
    }
    
    wpcm <- pcm*wtest
    
    if (sum(unique(data[,1])%in%c(-1,1))<2){
      yfac <- as.factor(data[,1])
      levels(yfac) <- c("-1", "1")
      data[,1] <- as.numeric(as.character(yfac))
    }
    ind1 <- which(data[,1]==1)
    ind2 <- which(data[,1]==-1)
    
    mwmcr <- apply(wpcm, 2, function(x) ((sum(abs(x[intersect(ind1, which(x<0))])))/sum(abs(x[ind1])) + sum(x[intersect(ind2, which(x>0))])/sum(abs(x[ind2])))/2)
    
    optind <- which.min(mwmcr)
    a <- as[optind[1]]
    
    prmspfit <- prmda(formula=formula, data=data, a=a, fun, probp1,
                      hampelp2, hampelp3, probp4, yweights=yweights, class=class, prior=prior, center, scale, numit, prec)
    
    if (plot==TRUE){
      plotdat <- data.frame(as=as, mwmcr=mwmcr)
      plotmwmcr <- ggplot(plotdat, aes(x=as, y=mwmcr)) + geom_line() + geom_point() + 
        geom_point(data=subset(plotdat,as==a), aes(x=as, y=mwmcr), color="red") +
        scale_x_continuous(breaks=as)
      print(plotmwmcr)
    }
    return(list(opt.mod=prmspfit, pcm=pcm))
  }
