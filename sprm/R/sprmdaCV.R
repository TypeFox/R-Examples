sprmdaCV <-
  function (formula, data, as, etas, nfold=10, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999, probp4=0.01, yweights=TRUE, class=c("regfit", "lda"), prior= c(0.5,0.5), center = "median", scale = "qn", print=FALSE, plot=TRUE, numit=100, prec=0.01)
    ## 14/10/09 IH
    ## Cross validation for partial robust M-regression for discriminanat analysis
    ##
    ## uses prmsDA, predict.prmda, weig
    ##
    ## formula .... an object of class formula.
    ## data ....... a data frame or list which contains the variables given in formula.
    ## as ......... a vector with positive integers, which are the number of PRMS components to be estimated in the models.
    ## etas ....... a vector of sparsity parameters; values between 0 and 1
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
    ## usesvd ..... logical, default is FALSE. If TRUE singular value decomposition is performed.
    ## plot ....... logical, default is TRUE. If TRUE a plot is generated with a measure of the prediction accuracy for each model (see Details).
    ## numit ...... the number of maximal iterations for the convergence of the coefficient estimates.
    ## prec ....... a value for the precision of estimation of the coefficients.
 
  {
    # library(pcaPP)
    # library(reshape2)
    #source("/home/hoffmann/Documents/Projects/Projekt PRM-DA/code/sprmsDA.R")
    #source("/home/hoffmann/Documents/Projects/Projekt PRM-DA/code/predict.sprmda.R")
    #source("/home/hoffmann/Documents/Projects/Projekt PRM-DA/code/weig.R")
    
    MCR <- number <- NULL
    
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
    pcm <- array(NA, dim = c(n, length(as), length(etas)))
    dimnames(pcm) <- list(1:n, 1:length(as), 1:length(etas))
    wtest <- array(NA, dim = c(n, length(as), length(etas)))
    dimnames(wtest) <- list(1:n, 1:length(as), 1:length(etas))
    nzcoef <- array(NA, dim = c(nfold, length(as), length(etas)))
    dimnames(nzcoef) <- list(1:nfold, 1:length(as), 1:length(etas))
    for (f in 1:nfold){
      dtrain <- data[-folds[[f]],]
      dtest <- data[folds[[f]],]
      for (i in 1:length(as)){
        a <- as[i]
        for (j in 1:length(etas)){
          e <- etas[j]
          res <- tryCatch({     
            trainmod <- sprmda(formula=formula, data=dtrain, a=a, eta=e, fun="Hampel", probp1 = probp1, hampelp2 = hampelp2, hampelp3 = hampelp3, 
								probp4=probp4, yweights=yweights, class=class, prior=prior, center = center, scale = scale, print=FALSE, numit=numit, prec=prec)
          }, error=function(err){
            error <- TRUE
            return(error)
          })          
          if (is.logical(res)){
            print(paste("sprmda broke off in CV loope for a=", a,"and eta=", e))
          } else{
            trainmod <- res
            pcm[folds[[f]],i,j] <- (((predict(trainmod, newdata=dtest)>0)-0.5)*2)
            wtest[folds[[f]],i,j] <- weig(x=trainmod, data=dtest, optcomp=a, yweights=probp4)
            nzcoef[f,i,j] <- sum(trainmod$coefficients!=0)
          }
        } 
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
    
    mwmcr <- apply(wpcm, c(2,3), function(x) ((-1*(sum(x[intersect(ind1, which(x<0))])))/sum(abs(x[ind1])) + sum(x[intersect(ind2, which(x>0))])/sum(abs(x[ind2])))/2)
    mnzcoef <- apply(nzcoef, c(2,3), mean)
    
    optind <- which(mwmcr==min(mwmcr, na.rm=TRUE), arr.ind=TRUE) 

    
    if (dim(optind)[1]>1){
    	optind <- matrix(optind[which(as[optind[,1]]==min(as[optind[,1]])),], byrow=FALSE, ncol=2)     # then lowest number of components
		if (dim(optind)[1]>1)
			optind <- matrix(optind[which(etas[optind[,2]]==max(etas[optind[,2]])),], byrow=TRUE, ncol=2) # first choose highest sparsity parameter
    }
    eta <- etas[optind[,2]]
    a <- as[optind[,1]]
    
    if (plot==TRUE){
      print(paste("optimal model: a =", a, "eta =", eta))
      mycol1 <- colorRampPalette(c("blue2","red", "orange", "yellow"))(64)
      mycol2 <- colorRampPalette(c("darkgreen", "yellow"))(64)
      ggmspe <- mwmcr
      rownames(ggmspe) <- as
      colnames(ggmspe) <- etas
      ggmspe <- melt(ggmspe)
      names(ggmspe) <- c("a", "eta", "MCR")
      mspeplot <- ggplot(ggmspe, aes(x=as.factor(a), y=as.factor(eta), fill=MCR)) +
        geom_tile() +
        scale_fill_gradientn(colours=mycol1) +
        ggtitle(paste0("Weighted Mean MCR (min at a=", a, ", eta=", eta, ")")) +
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
    
    sprmsDAFit <- sprmda(formula=formula, data=data, a=a, eta=eta, fun, probp1,
                      hampelp2, hampelp3,  probp4, yweights, class, prior, center, scale, print, numit, prec)
    
    return(list(opt.mod=sprmsDAFit, pcm=pcm, nzcoef=nzcoef))
    
  } 
    
