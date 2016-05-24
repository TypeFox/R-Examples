ANOM <- function(mc, xlabel=NULL, ylabel=NULL, printn=T, printp=T,
                 stdep=NULL, stind=NULL, pst=NULL, pbin=NULL, bg="white"){
  
  if(!(class(mc)[1] %in% c("glht", "SimCi", "mctp", "binomRDci"))){
    stop("Please insert an object of class 'glht', 'SimCi', 'mctp', or 'binomRDci' for mc!")
  }
  
  if(class(mc)[1]=="glht"){
    
    if(mc$type!="GrandMean"){
      stop("For ANOM you need a 'GrandMean' contrast matrix!")
    }
    
    modclass <- attr(mc$model, "class")[1]
    
    if(modclass=="aov" | modclass=="lm" | modclass=="glm"){
      dd <- attr(attr(mc$model$terms, "factors"), "dimnames")[[1]][1]
      ii <- attr(attr(mc$model$terms, "factors"), "dimnames")[[1]][2]
      dep <- mc$model$model[dd][, 1]
      ind <- mc$model$model[ii][, 1]    
    }
    
    if(modclass=="lme"){
      dd <- attr(attr(mc$model$terms, "factors"), "dimnames")[[1]][1]
      ii <- attr(attr(mc$model$terms, "factors"), "dimnames")[[1]][2]
      dep <- mc$model$data[, dd]
      ind <- mc$model$data[, ii]
    }
    
    if(modclass=="lmerMod"){
      dd <- rownames(attr(attr(mc$model@frame, "terms"), "factors"))[1]
      ii <- rownames(attr(attr(mc$model@frame, "terms"), "factors"))[2]
      dep <- mc$model@frame[, dd]
      ind <- mc$model@frame[, ii]
    }
    
    ss <- as.vector(tapply(ind, ind, length))
    ci <- confint(mc)
    
    if(modclass!="glm"){
      means <- as.vector(tapply(dep, ind, mean))
      mea <- mean(dep)
    }
    
    if(modclass=="glm"){
      if(mc$model$family$family=="poisson"){
        cit <- apply(ci$confint, 2, poisson()$linkinv)
        m <- mean(poisson()$linkinv(predict(mc$model)))
      }
      if(mc$model$family$family=="binomial"){
        cit <- apply(ci$confint, 2, binomial()$linkinv)
        m <- mean(binomial()$linkinv(predict(mc$model)))
      }
    }
    
    if(is.null(xlabel)==T){
      xlabel <- ii
    }else{
      xlabel <- xlabel
    }
    
    if(is.null(ylabel)==T){
      ylabel <- dd
    }else{
      ylabel <- ylabel
    }
    
    if(printp==T){
      pvals <- summary(mc)$test$pvalues[1:length(ss)]
    }else{
      pvals <- NULL
    }
    
    bg <- match.arg(bg, choices=c("gray", "grey", "white"))
    
    if(modclass!="glm"){
      if(modclass=="lmerMod"){
        ANOMgen(mu=means, n=ss, lo=ci$confint[, "lwr"], up=ci$confint[, "upr"],
                names=levels(mc$model@frame[, ii]), alternative=mc$alternative,
                xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg)
      }else{
        ANOMgen(mu=means, n=ss, lo=ci$confint[, "lwr"], up=ci$confint[, "upr"],
                names=mc$model$xlevels[[1]], alternative=mc$alternative,
                xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg)
      }
    }
    
    if(modclass=="glm"){
      
      if(mc$model$family$family=="poisson"){
        if(mc$alternative=="two.sided"){
          ANOMintern(mu=m*cit[, 1], n=ss, gm=m, lo=m*cit[, 2]-m, up=m*cit[, 3]-m,
                     names=mc$model$xlevels[[1]], alternative=mc$alternative,
                     xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg, whichone="glm")
        }
        if(mc$alternative=="greater"){
          ANOMintern(mu=m*cit[, 1], n=ss, gm=m, lo=m*cit[, 2]-m, up=Inf,
                     names=mc$model$xlevels[[1]], alternative=mc$alternative,
                     xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg, whichone="glm")
        }
        if(mc$alternative=="less"){
          ANOMintern(mu=m*cit[, 1], n=ss, gm=m, lo=-Inf, up=m*cit[, 3]-m,
                     names=mc$model$xlevels[[1]], alternative=mc$alternative,
                     xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg, whichone="glm")
        }
      }
      if(mc$model$family$family=="binomial"){
        if(mc$alternative=="two.sided"){
          ANOMintern(mu=cit[, 1], n=ss, gm=m, lo=cit[ ,2]-m, up=cit[ ,3]-m,
                     names=mc$model$xlevels[[1]], alternative=mc$alternative,
                     xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg, whichone="glm")
        }
        if(mc$alternative=="less"){
          ANOMintern(mu=cit[, 1], n=ss, gm=m, lo=-Inf, up=cit[ ,3]-m,
                     names=mc$model$xlevels[[1]], alternative=mc$alternative,
                     xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg, whichone="glm")
        }
        if(mc$alternative=="greater"){
          ANOMintern(mu=cit[, 1], n=ss, gm=m, lo=cit[ ,2]-m, up=Inf,
                     names=mc$model$xlevels[[1]], alternative=mc$alternative,
                     xlabel=xlabel, ylabel=ylabel, printn=printn, p=pvals, bg=bg, whichone="glm")
        }
      }
      
    }
    
  }
  
  if(class(mc)[1]=="SimCi"){
    
    if(mc$type!="GrandMean"){
      stop("For ANOM you need a 'GrandMean' contrast matrix!")
    }
    
    if(is.null(stdep)==T){
      stop("Please insert a vector giving the values\n
           of the dependent variable for stdep!")
    }
    
    if(is.null(stind)==T){
      stop("Please insert a numeric vector giving the values\n
           of the independent variable for stind!")
    }
    
    if(is.numeric(stdep)==F){stop("The dependent variable must be numeric.")}
    
    if((length(stdep)==length(stind))==F){
      stop("Dependent and independent variable must be vectors of equal length.")
    }
    
    stind <- as.factor(stind)
    ss <- as.vector(tapply(stdep, stind, length))
    
    if(mc$test.class=="ratios"){
      grame <- 100
    }else{
      grame <- mean(stdep)
    }
    
    if(is.null(xlabel)==T){
      xlabel <- "group"
    }else{
      xlabel <- xlabel
    }
    
    if(is.null(ylabel)==T){
      ylabel <- mc$resp
    }else{
      ylabel <- ylabel
    }
    
    if(printp==T){
      
      if(class(pst)!="SimTest"){
        stop("Please insert an object of class 'SimTest' for pst,\n
             or set printp to 'FALSE'.")
      }
      
      ppp <- pst$p.val.adj[1:length(ss)]
      
      }else{
        
        ppp <- NULL
        
      }
    
    if(mc$test.class=="ratios"){
      ANOMintern(mu=100*as.vector(mc$estimate), n=ss, lo=100*as.vector(mc$lower),
                 up=100*as.vector(mc$upper), alternative=mc$alternative,
                 xlabel=xlabel, ylabel=ylabel, printn=printn, p=ppp, bg=bg, whichone="ratio")
    }else{
      ANOMgen(mu=as.vector(mc$estimate)+grame, n=ss, lo=as.vector(mc$lower),
              up=as.vector(mc$upper), alternative=mc$alternative,
              xlabel=xlabel, ylabel=ylabel, printn=printn, p=ppp, bg=bg)
    }
    
    }
  
  if(class(mc)[1]=="mctp"){
    
    if(mc$input$type!="UserDefined"){
      stop("For ANOM you need a 'UserDefined' grand-mean-type contrast matrix!")
    }
    
    if(is.null(mc$Correlation)==T){
      stop("Set the argument 'correlation' in function 'mctp()' to 'TRUE'.")
    }
    
    if(is.null(xlabel)==T){
      xlabel <- as.character(mc$input$formula[3])
    }else{
      xlabel <- xlabel
    }
    
    if(is.null(ylabel)==T){
      ylabel <- as.character(mc$input$formula[2])
    }else{
      ylabel <- ylabel
    }
    
    if(printp==T){
      pvals <- mc$Analysis$p.Value
    }
    
    if(mc$input$alternative=="two.sided"){
      ANOMgen(mu=mc$Data.Info$Effect, n=mc$Data.Info$Size,
              lo=mc$Analysis$Lower, #abs(mc$Analysis$Estimator - mc$Analysis$Lower),
              up=mc$Analysis$Upper, #abs(mc$Analysis$Estimator - mc$Analysis$Upper),
              names=colnames(mc$Contrast), alternative=mc$input$alternative,
              xlabel=xlabel, ylabel=ylabel,
              printn=printn, p=pvals, bg=bg)
    }
    
    if(mc$input$alternative=="greater"){
      ANOMgen(mu=mc$Data.Info$Effect, n=mc$Data.Info$Size,
              lo=mc$Analysis$Lower, #abs(mc$Analysis$Estimator - mc$Analysis$Lower),
              up=Inf,
              names=colnames(mc$Contrast), alternative=mc$input$alternative,
              xlabel=xlabel, ylabel=ylabel,
              printn=printn, p=pvals, bg=bg)
    }
    
    if(mc$input$alternative=="less"){
      ANOMgen(mu=mc$Data.Info$Effect, n=mc$Data.Info$Size,
              lo=Inf,
              up=mc$Analysis$Upper, #abs(mc$Analysis$Estimator - mc$Analysis$Upper),
              names=colnames(mc$Contrast), alternative=mc$input$alternative,
              xlabel=xlabel, ylabel=ylabel,
              printn=printn, p=pvals, bg=bg)
    }
    
  }
  
  if(class(mc)[1]=="binomRDci"){
    
    if(attr(mc$cmat, "type")!="GrandMean"){
      stop("For ANOM you need a 'GrandMean' contrast matrix!")
    }
    
    if(is.null(xlabel)==T){
      xlabel <- "group"
    }else{
      xlabel <- xlabel
    }
    
    if(is.null(ylabel)==T){
      ylabel <- "probability of success"
    }else{
      ylabel <- ylabel
    }
    
    if(printp==T){
      
      if(class(pbin)!="binomRDtest"){
        stop("Please insert an object of class 'binomRDtest' for pbin,\n
             or set printp to 'FALSE'.")
      }
      
      ppp <- pbin$p.val.adj
      
      }else{
        
        ppp <- NULL
        
      }
    
    ANOMgen(mu=mc$p, n=mc$n, lo=mc$conf.int[, "lower"],
            up=mc$conf.int[, "upper"], names=mc$names,
            alternative=mc$alternative, xlabel=xlabel, ylabel=ylabel,
            printn=printn, p=ppp, bg=bg)
    
  }
  
}