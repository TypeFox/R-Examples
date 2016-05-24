GEE.var.fg <-
function(formula,id,family=gaussian,data,corstr="independence",b=0.75){
  #########################################################################
  # Arguments:
  # formula  specify the model of interest
  # family   "gaussian", "binomial" or "poisson"
  # data     data frame
  # corstr   Working correlation structure: "independence", "AR-M", "exchangeable", "unstructured".
  # value:   GEE returns a list containing the following elements
  #          cov.beta     estimate of robust variance for \hat{\beta}
  #          cov.var      estimate of the variance-covariance matrix for robust variance.
  #########################################################################
  # Delete the records with missing data in predictors or outcomes;
  if (is.null(data$id)){
        index <- which(names(data)==id)
        data$id <- data[,index]}
        
  ### na.action: only na.omit is used for gee;
  init <- model.frame(formula, data)
  init$num <- 1:length(init[,1])
  if(any(is.na(init))){
    index <- na.omit(init)$num
    data <- data[index,]
    ### Get the design matrix;
    m <- model.frame(formula, data)
    mt <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }else{
  	### Get the design matrix;
    m <- model.frame(formula, data)
    mt <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  
  ### Fit the GEE model to get the estimate of parameters \hat{\beta};
  gee.fit <- gee(formula,data=data,id=id,family=family,corstr=corstr)
  beta_est <- gee.fit$coefficient
  alpha <- gee.fit$working.correlation[1,2]
  len <- length(beta_est)
  len_vec <- len^2
  
  ### Estimate the robust variance for \hat{\beta}
  data$id <- gee.fit$id
  cluster<-cluster.size(data$id)
  ncluster<-max(cluster$n)
  size<-cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if(is.character(corstr)){
    var <- switch(corstr,
               "independence"=cormax.ind(ncluster),
               "exchangeable"=cormax.exch(ncluster, alpha),
               "AR-M"=cormax.ar1(ncluster, alpha),
               "unstructured"=summary(gee.fit)$working.correlation,)
  }else{
       print(corstr)
       stop("'working correlation structure' not recognized")
  }   
  if(is.character(family)){
    family <- switch(family,
                  "gaussian"="gaussian",
                  "binomial"="binomial",
                  "poisson"="poisson")
  }else{ 
    if(is.function(family)){
      family <- family()[[1]]
    }else{
       print(family)
       stop("'family' not recognized")
    }    
  }

  cov.beta<-unstr<-matrix(0,nrow=len,ncol=len)
  step11<-matrix(0, nrow=len, ncol=len)
  for (i in 1:size){
  	y<-as.matrix(data$response[data$id==unique(data$id)[i]])
  	covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      ncluster=cluster$n[i]
      var1=var[1:ncluster,1:ncluster] 
    if (family=="gaussian"){ 
       Vi=gee.fit$scale*var
       xx<-t(covariate)%*%solve(Vi)%*%covariate
       step11<-step11+xx  
    }else if (family=="poisson") {
       D<-mat.prod(covariate, exp(covariate%*%beta_est))
       Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)
       xx<-t(D)%*%solve(Vi)%*%D
       step11<-step11+xx
    }else if (family=="binomial"){
    	D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)
        xx<-t(D)%*%solve(Vi)%*%D
        step11<-step11+xx 
    }
  }
  step12<-matrix(0,nrow=len,ncol=len)
  step13<-matrix(0,nrow=len_vec,ncol=1)
  step14<-matrix(0,nrow=len_vec,ncol=len_vec)
  p<-matrix(0,nrow=len_vec,ncol=size)
  for (i in 1:size){
     y<-as.matrix(data$response[data$id==unique(data$id)[i]])
     covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
     ncluster=cluster$n[i]
      var1=var[1:ncluster,1:ncluster]
    if (family=="gaussian"){ 
    	## set up the scale parameter;
      Vi=gee.fit$scale*var
    	xx<-t(covariate)%*%solve(Vi)%*%covariate
    	Qi <- xx%*%solve(step11)
    	Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
       xy<-Ai%*%t(covariate)%*%solve(Vi)%*%(y-covariate%*%beta_est)
        step12<-step12+xy%*%t(xy)
        step13<-step13+vec(xy%*%t(xy))
        p[,i]<-vec(xy%*%t(xy))
     }else if (family=="poisson") {
     	## set up the scale parameter;
    	D<-mat.prod(covariate, exp(covariate%*%beta_est))
    	Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)
        xx<-t(D)%*%solve(Vi)%*%D
        Qi <- xx%*%solve(step11)
    	  Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
        xy<-Ai%*%t(D)%*%solve(Vi)%*%(y-exp(covariate%*%beta_est))
        step12<-step12+xy%*%t(xy)
        step13<-step13+vec(xy%*%t(xy))
        p[,i]<-vec(xy%*%t(xy))
     }else if (family=="binomial"){
       ## set up the scale parameter;
       D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
       Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)
       xx<-t(D)%*%solve(Vi)%*%D
       Qi <- xx%*%solve(step11)
    	 Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
       xy<-Ai%*%t(D)%*%solve(Vi)%*%(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))
       step12<-step12+xy%*%t(xy)
       step13<-step13+vec(xy%*%t(xy))
       p[,i]<-vec(xy%*%t(xy)) 
    }    
  }
  for (i in 1:size){
    dif<-(p[,i]-step13/size)%*%t(p[,i]-step13/size)
    step14<-step14+dif
  }
  cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
  cov.var<-size/(size-1)*kronecker(solve(step11), solve(step11))%*%step14%*%kronecker(solve(step11), solve(step11))
  return(list(cov.beta=diag(cov.beta), cov.var=cov.var))
}
