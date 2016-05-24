GEE.var.pan <-
function(formula,id,family=gaussian,data,corstr="independence"){
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
  step<-matrix(0, nrow=cluster$n[1], ncol=cluster$n[1])
  for (i in 1:size){
  	y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    if (family=="gaussian"){
       resid<-(y-covariate%*%beta_est)%*%t(y-covariate%*%beta_est) 
       step<-step+resid
    }else if (family=="poisson"){
       resid<-(y-exp(covariate%*%beta_est))%*%t(y-exp(covariate%*%beta_est))
       B<-matrix(0,nrow=cluster$n[i],ncol=cluster$n[i])
       diag(B)<-1/sqrt(exp(covariate%*%beta_est))
       step<-step+B%*%resid%*%B	
    }else if (family=="binomial"){
       resid<-(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))%*%t(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))
       B<-matrix(0,nrow=cluster$n[i],ncol=cluster$n[i])
       diag(B)<-1/sqrt(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)
       step<-step+B%*%resid%*%B	
    }
  }
  unstr<-step/size
  #diag(unstr)<-rep(1, cluster$n[i])
  step11<-matrix(0,nrow=len,ncol=len)
  step12<-matrix(0,nrow=len,ncol=len)
  step13<-matrix(0,nrow=len_vec,ncol=1)
  step14<-matrix(0,nrow=len_vec,ncol=len_vec)
  p<-matrix(0,nrow=len_vec,ncol=size)
  for (i in 1:size){
  	y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    var_i=var[1:cluster$n[i],1:cluster$n[i]]
    if (family=="gaussian"){
       xy<-t(covariate)%*%solve(var_i)%*%unstr%*%solve(var)%*%covariate
       xx<-t(covariate)%*%solve(var_i)%*%covariate
       step11<-step11+xx
       step12<-step12+xy
       step13<-step13+vec(xy)
       p[,i]<-vec(xy) 
    }else if (family=="poisson"){
       A<-matrix(0,nrow=cluster$n[i],ncol=cluster$n[i])
       diag(A)<-exp(covariate%*%beta_est)
       D<-mat.prod(covariate, exp(covariate%*%beta_est))
       Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])
       xy<-t(D)%*%solve(Vi)%*%sqrt(A)%*%unstr%*%sqrt(A)%*%solve(Vi)%*%D
       xx<-t(D)%*%solve(Vi)%*%D
       step12<-step12+xy
       step11<-step11+xx
       step13<-step13+vec(xy)
       p[,i]<-vec(xy) 
    }else if (family=="binomial"){
       A<-matrix(0,nrow=cluster$n[i],ncol=cluster$n[i])
       diag(A)<-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2
       D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
       Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])
       xy<-t(D)%*%solve(Vi)%*%sqrt(A)%*%unstr%*%sqrt(A)%*%solve(Vi)%*%D
       xx<-t(D)%*%solve(Vi)%*%D
       step12<-step12+xy
       step11<-step11+xx
       step13<-step13+vec(xy)
       p[,i]<-vec(xy) 
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
