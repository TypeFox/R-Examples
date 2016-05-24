DLR <- function(basemodel, augmentedmodel, diseasestatus,
                dataset, clustervar=NULL, alpha=0.05) {
  
  #require(gee)
  
  ### modelSettings
  d       <- diseasestatus   
  x       <- formula(basemodel)
  v       <- formula(augmentedmodel)
  dat     <- dataset
  dat[,d] <- as.numeric(as.character(dat[,d])) # needs to be numeric

  
  # make sure formula object v is ordered properly in order to correctly calculate difference in coefficients later on
  xterms <- paste(attr(terms(x),"term.labels"),collapse="+")
  yterms <- paste(setdiff(attr(terms(v),"term.labels"),attr(terms(x),"term.labels")),collapse="+")
  v      <- formula(paste("~",paste(xterms,yterms,sep="+")))

  # subject id
  # in case optional cluster variable is specified
  if (!is.null(clustervar)) dat$subjid <- dat[,clustervar] else dat$subjid <- 1:nrow(dat)  
  
  ### Gu and Pepe, Biostatistics, 2009, appendix
  # step 1: duplicate records 
  datdup <- rbind(dat,dat)
  
  # step 2: define indicator variable
  datdup$indicator <- rep(1:0,each=nrow(dat))   
  
  # step 3: define Z
  # first define model matrix to handle factors/interactions properly
  x.modelmatrix           <- model.matrix(x, datdup)
  colnames(x.modelmatrix) <- make.names(gsub(".Intercept.","intercept",paste("x",colnames(x.modelmatrix),sep=".")))
  dim.x                   <- ncol(x.modelmatrix) 
  v.modelmatrix           <- model.matrix(v, datdup)
  colnames(v.modelmatrix) <- make.names(gsub(".Intercept.","intercept",paste("v",colnames(v.modelmatrix),sep=".")))
  dim.v                   <- ncol(v.modelmatrix)
  # define and append Z
  z                       <- data.frame((1 - datdup$indicator) * x.modelmatrix, datdup$indicator * v.modelmatrix)
  datdup                  <- cbind(datdup,z)
  
  # step 4/5: fit simultaneously two logistic regression models in a GEE model setup for disease status d and predictors Z
  # prior to fit, data has to be sorted by cluster/subject
  datdup   <- datdup[order(datdup$subjid),]
  gee_xv   <- gee(as.formula(paste("d ~ ", paste(colnames(z),collapse="+"),"-1")), # no common intercept due to x an v specific ones
                  data=datdup, corstr="independence",family=binomial, id=datdup$subjid, maxiter=100, tol=1e-5)
  # model estimates
  coefs    <- coef(gee_xv)
  robvcov  <- gee_xv$robust.variance # sandwich estimators for all x and v coefficients, this is actually robust.SE^2
  
  ### inference
  # get coefficients in the logDLR model, difference between augmented and base model
  coef_logDLR          <- coefs[colnames(v.modelmatrix)]
  coef_logDLR[1:dim.x] <- coef_logDLR[1:dim.x] - coefs[colnames(x.modelmatrix)]
  names(coef_logDLR)   <- substr(names(coef_logDLR),3,100)
  
  # get covariance matrix for coefficients in the logDLR model
  # see chapter 3 in Gu and Pepe, 2009
  # the definition of I and Zero matrices is slightly more flexible/complex than in chapter 3 in order to allow for models w/o interaction 
  I.x.matrix <- diag(1,dim.x)                            # I matrix for x covariates
  Z.x.matrix <- diag(0,nrow=dim.x, ncol=dim.v - dim.x)   # Zero matrix for x covariates
  I.v.matrix <- diag(1,dim.v - dim.x)                    # I matrix for v covariates
  Z.v.matrix <- diag(0,ncol=dim.x, nrow=dim.v - dim.x)   # Zero matrix for v covariates
  A.matrix   <- rbind(cbind(-I.x.matrix,I.x.matrix,Z.x.matrix),cbind(Z.v.matrix,Z.v.matrix,I.v.matrix)) # A matrix

  # robust SE of logDLR coefficients
  robSE_logDLR    <- sqrt(diag(A.matrix %*% robvcov %*% t(A.matrix)))
  
  # corresponding z-statistic
  robZstat_logDLR <- coef_logDLR/robSE_logDLR
  
  # and two-sided p-values
  pvalues         <- 2 * pnorm(abs(robZstat_logDLR), lower.tail = FALSE)   
  
  # confidence intervals using normal approximation
  lowerci_logDLR  <- coef_logDLR - qnorm(1-alpha/2) * robSE_logDLR
  upperci_logDLR  <- coef_logDLR + qnorm(1-alpha/2) * robSE_logDLR
  
  ### results for logDLR model
  logPreTestResults  <- summary(gee_xv)$coefficients[colnames(x.modelmatrix),,drop=F]
  logPostTestResults <- summary(gee_xv)$coefficients[colnames(v.modelmatrix),,drop=F]
  logDLRResults      <- data.frame(coefs=coef_logDLR,robustSE=robSE_logDLR,Zstat=robZstat_logDLR,pvalue=pvalues,lowerCI=lowerci_logDLR,upperCI=upperci_logDLR)
  
  ###  DLR(Y|X) 
  #    DLR+, ie. DLR(Y=1|X=1), and DLR-, ie DLR(Y=0,X=1)
  xoneyone   <- rep(1,ncol(v.modelmatrix)) # X=1, Y=1
  xoneyzero  <- c(rep(1,ncol(x.modelmatrix)),rep(0,ncol(v.modelmatrix)-ncol(x.modelmatrix))) # X=1, Y=0
  logDLR_pos <- coef_logDLR  %*% xoneyone
  logDLR_neg <- coef_logDLR  %*% xoneyzero
  DLRy_pos   <- exp(logDLR_pos)
  DLRy_neg   <- exp(logDLR_neg)
  
  # corresponding robust SE (chapter 3, Gu and Pepe, 2009) at X=1, and confidence intervals
  SE_logDLR_pos <- sqrt(xoneyone  %*% A.matrix %*% robvcov %*% t(A.matrix) %*% xoneyone)
  SE_logDLR_neg <- sqrt(xoneyzero %*% A.matrix %*% robvcov %*% t(A.matrix) %*% xoneyzero)

  lowerCI_DLRy_pos <- exp(logDLR_pos - qnorm(1-alpha/2) * SE_logDLR_pos)
  upperCI_DLRy_pos <- exp(logDLR_pos + qnorm(1-alpha/2) * SE_logDLR_pos)
  lowerCI_DLRy_neg <- exp(logDLR_neg - qnorm(1-alpha/2) * SE_logDLR_neg)
  upperCI_DLRy_neg <- exp(logDLR_neg + qnorm(1-alpha/2) * SE_logDLR_neg)
  
  DLR_results <- data.frame(rbind(c(DLRy_pos,lowerCI_DLRy_pos,upperCI_DLRy_pos),
                                  c(DLRy_neg,lowerCI_DLRy_neg,upperCI_DLRy_neg)))
  colnames(DLR_results) <- c("DLR","lowerCI","upperCI")
  rownames(DLR_results) <- c("DLR(Y=1|X=1)","DLR(Y=0|X=1)")
  
  ### output
  return(list("logPreTestModel"  = logPreTestResults,
              "logPostTestModel" = logPostTestResults,
              "logDLRModel"      = logDLRResults,
              "DLR"              = DLR_results))
}
