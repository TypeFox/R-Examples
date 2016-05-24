ber_bg<-function(Y, b, covariates=NULL,partial=TRUE,nSim=150){

################################################################################

if(missing(Y)){stop("Argument 'Y' missing, with no default\n")}

if(missing(b)){stop("Argument 'b' missing, with no default\n")}

if(class(Y)!='matrix'){stop("'Y' must be of class 'matrix'\n")}

if(class(b)!='factor'){stop("'b' must be of class 'factor'\n")}

if(any(is.na(Y))){stop("NA values are not allowed in 'Y'\n")}

if(any(is.na(b))){stop("NA values are not allowed in 'b'\n")}

if(length(b)!=nrow(Y)){stop("length(b) is different from nrow(Y)\n")}

if(any(apply(Y,2,mode)!='numeric')){stop('Array expression columns contain non-numeric values!\n')}

################################################################################

if(!is.null(covariates)){
if(class(covariates)!="data.frame"){stop("'covariates' must be of class 'data.frame'\n")}

col.cov<-ncol(covariates)
for(i in 1:col.cov){if(class(covariates[,i])!="numeric" & class(covariates[,i])!="factor"){
stop("column ", i, " of 'covariates' must be of class 'factor' or 'numeric'\n")}}

if(any(is.na(covariates))){stop("NA values are not allowed in 'covariates'\n")}
}

################################################################################

cnames<-colnames(Y)
rnames<-rownames(Y)

################################################################################
################################################################################
if(is.null(covariates) & partial==TRUE)
{
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m <- nlevels(b)
       
#######################################    
# B bagging estimator     

B <- matrix(rep(0,m*g), ncol = g, byrow = TRUE)

    for(i in 1:nSim)
    {
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
             
    X_bag <- model.matrix(~b_bag-1)
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
     
    B<-B+B_bag   
    }
    B<-B/nSim
    
#####################################
# data adjustment

    X <- model.matrix(~b-1)
    Xinv <- ginv(X)
    
    XB_hat<-X %*% B
    genes_mean <-  apply(XB_hat,2,mean)     
    genes_mean_matrix <- matrix(rep(genes_mean, n), byrow = TRUE,nrow = n)

    adjY <- Y - XB_hat
    res_squared <- adjY^2
    
    Bdouble <- Xinv %*% res_squared                                             
    XD_hat<-X%*%Bdouble
    genes_var <-  apply(XD_hat,2,mean)     
    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE,nrow = n)
    ScaleFactors <- sqrt(genes_var_matrix/(XD_hat))
    adjYdouble <- ScaleFactors * adjY
    adjYdouble <- adjYdouble +  genes_mean_matrix
    
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
}

################################################################################
################################################################################
if(is.null(covariates) & partial==FALSE)
{
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m <- nlevels(b)
       
#######################################    
# bagging estimators     

B <- matrix(rep(0,m*g), ncol = g, byrow = TRUE)
Bdouble <- matrix(rep(0,m*g), ncol = g, byrow = TRUE)

    for(i in 1:nSim)
    {
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
             
    X_bag <- model.matrix(~b_bag-1)
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
    
    adjY_bag <- Y_bag - X_bag %*% B_bag
    res_bag <- Y_bag - X_bag %*% B_bag
    res_squared_bag <- res_bag^2
    Bdouble_bag <- Xinv_bag %*% res_squared_bag
     
    B<-B+B_bag
    Bdouble<-Bdouble+Bdouble_bag   
    }
    B<-B/nSim
    Bdouble<-Bdouble/nSim
    
#####################################
# data adjustment

    X <- model.matrix(~b-1)
    Xinv <- ginv(X)
        
    XB_hat<-X %*% B
    genes_mean <-  apply(XB_hat,2,mean)     
    genes_mean_matrix <- matrix(rep(genes_mean, n), byrow = TRUE,nrow = n) 
    
    adjY <- Y - XB_hat             

    XD_hat<-X%*%Bdouble
    genes_var <-  apply(XD_hat,2,mean)     
    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE,nrow = n)
    ScaleFactors <- sqrt(genes_var_matrix/(XD_hat))
    adjYdouble <- ScaleFactors * adjY
    adjYdouble <- adjYdouble +  genes_mean_matrix
    
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
}


################################################################################
################################################################################
if(!is.null(covariates) & partial==TRUE)
{
    library(MASS)
   
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m1 <- nlevels(b)
    X1 <- model.matrix(~b-1)
    colnames(covariates) <- paste("col",1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[,-1])
    X <- cbind(X1, X2)
    m2 <- dim(X2)[2]
    
#####################################
# bagging estimators

    B_Aggregating<-matrix(rep(0,g*(m1+m2)),ncol=g,nrow=m1+m2)
   
for(i in 1:nSim){
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
    covariates_bag <- covariates[sample_labels,]
    covariates_bag <- as.data.frame(covariates_bag)
           
    X1_bag <- model.matrix(~b_bag-1)
    colnames(covariates_bag) <- paste("col",1:ncol(covariates_bag),sep = "")
    fmla <- as.formula(paste("~",paste(colnames(covariates_bag),collapse = "+")))
    X2_bag <- model.matrix(fmla,covariates_bag)
    X2_bag <- as.matrix(X2_bag[,-1])
    X_bag <- cbind(X1_bag, X2_bag)
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag

    
   B_Aggregating <- B_Aggregating+B_bag   
   }

   B_Aggregating <- B_Aggregating/nSim 
   res_squared <- (Y - X %*% B_Aggregating)^2

#####################################
# data adjustment

    X1inv <-  ginv(X1)
    Bdouble_Aggregating <- X1inv %*% res_squared

    B1 <- B_Aggregating[1:m1,]
    B2 <- B_Aggregating[(m1 + 1):(m1 + m2),]
    
    genes_mean <-  apply(X1 %*% B1,2,mean)     
    genes_mean_matrix <- matrix(rep(genes_mean, n), byrow = TRUE,nrow = n) 

    adjY <- Y - X %*% B_Aggregating
    res_squared <- adjY^2

    X1D1_hat<-X1%*%Bdouble_Aggregating
    genes_var <- apply(X1D1_hat,2,mean)
    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE, nrow =n)
    ScaleFactors <- sqrt(genes_var_matrix/(X1D1_hat))
    adjYdouble <- ScaleFactors * adjY
    adjYdouble <- adjYdouble + genes_mean_matrix  + X2%*%B2

    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)

}

################################################################################
################################################################################
if(!is.null(covariates) & partial==FALSE)
{
    library(MASS)
   
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m1 <- nlevels(b)
    X1 <- model.matrix(~b-1)
    colnames(covariates) <- paste("col", 1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[, -1])
    m2 <- dim(X2)[2]

    B_Aggregating_bag<-matrix(rep(0,g*(m1+m2)),ncol=g,nrow=m1+m2)
    Bdouble_Aggregating_bag<-matrix(rep(0,g*m1),ncol=g,nrow=m1)

#####################################
# bagging estimators
   
for(i in 1:nSim){
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
    covariates_bag <- covariates[sample_labels,]
    covariates_bag <- as.data.frame(covariates_bag)
           
    X1_bag <- model.matrix(~b_bag-1)
    colnames(covariates_bag) <- paste("col",1:ncol(covariates_bag),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates_bag),collapse = "+")))
    X2_bag <- model.matrix(fmla, covariates_bag)
    X2_bag <- as.matrix(X2_bag[,-1])
    X_bag <- cbind(X1_bag,X2_bag)
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
    adjY_bag <- Y_bag - X_bag %*% B_bag
    res_squared_bag <- adjY_bag^2


    X1inv_bag <- ginv(X1_bag)
    Bdouble_bag <- X1inv_bag %*% res_squared_bag


    B_Aggregating_bag <- B_Aggregating_bag+B_bag 
    Bdouble_Aggregating_bag <- Bdouble_Aggregating_bag+Bdouble_bag  
}

   B_Aggregating <- B_Aggregating_bag/nSim
   Bdouble_Aggregating <- Bdouble_Aggregating_bag/nSim

#####################################
# data adjustment

    X <- cbind(X1, X2)
    Xinv <- ginv(X)
    B1 <- B_Aggregating[1:m1,]
    B2 <- B_Aggregating[(m1 + 1):(m1 + m2),]
    adjY <- Y - X %*% B_Aggregating

    genes_mean <-  apply(X1 %*% B1,2,mean)     
    genes_mean_matrix <- matrix(rep(genes_mean, n), byrow = TRUE,nrow = n) 

    X1D1_hat<-X1%*%Bdouble_Aggregating
    genes_var <- apply(X1D1_hat,2,mean)
    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE, nrow =n)
    ScaleFactors <- sqrt(genes_var_matrix/(X1D1_hat))
    adjYdouble <- ScaleFactors * adjY
    adjYdouble <- adjYdouble + genes_mean_matrix  + X2%*%B2

    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)    
}

################################################################################
################################################################################

}
