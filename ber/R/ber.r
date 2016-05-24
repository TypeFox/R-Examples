ber<-function(Y, b, covariates=NULL){

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

if(is.null(covariates)){
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m <- nlevels(b)
    X <- model.matrix(~b-1)
    Xinv <- ginv(X)
#    e_1 <- c(1, rep(0, m - 1))
#    e_1_matrix <- matrix(rep(e_1, n), ncol = m, byrow = TRUE)
    B <- Xinv %*% Y
    XB_hat<-X %*% B
    adjY <- Y - XB_hat
    res_squared <- adjY^2
    
    genes_mean <-  apply(XB_hat,2,mean)     
    genes_mean_matrix <- matrix(rep(genes_mean, n), byrow = TRUE,nrow = n)
    
    
    Bdouble <- Xinv %*% res_squared
    XD_hat<-X%*%Bdouble
    genes_var <-  apply(XD_hat,2,mean)     
    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE,nrow = n)
    ScaleFactors <- sqrt(genes_var_matrix/(XD_hat))
    adjYdouble <- ScaleFactors * adjY
#    adjYdouble <- adjYdouble +  e_1_matrix %*% B
    adjYdouble <- ScaleFactors * adjY +  genes_mean_matrix
    
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)      
    }else{
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m1 <- nlevels(b)
    X1 <- model.matrix(~b - 1)
    colnames(covariates) <- paste("col", 1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[,-1])
    m2 <- dim(X2)[2]
    X <- cbind(X1, X2)
    Xinv <- ginv(X)
    B <- Xinv %*% Y
    B1 <- B[1:m1,]
    B2 <- B[(m1 + 1):(m1 + m2),]
#    e_1 <- c(1, rep(0, m1 + m2 -1))
#    e_1_matrix <- matrix(rep(e_1, n), ncol = m1 + m2, byrow = TRUE)
    XB_hat<-X %*% B
    adjY <- Y - XB_hat
    res_squared <- adjY^2

    genes_mean <-  apply(X1 %*% B1,2,mean)     
    genes_mean_matrix <- matrix(rep(genes_mean, n), byrow = TRUE,nrow = n)

    
    X1inv <- ginv(X1)
    Bdouble1 <- X1inv %*% res_squared
    X1D1_hat<-X1%*%Bdouble1
    genes_var <- apply(X1D1_hat,2,mean)
    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE, nrow =n)

    
#    dd <- t(model.matrix(~factor(b) - 1))
#    matrix1 <- matrix(rep(1, n*g), ncol = n, byrow = TRUE)
#    nsamples <- dd %*% matrix1
#    Y1sum <- dd%*%res_squared
#    Bdouble <- Y1sum/nsamples
#    genes_var <- apply(Bdouble,2,mean)
#    genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE, nrow =n)
    

    ScaleFactors <- sqrt(genes_var_matrix/(X1D1_hat))
    adjYdouble <- ScaleFactors * adjY
#    adjYdouble <- adjYdouble + e_1_matrix %*% B  + X2%*%B2
    adjYdouble <- adjYdouble +  genes_mean_matrix  + X2%*%B2
    
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
    }
}
