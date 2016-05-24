#' Bootstrap ranking LASSO model.
#'
#' This function performs a LASSO logistic regression model using a bootstrap ranking procedure.
#'
#' @param x the predictor matrix
#' @param y the response variable, a factor object with values of 0 and 1 
#' @param B the external loop for intersection operation, with the default value 5
#' @param Boots the internal loop for bootstrap sampling, with the default value 100
#' @param kfold the K-fold cross validation, with the default value 10
#' @export
#' @import glmnet
#' @import SiZer
#' @references
#' Guo, P., Zeng, F., Hu, X., Zhang, D., Zhu, S., Deng, Y., & Hao, Y. (2015). Improved Variable 
#' Selection Algorithm Using a LASSO-Type Penalty, with an Application to Assessing Hepatitis B 
#' Infection Relevant Factors in Community Residents. PLoS One, 27;10(7):e0134151.
#' @examples
#' library(datasets)
#' head(iris)
#' X <- as.matrix(subset(iris,iris$Species!="setosa")[,-5])
#' Y <- as.factor(ifelse(subset(iris,iris$Species!="setosa")[,5]=='versicolor',0,1))
#' # Fitting a bootstrap ranking LASSO (BRLASSO) logistic regression model
#' BRLasso.fit <- BRLasso(x=X, y=Y, B=2, Boots=10, kfold=10)
#' # Variables selected by the BRLASSO model
#' BRLasso.fit$var.selected
#' # Coefficients of the selected variables
#' BRLasso.fit$var.coef  
BRLasso=function(x, y, B=5, Boots=100, kfold=10){
    varx <- colnames(x)
    varx0 <- varx
    nvar <- length(varx)
    rowx <- nrow(x)
    n <- length(y)
    var_list_all <- vector("list",B)
    if (rowx!=n){
      stop("The number of rows in x is not equal to the length of y!")
        }
    if (nvar==1){
      stop("The x matrix has only one variable!")
        }
    RecordM2=matrix(0, nrow=B, ncol=nvar)
    rownames(RecordM2)=paste("B", 1:B, sep="")
    colnames(RecordM2)=varx
for (ij in 1:B){
    Nboot <- n
    RecordM=matrix(0, nrow=Boots, ncol=nvar)
    rownames(RecordM)=paste("BSample", 1:Boots, sep="")
    colnames(RecordM)=varx 
    for(i in 1:Boots){
      repeat{ 
        s <- sample(Nboot, replace=TRUE)  
        if(length(table(y[s])) >= 2 & length(table(y[-s])) >= 2)
        break
             }
      BoostrapX <- x[s, ]
      BoostrapY <- y[s]
      cvfit <- cv.glmnet(x=BoostrapX, y=BoostrapY, type.measure="deviance", 
                         nfolds=kfold, family="binomial")
      model.final <- cvfit$glmnet.fit
      nzero <- as.matrix(coef(model.final, s=cvfit$lambda.min))
      var_nz <- names(nzero[nzero[,1]!=0,])
      for(i0 in 1:dim(nzero)[1]){
        for(j0 in 1:nvar){
           if (dimnames(nzero)[[1]][i0]==colnames(RecordM)[j0]){
               RecordM[i, j0]=nzero[i0, ]}
                          }
                                 }
      cat("Step 1-Current iteration: ", i, "\n")
                      }
    Impor <- abs(apply(RecordM, 2, mean))
    score <- sort(Impor)
    index <- seq(1:length(score))
    model.pwl <- piecewise.linear(x=index, y=score, middle=1, CI=FALSE, 
                                  bootstrap.samples=1000, sig.level=0.05)
    break_p=round(model.pwl$change.point)
    var_list_all[[ij]]=names(score[break_p:length(index)])
    if (length(var_list_all[[ij]])!=0){
      varx0 <- intersect(var_list_all[[ij]], varx0)}
    Impor0 <- apply(RecordM, 2, mean)[var_list_all[[ij]]]
    for(gp in 1:length(Impor0)){
      RecordM2[ij, colnames(RecordM2)==names(Impor0)[gp]]=Impor0[gp]
                                }
    cat("Step 2-Current iteration: ", ij, "\n")
} 
    if(length(varx0)==1){
      RecordM2 <- as.matrix(RecordM2[, varx0],ncol=1)
      colnames(RecordM2) <- varx0
                         } else {
      RecordM2 <- RecordM2[, varx0]
                                }
    Myresult <- list(var.selected=colnames(RecordM2), var.coef=apply(RecordM2, 2, mean))
    return(Myresult)
}