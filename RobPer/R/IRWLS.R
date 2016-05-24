IRWLS <- function(yy, matrix_, W, residuals_, scale_, tol){
    # Performs IRWLS
    # Used by singlePerM
    # Arguments not checked as it is assumed that the function singlePer gives the right arguments
    res <- residuals_
    repeat {
        res_old<- res
        yy_hat <- yy * sqrt(W(res_old/scale_))
        mat_hat <-  matrix_ * sqrt(W(res_old/scale_))
        if(any(apply(mat_hat,2, sum)==0)) {
            bad_col<- which(apply(mat_hat,2, sum)==0)
            for(bc in bad_col){
                mat_hat[which(matrix_[, bc]>0), bc] <- 1
            }
        }
        tempIRWLS  <- lm(yy_hat~0+ mat_hat)
        columns<- 1:dim(matrix_)[2]
        if( any(is.na(tempIRWLS$coeff))) {
            nas<-which(is.na(tempIRWLS$coeff))
            columns<- columns[-nas]
        }
        res <- as.vector(yy-matrix_[,columns]%*%as.matrix(tempIRWLS$coeff[columns]))
        if(max(abs((res_old -res)/scale_)) < tol) break
    }
    return(tempIRWLS$coeff)
}
