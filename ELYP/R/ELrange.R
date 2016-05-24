ELrange <- function(mle, loglik, step = rep(1,length(mle)), DataMat){
################################################
#######  the DataMat is for loglik( ) use.
#######  the loglik fun has two input, mle and DataMat
####################################################
MaxValue <- loglik(mle, DataMat)
BorderV <- MaxValue - 4
k <- length(mle)
step <- rep(NA, k)
temp <- matrix(NA, nrow=k, ncol=2)

for(J in 1:k)  { 
    step[J] <- 1      # any better initial value?
    para <- mle
    for( i in 1:10 ) {
    para[J] <- mle[J]+step[J]
    temp[J,1] <- loglik(para, DataMat)
    para[J] <- mle[J]-step[J]
    temp[J,2] <- loglik(para, DataMat)
    if ( (temp[J,1]+temp[J,2]) < BorderV )  {step[J] <- step[J]/2} ##??
        else  {step[J] <- 2*step[J]}
    }
}

#step2 <- 1
#for( i in 1:10 ) {
#temp21 <- loglik(c(mle[1], mle[2]+step2), DataMat)
#temp22 <- loglik(c(mle[1], mle[2]-step2), DataMat)
#if( (temp21+temp22) < BorderV ) {step2 <- step2/2}
#  else {step2 <- 2*step2}
#}

list( Steps = step, TempV = temp )
}
