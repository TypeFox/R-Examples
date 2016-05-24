################################################################################
# a function that creates the response matrix and it's structure attribute     #
# for a multivariate glm with a binomial family                                #
# the returned matrix can be used in a formula for a manyglm                    #
################################################################################

# Y: a matrix with numeric values, the response matrix for a glm with a binomial
#    family
# ORDER: a numeric value, either
#
# "1" = y1 is a matrix of no. of successes, y2 the corresponding matrix of no. of failures
# "2" = y1 is a matrix of no. of successes, y2 the corresponding matrix of no. of trials
# or a a vector of no. of trials (the same for each row in y1)


binstruc <- function(y1, y2, order =1){

if(any(y1<0) | any(y2<0)) stop("binomial data must consist of integers only")

y1 <- as.matrix(y1)

if(nrow(y1)!=NROW(y2)) stop("'y1' must have the same number of rows as 'y2'")

if(order==1){

   nvars <- ncol(y1)
   y2 <- as.matrix(y2)
   if(is.null(colnames(y2))) { colnames(y2) <- colnames(y1)
   } else if(is.null(colnames(y1))) colnames(y1) <- colnames(y2)
   if(nvars!=ncol(y2)) stop("'y1' must have the same number of columns as 'y2'")
   y <- cbind(y1,y2)
   colnam <- colnames(y)

   if(is.null(colnam)) {
     colnames(y) <- c(paste("succ",1:nvars, sep=""), paste("fail",1:nvars, sep=""))
     struc <-  c( rep("succ", nvars), rep("fail", nvars) )
   } else {
       struc <- c( rep("succ",times = nvars), rep("fail", times=nvars) )
       colnames(y) <- paste(struc , ".", colnam, sep = "")
   }

   n <- y[ ,  struc=="succ", drop=FALSE ] + y[, struc=="fail", drop=FALSE]
   if(any(n - n[,1] != 0))
        stop("the number of trials must be the same for each variable")
   
} else  if(order==2){

   nvars <- ncol(y1)
   if(NCOL(y2)==1){
     y2 <- c(y2)
   } else if(nvars!=ncol(y2)) {
      stop("'y2' must either have the same number of columns as 'y1' or be a vector")
   } else if(any(y2- y2[,1]!=0))
      stop("the number of trials must be the same for each variable")
   
   if(any(y2-y1<0)) stop("the number of trials must be larger than the number of successses")
   y <- cbind(y1,y2-y1)
   colnam <- colnames(y)

   if(is.null(colnam)) {
      colnames(y) <- c(paste("succ",1:nvars, sep=""), paste("fail",1:nvars, sep=""))
   } else {
       struc <- c( rep("succ",times = nvars), rep("fail", times=nvars) )
       colnames(y) <- paste(struc , ".", colnam, sep = "")
   }

   
} else  stop("'order' must be in either '1' or '2'")

return(y)

}


