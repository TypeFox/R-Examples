get.svmstep<-function(lambda,lambda.old){
      lambda.range<-range(c(lambda,lambda.old))
      lambda.cuts<-c(lambda.range[2]+1,lambda.old,lambda.range[1]-1)
      as.numeric(cut(-lambda,-lambda.cuts))-1
    }
      
