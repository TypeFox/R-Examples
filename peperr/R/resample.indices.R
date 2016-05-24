`resample.indices` <-
function(n, sample.n=100, method=c("no", "cv", "boot", "sub632"))
{
   method <- match.arg(method)
   if (method=="no"){
      sample.index <- list(1:n)
      not.in.sample <- list(1:n)
   } 

   if (method=="boot"){
         indices <- matrix(sapply(1:sample.n,function(b){sample(1:n,replace=TRUE)}), byrow=TRUE, nrow=sample.n, ncol=n)
         sample.index <- list()
         not.in.sample <- list()
         for (i in 1:sample.n){
            not.in.sample[[i]] <- (1:n)[-unique(indices[i,])]
            sample.index[[i]] <- indices[i,]
         }
      }
   if (method=="sub632") {
         sample.index <- list()
         not.in.sample <- list()
         for (i in 1:sample.n) {
            sample.index[[i]] <- sample(n, round(n*.632), replace=FALSE)
            not.in.sample[[i]] <- (1:n)[-unique(sample.index[[i]])]
         }
      }

   if (method=="cv"){
      if (n<sample.n) stop("Number of observations 'n' has to be larger than sample size 'sample.n'")

      not.in.sample <- split(sample(1:n), rep(1:sample.n, length = n))
      sample.index <- list()
      for (i in 1:sample.n) {
         sample.index[[i]] <- (1:n)[-unique(not.in.sample[[i]])]
      }
   }
   list(sample.index=sample.index, not.in.sample=not.in.sample)
}

