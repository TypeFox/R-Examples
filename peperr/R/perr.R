`perr` <-
function(peperrobject, type= c("632p", "632", "apparent", "NoInf", "resample", "nullmodel")){
   if (is.list(peperrobject$selected.complexity)){
   length.cplx <- length(peperrobject$sample.complexity)/length(peperrobject$selected.complexity)/length(peperrobject$indices$sample.index)
   } else {
   length.cplx <- length(peperrobject$selected.complexity)
   }
   sample.n <- length(peperrobject$indices$sample.index)
   full.apparent <- peperrobject$full.apparent
   noinf.error <- peperrobject$noinf.error
   sample.error.l <- peperrobject$sample.error
   type <- match.arg(type)

   if (!is.null(peperrobject$noinf.error)){
      eval.times <- peperrobject$attribute
      if (is.null(eval.times)) eval.times <- 1
      prediction.error <- matrix(ncol=length(eval.times), nrow=length.cplx)
      mean.sample.error <- matrix(ncol=length(eval.times), nrow=length.cplx)

      for (i in 1:length.cplx){
         sample.error <- matrix(ncol=length(eval.times), nrow=sample.n)
         for (j in 1:sample.n){
	   sample.error[j,1:length(sample.error.l[[j]][i,])] <- sample.error.l[[j]][i,]
         }
         mean.sample.error[i,] <- apply(sample.error, 2, mean, na.rm=TRUE)
         if (type == "632p"){
            relative.overfit <- ifelse(noinf.error[i,] > full.apparent[i,],
                                        (ifelse(mean.sample.error[i,] < noinf.error[i,],
                                         mean.sample.error[i,],noinf.error[i,]) - full.apparent[i,])/
                                        (noinf.error[i,] - full.apparent[i,]), 0)
            weights <- .632/(1-.368*relative.overfit)

            prediction.error[i,] <- (1-weights)*full.apparent[i,] + 
                                       weights*ifelse(mean.sample.error[i,] < noinf.error[i,], 
                                       mean.sample.error[i,], noinf.error[i,])
         } 
         if (type == "resample"){
            prediction.error <- mean.sample.error
         }
         if (type == "632"){
            prediction.error[i,] <- (1 - 0.632)*full.apparent[i,] + 0.632*mean.sample.error[i,]
         }
      }
#       if(nrow(prediction.error)==1){
#         prediction.error <- as.vector(prediction.error)
#       }
       if (type == "apparent"){
         prediction.error <- full.apparent
      }
      if (type == "NoInf"){
         prediction.error <- noinf.error
      }
      if (type == "nullmodel"){
         prediction.error <- peperrobject$null.model
         attr(prediction.error, "addattr") <- NULL
      }

   } else {
 
      mean.sample.error <- matrix(ncol=1, nrow=length.cplx)

      for (i in 1:length.cplx){
         sample.error <- matrix(ncol=1, nrow=sample.n)
         for (j in 1:sample.n){
	   sample.error[j,1:length(sample.error.l[[j]][i,])] <- sample.error.l[[j]][i,]
         }
         mean.sample.error[i,] <- apply(sample.error, 2, mean, na.rm=TRUE)

      if (type== "632"){
         prediction.error <- matrix(ncol=1, nrow=length.cplx)
         for (i in 1:length.cplx){
            prediction.error[i,] <- (1 - 0.632)*full.apparent[i,] + 0.632*mean.sample.error[i,]
         }
      }
      if (type== "632p"){
         stop("Type ''632p'' only available if no-information error available.")
      }
      if (type == "apparent"){
         prediction.error <- peperrobject$full.apparent
      }
      if (type == "NoInf"){
         prediction.error <- peperrobject$noinf.error
      }
      if (type == "resample"){
         prediction.error <- mean.sample.error
      }
      if (type == "nullmodel"){
         prediction.error <- peperrobject$null.model
      }
   }
}
   prediction.error
}

