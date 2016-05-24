########## S function: random.read ##########

# For reading in random effect information

# Last changed: 28 JUN 2002

random.read <- function(random,group)
{
   rhs <- as.character(random)[-1]

   rhs <- rm.char(rhs,"\n")
   rhs <- rm.char(rhs,"\t")

   rhs <- break.string(rhs,"+")

   # Build the random design matrix
   # denoted by X^R in the Zhao et al
   # paper.

   if (rhs[1]=="1")
      X.random <- matrix(1,length(group),1)

   if (rhs[1]!="1")
      X.random <- eval(parse(text=rhs[1]))

   if (length(rhs)>1)
      for (i in 2:length(rhs))
         X.random <- cbind(X.random,eval(parse(text=rhs[i])))

   num.meas <- length(group)

   num.subjects <- length(unique(group))

   sample.sizes <- unlist(lapply(split(1:num.meas,group),length))
   csum <- cumsum(sample.sizes)

   group.inds <- list()
   group.inds[[1]] <- 1:csum[1]
   for (i in 2:num.subjects)
      group.inds[[i]] <- (csum[i-1]+1):csum[i]

   random.info <- list(X.random=X.random,group.inds=group.inds)

   return(random.info)

}

########## End of random.read ##########

