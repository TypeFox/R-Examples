BigQuic.select = function(BigQuic_result = NULL, stars.thresh = 0.1, 
                          stars.subsample.ratio = NULL, rep.num = 20, 
                          verbose = TRUE, verbose2 = 0)
{
  if(length(BigQuic_result$lambda) == 1){
    BigQuic_result$setOptLambda(BigQuic_result$lambda)
    return(BigQuic_result)    
  }
  X <- read.table(file = BigQuic_result$inputFileName, skip = 1, )
  gcinfo(FALSE)
	n = nrow(X)
	d = ncol(X)
	nlambda = length(BigQuic_result$lambda)

	if(is.null(stars.subsample.ratio))
	{
		if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
		if(n<=144) stars.subsample.ratio = 0.8
	} 
	
	merge <- list()
	for(i in 1:nlambda)
  {
    merge[[i]] <- Matrix(0,d,d)
	}

	for(i in 1:rep.num)
	{
		if(verbose)
		{
		  mes <- paste(c("Conducting Subsampling....in progress:", 
                     floor(100*i/rep.num), "%"), collapse="")
			cat(mes, "\r")
      flush.console()
		}
		ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)		
    #HERE RUN BIGQUIC lambda times and store the resulting matrices as a list in tmp
    path <- list()
    for (j in 1:nlambda)
    {
        temp <- BigQuic(X = as.matrix(X[ind.sample,]), 
                        lambda = BigQuic_result$lambda[j], 
                        numthreads = BigQuic_result$numthreads, 
                        maxit = BigQuic_result$maxit, 
                        epsilon = BigQuic_result$epsilon, 
                        k = BigQuic_result$k, 
                        memory_size = BigQuic_result$memory_size, 
                        verbose = verbose2, 
                        isnormalized = BigQuic_result$isnormalized, 
                        seed = BigQuic_result$seed, use_ram = TRUE)
        path[[j]] <- temp$precision_matrices[[1]]
        #Probably no longer needed because they are cleaned by the garbage 
        #collection now.  
        #temp$cleanFiles(verbose = FALSE)
    }
    			
		for(k in 1:nlambda)
		{
			merge[[k]] = merge[[k]] + path[[k]]
		}
		rm(ind.sample,path)
		gc()  #Should there really be an explicit call to the garbage collector here?
	}
		
	if(verbose)
	{
			mes = "Conducting Subsampling....done.                 "
  		cat(mes, "\r")
  		cat("\n")
  		flush.console()
	}
        
  variability = rep(0,nlambda)
	for(i in 1:nlambda)
  {
		merge[[i]] = merge[[i]]/rep.num
  	variability[i] = 4*sum(merge[[i]]*(1-merge[[i]]))/(d*(d-1))
  }
    
  #find an index, use max to take 1 if index doesn't make sense
  #The index is where the smallest lambda is that has maximum low enough variability
  #opt.index = max(which.max(variability <= stars.thresh)[1]-1,1)
  variability <- rev(variability)
  monotonocity_index <- 0
  for (i in 1:(nlambda - 1))
  {
    if (variability[i] > variability[i + 1])
    {
      monotonocity_index <- i
      break
    }
  }
  variability = variability[1:monotonocity_index] # monotonocity and hence sparsity
  variability = variability[variability <= stars.thresh]     # reproducibility
  if (is.null(variability))
  {
    #Should report to the user that no lambdas are suitable.  
    return("No lambda's are suitable for this method because they have too much variability and are therefore not reproducible!")
  } else
  {
    opt.index = which.max(variability)
  }

  opt.lambda = BigQuic_result$lambda[nlambda - opt.index + 1]

  BigQuic_result$setOptLambda(opt.lambda)
#   if (BigQuic_result$use_ram == TRUE)
#   {
#     BigQuic_result <- BigQuic_result$precision_matrices[[which(BigQuic_result$lambda == BigQuic_result$opt.lambda)]]
#   }
#   else
#   {
#     #BigQuic_result <- 
#   }
  return(BigQuic_result)
}

plot.BigQuic_object = function(x, ...){

  if(x$use_ram == FALSE){
    format_Check <- read.table(file = inputFileName, nrows = 1)
    if (!is.integer(format_Check[[1]]) || !is.integer(format_Check[[1]]))
    {
      stop("The file is not formatted correctly for BigQuic, the first line 
           should be p (the number of attributes) then n (the number of 
           samples).  Then the rest of the file should contain the matrix, 
           e.g. 
           4 2
           1 2 3 4
           4 3 2 1")
    }
    M <- read.table(file = x$output_file_names[which(x$lambda == x$opt.lambda, arr.ind = TRUE)], skip = 1, )
    x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]] <- sparseMatrix(i = M[,1], j = M[,2], x = M[,3], dims = c(format_Check[1],format_Check[1]), symmetric = FALSE)
  }
  
  col_vec1 <- unlist(summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[3])
  col_vec1[col_vec1 < 0] <- 0
  #col_vec1 excludes color from the diagonal as not interesting
  col_vec1[summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[1] == summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[2]] <- 0
  
  col_vec2 <- unlist(summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[3])
  col_vec2[col_vec2 > 0] <- 0
  col_vec2 <- col_vec2*-1
  #col_vec2 excludes color from the diagonal as not interesting
  col_vec2[summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[1] == summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[2]] <- 0
  
  #Make anything green very green to see it better.  
  #col_vec2[col_vec2 != 0] <- 1
  
  if(max(col_vec1) > 0){
    col_vec1 <- col_vec1/max(col_vec1)
  } else{
    print("No positive associations.")
  }
  if(max(col_vec2) > 0){
    col_vec2 <- col_vec2/max(col_vec2)
  } else{
    print("No negative associations.")
  }
  
  plot(as.vector(summary(x$precision_matrices[[which(x$lambda == x$opt.lambda, arr.ind = TRUE)]])[1:2]), col=rgb(col_vec1,col_vec2,0,1))
}