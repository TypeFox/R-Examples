BigQuic <- function(X = NULL, inputFileName = NULL, outputFileName = NULL, lambda = 0.5, numthreads = 4, maxit = 5, epsilon = 1e-3, k = 0, memory_size = 8000, verbose = 0, isnormalized = 1, seed = NULL, use_ram = FALSE){
  outputFileNames <- vector(length = length(lambda))
  outFlag <- TRUE
  inFlag <- TRUE
  if (use_ram)
  {
    #precMatrices <- vector(length = length(lambda))
    precMatrices <- vector(length = 0)    
  }
  if (!is.null(seed) && length(seed) != 0)
  {
    set.seed(seed)
  }
  if (!is.null(X))
  {
    inFlag <- FALSE
    #WRITE INPUT FILE FOR BIGQUIC
    inputFileName <- tempfile(pattern = "BigQuic_input_matrix", fileext = ".Bmat")
    write.table(x = cbind(dim(X)[2], dim(X)[1]), file = inputFileName, row.names = FALSE, col.names = FALSE)
    write.table(x = X, file = inputFileName, append = TRUE, row.names = FALSE, col.names = FALSE)
  }
  if (!is.null(outputFileName))
  {
    tempFileName <- outputFileName
  }
  else
  {
    tempFileName <- inputFileName
    outFlag <- FALSE
  }
  for (i in 1:length(lambda))
  {
    j = i
    while (file.exists(paste0(tempFileName, ".", j, ".output")))
    {
      j <- j + 1
    }
    outputFileName <- paste0(tempFileName, ".", j, ".output")
    #Check input file is at least kinda valid
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
    BigQuicHelper(argvPasser = c("-l", lambda[i], "-n", numthreads, "-t", maxit, "-e", epsilon,  "-k", k, "-m", memory_size, "-q", verbose,  "-r", isnormalized,  inputFileName, outputFileName))
    outputFileNames[i] <- outputFileName
    if (use_ram)
    {
      M <- read.table(file = outputFileName, skip = 1, )
      #########GET FROM FILE DIM(X)
      precMatrices <- c(precMatrices, sparseMatrix(i = M[,1], j = M[,2], x = M[,3], dims = c(format_Check[1],format_Check[1]), symmetric = FALSE))  # Inverse/thetahat
    }
    if (outFlag == FALSE && use_ram == TRUE)
    {
      file.remove(outputFileName)
    }
  }
  if(use_ram && outFlag == FALSE)
  {
    output <- BigQuic_object_builder$new(precision_matrices = precMatrices, inputFileName = inputFileName, lambda = lambda, numthreads = numthreads, maxit = maxit, epsilon = epsilon, k = k, memory_size = memory_size, verbose = verbose, isnormalized = isnormalized, use_ram = use_ram, clean = FALSE, inFlag = inFlag, outFlag = outFlag)
  } else if (use_ram && outFlag == TRUE)
  {
    output <- BigQuic_object_builder$new(precision_matrices = precMatrices, output_file_names = outputFileNames, inputFileName = inputFileName, lambda = lambda, numthreads = numthreads, maxit = maxit, epsilon = epsilon, k = k, memory_size = memory_size, verbose = verbose, isnormalized = isnormalized, use_ram = use_ram, clean = FALSE, inFlag = inFlag, outFlag = outFlag)    
  } else 
  {
    output <- BigQuic_object_builder$new(output_file_names = outputFileNames, inputFileName = inputFileName, lambda = lambda, numthreads = numthreads, maxit = maxit, epsilon = epsilon, k = k, memory_size = memory_size, verbose = verbose, isnormalized = isnormalized, use_ram = use_ram, clean = FALSE, inFlag = inFlag, outFlag = outFlag)
  }
  if("AsIs" %in% class(X)){
    class(X) <- class(X)[-match("AsIs", class(X))]
    if(class(X) != "matrix"){
      X <- as.matrix(X)
      if(class(X) != "matrix")
      {
        stop("X is not a matrix, nor a matrix protected with AsIs")
      }
    }
  }

  output$setX(X)
  output$setSeed(seed)
  return(output)
}


BigQuic_object_builder <- setRefClass(Class = "BigQuic_object", 
            fields = list(precision_matrices = "list", X = "matrix", inputFileName = "character", lambda = "numeric", 
                                      numthreads = "numeric", maxit = "numeric", epsilon = "numeric", k = "numeric", memory_size = "numeric", 
                                      verbose = "numeric", isnormalized = "numeric", seed = "numeric", use_ram = "logical", 
                                      clean = "logical", output_file_names = "character", opt.lambda = "numeric", inFlag = "logical", 
                                      outFlag = "logical"), 
                          methods = list(cleanFiles = function(verbose = FALSE)
                            {
                              if (clean == FALSE)
                              {
                                if (outFlag == TRUE || use_ram == FALSE)
                                {
                                  for (i in 1:length(output_file_names))
                                  {
                                    file.remove(output_file_names[i])
                                    if (verbose)
                                    {
                                      print(c("Deleted file: ", output_file_names[i]))
                                    }
                                  }
                                }
                                if (inFlag == FALSE)
                                {
                                  file.remove(inputFileName)
                                  if (verbose)
                                  {
                                    print(c("Deleted file: ", inputFileName))
                                  }
                                }
                                clean <<- TRUE
                              }
                              else
                              {
                                print("Files were already cleaned up.  ")
                              }
                            },
                          setSeed = function(inputSeed)
                          {
                            if (!is.null(inputSeed))
                            {
                              seed <<- inputSeed
                            }
                          }, 
                          setOptLambda = function(optLambda)
                          {
                            if(!is.null(optLambda))
                            {
                              opt.lambda <<- optLambda
                            }
                          },
                          setX = function(inputX)
                          {
                            if(!is.null(inputX))
                            {
                              X <<- inputX
                            }
                          }, 
                          finalize = function()
                            {
                              cleanFiles()
                            }
                          ))

generate_sample <- function(n = 200, p = 150, seed = NULL)
{
  if (!is.null(seed))
  {
    set.seed(seed)
  } else
  {
    set.seed('1') 
  }
  
  X <- rbinom(p*n,1,prob=0.15);
  dim(X) <- c(n,p);
  X <- X %*% diag(1+9*runif(p))
  return(X)
}

