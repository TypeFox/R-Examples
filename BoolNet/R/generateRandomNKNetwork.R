generateRandomNKNetwork <- function(n,k,topology=c("fixed","homogeneous","scale_free"),
          linkage=c("uniform","lattice"),
          functionGeneration=c("uniform","biased"),
          validationFunction, failureIterations=10000,
          simplify=FALSE, noIrrelevantGenes=TRUE, 
          readableFunctions=FALSE,
          d_lattice=1, zeroBias=0.5, gamma=2.5, approx_cutoff=100)
{
  k_i_vec <- switch(match.arg(topology),
      fixed = {
          if (length(k) == n)
            k
          else
          if (length(k) == 1)
            rep(k,n)
          else
            stop("k must have 1 or n element(s)!")
        },
      homogeneous = round(rpois(n,k)),
      scale_free = rzeta(n,k,gamma=gamma,approx_cutoff=approx_cutoff),
      stop("'topology' must be one of \"fixed\",\"homogeneous\"")
      )
      
  k_i_vec[k_i_vec > n] <- n
  
  if ((zeroBias == 0 | zeroBias == 1) && noIrrelevantGenes && any(k_i_vec > 0))
    stop("If setting 'zeroBias' to 0 or 1, you must set 'noIrrelevantGenes' to FALSE!")
  
  geneNames <- paste("Gene",seq_len(n),sep="")
  
  if (!missing(validationFunction) && !is.null(validationFunction))
    validationFunction <- match.fun(validationFunction)
  else
    validationFunction <- NULL    
      
  interactions <- mapply(function(i,k_i)
      {
        if (k_i == 0)
        {
          genes <- 0
          func <- round(runif(min=0,max=1,n=1))
        }
        else
        {
          table <- allcombn(2,k_i) - 1
          genes <- switch(match.arg(linkage,c("uniform","lattice")),
            uniform = sample(seq_len(n),k_i,replace=FALSE),
            lattice = {
                region <- c(max(1,round(i - k_i*d_lattice)):max(1,i-1),
                      min(n,i+1):min(n,round(i + k_i*d_lattice)))
                
                sample(region,k_i,replace=FALSE)
              },
            stop("'linkage' must be one of \"uniform\",\"lattice\""))

            containsIrrelevant <- TRUE
            validFunction <- TRUE
            
            counter <- 0
            while((!validFunction) || containsIrrelevant)
            {
                tryCatch(
                {
                  functionGeneration <- match.fun(functionGeneration)
                },
                error = function(e){}
                )
                if (is.function(functionGeneration))
                  func <- functionGeneration(input=genes)
                else
                  func <- switch(match.arg(functionGeneration,c("uniform","biased")),
                      uniform = round(runif(min=0,max=1,n=2^k_i)),
                      biased = as.integer(runif(min=0,max=1,n=2^k_i) > zeroBias),
                      stop("'functionGeneration' must be one of \"uniform\",\"biased\" or a function!"))

                if (noIrrelevantGenes)
                {
                  dropGenes <- apply(table,2,function(gene)
                  # determine all genes that have no influence on the results,
                  # i.e. the result column is equal for 0 values and 1 values
                          {
                            (identical(func[gene==1],
                                  func[gene==0]))
                          })                                    
                  containsIrrelevant <- (sum(dropGenes) > 0)                  
                }
                else
                  containsIrrelevant <- FALSE
                
                if (!is.null(validationFunction))
                  validFunction <- validationFunction(table, func)
                  
                counter <- counter + 1
                if (counter > failureIterations)
                  stop(paste("Could not find a transition function matching the restrictions of",
                             "\"validationFunction\" in", failureIterations, " runs.",
                             "Change \"validationFunction\" or increase \"failureIterations\"!"))
            }
        }
        return(list(input=genes,func=func,
          expression=getInteractionString(readableFunctions,
                         func,geneNames[genes])))
      },seq_len(length(k_i_vec)),k_i_vec,SIMPLIFY=FALSE)
  fixed <- sapply(interactions,function(i)
      {
        if (i$input[1] == 0)
          i$func[1]
        else
          -1
      })

  names(fixed) <- geneNames
  names(interactions) <- geneNames
  
  net <- list(genes=geneNames,interactions=interactions,fixed=fixed)
  class(net) <- "BooleanNetwork"
  if (simplify)
    net <- simplifyNetwork(net,readableFunctions)
  return(net)  
}

generateCanalyzing <- function(input) 
{
 k <- length(input)
 table <- allcombn(2, k) - 1
 canalyzingInput <- sample(1:k,size=1)
 res <- round(runif(min = 0, 
                    max = 1, n = 2^k))
 
 res[table[,canalyzingInput] == sample(c(0,1),size=1)] <- sample(c(0,1), size=1)
 return(res)
}

generateNestedCanalyzing <- function(input) 
{
 k <- length(input)
 table <- allcombn(2, k) - 1
 remainingIndices <- 1:(2^k)
 res <- rep(sample(c(0,1), size=1),2^k)
 
 for (canalyzingInput in sample(1:k,size=k,replace=FALSE))
 {
   newIndices <- which(table[remainingIndices,canalyzingInput] == sample(c(0,1),size=1))
   res[remainingIndices[newIndices]] <- sample(c(0,1), size=1)
   remainingIndices <- setdiff(remainingIndices, remainingIndices[newIndices])
 }
 return(res)
}
