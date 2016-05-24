ndlCrossvalidate <- function(formula, data, frequency=NA, k=10, folds=NULL, ...)
{ call <- match.call()
  N <- NROW(data)

  if(class(formula)!="formula")
    if(is.na(formula))
      { if(!all(colnames(data) %in% c("Frequency","Cues","Outcomes")))
          stop("Data does not have proper structure with three columns: Frequency, Cues, Outcomes.\n")
      }
    else
      stop("Incorrect specification of argument 'formula'.\n")

  if(is.null(folds))
    { if(is.na(frequency[1]))
        { n.data=N
          sampling.space <- 1:N
        }
      else
        { if(is.character(frequency) & length(frequency)==1)
            { n.data <- sum(data[,frequency])
              sampling.space <- rep(1:N, times=data[,frequency])
            }
          if(is.numeric(frequency) & length(frequency)==N)
            { n.data <- sum(frequency)
              sampling.space <- rep(1:N, times=frequency)
            }
        }
      order.random <- sample(sampling.space, n.data, FALSE)
      if(n.data/k != round(n.data/k))
        order.random <- c(order.random, order.random[1:(ceiling(n.data/k)*k-n.data)])
      folds <- lapply(1:k, function(x) matrix(order.random,k,byrow=T)[x,])
    }
  else
    { k = length(folds)
      if(sd(sapply(folds, length))!=0)
        stop("Folds do not have same length")
      sum.folds = sum(sapply(folds, length))
      if(is.na(frequency[1]))
        n.data=NROW(data)
      else
        { if(is.character(frequency) & length(frequency)==1)
            n.data=sum(data[,frequency])
          if(is.numeric(frequency) & length(frequency)==N)
            n.data=sum(frequency)
        }
      if(sum.folds < floor(n.data/k)*k | sum.folds > ceiling(n.data/k)*k)
        stop("Folds together do not match size of data")
    }
  n.test = mean(sapply(folds, length))
  n.train = n.data - n.test

# One might want to finetune more exactly what happens if n.data/k is not not a natural number

options(ndlCrossvalidate.counter=1)

fits <- lapply(folds,
   function(test)
      { 
        cat(paste("[",getOption("ndlCrossvalidate.counter"),"]",sep=""))
        options(ndlCrossvalidate.counter=getOption("ndlCrossvalidate.counter")+1)

        if(class(formula)=="formula")
          cuesOutcomes.teach <- ndlCuesOutcomes(formula, data[-test,])
        else
          cuesOutcomes.teach <- data[-test,]

        weightMatrix.teach = estimateWeights(cuesOutcomes.teach, ...)

        if(class(formula)=="formula")
          cuesOutcomes.test <- ndlCuesOutcomes(formula, data[test,])
        else
          cuesOutcomes.test <- data[test,]
        activationMatrix.test = estimateActivations(cuesOutcomes.test, weightMatrix.teach, ...)$activationMatrix 

        ndl.test <- list(activationMatrix = activationMatrix.test,  weightMatrix = weightMatrix.teach, cuesOutcomes = cuesOutcomes.test)
        test.result <- ndlStatistics(ndl.test, ...)
        names(test.result)[which(names(test.result)=="n.data")]="n.test"

        return(test.result)
      })

     options(ndlCrossvalidate.counter=NULL)
     cat("\n")

result = list(call=call, formula = formula, fits = fits, k = k, n.total=n.data, n.train=n.train, n.test=n.test, folds = folds)
class(result) <- "ndlCrossvalidate"

return(result)

}

