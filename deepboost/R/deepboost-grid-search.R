#' Returns optimised parameter list for deepboost model on given data
#' @param formula A R Formula object see : ?formula
#' @param data input data.frame as training for model
#' @param k number of folds (default = 10) for cross validation optimisation
#' @param seed for random split to train / test (default 666)
#' @param logging_level print extra data while training 0 - no data, 1 - gridSearch data (default), 2 - all data
#' @details Finds optimised parameters for deepboost training.
#'  using grid search techniques over:
#'  - predefined, battle tested parameter possible values
#'  - cross validation over k folds
#' @return vector with average accuracy for chosen parameters, and a list of the best parameter combination: (accuracy, (num_iter, beta, lambda, loss_type))
#' @examples
#' deepboost.gridSearch(y ~ .,
#'  data.frame(x1=rep(c(0,0,1,1),2),x2=rep(c(0,1,0,1),2),y=factor(rep(c(0,0,0,1),2))), k=2)
#' @export
deepboost.gridSearch <- function(formula, data, k=10, seed=666, logging_level=1) {

  if (!(is.numeric(k)) || k <= 1 || !(k%%1==0))
  {
    stop("ERROR_paramter_setting : k must be >= 2 and integer (Default : 10)" )
  }

  if (!(is.numeric(logging_level)) || logging_level < 0 || logging_level > 2 || !(k%%1==0))
  {
    stop("ERROR_paramter_setting : logging_level must be integer (0 / 1 / 2) (Default : 1)" )
  }

  verbose <- ifelse(logging_level>1,TRUE,FALSE)

  num_iter_vals = c(5,10,25,50)
  beta_vals = c(2^-0, 2^-1, 2^-2, 2^-3, 2^-4, 2^-5, 2^-6)
  lambda_vals = c(0.0001, 0.005, 0.01, 0.05, 0.1, 0.5)
  loss_type_vals = c("l","e")
  dpbGrid <-  expand.grid(num_iter = num_iter_vals,
                          beta = beta_vals,
                          lambda = lambda_vals,
                          loss_type = loss_type_vals)

  set.seed(seed)

  #Randomly shuffle the data
  data<-data[sample(nrow(data)),]

  folds <- cut(seq(1,nrow(data)),breaks=k,labels=FALSE)
  best_acc <- -Inf
  avg_acc <- 0

  for(combination in 1:nrow(dpbGrid)){
    num_iter <- dpbGrid[combination,"num_iter"]
    beta <- dpbGrid[combination,"beta"]
    lambda <- dpbGrid[combination,"lambda"]
    loss_type <- as.character(dpbGrid[combination,"loss_type"])
    acc <- 0

    for(fold in 1:k){
      testIndexes <- which(folds==fold,arr.ind=TRUE)
      testData <- data[testIndexes, ]
      trainData <- data[-testIndexes, ]

      eval_model <- deepboost.formula(formula, trainData, num_iter = num_iter, beta = beta, lambda = lambda, loss_type = loss_type, verbose=verbose)
      acc <-  acc + sum(predict(eval_model, testData) == testData[,length(testData)]) / nrow(testData)
    }
    acc <- acc / k
    if(acc > best_acc){
      best_acc <- acc
      best_num_iter <- num_iter
      best_lambda <- lambda
      best_beta <- beta
      best_loss_type <- loss_type
    }
    avg_acc <- avg_acc + acc

  }
  avg_acc <- avg_acc / nrow(dpbGrid)

  if(logging_level > 0)
  {
    print(paste0("average accuracy : ", avg_acc))
    print(paste0("accuracy: ", best_acc, ", num_iter: ", best_num_iter, ", beta: ", best_beta, ", lambda: ", best_lambda, ", loss_type: ", best_loss_type))
  }

  RET <-
    c(avg_acc,
      list(best_num_iter,
           best_lambda,
           best_beta,
           best_loss_type))

  return(RET)
}
