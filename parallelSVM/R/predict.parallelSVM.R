#' @importFrom parallel detectCores
#' @importFrom foreach foreach %do%
#' @export
predict.parallelSVM <- function(object, newdata, compute = TRUE, probability = FALSE, ...){
  # Predict for each model stored in object
  
  # Use the amount of cores the model was build with
  registerCores(length(object))
  
  # Get the function call and its arguments
  fcall <- match.call()
  
  # Remove compute from the function call
  fcall$compute <- NULL
  
  # remember the call used
  call <- gsub('.parallelSVM','',deparse(fcall))
  call <- gsub('    ','',call)
  
  # Construct the new call expression
  fcall[[1]] <- predict
  
  # Create copies with the correct data
  function_call <- list()
  for (i in 1:length(object)){
    function_call[[i]]              <- fcall
    function_call[[i]]$object       <- object[[i]]
  }
  
  # parallel prediction
  predictDataSvm <- foreach(i = 1:length(object))  %do% {
    # Do the call
    eval(function_call[[i]], parent.frame())
  }
  
  # Calculate the average
  if (compute){
    # the factors
    result <- as.factor(apply(data.frame(predictDataSvm),
                              1,function(vec){names(sort(table(vec),decreasing=TRUE))[1]}))
    
    if (probability){
      # probability attributes
      column_names          <- colnames(attr(predictDataSvm[[1]],"probabilities"))
      row_names             <- rownames(attr(predictDataSvm[[1]],"probabilities"))
      attr_matrix           <- matrix(NA,length(predictDataSvm[[1]]),length(column_names))
      colnames(attr_matrix) <-column_names
      
      for (this_row in column_names){
        temp_frame <- matrix(NA,length(predictDataSvm[[1]]),length(predictDataSvm))
        for (j in 1:length(predictDataSvm)){
          temp_frame[,j] <- attr(predictDataSvm[[j]],"probabilities")[,this_row]
        }
        attr_matrix[,this_row] <- apply(temp_frame,1,mean)
      }
      rownames(attr_matrix)        <- row_names
      attr(result,"probabilities") <- attr_matrix
      return(result)
    } else{
        return(result)
    }
  } else{
    return(predictDataSvm)
  }
}

