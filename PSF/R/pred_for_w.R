#' It predict the values for given data, Window size (W) and Cluster size (k).
#'
#' Takes Time Series data and Window size as Input
#' @param data_in as Input data, in any format (data matrix data frame or vector). All variables should be numeric and NA values will get removed while execution.
#' @param w as window size (Can be obtained with function "optimum_w")
#' @param k as cluster size for Kmeans (Can be obtained with function "optimum_k")
#' @param next_val as Integer number. It states the number of predicted values to be obtained.
#' @return deno_prediction as predicted value
#' @export
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)      #data_in = x
#' ## Consider, w = 3, K = 4, next_val = 5
#' pred_for_w(x, 3, 4, 5)

#----------------------------------------------------------------------------------------------
# a Function "pred_for_w" to predict the future values for given input data, w and number of future values
#==============================================================================================
pred_for_w <- function(data_in, w, k, next_val)
{
  data_x <- data_in
  prediction <- NULL
  w1 <- w
  set.seed(1)
  data_in <- data_norm(data_in)$x
  result<- kmeans(data_in,k,nstart=1)      # Asumming k = 4
  data_in <- result$cluster

  # Loop for multiple number of predictions
  while(next_val != 0)
  {
    # Window Selection
    y <- length(data_in)
    # n1 contains last W digits
    n1 <- NULL
    while(w != 0)
    {
      n1 <- append(n1, data_in[y-w+1])
      w <- w - 1
    }

    #Searching
    #sapply(x, match, n1[1], nomatch=0)------- it is showing that digit is present
    #n2 will contain list of all next possible values (need to avg or nedian it)
    #count1 will contain number of such sequence
    count1 <- 0
    n3 <- NULL
    z <- w1
    for (i in 1:y)
    {
      v <- 1
      j <- i
      while(w1 != 0)
      {
        if(!(is.na(data_in[j] == n1[v])) && data_in[j] == n1[v])
        {
          if(z == v)
          {
            n3 <- append(n3, data_in[j+1])
            count1 <- count1 + 1
          }
          else
          {
            j <- j + 1   #original place
            v <- v + 1
          }
        }
        w1 <- w1 - 1
      }
      w1 <- z
      w <- z
    }


    n3 <- na.omit(n3)
    b <- round(mean(n3))         # b is the next predicted value
    #prediction <- append(prediction,b)
    prediction <- append(prediction,result$centers[b])
    data_in <- append(data_in, b)
    next_val <- next_val - 1
  }
  w <- w + 1

  deno_prediction <- data_denorm(prediction,data_norm(data_x)$dmax,data_norm(data_x)$dmin)

  #jk <- plot_PSF(data_x,deno_prediction)

  options(warn=-1)
  #output <- list("Predicted_Values"=deno_prediction,"Plot"=jk)
  return(deno_prediction)
  #return(output)
}
### Function (pred_for_w) ends here=============================================================
