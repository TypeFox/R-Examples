#' Finds optimum value of Window Size (w)
#'
#' Takes time series data and number of values to predict as input
#' @param data_in as Input data, in any format (data matrix data frame or vector). All variables should be numeric and NA values will get removed while execution.
#' @param next_val as Integer number. It states the number of predicted values to be obtained.
#' @return Optimum_W as optimum value of Window size
#' @return RMSE_val as RMSE values corresponding to each window size (W)
#' @return prediction as "next_val" predicted values corresponding to each window size (W)
#' @export
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)    # data_in = x
#' y <- optimum_w (x, 2)        # next_val = 2
#' y$Optimum_W
#' y$RMSE_Values
#' y$Prediction
#' y
#' ## Note: Sometime this function may display a error,
#' ## in such case, kindly minimize the value of "next_val" and proceed.

#-----------------------------------------------------------------------------------------------
#### Function (optimum_w) for finding optimum size of Window====================================
####============================================================================================

optimum_w <- function(data_in,next_val)
{
  k <- optimum_k(data_in)
  #k <- 4
  set.seed(1)
  #----- Delete Capture last "next_val" digits from Data compare with forecasted walues with various "win" sizes and Calculate RMSE/MAE------##
  #data_in <- data.frame(data_in)
  if(!is.vector(data_in))
  {
    data_in <- data_in[, 1]
  }

  win1 <- 1
  Original_Data <- NULL
  y <- max(c(length(data_in),nrow(data_in)))
  y1 <- next_val
  data_in_cut <- data_in[1:(y - y1)]
  while(y1 != 0)
  {
    Original_Data  <- append(Original_Data,data_in[y-y1+1])
    y1 <- y1 - 1
    na.rm = TRUE
  }
  ## "Original_Data" is actual data string
  ## "nb"="W_size" are predicted data string

  W_size <- pred_for_w(data_in_cut, win1, k, next_val)[1]
  W_size <- unlist(W_size)
  sv <- data.frame(Original_Data)
  sv[,1] <- data.frame(Original_Data)

  ## This while loop will generate a table for predicted values

  #while(!is.na(nb[1]))
  #while(!is.na(sv[1,win1]))
  while(!(is.na(W_size[1] == 0)) && W_size[1] != 0)
    #while(!(is.na(nb[1])) && nb[1] != 0)
  {
    W_size <- pred_for_w(data_in_cut, win1, k, next_val)[1]
    W_size <- unlist(W_size)
    #dg <- (unlist(W_size[1]))
    # sv[,win1 + 1] <- data.frame(as.numeric(unlist(W_size)))
    #sv[,win1 + 1] <- (as.double(unlist(W_size[1])))
    sv[,win1 + 1] <- data.frame(W_size)
    #sv[,win1 + 1] <- (dg)
    win1 <- win1 + 1
  }
  sv[,2] <- NA
  sv <- sv[ , colSums(is.na(sv)) == 0]

  #sv <- sv[ , colSums(sv[,2] == 0]

  rmse_val <- NULL
  mae_val <- NULL
  #for(i in 3:ncol(sv)-1 )
  for(i in 2:ncol(sv) )
  {
    rmse_t <- rmse(sv[,i]-sv[,1])
    rmse_val <- append(rmse_val,rmse_t)

    mae_t <- mae(sv[,i]-sv[,1])
    mae_val <- append(mae_val,mae_t)
  }
  rmse_val
  mae_val

  #perfect_w <- which.min(rmse_val)
  perfect_w <- max(which(rmse_val == min(rmse_val)))
  output <- list("Optimum_W"= perfect_w,"RMSE_Values"=rmse_val,"Prediction"=sv)
  #return(perfect_w)
  return(output)
}
#### End of the Function (optimum_w)===========================================================
