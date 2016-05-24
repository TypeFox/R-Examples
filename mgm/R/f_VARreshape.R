



f_VARreshape <- function(data, 
                         lag=1) 
{
  n_col <- ncol(data)
  n_row <- nrow(data)
  d_out <- matrix(NA, nrow = n_row, ncol=n_col*2)
  
  d_out[,1:n_col] <- data
  d_out[-n_row,(n_col+1):(n_col*2)] <- data[-1,]
  
  d_out <- d_out[-n_row, ] # one case is omitted, because the first time point cannot be predicted
  
  return(d_out)

}