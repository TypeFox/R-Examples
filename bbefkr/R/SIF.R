SIF <-
function(BAND_MATRIX, NUM_ITERATIONS, NUM_BATCH) 
{
    size_batch = NUM_ITERATIONS/NUM_BATCH
    h_mean <- colMeans(BAND_MATRIX)
    h_mean <-  matrix(rep(h_mean,times= dim(BAND_MATRIX)[1]), nrow = dim(BAND_MATRIX)[1], 
                      ncol = dim(BAND_MATRIX)[2], byrow = TRUE)
    sigma_square_tilde <- (1/(NUM_ITERATIONS - 1)) * colSums((BAND_MATRIX - h_mean)^2) 
    sum_par_mean = rep(0, times = dim(BAND_MATRIX)[2])        
    for (i in 1: NUM_BATCH)
    {                
        sum_par_mean = sum_par_mean + (colMeans(BAND_MATRIX[((i-1) * size_batch + 1):(size_batch * i),]) - h_mean[1,] )^2            
    }
    sigma_square_hat = size_batch/(NUM_BATCH - 1) * sum_par_mean
    var_h = sqrt(sum_par_mean/(NUM_BATCH^2 - NUM_BATCH))
    return(list(batch_se = round(sqrt(sigma_square_hat),4), 
                total_se = round(sqrt(sigma_square_tilde),4), 
                SIF = round(sigma_square_hat/sigma_square_tilde,4), VAR_H = round(var_h,4)))
}
