entropy <- function(x) {
  x = x/sum(x)
  prod((1/x)^x)
}

columnEntropy <- function(m){
  m = methods::as(m, 'dgTMatrix')
  m@x = m@x / Matrix::colSums(m)[m@j+1]
  m@x = (1/m@x)^m@x
  tapply(m@x, m@j, prod)
}

# get 1-day halflife decay
exp.decay <- function(r, decay_constant=0.6931472, halflife=NULL) {
  if(!is.null(halflife)) decay_constant = -(log(0.5) / halflife) 
  exp(-r*decay_constant)
}

matrix.autocor <- function(m, rows_lag=1){
  ## a sparse matrix solution for calculating the autocorrelation for each column vector
  ## (which, naturally, makes sense only if column vectors are time-series)
  m_lag = m[1:(nrow(m) - rows_lag),]
  if(rows_lag > 1){
    for(i in 2:rows_lag){
      m_lag = m_lag + m[i:(nrow(m) - rows_lag + (i-1)),]
    }
  } 
  m = m[(rows_lag+1):nrow(m),]
 
  mean_prod = Matrix::colMeans(m) * Matrix::colMeans(m_lag) * nrow(m)
  ss = Matrix::colSums(m * m_lag) - mean_prod
  
  l_mean_m = Matrix::colMeans(m)^2 * nrow(m)
  l_mean_mlag = Matrix::colMeans(m_lag)^2 * nrow(m_lag)
  l = sqrt(Matrix::colSums(m^2) - l_mean_m) * sqrt(Matrix::colSums(m_lag^2) - l_mean_mlag)
  cors = ss / l
  cors[is.na(cors)] = 0
  cors  
}