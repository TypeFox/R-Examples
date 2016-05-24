`outliers.thres.lrt` <- function(data, nb = 200, suav = 0.05, trim = 0.10,...){
  functions = t(data$y)
  n <- dim(functions)[1]
  m <- dim(functions)[2]
  if(is.null(n) && is.null(m)) 
     stop("I do not have a matrix")
     maximos <- c()
     aux <- c()
     resample.boot <- matrix(NA, nrow = n, ncol = m)
     for(i in 1:nb){
         bsample <- functions[sample(1:n, size = n, replace = T),]
         if(suav > 0){
            bsample <- bsample + mvrnorm(n = n, rep(0,m), var(functions) * suav)}
            bsample2 = fts(1:dim(bsample)[2], t(bsample))
            auxmean <- func.trim.mode(bsample2, trim = trim, ...)
            auxdt <- sqrt(as.vector(func.trimvar.mode(bsample2, trim = trim, ...)))
            d <- matrix(NA, nrow = n, ncol = m)
            for(j in 1:m){
                d[,j] <- 1 - abs(.5 - rank(bsample[,j], ties.method = "average") / n)
            }
            ans <- apply(d, 1, sum)
            rid <- rank(ans, ties.method = "first")
            bsample.trim <- bsample[rid >= floor(trim * n),]
            for(j in 1:(n - floor(trim * n))){
                aux[j] <- metri.p(bsample.trim[j,] / auxdt, auxmean / auxdt, ...)
            }
            maximos[i] <- as.numeric(max(aux))
     }
     as.numeric(quantile(maximos, .99))
}

