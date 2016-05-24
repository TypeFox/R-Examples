`outliers.depth.pond` <- function(data, dfunc = dfunc, nb = 200, suav = 0.05,...){
  functions = t(data$y)
  n <- dim(functions)[1]
  m <- dim(functions)[2]
  if(is.null(n) && is.null(m)) 
     stop("I do not have a matrix")
     cutoff <- median(quantile.outliers.pond(data, dfunc = dfunc, nb = nb, suav = suav,...))
     hay <- 1
     outliers <- c()
     prof.out <- c()
     functionsgood <- functions
     row.names(functionsgood) = 1:n
     functionsgood2 = fts(1:dim(functionsgood)[2], t(functionsgood))
     while(hay == 1){
           d = dfunc(functionsgood2, trim = 0.25, ...)$prof
           if(is.null(outliers)){
              dtotal <- d
           }
           fecha <- as.numeric(row.names(functionsgood)[d < cutoff])
           elim <- which(d < cutoff)
           if(length(elim) > 0){
              prof.out <- c(prof.out, d[d < cutoff])
              functionsgood <- functionsgood[-elim, ]
              outliers <- c(outliers, fecha)
           }
           if(length(elim) == 0||length(outliers) > n/5){
              hay <- 0
           }                    
     }
  return(list("outliers" = outliers, "cutoff" = cutoff, "depth.total" = dtotal, 
         "depth.out" = prof.out))
}

