`outliers.depth.trim` <- function(data, dfunc = depth.RP, nb = 200, suav = 0.05, trim,...){
       functions = t(data$y)
       n <- dim(functions)[1]
       m <- dim(functions)[2]
       cutoff <- median(quantile.outliers.trim(data, dfunc = dfunc, nb = nb,
                        suav = suav, trim = trim,...))
       hay <- 1
       outliers <- c()
       prof.out <- c()
       functionsgood <- functions
       rownames(functionsgood) = 1:n
       functionsgood2 = fts(1:dim(functionsgood)[2], t(functionsgood))
       while(hay == 1){
             d = dfunc(functionsgood2, trim = trim, ...)$prof
             if(is.null(outliers)){
                dtotal <- d
             }
             fecha <- as.numeric(rownames(functionsgood)[d < cutoff])
             elim <- which(d < cutoff)
             if(length(elim) > 0){
                prof.out <- c(prof.out, d[d<cutoff])
                functionsgood <- functionsgood[-elim, ]
                outliers <- c(outliers, fecha)
             }
             if(length(elim) == 0||length(outliers) > n / 5){
                hay <- 0
             }                    
       }
       return(list("outliers" = outliers, "cutoff" = cutoff, "depth.total" = dtotal, "depth.out" = prof.out))
}

