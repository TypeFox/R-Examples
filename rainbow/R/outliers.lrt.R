`outliers.lrt` <- function(data, nb = 200, suav = 0.05, trim = 0.10,...){
     functions = t(data$y)
     n <- dim(functions)[1]
     m <- dim(functions)[2]
     out.thres.lrt <- outliers.thres.lrt(data, nb = nb, suav = suav, trim = trim,...)
     hay <- 1
     outliers <- c()
     valor.estadistico <- c()
     nout <- 0
     ngood <- n - nout
     functionsgood <- functions
     functionsgood2 = fts(1:dim(functionsgood)[2], t(functionsgood))
     while(hay == 1){
           aux <- c()
           auxmean <- as.vector(func.trim.mode(functionsgood2, trim = trim,...))
           auxdt <- sqrt(as.vector(func.trimvar.mode(functionsgood2, trim = trim,...)))
           for(j in 1:ngood){
               aux[j] <- metri.p(functionsgood[j,] / auxdt, auxmean / auxdt, ...)
           }
           maximo <- as.numeric(max(aux))
           fecha <- as.numeric(row.names(functionsgood)[maximo == aux])
           elim <- which(maximo == aux)
           if(maximo > out.thres.lrt){
              functionsgood <- functionsgood[-elim, ]
              outliers <- c(outliers, fecha)
              valor.estadistico <- c(valor.estadistico, maximo)
              nout <- nout + 1
              ngood <- n - nout
           }
           if(maximo < out.thres.lrt){
              hay <- 0
           }                    
    }
    return(list("outliers" = outliers, "depth.out" = valor.estadistico, "cutoff" = out.thres.lrt))
}

