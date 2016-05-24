chiM <-
function (data,alpha=0.05){
        p <- dim(data)[2]
        discredata <- data
        cutp <- list()
          for(i in 1:(p-1)){
                 val <- value(i,data,alpha)
                 cutp[[i]] <- val$cuts
                 discredata[,i] <- val$disc[,i]
        }
        return(list(cutp=cutp,Disc.data=discredata))
}

