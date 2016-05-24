`pick.high.conc` <-
function(x,highest=("dilution"), sample.id=c("sample","sample.n")){

    xi <- create.ID.col(x, sample.id=sample.id) 

    ttt <- by(xi[[4]], factor(xi[[4]]$identifier, levels=unique(xi[[4]]$identifier)), function(xx) {
             
             maxConc <- max(xx[,highest])
             IDs <- xx[xx[,highest]==maxConc,"ID"]

             return(as.character(IDs))   
                
            })

    lines <- do.call("c", as.list(ttt))

    lines <- xi[[4]][,"ID"] %in% lines
    
    ret <- x

    ret[[1]] <- ret[[1]][lines,]
    ret[[2]] <- ret[[2]][lines,]
    ret[[4]] <- ret[[4]][lines,]
           
    return(ret)
}

