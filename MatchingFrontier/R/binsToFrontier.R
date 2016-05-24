binsToFrontier <-
function(strataholder){
    Ys <- c()
    drop.order <- list()
    num.treated <- sum(unlist(lapply(strataholder, function(x) sum(names(x) == 1))))
    num.control <- sum(unlist(lapply(strataholder, function(x) sum(names(x) == 0))))

    Ys <- c(Ys, .5 * sum(unlist(lapply(strataholder, function(x)
                                       abs(sum(names(x) == 0) / num.control - sum(names(x) == 1) / num.treated)))))

    while(1){
        diffs <- unlist(lapply(strataholder, function(x)
                               sum(names(x) == 0) / num.control - sum(names(x) == 1) / num.treated))
        drop.from <- which(diffs == max(diffs))[1]
        dropped.element.ind <- which(names(strataholder[[drop.from]]) == 0)[1]
        
        drop <- strataholder[[drop.from]][dropped.element.ind]
        strataholder[[drop.from]] <- strataholder[[drop.from]][-dropped.element.ind]
        new.L1 <- .5 * sum(unlist(lapply(strataholder, function(x)
                                         abs(sum(names(x) == 0) / num.control - sum(names(x) == 1) / num.treated)
                                         )
                                  )
                           )
        
        if(new.L1 > tail(Ys, 1)){
            break
        }
        Ys <- c(Ys, new.L1)
        drop.order[[length(drop.order) + 1]] <- drop
        num.control <- num.control - 1        
    }
    drop.order[[length(drop.order) + 1]] <- unlist(strataholder)
    
    Xs <- 1:length(Ys)
    return(list(drop.order = drop.order, Xs = Xs, Ys = Ys))
}
