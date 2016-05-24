selLmBoost <-
function(X,pars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),output = FALSE,k)
{
    if(output)
    {
        cat("Performing variable selection for variable", k, ": \n")
    }
    result <- list()
    p <- dim(as.matrix(X))
    if(p[2] > 1)
    {
        selVec <- rep(FALSE, p[2])
        a <- X[,-k]
        #names(a) <- LETTERS[1:(dim(X)[2]-1)]
        b <- X[,k]
        #names(b) <- "b"
        modfitLm <- train_LMboost(X[,-k],X[,k],pars)
        cc <- unique(modfitLm$model$xselect())
        if(output)
        {
            cat("The following variables \n")
            show(cc)
        }
        nstep <- length(modfitLm$model$xselect())
        howOftenSelected <- rep(NA,length(cc))
        for(i in 1:length(cc))
        {
            howOftenSelected[i] <- sum(modfitLm$model$xselect() == cc[i])/nstep
        }
        if(output)
        {
            cat("... have been selected that many times: \n")
            show(howOftenSelected)
        }
        howOftenSelectedSorted <- sort(howOftenSelected, decreasing = TRUE)
        if( sum(howOftenSelected>pars$atLeastThatMuchSelected) > pars$atMostThatManyNeighbors)
        {
            cc <- cc[howOftenSelected>howOftenSelectedSorted[pars$atMostThatManyNeighbors + 1]]
        } else
        {
            cc <- cc[howOftenSelected>pars$atLeastThatMuchSelected]
        }
        if(output)
        {
            cat("We finally choose as possible parents: \n")
            show(cc)
            cat("\n")
        }
        tmp <- rep(FALSE,p[2]-1)
        tmp[cc] <- TRUE
        selVec[-k] <- tmp
    } else
    {
        selVec <- list()
    }
    return(selVec)
}
