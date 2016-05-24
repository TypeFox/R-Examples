

getblanketboosting <- function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
    p <- ncol(X)
    
    if(p <=maxNoVariables)
    {
        usevar <- 1:p
    }else
    {
        usevar <- selLmBoost(X = cbind(X,Y), k = p + 1, output = FALSE, pars = list(atLeastThatMuchSelected = 0.00, atMostThatManyNeighbors = maxNoVariables))
        usevar <- which(usevar)
        
    }
    
    testsets <- list()
    if(length(usevar)>0)
    {
        for (ic in ((1:2^length(usevar))-1))
        {
            testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
        }
    }
    testsets <- unique(testsets)
    le <- sapply(testsets,length)
    testsets <- testsets[ keep <- which(le>0 & le <= maxNoVariablesSimult) ]
    testsets <- testsets[order(le[keep])]
    return(testsets)
}

selLmBoost <- function(X,k,pars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),output = FALSE)
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


train_LMboost <- function(X,y,pars = list()) #
{
   
    y <- y - rep( mean(y), length(y))
    
    ## begin old version
    #dat <- data.frame(cbind(yy,X))
    #gb <- glmboost(yy ~ .,data=dat)
    # EXPLANATION: surprisingly, it turned out that this cannot be applied to large p (private discussion with T. Hothorn in Sep 2013)
    yy <- as.vector(y)
    options(warn=-1)
    gb <- glmboost(X,yy, center = TRUE)
    options(warn=1)
    ## end old version
    
    ## begin new version
    #dat <- as.data.frame(X)
    #bl <- lapply(dat, bols)
    #gb <- mboost_fit(bl, y)
    ## end new version
    
    result <- list()
    result$Yfit <- gb$fitted()
    result$residuals <- gb$resid()
    result$model <- gb
    return(result)
}
