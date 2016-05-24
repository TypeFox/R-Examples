checks <-
function(trainset, traingroups, subset, resmatrix, restext, gamma, prior)
{
    if(!is.null(subset))
    {
        if(dim(cbind(subset))[2] != 1 || dim(cbind(subset))[1] > dim(trainset)[1])
        {
            cat("Invalid index vector.\n\n")
            return(output)
        }
        dataset <- cbind(traingroups, trainset)[subset, ]
        traingroups <- dataset[, 1]
        trainset <- dataset[, -1, drop = FALSE]
    }

    if(class(traingroups) == "factor")
    {
        if(sum(suppressWarnings(is.na(as.numeric(levels(traingroups))))) > 0)
        {
            cat("If grouping parameter is a factor, its levels must be numbers.\n\n")
            return(NULL)
        }
        traingroups <- as.numeric(levels(traingroups)[traingroups])
    }
    if(sum(!sapply(as.list(trainset), typeof) %in% 
       c("integer", "double", "complex", "logical")) > 0)
    {
        cat("Invalid data set.\n\n")
        return(NULL)
    }
    trainset <- cbind(as.matrix(trainset), traingroups)
    if(sum(is.na(trainset)) > 0)
    {
        trainset <- trainset[complete.cases(trainset), , drop = FALSE]
        cat("Missing values in the training or grouping set have been deleted.\n\n")
    }
    traingroups <- trainset[, dim(trainset)[2]]
    trainset <- trainset[, -dim(trainset)[2], drop = FALSE]

    numgroups <- length(levels(as.factor(traingroups)))
    dimension <- dim(trainset)[2]

    n <- tapply(trainset[, 1], traingroups, function(x) length(x))
    if(sum(abs(as.numeric(levels(as.factor(traingroups))) - {1:numgroups})) > 0)
    {
        cat("Classes must be identified by natural numbers varying from 1 to the number of classes.\n\n")
        return(NULL)
    }
    if(numgroups == 1)
    {
        cat("Training set has only one class.\n\n")
        return(NULL)
    }
    if(sum(n < dimension) > 0)
    {
        cat("Not enough training sample size for at least one class.\n\n")
        return(NULL)
    }
    if(!is.null(prior))
    {
        if(length(prior) != numgroups)
        {
            cat("Wrong number of classes in a priori probabilities.\n\n")
            return(NULL)
        }
        if(sum(prior > 1 | prior < 0) > 0 || abs(sum(prior) - 1) > 1e-12)
        {
            cat("prior values must be in the interval [0, 1] and sum 1.\n\n")
            return(NULL)
        }
    }
    if(is.null(gamma) || sum(gamma > 1 | gamma < 0) > 0)
    {
        cat("gamma must take values in the interval [0, 1].\n\n")
        return(NULL)
    }
    if(is.null(resmatrix) && is.null(restext))
    {
        cat("Either resmatrix or restext must be specified.\n\n")
        return(NULL)
    }
    if(!is.null(resmatrix) && dim(rbind(resmatrix))[2] != dimension * numgroups)
    {
        cat("restrictions matrix has wrong number of columns.\n\n")
        return(NULL)
    }
    if(is.null(resmatrix) && !is.null(restext))
    {
        resmatrix <- resmatrix(restext, numgroups, dimension)
        if(is.null(resmatrix))
            return(NULL)
    }

    m <- length(gamma)

    if(is.null(prior))
        prior <- n / sum(n)
    names(n) <- names(prior) <- paste("class", 1:numgroups, sep = "")

    output <- list()
    output$trainset <- trainset
    output$traingroups <- traingroups
    output$numgroups <- numgroups
    output$dimension <- dimension
    output$n <- n
    output$resmatrix <- resmatrix
    output$prior <- prior
    output$m <- m
    output
}
