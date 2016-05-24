rqda <- 
function(x, ...) UseMethod("rqda")


rqda.formula <- function(formula, data, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("rqda")
    output$call <- call
    if(missing(data))
    {
        cat("There is no data set.\n\n")
        return(output)
    }
    if(is.null(data))
    {
        cat("There is no data set.\n\n")
        return(output)
    }
    if(class(data) != "data.frame") 
    {
        cat("data parameter is not a data.frame.\n\n")
        return(output)
    }  
    dataset <- model.frame(formula = formula, data = data)
    grouping <- dataset[, 1]
    dataset <- dataset[, -1, drop = FALSE]
    output <- rqda.default(dataset, grouping, ...)
    output$call <- call
    output
} 


rqda.data.frame <-
function(x, grouping, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("rqda")
    output$call <- call
    if(missing(grouping))
    {
        cat("grouping vector is missing.\n\n")
        return(output)
    }   
    output <- rqda.default(x, grouping, ...)
    output$call <- call
    output
}


rqda.matrix <- 
function(x, ...)
{
    output <- rqda(as.data.frame(x), ...)
    call <- match.call()
    call[[1L]] <- as.name("rqda")   
    output$call <- call
    output
}


rqda.default <- 
function(x, grouping, subset = NULL, resmatrix = NULL, restext = NULL, 
         gamma = c(0, 1), prior = NULL, ...)
{
    output <- list()
    output$call <- match.call()

    if(missing(x))
    {
        cat("data set parameter is missing")
        return(output)
    }
    if(missing(grouping))
    {
        cat("grouping vector is missing.\n\n")
        return(output)
    }
    if(class(x) != "data.frame")
    {
        cat("data set parameter is not a data.frame.\n\n")
        return(output)
    }

    check <- checks(x, grouping, subset, resmatrix, restext, gamma, prior)
    if(is.null(check))
        return(output)
    trainset <- check$trainset
    traingroups <- check$traingroups
    numgroups <- check$numgroups
    dimension <- check$dimension
    sizes <- check$n
    resmatrix <- check$resmatrix
    resvector <- rep(0, dim(rbind(resmatrix))[1])
    prior <- check$prior
    m <- check$m
    rm(check)
   
    samplemeans <- getmeans(trainset, traingroups)
    rownames(samplemeans) <- paste("class", 1:numgroups, sep = "")

    meansVector <- c(array(t(samplemeans)))
    variances <- getvariances(trainset, traingroups)

    output$trainset <- cbind(trainset, traingroups)
    output$restrictions <- restrictions(resmatrix, dimension)
    output$resmatrix <- resmatrix
    output$prior <- prior
    output$counts <- sizes
    output$N <- sum(sizes)
    output$samplemeans <- samplemeans
    rownames(variances) <- colnames(variances) <- colnames(samplemeans)
    dimnames(variances)[[3]] <- rownames(samplemeans)
    output$samplevariances <- variances
    output$gamma <- gamma

    greatS <- array(0, c(dimension * numgroups, dimension * numgroups))
   
    for(i in 1:numgroups)
        variances[, , i]/sizes[i] -> greatS[1:dimension + dimension * {i - 1},
                                            1:dimension + dimension * {i - 1}]

    estimation <- list()
    estimatedmeans <- array(NA, c(numgroups, dimension, m))

    trainClassification <- array(0, c(sum(sizes), m + 1))
    trainClassification[, 1] <- traingroups

    for(j in 1:m)
    {
        estimation[[j]] <- rest.est(gamma[j], meansVector, resvector, resmatrix, greatS)
        estimatedmeans[, , j] <- t(array(estimation[[j]],c(dimension, numgroups)))
        trainClassification[, j + 1] <- classify(trainset, estimation[[j]], variances,
                                                 prior, numgroups, dimension) 
    }
    rownames(estimatedmeans) <- rownames(samplemeans)
    colnames(estimatedmeans) <- colnames(samplemeans)
    dimnames(estimatedmeans)[[3]] <- paste("gamma=", gamma, sep = "")
   
    output$estimatedmeans <- estimatedmeans

    apparent <- 100 * colSums(trainClassification[, -1] != traingroups) / sum(sizes)
    names(apparent) <- paste("gamma=", gamma, sep = "")

    output$apparent <- apparent
    class(output) <- "rqda"
    output
} 


print.rqda <- 
function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    if(!is.null(x$apparent))
    {
        cat("\nRestrictions:\n")
        cat(x$restrictions)
        cat("\nPrior probabilities of classes:\n")
        print(x$prior)
        cat("\nApparent error rate (%):\n")
        print(x$apparent)
    }
    cat("\n")
}


predict.rqda <-
function(object, newdata, prior=object$prior, gamma = object$gamma, grouping = NULL, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("predict")
    output$call <- call

    if(!inherits(object, "rqda"))
    {
        cat("object not of class \"rqda\".\n\n")
        return(output)
    }
    numgroups <- length(object$prior)
    dimension <- nrow(object$samplevariances)
    if(missing(newdata))
    {
        cat("newdata is missing.\n\n")
        return(output)
    }
    if(class(newdata) != "data.frame")
    {
        cat("newdata parameter is not a data.frame.\n\n")
        return(output)
    }
    if(sum(rownames(object$samplevariances) %in% colnames(newdata)) != dimension)
    {
        cat("Missing variables in newdata set.\n\n")
        return(output)
    }
    if(!is.null(grouping))
    {
        if(class(grouping) == "factor")
        {
            if(sum(suppressWarnings(is.na(as.numeric(levels(grouping))))) > 0)
            {
                cat("If grouping variable is a factor, its levels must be numbers.\n\n")
                return(NULL)
            }  
            grouping <- levels(grouping)[grouping]
        }
        grouping <- c(as.matrix(as.numeric(grouping)))
        if(dim(cbind(array(as.matrix(grouping)))) != c(dim(newdata)[1], 1) ||
            sum(grouping > floor(grouping)) > 0)
        {
            cat("Invalid grouping vector.\n\n")
            return(output)
        }
    }
    rown <- rownames(newdata)
    newdata <- newdata[, rownames(object$samplevariances), drop = FALSE]
    if(sum(!sapply(as.list(newdata), typeof) %in%
        c("integer", "double", "complex", "logical")) > 0)
    {
        cat("Invalid newdata.\n\n")
        return(NULL)
    }
    newdata <- as.matrix(newdata)
    rownames(newdata) <- rown
    if(sum(is.na(newdata)) > 0)
    {
        cat("Missing values in the newdata set have been deleted.\n\n")
        newdata <- newdata[complete.cases(newdata), , drop = FALSE]
    }
    if(sum(is.na(newdata)) > 0)
    {
        cat("There are missing values in the newdata set.\n\n")
        newdata <- newdata[complete.cases(newdata), , drop = FALSE]
        if(!is.null(grouping))
            grouping <- grouping(rownames(newdata))
    }
    if(is.null(prior))
    {
        cat("prior parameter is NULL.\n\n")
        return(output)
    }
    if(dim(rbind(prior))[1] != 1)
    {
        cat("Invalid prior parameter.\n\n")
        return(output)
    }
    if(length(prior) != numgroups)
    {
        cat("Wrong number of classes in a priori probabilities.\n\n")
        return(output)
    }
    if(sum(prior > 1 | prior < 0) > 0 || abs(sum(prior) - 1) > 1e-12)
    {
        cat("prior values must be in the interval [0,1] and sum 1.\n\n")
        return(output)
    }
    if(is.null(gamma))
    {
        cat("gamma parameter is NULL.\n\n")
        return(output)
    }
    if(sum(gamma %in% object$gamma) < length(gamma))
    {
        cat("Invalid gamma values.\n\n")
        return(output)
    }
    gammaindex <- c(outer(gamma, object$gamma, "==") %*% {1:length(object$gamma)})
    m <- length(gamma)
    classification <- array(NA, c(dim(rbind(newdata))[1], m))
    posteriorprob <- array(NA, c(dim(rbind(newdata))[1], numgroups, m))

    for(j in 1:m)
    {
        posteriorprob[, , j] <- posterior(newdata, array(t(object$estimatedmeans[, , gammaindex[j]])),
                                          object$samplevariances, prior, numgroups, dimension)
        maxClass <- {posteriorprob[, , j] == apply(posteriorprob[, , j, drop = FALSE], 1, max)} * rep(1:numgroups,
                                                   each = length(as.matrix(posteriorprob[, , j])) / numgroups)
        classification[, j] <- apply(maxClass, 1, function(x) x[x != 0][sample(length(x[x != 0]), 1)])
    }

    rownames(posteriorprob) <- rownames(newdata)
    colnames(posteriorprob) <- paste("class", 1:numgroups, sep = "")
    dimnames(posteriorprob)[[3]] <- paste("gamma=", gamma, sep = "")
    namesnewdata <- paste("gamma=", gamma, sep = "")
    namesnewdata[1] <- paste("Predicted:", namesnewdata[1])
    colnames(classification) <- namesnewdata
    rownames(classification) <- rownames(newdata)
    output$class <- classification
    output$prior <- prior
    output$posterior <- posteriorprob
    if(!is.null(grouping))
    {
        classif <- classification
        if(sum(is.na(grouping)) > 0)
        {
            classif <- classif[!is.na(c(grouping)), ]
            grouping <- c(grouping)[!is.na(c(grouping))]
            cat("Missing values in grouping set have been deleted and will not be considered on calculating the true error rate.\n\n")
        }
        error.rate <- array(0, c(1, m))
        error.rate[1,] <- 100 * colSums(classif != grouping) / dim(classif)[1]
        colnames(error.rate) <- paste("gamma=", gamma, sep = "")
        rownames(error.rate) <- c("True error rate (%):")
        output$error.rate <- error.rate
        classification <- cbind(grouping, classification)
        colnames(classification)[1] <- "Observed"
        output$class <- classification
    }
    cat("Call:\n")
    print(output$call)
    output
}


err.est.rqda <-
function(x, nboot = 50, gamma = x$gamma, prior = x$prior, ...)
{
    time <- proc.time()
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("err.est")
    output$call <- call

    if(!inherits(x, "rqda"))
    {
        cat("object not of class \"rqda\".\n\n")
        return(output)
    }
    numgroups <- length(x$prior)
    if(is.null(nboot) || length(array(as.array(nboot))) > 1 || 
        nboot > floor(nboot))
    {
        cat("Invalid nboot parameter.\n\n")
        return(output)
    }
    if(is.null(prior))
    {
        cat("prior parameter is NULL.\n\n")
        return(output)         
    }
    if(dim(rbind(prior))[1] != 1)
    {
        cat("Invalid prior parameter.\n\n")
        return(output)
    }
    if(length(prior) != numgroups)
    {
        cat("Wrong number of classes in a priori probabilities.\n\n")
        return(output)
    }
    if(sum(prior > 1 | prior < 0) > 0 || abs(sum(prior) - 1) > 1e-12)
    {
        cat("prior values must be in the interval [0,1] and sum 1.\n\n")
        return(output)
    }
    if(is.null(gamma))
    {
        cat("gamma parameter is NULL.\n\n")
        return(output)
    }
    if(sum(gamma %in% x$gamma) < length(gamma))
    {
        cat("Invalid gamma values.\n\n")
        return(output)
    }
    names(prior) <- names(x$prior)
    resmatrix <- x$resmatrix
    resvector <- rep(0, dim(rbind(resmatrix))[1])
    dimension <- nrow(x$samplevariances)
    m <- length(gamma)
    n <- x$counts
    traingroups <- x$trainset[, dim(x$trainset)[2]]
    trainset <- x$trainset[, -dim(x$trainset)[2], drop=FALSE]
    rm(x)
    if(sum(n < dimension / 0.632) > 0)
    {
        cat("Not enough training sample size for bootstraping for at least one class.\n\n")
        return(output)
    }

    datagroups <- restrictedMeans <- bootHowMany <- bootWho <- list()
    casesClassified <- casesClassifiedcv <- 0
    indexes <- c()
    for(i in 1:numgroups)
        for(j in 1:dimension)
            for(k in 1:dimension)
               indexes <- c(indexes, {i - 1}*dimension + j, {i - 1}*dimension + k)
    indexes <- matrix(indexes, ncol = 2, byrow = TRUE)

    variancesb <- array(NA, c(dimension, dimension, numgroups))
    greatS <- array(0, c(dimension * numgroups, dimension * numgroups))
    misclassifs2 <- misclassifs3 <- misclassifs2cv <- misclassifs3cv <- rep.int(0, m)
    first <- TRUE
    mean <- rep(NA, dimension)
    cases <- vector("list", numgroups)

    trainingMeans <- array(getmeans(trainset, traingroups), c(numgroups, dimension))
    trainingMeansVector <- c(array(t(trainingMeans)))

    A2 <- resmatrix
    signs <- resmatrix %*% cbind(trainingMeansVector) >= 0
    A2[signs,] <- (-1) * A2[signs, ]

    variances <- getvariances(trainset, traingroups)
    greatS[indexes] <- array(variances) / rep(n, each = dimension*dimension)

    for(group in 1:numgroups)
    {
        datagroups[[group]] <- trainset[traingroups == group, , drop = FALSE]
        mb <- bootstrapsamples(datagroups[[group]], 5 * nboot)
        valid <- rowSums(mb$howMany > 0) > dimension
        bootHowMany[[group]] <- mb$howMany[valid, , drop = FALSE]
        bootWho[[group]] <- mb$who[valid, , drop = FALSE]
    }

    rm(mb, valid, signs)
   
    for(j in 1:m) 
        restrictedMeans[[j]] <- rest.est(gamma[j], trainingMeansVector,
                                         resvector, resmatrix, greatS)
    
    time2 <- proc.time()
    meansVector <- rep(NA, dimension * numgroups)
    nClasses <- rep(NA, numgroups)
    for(k in 1:nboot)
    {
        for(group in 1:numgroups)
        {
            dd <- datagroups[[group]][bootWho[[group]][k, ], , drop = FALSE]
            meansVector[1:dimension + (group - 1) * dimension] <- colMeans(dd)
            variancesb[, , group] <- var(dd)
            cases[[group]] <- datagroups[[group]][bootHowMany[[group]][k, ] == 0, , drop = FALSE]
            nClasses[group] <- dim(cases[[group]])[1]
        }
        toClassSum <- sum(nClasses)
        meansVector3 <- meansVector - trainingMeansVector
        if(toClassSum>0)
        {
            casesClas <- do.call(rbind, cases)
            greatS[indexes] <- array(variancesb) / rep(n, each = dimension * dimension)
            casesClassified <- casesClassified + toClassSum
            casesClasses <- rep.int(1:numgroups, times = nClasses)
            for(j in 1:m)
            {
                misclassifs2[j] <- misclassifs2[j] + 
                                   sum(classify(casesClas, rest.est(gamma[j], meansVector, resvector, A2, greatS), 
                                   variancesb, prior, numgroups, dimension) != casesClasses)
                misclassifs3[j] <- misclassifs3[j] + sum(classify(casesClas - trainingMeans[casesClasses,] +
                                   t(array(restrictedMeans[[j]], c(dimension, numgroups)))[casesClasses, , drop = FALSE],
                                   rest.est(gamma[j], meansVector3 + restrictedMeans[[j]], resvector, resmatrix,
                                   greatS), variancesb, prior, numgroups, dimension) != casesClasses)
            }
        }
        for(group in 1:numgroups)
        {
            meansVectorcv <- meansVector
            meansVectorcv3 <- meansVector3
            howManySample <- bootHowMany[[group]][k, ]
            varb <- variancesb[, , group]
            for(kk in {1:n[group]}[howManySample > 0])
            {
                toClassify <- howManySample[kk]
                if(toClassify > 1 || sum(howManySample > 0) > dimension + 1)
                {
                    aux <- bootWho[[group]][k, ]
                    aux[match(kk,aux)] <- 0
                    datacv <- datagroups[[group]][aux, , drop = FALSE]
                    mean <- colMeans(datacv)
                    meansVectorcv[1:dimension + dimension * {group - 1}] <- mean
                    meansVectorcv3[1:dimension + dimension * {group - 1}] <- mean - trainingMeans[group,]
                    variancesb[, , group] <- var(datacv)
                    variancesb[, , group]/{n[group]-1} -> greatS[1:dimension + (group - 1) * dimension,
                                                                 1:dimension + (group - 1) * dimension]
                 
                    casesClassifiedcv <- casesClassifiedcv + toClassify
                    case <- datagroups[[group]][kk, , drop = FALSE]
                    for(j in 1:m)
                    {
                        misclassifs2cv[j] <- misclassifs2cv[j] + {classify(case, rest.est(gamma[j],
                                             meansVectorcv, resvector, A2, greatS), variancesb, prior,
                                             numgroups, dimension) != group}*toClassify
                        misclassifs3cv[j] <- misclassifs3cv[j] + {classify(case - trainingMeans[group, , drop = FALSE] +
                                             matrix(restrictedMeans[[j]][{1:dimension} + {group - 1}*dimension], nrow = 1),
                                             rest.est(gamma[j], meansVectorcv3 + restrictedMeans[[j]], resvector,
                                             resmatrix, greatS), variancesb, prior, numgroups, dimension) != group}*toClassify
                    }
                }  
            }
            variancesb[, , group] <- varb
            varb / n[group] -> greatS[1:dimension + {group - 1} * dimension,
                                      1:dimension + {group - 1} * dimension]
        }
        if(first)
        {
            cat(paste("The procedure will end in approximately ",
                      round({proc.time() - time2}[3] * {nboot - 1}, 2),
                      " seconds.\n", sep = ""))
            first <- FALSE
        }
    }
    result <- 100 * rbind(misclassifs2, misclassifs3, misclassifs2cv, misclassifs3cv) /
              c(casesClassified, casesClassified, casesClassifiedcv, casesClassifiedcv)
    colnames(result) <- paste("gamma=", gamma, sep = "")
    rownames(result) <- c("BT2", "BT3", "BT2CV", "BT3CV")
    output$restrictions <- restrictions(resmatrix, dimension)
    output$prior <- prior
    output$counts <- n
    output$N <- sum(n)
    output$estimationError <- result
    cat(paste("The procedure lasted ", round({proc.time() - time}[3], 2),
              " seconds.\n\n", sep = ""))
    cat("Call:\n")
    print(output$call)
    cat("\nRestrictions:\n")
    cat(output$restrictions)
    cat("\nPrior probabilities of classes:\n")
    print(output$prior)
    cat("\nTrue error rate estimation (%):\n")
    print(output$estimationError)
    invisible(output)
}
