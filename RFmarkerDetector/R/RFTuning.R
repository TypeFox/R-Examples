
#' Tuning of the mtry parameter for a Random Forest model
#' 
#' \code{tuneMTRY} tries to identify the 'optimal' value of the mtry parameter which indicates the 
#' number of input variables randomly chosen at each node
#' @param data the n x p dataframe used to build the models and to tune the parameter mtry. The first two columns
#' must represent respectively the sample names and the class labels related to each sample 
#' @param iterations the number of different random forest models built for each value of mtry
#' @param maxntree the maximum number of trees of each random forest model
#' @param mtry_length an integer value representing the number of mtry values to test.
#' @param changeTreeNum a logical value indicating whether or not to change
#' @param graph a logical value indicating whether to plot the OOB error as a function of the parameter mtry 
#' the number of trees during the tuning of mtry
#' @return a list of two elements:\itemize{
#'   \item a diagram of the average value of the OOB error as a function of the mtry with its 95\% confidence interval
#'   \item a n x p matrix of n iterations and p mtry values tested 
#' }
#' @author Piergiorgio Palla
#' @import randomForest randomForest
#' @details The function searches for the optimal value of mtry assigning to it a set of values and building different random forests (also with a different 
#' number of trees) for each value of the mtry. The number of models built for each mtry is defined by the \code{iteration} parameter. 
#' The oob errors of each random forest model, computed for each mtry value, are then arranged in a matrix
#' @examples
#' ## data(cachexiaData)
#' ## res <- tuneMTRY(cachexiaData, iterations = 10, maxntree = 600, mtry_length = 10, graph = FALSE)
#' @export

tuneMTRY <- function(data, iterations, maxntree, mtry_length, changeTreeNum = F, graph = T) {
    
    # numOfVariables = dim(data)[2] - 2
    input <- data[, 3:dim(data)[2]]
    numOfVariables <- dim(input)[2]
    DEFAULT_N_TREES <- 500
    DEFAULT_ITERATIONS <- 40
    DEFAULT_MTRY_VALUES_TO_TEST <- numOfVariables
    
    
    if (hasArg("iterations") && !is.null(iterations)) {
        iterat <- iterations
    } else {
        iterat <- DEFAULT_ITERATIONS
    }
    
    if (hasArg("maxntree") && !is.null(maxntree)) {
        nt <- maxntree
    } else {
        nt <- DEFAULT_N_TREES
    }
    
    if (hasArg("mtry_length") && !is.null(mtry_length) && mtry_length <= numOfVariables) {
        mtry_length <- mtry_length
    } else {
        mtry_length <- DEFAULT_MTRY_VALUES_TO_TEST
    }
    
    set.seed(1234)
    mtry_values <- seq(from = 1, to = numOfVariables, length = mtry_length)
    typical <- sqrt(numOfVariables)
    
    ### Here we add the optimal mtry
    if (length(mtry_values) < numOfVariables) {
        mtry_values <- c(mtry_values, typical)
    }
    mtry_values <- floor(sort(mtry_values))
    mtry_values <- unique(mtry_values)
    
    
    if (changeTreeNum == T) {
        ntree_values <- seq(from = 50, to = nt, length = iterat)
        ntree_values <- floor(ntree_values)
    } else {
        ntree_values <- rep(x = nt, times = iterat)
    }
    
    
    global_oob = c()
    
    START <- proc.time()
    
    # set.seed(115)
    for (mtr in mtry_values) {
        ### The sequence of models built for each mtry value starts from the same seed ####
        set.seed(1234)
        oob = c()
        for (i in 1:iterat) {
            ntree <- ntree_values[i]
            class <- data[, 2]
            rf <- randomForest::randomForest(class ~ ., data = input, mtry = mtr, ntree = ntree)
            print(rf)
            ####### We consider only the error rate after the completion of the model
            errors <- rf$err.rate[ntree, 1]
            oob <- append(oob, errors)
        }
        global_oob = cbind(global_oob, oob)
    }
    
    colnames(global_oob) <- mtry_values
    
    ELAPSED <- proc.time() - START
    cat(paste("Elapsed time:", round(ELAPSED[3], 2), "sec", sep = " "), "\n")
    
    
    cat("\nNumber of trees used to build each Random Forest model\n")
    cat(ntree_values, "\n")
    
    cat("\nMTRY values tested\n")
    cat(mtry_values)
    
    ### 
    g <- global_oob
    
    ## Here we set up the data needed to plot the 95% confidence interval of the mean of the OOB errors for each mtry value
    
    res <- optimizeMTRY(g)
    
    ### Then we plot with the package ggplot2
    mean_matrix <- res$mean_matrix
    ci_matrix <- res$ci_matrix
    
    if (graph == T) {
        diagram = plotOOBvsMTRY(mean_matrix, ci_matrix)
    } else {
        diagram = NULL
    }
    
    diagram
    cat("\n\nOOB matrix\n")
    return(list(oob = g, diagram = diagram))
    
    
}

#' Mtry Optimization
#' 
#' This function provides a 'population' estimate of the average OOB error computed for different mtry values,
#' starting from a sample of N models. These values will be used to compute the mtry associated to the minimum averaged OOB error,
#' that is the optimal parameter we are looking for.
#' 
#' 
#' @param oob_matrix a n x p of n OOB error values (one for each iteration) and p columns (one for each mtry value tested)
#' Each value of a column is the oob error of a model growth with a particular mtry. Typically for each mtry,
#' we will have N different models (N > 30), a sample large enough to provide an estimate of the average OOB
#' error for the corresponding population of models.
#' @return a list of two elements: \itemize{
#' \item mean_matrix a 1 x p matrix which contains the mean of each OOB errors sample (resulting from the training of N different Random Forest models growth
#' with N different mtry values) 
#' \item ci_matrix a 2 x p matrix in which each column represents the 95\% confidence interval of the mean of the population of the OOB errors for each 
#' mtry value 
#' \item sd_matrix a 1 x p matrix which contains the standard deviatiaon of each OOB error sample resulting from the training of N different models
#' built for each value of mtry}
#' 
#' @examples
#' ## data(cachexiaData)
#' ## res <- tuneMTRY(cachexiaData, iterations = 50, maxntree = 600, mtry_length = 10, graph = F)
#' ## l <- optimizeMTRY(res$oob)
#' @author Piergiorgio Palla
#' @import UsingR
#' @export

optimizeMTRY <- function(oob_matrix) {
    
    ## Here we initialize the vector that will contain the mean of each sample and the vector that will contain the interval
    ## estimate of the 'population' (from which each sample has been drawn) mean.
    
    mean_matrix = colMeans(oob_matrix)
    sd_matrix <- apply(oob_matrix, 2, sd)
    ci_matrix = c()
    
    for (i in 1:ncol(oob_matrix)) {
        sample <- oob_matrix[, i]
        
        ### sample standard deviation (denominator: n-1)
        stdev <- sd(sample)
        ### The use of the z test requires the knowledge of the variance / standard deviation of the population.  In this case we do
        ### not know those values, but for 'large' sample sizes (n > 30) the sample standard deviation is a good unbiased estimator
        ### of the corresponding parameter. This way we can use the z-test to estimate our confidence intervals of the population
        ### mean oob errors. We can also use a one sample t-test.
        
        ci <- UsingR::simple.z.test(sample, sigma = stdev, conf.level = 0.95)
        ci_matrix <- cbind(ci_matrix, ci)
        
        # ci <- t.test(sample) ci_matrix <- cbind(ci_matrix, ci)
        
    }
    
    colnames(ci_matrix) <- names(mean_matrix)
    res <- list(mean_matrix = mean_matrix, ci_matrix = ci_matrix, sd_matrix = sd_matrix)
    
}


#' Tuning of the ntree parameter (i.e. the number of trees) for a Random Forest model
#' 
#' This function tries to find the 'optimal' value for the parameter ntree which indicates 
#' the number of trees used to grow the ensemble of trees. To do that it will build several random forest models with a different 
#' number of trees for the mtry value considered. The number of models built for each ntree value will be equal to
#' the parameter iteration. The oob errors of each random forest model, computed 
#' for each ntree value will be arranged in a matrix.
#' @param data the n x p dataframe used to build the Random Forest models. The first two columns must represent respectively
#' the sample names and the class labels associated to each sample
#' @param mtry the chosen mtry value
#' @param iterations the number of Random Forest models to be built for each value of ntree
#' @param minNTREE the minimum number of trees of each random forest model.
#' @param pace the pace between each value of ntree to be tested
#' @param seq_length the number of ntree values to be tested
#' @return a n x p matrix in which n is the number of models considered and p is the number of ntree values tested.
#' Each column represents the oob errors resulting from each model and corresponding to the different ntree values
#' @examples
#' ## data(cachexiaData)
#' ## res <- tuneNTREE(cachexiaData, 8, iterations = 50, minNTREE = 600, pace = 100, seq_length = 10)
#' @author Piergiorgio Palla
#' @export

tuneNTREE <- function(data, mtry, iterations, minNTREE = 500, pace = 100, seq_length = 5) {
    
    ### input dataset ###
    input <- data[, 3:dim(data)[2]]
    
    DEFAULT_NTREE <- 500
    
    if (!hasArg(minNTREE)) {
        minNTREE = DEFAULT_NTREE
    }
    
    
    ntree_range <- seq(from = minNTREE, by = pace, length.out = seq_length)
    
    
    
    START <- proc.time()
    
    ### Output matrix initialization ###
    global_oob = c()
    
    for (t in ntree_range) {
        ### The sequence of models built for each mtry value starts from the same seed ####
        set.seed(1234)
        oob = c()
        for (i in 1:iterations) {
            class = data[, 2]
            rf <- randomForest::randomForest(class ~ ., data = input, mtry = mtry, ntree = t)
            print(rf)
            ####### We consider only the error rate after the completion of the model
            errors <- rf$err.rate[t, 1]
            oob <- append(oob, errors)
        }
        global_oob = cbind(global_oob, oob)
    }
    
    colnames(global_oob) <- ntree_range
    
    ELAPSED <- proc.time() - START
    cat(paste("Elapsed time:", round(ELAPSED[3], 2), "sec", sep = " "), "\n")
    
    
    cat("\n\nNumber of trees used to build each Random Forest model\n")
    cat(ntree_range, "\n")
    
    ### Matrix with the distribution of the oob errors for each ntree value considered ### Its number of rows is equal to the
    ### number of models built for each ntree value (i.e. iterations) ####
    g <- global_oob
    
}



#' Plotting the average OOB error and its 95\% confidence interval as a function of the mtry parameter 

#' @param mean_matrix a 1 x p matrix where p is the number of mtry values tested. Each value represents the average
#' OOB error obtained training multiple Random Forest models with a defined value of mtry
#' @param ci_matrix a 2 x p matrix containing the extremes of the confidence interval of the average OOB error. 
#' @return a graphical representation of the average OOB error as a function of the mtry parameter.
#' @author Piergiorgio Palla
#' @examples
#' ## data(cachexiaData)
#' ## res <- tuneMTRY(cachexiaData, iterations = 5, maxntree = 600, mtry_length = 10, graph = F)
#' ## l <- optimizeMTRY(res$oob)
#' ## plotOOBvsMTRY(l$mean_matrix, l$ci_matrix)
# @import ggplot2

plotOOBvsMTRY <- function(mean_matrix, ci_matrix) {
    
    # print(ci_matrix[2, ])
    min_mtry <- as.numeric(names(which.min(mean_matrix)))
    ### R CMD check workaround for ggplot2 ###
    mtry <- NULL
    mtry_vector <- as.numeric(names(mean_matrix))
    OOB.error <- NULL
    L <- NULL
    U <- NULL
    
    df <- data.frame(mtry = mtry_vector, OOB.error = mean_matrix, U = ci_matrix[1, ], L = ci_matrix[2, ], min_mtry = min_mtry)
    diagram <- ggplot2::ggplot(df, ggplot2::aes(mtry)) + ggplot2::geom_line(ggplot2::aes(y = OOB.error), colour = "blue") + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin = L, ymax = U), alpha = 0.2)
    
    
    # diagram <- diagram + #scale_x_discrete(limits=names(mean_matrix))
    
    diagram <- diagram + ggplot2::geom_vline(ggplot2::aes(xintercept = min_mtry), colour = "#BB0000", linetype = "dashed") + 
        ggplot2::ggtitle("Mtry Optimization")  #+ annotate('text', x=13, y=0, label=13)
    diagram
    # diagram <- diagram + theme(plot.title = element_text(lineheight=.8, face='bold')) plot(df$x, df$M, col='black')
    # polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = 'grey75', border = FALSE) lines(df$x, df$F, lwd = 2) #add red lines on
    # borders of polygon lines(df$x, df$U, col='red',lty=2) lines(df$x, df$L, col='red',lty=2)
    
} 
