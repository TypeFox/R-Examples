#' A function to calculate the shared nearest neighbors (SNN)
#' 
#' This function calculate the shared nearest neighbors (SNN). 
#' SNN is reported to be more robust than k nearest neighbors. 
#' Firstly, the k nearest neighbor distances for each observation is calculated. 
#' Then, the shared nearest neighbor similarity is calculated based on the result of k nearest neighbor.
#' Note that k.nn should be greater than k.sel.
#' 
#' @import plyr
#' 
#' @param data is the data frame containing the observations (should be numeric data). Each row represents an observation and each variable is stored in one column.
#' @param k.nn specifies the value used for calculating the shared nearest neighbors.
#' @param k.sel specifies the number of shared nearest neighbors
#'
#' @return The function returns the matrix containing the indices of top k shared nearest neighbors for each observation
#'
#' @examples
#' Func.SNN(data=TestData[,1:2], k.nn=5, k.sel=3)
#'
#' @export

Func.SNN <- function(data, k.nn, k.sel) {
    #Get the knn index
    mat.ind <- FNN::get.knn(data = data, k = k.nn)$nn.index
    #Define distance function
    func.dist <- function(x1, x2) {
        length(intersect(x1, x2))
    }
    #Count the distance using the customized function
    mat.count <- as.matrix(proxy::dist(x = mat.ind, method = func.dist, diag = T, upper = T))
    #Formulate the final matrix for use
    mat.final <- plyr::aaply(.data = mat.count, .margins = 1, .fun = function(x) {order(x, decreasing=T)[1:k.sel]})
    return(mat.final)
}

#' Angle-based outlier detection (ABOD) algorithm 
#' 
#' This function performs the basic and aprroximated version of angle-based outlier detection algorithm. 
#' The ABOD method is especially useful for high-dimensional data, as angle is a more robust measure than distance in high-dimensional space.
#' The basic version calculate the angle variance based on the whole data. The results obtained are more reliable. However, the speed can be very slow.
#' The approximated version calculate the angle variance based on a subset of data and thereby, increasing the calculation speed.
#' This function is based on the work of Krigel, H.P., Schubert, M., Zimek, A., Angle-based outlier detection in high dimensional data, 2008.
#' 
#' @import foreach
#' @import plyr
#' @import ggplot2
#' 
#' @param data is the data frame containing the observations. Each row represents an observation and each variable is stored in one column.
#' @param basic is a logical value, indicating whether the basic method is used. The speed of basic version can be very slow if the data size is large.
#' @param perc defines the percentage of data to use when calculating the angle variance. It is only needed when basic=F.  
#' 
#' @return The function returns the vector containing the angle variance for each observation
#' 
#' @examples
#' library(ggplot2)
#' res.ABOD <- Func.ABOD(data=TestData[,1:2], basic=FALSE, perc=0.2)
#' data.temp <- TestData[,1:2]
#' data.temp$Ind <- NA
#' data.temp[order(res.ABOD, decreasing = FALSE)[1:10],"Ind"] <- "Outlier"
#' data.temp[is.na(data.temp$Ind),"Ind"] <- "Inlier"
#' data.temp$Ind <- factor(data.temp$Ind)
#' ggplot(data = data.temp) + geom_point(aes(x = x, y = y, color=Ind, shape=Ind))
#' 
#' @export

Func.ABOD <- function(data, basic=FALSE, perc) {
    i=j=NULL
    if(basic==T) {
        res <- foreach(i = 1:(dim(data)[1]), .combine = c) %dopar% {
            obs <- data[i,]
            com <- t(combn(x = c(1:dim(data)[1])[-i], m = 2))
            
            cos.angles <- foreach(j = 1:(dim(com)[1]), .combine = c) %do% {
                vec.1 <- data[com[j,1],] - obs
                vec.2 <- data[com[j,2],] - obs    
                round(acos(sum(vec.1 * vec.2)/(sqrt(sum(vec.1^2))*sqrt(sum(vec.2^2)+0.01)))/(sqrt(sum(vec.1^2))*sqrt(sum(vec.2^2)+0.01)), digits = 2)
            }
            return(var(x = cos.angles))
        }
        return(res) 
    } else {
        nu <- round(dim(data)[1]*perc, digits = 0)
        res <- foreach(i = 1:(dim(data)[1]), .combine = c) %dopar% {
            obs <- data[i,]
            index.used <- sample(x = c(1:dim(data)[1])[-i], size = nu, replace = F)
            com <- t(combn(x = index.used, m = 2))
            
            cos.angles <- foreach(j = 1:(dim(com)[1]), .combine = c) %do% {
                vec.1 <- data[com[j,1],] - obs
                vec.2 <- data[com[j,2],] - obs    
                round(acos(sum(vec.1 * vec.2)/(sqrt(sum(vec.1^2))*sqrt(sum(vec.2^2)+0.01)))/(sqrt(sum(vec.1^2))*sqrt(sum(vec.2^2)+0.01)), digits = 2)
            }
            return(var(x = cos.angles))
        }
        return(res)    
    }
}

#' Subspace outlier detection (SOD) algorithm 
#' 
#' This function performs suspace outlier detection algorithm 
#' The implemented method is based on the work of Krigel, H.P., Kroger, P., Schubert, E., Zimek, A., Outlier detection in axis-parallel subspaces of high dimensional data, 2009.
#' 
#' @import foreach
#' @import plyr
#' @import ggplot2
#' 
#' @param data is the data frame containing the observations. Each row represents an observation and each variable is stored in one column.
#' @param k.nn specifies the value used for calculating the shared nearest neighbors. Note that k.nn should be greater than k.sel.
#' @param k.sel specifies the number shared nearest neighbors. It can be interpreted as the number of reference set for constructing the subspace hyperplane.
#' @param alpha specifies the lower limit for selecting subspace. 0.8 is set as default as suggested in the original paper.
#' 
#' @return The function returns a vector containing the SOD outlier scores for each observation
#' 
#' @examples
#' library(ggplot2)
#' res.SOD <- Func.SOD(data = TestData[,1:2], k.nn = 10, k.sel = 5, alpha = 0.8)
#' data.temp <- TestData[,1:2]
#' data.temp$Ind <- NA
#' data.temp[order(res.SOD, decreasing = TRUE)[1:10],"Ind"] <- "Outlier"
#' data.temp[is.na(data.temp$Ind),"Ind"] <- "Inlier"
#' data.temp$Ind <- factor(data.temp$Ind)
#' ggplot(data = data.temp) + geom_point(aes(x = x, y = y, color=Ind, shape=Ind))
#' 
#' @export

Func.SOD <- function(data, k.nn, k.sel, alpha=0.8) {
    i=j=NULL
    mat.ref <- Func.SNN(data = data, k.nn = k.nn, k.sel = k.sel)
    res <- foreach::foreach(i=1:dim(data)[1], .combine = c) %dopar% {
        obs <- data[i,]
        ref <- as.matrix(data[mat.ref[i,],])
        means <- colMeans(ref)
        var.total <- sum(aaply(.data = ref, .margins = 1, .fun = function(x) sum((x-means)^2)))/k.sel
        var.expect <- alpha*var.total/dim(data)[2]
        var.actual <- foreach(j = 1:dim(ref)[2], .combine = c) %dopar% {
            var(ref[,j])
        }
        var.ind <- ifelse(var.actual<var.expect, yes = 1, no = 0)
        res.hyper <- sqrt(sum(var.ind*(obs-means)^2))/length(which(var.ind==1))
        return(res.hyper)
    }
    return(res)
}

#' Feature-bagging outlier detection (FBOD) algorithm 
#' 
#' This function performs feature-bagging based outlier detection algorithm 
#' The implemented method is based on the work of "Lazarevic, A., Kumar, V., Feature bagging for outlier detection, 2005"
#' This method can be regarded as an ensemble method, which based on the results of local outlier factor (LOF). 
#' During each iteration, a random subset of variables, whose size is randomly chosen between d/2 to d (where d is the dimensionality of the input data), is selected.
#' The LOF method is applied to calculate the LOF scores based on the selected data subset.
#' The final score of FBOD is the cumulative sum of each iteration.
#'  
#' @import foreach
#' @import plyr
#' @import ggplot2
#' 
#' @param data is the data frame containing the observations. Each row represents an observation and each variable is stored in one column.
#' @param iter is the iteration used.
#' @param k.nn is the value used for calculating the LOF score
#' 
#' @return The function returns a vector containing the FBOD outlier scores for each observation
#' 
#' @examples
#' library(ggplot2)
#' res.FBOD <- Func.FBOD(data = TestData[,1:2], iter=10, k.nn=5)
#' data.temp <- TestData[,1:2]
#' data.temp$Ind <- NA
#' data.temp[order(res.FBOD, decreasing = TRUE)[1:10],"Ind"] <- "Outlier"
#' data.temp[is.na(data.temp$Ind),"Ind"] <- "Inlier"
#' data.temp$Ind <- factor(data.temp$Ind)
#' ggplot(data = data.temp) + geom_point(aes(x = x, y = y, color=Ind, shape=Ind))
#' 
#' @export

Func.FBOD <- function(data, iter, k.nn) {
    i=NULL
    res <- foreach(i=1:iter, .combine = cbind) %dopar% {
        d <- dim(data)[2]
        l <- sample(x = round(d/2, digits = 0):(d-1), size = 1)
        ind <- sample(x = 1:d, size = l, replace = F)
        data.use <- data[,ind]
        score <- DMwR::lofactor(data = data.use, k = k.nn)
        return(score)
    }
    res[is.nan(res)] <- NA
    res[is.infinite(res)] <- max(res[!is.infinite(res)], na.rm = T)
    res.final <- plyr::aaply(.data = res, .margins = 1, .fun = function(x) mean(x, na.rm = T))
    return(round(res.final, digits = 3))
}

#' Outlier score transformation
#' 
#' This function calculate the transformed outlier scores, with the aim of unifying the results from different methods.
#' The method is based on the work of Kriegel, H.P., Kroger, P., Schubert, E., Zimek, A., Interpreting and unifying outlier scores, 2011.
#' It consists of two steps, regularization and normalization. 
#' For the ABOD scores, logarithmic inversion is used for regularization
#' For the SOD scores, no action is taken to perform regularization
#' For the FBOD method, the basic regularization, i.e., score-1, is used for regularization
#' For the normalization step, the gaussian scaling method is used.
#' The final output can be interpreted as the outlier probability, ranging from 0 to 1.
#' 
#' @param raw.score is the scores returned by each method
#' @param method should be a character specifying  the method used to generate the raw score. It has 3 possible values, "ABOD", "SOD", and "FBOD".
#' @return The function returns the transformed outlier scores
#' @export

Func.trans <- function(raw.score, method) {
    #Regularization
    if(method=="FBOD") {score.reg <- raw.score-1}
    if(method=="ABOD") {score.reg <- -log(x = raw.score/max(raw.score), base = 10)}
    if(method=="SOD") {score.reg <- raw.score}
    
    #Normalization
    erf <- function(x) 2*pnorm(x*sqrt(2))-1
    if (sd(score.reg, na.rm = T)==0) {score.norm <- rep(0, length(raw.score))} else {
        score.norm <- erf(x=(score.reg-mean(score.reg, na.rm = T))/(sqrt(x=2)*sd(score.reg, na.rm = T)))
        score.norm[which(score.norm<0)] <- 0   
    }
    return(score.norm)
}

#' Testing data for testing the performance of different algorithms
#' 
#' The data are generated according to the example in the paper "Kriegel, Kroger, Schubert, and Arthur Zimek, 2009, Outlier Detection in Axis-Parallel Subspaces of High Dimensional Data".
#' The data has 60 rows and 3 variables. The first two variables are the x and y coordinates.
#' The third variable indicates the type of data, i.e., "Pattern_1", "Pattern_2", or "Outlier".
#' The first 25 observations belong to "Pattern_1". Another 25 observations represent "Pattern_2". 
#' The other 10 observations are "Outliers". 
#'  
#' @docType data
#' 
#' @usage data(TestData)
#' 
#' @format A data frame with 60 rows and 3 variables:
#' 
#' @examples
#' library(ggplot2)
#' data(TestData)
#' ggplot(data = TestData, aes(x = x, y = y, shape=Lab, color=Lab)) + geom_point()

"TestData"

#' The player statistics of the Golden States Warriors in the season 2013-2014
#' 
#' The data contains the statistics of the players in one NBA team, i.e., Golden States Warriors, during the season 2013-2014.
#' It can be obtained from the following website: http://www.basketball-reference.com.
#' The data has 18 rows since there were 18 players shown in the lineups for the Golden States Warriors during that season.
#' The data has 27 columns, including the player names, age, games played, games started, minutes played, field goal made, field goal attemps, field goal percentage, 3-pointers made, 3-pointer attemps, 3-pointer percentage, 2-pointers made, 2-pointer attemps, 2-pointer percentages, effective field goal percentage, free throws made, free throw attemps, free throw percentage, offensive rebounds, defensive rebounds, total rebounds, assists, steals, blocks, turnovers, personal fouls, and total points.  
#'  
#' @docType data
#' 
#' @usage data(TestData)
#' 
#' @format A data frame with 18 rows and 27 variables:

"GoldenStatesWarriors"