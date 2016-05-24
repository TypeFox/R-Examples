mgMEM <-
function(locD, truncation=NULL, transformation=NULL) {
    
    locD <- as.matrix(locD)
    
    ## Internal function to perform truncation and transformation
    .distTransform <- function (locD, transformation, truncProp) {
        locDTr <- locD
        if (is.null(truncProp) || !(truncProp)) {
            truncProp <- 1
        }
        if ((truncProp < 0) || (truncProp > 1)) {
            stop("memgene: distance truncation proportion must be between 0 and 1", call.=FALSE)
        }
        truncRange <- quantile(locD, truncProp)
        if (!is.null(transformation)) {
            transformation <- tolower(transformation)
            if (transformation=="exponential") {
                locDTr[locDTr > truncRange] <- 0
                locDTr <- exp(-locDTr/truncRange)
                locDTr <- 1-locDTr/max(locDTr)
                locDTr <- as.matrix(as.dist(locDTr))
                locDTr[locDTr == 0] <- 1
            } else if (transformation=="gaussian") {
                locDTr[locDTr > truncRange] <- 0
                locDTr <- exp(-(locDTr/truncRange)^2)
                locDTr <- 1-locDTr/max(locDTr)
                locDTr <- as.matrix(as.dist(locDTr))
                locDTr[locDTr == 0] <- 1
            }
        }
        else {
            ## CHECK **
            ## Is this correct way to truncate when there is no transformation?
            locDTr[locDTr > truncRange] <- 1
        }
        return(locDTr)
    }
    ## Truncation to a specified distance and transformation if required   
    #if (is.numeric(truncation) || !is.null(transformation)) {
        locD <- .distTransform(locD, transformation, truncation)
    #} 
  
    ## Truncate (or additionally truncate) the input distance matrix
    ## by determining the longest link on a minimum spanning tree
    ## of the points, and then setting the threshold as 4x this distance
    if (is.null(truncation) || is.numeric(truncation)) {
        ## Produce a minimum spanning tree using vegan
        spanning <- vegan::spantree(locD)
        threshh <- max(spanning$dist)
        locD_truncated <- locD
        locD_truncated[locD_truncated > threshh] <- 4 * threshh
        # note that locD_truncated has a main diagonal of zero so
        # that is congruent with the MEM framework as in Dray et al. (2006)
        # double-centred truncated matrix
        diag(locD_truncated) <- 4 * threshh
    }
    else {
        ## Truncation was FALSE or some other value
        locD_truncated <- locD
    }
    ## Double-centre truncated matrix
    row.wt = rep(1, nrow(locD_truncated))
    col.wt = rep(1, ncol(locD_truncated))
    st <- sum(col.wt)
    sr <- sum(row.wt)
    row.wt <- row.wt/sr
    col.wt <- col.wt/st
    Centered_Matrix <- -0.5*(locD_truncated*locD_truncated)
    row.mean <- apply(row.wt * Centered_Matrix, 2, sum)
    col.mean <- apply(col.wt *t(Centered_Matrix), 2, sum)
    col.mean <- col.mean - sum(row.mean * col.wt)
    Centered_Matrix <- sweep(Centered_Matrix, 2, row.mean)
    Centered_Matrix <- t(sweep(t(Centered_Matrix), 2, col.mean)) 
    
    ## Extract MEMs
    D.eig <- eigen(Centered_Matrix, symmetric = TRUE)
    ## Remove eigenvector with positive, near-zero variance due to centering
    Zero_eigenvalue <- which(D.eig$values==min(D.eig$values[D.eig$values > 0]))
    valuesMEM <- D.eig$values[-Zero_eigenvalue]
    vectorsMEM <- subset(D.eig$vectors, select=-Zero_eigenvalue)
    weight <- sqrt(abs(valuesMEM))
    
    ## Standardize
    vectorsMEM <- vectorsMEM %*% diag(weight)
    
    if (is.null(truncation)) {
        truncationStr <- "MST"
    } else if (!truncation) {
        truncationStr <- "None"
    } else if (is.numeric(truncation)) {
        truncationStr <- paste(round(truncation*100), "% + MST", sep="")
    }
    
    if (is.null(transformation)) {
        transformationStr <- "None"
    }
    else {
        transformationStr <- transformation
    }
    return(list(transformation=transformationStr, truncation=truncationStr,
                valuesMEM=valuesMEM, vectorsMEM=vectorsMEM))

}
