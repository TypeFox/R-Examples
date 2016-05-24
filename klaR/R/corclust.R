corclust <- function(x, cl = NULL, mincor = NULL, prnt = FALSE, method = "complete"){ 
    
    # x         : Data frame or matrix of variables to be grouped
    # cl        : optional vector containing class labels
    # mincor    : minimum required correlation between the variables of a cluster (a line is drawn)
    # prnt      : should distance matrix of absolute correlations be printed
    # method    : hierarchical clustering algorithm to be used (single, average, ward, default is complete linkage)

    if(is.matrix(x)) if(!is.numeric(x)) stop("Attributes of x must be numerical!")
    if(is.data.frame(x)) if(!all(sapply(x,is.numeric))) stop("Attributes of x must be numerical!")
    if(!is.null(cl)) if(!is.factor(cl))  stop("If classes are defined cl must be a factor!") 
    if(length(dim(x)) < 2)  stop("x must be either a matrix or a data frame of numeric attributes!")
    if(min(dim(x))==1)  stop("x must contain more than one variable and more than one observation!")
    if(!is.null(cl)){ 
        if(min(table(cl)) < 2) stop("At least one class consist of less than 2 observations!")
        if(min(table(cl)) < dim(x)[2]) warning("At least one of the classes consist of less observations than variables!")
        }

    
    if(is.null(cl)){
        kor <- 1-abs(cor(x,use="pairwise.complete.obs")) # not yet the matrix of distances!
        }
    else{
        kor <- matrix(1,ncol(x),ncol(x))
        for(k in seq(along=levels(cl))){
            cls <- levels(cl)[k]
            oldkor <- kor
            kor <- 1-abs(cor(x[which(cl==cls),],use="pairwise.complete.obs")) # not yet the matrix of distances!
            newkor <- (oldkor - kor) > 0
            kor[newkor] <- kor[newkor]  
            }    
        }
        
    distances <- NULL
    for (i in 1:(ncol(kor)-1)) distances <- c(distances,kor[(i+1):nrow(kor),i]) # only lower triangular matrix, without diagonal elements
    attr(distances,"Labels") <- colnames(x)
    attr(distances,"Size") <- ncol(x)
    attr(distances,"Metric") <- "absolute correlations"
    class(distances) <- "dissimilarity" 
        
    if(prnt) print(kor)
    
    variables <- distances
    res <- hclust(variables, method=method)
    ylb <- "1 - absolute correlation within cluster"
    if(method == "complete") ylb <- "1 - minimum absolute correlation within cluster"
    if(method == "average") ylb <- "1 - average absolute correlation within cluster"
    if(method == "single") ylb <- "1 - maximum absolute correlation within cluster"
    plot(res, ylab = ylb) # , main=paste("1 - absolute correlations between variables")
    if(!is.null(mincor)) abline(h=mincor, col=2)
    return(list(min.abs.cor=kor, clustering=res))
    }
