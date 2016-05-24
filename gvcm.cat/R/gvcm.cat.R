gvcm.cat <-
function(
formula,
data,
family = gaussian,
method = c("lqa", "AIC", "BIC"),
tuning = list(lambda=TRUE, specific=FALSE, phi=0.5, grouped.fused=0.5, elastic=0.5, vs=0.5, spl=0.5),
weights,
offset, 
start,
control,
model = FALSE,
x = FALSE,
y = FALSE,
plot=FALSE,
...
)

{

# check
    Call <- match.call()
    indx <- match(c("formula", "data"),
        names(Call), nomatch = 0)
    if (indx[1] == 0)
        stop("A formula argument is required. \n")
    if (indx[2] == 0)
        stop("A data argument is required. \n")
    if (missing(control))
        control <- cat_control(...)
    
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (family$family=="Gamma") family <- Gamma(link="log")

    if (!is.logical(model) || !is.logical(x) || !is.logical(y) || !is.logical(plot))
         stop ("Error in input arguments. \n")

    method <- match.arg(method)
    if (!(method %in% c("lqa", "AIC", "BIC"))) # check method!!
         stop ("method is incorrect. \n")

# na remove
    if (missing(data))
        data <- environment(formula)
    data <- na.omit(data)
 
# model.matrix
    dsgn <- design(formula,data)
    X <- dsgn$X
    n <- nrow(X)
    
# response
    Y <- model.extract(dsgn$m, "response")
    if (is.factor(Y)==TRUE){Y <- as.numeric(Y)-1}
    if (missing(weights))
        weights <- rep(1, times=n)
    if (length(weights)!=nrow(X) || !is.vector(weights) || !is.numeric(weights))
        stop("Error in input weights. ") 
    if (!is.null(dim(Y)[2]) && family$family=="binomial") {
        weights <- (Y[,1]+Y[,2])*weights 
        Y <- Y[,1]/(Y[,1]+Y[,2])
        } 

    if (family$family=="binomial" && (sum(Y>1) || sum(Y<0))) 
        stop("No binomial response. \n") 
    if (family$family=="Gamma" && (sum(Y<=0))) 
        stop("No Gamma-distributed response. \n") 
        
# definitions
    indices <- index(dsgn, data, formula)
    indices["index2b",1] <- if (control$assured.intercept) 0 else 1 

# standardize - splines ausgenommen!!
    not <- c()
    if (sum(X[,1])==n || sum(X[,1])==sum(weights)) {not <- 1}
    if (any(rowSums(indices)[c(7,10)]!=0)){   
        sm <- which( c(indices[7,]+ indices[10,]) != 0)
        for (i in 1:length(sm)) { not <- c(not, (sum(indices[1, 1:sm[i]])-indices[1,sm[i]]+1):(sum(indices[1, 1:sm[i]])) )}
    }
    not <- if (length(not)>0) {-not} else {1:ncol(X)}
    
    if(control$center){
       centering <- colSums(diag(weights)%*%X)/sum(weights) # colMeans(X[,-1])
       X[, not] <- scale(X[, not], center = centering[not], scale = FALSE)           
    }
    if(control$standardize){
       scaling <- sqrt(colSums(t((t(X) - colSums(diag(weights)%*%X)/sum(weights))^2*weights))/(sum(weights)-1))
       X[, not] <- scale(X[, not], center = FALSE, scale = scaling[not])
    } 


# default method        
    if (method %in% c("AIC", "BIC")) {
        output <-  abc(X, Y, indices, family, method, weights, offset, start, control, plot)
        } else {
        output <- pest(X, Y, indices, family, tuning, weights, offset, start, control, plot)
        }
    if (!exists("output"))
        stop("Error in argument 'method'")     

# boostrap
    if (control$bootstrap>0) {
    bootstrap.errors <- bootstrap(X, Y, indices, family, tuning=output$tuning, weights, 
                               offset, start, control=output$control, method)
    } else {
    bootstrap.errors <- NULL
    }

# output                                
    output$call <- Call
    output$formula <- dsgn$formula
    output$terms <- dsgn$Terms
    output$data <- data 
    output$x <- if(x==TRUE) X else NULL
    output$y <- if(y==TRUE) Y else NULL
    output$model <- if(model==TRUE) dsgn$m else NULL
    output$xlevels <- .getXlevels(dsgn$Terms, dsgn$m)
    output$bootstrap.errors <- bootstrap.errors
    output$method <- method
    output$scaling <- if(control$standardize == TRUE) scaling else NULL
    output$centering <- if(control$center == TRUE) centering else NULL
    class(output) <- c("gvcm.cat", "glm", "lm")
    output

}

