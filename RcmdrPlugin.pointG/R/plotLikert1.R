plotLikert1<-
function (X, y, tri = 0, scale=0, ...) 
{
    trellis.device(color = FALSE)
    p <- ncol(X)
    X <- as.matrix(X)
    Min <- min(X, na.rm = TRUE)
    Max <- max(X, na.rm = TRUE)
    facto <- factor(y)
    k <- length(levels(facto))
    XXX <- matrix(0, nrow = k, ncol = p)
    for (j in 1:p) {
        XXX[, j] <- tapply(X[, j], facto, mean, na.rm = TRUE)
    }
    XXXX <- as.vector(XXX)
    f1 <- factor(rep(levels(facto), p))
    f2 <- factor(rep(colnames(X), rep(k, p)))
    
    
if (scale==0){
    if (tri == 0) {
        dotplot(f2 ~ XXXX, groups = f1, key = simpleKey(levels(f1), 
            space = "right"), xlim = c(Min, Max),
xlab = paste("\U00E9","chelle de Likert",sep=""),  
            ...)
    }
    else {
        f2 <- reorder(f2, XXXX, mean)
        dotplot(f2 ~ XXXX, groups = f1, key = simpleKey(levels(f1), 
            space = "right"), xlim = c(Min, Max),
xlab = paste("\U00E9","chelle de Likert",sep=""),  
            ...)
    }
}
else{
if (tri == 0) {
        dotplot(f2 ~ XXXX, groups = f1, key = simpleKey(levels(f1), 
            space = "right"),
xlab = paste("\U00E9","chelle de Likert",sep=""),  ...)
    }
    else {
        f2 <- reorder(f2, XXXX, mean)
        dotplot(f2 ~ XXXX, groups = f1, key = simpleKey(levels(f1), 
            space = "right"), 
xlab = paste("\U00E9","chelle de Likert",sep=""),  ...)
    }


}


}