mv.pca <- function(X, shape.mat, scores= TRUE, na.action=na.fail, ...)
    {
    covmat <- match.arg(shape.mat,c("cov","sign","symmsign","tyler.shape","duembgen.shape"))
    
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    n<-dim(X)[1]
    p<-dim(X)[2]
    if (is.null(colnames(X))) colnames(X) <- paste(rep("X",p),1:p,sep=".")
    c.names <- colnames(X)
    
    cov.list <- switch(covmat,
              "cov" = {
                      covM <- cov(X)
                      center <- colMeans(X)
                      list(cov=covM, center=center, n.obs=n)
                      }
              ,
              "sign" = {
                       center <- spatial.median(X, ...)
                       covM <- SCov(X, center)
                       list(cov=covM, center=center, n.obs=n)
                       } 
              ,
              "symmsign" = {
                           center <- spatial.median(X, ...)
                           covM <- SSCov(X)
                           list(cov=covM, center=center, n.obs=n)
                           } 
              ,
              "tyler.shape" = {
                               Y <- HR.Mest(X, ...) 
                               center <- Y$center
                               covM <- Y$scatter
                               list(cov=covM / sum(diag(covM)), center=center, n.obs=n)
                               }
              ,
              "duembgen.shape" = {
                           center <- spatial.median(X, ...)
                           covM <- duembgen.shape(X)
                           list(cov=covM / sum(diag(covM)), center=center, n.obs=n)
                           }
                      )
    
    res <- princomp(covmat=cov.list)
    attr(res$loadings,"dimnames")[[1]] <- c.names
    if (scores) res$scores <- predict(res,X)
    return(res)
    }
