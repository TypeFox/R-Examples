intrinsic.pars <-
function(y = y, data = parent.frame(), id = id, repeated = NULL, rscale="ordinal")
{
    call <- match.call() 
    mcall <- match.call(expand.dots=FALSE)
    mf <- match(c("y", "data", "id", "repeated"), names(mcall), 0L)
    m <- mcall[c(1L, mf)]
    if (is.null(m$id)) 
        m$id <- as.name("id")
    m$formula <- y~1
    m[[1]] <- as.name("model.frame")
    m <- eval(m, envir = parent.frame())
    Terms <- attr(m, "terms") 
    if(attr(Terms,"intercept")!=1) 
       stop("an intercept must be included")
    Y <- as.numeric(factor(model.response(m)))
    if (is.null(Y)) {
        stop("response variable not found")
    }
    ncategories <- nlevels(factor(Y))
    if (ncategories <= 2) 
        stop("The response variable should have more than 2 categories")
    id <- model.extract(m, "id")
    if (is.null(id)) {
        stop("'id' variable not found")
    }
    if (length(id) != length(Y)) 
        stop("response variable and 'id' are not of same length")
    repeated <- model.extract(m, "repeated")
    if (is.null(repeated)) {
        index <- order(unlist(split(1:length(id),id)))
        repeated <- c(unlist(sapply(unlist(lapply(split(id, id), length)), function(x) 1:x)))
        repeated <- repeated[index]
    }
    if (length(repeated) != length(Y)) 
        stop("response variable and 'repeated' are not of same length")
    id <- as.numeric(factor(id))
    repeated <- as.numeric(factor(repeated))
    if(all(id==repeated)) 
         stop("'repeated' and 'id' must not be equal")
    dummy <- split(repeated, id)
    if (any(unlist(lapply(dummy, length)) != unlist(lapply(lapply(dummy, 
        unique), length)))) 
        stop("'repeated' does not have unique values per 'id'")
cdata <- datacounts(Y,id,repeated,ncategories)
if(rscale=="ordinal"){
cmod <- gnm(counts~(factor(x)+factor(y))*factor(tp)+factor(tp):x:y,
                                             family=poisson,data=cdata)
ans <- as.vector(coef(cmod)[pickCoef(cmod,"x:y")])
                     } else {
ans <- rep(0,max(cdata$tp))
cdata$x <- factor(cdata$x)
cdata$y <- factor(cdata$y)
for(i in 1:max(cdata$tp))
{
cmod <- gnm(counts~x+y+MultHomog(x,y),
       family=poisson,data=cdata[cdata$tp==i,])
cscores <- coef(cmod)[pickCoef(cmod,"MultHomog")]
cscores <- c(tcrossprod(normscores(cscores)))
cmod <- gnm(counts~factor(x)+factor(y)+cscores,
        family=poisson,data=cdata[cdata$tp==i,])
ans[i] <- coef(cmod)["cscores"]
}
}
ans
}
