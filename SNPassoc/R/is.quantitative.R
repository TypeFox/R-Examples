`is.quantitative` <-
function(formula, data){
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    ans<-length(table(unique(mf[,1])))>2
    ans
}

