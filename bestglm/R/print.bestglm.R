print.bestglm <-
function (x, ...)
{
ti<-x$Title
cat(ti, fill=TRUE)
if ((x$ModelReport$Bestk>0) || (x$ModelReport$IncludeInterceptQ)){
    cat("Best Model:", fill=TRUE)
    if (any(x$ModelReport$NumDF > 1))
        out<-summary(aov(x$BestModel))
    else
        out<-summary(x$BestModel)$coefficients
    print(out)
    }
else
    cat("Best Model is the null model with no parameters.",fill=TRUE)
}
