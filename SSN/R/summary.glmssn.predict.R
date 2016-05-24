summary.glmssn.predict <- function(object, ...)
{
    print(object)
}


print.glmssn.predict <- function(x, ...)
{
    cat("Object of class glmssn.predict\n")
    Preds <- getPreds(x)
    summary(Preds[,2])
}

