## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))
library(MASS)

dodata <- function(method) {

    options(digits = 5)
    set.seed(101)

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("===================================================\n")

    data(pottery);      show(lda <- LdaPP(origin~., data=pottery, method=method)); show(predict(lda))
    data(hemophilia);   show(lda <- LdaPP(as.factor(gr)~., data=hemophilia, method=method)); show(predict(lda))
    data(anorexia);     show(lda <- LdaPP(Treat~., data=anorexia, method=method)); show(predict(lda))
    data(Pima.tr);      show(lda <- LdaPP(type~., data=Pima.tr, method=method)); show(predict(lda))
    data(crabs);        show(lda <- LdaPP(sp~., data=crabs, method=method)); show(predict(lda))

    cat("===================================================\n")
}


## -- now do it:

## Commented out - still to slow
##dodata(method="huber")
dodata(method="mad")
##dodata(method="sest")
dodata(method="class")
