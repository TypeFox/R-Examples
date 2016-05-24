## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))
library(MASS)

dodata <- function(method) {

    options(digits = 5)
    set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("===================================================\n")

    cat("\nData: ", "anorexia\n");
    data(hemophilia);
    show(Linda(as.factor(gr)~., data=hemophilia, method=method))
    show(Linda(as.factor(gr)~., data=hemophilia))

    cat("\nData: ", "anorexia\n");
    data(anorexia);
    show(Linda(Treat~., data=anorexia, method=method))
    show(Linda(Treat~., data=anorexia))

    cat("\nData: ", "Pima\n");
    data(Pima.tr);
    show(Linda(type~., data=Pima.tr, method=method))
    show(Linda(type~., data=Pima.tr))

##    cat("\nData: ", "iris\n");
##    data(iris);
##    show(Linda(Species~., data=iris, method=method))
##    show(Linda(Species~., data=iris))

    cat("\nData: ", "crabs\n");
    data(crabs);
    show(Linda(sp~., data=crabs, method=method))
    show(Linda(sp~., data=crabs))

    cat("===================================================\n")
}


## -- now do it:
dodata(method="mcdA")
dodata(method="mcdB")
dodata(method="mcdC")
#dodata(method="fsa")
