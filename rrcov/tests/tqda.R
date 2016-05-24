## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))

dodata <- function(method) {

    options(digits = 5)
    set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("===================================================\n")

    data(hemophilia);   show(QdaCov(as.factor(gr)~., data=hemophilia, method=method))
    data(anorexia, package="MASS");     show(QdaCov(Treat~., data=anorexia, method=method))
    data(Pima.tr, package="MASS");      show(QdaCov(type~., data=Pima.tr, method=method))
    data(iris);                         # show(QdaCov(Species~., data=iris, method=method))
    data(crabs, package="MASS");        # show(QdaCov(sp~., data=crabs, method=method))

    show(QdaClassic(as.factor(gr)~., data=hemophilia))
    show(QdaClassic(Treat~., data=anorexia))
    show(QdaClassic(type~., data=Pima.tr))
    show(QdaClassic(Species~., data=iris))
##    show(QdaClassic(sp~., data=crabs))
    cat("===================================================\n")
}


## -- now do it:
dodata(method="mcd")
dodata(method="m")
dodata(method="ogk")
dodata(method="sde")
