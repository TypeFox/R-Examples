"predict.knncat" <- 
function(object, train, newdata, train.classcol = 1, newdata.classcol = 1, 
         return.classes = TRUE, more = FALSE, verbose = 0, ...)
{
#
# Perform knncat-type classification on new data.
#
#
# First ensure that all the predictor variables in the original appear in the
# newdata.
#
vars <- names(object$vars)
factor.vars <- object$vars == "factor"
not.found <- vars[!is.element (vars, names(newdata))]
if (length(not.found) > 0)
    stop (paste ("Variables ", paste (not.found, collapse=", "),
          "not found in newdata."))
#
# Ensure that the sets of levels match up. This might be stronger than
# we need.
#
train <- data.frame (train[,train.classcol], train[,vars,drop=FALSE])
if (is.factor (train[,1]))
    class.labels = levels(train[,1])
else
    class.labels = seq (0, max(train[,1]))
if (newdata.classcol <= 0)
    newdata <- data.frame (class = rep (0, nrow(newdata)), 
               newdata[,vars,drop=FALSE])
else
    newdata <- data.frame (newdata[,newdata.classcol], newdata[,vars,drop=FALSE])
newdata.true.class <- newdata[,1]
if (any (factor.vars))
{
    train.names <- sapply (object$phi[factor.vars,drop=FALSE], names)
    newdata.names <- sapply (newdata[,vars[factor.vars],drop=FALSE], levels)
    if (!identical (TRUE, all.equal (train.names, newdata.names)))
        stop ("Some level names differ in train and newdata.")
}
#
# Convert categoricals to numerics in the usual way
#
for (i in 1:ncol(train))
{
    if (is.factor (train[,i]))
    {
        train[,i] <- as.numeric (train[,i]) - 1
        newdata[,i] <- as.numeric (newdata[,i]) - 1
    }
}
nrow.train <- nrow(train)
nrow.newdata <- nrow(newdata)
ncol.train <- ncol(train)
ncol.newdata <- ncol(newdata)
train <- c(t(matrix (unlist (train), nrow.train, ncol.train)))
newdata <- c(t(matrix (unlist (newdata), nrow.newdata, ncol.newdata)))
cats.in.var <- sapply (object$phi, length)
cum.cats.this.subset <- c(0, cumsum (cats.in.var)[-length(cats.in.var)])
cdata <- rep (1, ncol.train - 1)
phidata <- unlist (object$phi)
prior.ind <- object$prior.ind
priordata <- object$prior
number.of.classes <- length(priordata)
if (any (names(object) == "knot.values"))
   knots <- c(object$knot.values)
else
    knots <- 0
error.rate <- 0
increase <- numeric (ncol.train - 1)
increase[!factor.vars] <- -5
increase[factor.vars] <- -1
if (return.classes)
    classes <- numeric (nrow.newdata)
else
    classes <- 0
status <- 0
thang <- .C ("donnwrap", 
    as.double (train), as.integer(nrow.train), as.integer(ncol.train),
    as.double (newdata), as.integer(nrow.newdata), as.integer(ncol.newdata),
    as.integer (cats.in.var), as.integer(cum.cats.this.subset),
    as.double (cdata), as.double (phidata), as.integer (prior.ind), 
    as.double (priordata), 
    as.integer (number.of.classes), as.integer (increase), 
    as.double (knots), rate = as.double (error.rate), 
    as.integer (object$best.k),
    as.integer (return.classes), classes = as.integer (classes), 
    as.integer (verbose), as.integer (status), PACKAGE="knncat")
if (more == TRUE)
    cat ("Test set error rate is ", 
        paste (signif (100* thang$rate, 3), "%", sep=""), "\n")
if (return.classes)
{
    preds <- factor (class.labels[thang$classes + 1])
    return (preds)
}
}
