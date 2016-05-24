"knncat" <-
function (train, test, k = c(1, 3, 5, 7, 9), xvals = 10, xval.ceil = -1, 
    knots = 10, 
    prior.ind = 4, prior,
    permute = 10, permute.tail = 1, improvement = .01, ridge = .003, 
    once.out.always.out = FALSE, classcol = 1, verbose = 0)
{
#
#  knncat: Create a knncat object from a training set, and optionally 
#          compute the predictions for a test set.
#
# Arguments: 
#        train: data frame of training data, classification in classcol column
#         test: data frame of test data (can be omitted). This should have
#               the correct classification in the classcol column, too.
#            k: vector of choices for # nn's. Default c(1, 3, 5, 7, 9)
#        xvals: number of cross-validations to use to find the best model
#               size and number of nn's. Default 10.
#    xval.ceil: Maximum number of variables to add. -1 = Use the smallest
#               number from any xval; 0 = use the smallest number from the
#               first xval; >= 0, use that.
#        knots: vector of number of knots for numeric variables. Reused if
#               necessary. Default: 10 for each.
#    prior.ind: Integer telling how to compute priors. 1 = estimated from
#               training set; 2 = all equal; 3 = supplied in "prior"; 4 = 
#               ignored. Default: 1.
#        prior: Numeric vector, one entry per unique element in the training
#               set's first column, giving prior probabilities. Ignored unless
#               prior.ind = 3; then they're normalized to sum to 1 and each
#               entry must be strictly > 0.
#      permute: Number of permutations for variable selection. Default: 10.
# permute.tail: A variable fails the permutation test if permute_tail or more
#               permutations do better than the original. Default: 1.
#  improvement: Minimum improvement for variable selection. Ignored unless
#               present and permute missing, or permute = 0; then default = .01.
#        ridge: Amount by which to "ridge" the W matrix for numerical 
#               stability. Default: .003.
# once.out.always.out: if TRUE, a variable that fails a permutation test
#               or doesn't improve by enough is excluded from further
#               consideration during that cross-validation run. Default FALSE.
#     classcol: Column with classification in it. Default: 1.
#      verbose: Controls level of diagnostic output. Higher numbers produce
#               more output, sometimes 'way too much. 0 produces no output;
#               1 gives progress report for xvals. Default: 1.
#
#
# Save numbers of rows and columns. If test is missing, pass zeros for
# that data. If not, test.classes will hold the predictions on return.
#
train.data.name <- deparse (substitute(train))
nrow.train <- nrow(train)
ncol.train <- ncol(train)
if (classcol != 1){
    classcol.name <- names(train)[classcol]
    train <- data.frame (train[,classcol], train[,-classcol])
    names(train)[1] <- classcol.name
}
factor.vars <- sapply (train[,-1,drop=FALSE], is.factor)
#
# "Missing.values" holds the values we'll use in case there are any
# missing values in the data. We compute these in C, but we may as well
# do it here. Whether there are missings in the training set or not, we 
# need to save them in case the test set has missings.
#
missing.values <-  train[1,-1]
missing.values[!factor.vars] <- sapply (train[,-1,drop=FALSE][,!factor.vars], mean,
                                na.rm=TRUE)
most.common <- function (x) {
    tbl <- table (x)
    names(tbl)[tbl == max(tbl)][1]
}
missing.values[factor.vars] <- sapply (train[,-1,drop=FALSE][,factor.vars], most.common)
empty.levels <- sapply (train, function (x) is.factor (x) && any (table (x) == 0))
if (any (empty.levels)) stop ("Some factor has empty levels")
vars.with.na <- sapply (train, function(x) any (is.na(x)))
if (vars.with.na[1] == TRUE)
    stop ("No missing values allowed in response variable")
if (any (vars.with.na))
    for (i in 1:(ncol.train - 1))
        if (vars.with.na[i] == TRUE)
            train[is.na (train[,i+1]),i + 1] <- missing.values[i]
    
train.levs <- lapply (train[,-1], levels)
if (missing (test))
{
    nrow.test <- ncol.test <- test.classes <- 0
}
else
{
    nrow.test <- nrow(test)
    ncol.test <- ncol (test)
    test.classes <- numeric (nrow(test))
#
# Ensure that all the columns of "train" are in "test," too. If they
# are, re-order the columns of test so they match up.
#
    if (any (!is.element (names(train), names(test))))
        stop ("Some train columns aren't present in test.")
    test <- test[,names(train)]
#
# Ensure that the sets of levels for the factor variables are the same.
#
    test.levs <- sapply (test[,-1,drop=FALSE], levels)
    if (!identical (TRUE, all.equal (train.levs, test.levs)))
        stop ("Sets of levels in train and test don't match.")
#
# Handle missing values in the test data. No na's allowed for now.
#
vars.with.na <- sapply (test, function(x) any (is.na(x)))
if (vars.with.na[1] == TRUE)
    stop ("No missing values allowed in test set response variable")
if (any (vars.with.na))
    for (i in 1:(ncol.test - 1))
        if (vars.with.na[i] == TRUE)
            test[is.na (test[,i+1]),i + 1] <- missing.values[i]
}
#
# Find the number in each class. Table() sorts by levels, so we don't
# need to worry about ordering the levels.
#
number.in.class <- table (train[,1])
if (any (number.in.class <= 1))
    stop ("Some classes have only one member. Check \"classcol\"")
nclass <- length(number.in.class)
misclass.mat <- numeric (nclass * nclass)
#
# Identify factor variables. Set "increase",  the vector that indicates
# which variables have which type, accordingly. Also set up cats.in.var
# (whose ith element tells us how many levels are in predictor variable
# i) and cdata (which says whether i is in or out). By the way, if the
# ith element of cats.in.var is -j, variable i is numeric with j knots.
#
increase <- cats.in.var <- cdata <- numeric (ncol(train) - 1)
increase[!factor.vars] <- -5
increase[factor.vars] <- -1
#
# Set up the knots vectors if there are any continuous variables. Even
# though it's a little wasteful, we'll send in a matrix of knots values.
# Each column will contain all the knots for a variable. If the numbers
# of knots differ among columns, some spots will be left as 0's.
#
if (sum (!factor.vars) == 0)
{
    knots.vec <- 0
    knot.values <- 0
}
else
{
    knots.vec <- numeric (sum (!factor.vars))
    knots.vec[] <- knots # re-use as necessary
    knot.values <- matrix (0, max(knots.vec), sum (!factor.vars))
    num.vars <- (1:(ncol(train) - 1))[!factor.vars]
    dimnames(knot.values) <- list (NULL, names(train)[-1][!factor.vars])
#
# Put knot values in columns, so we don't have to transpose this matrix.
# This is a good time to update "train.levs" for numeric.variables.
#
    num.ctr <- 1
    for (i in 1:length(factor.vars))
    {
        if (!factor.vars[i]) 
        {
            qq <- quantile (train[,num.vars[num.ctr]+1], 
                               seq (0, 1, length = knots.vec[num.ctr] + 1))
            qq <- qq[-length(qq)] # chop off the 100% point.
            knot.values[1:length(qq),num.ctr] <- qq
            train.levs[[i]] <- paste ("knot.", 1:knots.vec[num.ctr], sep="")
            num.ctr <- num.ctr + 1
        }
    }
}
#
# Count the number of levels in each factor. phi is the vector of coefficients;
# there's one for each level, plus the right number of knots for each numeric 
# variable.
#
factor.count <- sapply (train[,-1,drop=FALSE][,factor.vars], function(x) length(unique(x)))
if (length(factor.count) == 0)
    factor.count <- 0
cats.in.var[factor.vars] <- factor.count
cats.in.var[!factor.vars] <- knots.vec
phi <- numeric (sum (factor.count) + sum (knots.vec))
#
# Make sure that the k's are valid integers.
#
k <- round (k)
if (any (k <= 0))
{
    warning ("Some invalid k's excluded")
    k <- k[k > 0]
}
k.len <- length(k)
best.k <- 0
if (missing (permute) & !missing(improvement))
    permute <- 0
classif <- 1      # Required by "ords"
status <- 0
#
# Set up prior, if necessary
#
if (prior.ind == 3)
{
    if (missing (prior) || any (prior <= 0))
        stop ("Missing or invalid prior")
    if (length(prior) != nclass)
        stop ("Prior has wrong length")
    prior <- prior / sum (prior)
}
else
    prior <- numeric (nclass)
#
# Convert factors to zero-based numerics. Save the levels for later.
#
for (i in 1:ncol(train))
    if (is.factor (train[,i]))
        train[,i] <- as.numeric (train[,i]) - 1
if (nrow.test > 0)
    for (i in 1:ncol(test))
        if (is.factor (test[,i]))
            test[,i] <- as.numeric (test[,i]) - 1
#
# Right now R has to transpose the data. We should do this in C.
#
train.names <- names(train)
train <- c(t(matrix (unlist (train), nrow(train), ncol(train))))
if (nrow.test == 0)
    test <- 0
else
{
    test <- c(t(matrix (unlist (test), nrow(test), ncol(test))))
    names(test) <- NULL
}
thang <- .C("ord", 
    as.double (train), as.integer (nrow.train), as.integer (ncol.train), 
    as.double (test), as.integer (nrow.test), as.integer (ncol.test), 
    test.classes = as.integer (test.classes),
    cdata = as.double (cdata), phi = as.double (phi),
    as.integer (nclass), as.integer (xvals),
    as.integer (increase), as.integer (permute), as.integer (permute.tail),
    as.double (ridge), as.double (c(knot.values)),
    as.integer (k.len), as.integer (k), 
    best.k = as.integer (best.k),
    as.integer (classif), as.double (improvement), as.integer (cats.in.var),
    as.integer (number.in.class),
    misclass.mat = as.double (misclass.mat),
    as.integer (xval.ceil),
    as.integer (once.out.always.out),
    prior.ind = as.integer (prior.ind), prior = as.double (prior),
    as.integer (verbose), status = as.integer(status),
    PACKAGE = "knncat")
thang <- thang[names(thang) != ""]
thang$misclass.mat <- matrix (thang$misclass.mat, nclass, nclass)
if (nrow.test == 0)
{
    thang$misclass.type <- "train"
    thang$test.classes <- NULL
}
else
{
    thang$misclass.type <- "test"
}
#
# Store names of variables that were actually used.
#
thang$train <- train.data.name
vars <- train.names[-1][thang$cdata != 0]
factor.vars <- factor.vars[thang$cdata != 0]
thang$vars <- rep ("factor", length(factor.vars))
thang$vars[!factor.vars] <- "numeric"
names(thang$vars) <- vars
train.levs <- train.levs[thang$cdata != 0]

phi <- thang$phi[1:(sum (sapply (train.levs, length)))]
#
# Replace levels of numeric variables, which are currently NULL,
# with some useful text.
#
thang$knots.vec <- knots.vec
knots.vec <- knots.vec[thang$cdata != 0]
if (any (!factor.vars))
{
    num.ctr <- 1
    for (i in 1:length(factor.vars))
    {
        if (!factor.vars[i])
        {
            knot.levels <- paste ("knot.", 1:knots.vec[num.ctr], sep="")
            train.levs[[i]] <- knot.levels
            num.ctr <- num.ctr + 1
        }
    }
}
#
# Create list of phis by variable, then attach the names from train.levs.
#
phi.list <- split (phi, rep (names(train.levs), sapply (train.levs, length)))
phi.list <- phi.list[names(train.levs)]
for (i in 1:length(phi.list))
{
    names(phi.list[[i]]) <- train.levs[[i]]
}
thang$phi <- phi.list
thang$build <- numeric (5)
names(thang$build) <- c("Permute", "Improvement", "Once.Out", "Ridge", 
    "Xval")
thang$build[] <- c(permute, improvement, once.out.always.out, ridge, 
     xvals)
thang$k <- k
thang$missing <- missing.values[thang$cdata != 0]
if (any (!factor.vars))
{
    matchers <- is.element (dimnames(knot.values)[[2]], vars)
    if (any (matchers))
        thang$knot.values <- knot.values[,matchers, drop=FALSE]
}
oldClass (thang) <- "knncat"
invisible (return (thang))
}
