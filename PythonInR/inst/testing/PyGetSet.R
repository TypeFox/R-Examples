#' # pyGet
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

#' ## Define some auxilary functions
#' ### buildAscii
# The function build Ascii generates a charcter vector of length len with clen
# number of characters.
# example usage: buildAscii(3, 2)
buildAscii <- function(len, clen=100){
    fun <- function(x) paste(rawToChar(as.raw(sample(33:126, clen, replace=TRUE)), multiple=TRUE), collapse="")
    sapply(1:len, fun)
}

is.identical <- function(x, y){
    if (class(x) == "list" || class(y) == "list"){
        xv <- rapply(x, function(z) z)
        xc <- rapply(x, function(z) class(z))
        xt <- rapply(x, function(z) typeof(z))
        yv <- rapply(y, function(z) z)
        yc <- rapply(y, function(z) class(z))
        yt <- rapply(y, function(z) typeof(z))
    }else{
        xv <- x
        xc <- class(x)
        xt <- typeof(x)
        yv <- y
        yc <- class(y)
        yt <- typeof(y)
    }
    if ( all(xv == yv) & all(xc == yc) & all(xt == yt) ) return(TRUE)
    FALSE
}

# a function to generate test values
genVal <- function(i=-1){
    if (i<0) i <- sample(1:4, 1)
    n <- sample(1:100, 1)
    if (i == 1){
        x <- as.list(as.logical(sample(c(0,1), n, replace=TRUE)))
    }else if (i == 2) {
        x <- as.list(sample(1:1000000, n))
    }else if (i == 3) {
        x <- as.list(rnorm(n) * 1000000)
    }else {
        x <- as.list(sapply(1:n, function(x) intToUtf8(sample(1:1000, 100, replace=TRUE))))
    }
    return(x)
}

genNamedList <- function(difficult=TRUE){
    val <- genVal()
    if (difficult){
        fun <- function(x) intToUtf8(sample(32:256, 10))
    }else{
        fun <- function(x) paste(sample(LETTERS, 10), collapse="")
    }
    nam <- sapply(1:length(val), fun)
    names(val) <- nam
    val
}

#' ## Empty Elements
x <- list(logical(), numeric(), integer(), character(), list(), data.frame(), matrix())
y <- list(NULL, NULL, NULL, NULL, NULL, data.frame(), matrix())
z <- list()
for (i in 1:length(x)){
    pySet(sprintf("r%i", i), x[[i]])
    z[[i]] <- pyGet(sprintf("r%i", i))
}
expect_that(z, equals(y))

#' ## One dimensional elements
x <- list(logical=TRUE, integer=1, double=pi,
          ascii = paste(LETTERS, collapse=" "),
          utf8 = "Some text with some Utf8 characters at the end! äöüß")
for (i in 1:length(x)){
    ##print(sprintf("%s", names(x)[i]))
    pySet(sprintf("r%i", i), x[[i]])
    expect_that(pyGet(sprintf("r%i", i)), equals(x[[i]]))
}

#' ## Objects with length bigger than: 1
x <- list(logical=as.logical(sample(c(0,1), 1000, replace=TRUE)),
          integer=sample(1:1000000, 1000),
          double=rnorm(1000) * 1000000,
          ascii=buildAscii(1000),
          utf8=sapply(1:1000, function(x) intToUtf8(sample(1:1000, 1000, replace=TRUE)))
          )
for (i in 1:length(x)){
    ##print(names(x)[i])
    pySet(sprintf("r%i", i), x[[i]])
    expect_that(pyGet(sprintf("r%i", i)), equals(x[[i]]))
}

#' ## Lists
#' ### simple lists
x <- list(list1=list(1),
          list2=setNames(list(intToUtf8(sample(1:1000, 1000, replace=TRUE))), "charList")
          )
for (i in 1:length(x)){
    pySet(sprintf("r%i", i), x[[i]])
    expect_that(pyGet(sprintf("r%i", i), simplify=FALSE), is_identical_to(x[[i]]))
}

#' ### more complicated lists
# The comparison is a little bit complicated since Python dictonaries change
# the order of the elements therefore the order of the elements must be ignored.
x <- list(A=genVal(1),
          B=genVal(2),
          C=genVal(3),
          D=as.list(buildAscii(1000)),
          E=genVal(4)
         )
for (i in 1:length(x)){
    pySet(sprintf("r%i", i), x[[i]])
    expect_that(pyGet(sprintf("r%i", i), simplify=FALSE), is_identical_to(x[[i]]))
}

#' ### nested lists with names
d <- FALSE
x <- list(genNamedList(d), genNamedList(d), genNamedList(d), 
          genNamedList(d), 
          list(genNamedList(d), 
               list(genNamedList(d)),  genNamedList(d))
          )
expect_that(pySet("r", x), equals(0))
expect_that(is.identical(sort(names(unlist(x))),
                         sort(names(unlist(pyGet("r"))))), equals(TRUE))

#' ## Matrices
x <- matrix(1:8, 4, 2)
rownames(x) <- paste0("row", (1:dim(x)[1]))
colnames(x) <- paste0("col", (1:dim(x)[2]))
expect_that(pySet("r", x), equals(0))
expect_that(pyGet("r"), is_identical_to(x))
M <- as.matrix(cars)
rownames(M) <- paste0("row", (1:dim(M)[1]))
expect_that(pySet("rmatrix", M), equals(0))
expect_that(pyGet("rmatrix"), is_identical_to(M))

#' ## Data.Frame
rownames(cars) <- buildAscii(dim(cars)[1], 5)
expect_that(pySet("r", cars) , equals(0L))
# NOTE: since Python dict changes the order of the columns I can't translate it 1:1
expect_that(pyGet("r")[,colnames(cars)], is_identical_to(cars))
