#' Letter summary of similarities and differences
#' 
#' Convert a logical vector or a vector of p-values or a correlation or
#' distance matrix into a character-based display in which common characters
#' identify levels or groups that are not significantly different.  Designed
#' for use with the output of functions like TukeyHSD, diststats, simint,
#' simtest, csimint, csimtestmultcomp, friedmanmc, kruskalmcpgirmess.
#' 
#' Produces a "Letter-Based Representation of All- Pairwise Comparisons" as
#' described by Piepho (2004).  (The present algorithm does NOT perform his
#' "sweeping" step.)  \code{multcompLettersx} are wrapper of multcompLetters
#' that will reorder the levels of the data so that the letters appear in a
#' descending order of the mean. \code{mulcompletters3} is similar to
#' \code{multcompletters2} except that it uses vector names to separte and the
#' later has an formula interface. \code{multcompLetters4} will take a aov or
#' lm object and a comparison test and will produce all the letters for the
#' terms and interactions.
#' 
#' @aliases multcompLetters multcompLetters2 multcompLetters3 multcompLetters4
#' @param x One of the following: (1) A square, symmetric matrix with row
#' names.  (2) A vector with hyphenated names, which identify individual items
#' or factor levels after "strsplit".  (3) An object of class "dist".  If x (or
#' x[1]) is not already of class "logical", it is replaced with
#' do.call(compare, list(x, threshold)), which by default converts numbers
#' (typically p-values) less than 0.05 to TRUE and everything else to FALSE.
#' If x is a matrix, its diagonal must be or must convert to FALSE.
#' @param compare function or binary operator; not used if class(x) is
#' "logical".
#' @param threshold Second (reference) argument to "compare".
#' @param Letters Vector of distinct characters (or character strings) used to
#' connect levels that are not significantly different.  They should be
#' recogizable when concatonated.  The last element of "Letters" is used as a
#' prefix for a reuse of "Letters" if more are needed than are provided.  For
#' example, with the default "Letters", if 53 distinct connection colums are
#' required, they will be "a", ..., "z", "A", ..., "Z", and ".a".  If 54 are
#' required, the last one will be ".b".  If 105 are required, the last one will
#' be "..a", etc.  (If the algorithm generates that many distinct groups, the
#' display may be too busy to be useful, but the algorithm shouldn't break.)
#' @param reversed A logical value indicating whether the order of the letters
#' should be reversed. Defaults to FALSE.
#' @param formula The formula used to make the test (lm, aov, glm, etc.).  Like
#' y ~ x.
#' @param data Data used to make the test.
#' @param y Value of the response variable.
#' @param z Categorical variables used in the test.
#' @param object An object of class aov or lm for the time being.
#' @param comp A object with multiple comparison or a function name to perform
#' a multiple comparison.
#' @param ...  Extra arguments passed to multcompLetters.
#' @return An object of class 'multcompLetters', which is a list with the
#' following components: \item{Letters }{character vector with names = the
#' names of the levels or groups compared and with values = character strings
#' in which common values of the function argument "Letters" identify levels or
#' groups that are not significantly different (or more precisely for which the
#' corresponding element of "x" was FALSE or was converted to FALSE by
#' "compare").  } \item{monospacedLetters }{ Same as "Letters" but with spaces
#' so the individual grouping letters will line up with a monspaced type font.
#' } \item{LetterMatrix }{Logical matrix with one row for each level compared
#' and one column for each "Letter" in the "letter-based representation".  The
#' output component "Letters" is obtained by concatonating the column names of
#' all columns with TRUE in that row.  } multcompLetters4 will return a named
#' list with the terms containing a object of class 'multcompLetters' as
#' produced by \code{multcompLetters}.
#' @author Spencer Graves, Hans-Peter Piepho and Luciano Selzer
#' @seealso \code{\link{multcompBoxplot}} \code{\link{plot.multcompLetters}}
#' \code{\link{print.multcompLetters}} \code{\link{multcompTs}}
#' \code{\link{vec2mat}}
#' @references Piepho, Hans-Peter (2004) "An Algorithm for a Letter-Based
#' Representation of All-Pairwise Comparisons", Journal of Computational and
#' Graphical Statistics, 13(2)456-466.
#' @keywords dplot
#' @export
#' @examples
#' 
#' ##
#' ## 1.  a logical vector indicating signficant differences
#' ##
#' dif3 <- c(FALSE, FALSE, TRUE)
#' names(dif3) <- c("A-B", "A-C", "B-C")
#' dif3L <- multcompLetters(dif3)
#' dif3L
#' print(dif3L)
#' print(dif3L, TRUE)
#' 
#' ##
#' ## 2.  numeric vector indicating statistical significance
#' ##
#' dif4 <- c(.01, .02, .03, 1)
#' names(dif4) <- c("a-b", "a-c", "b-d", "a-d")
#' (diff4.T <- multcompLetters(dif4))
#' 
#' (dif4.L1 <- multcompLetters(dif4,
#'        Letters=c("*", ".")))
#' # "Letters" can be any character strings,
#' # but they should be recognizable when
#' # concatonated.
#' 
#' ##
#' ## 3.  distance matrix
#' ##
#' dJudge <- dist(USJudgeRatings)
#' dJl <- multcompLetters(dJudge, compare='>', threshold = median(dJudge))
#' # comparison of 43 judges;  compact but undecipherable:
#' dJl
#' 
#' x <- array(1:9, dim=c(3,3),
#'    dimnames=list(LETTERS[1:3], NULL) )
#' d3 <- dist(x)
#' dxLtrs <- multcompLetters(d3, compare=">", threshold=2)
#' 
#' d3d <- dist(x, diag=TRUE)
#' dxdLtrs <- multcompLetters(d3d, compare=">", threshold=2)
#' 
#' \dontshow{stopifnot(}
#' all.equal(dxLtrs, dxdLtrs)
#' \dontshow{)}
#' 
#' d3u <- dist(x, upper=TRUE)
#' dxuLtrs <- multcompLetters(d3d, compare=">", threshold=2)
#' 
#' \dontshow{stopifnot(}
#' all.equal(dxLtrs, dxuLtrs)
#' \dontshow{)}
#' 
#' ##
#' ## 4.  cor matrix
#' ##
#' set.seed(4)
#' x100 <- matrix(rnorm(100), ncol=5,
#'                dimnames=list(NULL, LETTERS[1:5]) )
#' cx <- cor(x100)
#' cxLtrs <- multcompLetters(abs(cx), threshold=.3)
#' 
#' 
#' ##
#' ##5. reversed
#' ##
#' 
#' dif3 <- c(FALSE, FALSE, TRUE)
#' names(dif3) <- c("A-B", "A-C", "B-C")
#' dif3L <- multcompLetters(dif3)
#' dif3L.R <- multcompLetters(dif3, rev = TRUE)
#' dif3L
#' dif3L.R
#' 
#' 
#' ##
#' ##6. multcompletters2 usage
#' 
#' experiment <- data.frame(treatments = gl(11, 20, labels = c("dtl", "ctrl", "treat1", 
#'               "treat2", "treatA2", "treatB", "treatB2",
#'               "treatC", "treatD", "treatA1", "treatX")),
#'               y = c(rnorm(20, 10, 5), rnorm(20, 20, 5), rnorm(20, 22, 5), rnorm(20, 24, 5),
#'                rnorm(20, 35, 5), rnorm(20, 37, 5), rnorm(20, 40, 5), rnorm(20, 43, 5),
#'                rnorm(20, 45, 5), rnorm(20, 60, 5), rnorm(20, 60, 5)))
#' exp_tukey <- TukeyHSD(exp_aov <- aov(y  ~ treatments, data = experiment))
#' exp_letters1 <- multcompLetters(exp_tukey$treatments[,4])
#' exp_letters1
#' #Notice lowest mean treatments gets a "e"
#' #Ordered letters
#' multcompLetters2(y ~ treatments, exp_tukey$treatments[,"p adj"], experiment)
#' multcompLetters2(y ~ treatments, exp_tukey$treatments[,"p adj"], experiment, reversed = TRUE)
#' 
#' ##7. multcompletters3 usage
#' 
#' multcompLetters3("treatments", "y", exp_tukey$treatments[,"p adj"], experiment)
#' 
#' ##8. multcompletters4 usage
#' 
#' 
#' multcompLetters4(exp_aov, exp_tukey)
#' 
#' 
"multcompLetters" <-
function(x, compare="<",
   threshold=0.05, Letters=c(letters, LETTERS, "."),
   reversed = FALSE){
##
## 1.  Covert to logical
##
  x.is <- deparse(substitute(x))
  if(class(x)=="dist")x <- as.matrix(x)  
  if(!is.logical(x))
    x <- do.call(compare, list(x, threshold))
##
## 2.  Create array of distinct pairs
##
  dimx <- dim(x)
  {
    if((length(dimx)==2) && (dimx[1]==dimx[2])){
      Lvls <- dimnames(x)[[1]]
      if(length(Lvls)!=dimx[1])
        stop("Names requred for ", x.is)
      else{
#       Create a matrix with 2 columns
#       with the names of all pairs         
        x2. <- t(outer(Lvls, Lvls, paste,
                     sep=""))
        x2.n <- outer(Lvls, Lvls,
           function(x1, x2)nchar(x2))
        x2.2 <- x2.[lower.tri(x2.)]
        x2.2n <- x2.n[lower.tri(x2.n)]
        x2a <- substring(x2.2, 1, x2.2n)
        x2b <- substring(x2.2, x2.2n+1)
        x2 <- cbind(x2a, x2b)
        x <- x[lower.tri(x)]        
      }
    }
    else{  
      namx <- names(x)
      if(length(namx)!=length(x))
        stop("Names required for ", x.is)
      x2 <- vec2mat2(namx)
      Lvls <- unique(as.vector(x2))
    }
  }
##
## 3.  Find the names of the levels 
##  
  n <- length(Lvls)
#   Generate an initial column
  LetMat <- array(TRUE, dim=c(n, 1),
               dimnames=list(Lvls, NULL))
##
## 4.  How many distinct pairs?  
##  
  k2 <- sum(x)
  if(k2==0){
    Ltrs <- rep(Letters[1], n)
    names(Ltrs) <- Lvls
    dimnames(LetMat)[[2]] <- Letters[1]
    return(list(Letters=Ltrs,
                LetterMatrix=LetMat))  
  }
##
## 4.  At last 2 levels are different: 
##     insert & absorb
##  
  distinct.pairs <- x2[x,,drop=FALSE]
  absorb <- function(A.){
#    Do the work in a recursive function:      
#    Delete any column for which the TRUE 
#    connections are a subset of another column
    k. <- dim(A.)[2]
    if(k.>1){ #i. <- 1; j. <- 2
      for(i. in 1:(k.-1))for(j. in (i.+1):k.){
        if(all(A.[A.[, j.], i.])){
#### drop a redundant column and recurse ###
          A. <- A.[, -j., drop=FALSE]
          return(absorb(A.))
        }
        else {
          if(all(A.[A.[, i.], j.])){
#### drop a redundant column and recurse ###
            A. <- A.[, -i., drop=FALSE]
            return(absorb(A.))
          }
        }          
      }
    }
#### end internal function absorb #######      
    A.
  }
# Now apply this function 
  for(i in 1:k2){ # i <- 1+i
#     Process the distinct differences one at a time       
#     Insert    i <- 1+i
#     Are (distinct) levels Td2[i, 1] and Td2[i,2]
#        connected in any columns of A?
    dpi <- distinct.pairs[i,]
    ijCols <- (LetMat[dpi[1],] & LetMat[dpi[2], ])
    if(any(ijCols)){
#     "Insert":  Break this connection 
      A1 <- LetMat[, ijCols, drop=FALSE]
      A1[dpi[1],] <- FALSE
      LetMat[dpi[2], ijCols] <- FALSE
      LetMat <- cbind(LetMat, A1)
#     Absorb   A. <- A
      LetMat <- absorb(LetMat)
    }
  }
##
## 5.  Sort the columns for visual appeal 
##  
  sortCols <- function(B){
    firstRow <- apply(B, 2, function(x)which(x)[1])
    B <- B[, order(firstRow)]
#     If ties, sort submatrices
    firstRow <- apply(B, 2, function(x)which(x)[1])
    reps <- (diff(firstRow)==0)
    if(any(reps)){
#     Break ties
      nrep <- table(which(reps))
      irep <- as.numeric(names(nrep))
      k <- dim(B)[1]
      for(i in irep){
        i. <- i:(i+nrep[as.character(i)])
        j. <- (firstRow[i]+1):k
        B[j., i.] <- sortCols(B[j., i., drop=FALSE])
      }
    }
#### end internal function sortCols #######      
    B
  }
  LetMat. <- sortCols(LetMat)
### Should the letters go in the reversed order?
  if(reversed) LetMat. <- LetMat.[ ,rev(1:ncol(LetMat.))]
# DON'T Sweep
    #...
##
## 6.  Create "Letters" for column names
##
  k.ltrs <- dim(LetMat.)[2]
  makeLtrs <- function(kl, ltrs=Letters){
    kL <- length(ltrs)
    if(kl<kL)return(ltrs[1:kl])
    ltrecurse <- c(paste(ltrs[kL], ltrs[-kL],
            sep=""), ltrs[kL])
    c(ltrs[-kL], makeLtrs(kl-kL+1,
                          ltrecurse))
  }
  Ltrs <- makeLtrs(k.ltrs, Letters)
  dimnames(LetMat.)[[2]] <- Ltrs
##
## 7.  Create simple summaries
##
  LetVec <- rep(NA, n)
  names(LetVec) <- Lvls
  for(i in 1:n)
    LetVec[i]<- paste(Ltrs[LetMat.[i, ]],
                    collapse="")
  nch.L <- nchar(Ltrs)
# To allow for multicharacter "Letters", create
# a vector of blanks with the right number
# of characters for each.  
  blk.L <- rep(NA, k.ltrs)
  for(i in 1:k.ltrs)
    blk.L[i] <- paste(rep(" ", nch.L[i]), collapse="")
# Now create monospacedLetters:    
  monoVec <- rep(NA, n)
  names(monoVec) <- Lvls
  for(j in 1:n){
    ch2 <- blk.L
    if(any(LetMat.[j,]))
      ch2[LetMat.[j,]] <- Ltrs[LetMat.[j,]]
    monoVec[j] <- paste(ch2, collapse="")
  }
##
## 8.  done
##
  InsertAbsorb <- list(Letters=LetVec,
        monospacedLetters=monoVec, 
        LetterMatrix=LetMat.)
  class(InsertAbsorb) <- "multcompLetters"
  InsertAbsorb  
}

#' @export
#' @describeIn multcompLetters 
"multcompLetters2" <- 
  function (formula, x, data, ...) {
    #Convert formula to character, get rid of "~"
    fm <- as.character(formula)
    fm <- fm[-1]
    #Split char vector with ":" as this points an
    #interaction and is not included in the data
    #per se
    fm <- strsplit(fm, ":", fixed = TRUE)
    y.z  <- tapply(data[,fm[[1]]], data[,fm[[2]]], 
                   function(x) do.call(mean, list(x=x)))
    oz <- order(y.z, decreasing= T )
    #This is to handle interactions
    if (length(fm[[2]] > 1)) {
      Lvls <- levels(interaction(data[,fm[[2]]], sep = ":"))[oz]
    } else {
      Lvls <- levels(data[,fm[[2]]])[oz]
    }
    value <- vec2mat(x)
    value <- value[Lvls, Lvls]
    multcompLetters(value, ...)
  }

#' @export
#' @describeIn multcompLetters
"multcompLetters3" <- 
  function (z, y , x, data, ...) {
    y.z  <- tapply(data[, y], data[, z], 
                   function(x) do.call(mean, list(x=x)))
    oz <- order(y.z, decreasing= T )
    #This is to handle interactions
    if (length(z > 1)) {
      Lvls <- levels(interaction(data[, z], sep = ":"))[oz]
    } else {
      Lvls <- levels(data[, z])[oz]
    }
    value <- vec2mat(x)
    value <- value[Lvls, Lvls]
    multcompLetters(value, ...)
  }
#' @export
#' @describeIn multcompLetters
#' 
"multcompLetters4" <- 
  function (object, comp, ...) {
    #Extract needed data from object
    formula <- terms(object)
    Terms <- colnames(attr(terms(object), "factors"))
    data <- model.frame(object)
    fm <- as.character(formula)
    fm <- fm[-1]
    fms <- list()
    for (i in 1:length(Terms)){
      fms[[i]] <- formula(paste(fm[1], "~", Terms[i]))
    }
    names(fms) <- Terms
    if(is.character(comp) | is.symbol(comp)) {
      comp <- match.fun(comp)
      comp <- comp(object)
    }
    comp <- extract_p(comp)
    ans <- list()
    for(i in 1:length(Terms)){
      ans[[i]] <- list(formula = fms[[i]], p = comp[[i]])
    }
    names(ans) <- Terms
    lapply(ans, function(x) multcompLetters2(x$formula, x$p, data, ...))
  }
