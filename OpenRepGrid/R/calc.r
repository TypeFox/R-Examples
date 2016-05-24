
#' Descriptive statistics for constructs and elements of a grid.
#'
#' @param x       \code{repgrid} object.
#' @param index   Whether to print the number of the element. 
#' @param trim    The number of characters an element or a construct
#'                is trimmed to (default is
#'                \code{20}). If \code{NA} no trimming occurs. Trimming
#'                simply saves space when displaying correlation of constructs
#'                or elements with long names.
#' @return        A dataframe containing the following measures is returned invisibly 
#'                (see \code{\link{describe}}): \cr
#'                item name \cr
#'                item number \cr 
#'                number of valid cases \cr  
#'                mean standard deviation \cr
#'                trimmed mean (with trim defaulting to .1) \cr
#'                median (standard or interpolated) \cr
#'                mad: median absolute deviation (from the median) \cr
#'                minimum \cr
#'                maximum \cr
#'                skew \cr
#'                kurtosis \cr
#'                standard error \cr
#'
#' @note          Note that standard deviation and variance are estimations, 
#'                i.e. including Bessel's correction. For more info type \code{?describe}.
#'
#' @author        Mark Heckmann
#' @export
#' @aliases statsElements statsConstructs 
#' @rdname stats
#' @examples \dontrun{
#'
#'    statsConstructs(fbb2003)
#'    statsConstructs(fbb2003, trim=10)
#'    statsConstructs(fbb2003, trim=10, index=F)
#'
#'    statsElements(fbb2003)
#'    statsElements(fbb2003, trim=10)
#'    statsElements(fbb2003, trim=10, index=F)
#'
#'    # save the access the results
#'    d <- statsElements(fbb2003)
#'    d
#'    d["mean"]
#'    d[2, "mean"]  # mean rating of 2nd element
#'  
#'    d <- statsConstructs(fbb2003)
#'    d
#'    d["sd"]
#'    d[1, "sd"]   # sd of ratings on first construct
#' }
#'
statsElements <- function(x, index=TRUE, trim=20)
{
  if (!inherits(x, "repgrid")) 							      # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  s <- getRatingLayer(x)
  res <- describe(s)                              # psych function
  enames <- getElementNames2(x, index=index, trim=trim) 
  ne <- getNoOfElements(x)
  if (length(unique(enames)) != ne){
    stop("please chose a longer value for 'trim' or set 'index' to TRUE", 
         "as the current value produces indentical rownames")
  }  
  rownames(res) <- enames
  class(res) <- c("statsElements", "data.frame")
  return(res)
}
# compare Bell, 1997, p. 7


#' Print method for class statsElements
#' 
#' @param x         Object of class statsElements.
#' @param digits    Numeric. Number of digits to round to (default is 
#'                  \code{1}).
#' @param ...       Not evaluated.
#' @export
#' @method          print statsElements
#' @keywords        internal
#'
print.statsElements <- function(x, digits=2, ...) 
{
  cat("\n##################################")
  cat("\nDesriptive statistics for elements")
  cat("\n##################################\n\n")
  x <- as.data.frame(x)
  print(round(x, 2))
}


#' Descriptive statistics for constructs and elements of a grid.
#'
#' @note          Note that standard deviation and variance are estimated ones, 
#'                i.e. including Bessel's correction. For more info type \code{?describe}.
#'
#' @author        Mark Heckmann
#' @export
#' @rdname stats
#'
statsConstructs <- function(x, index=T, trim=20){
  if (!inherits(x, "repgrid")) 							      # check if x is repgrid object
   	stop("Object x must be of class 'repgrid'")
  s <- getRatingLayer(x)
  res <- describe(t(s))
  cnames <- getConstructNames2(x, index=index, trim=trim)
  rownames(res) <- cnames
  class(res) <- c("statsConstructs", "data.frame")
  return(res)
}
# compare Bell, 1997, p. 9


#' Print method for class statsConstructs
#' 
#' @param x         Object of class statsConstructs.
#' @param digits    Numeric. Number of digits to round to (default is 
#'                  \code{1}).
#' @param ...       Not evaluated.
#' @export
#' @method          print statsConstructs
#' @keywords        internal
#'
print.statsConstructs <- function(x, digits=2, ...) 
{
  cat("\n####################################")
  cat("\nDesriptive statistics for constructs")
  cat("\n####################################\n\n")
  x <- as.data.frame(x)
  print(round(x, 2))
}


# TODO:
# Golden Section Statistics, p. 10
# Adams-Webber (1990) for dichotmous data
# Bell & McGorry for rating grids
#
# Adams-Webber, J.R. (1990) Some fundamental asymmetries in the structure of personal
# constructs. Advances in Personal Constructs, 1, 49-55.
#
# Bell, R. C., & McGorry P. (1992) The analysis of repertory 
# grids used to monitor the perceptions of recovering psychotic 
# patients. In A. Thomson & P. Cummins (Eds) European 
# Perspectives in Personal Construct Psychology. Lincoln UK: 
# European Personal Construct Association.
# --> ask Richard for paper
# statsGoldenSection


getScoreDataFrame <- function(x){
  sc <- x@ratings[,,1]
  rownames(sc) <- getConstructNames(x)$l
  colnames(sc) <- getElementNames(x)
  sc
}

# disc    element to be used as discrepancy
# remove  logical. remove element that was used as discrepancy
#
makeDiscrepancy <- function(x, disc, remove=TRUE){
  sc <- x@ratings[,,1]
  colnames(sc) <- getElementNames(x)
  scDiscElement <- sc[, disc]
  if (remove)
    sc[, -disc] - scDiscElement else 
    sc[,] - scDiscElement
}

statsDiscrepancy <- function(x, disc, sort=TRUE){
  a <- describe(makeDiscrepancy(x, disc))
  if (sort)
    a[order(a$mean),] else
    a 
}



# Measures comparing Constructs: Intensity, Cognitive Complexity 
# and other measures of Grid variation.

# Ward Hierarchical Clustering
# cluster <- function(){
#   
#   
# }
# a <- doubleEntry(x)
# sc <- getScoreDataFrame(x)
# d <- dist(sc, method = "euclidean") # distance matrix
# fit <- hclust(d, method="ward")
# plot(as.dendrogram(fit), horiz = F)
# plot(as.dendrogram(fit), horiz = F)


###############################################################################
# order elements and constructs by angles in first two dimensions from
# singular value decomposition approach (cf. Raeithel ???)
###############################################################################

#' Calculate angles for points in first two columns.
#'
#' The angles of the points given by the values in the first and second
#' column of a matrix seen from the origin are calulated in degrees.
#'
#' @param x           A matrix.
#' @param dim         Dimensions used for calculating angles.
#' @param clockwise   Logical. Positive angles are clockwise with x axis as 
#'                    basis.
#' @return  vector.   The angles of each row point with the origin as reference.
#'
#' @author       Mark Heckmann
#' @keywords internal
#'
#' @examples \dontrun{
#'
#'    calcAngles(matrix(rnorm(9), 3))
#' }
#'
calcAngles <- function(x, dim=c(1,2), clockwise=TRUE){
  angles <- atan2(x[ ,dim[2]], x[ ,dim[1]]) / (2 * pi / 360)
  angles <- angles * -1                    # positive angles are counted clockwise atan2 does anticlockwise by default
  angles[angles < 0] <- 360 + angles[angles < 0]  # map to 0 to 360 degrees i.e. only positive  values for ordering
  if (!clockwise)
    angles <- 360 - angles
  angles
}


#' Make indexes to order grid by angles in given dimensions.
#'
#' Reorder indexes for constructs and elements are calculated 
#' using the coordinates of the given dimensions.
#'
#' @param x           \code{repgrid} object that has been submitted to
#'                    \code{\link{calcBiplotCoords}}.
#' @param dim         Dimensions used to calculate angles for reordering grid.
#' @param clockwise   Logical. Positive angles are clockwise with x axis as 
#'                    basis.
#' @return  A list containing the indexes to reorder the grid. The 
#'          first list element for the constructs, the second for the elements indexes.
#'
#' @author        Mark Heckmann
#' @keywords internal
#' @examples \dontrun{
#'
#'    x <- randomGrid(15,30)      # make random grid
#'    i <- angleOrderIndexes2d(x) # make indexes for ordering
#'    x <- x[i[[1]], i[[2]]]      # reorder constructs and elements
#'    x                           # print grid
#' }
#'                  
angleOrderIndexes2d <- function(x, dim=c(1,2), clockwise=TRUE){
  if (!inherits(x, "repgrid")) 							      # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'") 
  E.mat <- x@calcs$biplot$elem
  C.mat <- x@calcs$biplot$con
  C.angles <- calcAngles(C.mat, dim=dim, clockwise=clockwise)
  E.angles <- calcAngles(E.mat, dim=dim, clockwise=clockwise)
  list(c.order=order(C.angles), 
       e.order=order(E.angles))
}


#' Order grid by angles between construct and/or elements in 2D.
#'
#' The approach is to reorder 
#' the grid matrix by their polar angles on the first two principal 
#' components from a data reduction technique 
#' (here the biplot, i.e. SVD). The function 
#' \code{reorder2d} reorders the grid according to the angles 
#' between the x-axis and the element (construct) 
#' vectors derived from a 2D biplot solution.
#' This approach is apt to identify circumplex 
#' structures in data indicated by the diagonal stripe 
#' in the display (see examples).
#'
#' @param x           \code{repgrid} object.
#' @param dim         Dimension of 2D solution used to calculate angles
#'                    (default \code{c(1,2)}).
#' @param center		  Numeric. The type of centering to be performed. 
#'                    0= no centering, 1= row mean centering (construct), 
#'                    2= column mean centering (elements), 3= double-centering (construct and element means),
#'                    4= midpoint centering of rows (constructs).
#'                    The default is \code{1} (row centering).
#' @param normalize   A numeric value indicating along what direction (rows, columns)
#'                    to normalize by standard deviations. \code{0 = none, 1= rows, 2 = columns}
#'                    (default is \code{0}).
#' @param g           Power of the singular value matrix assigned to the left singular 
#'                    vectors, i.e. the constructs.
#' @param h           Power of the singular value matrix assigned to the right singular 
#'                    vectors, i.e. the elements.
#' @param rc          Logical. Reorder constructs by similarity (default \code{TRUE}).
#' @param re          Logical. Reorder elements by similarity (default \code{TRUE}).
#' @param ...         Not evaluated.
#'
#' @return            Reordered \code{repgrid} object. 
#'
#' @author            Mark Heckmann
#' @export
#'
#' @examples \dontrun{
#'
#'    x <- feixas2004  
#'    reorder2d(x)            # reorder grid by angles in first two dimensions
#'    reorder2d(x, rc=F)      # reorder elements only
#'    reorder2d(x, re=F)      # reorder constructs only
#' }
#'
reorder2d <- function(x, dim=c(1,2), center=1, normalize=0, g=0, h=1-g, 
                    rc=TRUE, re=TRUE, ... ){
  if (!inherits(x, "repgrid")) 							  # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  x <- calcBiplotCoords(x, center=center, normalize=normalize, g=g, h=h, ...)
  i <- angleOrderIndexes2d(x, dim=dim)        # make indexes for ordering
  if(rc)
    x <- x[i[[1]], ,drop=FALSE]
  if(re)
    x <- x[ ,i[[2]], drop=FALSE]
  x
}



  
  
  
###############################################################################


#' Calculate the correlations between elements.
#' 
#' Note that simple element correlations as a measure of similarity
#' are flawed as they are not invariant to construct reflection (Mackay, 1992; 
#' Bell, 2010). A correlation index invariant to construct reflection is 
#' Cohen's rc measure (1969), which can be calculated using the argument 
#' \code{rc=TRUE} which is the default option.
#'
#' @param x         \code{repgrid} object.
#' @param rc        Use Cohen's rc which is invariant to construct 
#'                  reflection (see desciption above). It is used as the default.
#' @param method    A character string indicating which correlation coefficient 
#'                  is to be computed. One of \code{"pearson"} (default), 
#'                  \code{"kendall"} or \code{"spearman"}, can be abbreviated.
#'                  The default is \code{"pearson"}.
#' @param trim      The number of characters a construct is trimmed to (default is
#'                  \code{20}). If \code{NA} no trimming occurs. Trimming
#'                  simply saves space when displaying correlation of constructs
#'                  with long names.
#' @param index     Whether to print the number of the construct. 
#'
#' @references      Bell, R. C. (2010). A note on aligning constructs. 
#'                  \emph{Personal Construct Theory & Practice}, (7), 42-48.
#'
#'                  Cohen, J. (1969). rc: A profile similarity coefficient 
#'                  invariant over variable reflection. \emph{Psychological 
#'                  Bulletin, 71}(4), 281-284.
#'
#'                  Mackay, N. (1992). Identification, Reflection, 
#'                  and Correlation: Problems In The Bases Of Repertory 
#'                  Grid Measures. \emph{International Journal of Personal
#'                  Construct Psychology, 5}(1), 57-75. 
#'
#' @return  \code{matrix} of element correlations
#'
#' @author        Mark Heckmann
#' @export
#' @seealso \code{\link{constructCor}}
#'
#' @examples
#'
#'    elementCor(mackay1992)                      # Cohen's rc
#'    elementCor(mackay1992, rc=FALSE)            # PM correlation
#'    elementCor(mackay1992, rc=FALSE, method="spearman")  # Spearman correlation
#'
#'    # format output
#'    elementCor(mackay1992, trim=6)
#'    elementCor(mackay1992, index=FALSE, trim=6)
#'
#'    # save as object for further processing
#'    r <- elementCor(mackay1992)
#'    r
#'    
#'    # change output of object
#'    print(r, digits=5)
#'    print(r, col.index=FALSE)
#'    print(r, upper=FALSE)
#'    
#'    # accessing elements of the correlation matrix
#'    r[1,3]
#'    
elementCor <- function(x, rc=TRUE, method="pearson", 
                          trim=20, index=TRUE)
{
  method <- match.arg(method,  c("pearson", "kendall", "spearman"))
  if (!inherits(x, "repgrid")) 							# check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  if (rc) {
    x <- doubleEntry(x)                     # double entry to get rc correlation
    method <- "pearson"                     # Cohen's rc is only defined for pearson correlation
  }
  scores <- getRatingLayer(x)
  res <- cor(scores, method=method)
  res <- addNamesToMatrix2(x, res, index=index, trim=trim, along=2)
  class(res) <- c("elementCor", "matrix")
  attr(res, "arguments") <- list(method=method, rc=rc)
  res
}


#' Print method for class elementCor.
#' 
#' @param x           Object of class elementCor
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param col.index   Logical. Whether to add an extra index column so the 
#'                    column names are indexes instead of construct names. This option 
#'                    renders a neater output as long construct names will stretch 
#'                    the output (default is \code{TRUE}).
#' @param upper       Whether to display upper triangle of correlation matrix only 
#'                    (default is \code{TRUE}).
#' @param ...         Not evaluated.
#' @export
#' @method            print elementCor
#' @keywords          internal
#'
print.elementCor <- function(x, digits=2, col.index=TRUE, upper=TRUE, ...)
{
  args <- attr(x, "arguments")
  class(x) <- "matrix"
  x <- round(x, digits)
  d <- format(x, nsmall=digits)
  
  # console output
  if (upper)
    d[lower.tri(d, diag=TRUE)] <- paste(rep(" ", digits + 1), collapse="", sep="")
  if (col.index)                                   # make index column for neater colnames
    d <- addIndexColumnToMatrix(d) 
  d <- as.data.frame(d)
  cat("\n############################")
  cat("\nCorrelation between elements")
  cat("\n############################")
  if (args$rc)
    args$method <- "Cohens's rc (invariant to scale reflection)"
  cat("\n\nType of correlation: ", args$method, "\n")
  if (!args$rc)
    cat("Note: Standard correlations are not invariant to scale reflection.\n")
  cat("\n")
  print(d)
}


#' Root mean square (RMS) of inter-element correlations.
#'
#' The RMS is also known as 'quadratic mean' of 
#' the inter-element correlations. The RMS serves as a simplification of the 
#' correlation table. It reflects the average relation of one element with all 
#' other elements. Note that as the correlations are squared during its calculation, 
#' the RMS is not affected by the sign of the correlation (cf. Fransella, 
#' Bell & Bannister, 2003, p. 86).
#' 
#' Note that simple element correlations as a measure of similarity
#' are flawed as they are not invariant to construct reflection (Mackay, 1992; 
#' Bell, 2010). A correlation index invariant to construct reflection is 
#' Cohen's rc measure (1969), which can be calculated using the argument 
#' \code{rc=TRUE} which is the default option in this function.
#'                  
#' @param x       \code{repgrid} object.
#' @param rc      Whether to use Cohen's rc which is invariant to construct 
#'                reflection (see desciption above). It is used as the default.
#' @param method  A character string indicating which correlation coefficient 
#'                to be computed. One of \code{"pearson"} (default), 
#'                \code{"kendall"} or \code{"spearman"}, can be abbreviated.
#'                The default is \code{"pearson"}.
#' @param trim    The number of characters an element is trimmed to (default is
#'                \code{NA}). If \code{NA} no trimming occurs. Trimming
#'                simply saves space when displaying correlation of constructs
#'                with long names.
#' @return        \code{dataframe} of the RMS of inter-element correlations
#'
#' @author        Mark Heckmann
#' @export
#' @seealso  \code{\link{constructRmsCor}}, \code{\link{elementCor}}
#'
#' @references    Fransella, F., Bell, R. C., & Bannister, D. (2003). 
#'                A Manual for Repertory 
#'                Grid Technique (2. Ed.). Chichester: John Wiley & Sons.
#'
#' @examples 
#'
#'    # data from grid manual by Fransella, Bell and Bannister
#'    elementRmsCor(fbb2003)    
#'    elementRmsCor(fbb2003, trim=10)
#'    
#'    # modify output
#'    r <- elementRmsCor(fbb2003) 
#'    print(r, digits=5)
#'
#'    # access second row of calculation results
#'    r[2, "RMS"]
#'
elementRmsCor <- function(x, rc = TRUE, method = "pearson", trim = NA)
{
  method <- match.arg(method,  c("pearson", "kendall", "spearman"))
  res <- elementCor(x, rc = rc, method = method, trim = trim,
             index = TRUE)                              # calc correlations 
  diag(res) <- NA                                       # remove diagonal 
  res <- apply(res^2, 1, mean, na.rm=TRUE)              # mean of squared values 
  res <- data.frame(RMS = res^.5)                       # root of mean squares
  class(res) <- c("rmsCor", "data.frame")
  attr(res, "type") <- "elements"
  return(res)
}


###############################################################################

#' Calculate correlations between constructs. 
#'
#' Different types of correlations can be requested: 
#' PMC, Kendall tau rank correlation, Spearman rank correlation.
#'
#' @param x         \code{repgrid} object.
#' @param method    A character string indicating which correlation coefficient 
#'                  is to be computed. One of \code{"pearson"} (default), 
#'                  \code{"kendall"} or \code{"spearman"}, can be abbreviated.
#'                  The default is \code{"pearson"}.
#' @param trim      The number of characters a construct is trimmed to (default is
#'                  \code{20}). If \code{NA} no trimming occurs. Trimming
#'                  simply saves space when displaying correlation of constructs
#'                  with long names.
#' @param index     Whether to print the number of the construct. 
#' @return          Returns a matrix of construct correlations.
#'
#' @author          Mark Heckmann
#' @export
#' @seealso \code{\link{elementCor}}
#'
#' @examples 
#'
#'    # three different types of correlations
#'    constructCor(mackay1992)                
#'    constructCor(mackay1992, method="kendall")
#'    constructCor(mackay1992, method="spearman")
#'
#'    # format output
#'    constructCor(mackay1992, trim=6)
#'    constructCor(mackay1992, index=TRUE, trim=6)
#'    
#'    # save correlation matrix for further processing
#'    r <- constructCor(mackay1992)
#'    r
#'    print(r, digits=5)
#'    
#'    # accessing the correlation matrix
#'    r[1, 3]
#'
constructCor <- function(x, method = c("pearson", "kendall", "spearman"), 
                         trim=20, index=FALSE){
  if (!inherits(x, "repgrid")) 							      # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  method <- match.arg(method)
  scores <- getRatingLayer(x)
  res <- cor(t(scores), method=method)
  res <- addNamesToMatrix2(x, res, index=index, trim=trim)
  class(res) <- c("constructCor", "matrix") 
  attr(res, "arguments") <- list(method=method)  
  return(res)
}


#' Print method for class constructCor.
#' 
#' @param x           Object of class constructCor.
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param col.index   Logical. Whether to add an extra index column so the 
#'                    column names are indexes instead of construct names. This option 
#'                    renders a neater output as long construct names will stretch 
#'                    the output (default is \code{TRUE}).
#' @param upper       Whether to display upper triangle of correlation matrix only 
#'                    (default is \code{TRUE}).
#' @param header      Whether to print additional information in header.
#' @param ...         Not evaluated.
#' @export
#' @method            print constructCor
#' @keywords          internal
#'
print.constructCor <- function(x, digits=2, col.index=TRUE,
                               upper=TRUE, header=TRUE, ...)
{
  args <- attr(x, "arguments")
  d <- x
  class(d) <- "matrix"
  d <- round(d, digits) 
  d <- format(d, nsmall=digits)
  if (upper)
    d[lower.tri(d, diag=TRUE)] <- paste(rep(" ", digits + 1), collapse="", sep="")
  if (col.index)                                   # make index column for neater colnames
    d <- addIndexColumnToMatrix(d) else
      colnames(d) <- seq_len(ncol(d))
  d <- as.data.frame(d)
  if (header) {
    cat("\n##############################")
    cat("\nCorrelation between constructs")
    cat("\n##############################")
    cat("\n\nType of correlation: ", args$method, "\n\n")
  }
  print(d)
}



#' Root mean square (RMS) of inter-construct correlations.
#'
#' The RMS is also known as 'quadratic mean' of 
#' the inter-construct correlations. The RMS serves as a simplification of the 
#' correlation table. It reflects the average relation of one construct to all 
#' other constructs. Note that as the correlations are squared during its calculation, 
#' the RMS is not affected by the sign of the correlation (cf. Fransella, 
#' Bell & Bannister, 2003, p. 86).
#'
#' @param x       \code{repgrid} object
#' @param method  A character string indicating which correlation coefficient 
#'                is to be computed. One of \code{"pearson"} (default), 
#'                \code{"kendall"} or \code{"spearman"}, can be abbreviated.
#'                The default is \code{"pearson"}.
#' @param trim    The number of characters a construct is trimmed to (default is
#'                \code{NA}). If \code{NA} no trimming occurs. Trimming
#'                simply saves space when displaying correlation of constructs
#'                with long names.
#' @return        \code{dataframe} of the RMS of inter-construct correlations
#'
#' @author        Mark Heckmann
#' @export
#' @seealso    \code{\link{elementRmsCor}}, \code{\link{constructCor}}
#'
#' @references    Fransella, F., Bell, R. C., & Bannister, D. (2003). 
#'                A Manual for Repertory 
#'                Grid Technique (2. Ed.). Chichester: John Wiley & Sons.
#'
#' @examples 
#'
#'    # data from grid manual by Fransella, Bell and Bannister
#'    constructRmsCor(fbb2003)    
#'    constructRmsCor(fbb2003, trim=20)
#'    
#'    # modify output
#'    r <- constructRmsCor(fbb2003) 
#'    print(r, digits=5)

#'    # access calculation results
#'    r[2, 1]
#'
constructRmsCor <- function(x, method = "pearson", trim = NA)
{
  method <- match.arg(method,  c("pearson", "kendall", "spearman"))
  res <- constructCor(x, method = method, trim=trim, 
                      index=TRUE)                       # calc correlations
  diag(res) <- NA                                       # remove diagonal 
  res <- apply(res^2, 1, mean, na.rm=TRUE)              # mean of squared values
  res <- data.frame(RMS = res^.5)                       # root of mean squares
  class(res) <- c("rmsCor", "data.frame")
  attr(res, "type") <- "constructs" 
  return(res)
}


#' Print method for class rmsCor (RMS correlation for constructs or elements)
#' 
#' @param x           Object of class rmsCor.
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param ...         Not evaluated.
#' @export
#' @method            print rmsCor
#' @keywords          internal
#'
print.rmsCor <- function(x, digits=2, ...)
{
  d <- as.data.frame(x)
  d <- round(d, digits)  
  type <- attr(x, "type")
  cat("\n##########################################")
  cat("\nRoot-mean-square correlation of", type)
  cat("\n##########################################\n\n")
  print(d) 
  cat("\nAverage of statistic", round(mean(unlist(d), na.rm=TRUE), digits), "\n")
  cat("Standard deviation of statistic", round(sdpop(d, na.rm=TRUE), digits), "\n")
}


#' Calculate Somers' d for the constructs. 
#'
#' Somer'ss d is an  assymetric association measure as it depends on which 
#' variable is set as dependent and independent.
#' The direction of dependency needs to be specified.
#'
#' @param x           \code{repgrid} object
#' @param dependent   A string denoting the direction of dependency in the output 
#'                    table (as d is assymetrical). Possible values are \code{"columns"}
#'                    (the default) for setting the columns as dependent, \code{"rows"} 
#'                    for setting the rows as the dependent variable and 
#'                    \code{"symmetric"} for the 
#'                    symmetrical Somers' d measure (the mean of the two directional 
#'                    values for code{"columns"} and \code{"rows"}).
#' @param trim        The number of characters a construct is trimmed to (default is
#'                    \code{30}). If \code{NA} no trimming occurs. Trimming
#'                    simply saves space when displaying correlation of constructs
#'                    with long names.
#' @param index       Whether to print the number of the construct 
#'                    (default is \code{TRUE}). 
#' @return            \code{matrix} of construct correlations.
#' @note              Thanks to Marc Schwartz for supplying the code to calculate
#'                    Somers' d.
#' @references        Somers, R. H. (1962). A New Asymmetric Measure of Association
#'                    for Ordinal Variables. \emph{American Sociological Review, 27}(6),
#'                    799-811.
#'
#' @author            Mark Heckmann
#' @export
#'
#' @examples \dontrun{
#'
#'    constructD(fbb2003)       # columns as dependent (default)
#'    constructD(fbb2003, "c")  # row as dependent
#'    constructD(fbb2003, "s")  # symmetrical index
#'  
#'    # surpress printing
#'    d <- constructD(fbb2003, out=0, trim=5)
#'    d
#'    
#'    # more digits
#'    constructD(fbb2003, dig=3)
#'
#'    # add index column, no trimming
#'    constructD(fbb2003, col.index=TRUE, index=F, trim=NA)  
#'
#' }
#'
constructD <- function(x, dependent = "columns", trim=30, index=TRUE)
{
  if (!inherits(x, "repgrid")) 							    # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  dependent <- match.arg(dependent, c("columns", "rows", "symmetric")) 
  scores <- getRatingLayer(x)
  l <- lapply(as.data.frame(t(scores)),  I)     # put each row into a list
    
  somersd <- function(x, y, dependent, smin, smax){
    na.index <- is.na(x) | is.na(y)
    x <- x[!na.index]
    y <- y[!na.index]
    x <- factor(unlist(x), levels=seq(smin, smax))
    y <- factor(unlist(y), levels=seq(smin, smax))
    m <-  as.matrix(table(x,y))
    
    if (dependent == "rows")
      i <- 1 else 
    if (dependent == "columns")
      i <- 2 else 
    if (dependent == "symmetric")
      i <- 3
    calc.Sd(m)[[i]]
  }  
  
  nc <- length(l)
  smin <- x@scale$min
  smax <- x@scale$max
  sds <- mapply(somersd, rep(l,each=nc), rep(l, nc), 
                MoreArgs=list(dependent=dependent, 
                              smin=smin, smax=smax))
  res <- matrix(sds, nc)
  res <- addNamesToMatrix2(x, res, index=index, trim=trim, along=1)
  class(res) <- c("constructD", "matrix")
  attr(res, "arguments") <- list(dependent=dependent)
  res
}


#' Print method for class constructD.
#' 
#' @param x           Object of class constructD.
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param col.index   Logical. Whether to add an extra index column so the 
#'                    column names are indexes instead of construct names. This option 
#'                    renders a neater output as long construct names will stretch 
#'                    the output (default is \code{TRUE}).
#' @param ...         Not evaluated.
#' @export
#' @method            print constructD
#' @keywords          internal
#'
print.constructD <- function(x, digits=1, col.index=TRUE, ...)
{
  class(x) <- "matrix"
  x <- round(x, digits)
  args <- attr(x, "arguments")
  attr(x, "arguments") <- NULL
  cat("\n############################")
  cat("\nSomers' D between constructs")
  cat("\n############################\n\n")
  cat("Direction:", args$dependent, "are set as dependent\n")
  if (col.index)                                  # make index column for neater colnames
    x <- addIndexColumnToMatrix(x) else
      colnames(x) <- seq_len(ncol(x))
  print(x)
}


#' Principal component analysis (PCA) of inter-construct correlations.
#'
#' Various methods for rotation and methods for the calculation of the 
#' correlations are available. Note that the number of factors
#' has to be specified. For more information on the PCA function itself type 
#' \code{?principal}. 
#'
#' @param x           \code{repgrid} object.
#' @param nfactors    Number of components to extract (default is \code{3}).
#' @param rotate      \code{"none"}, \code{"varimax"}, \code{"promax"} and \code{"cluster"} 
#'                    are possible rotations (default is \code{none}).
#' @param method      A character string indicating which correlation coefficient 
#'                    is to be computed. One of \code{"pearson"} (default), 
#'                    \code{"kendall"} or \code{"spearman"}, can be abbreviated.
#'                    The default is \code{"pearson"}.
#' @param trim        The number of characters a construct is trimmed to (default is
#'                    \code{7}). If \code{NA} no trimming occurs. Trimming
#'                    simply saves space when displaying correlation of constructs
#'                    with long names.
#' @return            Returns an object of class \code{constructPca}.
#'
#' @seealso           To extract the PCA loadings for further processing see
#'                    \code{\link{constructPcaLoadings}}.
#' @author            Mark Heckmann
#' @export
#' 
#' @references        Fransella, F., Bell, R. & Bannister, D. (2003). \emph{A Manual for Repertory 
#'                    Grid Technique} (2. Ed.). Chichester: John Wiley & Sons.
#'
#' @examples \dontrun{
#'
#'    constructPca(bell2010)
#'    
#'    # data from grid manual by Fransella et al. (2003, p. 87)
#'    # note that the construct order is different
#'    constructPca(fbb2003, nfactors=2)
#'
#'    # no rotation
#'    constructPca(fbb2003, rotate="none")
#'    
#'    # use a different type of correlation (Spearman)
#'    constructPca(fbb2003, method="spearman")
#'    
#'    # save output to object           
#'    m <- constructPca(fbb2003, nfactors=2)
#'    m
#'    
#'    # different printing options
#'    print(m, digits=5)
#'    print(m, cutoff=.3)
#'    
#' }
#'
constructPca <- function(x, nfactors=3, rotate="varimax", method = "pearson" , 
                         trim=NA) {
  
  method <- match.arg(method,c("pearson", "kendall", "spearman")) 
  rotate <- match.arg(rotate, c("none", "varimax", "promax", "cluster"))
  if (!rotate %in% c("none", "varimax", "promax", "cluster"))
    stop('only "none", "varimax", "promax" and "cluster" are possible rotations')
  
  res <- constructCor(x, method=method, trim=trim)          # calc inter constructs correations                    
  pc <- principal(res, nfactors = nfactors, rotate=rotate)  # do PCA
  class(pc) <- c("constructPca", class(pc))
  attr(pc, "arguments") <- list(nfactors=nfactors, rotate=rotate, method=method)
  return(pc)
}


#' Extract loadings from PCA of constructs.
#' 
#' @param x       \code{repgrid} object. This object is returned by the 
#'                function \code{\link{constructPca}}.
#' @return        A matrix containing the factor loadings.
#' @export
#' @examples 
#'  
#'  p <- constructPca(bell2010)
#'  l <- constructPcaLoadings(p)
#'  l[1, ]
#'  l[, 1]
#'  l[1,1]
#'  
constructPcaLoadings <- function(x)
{
  if (!inherits(x, "constructPca"))
    stop("'x' must be an object of class 'constructPca'",
         "as returned by the function 'constructPca'")
  loadings(x)  
} 


#' Print method for class constructPca.
#' 
#' @param x           Object of class constructPca.
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param cutoff      Loadings smaller than cutoff are not printed.
#' @param ...         Not evaluated.
#' @export
#' @method            print constructPca
#' @keywords          internal
#'
print.constructPca <- function(x, digits=2, cutoff=0, ...)
{
  args <- attr(x, "arguments")
  
  cat("\n#################")
  cat("\nPCA of constructs")
  cat("\n#################\n")
  
  cat("\nNumber of components extracted:", args$nfactors)
  cat("\nType of rotation:", args$rotate, "\n")
  print(loadings(x), cutoff=cutoff, digits=digits)  
}


#' Align constructs by loadings on first pricipal component. 
#'
#' In case a construct loads negatively on the first principal
#' component, the function \code{\link{alignByLoadings}} will reverse it 
#' so that all constructs have positive loadings on the first 
#' principal component (see deatil section for more).
#'
#' The direction of the constructs in a grid is arbitrary and 
#' a reflection of a scale does not affect the information 
#' contained in the grid. Nonetheless, the direction of a scale 
#' has an effect on inter-element correlations (Mackay, 1992) 
#' and on the spatial representation and clustering of the grid 
#' (Bell, 2010). Hence, it is desirable to follow a protocol to align 
#' constructs that will render unique results. A common approach 
#' is to align constructs by pole preference, but this information 
#' is not always accessible. Bell (2010) proposed another solution for
#'  the problem of construct 
#' alignment. As a unique protocol he suggests to align constructs 
#' in a way so they all have positive loadings on the first 
#' component of a grid PCA.
#'
#' @param x         \code{repgrid} object.
#' @param trim      The number of characters a construct is trimmed to (default is
#'                  \code{10}). If \code{NA} no trimming is done. Trimming
#'                  simply saves space when displaying the output.
#' @param index     Whether to print the number of the construct (e.g. for correltion 
#'                  matrices). The default is \code{TRUE}.
#' @return          An object of class \code{alignByLoadings} containing a list 
#'                  of calculations with the following entries:
#'                  
#'                   \item{cor.before}{Construct correlation matrix before reversal}
#'                   \item{loadings.before}{Loadings on PCs before reversal}
#'                   \item{reversed}{Constructs that have been reversed}
#'                   \item{cor.after}{Construct correlation matrix after reversal}
#'                   \item{loadings.after}{Loadings on PCs after reversal}
#' 
#' @note            Bell (2010) proposed a solution for the problem of construct 
#'                  alignment. As construct reversal has an effect on element 
#'                  correlation and thus on any measure that based on element 
#'                  correlation (Mackay, 1992), it is desireable to have a 
#'                  standard method for 
#'                  construct alignment independently from its semantics (preferred 
#'                  pole etc.). Bell (2010) proposes to align constructs in a way
#'                  so they all have positive loadings on the first component of a 
#'                  grid PCA.
#' 
#' @references      Bell, R. C. (2010). A note on aligning constructs.
#'                  \emph{Personal Construct Theory & Practice, 7}, 42-48.
#' 
#'                  Mackay, N. (1992). Identification, Reflection, 
#'                  and Correlation: Problems in ihe bases of repertory 
#'                  grid measures. \emph{International Journal of Personal
#'                  Construct Psychology, 5}(1), 57-75. 
#'
#' @author          Mark Heckmann
#' @export
#'
#' @seealso \code{\link{alignByIdeal}}
#'
#' @examples 
#'
#'   # reproduction of the example in the Bell (2010)
#'   # constructs aligned by loadings on PC 1
#'   bell2010                    
#'   alignByLoadings(bell2010)   
#'
#'   # save results
#'   a <- alignByLoadings(bell2010)
#'   
#'   # modify printing of resukts
#'   print(a, digits=5)
#'   
#'   # access results for further processing  
#'   names(a)
#'   a$cor.before
#'   a$loadings.before
#'   a$reversed
#'   a$cor.after
#'   a$loadings.after
#'
alignByLoadings <- function(x, trim=20, index=TRUE){
  options(warn=1)                                      # surpress warnings (TODO sometimes error in SVD due to singularities in grid)
  ccor.old <- constructCor(x, trim=trim, index=index)  # construct correlation unreversed
  pc.old <- principal(ccor.old)                        # calc principal component (psych pkg)
  reverseIndex <- 
    which(pc.old$loadings[ ,1] < 0)                    # which constructs to reverse
  x2 <- swapPoles(x, reverseIndex)                     # reverse constructs
  ccor.new <- constructCor(x2, trim=trim, index=index) # correlation with reversed constructs
  pc.new <- principal(ccor.new)                        # 2nd principal comps
  options(warn=0)                                      # reset to do warnings
  
  res <- list(cor.before=ccor.old, 
              loadings.before=pc.old$loadings[ , 1, drop=FALSE],
              reversed=data.frame(index=reverseIndex),
              cor.after=ccor.new,
              loadings.after=pc.new$loadings[ , 1, drop=FALSE])
  class(res) <- "alignByLoadings"
  res
}


#' Print method for class alignByLoadings.
#' 
#' @param x           Object of class alignByLoadings.
#' @param digits      Numeric. Number of digits to round to (default is 
#'                    \code{2}).
#' @param col.index   Logical. Whether to add an extra index column so the 
#'                    column names are indexes instead of construct names (e.g. for 
#'                    the correlation matrices). This option 
#'                    renders a neater output as long construct names will stretch 
#'                    the output (default is \code{TRUE}).
#' @param ...         Not evaluated.
#' @export
#' @method            print alignByLoadings
#' @keywords          internal
#'
print.alignByLoadings <- function(x, digits=2, col.index=TRUE, ...)
{
    cat("\n###################################")
    cat("\nAlignment of constructs by loadings")
    cat("\n###################################\n")
    
    cat("\nConstruct correlations - before alignment\n\n")
    print(x$cor.before, digits=digits, col.index=col.index, header=FALSE)
    
    cat("\nConstruct factor loadiongs on PC1 - before alignment\n\n")
    print(x$loadings.before, digits=digits)
    
    cat("\nThe following constructs are reversed:\n\n")
    if (dim(x$reversed)[1] == 0) {
      cat("None. All constructs are already aligned accordingly.\n")
    } else {
      print(x$reversed)
    }
    
    cat("\nConstruct correlations - after alignment\n\n")
    print(x$cor.after, digits=digits, col.index=col.index, header=FALSE)
      
    cat("\nConstruct factor loadings on PC1 - after alignment\n\n")
    print(x$loadings.after, digits=digits)
    cat("\n\n")
}


#' Align constructs using the ideal element to gain pole preferences.
#'
#' The direction of the constructs in a grid is arbitrary and a reflection of 
#' a scale does not affect the information contained in the grid. 
#' Nonetheless, the direction of a scale has an effect on inter-element 
#' correlations (Mackay, 1992) and on the spatial representation and clustering
#' of the grid (Bell, 2010). Hence, it is desirable to follow a protocol to 
#' align constructs that will render unique results. A common approach is
#' to align constucts by pole preference, i. e. aligninig all positive and  
#' negative poles. This can e. g. be achieved using \code{\link{swapPoles}}.
#' If an ideal element is present, this element can be used to identify
#' the positive and negative pole. The function \code{alignByIdeal} will
#' align the constructs accordingly. Note that this approach does not always 
#' yield definite results as sometimes ratings do not show a clear 
#' preference for one pole (Winter, Bell & Watson, 2010).
#' If a preference cannot be determined definitely,
#' the construct direction remains unchanged (a warning is issued in that case).
#' 
#' @param x       \code{repgrid} object
#' @param ideal   Number of the element that is used for alignment 
#'                (the ideal).
#' @param high    Logical. Whether to align the constructs so the ideal 
#'                will have high
#'                ratings on the constructs (i.e. \code{TRUE}, default) or low
#'                ratings (\code{FALSE}). High scores will lead to the preference pole
#'                on the right side, low scores will align the preference pole
#'                on the left side.
#'
#' @return        \code{repgrid} object with aligned constructs.
#'
#' @references    Bell, R. C. (2010). A note on aligning constructs.
#'                \emph{Personal Construct Theory & Practice, 7}, 42-48.
#' 
#'                Mackay, N. (1992). Identification, Reflection, 
#'                and Correlation: Problems in ihe bases of repertory 
#'                grid measures. \emph{International Journal of Personal
#'                Construct Psychology, 5}(1), 57-75. 
#'          
#'                Winter, D. A., Bell, R. C., & Watson, S. (2010). Midpoint 
#'                ratings on personal constructs: Constriction or the middle 
#'                way? \emph{Journal of Constructivist Psychology, 23}(4), 337-356.
#'
#' @export
#' @author    Mark Heckmann
#' @seealso \code{\link{alignByLoadings}}
#'
#' @examples \dontrun{
#'
#'   feixas2004                             # original grid
#'   alignByIdeal(feixas2004, 13)           # aligned with preference pole on the right
#'
#'   raeithel                               # original grid
#'   alignByIdeal(raeithel, 3, high=FALSE)  # aligned with preference pole on the left
#'
#' }
#'
alignByIdeal <- function(x, ideal, high=TRUE){
  if (!inherits(x, "repgrid")) 							  # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  	
  idealRatings <- getRatingLayer(x)[, ideal]
  unclear <- which(idealRatings == getScaleMidpoint(x))
  positive <- which(idealRatings > getScaleMidpoint(x))
  negative <- which(idealRatings < getScaleMidpoint(x))
  
  if (high)     
    x <- swapPoles(x, negative) else   # align such that ratings on ideal are high
    x <- swapPoles(x, positive)        # align such that ratings on ideal are low

  if (length(unclear) != 0)
    warning("The following constructs do not show a preference for either pole",
        "and have thus not been aligned: ", paste(unclear, collapse=","))
  x
}


#### CLUSTER ANALYIS ####

#' Cluster analysis (of constructs or elements).
#'
#' \code{cluster} is a preliminary implementation of a cluster function. 
#' It supports
#' various distance measures as well as cluster methods. More is to come.
#'  
#' \bold{align}: Aligning will reverse constructs if necessary to yield a
#' maximal similarity between constructs. In a first step the constructs are
#' clustered including both directions. In a second step the direction of a
#' construct that yields smaller distances to the adjacent constructs is
#' preserved and used for the final clustering. As a result, every construct is
#' included once but with an orientation that guarantees optimal clustering.
#' This approach is akin to the procedure used in FOCUS (Jankowicz & Thomas,
#' 1982).
#' 
#' @param x          \code{repgrid} object.
#' @param along      Along which dimension to cluster. 1 = constructs only, 
#'                   2= elements only, 0=both (default).
#' @param dmethod    The distance measure to be used. This must be one of 
#'                   "euclidean", "maximum", "manhattan", "canberra", "binary" 
#'                   or "minkowski". Any unambiguous substring can be given. 
#'                   For additional information on the different types type
#'                   \code{?dist}. 
#' @param  cmethod   The agglomeration method to be used. This should be (an
#'                   unambiguous abbreviation of) one of \code{"ward"}, 
#'                   \code{"single"}, \code{"complete"}, \code{"average"}, 
#'                   \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
#' @param  p         The power of the Minkowski distance, in case \code{"minkowski"}
#'                   is used as argument for \code{dmethod}.
#' @param align      Whether the constructs should be aligned before clustering
#'                   (default is \code{TRUE}). If not, the grid matrix is clustered 
#'                   as is. See Details section for more information.
#' @param trim       the number of characters a construct is trimmed to (default is
#'                   \code{10}). If \code{NA} no trimming is done. Trimming
#'                   simply saves space when displaying the output.
#' @param main       Title of plot. The default is a name indicating the distance 
#'                   function and cluster method.
#' @param mar        Define the plot region (bottom, left, upper, right).
#' @param cex        Size parameter for the nodes. Usually not needed.              
#' @param lab.cex    Size parameter for the constructs on the right side.
#' @param cex.main   Size parameter for the plot title (default is \code{.9}).
#' @param print      Logical. Wether to print the dendrogram (default is \code{TRUE}).
#' @param ...        Additional parameters to be passed to plotting function from
#'                   \code{as.dendrogram}. Type \code{?as.dendrogram} for further 
#'                   information. This option is usually not needed, except if special
#'                   designs are needed.
#' @return           Reordered \code{repgrid} object.
#' @references 
#' 
#' Jankowicz, D., & Thomas, L. (1982). An Algorithm for the Cluster Analysis of 
#' Repertory Grids in Human Resource Development. \emph{Personnel Review,
#' 11}(4), 15-22. doi:10.1108/eb055464.
#' 
#' @author            Mark Heckmann
#' @export
#' @seealso           \code{\link{bertinCluster}}
#'
#' @examples \dontrun{
#'
#'   cluster(bell2010)
#'   cluster(bell2010, main="My cluster analysis")   # new title
#'   cluster(bell2010, type="t")                     # different drawing style
#'   cluster(bell2010, dmethod="manhattan")          # using manhattan metric
#'   cluster(bell2010, cmethod="single")             # do single linkage clustering
#'   cluster(bell2010, cex=1, lab.cex=1)             # change appearance
#'   cluster(bell2010, lab.cex=.7,                   # advanced appearance changes
#'           edgePar = list(lty=1:2, col=2:1))
#' }
#'
cluster <- function(x, along=0, dmethod="euclidean", cmethod="ward", p=2, 
                    align=TRUE, trim=NA, main=NULL,
                    mar=c(4, 2, 3, 15), cex=0, lab.cex=.8, cex.main=.9, 
                    print=TRUE, ...)
{
  dmethods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  dmethod <- match.arg(dmethod, dmethods)
  
  cmethods <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
  cmethod <- match.arg(cmethod, cmethods)
  
  if (is.null(main))
    main <- paste(dmethod, "distance and", cmethod, "clustering")
  if (! along %in% 0:2)
    stop("'along must take the value 0, 1 or 2.")
  if (align)
    x <- align(x, along = along, dmethod = dmethod, 
               cmethod = cmethod, p = p)  
  r <- getRatingLayer(x, trim=trim)             # get ratings

  # dendrogram for constructs
  if (along %in% 0:1){
    d <- dist(r, method = dmethod, p=p)           # make distance matrix for constructs
    fit.constructs <- hclust(d, method=cmethod)   # hclust object for constructs
    dend.con <- as.dendrogram(fit.constructs)
    con.ord <- order.dendrogram(rev(dend.con))
    x <- x[con.ord, ]                             # reorder repgrid object    
  }
  
  if (along %in% c(0,2)){
    # dendrogram for elements
    d <- dist(t(r), method = dmethod, p=p)        # make distance matrix for elements
    fit.elements <- hclust(d, method=cmethod)     # hclust object for elements
    dend.el <- as.dendrogram(fit.elements)
    el.ord <- order.dendrogram(dend.el)
    x <- x[ , el.ord]                               # reorder repgrid object
  }
   
  if (print){                                       # print dendrogram?
    op <- par(mar=mar)                              # change mar settings and save old mar settings
    if (along == 0)
      layout(matrix(1:2, ncol=1))  
    if (along %in% c(0,1)){
      plot(dend.con, horiz=TRUE, main=main, xlab="",  # print cluster solution
           nodePar=list(cex=cex, lab.cex=lab.cex), 
           cex.main=cex.main, ...)
    }
    if (along %in% c(0,2)){
      plot(dend.el, horiz=TRUE, main=main, xlab="",  # print cluster solution
           nodePar=list(cex=cex, lab.cex=lab.cex), 
           cex.main=cex.main, ...)
    } 
    par(op)                                         # reset to old mar settings  
  }
  invisible(x)                                      # return reordered repgrid object
}


# function calculates cluster dendrogram from doublebind grid matrix
# and reverses the constructs accoring to the upper big cluster
align <- function(x, along = 0, dmethod = "euclidean", 
                  cmethod = "ward", p = 2, ...) 
{
  x2 <- doubleEntry(x)
  xr <- cluster(x2, dmethod=dmethod, cmethod=cmethod, p=p, 
                align=FALSE, print=FALSE)
  nc <- getNoOfConstructs(xr) / 2
  xr[1:nc, ]
}


#' Multiscale bootstrap cluster analysis.
#' 
#' p-values are calculated for each branch of the cluster dendrogram to indicate
#' the stability of a specific partition. \code{clusterBoot} will yield the same
#' clusters as the \code{\link{cluster}} function (i.e. standard hierarchical 
#' clustering) with additional p-values. Two kindes of p-values are reported: 
#' bootstrap probabilities (BP) and approximately unbiased (AU) probabilities
#' (see Details section for more information).
#' 
#' In standard (hierarchical) cluster analysis the question arises which of the 
#' idientified structures are significant or just emerged by chance. Over the
#' last decade several methods have been developed to test structures for
#' robustness. One line of research in this area is based on resampling. The
#' idea is to resample the rows or columns of the data matrix and to build the
#' dendrogram for each bootstrap sample (Felsenstein, 1985). The p-values
#' indicates the pecentage of times a specific structure is identified across
#' the bootstrap samples. It was shown that the p-value is biased (Hillis &
#' Bull, 1993; Zharkikh & Li, 1995). In the literature several methods for bias
#' correction have been proposed. In \code{clusterBoot} a method based on the
#' \emph{multiscale bootstrap} is used to derive corrected (approximately 
#' unbiased) p-values (Shimodaira, 2002, 2004). In conventional bootstrap 
#' analysis the size of the bootstrap sample is identical to the orginal sample 
#' size. Multiscale bootstrap varies the bootstrap sample size in order to infer
#' a correction formula for the biased p-value on the basis of the variation of 
#' the results for the different sample sizes (Suzuki & Shimodaira, 2006).
#'
#' \bold{align}: Aligning will reverse constructs if necessary to yield a
#' maximal similarity between constructs. In a first step the constructs are
#' clustered including both directions. In a second step the direction of a
#' construct that yields smaller distances to the adjacent constructs is
#' preserved and used for the final clustering. As a result, every construct is
#' included once but with an orientation that guarantees optimal clustering.
#' This approach is akin to the procedure used in FOCUS (Jankowicz & Thomas,
#' 1982).
#' 
#' @references 
#' 
#' Felsenstein, J. (1985). Confidence Limits on Phylogenies: An Approach Using 
#' the Bootstrap. \emph{Evolution, 39}(4), 783. doi:10.2307/2408678
#' 
#' Hillis, D. M., & Bull, J. J. (1993). An Empirical Test of Bootstrapping as a
#' Method for Assessing Confidence in Phylogenetic Analysis. \emph{Systematic Biology,
#' 42}(2), 182-192.
#' 
#' Jankowicz, D., & Thomas, L. (1982). An Algorithm for the Cluster Analysis of 
#' Repertory Grids in Human Resource Development. \emph{Personnel Review,
#' 11}(4), 15-22. doi:10.1108/eb055464.
#' 
#' Shimodaira, H. (2002) An approximately unbiased test of phylogenetic tree
#' selection. \emph{Syst, Biol., 51}, 492-508.
#' 
#' Shimodaira,H. (2004) Approximately unbiased tests of regions using multistep-
#' multiscale bootstrap resampling. \emph{Ann. Stat., 32}, 2616-2614.
#'  
#' Suzuki, R., & Shimodaira, H. (2006). Pvclust: an R package for assessing the
#' uncertainty in hierarchical clustering. \emph{Bioinformatics,
#' 22}(12), 1540-1542. doi:10.1093/bioinformatics/btl117
#' 
#' Zharkikh, A., & Li, W.-H. (1995). Estimation of confidence in phylogeny: the
#' complete-and-partial bootstrap technique. \emph{Molecular Phylogenetic Evolution,
#' 4}(1), 44-63. 
#' 
#' @param x       \code{grid object}
#' @param align   Whether the constructs should be aligned before clustering
#'                (default is \code{TRUE}). If not, the grid matrix is clustered 
#'                as is. See Details section for more information.
#' @param along   Along which dimension to cluster. 1 = constructs, 2= elements.                 
#' @inheritParams cluster
#' @param p       Power of the Minkowski metric. Not yet passed on to pvclust!
#' @inheritParams pvclust::pvclust
#' @param seed    Random seed for bootstrapping. Can be set for reproducibility (see
#'                \code{\link{set.seed}}). Usually not needed.
#' @param ...     Arguments to pass on to \code{\link{pvclust}}.
#' @return        A pvclust object as returned by the function \code{\link{pvclust}}
#'  
#' @author        Mark Heckmann
#' @export
#'
#' @examples \dontrun{
#'
#'  # p-values for construct dendrogram
#'  s <- clusterBoot(boeker)
#'  plot(s)
#'  pvrect(s, max.only=FALSE)
#'  
#'  # p-values for element dendrogram
#'  s <- clusterBoot(boeker, along=2)
#'  plot(s)
#'  pvrect(s, max.only=FALSE)
#' }
#'
clusterBoot <- function(x, along=1, align=TRUE, dmethod = "euclidean", 
                        cmethod = "ward", p=2, nboot=1000,
                        r=seq(.8, 1.4, by=.1), seed=NULL, ...)
{
  if (! along %in% 1:2)
    stop("along must either be 1 for constructs (default) or 2 for element clustering", call.=FALSE)
  if (align)
    x <- align(x)
  xr <- getRatingLayer(x)
  if (!is.null(seed) & is.numeric(seed))
    set.seed(seed)  
  if (along ==1)
    xr <- t(xr)       
  pv.e <- pvclust(xr, method.hclust=cmethod, method.dist=dmethod, r=r, nboot=nboot, ...)
  pv.e
}







### using gridn graphics
# data(mtcars)
# x  <- t(as.matrix(scale(mtcars)))
# dd.row <- as.dendrogram(hclust(dist(x)))
# row.ord <- order.dendrogram(dd.row)
# 
# dd.col <- as.dendrogram(hclust(dist(t(x))))
# col.ord <- order.dendrogram(dd.col)
# library(latticeExtra)
# dendrogramGrob(dd.col)


# data(mtcars)
# x  <- t(as.matrix(scale(mtcars)))
# dd.row <- as.dendrogram(hclust(dist(x)))
# row.ord <- order.dendrogram(dd.row)
# 
# dd.col <- as.dendrogram(hclust(dist(t(x))))
# col.ord <- order.dendrogram(dd.col)
# 
# library(lattice)
# 
# levelplot(x[row.ord, col.ord],
#           aspect = "fill",
#           scales = list(x = list(rot = 90)),
#           colorkey = list(space = "left"),
#           legend =
#           list(right =
#                list(fun = dendrogramGrob,
#                     args =
#                     list(x = dd.col, ord = col.ord,
#                          side = "right",
#                          size = 10)),
#                top =
#                list(fun = dendrogramGrob,
#                     args =
#                     list(x = dd.row, 
#                          side = "top",
#                          type = "triangle"))))





#' Normalize rows or columns by its standard deviation.  
#'
#' @param  x          \code{matrix}
#' @param  normalize  A numeric value indicating along what direction (rows, columns)
#'                    to normalize by standard deviations. \code{0 = none, 1= rows, 2 = columns}
#'                    (default is \code{0}).
#' @param ...         Not evaluated.
#' @return            Not yet definde TODO!
#'
#' @author        Mark Heckmann
#' @export
#'
#' @examples \dontrun{
#'
#'  x <- matrix(sample(1:5, 20, rep=T), 4)
#'  normalize(x, 1)						      # normalizing rows
#'  normalize(x, 2)						      # normalizing columns
#' }
#'
normalize <- function(x, normalize=0, ...){
  if (!normalize %in% 0:2 )
    stop("along must take a numeric value:\n",
         "normalize, 0 = none, 1 = rows, 2=columns")
  if (normalize == 1){
    x <- t(x) 
    x <- scale(x, center=FALSE, scale=apply(x, 2, sd, na.rm=TRUE))  # see ? scale   
    x <- t(x)
  } else  if (normalize == 2){
    x <- scale(x, center=FALSE, scale=apply(x, 2, sd, na.rm=TRUE))  # see ? scale   
  }
  x        
}


#' Centering of rows (constructs) and/or columns (elements).
#'
#' @param  x          \code{repgrid} object.
#' @param  center		  Numeric. The type of centering to be performed. 
#'                    0= no centering, 1= row mean centering (construct), 
#'                    2= column mean centering (elements), 3= double-centering (construct and element means),
#'                    4= midpoint centering of rows (constructs).
#'                    of the scale(default \code{FALSE}). Default is \code{1} (row centering).
#' @param  ...		    Not evaluated.
#' @return  \code{matrix} containing the transformed values.
#'
#' @note  If scale midpoint centering is applied no row or column centering can be 
#'        applied simultaneously.
#'        TODO: After centering the standard representation mode does not work any more as 
#'         it remains unclear what color values to attach to the centered values.
#'
#' @author        Mark Heckmann
#' @export
#'
#' @examples \dontrun{
#'
#'  center(bell2010)						      # no centering
#'  center(bell2010, rows=T)				  # row centering of grid 
#'  center(bell2010, cols=T)				  # column centering of grid
#'  center(bell2010, rows=T, cols=T)	# row and column centering
#' }
#'
center <- function(x, center=1, ...){
	dat <- x@ratings[ , ,1]
	if (center == 0)        # no centering
	  res <- dat
	if (center == 1)        # center at row (construct) mean
	  res <- sweep(dat, 1, apply(dat, 1, mean, na.rm=TRUE))
	if (center == 2)        # center at column (element) mean
    res <- sweep(dat, 2, apply(dat, 2, mean, na.rm=TRUE))
	if (center == 3){       # center at row and column mean
    res <- sweep(dat, 1, apply(dat, 1, mean, na.rm=TRUE))
    res <- sweep(res, 2, apply(res, 2, mean, na.rm=TRUE))
  }
  if (center == 4)        # center at midpoint i.e. middle of scale
	  res <- dat - getScaleMidpoint(x)  
	res
}


#' Calculate SSQ (accuracy) of biplot representation for elements 
#' and constructs. 
#'
#' Each construct and element are vectors in a 
#' multidimensional space. When reducing the representation 
#' to a lower dimensional space, a loss
#' of information (sum-of-squares) will usually occur. The output of the function
#' shows the proportion of sum-of-squares (SSQ) explained for the elements 
#' (constructs) and the amount explained by each principal component. This 
#' allows to assess which elements (construct) are represented how well in the 
#' current representation. Also it shows how much of the total variation is
#' explained.
#'
#' @param x                 \code{repgrid} object.
#' @param along             Numeric. Table of sum-of-squares (SSQ) for 1=constructs, 
#'                          2=elements (default). 
#'                          Note that currently these calculations only make sense
#'                          for biplot reperesentations with \code{g=1} and \code{h=1}
#'                          respectively.
#' @param  center		        Numeric. The type of centering to be performed. 
#'                          0= no centering, 1= row mean centering (construct), 
#'                          2= column mean centering (elements), 3= double-centering (construct and element means),
#'                          4= midpoint centering of rows (constructs).
#'                          The default is \code{1} (row centering).
#' @param normalize         A numeric value indicating along what direction (rows, columns)
#'                          to normalize by standard deviations. \code{0 = none, 1= rows, 2 = columns}
#'                          (default is \code{0}).
#' @param g                 Power of the singular value matrix assigned to the left singular 
#'                          vectors, i.e. the constructs.
#' @param h                 Power of the singular value matrix assigned to the right singular 
#'                          vectors, i.e. the elements.
#' @param col.active        Columns (elements) that are no supplementary points, i.e. they are used
#'                          in the SVD to find principal components. default is to use all elements.
#' @param col.passive       Columns (elements) that are supplementary points, i.e. they are NOT used
#'                          in the SVD but projecte into the component space afterwards. They do not 
#'                          determine the solution. Default is \code{NA}, i.e. no elements are set 
#'                          supplementary.
#' @return                  A list containing three wo elements \cr
#'                          \item{ssq.table}{dataframe with sum-of-squares explained for 
#'                                           element/construct by each dimension}
#'                          \item{ssq.table.cumsum}{dataframe with cumulated 
#'                                                  sum-of-squares explained for 
#'                                                  element/construct number of dimensions}                          
#'                          \item{ssq.total}{total sum-of-squares after pre-transforming grid matrix}
#'                          
#' @note                    TODO: if g or h is not equal to 1 the SSQ does not measure
#'                          accuracy of representation as currently the ssq of each point
#'                          are set in constrast with the pre-transformed matrix.
#'
#' @author        Mark Heckmann
#' @keywords      internal
#' @export
#' @examples
#' 
#'  # explained sum-of-squares for elements
#'  ssq(bell2010)
#'  
#'  # explained sum-of-squares for constructs
#'  ssq(bell2010, along=1)
#'  
#'  # save results
#'  s <- ssq(bell2010)
#'  
#'  # printing options
#'  print(s)
#'  print(s, digits=4)
#'  print(s, dim=3)
#'  print(s, cumulated=FALSE)
#'  
#'  # access results
#'  names(s)
#'  s$ssq.table
#'  s$ssq.table.cumsum
#'  s$ssq.total
#'  
ssq <- function(x, along=2, center=1, normalize=0, 
                g=0, h=1-g, col.active=NA, col.passive=NA, ...){ 
  x <- calcBiplotCoords(x, center=center, normalize=normalize, 
                        g=g, h=h,  col.active= col.active, 
                        col.passive=col.passive)
  if (is.null(x@calcs$biplot))
    stop("biplot coordinates have not yet been calculated")
  
  E <- x@calcs$biplot$el
  C <- x@calcs$biplot$con
  X <- x@calcs$biplot$X
  
  enames <- getElementNames(x)
  cnames <- getConstructNames(x)$rightpole
  
  ssq.c <- diag(X %*% t(X))         # SSQ for each row (constructs)
  ssq.e <- diag(t(X) %*% X)         # SSQ for each column(element)
  ssq.c.prop <- C^2 / ssq.c         # ssq of row and dimension / total ssq element  
  ssq.e.prop <- E^2 / ssq.e         # ssq of element and dimension / total ssq element  
 
  rownames(ssq.c.prop) <- cnames    # add construct names to rows
  colnames(ssq.c.prop) <-           # add dimension labels to columns
    paste(1L:ncol(ssq.c.prop), "D", sep="") 
  rownames(ssq.e.prop) <- enames    # add element names to rows
  colnames(ssq.e.prop) <-           # add dimension labels to columns
    paste(1L:ncol(ssq.e.prop), "D", sep="")

  ssq.c.cumsum <- t(apply(ssq.c.prop, 1, cumsum))         # cumulated ssq of construct per dimension / total ssq element 
  ssq.e.cumsum <- t(apply(ssq.e.prop, 1, cumsum))         # cumulated ssq of elements per dimension / total ssq element 
  ssq.c.avg <- apply(C^2, 2, sum, na.rm=T) / sum(ssq.c)   # ssq per dimension / ssq total 
  ssq.e.avg <- apply(E^2, 2, sum, na.rm=T) / sum(ssq.e)   # ssq per dimension / ssq total
  ssq.c.avg.cumsum <- cumsum(ssq.c.avg)                     # cumulated ssq per dimension / ssq total
  ssq.e.avg.cumsum <- cumsum(ssq.e.avg)                     # cumulated ssq per dimension / ssq total
  ssq.c.table <- rbind(ssq.c.prop, TOTAL=ssq.c.avg)         # 
  ssq.e.table <- rbind(ssq.e.prop, TOTAL=ssq.e.avg)         # 
  ssq.c.table.cumsum <- rbind(ssq.c.cumsum, TOTAL=ssq.c.avg.cumsum)
  ssq.e.table.cumsum <- rbind(ssq.e.cumsum, TOTAL=ssq.e.avg.cumsum)
  ssq.c.table <- ssq.c.table * 100
  ssq.e.table <- ssq.e.table * 100
  ssq.c.table.cumsum <- ssq.c.table.cumsum * 100
  ssq.e.table.cumsum <- ssq.e.table.cumsum * 100
  
  if (along == 1){          # elements or constructs?
    ssq.table <- ssq.c.table
    ssq.table.cumsum <- ssq.c.table.cumsum
    along.text <- "constructs"
  } else {
    ssq.table <- ssq.e.table
    ssq.table.cumsum <- ssq.e.table.cumsum
    along.text <- "elements"
  }
  
  l <- list(ssq.table=ssq.table, 
            ssq.table.cumsum=ssq.table.cumsum, 
            ssq.total=sum(ssq.e))
  attr(l, "arguments") <- list(along.text=along.text)
  class(l) <- "ssq"
  return(l) 
}


#' Print method for class ssq.
#' 
#' @param x                 Object of class ssq.
#' @param digits            Number of digits to round the output to (default is \code{2}).
#' @param dim               The number of PCA dimensions to print. Default
#'                          is \code{5} dimensions. \code{NA} will print all 
#'                          dimensions.
#' @param cumulated         Logical (default is \code{TRUE}). 
#'                          Print a cumulated table of sum-of-squares?
#'                          If \code{FALSE} the uncumulated sum-of-squares are printed.                         
#'                          (default is \code{TRUE}).
#' @param ...               Not evaluated.
#' @export
#' @method            print ssq
#' @keywords          internal
#'
#'
print.ssq <- function(x, digits=2, dim=5, cumulated=TRUE, ...)
{
  # dimensions to print
  if (is.na(dim[1]))
    dim <- ncol(x$ssq.table)
  args <- attr(x, "arguments")
  if (cumulated) {                 # output cumulated table?
    ssq.out <- x$ssq.table.cumsum 
    cum.text <- "Cumulated proportion"
  } else {
    ssq.out <- x$ssq.table    
    cum.text <- "Proportion"
  }
  cat("\n", cum.text, "of explained sum-of-squares for ", 
      args$along.text, "\n\n", sep="")
  print(round(ssq.out[, 1:dim], digits))
  cat("\nTotal sum-of-squares of pre-transformed ", 
      "(i.e. centered and scaled) matrix:", x$ssq.total) 
}


