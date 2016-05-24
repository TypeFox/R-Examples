#' Spectral Map Analysis
#' Produces an object of class \code{mpm} that allows for exploratory
#' multivariate analysis of large data matrices, such as gene expression data
#' from microarray experiments.
#' 
#' The function \code{mpm} presents a unified approach to exploratory
#' multivariate analysis encompassing principal component analysis,
#' correspondence factor analysis, and spectral map analysis. The algorithm
#' computes projections of high dimensional data in an orthogonal space. The
#' resulting object can subsequently be used in the construction of biplots
#' (i.e. \code{plot.mpm}).
#' 
#' The projection of the pre-processed data matrix in the orthogonal space is
#' calculated using the \code{La.svd} function.
#' 
#' @param data a data frame with the row descriptors in the first column. For
#'   microarray data rows indicate genes and columns biological samples.
#' @param logtrans an optional logical value. If \code{TRUE}, data are first
#'   transformed to logarithms (base e) before the other operations.
#'   Non-positive numbers are replaced by \code{logrepl}. If \code{FALSE}, data
#'   are left unchanged. Defaults to \code{TRUE}.
#' @param logrepl an optional numeric value that replaces non-positive numbers
#'   in log-transformations. Defaults to \code{1e-9}.
#' @param closure optional character string specifying the closure operation
#'   that is carried out on the optionally log-transformed data matrix. If
#'   \kbd{"double"}, data are divided by row- and column-totals. If \kbd{"row"}
#'   data are divided by row-totals. If \kbd{"column"} data are divided by
#'   column-totals. If \kbd{"none"} no closure is carried out. Defaults to
#'   \kbd{"none"}.
#' @param center optional character string specifying the centering operation
#'   that is carried out on the optionally log-transformed, closed data matrix.
#'   If \kbd{"double"} both row- and column-means are subtracted. If
#'   \kbd{"row"} row-means are subtracted. If \kbd{"column"} column-means are
#'   subtracted. If \kbd{"none"} the data are left uncentered. Defaults to
#'   \kbd{"double"}.
#' @param normal optional character string specifying the normalization
#'   operation that is carried out on the optionally log-transformed, closed,
#'   and centered data matrix. If \kbd{"global"} the data are normalized using
#'   the global standard deviation. If \kbd{"row"} data are divided by the
#'   standard deviations of the respective row. If \kbd{"column"} data are
#'   divided by their respective column standard deviation. If \kbd{"none"} no
#'   normalization is carried out. Defaults to \kbd{"global"}.
#' @param row.weight optional character string specifying the weights of the
#'   different rows in the analysis. This can be \kbd{"constant"},
#'   \kbd{"mean"}, \kbd{"median"}, \kbd{"max"}, \kbd{"logmean"}, or \kbd{"RW"}.
#'   If \kbd{"RW"} is specified, weights must be supplied in the vector
#'   \kbd{RW}. In other cases weights are computed from the data. Defaults to
#'   \kbd{"constant"}, i.e. constant weighting.
#' @param col.weight optional character string specifying the weights of the
#'   different columns in the analysis. This can be \kbd{"constant"},
#'   \kbd{"mean"}, \kbd{"median"}, \kbd{"max"}, \kbd{"logmean"}, or \kbd{"CW"}.
#'   If \kbd{"CW"} is specified, weights must be supplied in the vector
#'   \code{CW}. In other cases weights are computed from the data. Defaults to
#'   \kbd{"constant"}, i.e. constant weighting.
#' @param CW optional numeric vector with external column weights. Defaults to
#'   1 (constant weights).
#' @param RW optional numeric vector with external row weights. Defaults to 1
#'   (constant weights).
#' @param pos.row logical vector indicating rows that are not to be included in
#'   the analysis but must be positioned on the projection obtained with the
#'   remaining rows. Defaults to \code{FALSE}.
#' @param pos.column logical vector indicating columns that are not to be
#'   included in the analysis but must be positioned on the projection obtained
#'   with the remaining columns. Defaults to \code{FALSE}.
#' @return An object of class \code{mpm} representing the projection of data
#'   after the different operations of transformation, closure, centering, and
#'   normalization in an orthogonal space. Generic functions \code{plot} and
#'   \code{summary} have methods to show the results of the analysis in more
#'   detail. The object consists of the following components:
#'   \item{TData}{matrix with the data after optional log-transformation,
#'   closure, centering and normalization.} \item{row.names}{character vector
#'   with names of the row elements as supplied in the first column of the
#'   original data matrix} \item{col.names}{character vector with the names of
#'   columns obtained from the column names from the original data matrix}
#'   \item{closure}{closure operation as specified in the function call}
#'   \item{center}{centering operation as specified in the function call}
#'   \item{normal}{normalization operation as specified in the function call}
#'   \item{row.weight}{type of weighting used for rows as specified in the
#'   function call} \item{col.weight}{type of weighting used for columns as
#'   specified in the function call} \item{Wn}{vector with calculated weights
#'   for rows} \item{Wp}{vector with calculated weights for columns}
#'   \item{RM}{vector with row means of original data} \item{CM}{vector with
#'   column means of original data} \item{pos.row}{logical vector indicating
#'   positioned rows as specified in the function call}
#'   \item{pos.column}{logical vector indicating positioned columns as
#'   specified in the function call} \item{SVD}{list with components returned
#'   by \code{La.svd}} \item{eigen}{eigenvalues for each orthogonal factor from
#'   obtained from the weighted singular value decomposition}
#'   \item{contrib}{contributions of each factor to the total variance of the
#'   pre-processed data, i.e. the eigenvalues as a fraction of the total
#'   eigenvalue.} \item{call}{the matched call.}
#' @note Principal component analysis is defined as the projection onto an
#'   orthogonal space of the column-centered and column-normalized data. In
#'   correspondence factor analysis the data are pre-processed by double
#'   closure, double centering, and global normalization. Orthogonal projection
#'   is carried out using the weighted singular value decomposition. Spectral
#'   map analysis is in essence a principal component analysis on the
#'   log-transformed, double centered and global normalized data. Weighted
#'   spectral map analysis has been proven to be successful in the detection of
#'   patterns in gene expression data (Wouters et al., 2003).
#' @author Luc Wouters, Rudi Verbeeck, Tobias Verbeke
#' @seealso \code{\link{plot.mpm}}, \code{\link{summary.mpm}}
#' @references Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U.,
#'   Molenberghs, G., Lewi, P.J. (2003). Graphical exploration of gene
#'   expression data: a comparative study of three multivariate methods.
#'   \emph{Biometrics} \bold{59}, 1131-1140.
#' @keywords multivariate
#' @examples
#' 
#'   data(Golub)
#'   # Principal component analysis
#'   r.pca <- mpm(Golub[,1:39], center = "column", normal = "column")
#'   # Correspondence factor analysis
#'   r.cfa <- mpm(Golub[,1:39],logtrans = FALSE, row.weight = "mean",
#'              col.weight = "mean", closure = "double")
#'   # Weighted spectral map analysis
#'   r.sma <- mpm(Golub[,1:39], row.weight = "mean", col.weight = "mean")
#'
#' @export
mpm <- function(data, 
	logtrans = TRUE,    
  logrepl = 1e-9,     
	center = c("double", "row", "column", "global", "none"),
	normal = c("global", "row", "column", "none"),
	closure = c("none", "row", "column", "global", "double"),
	row.weight = c("constant", "mean", "median", "max", "logmean", "RW"),
	col.weight = c("constant", "mean", "median", "max", "logmean", "CW"),
	CW = rep(1, ncol(data)-1),  
	RW = rep(1, nrow(data)),     
	pos.row = rep(FALSE, nrow(data)),            # Positioned rows and columns are not taken into account during calculation,
	pos.column = rep(FALSE, ncol(data) - 1)){    # but are still plotted at the correct location
  
  
  ### Error checking and argument matching
  if (missing(data))
      stop("Argument \"data\" is missing, with no default")
  NData <- as.matrix(data[, -1]) # drop column 1 with row names
  if (any(is.na(NData)) || !is.numeric(NData))
      stop("Data must be numeric without NA's")
  if (length(pos.row) != dim(NData)[1])
      stop("Length of pos.row argument not equal to number of rows in table")
  if (length(pos.column) != dim(NData)[2])
      stop("Length of pos.column argument not equal to number of columns in table")
  if (length(RW) != nrow(NData))
      stop("Length of RW not equal to number of rows in table")
  if (length(CW) != ncol(NData))
      stop ("Length of CW not equal to number of columns in table")
    
    
  # Parse other arguments, store in variables
  # If no value is specified, the first option in the list is taken as the default value  
  center <- match.arg(center)
  normal <- match.arg(normal)
  closure <- match.arg(closure)
  row.weight <- match.arg(row.weight)
  col.weight <- match.arg(col.weight)
  logrepl <- max(1e-9, logrepl)
  
  ### Get row-descriptors as character variable from column 1
  Row.Names <- as.character(data[, 1]) # TV: simplified
  
  ###################
  ### Positioning ###
  ################### 
  
  #### Determine positioned rows-columns, set to NA 
  RData <- NData
  RData[pos.row, ] <- NA
  RData[, pos.column] <- NA
  
  ####################
  ### Reexpression ###
  ####################
  
  ### Logarithmic transform with replacement of non-positive numbers
  
  ## Define logtransform function
  logtransf <- function(x, logrepl, comp)
  {
    if (any(x <= 0, na.rm = TRUE)) # Test if there is non-positive data
    # (RV: added na.rm=T: test fails when positioning is used)
    {
      warning(paste("Non-positive data replaced by",
              logrepl, "in computing", comp, "in: spectralmap.\n"),
          call. = FALSE)
      # Replace non-positive data by logrepl
      x[x <= 0] <- logrepl
    }
    return(log(x))
  }

  LData <- if (logtrans) logtransf(NData, logrepl, comp = "logarithms") else NData # just copy the data from the previous step
  
  ### means of original data matrix
  RM <- rowMeans(NData[, !pos.column]) ### transformation does not 
                                       ### have impact on this !!
  CM <- colMeans(NData[!pos.row, ])   
  
  ### define weights
  Wn <- pmax(0, switch(row.weight,
      constant = rep(1, length(RM)),
  	  mean = apply(RData, 1, mean, na.rm = TRUE),
  	  median = apply(RData, 1, median, na.rm = TRUE),
  	  max = apply(RData, 1, max,na.rm = TRUE),
  	  logmean = apply(logtransf(RData, logrepl, comp = "logmean weights"),
                    1, mean, na.rm = TRUE),
  	  RW = RW))

  Wp <- pmax(0,switch(col.weight,
  	  constant = rep(1, length(CM)),
  	  mean = apply(RData, 2, mean, na.rm = TRUE),
  	  median = apply(RData, 2, median, na.rm = TRUE),
  	  max = apply(RData,2,max,na.rm=TRUE),
  	  logmean = apply(logtransf(RData, logrepl, comp = "logmean weights"),
                      2, mean, na.rm = TRUE) ,
  	  CW = CW))

  Wn[pos.row] <- 0
  Wp[pos.column] <- 0
  
  Wn <- Wn / sum(Wn)   # normalize weights to unit sum
  Wp <- Wp / sum(Wp)   
  
  ###############
  ### Closure ###
  ###############
  
  if (closure != "none" && any(LData < 0)) 
    warning("Closure operation with non-positive data")
  
  Tn <- rowSums(LData[, !pos.column], na.rm = TRUE) # row totals
  Tp <- colSums(LData[!pos.row, ], na.rm = TRUE)    # column totals
  # Positioned rows (in Tn) / columns (in Tp) are not excluded to preserve the dimensions
  # Consequently, Tn <> Tp <> Tt
  Tt <- sum(LData[!pos.row, !pos.column], na.rm = TRUE)   # global total
  
  ClData <- switch(closure,
      none = LData,
      row = sweep(LData, 1, Tn, "/"),
      column = sweep(LData, 2, Tp, "/"),
      global = LData / Tt,
      double = Tt * sweep(sweep(LData, 1, Tn, "/"), 2, Tp, "/"))
  
  if (any(!is.finite(ClData))) # TV: a bit indirect ?, changed to is.finite 
    stop("Division by 0 in closure operation")
    
  #################
  ### Centering ### 
  #################
  
  ## using weighted means (note: positioned data have weight 0)
  
  Mp <- colSums(sweep(ClData, 1, Wn, "*"))  # weighted row means. Use "sum" as Wn is normalized.
  Mn <- rowSums(sweep(ClData, 2, Wp, "*"))  # weighted column means
  Mg <- sum(Mp * Wp) # weighted global mean (weigh each row by Wn, each column by Wp and add everything together)

  CData <- switch(center,
  	double = Mg + sweep(sweep(ClData, 2, Mp), 1, Mn),
  	row = sweep(ClData, 1, Mn),
  	column = sweep(ClData, 2, Mp), # RV: originally the code swept across LData instead of ClData
  	global = ClData - Mg,
  	none = ClData)

  #######################
  ### Standardization ###
  #######################
  
  Vp <- colSums(sweep(CData^2, 1, Wn, "*"))
  Vn <- rowSums(sweep(CData^2, 2, Wp, "*"))
  Vg <- sum(Vp * Wp)
  
  SData <- switch(normal,
  	  global = CData / sqrt(Vg),
  	  row = sweep(CData, 1, sqrt(Vn), "/"),
  	  column = sweep(CData, 2, sqrt(Vp), "/"), # RV: originally the code swept across LData instead of ClData
  	  none = CData)
  
  #####################
  ### Factorization ###
  #####################
  
  WData <- sweep(sweep(SData, 1, sqrt(Wn), "*"), 2, sqrt(Wp), "*") #  weighted data matrix
  svd.res <- La.svd(WData)       # Singular Value Decomposition
  eigen <- svd.res$d^2           # Eigenvalues
  contrib <- eigen / sum(eigen)  # Contributions (normalized singular values)
  
  ### return
  r <- list(TData = SData,
          	row.names = Row.Names,
          	col.names = names(data)[-1], # drop name of column one (row names)
          	closure = closure,
          	center = center,
          	normal = normal,
          	row.weight = row.weight,
          	col.weight = col.weight,
          	eigen = eigen,
          	contrib = contrib,
          	Rm = RM,
          	Cm = CM,
          	Wn = Wn,
          	Wp = Wp,
          	SVD = svd.res,
          	pos.column = pos.column,
          	pos.row = pos.row,
          	call = match.call())
  class(r) <- "mpm"
  return(r)
}
