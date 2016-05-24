#' Model all possible end-member scenarios.
#' 
#' This function takes a definition of weight transformation 
#' limits and corresponding minimum and maximum numbers of end-members to 
#' model all end-member scenarios in accordance with these parameters. Based 
#' on the output the user can decide on robust end-members.
#' 
#' The plot output is an overlay of several data. The coloured lines in the 
#' background are end-member loadings (number noted in the plot title), 
#' resulting from all possible model scenarios. If \code{col.q == TRUE} they
#' are coloured according to the number of end-members with which the model 
#' was generated. This colour scheme allows to depict end-members that emerge
#' for model realisations with specific number of end-members. The thick 
#' black line is a kernel density estimate curve, generated from the mode 
#' positions of all end-members. The kernel bandwidth is set to 1 percent of 
#' the number of grain-size classes of the input data set, which gave useful
#' results for most of our test data sets. The cumulaitve dot-line-plot is a
#' further visualisation of end-member mode positions. The function is a 
#' modified wrapper function for the function \code{test.robustness()}.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param q Numeric matrix, definitions of minimum and maximum number of
#' end-members (cf. \code{get.q()}), required.
#' @param l Numeric vector, weight transformation limit values, corresponding
#' to the matrix q, required.
#' @param plot Logical scalar, option to plot the results (cf. details for 
#' explanations), default is \code{TRUE}.
#' @param col.q Logical scalar, option to colour end-member loadings by the 
#' number of end-members which were used to create the model realisation,
#' default is \code{TRUE}.
#' @param bw Numeric scalar, optional manual setting of the kde bandwidth. 
#' By default, bw is calculated as 1 percent of the number of grain-size 
#' classes.
#' @param \dots Further arguments passed to the function.
#' @return \code{List} object with all modelled end-members, each described by
#' input parameters, mode position, quality measures and value distributions.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{test.l.max}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X, envir = environment())
#' 
#' ## define input parameters
#' l <- c(0, 0.05, 0.10)
#' q <- cbind(c(2, 2, 3), c(5, 6, 4))
#' 
#' ## infer l-vector
#' em_pot <- model.em(X = X, q = q, l = l)
#' 
#' @export model.em
model.em <- function(
  X,
  q,
  l,
  plot = TRUE,
  col.q = TRUE,
  bw,
  ...
) {
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## check input data
  if(length(l) != nrow(q)) {
    stop("l and q are not of identical length!")
  }
  
  ## check/set bw
  if(missing(bw) == TRUE) {
    bw <- ncol(X) / 100
  }
  
  ## create P-matrix
  P <- cbind(q, l)
  
  ## run test.robustness()
  em <- test.robustness(X = X, P = P, ...)
  
  ## assign plot colour
  if(col.q == TRUE) {
    plot_col <- em$q - min(em$q) + 1
  } else {
    plot_col <- rep(x = "grey70",
                    times = nrow(em$loadings))
  }
  
  ## optionally, plot output
  if(plot == TRUE) {
    
    ## create empty, scaled plot
    plot(NA, 
         xlim = c(1,ncol(em$loadings)), 
         ylim = c(min(em$loadings), max(em$loadings) * 1.1),
         xlab = "Class",
         ylab = "Contribution",
         main = paste("Loadings (n = ", nrow(em$loadings), ")", sep = ""))
    
    ## add all potential loadings
    for(i in 1:nrow(em$loadings)) {
      lines(x = 1:ncol(em$loadings),
            y = em$loadings[i,],
            col = adjustcolor(col = plot_col[i], alpha.f = 0.3))
    }
    
    ## add cumulate mode position plot
    par(new = TRUE)
    plot(x = sort(em$modes),
         y = 1:length(em$modes),
         xlim = c(1,ncol(em$loadings)),
         type = "b",
         lwd = 2,
         col = plot_col,
         ann = FALSE,
         axes = FALSE)
    
    ## add mode position KDE plot
    par(new = TRUE)
    kde <- density(x = em$modes,
                   from = 1, 
                   to = ncol(X),
                   bw = bw)
    plot(kde,
         ylim = c(0, max(kde$y) * 1.2),
         lwd = 2,
         ann = FALSE,
         axes = FALSE)
  }
  
  if(col.q == TRUE) {
    legend(x = "top",
           legend = paste("q = ", seq(from = min(em$q), 
                                      to = max(em$q))),
           text.col = sort(x = unique(x = plot_col)),
           horiz = TRUE,
           box.lty = 0)
  }

  ## return output
  return(em)
}