#' Dot plots for influence diagnostics
#'
#' This is a function that can be used to create (modified) dotplots for the
#' diagnostic measures.  The plot allows the user to understand the distribution
#' of the diagnostic measure and visually identify unusual cases.
#' 
#' @note
#' The resulting plot uses \code{coord_flip} to rotate the plot, so when
#' adding customized axis labels you will need to flip the usage of 
#' \code{xlab} and \code{ylab}.
#'
#' @param x values of the diagnostic of interest
#' @param index index (IDs) of \code{x} values
#' @param data data frame to use (optional)
#' @param cutoff value(s) specifying the boundary for unusual values of the 
#' diagnostic. The cutoff(s) can either be supplied by the user, or automatically
#' calculated using measures of internal scaling if \code{cutoff = "internal"}
#' @param name what diagnostic is being plotted 
#' (one of \code{"cooks.distance"}, \code{"mdffits"}, \code{"covratio"}, 
#' \code{"covtrace"}, \code{"rvc"}, or \code{"leverage"}).
#' this is used for the calculation of "internal" cutoffs
#' @param modify specifies the \code{geom} to be used to produce a 
#' space-saving modification: either \code{"dotplot"} or \code{"boxplot"}
#' @param ... other arguments to be passed to \code{qplot()}
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @examples 
#' data(sleepstudy, package = 'lme4')
#' fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'
#' # Subject level deletion and diagnostics
#' subject.del  <- case_delete(model = fm, group = "Subject", type = "both")
#' subject.diag <- diagnostics(subject.del)
#' 
#' dotplot_diag(x = COOKSD, index = IDS, data = subject.diag[["fixef_diag"]], 
#'              name = "cooks.distance", modify = FALSE, 
#'              xlab = "Subject", ylab = "Cook's Distance")
#' 
#' dotplot_diag(x = sigma2, index = IDS, data = subject.diag[["varcomp_diag"]], 
#'              name = "rvc", modify = "dotplot", cutoff = "internal", 
#'              xlab = "Subject", ylab = "Relative Variance Change")
#'              
#' dotplot_diag(x = sigma2, index = IDS, data = subject.diag[["varcomp_diag"]], 
#'              name = "rvc", modify = "boxplot", cutoff = "internal", 
#'              xlab = "Subject", ylab = "Relative Variance Change")
#' @export
#' @keywords hplot
dotplot_diag <- function(x, index, data, cutoff, 
                         name = c("cooks.distance", "mdffits", "covratio",
                                  "covtrace", "rvc", "leverage"),
                         modify = FALSE, ... ){

#   if(!is.logical(modify)) stop("modify should be either TRUE or FALSE")
  
  if(!modify %in% c(FALSE, "boxplot", "dotplot")) {
    stop("modify should be FALSE or either 'boxplot' or 'dotplot'")
  }
  
  if(modify != FALSE & missing(cutoff)){
    stop("a cutoff should be specified if a modified dotplot is requested")
  }
  
  if(modify != FALSE & missing(name)){
    stop("a name should be specified if a modified dotplot is requested")
  }
  
  if(!missing(cutoff)){
    if( !is.numeric(cutoff) && cutoff != "internal" ){
      stop("cutoff should be numeric or 'internal'")
    }
  } else{
    cutoff <- NULL
  }
  
  if(!missing(name)){
    if(!name %in% c("cooks.distance", "mdffits", "covratio",
                    "covtrace", "rvc", "leverage")) {
      stop("name should be one of 'cooks.distance', 'mdffits', 'covratio', 
         'covtrace', 'rvc', 'leverage'")
    }
  }
  
  if(missing(data)) {
    data <- data.frame() 
    if(missing(index)) index <- factor(seq(1, length(x)))
  } else{
    if(missing(index)) index <- factor(seq(1, nrow(data)))
  }
  
  name <- match.arg(name)
  x <- eval(substitute(x), data, parent.frame())
  index <- eval(substitute(index), data, parent.frame())
  
  if(class(x) %in% c("fixef.dd", "vcov.dd")) x <- as.numeric(x)
  
  if(modify == FALSE) {
    p <- qplot(x = reorder(index, x, identity), y = x, geom = "blank", ... )
  }
  
  if(!is.null(cutoff)){
    
    if(!is.numeric(cutoff)) cutoff <- internal_cutoff(x = x, name = name)

    if(is.numeric(cutoff)){
      if( !name %in% c("covratio", "rvc") ){
        
        extreme <-  x > cutoff

        if(modify != FALSE){
          levels(index)[levels(index) %in% index[which(extreme == FALSE)]] <- "within cutoff"
          p <- qplot(x = reorder(index, x, mean), y = x, geom = "blank", ... )
          
          if(modify == "boxplot"){
            p <- p + geom_boxplot(aes(x = index[which(extreme == FALSE)],
                                      y = x[which(extreme == FALSE)]), 
                                  inherit.aes = FALSE)
          }
          if(modify == "dotplot"){
            p <- p  + geom_point(aes(x = index[which(extreme == FALSE)],
                                     y = x[which(extreme == FALSE)]), 
                                 colour = I("blue"), inherit.aes = FALSE)
          }
        } else {
          p <- p  + geom_point(aes(x = index[which(extreme == FALSE)],
                                   y = x[which(extreme == FALSE)]), 
                               colour = I("blue"), inherit.aes = FALSE)
        }

        if( sum(extreme) > 0 ){
          p <- p +
            geom_point(aes( x = index[which(extreme == TRUE)],
                            y = x[which(extreme == TRUE)]),
                       colour = I("red"), shape = 17, inherit.aes = FALSE) +
            geom_text(aes(x = index[which(extreme == TRUE)],
                          y = x[which(extreme == TRUE)],
                          label = index[which(extreme == TRUE)], 
                          hjust=.5, vjust=1.5, size=3),
                      inherit.aes = FALSE)
        }


#         	p + geom_point(aes(x = index[which(extreme == FALSE)],
#         	               y = x[which(extreme == FALSE)]), 
#                          colour = I("blue"), inherit.aes = FALSE) + 
#           		geom_hline(aes(yintercept = cutoff), colour=I("red")) +
#               	theme(legend.position = "none") +
#                 coord_flip()
		  p + geom_hline(aes(yintercept = cutoff), colour=I("red")) +
		    theme(legend.position = "none") +
		    coord_flip()
    }

      else{
        extreme <- x < cutoff[1] | x > cutoff[2]

        if(modify != FALSE){
          levels(index)[levels(index) %in% index[which(extreme == FALSE)]] <- "within cutoff"
          p <- qplot(x = reorder(index, x, mean), y = x, geom = "blank", ... )
          
          if(modify == "boxplot"){
            p <- p + geom_boxplot(aes(x = index[which(extreme == FALSE)],
                                      y = x[which(extreme == FALSE)]), 
                                  inherit.aes = FALSE)
          }
          if(modify == "dotplot"){
            p <- p  + geom_point(aes(x = index[which(extreme == FALSE)],
                                     y = x[which(extreme == FALSE)]), 
                                 colour = I("blue"), inherit.aes = FALSE)
          }
        } else {
          p <- p  + geom_point(aes(x = index[which(extreme == FALSE)],
                                   y = x[which(extreme == FALSE)]), 
                               colour = I("blue"), inherit.aes = FALSE)
        }
        
          
          if(sum(extreme) > 0){
            p <- p + 
              geom_point(aes(x = index[which(extreme == TRUE)],
                             y = x[which(extreme == TRUE)]), 
                         colour = I("red"), shape = 17, inherit.aes = FALSE) +
              geom_text(aes(x = index[which(extreme == TRUE)],
                            y = x[which(extreme == TRUE)], 
                            label = index[which(extreme == TRUE)], 
                            hjust=.5, vjust=1.5, size=3),
                        inherit.aes = FALSE)
          } 
        

        	p + 
#             geom_point(aes(x = index[which(extreme == FALSE)],
#                              y = x[which(extreme == FALSE)]), 
#                          colour = I("blue"), inherit.aes = FALSE) + 
          		geom_hline(aes(yintercept = cutoff[1]), colour=I("red")) + 
            	geom_hline(aes(yintercept = cutoff[2]), colour=I("red")) + 
                theme(legend.position = "none") +
                coord_flip()	
        
      }
        
      }
    }

  else{
    p + geom_point() + coord_flip()
  }
}



# Calculating a cutoff value for diagnostic measures
#
# This function provides cutoff values using internal scaling. 
# In other words, a measure of relative standing is used (3*IQR)
# to specify unusual values of the diagnostic of interest relative
# to the vector of diagnostics.
#
# @param x a vector
# @param name specification of which diagnostic to plot
internal_cutoff <- function(x, name){
  q3 <- quantile(x, p=0.75)
  x.iqr <- IQR(x)
	
  if(name %in% c("covratio", "rvc")){
    q1 <- quantile(x, p=0.25)
    cutoff <- c(lower = q1 - 3 * x.iqr, upper = q3 + 3 * x.iqr)	
  }
	
  else{cutoff <- q3 + 3 * x.iqr}
	
  return(cutoff)	
}
