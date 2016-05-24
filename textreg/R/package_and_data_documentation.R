#' Sparse regression package for text that allows for multiple word phrases.
#'
#' Built on Georgiana Ifrim's work, but allowing for regularization of phrases, this 
#' package does sparse regression using greedy coordinate descent.
#' In a nutshell, the textreg package allows for regressing a vector of +1/-1 labels onto raw text. 
#' The textreg package takes care of converting the text to all of the possible related features, allowing
#' you to think of the more organic statement of regressing onto ``text'' in some broad sense.
#' 
#' Implementation-wise, it is a wrapper for a modified version of the C++ code written by Georgiana Ifrim to do this regression.
#' It is also designed to (somewhat) integrate with the tm package, a commonly used R package for dealing with text.
#'
#' One warning: this package uses tm, but does need to generate vectors of character strings to pass to the textreg call, which can be quite expensive.  
#' You can also pass a filename to the textreg call instead, which allows one to avoid loading a large corpus into memory and then copying it over.
#' You can use a prior build.corpus command before textreg to mitigate this cost, but it is an imperfect method.
#'
#' The n-gram package is documented, but it is research code, meaning gaps and errors are possible; the author would appreciate notification of anything that is out of order.
#'
#' The primary method in this package is the regression call textreg(). This method takes a corpus and a labeling vector and returns a textreg.result
#' object that contains the final regression result along with diagnostic information that can be of use.  
#'
#' Start by reading the ``bathtub'' vignette, which walks through most of the functionality of this package.
#'
#' Special thanks and acknowledgements to Pavel Logacev, who found some subtle bugs on the windows platform and gave excellent advice in general.
#' Also thanks to Kevin Wu, who wrote earlier versions of the stemming and cross-validation code.  And Georgiana Ifrim, of course, for the earlier
#' version of the C++ code.
#'
#' @useDynLib textreg
#' @references Ifrim, G., Bakir, G., & Weikum, G. (2008). Fast logistic regression for text categorization with variable-length n-grams. 14th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 354-362.
#' @references Ifrim, G., & Wiuf, C. (2011). Bounded coordinate-descent for biological sequence classification in high dimensional predictor space. 17th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 708-716.
#' @references Jia, J., Miratrix, L., Yu, B., Gawalt, B., Ghaoui, El, L., Barnesmoore, L., & Clavier, S. (2014). Concise Comparative Summaries (CCS) of Large Text Corpora with a Human Experiment. The Annals of Applied Statistics, 8(1), 499-529.
#' @references Miratrix, L., & Ackerman, R. (2014). A method for conducting text-based sparse feature selection for interpretability of selected phrases.
#' @import Rcpp
#' @import tm
#' @docType package
#' @name textreg-package
NULL






#' Some small, fake test corpora.
#' 
#' A list of several fake documents along with some labeling schemes primarily used by the unit testing code.
#' Also used in some examples.
#'
#' @docType data
#' @keywords datasets
#' @format A list of dataframes
#' @name testCorpora
NULL



#' Sample of cleaned OSHA accident summaries.
#' 
#' bathtub consists of several accident reports plus a labeling with a +1 for any report
#' that had been tagged as related to METHELYNE CHLORIDE.
#'
#' @docType data
#' @keywords datasets
#' @import tm
#' @format Corpus object from the \code{tm} package.  Has a meta info of the METHELYNE CHLORIDE labeling called "meth.chl"
#' @name bathtub
#' @family bathtub
#' @examples
#' library( tm )
#' data( bathtub )
#' meta( bathtub, "meth.chl" )
NULL




#' Sample of raw-text OSHA accident summaries.
#' 
#' dirtyBathtub consists of the (more) raw data from which the \code{bathtub} dataset is derived.
#'
#' @docType data
#' @keywords datasets
#' @format Dataframe.  Has a meta info of the METHELYNE CHLORIDE labeling, plus 100s of other labels.
#' @name dirtyBathtub
#' @family bathtub
#' @examples
#' data( dirtyBathtub )
#' table( dirtyBathtub$fatality )
NULL
