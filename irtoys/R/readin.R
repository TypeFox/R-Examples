
#' Read in parameter estimates
#' 
#' From BILOG output, read in estimates of item parameters. Invoked
#' automatically when the model is estimated via \code{irtoys}. If that is not
#' the case, \code{file} must be a file produced with the \code{>SAVE PARm
#' file;} command in BILOG.
#' 
#' 
#' @param file File name
#' @return A list with two elements, \code{est} and \code{se}, for the parameter 
#' estimates and their standard errors, correspondingly. Each element is a  
#' matrix with one row per item, and three columns: [,1] item
#' discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and [,3] asymptote
#' \eqn{c}. For the 1PL and 2PL models, all asymptotes are equal to 0; for the
#' 1PL, the discriminations are all equal but not necessarily equal to 1.
#' @author Ivailo Partchev
#' @export
#' @keywords data
read.ip.bilog = function(file) {
  p = read.fwf(file=file, widths=c(8,8,rep(10,13),4,1,1), skip=4, header=FALSE)
  return(list(est=cbind(p[,5],p[,7],p[,11]),se=cbind(p[,6],p[,8],p[,12])))
}

#' Read in parameter estimates
#' 
#' From ICL output, read in estimates of item parameters. Invoked automatically
#' when the model is estimated via \code{irtoys}. If that is not the case,
#' \code{file} must be a file produced with the \code{write_item_param file}
#' command in ICL.
#' 
#' 
#' @param file File name
#' @return A list with two elements, \code{est} and \code{se}, for the parameter 
#' estimates and their standard errors, correspondingly. Because ICL does not
#' compute standard errors, \code{se} will be NULL. \code{est} is a matrix with 
#' one row per item, and three columns: [,1] item
#' discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and [,3] asymptote
#' \eqn{c}. For the 1PL and 2PL models, all asymptotes are equal to 0; for the
#' 1PL, the discriminations are all equal but not necessarily equal to 1.
#' @author Ivailo Partchev
#' @export
#' @keywords models
read.ip.icl = function(file) {
  return(list(est=matrix(scan(file),ncol=4,byrow=TRUE)[,-1],se=NULL))
}


#' Read in quadrature
#' 
#' From ICL output, read in estimates of item parameters. \code{file} must be a
#' file produced with the \code{write_latent_dist file} command in ICL.
#' Quadrature objects are used when estimating abilities with \code{eap} and
#' for some of the scaling methods in \code{sca}.
#' 
#' 
#' @param file File name
#' @return A list of: \item{quad.points}{A vector of quadrature points}
#' \item{quad.weights}{A vector of the corresponding quadrature weights}
#' @author Ivailo Partchev
#' @export
#' @seealso \code{\link{normal.qu}}, \code{\link{eap}}, \code{\link{sca}}
#' @keywords models
read.qu.icl = function(file) {
  f=read.table(file,header=FALSE)
  return(list(quad.points=f$V1, quad.weights=f$V2))
}

#' Read responses from a file
#' 
#' Reads responses to a questionnaire from a text file
#' 
#' Included for those who are too faint-hearted to write
#' \code{as.matrix(read.table(file, head=F))}. Of course, data can be entered
#' into R in many other ways.
#' 
#' The data values in the \code{file} must be separated with blanks.
#' 
#' Responses are the empirical data used in IRT. Note that \code{irtoys} deals
#' with models for dichotomous data, and typically expects data consisting of
#' zeroes and ones, without any missing values (non-responses are considered as
#' wrong responses). In fact, there are only two commands in \code{irtoys} that
#' accept other kinds of data: \code{sco} and \code{tgp}.
#' 
#' \code{read.resp} does accept missing data and values other than 0 and 1. Use
#' \code{sco} and a key to score multiple choice responses to 0/1. If you have
#' dichotomous data that contains NAs, you can use \code{sco} without a key to
#' change all NA to 0.
#' 
#' @param file File name
#' @param na The symbol used to represent missing data
#' @return A matrix, typically of zeroes and ones, representing the correct or
#' wrong responses given by persons (rows) to items (columns).
#' @author Ivailo Partchev
#' @export
#' @seealso \code{\link{sco}}, \code{\link{tgp}},
#' @keywords IO
#' @export
#' @examples
#' 
#' \dontrun{
#'   r <- read.resp("c:/myfiles/irt.dat")
#' }
#' 
read.resp = function(file, na=".") {
  return(as.matrix(read.table(file=file, header=FALSE, na.strings=na)))
}


