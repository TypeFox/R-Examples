#' Apply a simple homogeneity test after Galbraith (2003)
#'
#' A simple homogeneity test for De estimates
#'
#' For details see Galbraith (2003).
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#' @param log \code{\link{logical}} (with default): peform the homogeniety test
#' with (un-)logged data
#' @param \dots further arguments (for internal compatibility only).
#' @return Returns a terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following element:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call}
#'
#' The output should be accessed using the function
#' \code{\link{get_RLum}}
#' @section Function version: 0.2
#' @author Christoph Burow, University of Cologne (Germany)
#' @seealso \code{\link{pchisq}}
#' @references Galbraith, R.F., 2003. A simple homogeneity test for estimates
#' of dose obtained using OSL. Ancient TL 21, 75-77.
#' @examples
#'
#' ## load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ## apply the homogeneity test
#' calc_HomogeneityTest(ExampleData.DeValues$BT998)
#'
#' @export
calc_HomogeneityTest <- function(
  data,
  log=TRUE,
  ...
){

  ##============================================================================##
  ## CONSISTENCY CHECK OF INPUT DATA
  ##============================================================================##

  if(missing(data)==FALSE){
    if(is(data, "data.frame") == FALSE & is(data, "RLum.Results") == FALSE){
      stop("[calc_FiniteMixture] Error: 'data' object has to be of type
           'data.frame' or 'RLum.Results'!")
    } else {
      if(is(data, "RLum.Results") == TRUE){
        data <- get_RLum(data, signature(object = "De.values"))

      }
    }
  }

  ##==========================================================================##
  ## ... ARGUMENTS
  ##==========================================================================##

  extraArgs <- list(...)

  ## set plot main title
  if("verbose" %in% names(extraArgs)) {
    verbose<- extraArgs$verbose
  } else {
    verbose<- TRUE
  }

  ##============================================================================##
  ## CALCULATIONS
  ##============================================================================##

  if(log==TRUE){
    dat<- log(data)
  } else {
    dat<- data
  }

  wi<- 1/dat[2]^2
  wizi<- wi*dat[1]
  mu<- sum(wizi)/sum(wi)
  gi<- wi*(dat[1]-mu)^2

  G<- sum(gi)
  df<- length(wi)-1
  n<- length(wi)
  P<- pchisq(G, df, lower.tail = FALSE)

  ##============================================================================##
  ## OUTPUT
  ##============================================================================##

  if(verbose == TRUE) {
    cat("\n [calc_HomogeneityTest]")
    cat(paste("\n\n ---------------------------------"))
    cat(paste("\n n:                 ", n))
    cat(paste("\n ---------------------------------"))
    cat(paste("\n mu:                ", round(mu,4)))
    cat(paste("\n G-value:           ", round(G,4)))
    cat(paste("\n Degrees of freedom:", df))
    cat(paste("\n P-value:           ", round(P,4)))
    cat(paste("\n ---------------------------------\n\n"))
  }

  ##============================================================================##
  ## RETURN VALUES
  ##============================================================================##

  summary<- data.frame(n=n,g.value=G,df=df,P.value=P)

  call<- sys.call()
  args<- list(log=log)

  newRLumResults.calc_HomogeneityTest <- set_RLum(
    class = "RLum.Results",
    data = list(
      summary=summary,
      data=data,
      args=args,
      call=call
    ))

  invisible(newRLumResults.calc_HomogeneityTest)

}
