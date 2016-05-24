#' Apply the model after Fuchs & Lang (2001) to a given De distribution.
#'
#' This function applies the method according to Fuchs & Lang (2001) for
#' heterogeneously bleached samples with a given coefficient of variation
#' threshold.
#'
#' \bold{Used values} \cr If the coefficient of variation (c[v]) of the first
#' two values is larger than the threshold c[v_threshold], the first value is
#' skipped.  Use the \code{startDeValue} argument to define a start value for
#' calculation (e.g. 2nd or 3rd value).\cr
#'
#' \bold{Basic steps of the approach} \cr
#'
#' (1) Estimate natural relative variation of the sample using a dose recovery
#' test\cr (2) Sort the input values ascendingly\cr (3) Calculate a running
#' mean, starting with the lowermost two values and add values iteratively.\cr
#' (4) Stop if the calculated c[v] exceeds the specified \code{cvThreshold}\cr
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#' @param cvThreshold \link{numeric} (with default): coefficient of variation
#' in percent, as threshold for the method, e.g. \code{cvThreshold = 3}. See
#' details.
#' @param startDeValue \link{numeric} (with default): number of the first
#' aliquot that is used for the calculations
#' @param plot \link{logical} (with default): plot output
#' \code{TRUE}/\code{FALSE}
#' @param \dots further arguments and graphical parameters passed to
#' \code{\link{plot}}
#' @return Returns a plot (optional) and terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following elements:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call}
#' \item{usedDeValues}{\link{data.frame} containing the used values for the
#' calculation}
#' @note Please consider the requirements and the constraints of this method
#' (see Fuchs & Lang, 2001)
#' @section Function version: 0.4.1
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France) Christoph Burow, University of Cologne (Germany)
#' @seealso \code{\link{plot}}, \code{\link{calc_MinDose}},
#' \code{\link{calc_FiniteMixture}}, \code{\link{calc_CentralDose}},
#' \code{\link{calc_CommonDose}}, \code{\linkS4class{RLum.Results}}
#' @references Fuchs, M. & Lang, A., 2001. OSL dating of coarse-grain fluvial
#' quartz using single-aliqout protocols on sediments from NE Peloponnese,
#' Greece. In: Quaternary Science Reviews 20, 783-787.
#'
#' Fuchs, M. & Wagner, G.A., 2003. Recognition of insufficient bleaching by
#' small aliquots of quartz for reconstructing soil erosion in Greece.
#' Quaternary Science Reviews 22, 1161-1167.
#' @keywords dplot
#' @examples
#'
#'
#' ##load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ##calculate De according to Fuchs & Lang (2001)
#' temp<- calc_FuchsLang2001(ExampleData.DeValues$BT998, cvThreshold = 5)
#'
#' @export
calc_FuchsLang2001 <- function(
  data,
  cvThreshold=5,
  startDeValue=1,
  plot=TRUE,
  ...
){

  # Integrity Tests ---------------------------------------------------------

  if(missing(data)==FALSE){
    if(is(data, "data.frame") == FALSE & is(data,"RLum.Results") == FALSE){
      stop("[calc_FuchsLang2001] 'data' has to be of type 'data.frame' or 'RLum.Results'!")
    } else {
      if(is(data, "RLum.Results") == TRUE){
        data <- get_RLum(data,data.object="De.values")
      }
    }
  }

  # Deal with extra arguments -----------------------------------------------
  ##deal with addition arguments
  extraArgs <- list(...)

  verbose <- if("verbose" %in% names(extraArgs)) {extraArgs$verbose} else {TRUE}


  ##============================================================================##
  ##PREPARE DATA
  ##============================================================================##

  ##1. order values in acending order write used D[e] values in data.frame
  o <- order(data[1]) # o is only an order parameter
  data_ordered <- data[o,] # sort values after o and write them into a new variable

  ##2. estimate D[e]

  # set variables
  usedDeValues<-data.frame(De=NA,De_Error=NA,cv=NA)
  endDeValue<-startDeValue

  # if the frist D[e] values are not used write this information in the data.frame
  if (startDeValue!=1) {

    n <- abs(1-startDeValue)

    #  write used D[e] values in data.frame
    usedDeValues[1:n,1]<-data_ordered[1:n,1]
    usedDeValues[1:n,2]<-data_ordered[1:n,2]
    usedDeValues[1:n,3]<-"skipped"
  }

  ##=================================================================================================##
  ##LOOP FOR MODEL
  ##=================================================================================================##

  # repeat loop (run at least one time)
  repeat {

    #calculate mean, sd and cv
    mean<-round(mean(data_ordered[startDeValue:endDeValue,1]),digits=2) #calculate mean from ordered D[e] values
    sd<-round(sd(data_ordered[startDeValue:endDeValue,1]),digits=2)		#calculate sd from ordered D[e] values
    cv<-round(sd/mean*100, digits=2) #calculate coefficent of variation


    # break if cv > cvThreshold
    if (cv>cvThreshold & endDeValue>startDeValue){

      # if the first two D[e] values give a cv > cvThreshold, than skip the first D[e] value
      if (endDeValue-startDeValue<2) {

        #  write used D[e] values in data.frame
        usedDeValues[endDeValue,1]<-data_ordered[endDeValue,1]
        usedDeValues[endDeValue,2]<-data_ordered[endDeValue,2]
        usedDeValues[endDeValue-1,3]<-"not used"

        # go to the next D[e] value
        startDeValue<-startDeValue+1

      } else {

        usedDeValues[endDeValue,1]<-data_ordered[endDeValue,1]
        usedDeValues[endDeValue,2]<-data_ordered[endDeValue,2]
        usedDeValues[endDeValue,3]<-paste("# ",cv," %",sep="")

        break #break loop
      }

    }#EndIf
    else {

      # write used D[e] values in data.frame
      usedDeValues[endDeValue,1]<-data_ordered[endDeValue,1]
      usedDeValues[endDeValue,2]<-data_ordered[endDeValue,2]

      # first cv values alway contains NA to ensure that NA% is not printed test
      if(is.na(cv)==TRUE) {
        usedDeValues[endDeValue,3]<-cv
      } else {
        usedDeValues[endDeValue,3]<-paste(cv," %",sep="")
      }
    }#EndElse

    # go the next D[e] value until the maximum number is reached
    if (endDeValue<length(data_ordered[,1])) {
      endDeValue<-endDeValue+1
    } else {break}

  }#EndRepeat

  ##=================================================================================================##
  ##ADDITIONAL CALCULATIONS and TERMINAL OUTPUT
  ##=================================================================================================##

  # additional calculate weighted mean
  w<-1/(data_ordered[startDeValue:endDeValue,2])^2 #weights for weighted mean
  weighted_mean <- round(weighted.mean(data_ordered[startDeValue:endDeValue,1], w), digits=2)
  weighted_sd<-round(sqrt(1/sum(w)),digits=2)
  n.usedDeValues<-endDeValue-startDeValue+1

  # standard error
  se<- round(sd/sqrt(endDeValue-startDeValue+1), digits=2)

  if(verbose==TRUE){
    cat("\n [calc_FuchsLang2001]")
    cat(paste("\n\n----------- meta data --------------"))
    cat(paste("\n cvThreshold:            ",cvThreshold,"%"))
    cat(paste("\n used values:            ",n.usedDeValues))
    cat(paste("\n----------- dose estimate ----------"))
    cat(paste("\n mean:                   ",mean))
    cat(paste("\n sd:                     ",sd))
    cat(paste("\n weighted mean:          ",weighted_mean))
    cat(paste("\n weighted sd:            ",weighted_sd))
    cat(paste("\n------------------------------------\n\n"))
  }

  ##=================================================================================================#
  ##RETURN  VALUES
  ##=================================================================================================##

  summary<- data.frame(de=mean,
                       de_err=sd,
                       de_weighted=weighted_mean,
                       de_weighted_err=weighted_sd,
                       n.usedDeValues=n.usedDeValues)

  call<- sys.call()
  args<- list(cvThreshold = cvThreshold, startDeValue = startDeValue)

  newRLumResults.calc_FuchsLang2001<- set_RLum(
    class = "RLum.Results",
    data = list(summary = summary,
                data = data,
                args = args,
                call = call,
                usedDeValues=usedDeValues))

  ##=========##
  ## PLOTTING
  if(plot==TRUE) {
    try(plot_RLum.Results(newRLumResults.calc_FuchsLang2001, ...))
  }#endif::plot

  invisible(newRLumResults.calc_FuchsLang2001)

}
