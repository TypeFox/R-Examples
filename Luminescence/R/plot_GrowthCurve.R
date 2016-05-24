#' Fit and plot a growth curve for luminescence data (Lx/Tx against dose)
#'
#' A dose response curve is produced for luminescence measurements using a
#' regenerative protocol.
#'
#' \bold{Fitting methods} \cr\cr For all options (except for the \code{LIN}, \code{QDR} and
#' the \code{EXP OR LIN}), the \code{\link[minpack.lm]{nlsLM}} function with the
#' \code{LM} (Levenberg-Marquardt algorithm) algorithm is used. Note: For historical reasons
#' for the Monte Carlo simulations partly  the function \code{\link{nls}} using the \code{port} algorithm.
#'
#' The solution is found by transforming the function or using \code{\link{uniroot}}. \cr
#'
#' \code{LIN}: fits a linear function to the data using
#' \link{lm}: \deqn{y = m*x+n}
#'
#' \code{QDR}: fits a linear function to the data using
#' \link{lm}: \deqn{y = a + b * x + c * x^2}
#'
#' \code{EXP}: try to fit a function of the form
#' \deqn{y = a*(1-exp(-(x+c)/b))} Parameters b and c are approximated by a
#' linear fit using \link{lm}. Note: b = D0\cr
#'
#' \code{EXP OR LIN}: works for some cases where an \code{EXP} fit fails. If
#' the \code{EXP} fit fails, a \code{LIN} fit is done instead. \cr
#'
#' \code{EXP+LIN}: tries to fit an exponential plus linear function of the
#' form: \deqn{y = a*(1-exp(-(x+c)/b)+(g*x))} The De is calculated by
#' iteration.\cr \bold{Note:} In the context of luminescence dating, this
#' function has no physical meaning. Therefore, no D0 value is returned.\cr
#'
#' \code{EXP+EXP}: tries to fit a double exponential function of the form
#' \deqn{y = (a1*(1-exp(-(x)/b1)))+(a2*(1-exp(-(x)/b2)))} This fitting
#' procedure is not robust against wrong start parameters and should be further
#' improved.\cr\cr
#'
#' \bold{Fit weighting}\cr
#'
#' If the option \code{fit.weights =  TRUE} is chosen, weights are calculated using
#' provided signal errors (Lx/Tx error): \deqn{fit.weights = 1/error/(sum(1/error))}\cr
#'
#' \bold{Error estimation using Monte Carlo simulation}\cr
#'
#' Error estimation is done using a Monte Carlo (MC) simulation approach. A set of Lx/Tx values is
#' constructed by randomly drawing curve data from samled from normal
#' distributions. The normal distribution is defined by the input values (mean
#' = value, sd = value.error). Then, a growth curve fit is attempted for each
#' dataset resulting in a new distribution of single De values. The \link{sd}
#' of this distribution is becomes then the error of the De. With increasing
#' iterations, the error value becomes more stable. \bold{Note:} It may take
#' some calculation time with increasing MC runs, especially for the composed
#' functions (\code{EXP+LIN} and \code{EXP+EXP}).\cr Each error estimation is
#' done with the function of the chosen fitting method. \cr
#'
#' \bold{Subtitle information}\cr
#'
#' To avoid plotting the subtitle information, provide an empty user mtext \code{mtext = ""}.
#' To plot any other subtitle text, use \code{mtext}.
#'
#' @param sample \code{\link{data.frame}} (\bold{required}): data frame with
#' three columns for x=Dose,y=LxTx,z=LxTx.Error, y1=TnTx. The column for the
#' test dose response is optional, but requires 'TnTx' as column name if used.
#'
#' @param na.rm \code{\link{logical}} (with default): excludes \code{NA} values
#' from the data set prior to any further operations.
#'
#' @param fit.method \code{\link{character}} (with default): function used for
#' fitting. Possible options are: \code{LIN}, \code{QDR}, \code{EXP}, \code{EXP OR LIN},
#' \code{EXP+LIN} or \code{EXP+EXP}. See details.
#'
#' @param fit.force_through_origin \code{\link{logical}} (with default) allow to force
#' the fitted function through the origin. For \code{method = "EXP+EXP"} the function will
#' go to the origin in either case, so this option will have no effect.
#'
#' @param fit.weights \code{\link{logical}} (with default): option whether the
#' fitting is done with or without weights. See details.
#'
#' @param fit.includingRepeatedRegPoints \code{\link{logical}} (with default):
#' includes repeated points for fitting (\code{TRUE}/\code{FALSE}).
#'
#' @param fit.NumberRegPoints \code{\link{integer}} (optional): set number of
#' regeneration points manually. By default the number of all (!) regeneration
#' points is used automatically.
#'
#' @param fit.NumberRegPointsReal \code{\link{integer}} (optional): if the
#' number of regeneration points is provided manually, the value of the real,
#' regeneration points = all points (repeated points) including reg 0, has to
#' be inserted.
#'
#' @param fit.bounds \code{\link{logical}} (with default): set lower fit bounds
#' for all fitting parameters to 0. Limited for the use with the fit methods
#' \code{EXP}, \code{EXP+LIN} and \code{EXP OR LIN}. Argument to be inserted
#' for experimental application only!
#'
#' @param NumberIterations.MC \code{\link{integer}} (with default): number of
#' Monte Carlo simulations for error estimation. See details.
#'
#' @param output.plot \code{\link{logical}} (with default): plot output
#' (\code{TRUE/FALSE}).
#'
#' @param output.plotExtended \code{\link{logical}} (with default): If
#' \code{TRUE}, 3 plots on one plot area are provided: (1) growth curve, (2)
#' histogram from Monte Carlo error simulation and (3) a test dose response
#' plot. If \code{FALSE}, just the growth curve will be plotted.
#' \bold{Requires:} \code{output.plot = TRUE}.
#'
#' @param output.plotExtended.single \code{\link{logical}} (with default):
#' single plot output (\code{TRUE/FALSE}) to allow for plotting the results in
#' single plot windows. Requires \code{output.plot = TRUE} and
#' \code{output.plotExtended = TRUE}.
#'
#' @param cex.global \code{\link{numeric}} (with default): global scaling
#' factor.
#'
#' @param txtProgressBar \code{\link{logical}} (with default): enables or disables txtProgressBar.
#' If \code{verbose = FALSE} also no txtProgressBar is shown.
#'
#' @param verbose \code{\link{logical}} (with default): enables or disables terminal feedback.
#'
#' @param \dots Further arguments and graphical parameters to be passed. Note:
#' Standard arguments will only be passed to the growth curve plot. Supported:
#' \code{xlim}, \code{ylim}, \code{main}, \code{xlab}, \code{ylab}
#'
#' @return Along with a plot (so far wanted) an \code{RLum.Results} object is returned containing,
#' the slot \code{data} contains the following elements:\cr
#'
#' \tabular{lll}{
#' \bold{DATA.OBJECT} \tab \bold{TYPE} \tab \bold{DESCRIPTION} \cr
#' \code{..$De} : \tab  \code{data.frame} \tab Table with De values \cr
#' \code{..$De.MC} : \tab \code{numeric} \tab Table with De values from MC runs \cr
#' \code{..$Fit} : \tab \code{\link{nls}} or \code{\link{lm}} \tab object from the fitting for \code{EXP},
#' \code{EXP+LIN} and \code{EXP+EXP}. In case of a resulting  linear fit when using \code{LIN}, \code{QDR} or
#' \code{EXP OR LIN} \cr
#' \code{..$Formula} : \tab \code{\link{expression}} \tab Fitting formula as R expression \cr
#' \code{..$call} : \tab \code{call} \tab The original function call\cr
#' }
#'
#' @section Function version: 1.8.1
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France), \cr Michael Dietze, GFZ Potsdam (Germany)
#'
#' @seealso \code{\link{nls}}, \code{\linkS4class{RLum.Results}},
#' \code{\link{get_RLum}}, \code{\link[minpack.lm]{nlsLM}}, \code{\link{lm}}, \code{uniroot}
#'
#' @examples
#'
#' ##(1) plot growth curve for a dummy data.set and show De value
#' data(ExampleData.LxTxData, envir = environment())
#' temp <- plot_GrowthCurve(LxTxData)
#' get_RLum(temp)
#'
#' ##(1a) to access the fitting value try
#' get_RLum(temp, data.object = "Fit")
#'
#' ##(2) plot the growth curve only - uncomment to use
#' ##pdf(file = "~/Desktop/Growth_Curve_Dummy.pdf", paper = "special")
#' plot_GrowthCurve(LxTxData)
#' ##dev.off()
#'
#' ##(3) plot growth curve with pdf output - uncomment to use, single output
#' ##pdf(file = "~/Desktop/Growth_Curve_Dummy.pdf", paper = "special")
#' plot_GrowthCurve(LxTxData, output.plotExtended.single = TRUE)
#' ##dev.off()
#'
#' @export
plot_GrowthCurve <- function(
  sample,
  na.rm = TRUE,
  fit.method = "EXP",
  fit.force_through_origin = FALSE,
  fit.weights = TRUE,
  fit.includingRepeatedRegPoints = TRUE,
  fit.NumberRegPoints,
  fit.NumberRegPointsReal,
  fit.bounds = TRUE,
  NumberIterations.MC = 100,
  output.plot = TRUE,
  output.plotExtended = TRUE,
  output.plotExtended.single = FALSE,
  cex.global = 1,
  txtProgressBar = TRUE,
  verbose = TRUE,
  ...
) {

  ##1. check if sample is data.frame
  if(is.data.frame(sample)==FALSE){
    stop("\n [plot_GrowthCurve()] Sample has to be of type data.fame!")
  }

  ##2. check if sample contains a least three rows
  if(length(sample[,1])<3){
    stop("\n [plot_GrowthCurve()] At least two regeneration points are needed!")
  }

  ## optionally, count and exclude NA values and print result
  if(na.rm == TRUE) {
    n.NA <- sum(!complete.cases(sample))

    if (n.NA == 1) {
      warning("[plot_GrowthCurve()] 1 NA value excluded.")
    } else if (n.NA > 1) {
      warning(paste(" [plot_GrowthCurve()]", n.NA, "NA values excluded."))
    }

    sample <- na.exclude(sample)

    ##Check if anything is left after removal
    if(nrow(sample) == 0){

      warning("[plot_GrowthCurve()] Sorry, after NA removal nothing is left from the data set! NULL returned")
      return(NULL)
    }

  }

  ##3. verbose mode
  if(!verbose){
    txtProgressBar <- FALSE
  }


  ##remove rownames from data.frame, as this could causes errors for the reg point calculation
  rownames(sample) <- NULL


  ##NULL values in the data.frame are not allowed for the y-column
  if(length(sample[sample[,2]==0,2])>0){
    sample[sample[,2]==0,2]<-0.0001
    warning(paste("[plot_GrowthCurve()]", length(sample[sample[,2]==0,2]), "values with 0 for Lx/Tx detected; replaced by 0.0001."))
  }

  ##1. INPUT

  #1.0.1 calculate number of reg points if not set
  if(missing(fit.NumberRegPoints)==TRUE){fit.NumberRegPoints<-length(sample[-1,1])}
  if(missing(fit.NumberRegPointsReal)==TRUE){

    fit.RegPointsReal <- as.integer(
      rownames(sample[-which(duplicated(sample[,1]) | sample[,1]==0),]))

    fit.NumberRegPointsReal <- length(fit.RegPointsReal)

  }

  #1.1 Produce dataframe from input values
  xy<-data.frame(x=sample[2:(fit.NumberRegPoints+1),1],y=sample[2:(fit.NumberRegPoints+1),2])
  y.Error<-sample[2:(fit.NumberRegPoints+1),3]

  ##1.1.1 produce weights for weighted fitting
  if(fit.weights){

    fit.weights <- 1 / y.Error / (sum(1 / y.Error))

    if(is.na(fit.weights[1])){

      fit.weights <- NULL
      warning("fit.weights set to NULL as error column is invalid or 0.")

    }

  }else{
    fit.weights <- rep(1, length(y.Error))

  }


  # Deal with extra arguments -----------------------------------------------
  ##deal with addition arguments
  extraArgs <- list(...)

  main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
  {"Growth curve"}

  xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
  {"Dose [s]"}

  ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
  {expression(L[x]/T[x])}

  if("cex" %in% names(extraArgs)) {cex.global <- extraArgs$cex}

  ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else
  {c(min(xy$y)-max(y.Error),(max(xy$y)+if(max(xy$y)*0.1>1.5){1.5}else{max(xy$y)*0.2}))}

  xlim <- if("xlim" %in% names(extraArgs)) {extraArgs$xlim} else
  {c(0,(max(xy$x)+if(max(xy$x)*0.4>50){50}else{max(xy$x)*0.4}))}

  fun   <- if("fun" %in% names(extraArgs)) {extraArgs$fun} else {FALSE}


  #1.2 Prepare data sets regeneration points for MC Simulation
  data.MC<-t(matrix(sapply(seq(2,fit.NumberRegPoints+1,by=1),
                           function(x){sample(rnorm(10000,mean=sample[x,2], sd=abs(sample[x,3])),
                                              NumberIterations.MC, replace=TRUE)}), nrow=NumberIterations.MC
  )#end matrix
  )#end transpose matrix

  #1.3 Do the same for the natural signal
  data.MC.De <- numeric(NumberIterations.MC)
  data.MC.De <- sample(rnorm(10000,mean=sample[1,2], sd=abs(sample[1,3])),
                       NumberIterations.MC, replace=TRUE)

  #1.3 set x.natural
  x.natural<-numeric(NumberIterations.MC)
  x.natural <- NA

  ##1.4 set initialise variables
  De <- NA
  De.Error <- NA


  ##============================================================================##
  # FITTING ----------------------------------------------------------------------
  ##============================================================================##
  ##3. Fitting values with nonlinear least-squares estimation of the parameters

  ##set functions for fitting

  #EXP
  fit.functionEXP <- function(a,b,c,x) {a*(1-exp(-(x+c)/b))}
  fit.formulaEXP <- y ~ a * (1 - exp(-(x+c)/b))

  #EXP+LIN
  fit.functionEXPLIN<-function(a,b,c,g,x) {a*(1-exp(-(x+c)/b)+(g*x))}
  fit.formulaEXPLIN <- y ~ a*(1-exp(-(x+c)/b)+(g*x))

  #EXP+EXP
  fit.functionEXPEXP<-function(a1,a2,b1,b2,x){(a1*(1-exp(-(x)/b1)))+(a2*(1-exp(-(x)/b2)))}
  fit.formulaEXPEXP <- y ~ (a1*(1-exp(-(x)/b1)))+(a2*(1-exp(-(x)/b2)))

  ##input data for fitting; exclude repeated RegPoints
  if (fit.includingRepeatedRegPoints == FALSE) {
    data <-
      data.frame(x = xy[-which(duplicated(xy[,1])),1],y = xy[-which(duplicated(xy[,1])),2])
  }else{
    data <- data.frame(xy)
  }


  ##START PARAMETER ESTIMATION
  ##--------------------------------------------------------------------------##
  ##general setting of start parameters for fitting

  ##a - estimation for a a the maxium of the y-values (Lx/Tx)
  a <- max(data[,2])

  ##b - get start parameters from a linear fit of the log(y) data
  ##    (suppress the warning in case one parameter is negative)

  fit.lm <- lm(suppressWarnings(log(data$y))~data$x)
  b <- as.numeric(1/fit.lm$coefficients[2])

  ##c - get start parameters from a linear fit - offset on x-axis
  fit.lm<-lm(data$y~data$x)
  c <- as.numeric(abs(fit.lm$coefficients[1]/fit.lm$coefficients[2]))

  #take slope from x - y scaling
  g <- max(data[,2]/max(data[,1]))

  #set D01 and D02 (in case of EXp+EXP)
  D01 <- NA; D02 <- NA

  ##--------------------------------------------------------------------------##
  ##to be a little bit more flexible the start parameters varries within a normal distribution

  ##draw 50 start values from a normal distribution a start values
  a.MC<-rnorm(50,mean=a,sd=a/100)

  if(!is.na(b)) {
    b.MC <- rnorm(50,mean = b,sd = b / 100)
  }else{
    b <- NA

  }

  c.MC<-rnorm(50,mean=c,sd=c/100)
  g.MC<-rnorm(50,mean=g,sd=g/1)

  ##set start vector (to avoid errors witin the loop)
  a.start <- NA; b.start <- NA; c.start <- NA; g.start <- NA



  ##--------------------------------------------------------------------------##
  #===========================================================================##
  #QDR#
  if (fit.method == "QDR"){

    ##Do fitting with option to force curve through the origin
    if(fit.force_through_origin){

      ##linear fitting ... polynomial
      fit  <- lm(data$y ~  0 + I(data$x) + I(data$x^2), weights = fit.weights)

      ##give function for uniroot
      De.fs <- function(x, y) {
        0 + coef(fit)[1] * x + coef(fit)[2] * x ^ 2 - y

      }


    }else{


      ##linear fitting ... polynomial
      fit  <- lm(data$y ~  I(data$x) + I(data$x^2), weights = fit.weights)

      ##give function for uniroot
      De.fs <- function(x, y) {
        coef(fit)[1] + coef(fit)[2] * x + coef(fit)[3] * x ^ 2 - y

      }

    }

    ##solve and get De
    De.uniroot <- try(
      uniroot(De.fs, y = sample[1,2], lower = 0, upper = max(sample[,1]) * 1.5), silent = TRUE)

    if(!inherits(De.uniroot, "try-error")){
      De <- round(De.uniroot$root, digits = 2)
      if(verbose) writeLines(paste0("[plot_GrowthCurve()] >> De = ", De))

    }else{
      if(verbose) writeLines("[plot_GrowthCurve()] >> no solution found for QDR fit")
      De <- NA

    }


    # +++++++++++++++++++++++++++++++++++++++++

    ##terminal output fo MC
    if(verbose){
      cat("\n\t Run Monte Carlo loops for error estimation of the QDR fit\n")
    }

    ##set progressbar
    if(txtProgressBar){
      pb<-txtProgressBar(min=0,max=NumberIterations.MC, char="=", style=3)
    }


    #start loop for Monte Carlo Error estimation
    fit.MC <- sapply(1:NumberIterations.MC, function(i){

      data <- data.frame(x=xy$x, y=data.MC[,i])

      if(fit.force_through_origin){

        ##linear fitting ... polynomial
        fit.MC  <- lm(data$y ~  0 + I(data$x) + I(data$x^2), weights = fit.weights)

        ##give function for uniroot
        De.fs.MC <- function(x, y) {
          0 + coef(fit.MC)[1] * x + coef(fit.MC)[2] * x ^ 2 - y

        }


      }else{


        ##linear fitting ... polynomial
        fit.MC  <- lm(data$y ~  I(data$x) + I(data$x^2), weights = fit.weights)

        ##give function for uniroot
        De.fs.MC <- function(x, y) {
          coef(fit.MC)[1] + coef(fit.MC)[2] * x + coef(fit.MC)[3] * x ^ 2 - y

        }

      }

      ##solve and get De
      De.uniroot.MC <- try(uniroot(
        De.fs.MC,
        y = data.MC.De[i],
        lower = 0,
        upper = max(sample[, 1]) * 1.5
      ), silent = TRUE)

      if(!inherits(De.uniroot.MC, "try-error")){
        De.MC <- round(De.uniroot.MC$root, digits = 2)

      }else{
        De.MC <- NA

      }

      ##update progress bar
      if(txtProgressBar) setTxtProgressBar(pb, i)

      return(De.MC)

    })

    if(txtProgressBar) close(pb)

    x.natural<- fit.MC
  }
  #===========================================================================##
  #EXP#

  if (fit.method=="EXP" | fit.method=="EXP OR LIN" | fit.method=="LIN"){

    if(fit.method!="LIN" & length(data[,1])>=2){

      ##FITTING on GIVEN VALUES##
      #	--use classic R fitting routine to fit the curve

      ##try to create some start parameters from the input values to make
      ## the fitting more stable
      for(i in 1:50){

        a<-a.MC[i];b<-b.MC[i];c<-c.MC[i]

        fit.initial<-try(nls(y~fit.functionEXP(a,b,c,x),
                             data=data,
                             start=c(a=a,b=b,c=c),
                             trace=FALSE,
                             algorithm="port",
                             lower=c(a=0,b>0,c=0),
                             nls.control(maxiter=100,warnOnly=FALSE,minFactor=1/2048) #increase max. iterations
        ),silent=TRUE)

        if(class(fit.initial)!="try-error"){
          #get parameters out of it
          parameters<-(coef(fit.initial))
          b.start[i]<-as.vector((parameters["b"]))
          a.start[i]<-as.vector((parameters["a"]))
          c.start[i]<-as.vector((parameters["c"]))
        }
      }


      ##used median as start parameters for the final fitting
      a <- median(na.exclude(a.start))
      b <- median(na.exclude(b.start))
      c <- median(na.exclude(c.start))

      #FINAL Fit curve on given values
      fit <- try(minpack.lm::nlsLM(
        formula = fit.formulaEXP,
        data = data,
        start = list(a = a, b = b,c = c),
        weights = fit.weights,
        trace = FALSE,
        algorithm = "LM",
        lower = if (fit.bounds) {
          c(0,0,0)
        }else{
          c(-Inf,-Inf,-Inf)
        },
        upper = if (fit.force_through_origin) {
          c(Inf, Inf, 0)
        }else{
          c(Inf, Inf, Inf)
        },
        control = minpack.lm::nls.lm.control(maxiter = 500)
      ), silent = TRUE
      )

      if (inherits(fit, "try-error") & inherits(fit.initial, "try-error")){

        if(verbose) writeLines("[plot_GrowthCurve()] >> try-error for EXP fit")

      }else{

        ##this is to avoid the singular convergence failure due to a perfect fit at the beginning
        ##this may happen especially for simulated data
        if(inherits(fit, "try-error") & !inherits(fit.initial, "try-error")){
          fit <- fit.initial
          rm(fit.initial)

        }

        #get parameters out of it
        parameters<-(coef(fit))
        b<-as.vector((parameters["b"]))
        a<-as.vector((parameters["a"]))
        c<-as.vector((parameters["c"]))


        #calculate De
        De<-suppressWarnings(round(-c-b*log(1-sample[1,2]/a), digits=2))

        #print D01 value
        D01<-round(b,digits=2)
        if(verbose) writeLines(paste0("[plot_GrowthCurve()] >> D01 = ",D01, " | De = ", De))


        ##Monte Carlo Simulation
        #	--Fit many curves and calculate a new De +/- De_Error
        #	--take De_Error

        #set variables
        var.b<-vector(mode="numeric", length=NumberIterations.MC)
        var.a<-vector(mode="numeric", length=NumberIterations.MC)
        var.c<-vector(mode="numeric", length=NumberIterations.MC)

        #start loop
        for (i in 1:NumberIterations.MC) {

          ##set data set
          data <- data.frame(x = xy$x,y = data.MC[,i])

          fit.MC <- try(minpack.lm::nlsLM(
            formula = fit.formulaEXP,
            data = data,
            start = list(a = a, b = b,c = c),
            weights = fit.weights,
            trace = FALSE,
            algorithm = "LM",
            lower = if (fit.bounds) {
              c(0,0,0)
            }else{
              c(-Inf,-Inf,-Inf)
            },
            upper = if (fit.force_through_origin) {
              c(Inf, Inf, 0)
            }else{
              c(Inf, Inf, Inf)
            },
            control = minpack.lm::nls.lm.control(maxiter = 500)
          ), silent = TRUE
          )

          #get parameters out of it including error handling
          if (class(fit.MC)=="try-error") {

            x.natural[i]<-NA

          }else {

            #get parameters out
            parameters<-coef(fit.MC)
            var.b[i]<-as.vector((parameters["b"]))
            var.a[i]<-as.vector((parameters["a"])) #Imax
            var.c[i]<-as.vector((parameters["c"]))

            #calculate x.natural for error calculatio
            x.natural[i]<-suppressWarnings(
              round(-var.c[i]-var.b[i]*log(1-data.MC.De[i]/var.a[i]), digits=2))

          }

        }#end for loop
      }#endif::try-error fit
    }#endif:fit.method!="LIN"
    #========================================================================
    #LIN#
    ##two options: just linear fit or LIN fit after the EXP fit failed

    #set fit object, if fit objekt was not set before
    if(exists("fit")==FALSE){fit<-NA}

    if ((fit.method=="EXP OR LIN" & class(fit)=="try-error") |
        fit.method=="LIN" | length(data[,1])<3) {

      ##Do fitting again as just allows fitting through the origin
      if(fit.force_through_origin){

        fit.lm<-lm(data$y ~ 0 + data$x, weights = fit.weights)

        #calculate De
        De <- round((sample[1,2]/fit.lm$coefficients[1]), digits=2)


      }else{

        fit.lm<-lm(data$y ~ data$x, weights = fit.weights)

        #calculate De
        De <- round((sample[1,2]-fit.lm$coefficients[1])/fit.lm$coefficients[2], digits=2)

      }


      ##remove vector labels
      De <- as.numeric(as.character(De))


      #start loop for Monte Carlo Error estimation
      for (i in 1:NumberIterations.MC) {

        data <- data.frame(x=xy$x, y=data.MC[,i])

        if(fit.force_through_origin){

          ##do fitting
          fit.lmMC <- lm(data$y ~ 0 + data$x, weights=fit.weights)

          #calculate x.natural
          x.natural[i]<-round((data.MC.De[i]/fit.lmMC$coefficients[1]), digits=2)


        }else{

          ##do fitting
          fit.lmMC <- lm(data$y~ data$x, weights=fit.weights)


          #calculate x.natural
          x.natural[i]<-round((data.MC.De[i]-fit.lmMC$coefficients[1])/
                                fit.lmMC$coefficients[2], digits=2)

        }




      }#endfor::loop for MC

      #correct for fit.method
      fit.method<-"LIN"

      ##set fit object
      if(fit.method=="LIN"){fit<-fit.lm}

    }else{fit.method<-"EXP"}#endif::LIN
  }#end if EXP (this includes the LIN fit option)
  #===========================================================================
  #===========================================================================
  #EXP+LIN#
  else if (fit.method=="EXP+LIN") {


    ##try some start parameters from the input values to makes the fitting more stable
    for(i in 1:length(a.MC)){

      a<-a.MC[i];b<-b.MC[i];c<-c.MC[i];g<-g.MC[i]

      ##---------------------------------------------------------##
      ##start: with EXP function
      fit.EXP<-try(nls(y~fit.functionEXP(a,b,c,x),
                       data=data,
                       start=c(a=a,b=b,c=c),
                       trace=FALSE,
                       algorithm="port",
                       lower=c(a=0,b>10,c=0),
                       nls.control(maxiter=100,warnOnly=FALSE,minFactor=1/1048)
      ),silent=TRUE)


      if(class(fit.EXP)!="try-error"){
        #get parameters out of it
        parameters<-(coef(fit.EXP))
        b<-as.vector((parameters["b"]))
        a<-as.vector((parameters["a"]))
        c<-as.vector((parameters["c"]))

        ##end: with EXP function
        ##---------------------------------------------------------##
      }


      fit<-try(nls(y~fit.functionEXPLIN(a,b,c,g,x),
                   data=data,
                   start=c(a=a,b=b,c=c,g=g),
                   trace=FALSE,
                   algorithm="port",
                   lower = if(fit.bounds==TRUE){lower=c(a=0,b>10,c=0,g=0)}else{c()},
                   nls.control(maxiter=500,warnOnly=FALSE,minFactor=1/2048) #increase max. iterations
      ),silent=TRUE)

      if(class(fit)!="try-error"){
        #get parameters out of it
        parameters<-(coef(fit))
        b.start[i]<-as.vector((parameters["b"]))
        a.start[i]<-as.vector((parameters["a"]))
        c.start[i]<-as.vector((parameters["c"]))
        g.start[i]<-as.vector((parameters["g"]))
      }



    }##end for loop


    ##used mean as start parameters for the final fitting
    a<-median(na.exclude(a.start))
    b<-median(na.exclude(b.start))
    c<-median(na.exclude(c.start))
    g<-median(na.exclude(g.start))

    ##perform final fitting
    fit <- try(minpack.lm::nlsLM(
      formula = fit.formulaEXPLIN,
      data = data,
      start = list(a = a, b = b,c = c, g = g),
      weights = fit.weights,
      trace = FALSE,
      algorithm = "LM",
      lower = if (fit.bounds) {
        c(0,10,0,0)
      }else{
        c(-Inf,-Inf,-Inf,-Inf)
      },
      upper = if (fit.force_through_origin) {
        c(Inf, Inf, 0, Inf)
      }else{
        c(Inf, Inf, Inf, Inf)
      },
      control = minpack.lm::nls.lm.control(maxiter = 500)
    ), silent = TRUE
    )


    #if try error stop calculation
    if(class(fit)!="try-error"){

      #get parameters out of it
      parameters<-(coef(fit))
      b<-as.vector((parameters["b"]))
      a<-as.vector((parameters["a"]))
      c<-as.vector((parameters["c"]))
      g<-as.vector((parameters["g"]))

      #problem: analytically it is not easy to calculate x,
      #use uniroot to solve that problem ... readjust function first
      f.unirootEXPLIN <- function(a,b,c,g,x,LnTn){fit.functionEXPLIN(a,b,c,g,x)-LnTn}

      temp.De <-  try(uniroot(f = f.unirootEXPLIN,
                              interval = c(0,max(xy$x)*1.5),
                              tol = 0.001,
                              a = a,
                              b = b,
                              c = c,
                              g = g,
                              LnTn = sample[1,2],
                              extendInt = "yes",
                              maxiter = 3000), silent = TRUE)



      if (class(temp.De) != "try-error") {
        De <- round(temp.De$root, digits = 2)
      }else{
        De <- NA
      }


      if(verbose) writeLines(paste0("[plot_GrowthCurve()] >> De = ", De))

      ##Monte Carlo Simulation for error estimation
      #	--Fit many curves and calculate a new De +/- De_Error
      #	--take De_Error

      #set variables
      var.b <- vector(mode="numeric", length=NumberIterations.MC)
      var.a <- vector(mode="numeric", length=NumberIterations.MC)
      var.c <- vector(mode="numeric", length=NumberIterations.MC)
      var.g <- vector(mode="numeric", length=NumberIterations.MC)

      ##terminal output fo MC
      if(verbose){
        cat("\n\t Run Monte Carlo loops for error estimation of the EXP+LIN fit\n")
      }

      ##set progressbar
      if(txtProgressBar){
        pb<-txtProgressBar(min=0,max=NumberIterations.MC, char="=", style=3)
      }


      #start Monto Carlo loops
      for(i in  1:NumberIterations.MC){

        data <- data.frame(x=xy$x,y=data.MC[,i])

        ##perform MC fitting
        fit.MC <- try(minpack.lm::nlsLM(
          formula = fit.formulaEXPLIN,
          data = data,
          start = list(a = a, b = b,c = c, g = g),
          weights = fit.weights,
          trace = FALSE,
          algorithm = "LM",
          lower = if (fit.bounds) {
            c(0,10,0,0)
          }else{
            c(-Inf,-Inf,-Inf, -Inf)
          },
          control = minpack.lm::nls.lm.control(maxiter = 500)
        ), silent = TRUE
        )

        #get parameters out of it including error handling
        if (class(fit.MC)=="try-error") {

          x.natural[i]<-NA

        }else {
          parameters <- coef(fit.MC)
          var.b[i] <- as.vector((parameters["b"]))
          var.a[i] <- as.vector((parameters["a"]))
          var.c[i] <- as.vector((parameters["c"]))
          var.g[i] <- as.vector((parameters["g"]))

          #problem: analytical it is not easy to calculate x,
          #use uniroot to solve this problem

          temp.De.MC <-  try(uniroot(
            f = f.unirootEXPLIN,
            interval = c(0,max(xy$x) * 1.5),
            tol = 0.001,
            a = var.a[i],
            b = var.b[i],
            c = var.c[i],
            g = var.g[i],
            LnTn = data.MC.De[i]
          ), silent = TRUE)

          if (class(temp.De.MC) != "try-error") {
            x.natural[i] <- temp.De.MC$root
          }else{
            x.natural[i] <- NA
          }


        }
        ##update progress bar
        if(txtProgressBar) setTxtProgressBar(pb, i)

      }#end for loop

      ##close
      if(txtProgressBar) close(pb)

    }else{

      #print message
      if(verbose) writeLines(paste0("[plot_GrowthCurve()] >> FITTING FAILED"))

    } #end if "try-error" Fit Method

  } #End if EXP+LIN
  #==========================================================================
  #===========================================================================
  #EXP+EXP#
  else if (fit.method=="EXP+EXP") {

    a1.start <- NA
    a2.start <- NA
    b1.start <- NA
    b2.start <- NA

    ## try to create some start parameters from the input values to make the fitting more stable
    for(i in 1:50) {
      a1 <- a.MC[i];b1 <- b.MC[i];
      a2 <- a.MC[i] / 2; b2 <- b.MC[i] / 2

      fit.start <- try(nls(
        y ~ fit.functionEXPEXP(a1,a2,b1,b2,x),
        data = data,
        start = c(
          a1 = a1,a2 = a2,b1 = b1,b2 = b2
        ),
        trace = FALSE,
        algorithm = "port",
        lower = c(a1 > 0,a2 > 0,b1 > 0,b2 > 0),
        nls.control(
          maxiter = 500,warnOnly = FALSE,minFactor = 1 / 2048
        ) #increase max. iterations
      ),silent = TRUE)


      if (class(fit.start) != "try-error") {
        #get parameters out of it
        parameters <- coef(fit.start)
        a1.start[i] <- as.vector((parameters["a1"]))
        b1.start[i] <- as.vector((parameters["b1"]))
        a2.start[i] <- as.vector((parameters["a2"]))
        b2.start[i] <- as.vector((parameters["b2"]))
      }
    }

    ##use obtained parameters for fit input
    a1.start <- median(a1.start, na.rm = TRUE)
    b1.start <- median(b1.start, na.rm = TRUE)
    a2.start <- median(a2.start, na.rm = TRUE)
    b2.start <- median(b2.start, na.rm = TRUE)

    ##perform final fitting
    fit <- try(minpack.lm::nlsLM(
      formula = fit.formulaEXPEXP,
      data = data,
      start = list(a1 = a1, b1 = b1, a2 = a2, b2 = b2),
      weights = fit.weights,
      trace = FALSE,
      algorithm = "LM",
      lower = if (fit.bounds) {
        c(0,0,0,0)
      }else{
        c(-Inf,-Inf,-Inf, -Inf)
      },
      control = minpack.lm::nls.lm.control(maxiter = 500)
    ), silent = TRUE
    )


    ##insert if for try-error
    if (class(fit)!="try-error") {

      #get parameters out of it
      parameters <- (coef(fit))
      b1 <- as.vector((parameters["b1"]))
      b2 <- as.vector((parameters["b2"]))
      a1 <- as.vector((parameters["a1"]))
      a2 <- as.vector((parameters["a2"]))

      ##set D0 values
      D01 <- round(b1,digits = 2)
      D02 <- round(b2,digits = 2)


      #problem: analytically it is not easy to calculate x, use uniroot
      f.unirootEXPEXP <- function(a1,a2,b1,b2,x,LnTn){fit.functionEXPEXP(a1,a2,b1,b2,x)-LnTn}

      temp.De <-  try(uniroot(f = f.unirootEXPEXP,
                              interval = c(0,max(xy$x)*1.5),
                              tol = 0.001,
                              a1 = a1,
                              a2 = a2,
                              b1 = b1,
                              b2 = b2,
                              LnTn = sample[1,2],
                              extendInt = "yes",
                              maxiter = 3000), silent = TRUE)


      if (class(temp.De) != "try-error") {
        De <- round(temp.De$root, digits = 2)
      }else{
        De <- NA
      }

      ##remove object
      rm(temp.De)

      #print D0 and De value values
      if(verbose){
        writeLines(paste0("\n [plot_GrowthCurve()] >> D01 = ",D01, " | D02 = ",D02, " | De = ", De))
      }


      ##Monte Carlo Simulation for error estimation
      #	--Fit many curves and calculate a new De +/- De_Error
      #	--take De_Error from the simulation
      # --comparison of De from the MC and original fitted De gives a value for quality

      #set variables
      var.b1<-vector(mode="numeric", length=NumberIterations.MC)
      var.b2<-vector(mode="numeric", length=NumberIterations.MC)
      var.a1<-vector(mode="numeric", length=NumberIterations.MC)
      var.a2<-vector(mode="numeric", length=NumberIterations.MC)

      ##terminal output fo MC
      if(verbose){
        cat("\n\t Run Monte Carlo loops for error estimation of the EXP+EXP fit\n")
      }

      ##progress bar
      if(txtProgressBar){
        pb<-txtProgressBar(min=0,max=NumberIterations.MC, initial=0, char="=", style=3)
      }

      #start Monto Carlo loops
      for (i in 1:NumberIterations.MC) {

        #update progress bar
        if(txtProgressBar) setTxtProgressBar(pb,i)

        data<-data.frame(x=xy$x,y=data.MC[,i])

        ##perform final fitting
        fit.MC <- try(minpack.lm::nlsLM(
          formula = fit.formulaEXPEXP,
          data = data,
          start = list(a1 = a1, b1 = b1, a2 = a2, b2 = b2),
          weights = fit.weights,
          trace = FALSE,
          algorithm = "LM",
          lower = if (fit.bounds) {
            c(0,0,0,0)
          }else{
            c(-Inf,-Inf,-Inf, -Inf)
          },
          control = minpack.lm::nls.lm.control(maxiter = 500)
        ), silent = TRUE
        )

        #get parameters out of it including error handling
        if (class(fit.MC)=="try-error") {

          x.natural[i]<-NA

        }else {
          parameters <- (coef(fit.MC))
          var.b1[i] <- as.vector((parameters["b1"]))
          var.b2[i] <- as.vector((parameters["b2"]))
          var.a1[i] <- as.vector((parameters["a1"]))
          var.a2[i] <- as.vector((parameters["a2"]))

          #problem: analytically it is not easy to calculate x, here an simple approximation is made

          temp.De.MC <-  try(uniroot(
            f = f.unirootEXPEXP,
            interval = c(0,max(xy$x) * 1.5),
            tol = 0.001,
            a1 = var.a1[i],
            a2 = var.a2[i],
            b1 = var.b1[i],
            b2 = var.b2[i],
            LnTn = data.MC.De[i]
          ), silent = TRUE)

          if (class(temp.De.MC) != "try-error") {
            x.natural[i] <- temp.De.MC$root
          }else{
            x.natural[i] <- NA
          }

        } #end if "try-error" MC simulation

      } #end for loop

    }else{

      #print message
      if(verbose) writeLines(paste0("[plot_GrowthCurve()] >> FITTING FAILED"))

    } #end if "try-error" Fit Method


    ##close
    if(txtProgressBar) if(exists("pb")){close(pb)}
    #===========================================================================
  } #End if Fit Method


  #Get De values from Monto Carlo simulation

  #calculate mean and sd (ignore NaN values)
  De.MonteCarlo<-round(mean(na.exclude(x.natural)),digits=2)

  #De.Error is Error of the whole De (ignore NaN values)
  De.Error <- sd(na.exclude(x.natural))

  ##choose format in dependency of the size of the error
  De.Error <- ifelse(De.Error <= 0.01,
                     format(De.Error, scientific = TRUE, digits = 2),
                     round(De.Error, digits = 2))


  # Formula creation --------------------------------------------------------
  if(!is(fit,"try-error") & !is.na(fit[1])){

    if(fit.method == "EXP") {
      f <- parse(text = paste0(format(coef(fit)[1], scientific = TRUE), " * (1 - exp( - ( x + ",
                               format(coef(fit)[3], scientific = TRUE), ") / ",
                               format(coef(fit)[2], scientific = TRUE), "))"))

    }

    if(fit.method == "EXP+LIN") {
      f <- parse(text = paste0(format(coef(fit)[1], scientific = TRUE), " * (1-exp(-(x+",
                               format(coef(fit)[3], scientific = TRUE), ") / ",
                               format(coef(fit)[2], scientific = TRUE), ")+(",
                               format(coef(fit)[4], scientific = TRUE), " * x))"))
    }

    if(fit.method == "EXP+EXP") {
      f <- parse(text = paste0(format(coef(fit)[1], scientific = TRUE), " * (1 - exp( -x / ",
                               format(coef(fit)[3], scientific = TRUE), ")) + ",
                               format(coef(fit)[2], scientific = TRUE), " * (1 - exp( -x / ",
                               format(coef(fit)[4], scientific = TRUE), "))"))
    }

    if(fit.method == "LIN" &  fit.force_through_origin) {
      f <- parse(text = paste0(format(fit.lm$coefficients[1], scientific = TRUE), " * x"))

    }

    if(fit.method == "LIN" &  !fit.force_through_origin) {
      f <- parse(text = paste0(format(fit.lm$coefficients[2], scientific = TRUE),
                               "* x + ", format(fit.lm$coefficients[1], scientific = TRUE)))

    }

    if(fit.method == "QDR"  &  fit.force_through_origin) {
      f <- parse(text = paste0(format(coef(fit)[1], scientific = TRUE), " * x ",
                               " + ", format(coef(fit)[2], scientific = TRUE), " * x^2"
      ))

    }

    if(fit.method == "QDR" & !fit.force_through_origin) {
      f <- parse(text = paste0(format(coef(fit)[1], scientific = TRUE),
                               " + ", format(coef(fit)[2], scientific = TRUE), " * x ",
                               " + ", format(coef(fit)[3], scientific = TRUE), " * x^2"
      ))

    }

  }else{

    f <- NA

  }

  ##============================================================================##
  # PLOTTING ---------------------------------------------------------------------
  ##============================================================================##

  ##5. Plotting if plotOutput==TRUE
  if(output.plot==TRUE) {


    ##cheat the R check
    x<-NULL; rm(x)

    #PAR	#open plot area
    if(output.plot== TRUE &
       output.plotExtended== TRUE &
       output.plotExtended.single == FALSE ){

      ####grep recent plot parameter for later reset
      par.default.complex <- par(no.readonly = TRUE)

      ##set new parameter
      layout(matrix(c(1,1,1,1,2,3), 3, 2, byrow=TRUE), respect=TRUE)
      par(cex=0.8*cex.global)

    }else{

      par.default.single <- par(no.readonly = TRUE)$cex
      par(cex=cex.global)

    }

    #PLOT		#Plot input values

    ##Make selection to support manual number of reg points input
    if(exists("fit.RegPointsReal")==TRUE){

      ##here the object sample has to be used otherwise the first regeneration point is not plotted.
      temp.xy.plot  <- sample[fit.RegPointsReal,]

    }else{

      temp.xy.plot  <- xy[1:fit.NumberRegPointsReal,]

    }

    plot(
      temp.xy.plot[, 1:2],
      ylim = ylim,
      xlim = xlim,
      pch = 19,
      xlab = xlab,
      ylab = ylab
    )

    #ADD HEADER
    title(main=main,line=3)

    #CURVE	#plot fitted curve
    if (fit.method=="EXP+LIN") {
      try(curve(a * (1 - exp(-(x + c) / b) + (g * x)), lwd = 1.5, add = TRUE))
    }
    else
      if (fit.method  ==  "LIN" &
          fit.force_through_origin)
      {
        curve(
          fit.lm$coefficients[1]  *  x, lwd  =  1.5,
          add  =  TRUE)
      }
    else if (fit.method == "LIN") {
      curve(fit.lm$coefficients[2] * x + fit.lm$coefficients[1],
            lwd = 1.5,
            add = TRUE)
    }
    else if (fit.method == "QDR" & fit.force_through_origin) {
      curve(coef(fit)[1] * x + coef(fit)[2] * x ^ 2,
            lwd = 1.5,
            add = TRUE)
    }
    else if (fit.method == "QDR") {
      curve(coef(fit)[1] + coef(fit)[2] * x + coef(fit)[3] * x ^ 2,
            lwd = 1.5,
            add = TRUE)
    }
    else if (fit.method == "EXP") {
      try(curve(fit.functionEXP(a, b, c, x), lwd = 1.5, add = TRUE))
    }
    else
      if (fit.method  ==  "EXP+EXP")
      {
        try(curve(fit.functionEXPEXP(a1, a2, b1, b2, x),
                  lwd  =  1.5,
                  add  =  TRUE))
      }

    ##POINTS	#Plot Reg0 and Repeated Points

    #Natural value
    points(sample[1,1:2], col = "red")
    segments(sample[1,1],sample[1,2]-sample[1,3],
             sample[1,1],sample[1,2]+sample[1,3], col = "red")

    #Repeated Point
    points(xy[which(duplicated(xy[,1])),1],xy[which(duplicated(xy[,1])),2],
           pch=2)

    #Reg Point 0
    points(xy[which(xy==0),1],xy[which(xy==0),2], pch=1, cex = 1.5*cex.global)

    ##ARROWS	#y-error Bars

    segments(xy$x,xy$y-y.Error,xy$x,xy$y+y.Error)

    ##LINES	#Insert Ln/Tn
    if(is.na(De)){

      lines(c(0,max(sample[,1])*2),c(sample[1,2],sample[1,2]), col="red", lty=2,lwd=1.25)

    }else{

      try(lines(c(0,De),c(sample[1,2],sample[1,2]), col="red", lty=2,lwd=1.25),silent=TRUE)

    }

    try(lines(c(De,De),c(0,sample[1,2]), col="red", lty=2, lwd=1.25),silent=TRUE)
    try(points(De,sample[1,2], col="red", pch=19),silent=TRUE)

    ## check/set mtext
    mtext <- if("mtext" %in% names(list(...))) {
      list(...)$mtext
    } else {
      substitute(D[e] == De,
                 list(De=paste(De,"\u00B1",De.Error, " | fit: ",fit.method)))
    }



    ##TEXT		#Insert fit and result
    try(mtext(side=3, mtext, line=0.5, cex=0.8*cex.global),silent=TRUE)

    #write error message in plot if De is NaN
    try(if (De=="NaN") {
      text(sample[2,1],0,"Error: De could not be calculated!",
           adj=c(0,0), cex=0.8, col="red")
    },silent=TRUE)

    ##LEGEND	#plot legend

    legend("topleft", c("REG points", "REG point repeated", "REG point 0"),
           pch=c(19,2,1), cex=0.8*cex.global, bty="n")

    ##plot only if wanted
    if(output.plot==TRUE & output.plotExtended==TRUE){

      ##HIST		#try to plot histogramm of De values from the Monte Carlo simulation

      if(output.plotExtended.single != TRUE){

        par(cex=0.7*cex.global)

      }


      ##(A) Calculate histogram data
      try(histogram <- hist(x.natural, plot = FALSE), silent = TRUE)

      #to avoid errors plot only if histogram exists
      if (exists("histogram") && length(histogram$counts) > 2) {

        ##calculate normal distribution curves for overlay
        norm.curve.x <- seq(min(x.natural, na.rm = TRUE),
                            max(x.natural, na.rm = TRUE),
                            length = 101)

        norm.curve.y <- dnorm(norm.curve.x,
                              mean=mean(x.natural, na.rm = TRUE),
                              sd=sd(x.natural, na.rm = TRUE))

        ##plot histogram
        histogram <- hist(x.natural,
                          xlab = xlab,
                          ylab = "Frequency",
                          main=expression(paste(D[e], " from MC simulation")),
                          freq=FALSE,
                          border = "white",
                          axes = FALSE,
                          ylim = c(0,max(norm.curve.y)),
                          sub =
                            paste("n = ", NumberIterations.MC, ", valid fits =", length(na.exclude(x.natural))),
                          col="grey")

        ##add axes
        axis(side = 1)
        axis(side = 2,
             at = seq(min(histogram$density),max(histogram$density), length = 5),
             labels = round(
               seq(min(histogram$counts),max(histogram$counts), length = 5),
               digits = 0))

        ##add norm curve
        lines(norm.curve.x, norm.curve.y, col = "red")

        ##add rug
        rug(x.natural)

        ##write De + Error from Monte Carlo simulation + write quality of error estimation
        try(mtext(side=3,substitute(D[e[MC]] == De,
                                    list(De=paste(De.MonteCarlo,"\u00B1",De.Error,
                                                  " | quality = ",round((1-abs(De-De.MonteCarlo)/De)*100,
                                                                        digits=1),"%"))),cex=0.6*cex.global),silent=TRUE)

      } else {

        plot(NA,NA,
             xlim=c(0,10),
             ylim=c(0,10),
             main=expression(paste(D[e], " from Monte Carlo simulation")))
        text(5,5,"not available")

      }#end ifelse


      ##PLOT		#PLOT test dose response curve if available if not plot not available
      #plot Tx/Tn value for sensitiviy change

      if ("TnTx" %in% colnames(sample)==TRUE) {

        plot(1:length(sample[,"TnTx"]),sample[1:(length(sample[,"TnTx"])),"TnTx"]/sample[1,"TnTx"],
             xlab="SAR cycle",
             ylab=expression(paste(T[n]/T[x])),
             main="Test dose response",
             type="o",
             pch=20,
        )

        ##LINES		#plot 1 line
        lines(c(1,length(sample[,"TnTx"])),c(1,1), lty=2, col="gray")
      } else {

        plot(NA,NA,xlim=c(0,10), ylim=c(0,10), main="Test dose response")
        text(5,5,"not available\n no TnTx column")
      }#end if else


      ## FUN by R Luminescence Team
      if(fun==TRUE){sTeve()}

      ##END lines
    }#endif::output.plotExtended


    ##reset only the parameter that have been changed!
    if(exists("par.default.single")){

      par(cex = par.default.single)
      rm(par.default.single)

    }

    if(exists("par.default.complex")){

      par(par.default.complex)
      rm(par.default.complex)

    }

  }#end if plotOutput

  ##RETURN - return De values and parameter

  output <- try(data.frame(
    De = De,
    De.Error = De.Error,
    D01 = D01,
    D02 = D02,
    De.MC = De.MonteCarlo,
    Fit = fit.method
  ),
  silent = TRUE
  )

  ##make RLum.Results object
  output.final <- set_RLum(
    class = "RLum.Results",
    data = list(
      De = output,
      De.MC = x.natural,
      Fit = fit,
      Formula = f,
      call = sys.call()
    )
  )
  invisible(output.final)

}
