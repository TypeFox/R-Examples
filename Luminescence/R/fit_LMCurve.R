#' Nonlinear Least Squares Fit for LM-OSL curves
#'
#' The function determines weighted nonlinear least-squares estimates of the
#' component parameters of an LM-OSL curve (Bulur 1996) for a given number of
#' components and returns various component parameters. The fitting procedure
#' uses the function \code{\link{nls}} with the \code{port} algorithm.
#'
#' \bold{Fitting function}\cr\cr The function for the fitting has the general
#' form: \deqn{y = (exp(0.5)*Im_1*x/xm_1)*exp(-x^2/(2*xm_1^2)) + ,\ldots, +
#' exp(0.5)*Im_i*x/xm_i)*exp(-x^2/(2*xm_i^2))} where \eqn{1 < i < 8}\cr This
#' function and the equations for the conversion to b (detrapping probability)
#' and n0 (proportional to initially trapped charge) have been taken from Kitis
#' et al. (2008): \deqn{xm_i=\sqrt{max(t)/b_i}} \deqn{Im_i=exp(-0.5)n0/xm_i}\cr
#' \bold{Background subtraction}\cr\cr Three methods for background subtraction
#' are provided for a given background signal (\code{values.bg}).\cr
#' \code{polynomial}: default method. A polynomial function is fitted using
#' \link{glm} and the resulting function is used for background subtraction:
#' \deqn{y = a*x^4 + b*x^3 + c*x^2 + d*x + e}\cr \code{linear}: a linear
#' function is fitted using \link{glm} and the resulting function is used for
#' background subtraction: \deqn{y = a*x + b}\cr \code{channel}: the measured
#' background signal is subtracted channelwise from the measured signal.\cr\cr
#' \bold{Start values}\cr
#'
#' The choice of the initial parameters for the \code{nls}-fitting is a crucial
#' point and the fitting procedure may mainly fail due to ill chosen start
#' parameters. Here, three options are provided:\cr\cr \bold{(a)} If no start
#' values (\code{start_values}) are provided by the user, a cheap guess is made
#' by using the detrapping values found by Jain et al. (2003) for quartz for a
#' maximum of 7 components. Based on these values, the pseudo start parameters
#' xm and Im are recalculated for the given data set. In all cases, the fitting
#' starts with the ultra-fast component and (depending on \code{n.components})
#' steps through the following values. If no fit could be achieved, an error
#' plot (for \code{plot = TRUE}) with the pseudo curve (based on the
#' pseudo start parameters) is provided. This may give the opportunity to
#' identify appropriate start parameters visually.\cr\cr \bold{(b)} If start
#' values are provided, the function works like a simple \code{\link{nls}}
#' fitting approach.\cr\cr \bold{(c)} If no start parameters are provided and
#' the option \code{fit.advanced = TRUE} is chosen, an advanced start paramter
#' estimation is applied using a stochastical attempt. Therefore, the
#' recalculated start parameters \bold{(a)} are used to construct a normal
#' distribution. The start parameters are then sampled randomly from this
#' distribution. A maximum of 100 attempts will be made. \bold{Note:} This
#' process may be time consuming. \cr\cr \bold{Goodness of fit}\cr\cr The
#' goodness of the fit is given by a pseudoR^2 value (pseudo coefficient of
#' determination). According to Lave (1970), the value is calculated as:
#' \deqn{pseudoR^2 = 1 - RSS/TSS} where \eqn{RSS = Residual~Sum~of~Squares} \cr
#' and \eqn{TSS = Total~Sum~of~Squares}\cr\cr \bold{Error of fitted component
#' parameters}\cr\cr The 1-sigma error for the components is calculated using
#' the function \link{confint}. Due to considerable calculation time, this
#' option is deactived by default. In addition, the error for the components
#' can be estimated by using internal R functions like \link{summary}. See the
#' \link{nls} help page for more information.\cr \emph{For more details on the
#' nonlinear regression in R, see Ritz & Streibig (2008).}
#'
#' @param values \code{\linkS4class{RLum.Data.Curve}} or \link{data.frame}
#' (\bold{required}): x,y data of measured values (time and counts). See
#' examples.
#'
#' @param values.bg \code{\linkS4class{RLum.Data.Curve}} or \link{data.frame}
#' (optional): x,y data of measured values (time and counts) for background
#' subtraction.
#'
#' @param n.components \link{integer} (with default): fixed number of
#' components that are to be recognised during fitting (min = 1, max = 7).
#'
#' @param start_values \link{data.frame} (optional): start parameters for lm
#' and xm data for the fit. If no start values are given, an automatic start
#' value estimation is attempted (see details).
#'
#' @param input.dataType \link{character} (with default): alter the plot output
#' depending on the input data: "LM" or "pLM" (pseudo-LM). See: \link{CW2pLM}
#'
#' @param fit.method \code{\link{character}} (with default): select fit method,
#' allowed values: \code{'port'} and \code{'LM'}. \code{'port'} uses the 'port'
#' routine usint the funtion \code{\link{nls}} \code{'LM'} utilises the
#' function \code{nlsLM} from the package \code{minpack.lm} and with that the
#' Levenberg-Marquardt algorithm.
#'
#' @param sample_code \link{character} (optional): sample code used for the
#' plot and the optional output table (mtext).
#'
#' @param sample_ID \link{character} (optional): additional identifier used as
#' column header for the table output.
#'
#' @param LED.power \link{numeric} (with default): LED power (max.) used for
#' intensity ramping in mW/cm^2. \bold{Note:} This value is used for the
#' calculation of the absolute photoionisation cross section.
#'
#' @param LED.wavelength \link{numeric} (with default): LED wavelength in nm
#' used for stimulation. \bold{Note:} This value is used for the calculation of
#' the absolute photoionisation cross section.
#'
#' @param fit.trace \link{logical} (with default): traces the fitting process
#' on the terminal.
#'
#' @param fit.advanced \link{logical} (with default): enables advanced fitting
#' attempt for automatic start parameter recognition. Works only if no start
#' parameters are provided. \bold{Note:} It may take a while and it is not
#' compatible with \code{fit.method = "LM"}.
#'
#' @param fit.calcError \link{logical} (with default): calculate 1-sigma error
#' range of components using \link{confint}.
#'
#' @param bg.subtraction \link{character} (with default): specifies method for
#' background subtraction (\code{polynomial}, \code{linear}, \code{channel},
#' see Details). \bold{Note:} requires input for \code{values.bg}.
#'
#' @param verbose \link{logical} (with default): terminal output with
#' fitting results.
#'
#' @param plot \link{logical} (with default): returns a plot of the
#' fitted curves.
#'
#' @param plot.BG \link{logical} (with default): returns a plot of the
#' background values with the fit used for the background subtraction.
#'
#' @param \dots Further arguments that may be passed to the plot output, e.g.
#' \code{xlab}, \code{xlab}, \code{main}, \code{log}.
#'
#' @return
#' Various types of plots are returned. For details see above.\cr
#' Furthermore an \code{RLum.Results} object is returned with the following structure:\cr
#'
#' data:\cr
#' .. $fit : \code{nls} (nls object)\cr
#' .. $output.table : \code{data.frame} with fitting results\cr
#' .. $component.contribution.matrix : \code{list} component distribution matrix\cr
#' .. $call : \code{call} the original function call
#'
#' Matrix structure for the distribution matrix:\cr
#'
#' Column 1 and 2: time and \code{rev(time)} values\cr
#' Additional columns are used for the components, two for each component,
#' containing I0 and n0. The last columns \code{cont.} provide information on
#' the relative component contribution for each time interval including the row
#' sum for this values.
#'
#' @note The pseudo-R^2 may not be the best parameter to describe the goodness
#' of the fit. The trade off between the \code{n.components} and the pseudo-R^2
#' value currently remains unconsidered. \cr
#'
#' The function \bold{does not} ensure that the fitting procedure has reached a
#' global minimum rather than a local minimum! In any case of doubt, the use of
#' manual start values is highly recommended.
#'
#' @section Function version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\link{fit_CWCurve}}, \code{\link{plot}}, \code{\link{nls}},
#' \code{\link[minpack.lm]{nlsLM}}, \code{\link{get_RLum}}
#'
#' @references Bulur, E., 1996. An Alternative Technique For Optically
#' Stimulated Luminescence (OSL) Experiment. Radiation Measurements, 26, 5,
#' 701-709.
#'
#' Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of
#' blue-light stimulated luminescence components in different quartz samples:
#' implications for dose measurement. Radiation Measurements, 37 (4-5),
#' 441-449.
#'
#' Kitis, G. & Pagonis, V., 2008. Computerized curve deconvolution analysis for
#' LM-OSL. Radiation Measurements, 43, 737-741.
#'
#' Lave, C.A.T., 1970. The Demand for Urban Mass Transportation. The Review of
#' Economics and Statistics, 52 (3), 320-323.
#'
#' Ritz, C. & Streibig, J.C., 2008. Nonlinear Regression with R. R. Gentleman,
#' K. Hornik, & G. Parmigiani, eds., Springer, p. 150.
#'
#' @keywords dplot models
#'
#' @examples
#'
#'
#' ##(1) fit LM data without background subtraction
#' data(ExampleData.FittingLM, envir = environment())
#' fit_LMCurve(values = values.curve, n.components = 3, log = "x")
#'
#' ##(2) fit LM data with background subtraction and export as JPEG
#' ## -alter file path for your preferred system
#' ##jpeg(file = "~/Desktop/Fit_Output\%03d.jpg", quality = 100,
#' ## height = 3000, width = 3000, res = 300)
#' data(ExampleData.FittingLM, envir = environment())
#' fit_LMCurve(values = values.curve, values.bg = values.curveBG,
#'             n.components = 2, log = "x", plot.BG = TRUE)
#' ##dev.off()
#'
#' ##(3) fit LM data with manual start parameters
#' data(ExampleData.FittingLM, envir = environment())
#' fit_LMCurve(values = values.curve,
#'             values.bg = values.curveBG,
#'             n.components = 3,
#'             log = "x",
#'             start_values = data.frame(Im = c(170,25,400), xm = c(56,200,1500)))
#'
#' @export
fit_LMCurve<- function(
  values,
  values.bg,
  n.components = 3,
  start_values,
  input.dataType = "LM",
  fit.method = "port",
  sample_code = "",
  sample_ID = "",
  LED.power = 36,
  LED.wavelength = 470,
  fit.trace = FALSE,
  fit.advanced = FALSE,
  fit.calcError = FALSE,
  bg.subtraction = "polynomial",
  verbose = TRUE,
  plot = TRUE,
  plot.BG = FALSE,
  ...
){

  # (0) Integrity checks -------------------------------------------------------

  ##(1) data.frame or RLum.Data.Curve object?
  if(is(values, "data.frame") == FALSE & is(values, "RLum.Data.Curve") == FALSE){

    stop("[fit_LMCurve()] 'values' object has to be of type
         'data.frame' or 'RLum.Data.Curve'!")

  }else{

    if(is(values, "RLum.Data.Curve") == TRUE && (
      values@recordType!="RBR" & values@recordType!="LM-OSL")){

      stop("[fit_LMCurve()] recordType should be 'RBR' or 'LM-OSL'!
           Consider as(object,'data.frame') if you had used the pseudo transformation functions.")

    }else if(is(values, "RLum.Data.Curve") == TRUE){

      values <- as(values,"data.frame")

    }
  }

  ##(2) data.frame or RLum.Data.Curve object?
  if(missing(values.bg)==FALSE){

    if(is(values.bg, "data.frame") == FALSE & is(values.bg,
                                                 "RLum.Data.Curve") == FALSE){

      stop("[fit_LMCurve()] 'values.bg' object has to be of type 'data.frame' or 'RLum.Data.Curve'!")

    }else{

      if(is(values, "RLum.Data.Curve") == TRUE && values@recordType!="RBR"){

        stop("[fit_LMCurve()] recordType should be 'RBR'!")

      }else if(is(values.bg, "RLum.Data.Curve") == TRUE){

        values.bg <- as(values.bg,"data.frame")

      }
    }
  }

  ## Set plot format parameters -----------------------------------------------
  extraArgs <- list(...) # read out additional arguments list

  log       <- if("log" %in% names(extraArgs)) {extraArgs$log}
  else {""}

  xlim      <- if("xlim" %in% names(extraArgs)) {extraArgs$xlim}
  else {c(min(values[,1]),max(values[,1]))}

  ylim      <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim}
  else {

    if(input.dataType=="pLM"){
      c(0,max(values[,2]*1.1))
    }else{
      c(min(values[,2]),max(values[,2]*1.1))
    }

  }

  xlab      <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab}
  else {

    if(input.dataType=="LM"){"Time [s]"}else{"u [s]"}

  }

  ylab     <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab}
  else {

    if(input.dataType=="LM"){
      paste("LM-OSL [cts/",round(max(values[,1])/length(values[,1]),digits=2)," s]",sep="")
    }else{"pLM-OSL [a.u.]"}
  }


  main      <- if("main" %in% names(extraArgs)) {extraArgs$main}
  else {"Default"}

  cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex}
  else {0.8}


  fun       <- if("fun" %in% names(extraArgs)) {extraArgs$fun} else {FALSE}


  ##============================================================================##
  ##  BACKGROUND SUBTRACTION
  ##============================================================================##

  #   ##perform background subtraction if background LM measurment exists

  if(missing(values.bg)==FALSE){

    #set graphical parameters
    par.default <- par(mfrow=c(1,1), cex=1.5*cex)

    ##check if length of bg and signal is consistent
    if(length(values[,2])!=length(values.bg[,2])){stop("[fit_LMCurve] Length of values and values.bg differs!")}

    if(bg.subtraction=="polynomial"){

      #fit polynom function to background
      glm.fit<-glm(values.bg[,2] ~ values.bg[,1]+I(values.bg[,1]^2)+I(values.bg[,1]^3))
      glm.coef<-coef(glm.fit)

      #subtract background with fitted function
      values[,2]<-values[,2]-
        (glm.coef[4]*values[,1]^3+glm.coef[3]*values[,1]^2+glm.coef[2]*values[,1]+glm.coef[1])
      writeLines("[fit_LMCurve] >> Background subtracted (method=\"polynomial\")!")

      ##plot Background measurement if needed
      if(plot.BG==TRUE){

        plot(values.bg, ylab="LM-OSL [a.u.]", xlab="Time [s]", main="Background")
        curve((glm.coef[4]*x^3+glm.coef[3]*x^2+glm.coef[2]*x+glm.coef[1]),add=TRUE,col="red",lwd=2)
        text(0,max(values.bg[,2]),paste("y = ", round(glm.coef[4],digits=2),
                                        "*x^3+",
                                        round(glm.coef[3],digits=2),
                                        "*x^2+",
                                        round(glm.coef[2],digits=2),
                                        "*x+",
                                        round(glm.coef[1],digits=2),
                                        sep=""),pos=4)
        mtext(side=3,sample_code,cex=.8*cex)
      }

    }else if(bg.subtraction=="linear"){

      #fit linear function to background
      glm.fit<-glm(values.bg[,2] ~ values.bg[,1])
      glm.coef<-coef(glm.fit)

      ##substract bg
      values[,2]<-values[,2]-(glm.coef[2]*values[,1]+glm.coef[1])
      writeLines("[fit_LMCurve.R] >> Background subtracted (method=\"linear\")!")

      ##plot Background measurement if needed
      if(plot.BG==TRUE){

        plot(values.bg, ylab="LM-OSL [a.u.]", xlab="Time [s]", main="Background")
        curve((glm.coef[2]*x+glm.coef[1]),add=TRUE,col="red",lwd=1.5)
        text(0,max(values.bg[,2]),paste("y = ",
                                        round(glm.coef[2],digits=2),
                                        "*x+",
                                        round(glm.coef[1],digits=2),
                                        sep=""),pos=4)
        mtext(side=3,sample_code,cex=.8*cex)

      }#endif::plot BG

    }else if(bg.subtraction=="channel"){

      values[,2]<-values[,2]-values.bg[,2]
      writeLines("[fit_LMCurve.R] >> Background subtracted (method=\"channel\")!")

      if(plot.BG==TRUE){

        plot(values.bg, ylab="LM-OSL [a.u.]", xlab="Time [s]", main="Background")
        mtext(side=3,sample_code,cex=.8*cex)
      }

    }else{stop("Error: Invalid method for background subtraction")}

    ##reset par values
    par(par.default)
    rm(par.default)
  }


  ##============================================================================##
  ##  FITTING
  ##============================================================================##

  ##------------------------------------------------------------------------##
  ##set function for fit equation (according Kitis and Pagonis, 2008)
  ##////equation used for fitting////(start)
  fit.equation<-function(Im.i,xm.i){
    equation<-parse(
      text=paste("exp(0.5)*Im[",Im.i,"]*(values[,1]/xm[",xm.i,"])*exp(-values[,1]^2/(2*xm[",xm.i,"]^2))",
                 collapse="+",sep=""))
    return(equation)
  }
  ##////equation used for fitting///(end)
  ##------------------------------------------------------------------------##

  ##set formula elements for fitting functions
  ## the upper two funtions should be removed ... but chances are needed ... TODO
  ##////equation used for fitting////(start)
  fit.formula <- function(n.components){

    Im <- paste0("Im.",1:n.components)
    xm <- paste0("xm.",1:n.components)

    as.formula(paste0("y ~ ", paste("(exp(0.5) * ", Im, "* x/", xm, ") * exp(-x^2/(2 *",xm,"^2))", collapse=" + ")))

  }
  ##////equation used for fitting///(end)

  ##------------------------------------------------------------------------##
  ##automatic start parameter estimation

  ##set fit function
  fit.function<-fit.equation(Im.i=1:n.components,xm.i=1:n.components)

  if(missing(start_values)){

    ##set b (detrapping) values for a 7-component function taken from Jain et al. (2003)
    b.pseudo<-c(32,2.5,0.65,0.15,0.025,0.0025,0.00030)

    ##calculate xm parameters from values set based on the pseudo curves
    xm.pseudo<-sqrt(max(values[,1])/b.pseudo)

    ##the Im values obtaind by calculating residuals
    xm.residual<-sapply(1:length(b.pseudo),function(x){abs(values[,1]-xm.pseudo[x])})
    xm.residual<-cbind(xm.residual,values[,1])
    Im.pseudo<-sapply(1:length(xm.pseudo),function(x){
      min(xm.residual[which(xm.residual[,x]==min(xm.residual[,x])),8])#8 is time index
    })

    ##set additional variables
    b.pseudo_start<-1
    b.pseudo_end<-0
    fit.trigger<-FALSE

    while(fit.trigger==FALSE){


      xm<-xm.pseudo[b.pseudo_start:(n.components+b.pseudo_end)]
      Im<-Im.pseudo[b.pseudo_start:(n.components+b.pseudo_end)]

      if(fit.advanced){
        ##---------------------------------------------------------------##
        ##MC for fitting parameter
        ##make the fitting more stable by small variations of the parameters

        ##sample input parameters values from a normal distribution
        xm.MC<-sapply(1:length(xm),function(x){
          xm.MC<-sample(rnorm(30,mean=xm[x],sd=xm[x]/10), replace=TRUE)
        })

        Im.MC<-sapply(1:length(xm),function(x){
          Im.MC<-sample(rnorm(30,mean=Im[x],sd=Im[x]/10), replace=TRUE)

        })
        ##---------------------------------------------------------------##

        for(i in 1:length(xm.MC[,1])){

          ##NLS          ##try fit
          fit<-try(nls(y~eval(fit.function),
                       trace=fit.trace,
                       data=data.frame(x=values[,1],y=values[,2]),
                       algorithm="port",
                       start=list(Im=Im.MC[i,],xm=xm.MC[i,]),#end start values input
                       nls.control(
                         maxiter=500
                       ),#end nls control
                       lower=c(xm=min(values[,1]),Im=0),
                       upper=c(xm=max(values[,1]),Im=max(values[,2]*1.1))
          ),# nls
          silent=TRUE)# end try
          ##graphical output
          if(i==1){cat(paste("[fit_LMCurve()] >> advanced fitting attempt (#",
                             b.pseudo_start,"): ",sep=""))}
          cat("*")

          if(inherits(fit,"try-error") == FALSE){break}
        }#end::forloop

        cat("\n")

      }else{


        if(fit.method == "port") {
          fit <- try(nls(
            y ~ eval(fit.function),
            trace = fit.trace,
            data = data.frame(x = values[,1],y = values[,2]),
            algorithm = "port",
            start = list(Im = Im,xm = xm),#end start values input
            nls.control(maxiter = 500),#end nls control
            lower = c(xm = 0,Im = 0)
          ),# nls
          silent = TRUE)
          # end try

        }else if (fit.method == "LM") {
          ##re-name for method == "LM"
          names(Im) <- paste0("Im.", 1:n.components)
          names(xm) <- paste0("xm.", 1:n.components)
          start.list <- c(as.list(Im), as.list(xm))
          lower <-
            sapply(start.list, function(x) {
              start.list[[x]] <- 0
            })

          fit <- try(minpack.lm::nlsLM(
            fit.formula(n.components),
            data = data.frame(x = values[,1],
                              y = values[,2]),
            start = start.list,
            lower = lower,
            trace = fit.trace,
            control = minpack.lm::nls.lm.control(maxiter = 500)
          ), silent = TRUE)

        }else{

          stop("[fit_LMCurve()] unknow method for 'fit.method'")

        }


      }#endifelse::fit.advanced


      if(inherits(fit,"try-error")==FALSE){fit.trigger<-TRUE}
      else{

        if((n.components+b.pseudo_end)==7){fit.trigger<-TRUE
        }else{
          b.pseudo_start<-b.pseudo_start+1
          b.pseudo_end<-b.pseudo_end+1
        }#endif::maximum loops
      }#endif::try-error
    }#end:whileloop fit trigger

  }else{#endif::missing start values
    ##------------------------------------------------------------------------##

    fit<-try(nls(y~eval(fit.function),
                 trace=fit.trace, data.frame(x=values[,1],y=values[,2]),
                 algorithm="port", start=list(Im=start_values[,1],xm=start_values[,2]),#end start values input
                 nls.control(maxiter=500),
                 lower=c(xm=0,Im=0),
                 #upper=c(xm=max(x),Im=max(y)*1.1)# set lower boundaries for components
    )# nls
    )# end try
  }#endif::startparameter

  ##------------------------------------------------------------------------##

  ##grep parameters
  if(inherits(fit,"try-error")==FALSE){
    parameters<-coef(fit)

    ##write parameters in vectors and order parameters
    Im<-parameters[1:(length(parameters)/2)]
    Im.names <- names(Im)
    xm<-parameters[(1+(length(parameters)/2)):length(parameters)]
    xm.names <- names(xm)

    ##order parameters
    o <- order(xm)
    xm <- xm[o]
    names(xm) <- xm.names
    Im <- Im[o]
    names(Im) <- Im.names

    if (verbose){
      ##print rough fitting information - use the nls() control for more information
      writeLines("\n[fit_LMCurve()]")
      writeLines(paste("\nFitting was done using a ",n.components, "-component function:\n",sep=""))

      ##print parameters
      print(c(xm, Im))

      #print some additional information
      writeLines("\n(equation used for fitting according Kitis & Pagonis, 2008)")
    }#end if

    ##============================================================================##
    ##  Additional Calculations
    ##============================================================================##

    ##calculate stimulation intensity Schmidt (2008)

    ##Energy - E = h*v
    h<-6.62606957e-34 #in W*s^2 - Planck constant
    ny<-299792458/(LED.wavelength/10^9) #frequency of light
    E<-h*ny

    ##transform LED.power in W/cm^2
    LED.power<-LED.power/1000

    stimulation_intensity<-LED.power/E


    ##calculate b and n from the equation of Bulur(1996) to compare results
    ##Using Equation 5 and 6 from Kitis (2008)
    b<-as.vector(max(values[,1])/xm^2) #detrapping probability
    n0<-as.vector((Im/exp(-0.5))*xm)


    ##CALCULATE 1- sigma CONFIDENCE INTERVALL
    ##------------------------------------------------------------------------##
    b.error<-rep(NA, n.components)
    n0.error<-rep(NA, n.components)

    if(fit.calcError==TRUE){
      ##option for confidence interval
      values.confint<-confint(fit, level=0.68)
      Im.confint<-values.confint[1:(length(values.confint[,1])/2),]
      xm.confint<-values.confint[((length(values.confint[,1])/2)+1):length(values.confint[,1]),]

      ##error calculation
      b.error<-as.vector(abs((max(values[,1])/xm.confint[,1]^2)-(max(values[,1])/xm.confint[,2]^2)))
      n0.error<-as.vector(abs(((Im.confint[,1]/exp(-0.5))*xm.confint[,1]) - ((Im.confint[,2]/exp(-0.5))*xm.confint[,2])))
    }
    ##------------------------------------------------------------------------##


    ##calculate photoionisation cross section and print on terminal
    ##using EQ (5) in Kitis
    cs<-as.vector((max(values[,1])/xm^2)/stimulation_intensity)
    rel_cs<-round(cs/cs[1],digits=4)

    ##coefficient of determination after law
    RSS <- sum(residuals(fit)^2) #residual sum of squares
    TSS <- sum((values[,2] - mean(values[,2]))^2) #total sum of squares
    pR<-round(1-RSS/TSS,digits=4)

    ##============================================================================##
    ## COMPONENT TO SUM CONTRIBUTION MATRIX
    ##============================================================================##

    ##+++++++++++++++++++++++++++++++
    ##set matrix
    ##set polygon matrix for optional plot output
    component.contribution.matrix <- matrix(NA,
                                            nrow = length(values[,1]),
                                            ncol = (2*length(xm)) + 2)

    ##set x-values
    component.contribution.matrix[,1] <- values[,1]
    component.contribution.matrix[,2] <- rev(values[,1])

    ##+++++++++++++++++++++++++++++++
    ##set 1st polygon
    ##1st polygon (calculation)
    y.contribution_first <- (exp(0.5)*Im[1]*values[,1]/
                               xm[1]*exp(-values[,1]^2/(2*xm[1]^2))/
                               (eval(fit.function))*100)

    ##avoid NaN values (might happen with synthetic curves)
    y.contribution_first[is.nan(y.contribution_first)==TRUE] <- 0

    ##set values in matrix
    component.contribution.matrix[,3] <- 100
    component.contribution.matrix[,4] <- 100-rev(y.contribution_first)

    ##+++++++++++++++++++++++++++++++
    ##set polygons in between
    ##polygons in between (calculate and plot)
    if (length(xm)>2){

      y.contribution_prev <- y.contribution_first
      i<-2

      ##matrix stepping
      k <- seq(3, ncol(component.contribution.matrix), by=2)

      while (i<=length(xm)-1) {
        y.contribution_next<-(exp(0.5)*Im[i]*values[,1]/
                                xm[i]*exp(-values[,1]^2/(2*xm[i]^2))/
                                (eval(fit.function))*100)

        ##avoid NaN values
        y.contribution_next[is.nan(y.contribution_next)==TRUE] <- 0

        ##set values in matrix
        component.contribution.matrix[, k[i]] <- 100-y.contribution_prev
        component.contribution.matrix[, k[i]+1] <- rev(100-y.contribution_prev-
                                                         y.contribution_next)

        y.contribution_prev <- y.contribution_prev + y.contribution_next

        i<-i+1
      }#end while loop
    }#end if

    ##+++++++++++++++++++++++++++++++
    ##set last polygon

    ##last polygon (calculation)
    y.contribution_last<-(exp(0.5)*Im[length(xm)]*values[,1]/
                            xm[length(xm)]*exp(-values[,1]^2/
                                                 (2*xm[length(xm)]^2))/
                            (eval(fit.function))*100)

    ##avoid NaN values
    y.contribution_last[is.nan(y.contribution_last)==TRUE]<-0

    component.contribution.matrix[,((2*length(xm))+1)] <- y.contribution_last
    component.contribution.matrix[,((2*length(xm))+2)] <- 0

    ##change names of matrix to make more easy to understand
    component.contribution.matrix.names <- c("x", "rev.x",
                                             paste(c("y.c","rev.y.c"),rep(1:n.components,each=2), sep=""))


    ##calculate area for each component, for each time interval
    component.contribution.matrix.area <- sapply(
      seq(3,ncol(component.contribution.matrix),by=2),
      function(x){

        matrixStats::rowDiffs(cbind(rev(component.contribution.matrix[,(x+1)]),
                       component.contribution.matrix[,x]))

      })

    ##append to existing matrix
    component.contribution.matrix <- cbind(
      component.contribution.matrix,
      component.contribution.matrix.area,
      rowSums(component.contribution.matrix.area)
    )

    ##set final column names
    colnames(component.contribution.matrix) <- c(
      component.contribution.matrix.names,
      paste(c("cont.c"),rep(1:n.components,each=1), sep=""),
      "cont.sum")

    ##============================================================================##
    ##  Terminal Output (advanced)
    ##============================================================================##
    if (verbose){
      ##write fill lines
      writeLines("------------------------------------------------------------------------------")
      writeLines("(1) Corresponding values according the equation in Bulur, 1996 for b and n0:\n")
      for (i in 1:length(b)){
        writeLines(paste("b",i," = ",format(b[i],scientific=TRUE)," +/- ",format(b.error[i],scientific=TRUE),sep=""))
        writeLines(paste("n0",i," = ",format(n0[i],scientific=TRUE)," +/- ",format(n0.error[i],scientific=TRUE),"\n",sep=""))
      }#end for loop

      ##write photoionisation cross section on terminal
      for (i in 1:length(cs)){
        writeLines(paste("cs from component.",i," = ",format(cs[i],scientific=TRUE, digits=4), " cm^2",
                         "\t >> relative: ",round(cs[i]/cs[1],digits=4),sep=""))

      }#end for loop

      writeLines(paste(
        "\n(stimulation intensity value used for calculation: ",format(stimulation_intensity,scientific=TRUE)," 1/s 1/cm^2)",sep=""))
      writeLines("(errors quoted as 1-sigma uncertainties)")
      writeLines("------------------------------------------------------------------------------\n")

      #sum of squares
      writeLines(paste("pseudo-R^2 = ",pR,sep=""))
    }#end if

    ##============================================================================##
    ##  COMPOSE RETURN VALUES (data.frame)
    ##============================================================================##

    ##write output table if values exists
    if (exists("fit")){

      ##set data.frame for a max value of 7 components
      output.table <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,
                                 NA,NA,NA,NA,NA,NA,NA,NA,
                                 NA,NA,NA,NA,NA,NA,NA,NA,
                                 NA,NA,NA,NA,NA,NA,NA,NA,
                                 NA,NA,NA,NA,NA,NA,NA,NA,
                                 NA,NA,NA,NA,NA,NA,NA,NA,
                                 NA,NA,NA,NA,NA,NA,NA,NA)

      output.tableColNames<-c("Im1","xm1",
                              "b1","b1.error","n01","n01.error",
                              "cs1","rel_cs1",
                              "Im2","xm2",
                              "b2","b2.error","n02","n02.error",
                              "cs2","rel_cs2",
                              "Im3","xm3",
                              "b3","b3.error","n03","n03.error",
                              "cs3","rel_cs3",
                              "Im4","xm4",
                              "b4","b4.error","n04","n04.error",
                              "cs4","rel_cs4",
                              "Im5","xm5",
                              "b5","b5.error","n05","n05.error",
                              "cs5","rel_cs5",
                              "Im6","xm6",
                              "b6","b6.error","n06","n06.error",
                              "cs6","rel_cs6",
                              "Im7","xm7",
                              "b7","b7.error","n07","n07.error",
                              "cs7","rel_cs7")


      ##write components in output table
      i<-0
      k<-1
      while(i<=n.components*8){
        output.table[1,i+1]<-Im[k]
        output.table[1,i+2]<-xm[k]
        output.table[1,i+3]<-b[k]
        output.table[1,i+4]<-b.error[k]
        output.table[1,i+5]<-n0[k]
        output.table[1,i+6]<-n0.error[k]
        output.table[1,i+7]<-cs[k]
        output.table[1,i+8]<-rel_cs[k]
        i<-i+8
        k<-k+1
      }

      ##add pR and n.components
      output.table<-cbind(sample_ID,sample_code,n.components,output.table,pR)

      ###alter column names
      colnames(output.table)<-c("ID","sample_code","n.components",output.tableColNames,"pseudo-R^2")

      ##----------------------------------------------------------------------------##
    }#endif::exists fit
  }else{

    output.table <- NA
    component.contribution.matrix <- NA
    writeLines("[fit_LMCurve] Fitting Error: Plot without fit produced!")

  }
  ##============================================================================##
  ##  PLOTTING
  ##============================================================================##
  if(plot){

    ##cheat the R check routine
    x <- NULL; rm(x)

    ##grep package colour gallery
    col <- get("col", pos = .LuminescenceEnv)

    ##set plot frame
    par.default <- par(no.readonly = TRUE)

    layout(matrix(c(1,2,3),3,1,byrow=TRUE),c(1.6,1,1), c(1,0.3,0.4),TRUE)
    par(oma=c(1,1,1,1),mar=c(0,4,3,0), cex=cex)

    ##==uppper plot==##
    ##open plot area
    plot(NA,NA,
         xlim=xlim,
         ylim=ylim,
         xlab="",
         xaxt="n",
         main=main,
         log=log,
         ylab=ylab
    )#endplot

    mtext(side=3,sample_code,cex=0.8*cex)

    ##plotting measured signal
    points(values[,1],values[,2],pch=20, col=rgb(0.4,0.4,0.4,0.5))

    ##==pseudo curve==##------------------------------------------------------#

    ##curve for used pseudo values
    if(inherits(fit,"try-error")==TRUE & missing(start_values)==TRUE){
      fit.function<-fit.equation(Im.i=1:n.components,xm.i=1:n.components)
      Im<-Im.pseudo[1:n.components]
      xm<-xm.pseudo[1:n.components]

      ##draw pseudo curve
      lines(values[,1],eval(fit.function), lwd=2, col="red", lty=2)

      axis(side=1)
      mtext(side=1,xlab, cex=.9*cex,line=2)

      mtext(side=4,paste(n.components, " component pseduo function is shown",sep=""),cex=0.7, col="blue")

      ##draw information text on plot
      text(min(values[,1]),max(values[,2]),"FITTING ERROR!",pos=4)

      ##additional legend
      legend("topright",c("pseudo sum function"),lty=2,lwd=2,col="red",bty="n")

    }
    ##==pseudo curve==##------------------------------------------------------##

    ##plot sum function
    if(inherits(fit,"try-error")==FALSE){
      lines(values[,1],eval(fit.function), lwd=2, col="black")
      legend.caption<-"sum curve"
      curve.col<-1

      ##plot signal curves

      ##plot curve for additional parameters
      for (i in 1:length(xm)) {
        curve(exp(0.5)*Im[i]*x/xm[i]*exp(-x^2/(2*xm[i]^2)),col=col[i+1], lwd=2,add=TRUE)
        legend.caption<-c(legend.caption,paste("component ",i,sep=""))
        curve.col<-c(curve.col,i+1)
      }
      ##plot legend
      legend(if(log=="x"| log=="xy"){
        if(input.dataType=="pLM"){"topright"}else{"topleft"}}else{"topright"},
        legend.caption,lty=1,lwd=2,col=col[curve.col], bty="n")


      ##==lower plot==##
      ##plot residuals
      par(mar=c(4.2,4,0,0))
      plot(values[,1],residuals(fit),
           xlim=xlim,
           xlab=xlab,
           type="l",
           col="grey",
           ylab="Residual",
           lwd=2,
           log=log)

      ##ad 0 line
      abline(h=0)


      ##------------------------------------------------------------------------#
      ##++component to sum contribution plot ++##
      ##------------------------------------------------------------------------#

      ##plot component contribution to the whole signal
      #open plot area
      par(mar=c(4,4,3.2,0))
      plot(NA,NA,
           xlim=xlim,
           ylim=c(0,100),
           ylab="Contribution [%]",
           xlab=xlab,
           main="Component contribution to sum curve",
           log=if(log=="xy"){"x"}else{log})

      stepping <- seq(3,length(component.contribution.matrix),2)

      for(i in 1:length(xm)){

        polygon(c(component.contribution.matrix[,1],
                  component.contribution.matrix[,2]),
                c(component.contribution.matrix[,stepping[i]],
                  component.contribution.matrix[,stepping[i]+1]),
                col = col[i+1])
      }
      rm(stepping)

      ##reset par
      par(par.default)
      rm(par.default)
      ##------------------------------------------------------------------------##
    }#end if try-error for fit

    if(fun==TRUE){sTeve()}
  }
  ##-----------------------------------------------------------------------------
  ##remove objects
  try(unlist("parameters"))

  ##============================================================================#
  ## Return Values
  ##============================================================================#

  newRLumResults.fit_LMCurve <- set_RLum(
    class = "RLum.Results",
    data = list(
      fit = fit,
      output.table = output.table,
      component.contribution.matrix = list(component.contribution.matrix),
      call = sys.call()
    )
  )

  invisible(newRLumResults.fit_LMCurve)

}
