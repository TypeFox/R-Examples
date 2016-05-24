#' Model Luminescence Signals
#'
#' This function models luminescence signals for quartz based on published physical models.
#' It is possible to simulate TL, (CW-) OSL, RF measurements in a arbitrary sequence. This
#' sequence is definded as a \code{\link{list}} of certain abrivations. Furthermore it is possible to
#' load a sequence direct from the Riso Sequence Editor.
#' The output is an \code{\linkS4class{RLum.Analysis}}object and so the plots are done
#' by the \code{\link{plot_RLum.Analysis}} function. If a SAR sequence is simulated the plot output can be disabled and SAR analyse functions
#' can be used.
#'
#'
#' Defining a \bold{sequence}\cr
#'
#' \tabular{lll}{
#' \bold{Arguments} \tab \bold{Description} \tab \bold{Sub-arguments}\cr
#' TL \tab thermally stimulated luminescence \tab 'temp begin', 'temp end', 'heating rate'\cr
#' OSL\tab optically stimulated luminescence \tab 'temp', 'duration','optical_power'\cr
#' ILL\tab illumination \tab 'temp', 'duration','optical_power'\cr
#' LM_OSL\tab linear modulated OSL \tab 'temp', 'duration', optional: 'start_power', 'end_power'\cr
#' RL/RF\tab radioluminescence\tab 'temp','dose', 'dose_rate' \cr
#' IRR\tab irradiation \tab 'temp','dose', 'dose_rate' \cr
#' CH \tab cutheat \tab 'temp', optional: 'duration', 'heating_rate' \cr
#' PH  \tab preheat \tab 'temp', 'duration' optional: 'heating_rate' \cr
#' PAUSE \tab pause \tab 'temp', 'duration'
#' }
#'
#' Note: 100 \% illumination power equates 20 mW/cm^2
#'
#'
#' Defining a \bold{SAR-sequence}\cr
#'
#' \tabular{lll}{
#' \bold{Abrivation} \tab \bold{Description} \tab \bold{examples} \cr
#' RegDose \tab Dose points of the regenerative cycles\tab c(0, 80, 140, 260, 320, 0, 80)\cr
#' TestDose\tab Test dose for the SAR cycles  \tab 50 \cr
#' PH\tab Temperature of the preheat \tab 240 \cr
#' CH\tab Temperature of the cutheat \tab 200 \cr
#' OSL_temp\tab Temperature of OSL read out\tab  125 \cr
#' OSL_duration\tab  Duration of OSL read out\tab default: 40 \cr
#' Irr_temp \tab Temperature of irradiation \tab default: 20\cr
#' PH_duration  \tab Duration of the preheat \tab default: 10 \cr
#' dose_rate \tab Dose rate of the laboratory irradiation source \tab default: 1 \cr
#' optical_power \tab Percentage of the full illumination power \tab default: 90 \cr
#' Irr_2recover \tab Dose to be recovered in a dose-recovery-test \tab 20
#' }
#'
#' @param sequence \code{\link{list}} (\bold{required}): set sequence to model as \code{\link{list}} or as *.seq file from the
#' Riso sequence editor. To simulate SAR measurements there is an extra option to set the sequence list (cf. details).
#
#' @param model \code{\link{character}} (\bold{required}): set model to be used. Available models are:
#' "Bailey2001", "Bailey2002", "Bailey2004", "Pagonis2007", "Pagonis2008"
#'
#' @param lab.dose_rate \code{\link{numeric}} (with default): laboratory dose rate in XXX
#' Gy/s for calculating seconds into Gray in the *.seq file.
#'
#' @param simulate_sample_history \code{\link{logical}} (with default): FALSE (with default): simulation begins at laboratory conditions, TRUE: simulations begins at crystallization (all levels 0)
#' process
#'
#' @param plot \code{\link{logical}} (with default): Enables or disables plot output
#'
#' @param verbose \code{\link{logical}} (with default): Verbose mode on/off
#'
#' @param show.structure \code{\link{logical}} (with default): Shows the structure of the result.
#' Recommended to show record.id to analyse concentrations.
#'
#' @param \dots further arguments and graphical parameters passed to
#' \code{\link{plot.default}}. See details for further information.
#'
#' @return This function returns an \code{\linkS4class{RLum.Analysis}} object with all TL, (LM-) OSL and RF/RL steps
#' in the sequence. Every entry is an \code{\linkS4class{RLum.Data.Curve}} object and can be plotted, analysed etc. with
#' further \code{RLum}-functions.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany),
#' Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @references
#'
#' Bailey, R.M., 2001. Towards a general kinetic model for optically and thermally stimulated
#' luminescence of quartz. Radiation Measurements 33, 17-45.
#'
#' Bailey, R.M., 2002. Simulations of variability in the luminescence characteristics of natural
#' quartz and its implications for estimates of absorbed dose.
#' Radiation Protection Dosimetry 100, 33-38.
#'
#' Bailey, R.M., 2004. Paper I-simulation of dose absorption in quartz over geological timescales
#' and it simplications for the precision and accuracy of optical dating.
#' Radiation Measurements 38, 299-310.
#'
#' Pagonis, V., Chen, R., Wintle, A.G., 2007: Modelling thermal transfer in optically
#' stimulated luminescence of quartz. Journal of Physics D: Applied Physics 40, 998-1006.
#'
#' Pagonis, V., Wintle, A.G., Chen, R., Wang, X.L., 2008. A theoretical model for a new dating protocol
#' for quartz based on thermally transferred OSL (TT-OSL).
#' Radiation Measurements 43, 704-708.
#'
#' Soetaert, K., Cash, J., Mazzia, F., 2012. Solving differential equations in R.
#' Springer Science & Business Media.
#'
#' @seealso \code{\link{plot}}, \code{\linkS4class{RLum}},
#' \code{\link{read_SEQ2R}}
#'
#' @examples
#'
#'
#' ##================================================================##
#' ## Example 1: Simulate sample history of Bailey2001
#' ## (cf. Bailey, 2001, Fig. 1)
#' ##================================================================##
#'
#' ##set sequence with the following steps
#' ## (1) Irradiation at 20 deg. C with a dose of 10 Gy and a dose rate of 1 Gy/s
#' ## (2) TL from 20-400 deg. C with a rate of 5 K/s
#'
#' sequence <-
#'   list(
#'     IRR = c(20, 10, 1),
#'     TL = c(20, 400, 5)
#'   )
#'
#' ##model sequence
#' model.output <- model_LuminescenceSignals(
#'   sequence = sequence,
#'   model = "Bailey2001"
#' )
#'
#' \dontrun{
#' ##============================================================================##
#' ## Example 2: Simulate sequence at labour without sample history
#' ##============================================================================##
#'
#' ##set sequence with the following steps
#' ## (1) Irraditation at 20 deg. C with a dose of 100 Gy and a dose rate of 1 Gy/s
#' ## (2) Preheat to 200 deg. C and hold for 10 s
#' ## (3) LM-OSL at 125 deg. C. for 100 s
#' ## (4) Cutheat at 200 dec. C.
#' ## (5) Irraditation at 20 deg. C with a dose of 10 Gy and a dose rate of 1 Gy/s
#' ## (6) Pause at 200 de. C. for 100 s
#' ## (7) OSL at 125 deg. C for 100 s with 90 % optical power
#' ## (8) Pause at 200 deg. C for 100 s
#' ## (9) TL from 20-400 deg. C with a heat rate of 5 K/s
#' ## (10) Radiofluorescence at 20 deg. C with a dose of 200 Gy and a dose rate of 0.01 Gy/s
#'
#' sequence <-
#'  list(
#'    IRR = c(20, 100, 1),
#'    PH = c(200, 10),
#'    LM_OSL = c(125, 100),
#'    CH = c(200),
#'    IRR = c(20, 10, 1),
#'    PAUSE = c(200, 100),
#'    OSL = c(125, 100, 90),
#'    PAUSE = c(200, 100),
#'    TL = c(20, 400, 5),
#'    RF = c(20, 200, 0.01)
#' )
#'
#' # call function "model_LuminescenceSignals", set sequence = sequence,
#' # model = "Pagonis2008" (palaeodose = 200 Gy) and simulate_sample_history = FALSE (default),
#' # because the sample history is not part of the sequence
#'
#' model.output <- model_LuminescenceSignals(
#'    sequence = sequence,
#'    model = "Pagonis2008"
#'    )
#'
#'
#'
#' ##============================================================================##
#' ## Example 3: Simulate SAR sequence
#' ##============================================================================##
#'
#' ##set SAR sequence with the following steps
#' ## (1) RegDose: set regenerative dose [Gy] as vector
#' ## (2) TestDose: set test dose [Gy]
#' ## (3) PH: set preheat temperature in deg. C
#' ## (4) CH: Set cutheat temperature in deg. C
#' ## (5) OSL_temp: set OSL reading temperature in deg. C
#' ## (6) OSL_duration: set OSL reading duration in s
#'
#' sequence <- list(
#'  RegDose = c(0,10,20,50,90,0,10),
#'  TestDose = 5,
#'  PH = 240,
#'  CH = 200,
#'  OSL_temp = 125,
#'  OSL_duration = 70)
#'
#' # call function "model_LuminescenceSignals", set sequence = sequence,
#' # model = "Pagonis2007" (palaeodose = 20 Gy) and simulate_sample_history = FALSE (default),
#' # because the sample history is not part of the sequence
#'
#'  model.output <- model_LuminescenceSignals(
#'
#'  sequence = sequence,
#'  model = "Pagonis2007",
#'  plot = FALSE
#'  )
#'
#' # in environment is a new object "model.output" with the results of
#' # every step of the given sequence.
#' # Plots are done at OSL and TL steps and the growth curve
#'
#' # call "analyse_SAR.CWOSL" from RLum package
#'  results <- analyse_SAR.CWOSL(model.output,
#'                             signal.integral.min = 1,
#'                             signal.integral.max = 15,
#'                             background.integral.min = 601,
#'                             background.integral.max = 701,
#'                             fit.method = "EXP",
#'                             dose.points = c(0,10,20,50,90,0,10))
#'
#'
#' ##============================================================================##
#' ## Example 4: generate sequence from *.seq file and run SAR simulation
#' ##============================================================================##
#'
#' # load example *.SEQ file and construct a sequence.
#' # call function "model_LuminescenceSignals", load created sequence for sequence,
#' # set model = "Bailey2002" (palaeodose = 10 Gy)
#' # and simulate_sample_history = FALSE (default),
#' # because the sample history is not part of the sequence
#'
#' path <- system.file("extdata", "example_SAR_cycle.SEQ", package="RLumModel")
#'
#' sequence <- read_SEQ2R(file = path)
#'
#' model.output <- model_LuminescenceSignals(
#'   sequence = sequence,
#'   model = "Bailey2001",
#'   plot = FALSE
#' )
#'
#'
#' ## call RLum package function "analyse_SAR.CWOSL" to analyse the simulated SAR cycle
#'
#' results <- analyse_SAR.CWOSL(model.output,
#'                              signal.integral.min = 1,
#'                              signal.integral.max = 10,
#'                              background.integral.min = 301,
#'                              background.integral.max = 401,
#'                              dose.points = c(0,8,14,26,32,0,8),
#'                              fit.method = "EXP")
#'
#' print(get_RLum(results))
#'
#'
#' ##============================================================================##
#' ## Example 5: compare different optical powers of stimulation light
#' ##============================================================================##
#'
#' # call function "model_LuminescenceSignals", model = "Bailey2004"
#' # and simulate_sample_history = FALSE (default),
#' # because the sample history is not part of the sequence
#' # the optical_power of the LED is varied and then compared.
#'
#' optical_power <- seq(from = 0,to = 100,by = 20)
#'
#' model.output <- lapply(1:length(optical_power), function(x){
#'
#'  sequence <- list(IRR = c(20, 50, 1),
#'                   PH = c(220, 10, 5),
#'                   OSL = c(125, 50, optical_power[x])
#'                   )
#'
#'  data <- model_LuminescenceSignals(
#'            sequence = sequence,
#'            model = "Bailey2004",
#'            plot = FALSE
#'            )
#'
#'  return(get_RLum(data, recordType = "OSL$", drop = FALSE))
#' })
#'
#' ##combine output curves
#' model.output.merged <- merge_RLum(model.output)
#'
#' ##plot
#' plot_RLum(
#'  object = model.output.merged,
#'  xlab = "Illumination time [s]",
#'  ylab = "OSL signal [a.u.]",
#'  main = "OSL signal dependency on optical power of stimulation light",
#'  legend.text = paste("Optical power density", 20*optical_power/100, "mW/cm^2"),
#'  combine = TRUE)
#'
#'}
#' @export
model_LuminescenceSignals <- function(
  model,
  sequence,
  lab.dose_rate = 1,
  simulate_sample_history = FALSE,
  plot = TRUE,
  verbose = TRUE,
  show.structure = FALSE,
  ...
) {


# Integrity tests and conversion --------------------------------------------------------------

  #Check if model is supported
  model.allowed_keywords <- c("Bailey2001", "Bailey2004", "Pagonis2008", "Pagonis2007", "Bailey2002")

  if(!model%in%model.allowed_keywords){
    stop(paste0("[model_LuminescenceSignals()] Model not supported. Supported models are: ", paste(model.allowed_keywords, collapse = ", ")))

  }

  #Check sequence
  if(is(sequence,"character")){

    sequence <- read_SEQ2R(
      file = sequence,
      lab.dose_rate = lab.dose_rate,
      txtProgressBar = ifelse(verbose, TRUE, FALSE)
    )

  }

  else if(is.list(sequence)){

    if(!is.numeric(unlist(sequence))){
      stop("[model_LuminescenceSignals()] Sequence comprises non-numeric arguments!")
    }

    if("RegDose"%in%names(sequence)){# test if .create_SAR.sequence is requiered

      RegDose = sequence$RegDose
      TestDose = sequence$TestDose
      PH = sequence$PH
      CH = sequence$CH
      OSL_temp = sequence$OSL_temp

      Irr_temp = sequence$Irr_temp
      if(is.null(Irr_temp)){

        Irr_temp <- 20
      }

      OSL_duration = sequence$OSL_duration
      if(is.null(sequence$OSL_duration)){

        OSL_duration <- 40
      }

      PH_duration = sequence$PH_duration
      if(is.null(PH_duration)){

        PH_duration <- 10
      }

      dose_rate = sequence$dose_rate
      if(is.null(dose_rate)){

        dose_rate <- lab.dose_rate
      }

      optical_power = sequence$optical_power
      if(is.null(optical_power)){

        optical_power <- 90
      }


      if(!"Irr_2recover"%in%names(sequence)){# SAR sequence

      sequence <- .create_SAR.sequence(
        RegDose = RegDose,
        TestDose = TestDose,
        PH = PH,
        CH = CH,
        OSL_temp = OSL_temp,
        Irr_temp = Irr_temp,
        OSL_duration = OSL_duration,
        PH_duration = PH_duration,
        dose_rate = dose_rate,
        optical_power = optical_power
      )}
      else{# DRT sequence

        sequence <- .create_DRT.sequence(
          RegDose = RegDose,
          TestDose = TestDose,
          PH = PH,
          CH = CH,
          OSL_temp = OSL_temp,
          Irr_temp = Irr_temp,
          OSL_duration = OSL_duration,
          PH_duration = PH_duration,
          dose_rate = dose_rate,
          optical_power = optical_power,
          Irr_2recover = sequence$Irr_2recover
          )}


    }else{

      sequence <- sequence
    }
  }

  else{

    stop("[model_LuminescenceSignals()] Sequence has to be of class list or a *.seq file")
  }


  #check for wrong elements in the sequence

    ##allowed keywords
    sequence.allowed_keywords <- c("IRR","PH", "CH", "TL", "OSL", "PAUSE", "LM_OSL", "RL", "RF", "ILL")

    ##check
    if(!all(names(sequence)%in%sequence.allowed_keywords)){
      stop(paste0("[model_LuminescenceSignals()] Unknow sequence arguments: Allowed arguments are: ", paste(sequence.allowed_keywords, collapse = ", ")))

    }

  #check if lab.dose_rate > 0
    if(lab.dose_rate <= 0){

      stop("[model_LuminescenceSignals()] lab.dose_rate has to be a positive number! ")

    }

# Load model parameters ------------------------------------------------------------------------------------

    parms <- .set_pars(model)
    if(simulate_sample_history == TRUE){
      n <- Luminescence::set_RLum(class = "RLum.Results",
                                  data = list(n = rep(0,length(parms$N)+2),
                                              temp = 20,
                                              model = model))
    } else {
      n <- parms$n
    }


# sequence ------------------------------------------------------------------------------------

  #sequence, n and parms as arguments for the SequenceTranslator, who translates the sequence to different model steps
    model.output <-
      .translate_sequence(
        sequence = sequence,
        n = n,
        model = model,
        parms = parms,
        verbose = verbose
      )


# Plot settings -------------------------------------------------------------------------------

  if(plot){

    plot.data <- get_RLum(model.output, recordType = c("RF$", "TL$", "OSL$", "LM-OSL$"), drop = FALSE)
    Luminescence::plot_RLum(plot.data, ...)
  }

# model.output structure --------------------------------------------------

  if(show.structure){
    cat("[model_LuminescenceSignals()] \n\t>> Structure of output from 'model_LuminescenceSignals()' \n\n")
    print(Luminescence::structure_RLum(model.output))
  }

# return model.output -------------------------------------------------------------------------
  return(model.output)

}#end function
