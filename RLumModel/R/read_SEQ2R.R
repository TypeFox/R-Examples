#' Parse a Risoe SEQ-file to a sequence neccessary for simulating quartz luminescence
#'
#' A SEQ-file created by the Risoe Sequence Editor can be imported to simulate the sequence written
#' in the sequence editor.
#'
#' \bold{Supported versions}\cr
#'
#' Supppored and tested: version 4.36.\cr
#'
#' @param file \code{\link{character}} (\bold{required}): a *.seq file created by the Risoe Sequence Editor
#'
#' @param lab.dose_rate \code{\link{character}} (with default): set the dose rate of the radiation source
#' in the laboratory [Gy/s]. Default: 1 Gy/s
#'
#' @param txtProgressBar \code{\link{logical}} (with default): enables or disables the txtProgressBar for a visuell
#' control of the progress. Default: txtProgressBar = TRUE
#'
#' @return This function returns a \code{\link{list}} with the parsed *.seq file and the required steps for
#' \code{\link{model_LuminescenceSignals}}.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany),
#'
#' @references
#'
#' Riso: Sequence Editor User Manual.
#' Available at: \url{http://www.nutech.dtu.dk/english/-/media/Andre_Universitetsenheder/Nutech/Produkter
#' %20og%20services/Dosimetri/radiation_measurement_instruments/tl_osl_reader/Manuals/
#' SequenceEditor.ashx?la=da}
#'
#' @seealso \code{\link{model_LuminescenceSignals}}, \code{\link{readLines}}
#'
#' @examples
#' ##search "example_SAR_cycle.SEQ" in "extdata" in package "RLumModel"
#' path <- system.file("extdata", "example_SAR_cycle.SEQ", package="RLumModel")
#'
#' sequence <- read_SEQ2R(file = path)
#'
#' @export
read_SEQ2R <- function(
  file,
  lab.dose_rate = 1,
  txtProgressBar = TRUE
){

# Integrity tests and conversion --------------------------------------------------------------

if(class(file)!= "character"){
  stop("[read_SEQ2R()] class of file has to be a character.")
}

if(!file.exists(file)){
  stop("[read_SEQ2R()] file name doesn't seem to exist.")
}

if(lab.dose_rate < 0){
  stop("[read_SEQ2R()] Argument 'lab.dose_rate' has to be positiv.")
}

# parse *.SEQ file --------------------------------------------------------

file2read <- readLines(file, encoding = "UTF-8")


##(1)
##get all rows with the term "[Cell(...)]" - that's what we are interested in and it defines
##the number of elements we need
records.row_number <- grep(pattern = "\\[Cell\\(", x = file2read)

##(2)
##make a list with data of each sequence step
data.list <- lapply(1:length(records.row_number), function(x) {

  ##grep each element
  if (!is.na(records.row_number[x + 1])) {
    return(file2read[records.row_number[x]:(records.row_number[x + 1] - 1)])
  }
  else{
    return(file2read[records.row_number[x]:length(file2read)])

  }

})

##PROGRESS BAR
if(txtProgressBar){
  cat("\n [read_SEQ2R()] \n\t Parse *.seq file to sequence for RLumModel\n")
  pb <- txtProgressBar(min=0,max=length(data.list), char = "=", style=3)
}

sequence <- list()
names <- character()
for(x in 1:length(data.list)){

  ##get length of record
  sequence.ID <- as.numeric(gsub("Command=", "", data.list[[x]][2]))

  #identify ID with sequence step
  if(sequence.ID >= 1 && sequence.ID <=5){

    stop(paste("[read_SEQ2R()] Step ",x," of your sequence is not supported!", sep = ""))
  }

  if(sequence.ID == 6){#TL

    temp_begin <- as.numeric(gsub("Low=", "", data.list[[x]][grep(pattern = "^Low=",x = data.list[[x]])]))
    temp_end <- as.numeric(gsub("High=", "", data.list[[x]][grep(pattern = "^High=",x = data.list[[x]])]))
    b <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate",x = data.list[[x]])]))
    PH_time <- as.numeric(gsub("PhTime=", "", data.list[[x]][grep(pattern = "^PhTime",x = data.list[[x]])]))
    PH_temp <- as.numeric(gsub("PhTemp=", "", data.list[[x]][grep(pattern = "^PhTemp",x = data.list[[x]])]))

    names[x] <- "TL"
    sequence[[x]] <- c(temp_begin,temp_end,b,PH_time,PH_temp)


  }

  if(sequence.ID == 7){#OSL

    temp <- as.numeric(gsub("Temperature=", "", data.list[[x]][grep(pattern = "^Temperature=",x = data.list[[x]])]))
    duration <- as.numeric(gsub("High=", "", data.list[[x]][grep(pattern = "^High=",x = data.list[[x]])]))
    optical_power <- as.numeric(gsub("Current=", "", data.list[[x]][grep(pattern = "^Current=",x = data.list[[x]])]))
    b <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate=",x = data.list[[x]])]))
    lightsource <- as.numeric(gsub("Lightsource=", "", data.list[[x]][grep(pattern = "^Lightsource=",x = data.list[[x]])]))

    if(lightsource == 2){
      names[x] <- "IR_OSL"
    }
    if(lightsource == 4){
      names[x] <- "OSL"

    }
    sequence[[x]] <- c(temp,duration,optical_power,b)


  }
  if(sequence.ID == 8){#TOL

    stop(paste("[read_SEQ2R()] Step ",x," \"TOL\" of your sequence is not supported!", sep = ""))

  }
  if(sequence.ID == 9){#TR-POSL

    stop(paste("[read_SEQ2R()] Step ",x," \"TR-POSL\" of your sequence is not supportel!", sep = ""))
  }

  if(sequence.ID == 10){# Irradiation


    temp <- as.numeric(gsub("An_Temp=", "", data.list[[x]][grep(pattern = "^An_Temp=",x = data.list[[x]])]))
    duration <- as.numeric(gsub("Irr_Time=", "", data.list[[x]][grep(pattern = "^Irr_Time",x = data.list[[x]])]))
    if(duration == 0){

      duration <- 1e-13
    }

    dose <- duration*lab.dose_rate

    b <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate",x = data.list[[x]])]))
    PH_time <- as.numeric(gsub("An_Time=", "", data.list[[x]][grep(pattern = "^An_Time",x = data.list[[x]])]))

    names[x] <- "IRR"
    sequence[[x]] <- c(temp,dose,lab.dose_rate,b,PH_time)

  }

  if(sequence.ID == 11){# illumination

    temp <- as.numeric(gsub("Temperature=", "", data.list[[x]][grep(pattern = "^Temperature=",x = data.list[[x]])]))
    duration <- as.numeric(gsub("^Bl_Time=", "", data.list[[x]][grep(pattern = "^Bl_Time",x = data.list[[x]])]))
    optical_power <- as.numeric(gsub("Power=", "", data.list[[x]][grep(pattern = "^Power",x = data.list[[x]])]))

    b <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate",x = data.list[[x]])]))
    PH_time <- as.numeric(gsub("An_Time=", "", data.list[[x]][grep(pattern = "^An_Time",x = data.list[[x]])]))

    names[x] <- "ILL"
    sequence[[x]] <- c(temp,duration,optical_power,b,PH_time)


  }

  if(sequence.ID == 12){# PH

    PH_temp <- as.numeric(gsub("An_Temp=", "", data.list[[x]][grep(pattern = "^An_Temp",x = data.list[[x]])]))
    PH_time <- as.numeric(gsub("An_Time=", "", data.list[[x]][grep(pattern = "^An_Time",x = data.list[[x]])]))
    b <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate",x = data.list[[x]])]))

    names[x] <- "PH"
    sequence[[x]] <- c(PH_temp,PH_time,b)

  }

  if(sequence.ID == 13){# set temp

    PH_temp <- as.numeric(gsub("An_Temp=", "", data.list[[x]][grep(pattern = "^An_Temp",x = data.list[[x]])]))

    names[x] <- "PH"
    sequence[[x]] <- c(PH_temp)

  }

  if(sequence.ID == 14){# Pause

    Pause_time <- as.numeric(gsub("An_Time=", "", data.list[[x]][grep(pattern = "^An_Time",x = data.list[[x]])]))

    names[x] <- "PAUSE"
    sequence[[x]] <- c(Pause_time)

  }

  if(sequence.ID == 15){

    stop(paste("[read_SEQ2R()] Step ",x," of your sequence is not supported!", sep = ""))
  }

  if(sequence.ID == 16){

    stop(paste("[read_SEQ2R()] Step ",x," of your sequence is not supported!", sep = ""))
  }

  if(sequence.ID == 17){

    stop("[read_seq2R()] Low-Level: Not supported!")
  }

  if(sequence.ID == 18){# LM-OSL

    temp <- as.numeric(gsub("Temperature=", "", data.list[[x]][grep(pattern = "^Temperature=",x = data.list[[x]])]))
    duration <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate",x = data.list[[x]])]))
    start_power <- as.numeric(gsub("Low=", "", data.list[[x]][grep(pattern = "^Low",x = data.list[[x]])]))
    end_power <- as.numeric(gsub("High=", "", data.list[[x]][grep(pattern = "^High",x = data.list[[x]])]))

    b <- as.numeric(gsub("HeatRate=", "", data.list[[x]][grep(pattern = "^HeatRate",x = data.list[[x]])]))
    PH_time <- as.numeric(gsub("PreHeatTime=", "", data.list[[x]][grep(pattern = "^PreHeatTime",x = data.list[[x]])]))

    names[x] <- "LM_OSL"
    sequence[[x]] <- c(temp,duration,start_power,end_power,b,PH_time)


  }

  if(sequence.ID == 19){

    stop(paste("[read_SEQ2R()] Step ",x," \"SG OSL\" of your sequence is not supported!", sep = ""))

  }

  if(sequence.ID == 20){

    stop(paste("[read_SEQ2R()] Step ",x," \"User defined\" of your sequence is not supported!", sep = ""))

  }


  if(sequence.ID == 21){

    stop(paste("[read_SEQ2R()] Step ",x," \"SG LM-OSL\" of your sequence is not supported!", sep = ""))

  }

  if(sequence.ID == 22){

    stop(paste("[read_SEQ2R()] Step ",x," of your sequence is not supported!", sep = ""))

  }

  if(sequence.ID == 23){

    stop(paste("[read_SEQ2R()] Step ",x," \"POSL\" of your sequence is not supported!", sep = ""))

  }


  if(sequence.ID == 28){#RL

    temp <- as.numeric(gsub("Temperature=", "", data.list[[x]][grep(pattern = "^Temperature=",x = data.list[[x]])]))
    duration <- as.numeric(gsub("High=", "", data.list[[x]][grep(pattern = "^High",x = data.list[[x]])]))
    lab.dose_rate <- lab.dose_rate

    dose <- duration*lab.dose_rate
    b <- as.numeric(gsub("Rate=", "", data.list[[x]][grep(pattern = "^Rate",x = data.list[[x]])]))

    names[x] <- "RL"
    sequence[[x]] <- c(temp,dose,lab.dose_rate,b)
  }

  if(sequence.ID == 29){

    stop(paste("[read_SEQ2R()] Step ",x," \"XRF\" of your sequence is not supported!", sep = ""))

  }

  if(sequence.ID == 30){

    stop(paste("[read_SEQ2R()] Step ",x," \"Photo\" of your sequence is not supported!", sep = ""))

  }



  ##update progress bar
  if (txtProgressBar) {
    setTxtProgressBar(pb, x)
  }


  }
  ##close ProgressBar
  if(txtProgressBar){close(pb)}

  return(sequence <- setNames(object = sequence,nm = names))
}
