# This R-script performs diagnostic tests on PK/PD models using non-compartmental analysis parameters
# Chayan, 11/2014

#  roxygen comments
#' Performs NCA calculations and population PK model diagnosis.
#'
#' \pkg{ncappc} is a flexible tool, to
#' \enumerate{
#'  \item perform a traditional NCA
#'  \item perform simulation-based posterior predictive checks for a
#' population PK model using NCA metrics.
#' }
#'
#' Non-compartmental analysis (NCA) calculates pharmacokinetic (PK) metrics 
#' related to the systemic exposure to a drug following administration, e.g. 
#' area under the concentration-time curve and peak concentration. \pkg{ncappc} 
#' performs a traditional NCA using the observed plasma concentration-time data.
#' In the presence of simulated plasma concentration-time data, \pkg{ncappc} 
#' also performs simulation-based posterior predictive checks (ppc) using NCA 
#' metrics for the corresponding population PK (PopPK) model used to generate 
#' the simulated data. The diagnostic analysis is performed at the population as
#' well as the individual level. The distribution of the simulated population 
#' means of each NCA metric is compared with the corresponding observed 
#' population mean. The individual level comparison is performed based on the 
#' deviation of the mean of any NCA metric based on simulations for an 
#' individual from the corresponding NCA metric obtained from the observed data.
#' Additionaly, \pkg{ncappc} reports the normalized prediction distribution 
#' error (NPDE) of the simulated NCA metrics for each individual and their 
#' distribution within a population. \pkg{ncappc} produces two default outputs 
#' depending on the type of analysis performed, i.e., traditional NCA and PopPK 
#' diagnosis. The PopPK diagnosis feature of \pkg{ncappc} produces 7 sets of 
#' graphical outputs to assess the ability of a population model to simulate the
#' concentration-time profile of a drug and thereby identify model 
#' misspecification. In addition, tabular outputs are generated showing the 
#' values of the NCA metrics estimated from the observed and the simulated data,
#' along with the deviation, NPDE, regression parameters used to estimate the 
#' elimination rate constant and the related population statistics. The default 
#' values of the arguments used in \pkg{ncappc} are shown in the \strong{Useage}
#' section of this document and/or in \strong{bold} in the \strong{Arguments} 
#' section.
#'
#' @param obsFile Observed concentration-time data from an internal data frame
#'   or an external table with comma, tab or space as separator
#'   (\strong{"nca_original.npctab.dta"})
#' @param simFile NONMEM simulation output with the simulated concentration-time
#'   data from an internal data frame or an external table 
#'   (\strong{"nca_simulation.1.npctab.dta"}). \code{NULL} produces just the
#'   traditional NCA output, a filename or data frame prduces the NCA output as
#'   well as the PopPK diagnosis. If \code{new_data_method=TRUE} then this can
#'   be a compressed file as well.
#' @param str1Nm Column name for 1st level population stratifier
#'   (\strong{\code{NULL}})
#' @param str1 Stratification ID of the members within 1st level stratification 
#'   (e.g c(1,2)) (\strong{\code{NULL}})
#' @param str2Nm Column name for 2nd level population stratifier
#'   (\strong{\code{NULL}})
#' @param str2 Stratification ID of the members within 2nd level stratification
#'   (e.g c(1,2)) (\strong{\code{NULL}})
#' @param str3Nm Column name for 3rd level population stratifier
#'   (\strong{\code{NULL}})
#' @param str3 Stratification ID of the members within 3rd level stratification 
#'   (e.g c(1,2)) (\strong{\code{NULL}})
#' @param concUnit Unit of concentration (e.g. "ng/mL") (\strong{"M.L^-3"})
#' @param timeUnit Unit of time (e.g. "h") (\strong{"T"})
#' @param doseUnit Unit of dose amount (e.g. "ng") (\strong{"M"})
#' @param doseNormUnit Normalization factor of dose amount if used (e.g. "kg") 
#'   (\strong{\code{NULL}})
#' @param obsLog Concentration in observed data in logarithmic form
#'   (\code{TRUE}, \code{FALSE}) (\strong{\code{FALSE}})
#' @param simLog Concentration in simulated data in logarithmic form
#'   (\code{TRUE}, \code{FALSE}) (\strong{\code{FALSE}})
#' @param psnOut observed data is an output from PsN or in NONMEM output format 
#'   (\code{TRUE}, \code{FALSE}) (\strong{\code{TRUE}})
#' @param idNmObs Column name for ID in observed data (\strong{"ID"})
#' @param timeNmObs Column name for time in observed data (\strong{"TIME"})
#' @param concNmObs Column name for concentration in observed data
#'   (\strong{"DV"})
#' @param idNmSim Column name for ID in simulated data (\strong{"ID"})
#' @param timeNmSim Column name for time in simulated data (\strong{"TIME"})
#' @param concNmSim Column name for concentration in simulated data
#'   (\strong{"DV"})
#' @param AUCTimeRange User-defined window of time used to estimate AUC
#'   (\strong{\code{NULL}})
#' @param backExtrp If back-extrapolation is needed for AUC (\code{TRUE} or
#'   \code{FALSE}) (\strong{\code{FALSE}})
#' @param LambdaTimeRange User-defined window of time to estimate elimination 
#'   rate-constant. This argument lets the user to choose a specific window of
#'   time to be used to estimate the elimination rate constant (Lambda) in the
#'   elimination phase. The accepted format for the input to this argument is a 
#'   numeric array of two elements; \code{c(14,24)} will estimate the Lambda
#'   using the data within the time units 14 to 24.  If \code{NULL} then all
#'   times are considered.
#' @param LambdaExclude User-defined excluded observation time points for
#'   estimation of elimination rate-constant (\strong{\code{NULL}})
#' @param doseAmtNm Column name to specify dose amount (\strong{\code{NULL}})
#' @param adminType Route of administration
#'   (iv-bolus,iv-infusion,extravascular) (\strong{"extravascular"})
#' @param doseType Steady-state (ss) or nonsteady-state (ns) dose
#'   (\strong{"ns"})
#' @param doseTime Dose time prior to the first observation for steady-state data (\strong{\code{NULL}})
#' @param Tau Dosing interval for steady-state data (\strong{\code{NULL}})
#' @param TI Infusion duration (\strong{\code{NULL}})
#' @param method linear, log or linear-log (\strong{"linear-log"})
#' @param blqNm Name of BLQ column if used (\strong{\code{NULL}})
#' @param blqExcl Excluded BLQ value or logical condition (e.g. 1 or ">=1" or 
#'   c(1,">3")) (\strong{"1"})
#' @param evid Use EVID (\code{TRUE}, \code{FALSE}) (\strong{\code{TRUE}})
#' @param evidIncl Included EVID (\strong{"0"})
#' @param mdv Use MDV (\code{TRUE}(includes data for MDV==0), \code{FALSE})
#'   (\strong{\code{FALSE}})
#' @param filterNm Column name for filter (\strong{\code{NULL}})
#' @param filterExcl Filter identifier or logical condition used for row
#'   exclusion (e.g. c(1, 2, "<20", ">=100", "!=100")) (\strong{\code{NULL}})
#' @param negConcExcl Exclude -ve conc (\strong{\code{FALSE}})
#' @param param NCA parameters (AUClast, AUClower_upper, AUCINF_obs, 
#'   AUCINF_pred, AUMClast, Cmax, Tmax, HL_Lambda_z) (c(\strong{"AUClast",
#'   "Cmax"}))
#' @param timeFormat time format (number, H:M, H:M:S) (\strong{"number"})
#' @param dateColNm colunm name for date if used (Date, DATE)
#'   (\strong{\code{NULL}})
#' @param dateFormat date format (D-M-Y, D/M/Y or any other combination of
#'   D,M,Y) (\strong{\code{NULL}})
#' @param spread Measure of the spread of simulated data (\code{"ppi"} (95\%
#'   parametric prediction interval) or \code{"npi"} (95\% nonparametric
#'   prediction interval))
#' @param tabCol Output columns to be printed in the report in addition to ID, 
#'   dose and population strata information (list of NCA metrics in a string 
#'   array) (\strong{c("AUClast", "Cmax", "Tmax", "AUCINF_obs", "Vz_obs",
#'   "Cl_obs", "HL_Lambda_z")})
#' @param figFormat format of the produced figures (bmp, jpeg, tiff, png)
#'   (\strong{"tiff"})
#' @param noPlot Perform only NCA calculations without any plot generation
#'   (\code{TRUE} or \code{FALSE}) 
#' @param printOut Write/print output on the disk. No plot will be saved if 
#'   noPlot is set to \code{TRUE} (\code{TRUE}, \code{FALSE})
#'   (\strong{\code{TRUE}})
#' @param studyName Name of the study to be added as a description in the report
#'   (\strong{\code{NULL}})
#' @param new_data_method \code{TRUE} or \code{FALSE}. For testing a faster
#'   method of reading data. (\strong{\code{TRUE}})
#' @param overwrite_SIMDATA Can be \code{TRUE}, to create new information in the
#'   SIMDATA directory, \code{FALSE}, to use the information in the SIMDATA 
#'   directory or \code{NULL} to have a dialog come up to ask the user what to 
#'   do. (\strong{\code{NULL}})
#' @param outFileNm Additional tag to the name of the output html and pdf output
#'   file hyphenated to the standard ncappc report file name standard ncappc 
#'   report file name (\strong{Name of the observed data file})
#'
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import scales
#' @import gtable
#' @import knitr
#' @import xtable
#' @import reshape2
#' @import testthat
#' 
#' @return NCA results and diagnostic test results
#' @export
#' @examples
#' ncappc(obsFile=system.file("extdata","pkdata.csv",package="ncappc"), 
#' psnOut=FALSE,noPlot=TRUE,printOut=FALSE)
#'

ncappc <- function(obsFile="nca_original.npctab.dta",
                   simFile="nca_simulation.1.npctab.dta",
                   str1Nm=NULL,str1=NULL,
                   str2Nm=NULL,str2=NULL,str3Nm=NULL,str3=NULL,
                   concUnit=NULL,timeUnit=NULL,doseUnit=NULL,
                   doseNormUnit=NULL,obsLog=FALSE,simLog=FALSE,
                   psnOut=TRUE,idNmObs="ID",timeNmObs="TIME",
                   concNmObs="DV",idNmSim="ID",timeNmSim="TIME",
                   concNmSim="DV",AUCTimeRange=NULL,backExtrp=FALSE,
                   LambdaTimeRange=NULL,LambdaExclude=NULL,doseAmtNm=NULL,
                   adminType="extravascular",doseType="ns",doseTime=NULL,Tau=NULL,
                   TI=NULL,method="linear-log",blqNm=NULL,blqExcl=1,evid=TRUE,evidIncl=0,
                   mdv=FALSE,filterNm=NULL,filterExcl=NULL,negConcExcl=FALSE,
                   param=c("AUClast","Cmax"),timeFormat="number",dateColNm=NULL,
                   dateFormat=NULL,spread="npi",
                   tabCol=c("AUClast","Cmax","Tmax","AUCINF_obs","Vz_obs","Cl_obs","HL_Lambda_z"),
                   figFormat="tiff",noPlot=FALSE,printOut=TRUE,studyName=NULL,new_data_method=TRUE,
                   overwrite_SIMDATA=NULL,outFileNm=NULL){
  
  "..density.." <- "meanObs" <- "sprlow" <- "sprhgh" <- "AUClast" <- "AUCINF_obs" <- "Cmax" <- "Tmax" <- "FCT" <- "ID" <- "STR1" <- "STR2" <- "STR3" <- "NPDE" <- "mcil" <- "mciu" <- "sdu" <- "sducil" <- "sduciu" <- "scale_linetype_manual" <- "scale_color_manual" <- "xlab" <- "ylab" <- "guides" <- "guide_legend" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "grid.arrange" <- "unit.c" <- "grid.grab" <- "ggsave" <- "facet_wrap" <- "ggplot" <- "labs" <- "geom_point" <- "geom_errorbarh" <- "knit2html" <- "knit2pdf" <- "knit" <- "file_test" <- "tail" <- "read.csv" <- "read.table" <- "dev.off" <- "write.table" <- "head" <- "write.csv" <- "coef" <- "dist" <- "lm" <- "median" <- "na.omit" <- "percent" <- "qchisq" <- "qnorm" <- "qt" <- "quantile" <- "scale_y_continuous" <- "sd" <- "STRAT1" <- "STRAT2" <- "STRAT3" <- "sdcil" <- "sdciu" <- "str" <- NULL
  rm(list=c("..density..","meanObs","sprlow","sprhgh","AUClast","AUCINF_obs","Cmax","Tmax","FCT","ID","STR1","STR2","STR3","NPDE","mcil","mciu","sdu","sducil","sduciu","scale_linetype_manual","scale_color_manual","xlab","ylab","guides","guide_legend","theme","element_text","unit","element_rect","geom_histogram","aes","geom_vline","grid.arrange","unit.c","grid.grab","ggsave","facet_wrap","ggplot","labs","geom_point","geom_errorbarh","knit2html","knit2pdf","knit","file_test","tail","read.csv","read.table","dev.off","write.table","head","write.csv","coef","dist","lm","median","na.omit","percent","qchisq","qnorm","qt","quantile","scale_y_continuous","sd","STRAT1","STRAT2","STRAT3","sdcil","sdciu","str"))
  
  options(warning.length=5000)
  options(scipen=999)
  usrdir <- getwd()
  
  # Observed data
  if (is.null(obsFile)){stop("Name of the file with observed data is required\n")}
  if (!is.data.frame(obsFile)){
    if (!file_test("-f", obsFile)){stop("File for the observed data does not exist\n")}
    # read observed data file
    if (psnOut == FALSE){
      extn <- tail(unlist(strsplit(obsFile, ".", fixed=T)), n=1)
      if(extn=="csv"){indf <- read.csv(obsFile)}else{indf <- read.table(obsFile, header=T)}
    }else{
      print("Note: The observed data file is expected to be generated by PsN since psnOut is set to TRUE.")
      indf <- read.table(obsFile, header=T, skip=1)
    }
  }else{
    indf <- obsFile
  }
  
  # Simulated data
  # Simulated data is not supplied
  if(is.null(simFile)){print("Note: Simulated data file is not supplied. Only NCA module will be executed.")}
  
  # Simulated data is not found in the given path
  if ((!is.null(simFile)) && (!is.data.frame(simFile)) && (!file_test("-f", simFile))){print(paste0("Note: Simulated data file, ",simFile,", is not found in the working directory. Only NCA module will be executed.")); simFile <- NULL}
  
  # Check steady state dosing interval
  if (doseType == "ss"){
    if(is.null(Tau)){setwd(usrdir);stop("Time for dosing interval is required for steady-state data\n")}
    if(is.null(doseTime)){print("Note: Dose time prior to the first observation for steady-state data is missing. Steady state observation period will be estimated based on the first non-zero obverved concentration and Tau.\n")}
  }
  
  # Check for column names
  if (idNmObs%in%colnames(indf)==F | timeNmObs%in%colnames(indf)==F | concNmObs%in%colnames(indf)==F){
    setwd(usrdir);stop("Incorrect column names of ID, TIME and/or DV\n")
  }else{
    idCol   <- which(colnames(indf) == idNmObs)[1]; id <- unique(indf[,idCol])
    timeCol <- which(colnames(indf) == timeNmObs)[1]
    concCol <- which(colnames(indf) == concNmObs)[1]
  }
  
  # exclude data based on specific values on filter column (optional)
  if (!is.null(filterNm)){
    if(filterNm%in%colnames(indf)==T & !is.null(filterExcl)){
      # filterExcl  == values to be excluded
      filterCol <- which(colnames(indf) == filterNm)[1]
      for (i in 1:length(filterExcl)){
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i])){
          indf <- indf[indf[,filterCol] != filterExcl[i],]
        }else if(!grepl("[<>!=]", filterExcl[i])){
          indf <- indf[indf[,filterCol] != filterExcl[i],]
        }else{
          indf <- eval(parse(text=paste0("subset(indf, !",filterNm,"%in% indf[indf[,",filterCol,"]",filterExcl[i],",filterCol])")))
        }
      }
    }else{
      print("Note: Incorrect filterNm or filterExcl specification. filterNm will not be used to process the observed data.\n")
    }
  }
  
  # copy the input data before processing
  refdf <- indf
  
  # 1st level population stratification
  if (!is.null(str1Nm)){
    if (str1Nm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the 1st level stratification column\n")}
    if (is.null(str1)){
      str1 <- unique(sort(indf[,str1Nm]))
      #if(is.factor(indf[,str1Nm]))  str1 <- sort(unique(as.character(indf[,str2Nm])))
      #if(is.numeric(indf[,str1Nm]) && anyNA(indf[,str1Nm])) str1 <- c(unique(sort(indf[,str1Nm])),NA)
      #if(is.numeric(indf[,str1Nm]) && !anyNA(indf[,str1Nm])) str1 <- unique(sort(indf[,str1Nm]))
    }
    #for (i in 1:length(str1)){
    #  if (nrow(indf[indf[,str1Nm]==str1[i],]) == 0){setwd(usrdir);stop("1st level stratification ID does not match the values within 1st level stratification column\n")}
    #}
  }
  
  # 2nd level population stratification
  if (!is.null(str2Nm)){
    if (str2Nm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the 2nd level stratification column\n")}
    if (is.null(str2)){
      str2 <- unique(sort(indf[,str2Nm]))
      #if(is.factor(indf[,str2Nm]))  str2 <- sort(unique(as.character(indf[,str2Nm])))
      #if(is.numeric(indf[,str2Nm]) && anyNA(indf[,str2Nm])) str2 <- c(unique(sort(indf[,str2Nm])),NA)
      #if(is.numeric(indf[,str2Nm]) && !anyNA(indf[,str2Nm])) str2 <- unique(sort(indf[,str2Nm]))
    }
    #for (i in 1:length(str2)){
    #  if (nrow(indf[indf[,str2Nm]==str2[i],]) == 0){setwd(usrdir);stop("1st level stratification ID does not match the values within 1st level stratification column\n")}
    #}
  }
  
  # 3rd level population stratification
  if (!is.null(str3Nm)){
    if (str3Nm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the 3rd level stratification column\n")}
    if (is.null(str3)){
      str3 <- unique(sort(indf[,str3Nm]))
      #if(is.factor(indf[,str3Nm]))  str3 <- sort(unique(as.character(indf[,str2Nm])))
      #if(is.numeric(indf[,str3Nm]) && anyNA(indf[,str3Nm])) str3 <- c(unique(sort(indf[,str3Nm])),NA)
      #if(is.numeric(indf[,str3Nm]) && !anyNA(indf[,str3Nm])) str3 <- unique(sort(indf[,str3Nm]))
    }
    #for (i in 1:length(str3)){
    #  if (nrow(indf[indf[,str3Nm]==str3[i],]) == 0){setwd(usrdir);stop("1st level stratification ID does not match the values within 1st level stratification column\n")}
    #}
  }
  
  # check time range, if any
  if ((!is.null(AUCTimeRange)) && (length(AUCTimeRange) != 2 | class(AUCTimeRange) != "numeric")){
    print("Note: Incorrect time range for AUC calculation. AUCTimeRange will not be used.\n")
  }
  
  if ((!is.null(LambdaTimeRange)) && (length(LambdaTimeRange) != 2 | class(LambdaTimeRange) != "numeric")){
    print("Note: Incorrect time range for Lambda calculation. LambdaTimeRange will not be used.\n")
  }
  
  # check requirements for infusion data
  if (adminType == "iv-infusion" & is.null(TI) & ("AMT"%in%colnames(indf)==F | "RATE"%in%colnames(indf)==F)){setwd(usrdir);stop("Duration of the infusion time is needed if AMT and RATE are absent in the input data\n")}
  
  # Set backExtrp to FALSE in the presence of simulated data
  if (!is.null(simFile)){backExtrp <- FALSE}
  
  
  # Dose amount is extracted from doseAmtNm column
  if (!is.null(doseAmtNm)){
    if (doseAmtNm%in%colnames(indf)==T){
      doseAmtNm <- doseAmtNm
    }else if ("AMT"%in%colnames(indf)){
      doseAmtNm <- "AMT"
    }else{
      doseAmtNm <- NULL
      print("Note: Dose amount column name provided in doseAmtNm or AMT column does not exist in the observed data file. Dose related NCA metrics will not be estimated for the observed data.\n")
    }
  }else{
    if ("AMT"%in%colnames(indf)) doseAmtNm <- "AMT"
  }
  
  # Dose unit
  if (is.null(doseUnit)) doseUnit <- "M"
  
  # Preliminary description
  # Units for dose, time and conc
  dunit    <- ifelse(is.null(doseNormUnit), doseUnit, paste0(doseUnit,"/",doseNormUnit))
  tunit    <- ifelse(is.null(timeUnit), "T", timeUnit)
  cunit    <- ifelse(is.null(concUnit), "M.L^-3", concUnit)
  aucunit  <- paste0(tunit,"*",cunit)
  aumcunit <- paste0("(",tunit,"^2)*",cunit)
  clunit   <- paste0(dunit,"/(",aucunit,")")
  vlunit   <- paste0(dunit,"/",cunit)
  

  # ignore data with BLQ = 1 or user specified value (optional)
  if (!is.null(blqNm)){
    if(blqNm%in%colnames(indf) == T){
      blqCol <- which(colnames(indf) == blqNm)[1]
      for (i in 1:length(blqExcl)){
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i]) == T) {indf <- indf[indf[,blqNm] != blqExcl[i],]}
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i]) == F) {indf <- eval(parse(text=paste0("subset(indf, !",blqNm,"%in% indf[as.numeric(as.character(indf[,",blqCol,"]))",blqExcl[i],",blqCol])")))}
      }
    }else{
      print("Note: Incorrect BLQ column name. BLQ will not be used to process the observed data.\n")
    }
  }
  
  # include data based on specific values on EVID column (optional) but keep rows with TIME == 0
  if (evid == TRUE){
    if("EVID"%in%colnames(indf) == T){
      # uevid == unique values in EVID column
      # evidIncl == EVID values to be included
      # ievid == EVID values to be ignored
      uevid <- unique(as.numeric(as.character(indf$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
      if (length(ievid) != 0){for (i in 1:length(ievid)){indf <- indf[indf$EVID != ievid[i],]}}
    }else{
      print("Note: EVID column is not present. EVID will not be used to process the observed data.\n")
    }
  }
  
  # if MDV fiter is present, exclude data for MDV == 1 but keep rows with TIME == 0
  if (mdv == TRUE){
    if("MDV"%in%colnames(indf) == T){
      indf <- indf[indf$MDV == 0,]
    }else{
      print("Note: MDV column is not present. MDV will not be used to process the observed data.\n")
    }
  }
  
  # Appropriate date format
  if (!is.null(dateColNm)){
    if (dateColNm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the date column\n")}
    if (dateFormat == "D-M-Y"){
      dateFormat <- "%d-%m-%Y"
    }else if (dateFormat == "D-Y-M"){
      dateFormat <- "%d-%Y-%m"
    }else if (dateFormat == "M-D-Y"){
      dateFormat <- "%m-%d-%Y"
    }else if (dateFormat == "M-Y-D"){
      dateFormat <- "%m-%Y-%d"
    }else if (dateFormat == "Y-M-D"){
      dateFormat <- "%Y-%m-%d"
    }else if (dateFormat == "Y-D-M"){
      dateFormat <- "%Y-%d-%m"
    }else if (dateFormat == "D/M/Y"){
      dateFormat <- "%d/%m/%Y"
    }else if (dateFormat == "D/Y/M"){
      dateFormat <- "%d/%Y/%m"
    }else if (dateFormat == "M/D/Y"){
      dateFormat <- "%m/%d/%Y"
    }else if (dateFormat == "M/Y/D"){
      dateFormat <- "%m/%Y/%d"
    }else if (dateFormat == "Y/M/D"){
      dateFormat <- "%Y/%m/%d"
    }else if (dateFormat == "Y/D/M"){
      dateFormat <- "%Y/%d/%m"
    }else{setwd(usrdir);stop("Incorrect date format. Currently allowed date formats are \"D-M-Y\", \"D/M/Y\", or any combination of D, M and Y separated by either - or /\n")}
  }
  
  # Appropriate time format
  if (timeFormat != "number"){
    if(timeFormat == "H:M"){
      timeFormat <- "%H:%M"
    }else if(timeFormat == "H:M:S"){
      timeFormat <- "%H:%M:%S"
    }else{setwd(usrdir);stop("Incorrect time format. Currently allowed date formats are \"number\", \"H:M\", \"H:M:S\"\n")}
  }
  
  # Create empty data frame for output
  outData <- data.frame()
  
  # cpopStrNm = Combined stratifyting column names
  # npopStr   = Number of stratification levels
  # popStrNm1 = 1st level stratifying column name
  # popStr1   = 1st level stratification ID names
  # npopStr1  = Number of 1st level stratification ID names
  
  if(is.null(str1Nm)  & is.null(str2Nm)  & is.null(str3Nm)) {case<-1; npopStr<-0} # No stratification
  if(!is.null(str1Nm) & is.null(str2Nm)  & is.null(str3Nm)) {case<-2; cpopStrNm<-str1Nm; npopStr<-1; popStrNm1<-str1Nm; popStr1<-str1; npopStr1<-length(str1)} # Str1
  if(is.null(str1Nm)  & !is.null(str2Nm) & is.null(str3Nm)) {case<-2; cpopStrNm<-str2Nm; npopStr<-1; popStrNm1<-str2Nm; popStr1<-str2; npopStr1<-length(str2)} # Str2
  if(is.null(str1Nm)  & is.null(str2Nm)  & !is.null(str3Nm)){case<-2; cpopStrNm<-str3Nm; npopStr<-1; popStrNm1<-str3Nm; popStr1<-str3; npopStr1<-length(str3)} # Str3
  
  # Str1 & Str2
  if(!is.null(str1Nm) & !is.null(str2Nm) & is.null(str3Nm)){
    case<-3; cpopStrNm<-paste(str1Nm,str2Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStr1<-str1; popStr2<-str2; npopStr1<-length(str1); npopStr2<-length(str2)
  }
  # Str1 & Str3
  if(!is.null(str1Nm) & is.null(str2Nm) & !is.null(str3Nm)){
    case<-3; cpopStrNm<-paste(str1Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str3Nm; popStr1<-str1; popStr2<-str3; npopStr1<-length(str1); npopStr2<-length(str3)
  }
  # Str2 & Str3
  if(is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
    case<-3; cpopStrNm<-paste(str2Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str2Nm; popStrNm2<-str3Nm; popStr1<-str2; popStr2<-str3; npopStr1<-length(str2); npopStr2<-length(str3)
  }
  # Str1 & Str2 & Str3
  if(!is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
    case<-4; cpopStrNm<-paste(str1Nm,str2Nm,str3Nm,sep=", ")
    npopStr<-3
    popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStrNm3<-str3Nm
    popStr1<-str1; popStr2<-str2; popStr3<-str3
    npopStr1<-length(str1); npopStr2<-length(str2); npopStr3<-length(str3)
  }
  
  
  # Allowed NCA parameters
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  npr    <- length(param)
  fctNm  <- data.frame()
  for (p in 1:npr){
    if (param[p]%in%alwprm == F){setwd(usrdir);stop("Incorrect NCA parameters. Please select NCA parameters from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\"\n")}
    if (param[p] == "AUClast" | param[p] == "AUClower_upper" | param[p] == "AUCINF_obs" | param[p] == "AUCINF_pred"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste0(param[p]," (",cunit,"*",tunit,")",sep="")))
    }else if (param[p] == "AUMClast"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste0(param[p]," (",cunit,"*",tunit,"^2)")))
    }else if (param[p] == "Cmax"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste0(param[p]," (",cunit,")")))
    }else if (param[p] == "Tmax"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste0(param[p]," (",tunit,")")))
    }else if (param[p] == "HL_Lambda_z"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste0(param[p]," (",tunit,")")))
    }
  }
  
  obsFileNm <- ifelse(is.data.frame(obsFile), deparse(substitute(obsFile)), obsFile)
  
  if(!is.null(studyName)){
    txt <- paste0("Name of the study: \"",studyName,"\"")
    txt <- paste(txt,paste0("Name of the file with the observed data: \"",obsFileNm,"\""),sep="\n")
  }else{
    txt <- paste0("Name of the file with the observed data: \"",obsFileNm,"\"")
  }
  txt <- paste(txt,paste0("Route of administration: ",adminType),sep="\n")
  if (doseType == "ss"){
    txt <- paste(txt,paste0("Dose type: steady-state with dosing interval (Tau) of",Tau),sep="\n")
  }else{
    txt <- paste(txt,"Dose type: non-steady-state",sep="\n")
  }
  txt <- paste(txt,paste0("No. of population stratification level: ",npopStr),sep="\n")
  if(case==2){
    txt <- paste(txt,paste0("Population stratification column: ",cpopStrNm),sep="\n")
    txt <- paste(txt,paste0("Population stratification ID within ",popStrNm1,": ",paste(popStr1,collapse=", ")),sep="\n")
  }else if(case==3){
    txt <- paste(txt,paste0("Population stratification columns: ",cpopStrNm),sep="\n")
    txt <- paste(txt,paste0("1st level population stratification ID within ",popStrNm1,": ",paste(popStr1,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("2nd level population stratification ID within ",popStrNm2,": ",paste(popStr2,collapse=", ")),sep="\n")
  }else if(case==4){
    txt <- paste(txt,paste0("Population stratification columns: ",cpopStrNm),sep="\n")
    txt <- paste(txt,paste0("1st level population stratification ID within ",popStrNm1,": ",paste(popStr1,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("2nd level population stratification ID within ",popStrNm2,": ",paste(popStr2,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("3rd level population stratification ID within ",popStrNm3,": ",paste(popStr3,collapse=", ")),sep="\n")
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Set plot output dimensions
  if (npr<=2){ncol<-2;hth<-12;wth<-16;phth<-6;pwth<-7}else if(npr>2 & npr<=4){ncol<-2;hth<-18;wth<-18;phth<-9;pwth<-7}else if(npr>4 & npr<=6){ncol<-3;hth<-18;wth<-25;phth<-9;pwth<-9}else if(npr>6){ncol<-3;hth<-26;wth<-25;phth<-10;pwth<-9}
  
  # Initiate plot lists for pdf output
  concplot <- list(); histobsplot=list(); popplot <- list(); devplot <- list(); outlierplot <- list(); forestplot <- list(); npdeplot <- list(); histnpdeplot <- list()
  
  # calculate the NCA parameters for the observed data
  dset = "obs"
  
  # Function to generate time and conc data to calculate NCA metrics
  ncaId <- function(ifdf,ID){
    if(adminType == "iv-infusion" & is.null(TI)){
      amt  <- ifdf[ifdf[,idCol]==ID & ifdf$AMT > 0,"AMT"][1]
      rate <- ifdf[ifdf[,idCol]==ID & ifdf$RATE > 0,"RATE"][1]
      if (is.na(amt) | is.na(rate) | rate==0){setwd(usrdir);stop(paste0("Incorrect AMT and/or RATE value for ",ID))}else{TI <- amt/rate}
    }else{
      TI <- "NaN"
    }
    if (timeFormat != "number"){
      time <- numeric(0)
      if (!is.null(dateColNm)){
        tm <- as.POSIXct(paste(ifdf[ifdf[,idCol]==ID,dateColNm],ifdf[ifdf[,idCol]==ID,timeCol]),format=paste(dateFormat,timeFormat,sep=" "))
      }else{
        tm <- ifdf[ifdf[,idCol]==ID,timeCol]
      }
      for (j in 1:length(tm)){
        time[j] <- ifelse ((is.null(dateColNm)), as.numeric(difftime(strptime(tm[j], format=timeFormat), strptime(tm[1], format=timeFormat), units='hours')), as.numeric(difftime(strptime(tm[j], format="%Y-%m-%d %H:%M:%S"), strptime(tm[1], format="%Y-%m-%d %H:%M:%S"), units='hours')))
      }
    }else{
      time <- suppressWarnings(as.numeric(as.character(ifdf[ifdf[,idCol]==ID,timeCol])))
    }
    tconc <- ifdf[ifdf[,idCol]==ID, concCol]; conc <- numeric(0)
    for (c in 1:length(tconc)){conc <- c(conc, ifelse ((tconc[c]=="."), 0, as.numeric(as.character(tconc[c]))))}; rm(tconc)
    if (length(which(is.na(time)))!=0){
      zidx <- which(is.na(time))
      time <- time[-zidx]
      conc <- conc[-zidx]
    }
    if (length(which(is.na(conc)))!=0){
      zidx <- which(is.na(conc))
      time <- time[-zidx]
      conc <- conc[-zidx]
    }
    if (obsLog == TRUE){
      for (c in 1:length(conc)){conc[c] <- ifelse ((conc[c] != 0), exp(conc[c]), 0)}
    }
    
    tc <- data.frame(time,conc)
    if(nrow(tc)>0){
      tc <- tc[order(tc$time),]
      #if (tc$time[1]<0) tc$time <- tc$time + abs(min(tc$time))
    }
    return(tc)
  }
  
  # Estimate NCA metrics for case = 1
  pddf <- data.frame()
  
  if (case == 1){
    ifdf <- indf
    if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
    if (nrow(ifdf) == 0){next}
    idd        <- unique(as.character(ifdf[,idCol]))
    if(is.null(doseAmtNm)){
      doseAmount <- NA
    }else{
      doseAmount <- paste(as.numeric(unique(refdf[refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
    }
    
    # Description
    pddf  <- rbind(pddf, data.frame(a=doseAmount, b=length(idd)))
    cdata <- data.frame()
    for (i in 1:length(idd)){
      if (!is.null(doseAmtNm)){
        idzAmt <- as.numeric(refdf[refdf[,idNmObs]==as.character(idd[i]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
      }else{
        idzAmt <- NA
      }
      tc   <- ncaId(ifdf,as.character(idd[i]))
      if(nrow(tc)==0) next
      time <- as.numeric(tc$time)
      conc <- as.numeric(tc$conc)
      cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i])))
      NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
      outData <- rbind(outData, data.frame(ID=as.character(idd[i]),Dose=idzAmt,t(NCAprm)))
    }
    
    if(noPlot==FALSE){
      # DV plot
      figlbl    <- "All-data"
      cdata$FCT <- figlbl
      gdr       <- dv.plot(pdata=cdata,cunit=cunit,tunit=tunit)
      suppressMessages(suppressWarnings(grid.arrange(gdr)))
      ggr <- grid.grab()
      concplot[[length(concplot)+1]] <- ggr
      if (printOut==TRUE){
        fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
        eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=8,width=16,units=\"cm\",res=200)")))
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        dev.off()
      }
      
      # Obs hist plot
      if(nrow(outData)>=5){
        plotData    <- subset(outData, select=c(AUClast,AUCINF_obs,Cmax,Tmax))
        pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
        for (p in 1:length(pltPrm)){if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
        if (length(pltPrm) == 0) next
        figlbl      <- NULL
        histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=pltPrm,cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histobsgrob$gdr
        mylegend    <- histobsgrob$legend
        lheight     <- histobsgrob$lheight
        if (printOut==TRUE){
          fl <- paste0(usrdir,"/HistObs")
          eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)")))
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
        ggr <- grid.grab()
        histobsplot[[length(histobsplot)+1]] <- ggr
      }
    }
    cnm         <- c(paste0("Dose (",dunit,")"),"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 2
  if (case == 2){
    for (s1 in 1:npopStr1){
      cdata <- data.frame()
      ifdf <- indf[indf[,popStrNm1]==as.character(popStr1[s1]),]
      if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
      if (nrow(ifdf) == 0){next}
      idd <- unique(as.character(ifdf[,idCol]))
      if (!is.null(doseAmtNm)){
        doseAmount <- paste(as.numeric(unique(refdf[refdf[,popStrNm1]==as.character(popStr1[s1]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
      }else{
        doseAmount <- NA
      }
      # Description
      pddf <- rbind(pddf, data.frame(a=as.character(popStr1[s1]), b=doseAmount, c=length(idd)))
      for (i in 1:length(idd)){
        if (!is.null(doseAmtNm)){
          idzAmt <- as.numeric(refdf[refdf[,popStrNm1]==as.character(popStr1[s1]) & refdf[,idNmObs]==as.character(idd[i]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
        }else{
          idzAmt <- NA
        }
        tc   <- ncaId(ifdf,as.character(idd[i]))
        if(nrow(tc)==0) next
        time <- as.numeric(tc$time)
        conc <- as.numeric(tc$conc)
        
        cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i])))
        NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
        outData <- rbind(outData, data.frame(ID=as.character(idd[i]),STRAT1=popStr1[s1],Dose=idzAmt,t(NCAprm)))
      }
      
      if(noPlot==FALSE){
        # DV plot
        figlbl    <- paste0(popStrNm1,"-",as.character(popStr1[s1]))
        cdata$FCT <- figlbl
        gdr       <- dv.plot(pdata=cdata,cunit=cunit,tunit=tunit)
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        ggr <- grid.grab()
        concplot[[length(concplot)+1]] <- ggr
        if (printOut==TRUE){
          fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
          eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=8,width=16,units=\"cm\",res=200)")))
          suppressMessages(suppressWarnings(grid.arrange(gdr)))
          dev.off()
        }
        
        # Obs hist plot
        plotData    <- subset(outData, STRAT1==popStr1[s1], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
        if (nrow(plotData)<5) next
        pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
        for (p in 1:length(pltPrm)){if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
        if (length(pltPrm) == 0) next
        figlbl      <- paste0(popStrNm1,"-",as.character(popStr1[s1]))
        histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histobsgrob$gdr
        mylegend    <- histobsgrob$legend
        lheight     <- histobsgrob$lheight
        if (printOut==TRUE){
          fl <- paste0(usrdir,"/HistObs_",figlbl)
          eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)")))
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
        ggr <- grid.grab()
        histobsplot[[length(histobsplot)+1]] <- ggr
      }
    }
    cnm         <- c(popStrNm1,paste0("Dose (",dunit,")"),"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 3
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        cdata <- data.frame()
        ifdf <- indf[indf[,popStrNm1]==as.character(popStr1[s1]) & indf[,popStrNm2]==as.character(popStr2[s2]),]
        if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
        if (nrow(ifdf) == 0){next}
        idd <- unique(as.character(ifdf[,idCol]))
        if (!is.null(doseAmtNm)){
          doseAmount <- paste(as.numeric(unique(refdf[refdf[,popStrNm1]==as.character(popStr1[s1]) & refdf[,popStrNm2]==as.character(popStr2[s2]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
        }else{
          doseAmount <- NA
        }
        # Description
        pddf   <- rbind(pddf, data.frame(a=as.character(popStr1[s1]), b=as.character(popStr2[s2]), c=doseAmount, d=length(idd)))
        for (i in 1:length(idd)){
          if (!is.null(doseAmtNm)){
            idzAmt <- as.numeric(refdf[refdf[,popStrNm1]==as.character(popStr1[s1]) & refdf[,popStrNm2]==as.character(popStr2[s2]) & refdf[,idNmObs]==as.character(idd[i]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
          }else{
            idzAmt <- NA
          }
          tc   <- ncaId(ifdf,as.character(idd[i]))
          if(nrow(tc)==0) next
          time <- as.numeric(tc$time)
          conc <- as.numeric(tc$conc)
          cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i])))
          NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
          outData <- rbind(outData, data.frame(ID=as.character(idd[i]),STRAT1=popStr1[s1],STRAT2=popStr2[s2],Dose=idzAmt,t(NCAprm)))
        }
        
        if(noPlot==FALSE){
          # DV plot
          figlbl    <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]))
          cdata$FCT <- figlbl
          gdr <- dv.plot(pdata=cdata,cunit=cunit,tunit=tunit)
          suppressMessages(suppressWarnings(grid.arrange(gdr)))
          ggr <- grid.grab()
          concplot[[length(concplot)+1]] <- ggr
          if (printOut==TRUE){
            fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
            eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=8,width=16,units=\"cm\",res=200)")))
            suppressMessages(suppressWarnings(grid.arrange(gdr)))
            dev.off()
          }
          
          # Obs hist plot
          plotData    <- subset(outData, STRAT1==popStr1[s1] & STRAT2==popStr2[s2], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
          if (nrow(plotData)<5) next
          pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
          for (p in 1:length(pltPrm)){ if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
          if (length(pltPrm) == 0) next
          figlbl      <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]))
          histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histobsgrob$gdr
          mylegend    <- histobsgrob$legend
          lheight     <- histobsgrob$lheight
          if (printOut==TRUE){
            fl <- paste0(usrdir,"/HistObs_",figlbl)
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)")))
            suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          ggr <- grid.grab()
          histobsplot[[length(histobsplot)+1]] <- ggr
        }
      }
    }
    cnm         <- c(popStrNm1,popStrNm2,paste0("Dose (",dunit,")"),"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 4
  if (case == 4){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        for (s3 in 1:npopStr3){
          cdata <- data.frame()
          ifdf <- indf[(indf[,popStrNm1]==as.character(popStr1[s1]) & indf[,popStrNm2]==as.character(popStr2[s2]) & indf[,popStrNm3]==as.character(popStr3[s3])),]
          if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
          if (nrow(ifdf) == 0){next}
          idd        <- unique(as.character(ifdf[,idCol]))
          if (!is.null(doseAmtNm)){
            doseAmount <- paste(as.numeric(unique(refdf[refdf[,popStrNm1]==as.character(popStr1[s1]) & refdf[,popStrNm2]==as.character(popStr2[s2]) & refdf[,popStrNm3]==as.character(popStr3[s3]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
          }else{
            doseAmount <- NA
          }
          # Description
          pddf  <- rbind(pddf, data.frame(a=as.character(popStr1[s1]), b=as.character(popStr2[s2]), c=as.character(popStr3[s3]), d=doseAmount, e=length(idd)))
          for (i in 1:length(idd)){
            if (!is.null(doseAmtNm)){
              idzAmt <- as.numeric(refdf[refdf[,popStrNm1]==as.character(popStr1[s1]) & refdf[,popStrNm2]==as.character(popStr2[s2]) & refdf[,popStrNm3]==as.character(popStr3[s3]) & refdf[,idNmObs]==as.character(idd[i]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
            }else{
              idzAmt <- NA
            }
            tc     <- ncaId(ifdf,as.character(idd[i]))
            if(nrow(tc)==0) next
            time   <- as.numeric(tc$time)
            conc   <- as.numeric(tc$conc)
            cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i])))
            NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
            outData <- rbind(outData, data.frame(ID=as.character(idd[i]),STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Dose=idzAmt,t(NCAprm)))
          }
          
          if(noPlot==FALSE){
            # DV plot
            figlbl    <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]),"_",popStrNm3,"-",as.character(popStr3[s3]))
            cdata$FCT <- figlbl
            gdr <- dv.plot(pdata=cdata,cunit=cunit,tunit=tunit)
            suppressMessages(suppressWarnings(grid.arrange(gdr)))
            ggr <- grid.grab()
            concplot[[length(concplot)+1]] <- ggr
            if (printOut==TRUE){
              fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
              eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=8,width=16,units=\"cm\",res=200)")))
              suppressMessages(suppressWarnings(grid.arrange(gdr)))
              dev.off()
            }
            
            # Obs hist plot
            plotData    <- subset(outData, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
            if (nrow(plotData)<5) next
            pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
            for (p in 1:length(pltPrm)){ if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
            if (length(pltPrm) == 0) next
            figlbl      <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]),"_",popStrNm3,"-",as.character(popStr3[s3]))
            histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histobsgrob$gdr
            mylegend    <- histobsgrob$legend
            lheight     <- histobsgrob$lheight
            if (printOut==TRUE){
              fl <- paste0(usrdir,"/HistObs_",figlbl)
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)")))
              suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
            ggr <- grid.grab()
            histobsplot[[length(histobsplot)+1]] <- ggr
          }
        }
      }
    }
    cnm         <- c(popStrNm1,popStrNm2,popStrNm3,paste0("Dose (",dunit,")"),"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Statistical analysis for each patient group
  stNm <- c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","No_points_Lambda_z","AUC_pBack_Ext_obs","AUC_pBack_Ext_pred","AUClower_upper","Rsq","Rsq_adjusted","Corr_XY","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","AUC_pExtrap_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred","AUC_pExtrap_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","MRTINF_obs","MRTINF_pred","Vss_obs","Vss_pred","Tau","Tmin","Cmin","Cavg","p_Fluctuation","Accumulation_Index","Clss")
  grStat <- data.frame()
  if (case == 1){
    pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
    nm <- data.frame(character(0))
    counter <- 1
    for (i in 1:length(stNm)){
      Nm <- stNm[i]
      tdf <- as.numeric(as.character(outData[outData[,Nm]!="NaN",Nm]))
      if (length(tdf) < 2){
        nm <- rbind(nm, data.frame(Nm))
        pm[counter,] <- rep(NA,13)
      }else{
        nm <- rbind(nm, data.frame(Nm))
        stPrm <- calc.stat(x=tdf)        # calls calc.stat function
        stPrm <- unname(stPrm)
        pm[counter,] <- stPrm
      }
      counter <- counter + 1
    }
    if (nrow(pm) == 0) next
    pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
    tmpStat <- t(cbind(nm, pm))
    rownames(tmpStat)[1] <- "Name"
    tmpStat <- cbind(Stat=rownames(tmpStat),tmpStat)
    for (cnum in 2:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
    tmpStat <- tmpStat[-1,]
    grStat <- rbind(grStat,tmpStat)
  }
  if (case == 2){
    for (s1 in 1:npopStr1){
      pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
      nm <- data.frame(character(0))
      counter <- 1
      for (i in 1:length(stNm)){
        Nm <- stNm[i]
        tdf <- as.numeric(as.character(outData[(outData$STRAT1==as.character(popStr1[s1]) & outData[,Nm]!="NaN"),Nm]))
        if (length(tdf) < 2){
          nm <- rbind(nm, data.frame(Nm))
          pm[counter,] <- rep(NA,13)
        }else{
          nm <- rbind(nm, data.frame(Nm))
          stPrm <- calc.stat(x=tdf)        # calls calc.stat function
          stPrm <- unname(stPrm)
          pm[counter,] <- stPrm
        }
        counter <- counter + 1
      }
      if (nrow(pm) == 0) next
      pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
      tmpStat <- t(cbind(nm, pm))
      rownames(tmpStat)[1] <- "Name"
      tmpStat <- cbind(STRAT1=as.character(popStr1[s1]),Stat=rownames(tmpStat),tmpStat)
      for (cnum in 3:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
      tmpStat <- tmpStat[-1,]
      grStat <- rbind(grStat,tmpStat)
    }
    names(grStat)[names(grStat)%in%"STRAT1"] <- popStrNm1
  }
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
        nm <- data.frame(character(0))
        counter <- 1
        for (i in 1:length(stNm)){
          Nm <- stNm[i]
          tdf <- as.numeric(as.character(outData[(outData$STRAT1==as.character(popStr1[s1]) & outData$STRAT2==as.character(popStr2[s2]) & outData[,Nm]!="NaN"),Nm]))
          if (length(tdf) < 2){
            nm <- rbind(nm, data.frame(Nm))
            pm[counter,] <- rep(NA,13)
          }else{
            nm <- rbind(nm, data.frame(Nm))
            stPrm <- calc.stat(x=tdf)        # calls calc.stat function
            stPrm <- unname(stPrm)
            pm[counter,] <- stPrm
          }
          counter <- counter + 1
        }
        if (nrow(pm) == 0) next
        pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
        tmpStat <- t(cbind(nm, pm))
        rownames(tmpStat)[1] <- "Name"
        tmpStat <- cbind(STRAT1=as.character(popStr1[s1]),STRAT2=as.character(popStr2[s2]),Stat=rownames(tmpStat),tmpStat)
        for (cnum in 4:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
        tmpStat <- tmpStat[-1,]
        grStat <- rbind(grStat,tmpStat)
      }
    }
    names(grStat)[names(grStat)%in%c("STRAT1","STRAT2")] <- c(popStrNm1,popStrNm2)
  }
  if (case == 4){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        for (s3 in 1:npopStr3){
          pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
          nm <- data.frame(character(0))
          counter <- 1
          for (i in 1:length(stNm)){
            Nm <- stNm[i]
            tdf <- as.numeric(as.character(outData[(outData$STRAT1==as.character(popStr1[s1]) & outData$STRAT2==as.character(popStr2[s2]) & outData$STRAT3==as.character(popStr3[s3]) & outData[,Nm]!="NaN"),Nm]))
            if (length(tdf) < 2){
              nm <- rbind(nm, data.frame(Nm))
              pm[counter,] <- rep(NA,13)
            }else{
              nm <- rbind(nm, data.frame(Nm))
              stPrm <- calc.stat(x=tdf)        # calls calc.stat function
              stPrm <- unname(stPrm)
              pm[counter,] <- stPrm
            }
            counter <- counter + 1
          }
          if (nrow(pm) == 0) next
          pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
          tmpStat <- t(cbind(nm, pm))
          rownames(tmpStat)[1] <- "Name"
          tmpStat <- cbind(STRAT1=as.character(popStr1[s1]),STRAT2=as.character(popStr2[s2]),STRAT3=as.character(popStr3[s3]),Stat=rownames(tmpStat),tmpStat)
          for (cnum in 5:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
          tmpStat <- tmpStat[-1,]
          grStat <- rbind(grStat,tmpStat)
        }
      }
    }
    names(grStat)[names(grStat)%in%c("STRAT1","STRAT2","STRAT3")] <- c(popStrNm1,popStrNm2,popStrNm3)
  }
  
  if(printOut==TRUE) write.table(grStat, file=paste0(usrdir,"/ObsStat-",outFileNm,".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  
  # Print output files if simulated data does not exist
  if (is.null(simFile)){
    # Raname ID and stratifier columns and format output table sigfig
    if(case == 1){
      names(outData)[names(outData)%in%c("ID")] <- c(idNmObs)
      outData[,c(2:ncol(outData))] <- as.data.frame(lapply(outData[,c(2:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    if(case == 2){
      names(outData)[names(outData)%in%c("ID","STRAT1")] <- c(idNmObs,popStrNm1)
      outData[,c(3:ncol(outData))] <- as.data.frame(lapply(outData[,c(3:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    if(case == 3){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2")] <- c(idNmObs,popStrNm1,popStrNm2)
      outData[,c(4:ncol(outData))] <- as.data.frame(lapply(outData[,c(4:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    if(case == 4){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2","STRAT3")] <- c(idNmObs,popStrNm1,popStrNm2,popStrNm3)
      outData[,c(5:ncol(outData))] <- as.data.frame(lapply(outData[,c(5:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    
    # Subset table to print in the report
    if(case == 1) {prnTab1 <- head(cbind(outData[,1:2], subset(outData, select = tabCol)),100); names(prnTab1)[1:2] <- names(outData)[1:2]}
    if(case == 2) {prnTab1 <- head(cbind(outData[,1:3], subset(outData, select = tabCol)),100); names(prnTab1)[1:3] <- names(outData)[1:3]}
    if(case == 3) {prnTab1 <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),100); names(prnTab1)[1:4] <- names(outData)[1:4]}
    if(case == 4) {prnTab1 <- head(cbind(outData[,1:5], subset(outData, select = tabCol)),100); names(prnTab1)[1:5] <- names(outData)[1:5]}
    
    # Add unit to report table header
    tabUnit1 <- data.frame(NAME=c("Dose"),
                           UNIT1=c(paste0(names(outData)[names(outData)%in%"Dose"],"\n(",dunit,")")),
                           UNIT2=c(paste0(names(outData)[names(outData)%in%"Dose"]," (",dunit,")")))
    
    tabUnit2 <- data.frame(NAME=c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","AUClower_upper","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMCINF_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred"),
                           UNIT1=c(paste0("C0\n(",cunit,")"),paste0("Tmax\n(",tunit,")"),paste0("Cmax\n(",cunit,")"),paste0("Cmax_D\n(",cunit,"/",dunit,")"),paste0("Tlast\n(",tunit,")"),paste0("Clast\n(",cunit,")"),paste0("AUClast\n(",aucunit,")"),paste0("AUMClast\n(",aumcunit,")"),paste0("MRTlast\n(",tunit,")"),paste0("AUClower_upper\n(",aucunit,")"),paste0("Lambda_z (/",tunit,")"),paste0("Lambda_z_lower\n(",tunit,")"),paste0("Lambda_z_upper\n(",tunit,")"),paste0("HL_Lambda_z\n(",tunit,")"),paste0("AUCINF_obs\n(",aucunit,")"),paste0("AUCINF_D_obs\n(",aucunit,"/",dunit,")"),paste0("Vz_obs\n(",vlunit,")"),paste0("Cl_obs\n(",clunit,")"),paste0("AUCINF_pred\n(",aucunit,")"),paste0("AUCINF_D_pred\n(",aucunit,"/",dunit,")"),paste0("Vz_pred\n(",vlunit,")"),paste0("Cl_pred\n(",clunit,")"),paste0("AUMCINF_obs\n(",aumcunit,")"),paste0("AUMCINF_pred\n(",aumcunit,")"),paste0("MRTINF_obs\n(",tunit,")"),paste0("MRTINF_pred\n(",tunit,")"),paste0("Tau\n(",tunit,")"),paste0("Tmin\n(",tunit,")"),paste0("Cmin\n(",cunit,")"),paste0("Cavg\n(",cunit,")"),paste0("AUCtau\n(",aucunit,")"),paste0("AUMCtau\n(",aumcunit,")"),paste0("Clss\n(",clunit,")"),paste0("Vss_obs\n(",vlunit,")"),paste0("Vss_pred\n(",vlunit,")")),
                           UNIT2=c(paste0("C0 (",cunit,")"),paste0("Tmax (",tunit,")"),paste0("Cmax (",cunit,")"),paste0("Cmax_D (",cunit,"/",dunit,")"),paste0("Tlast (",tunit,")"),paste0("Clast (",cunit,")"),paste0("AUClast (",aucunit,")"),paste0("AUMClast (",aumcunit,")"),paste0("MRTlast (",tunit,")"),paste0("AUClower_upper (",aucunit,")"),paste0("Lambda_z (/",tunit,")"),paste0("Lambda_z_lower (",tunit,")"),paste0("Lambda_z_upper (",tunit,")"),paste0("HL_Lambda_z (",tunit,")"),paste0("AUCINF_obs (",aucunit,")"),paste0("AUCINF_D_obs (",aucunit,"/",dunit,")"),paste0("Vz_obs (",vlunit,")"),paste0("Cl_obs (",clunit,")"),paste0("AUCINF_pred (",aucunit,")"),paste0("AUCINF_D_pred (",aucunit,"/",dunit,")"),paste0("Vz_pred (",vlunit,")"),paste0("Cl_pred (",clunit,")"),paste0("AUMCINF_obs (",aumcunit,")"),paste0("AUMCINF_pred (",aumcunit,")"),paste0("MRTINF_obs (",tunit,")"),paste0("MRTINF_pred (",tunit,")"),paste0("Tau (",tunit,")"),paste0("Tmin (",tunit,")"),paste0("Cmin (",cunit,")"),paste0("Cavg (",cunit,")"),paste0("AUCtau (",aucunit,")"),paste0("AUMCtau (",aumcunit,")"),paste0("Clss (",clunit,")"),paste0("Vss_obs (",vlunit,")"),paste0("Vss_pred (",vlunit,")")))
    tabUnit <- rbind(tabUnit1,tabUnit2)
    
    prnTab2 <- prnTab1
    names(prnTab1) <- unlist(lapply(names(prnTab1), FUN=function(x){if(x%in%tabUnit$NAME){x <- as.character(tabUnit[tabUnit$NAME==x,"UNIT1"])[1]}else{x}}))
    names(prnTab2) <- unlist(lapply(names(prnTab2), FUN=function(x){if(x%in%tabUnit$NAME){x <- as.character(tabUnit[tabUnit$NAME==x,"UNIT2"])[1]}else{x}}))
    
    
    # Add unit to output table header
    names(outData)[names(outData)%in%"Dose"] <- paste0(names(outData)[names(outData)%in%"Dose"]," (",dunit,")")
    names(outData)[names(outData)%in%c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","AUClower_upper","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMCINF_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred")] <- c(paste0("C0 (",cunit,")"),paste0("Tmax (",tunit,")"),paste0("Cmax (",cunit,")"),paste0("Cmax_D (",cunit,"/",dunit,")"),paste0("Tlast (",tunit,")"),paste0("Clast (",cunit,")"),paste0("AUClast (",aucunit,")"),paste0("AUMClast (",aumcunit,")"),paste0("MRTlast (",tunit,")"),paste0("AUClower_upper (",aucunit,")"),paste0("Lambda_z (/",tunit,")"),paste0("Lambda_z_lower (",tunit,")"),paste0("Lambda_z_upper (",tunit,")"),paste0("HL_Lambda_z (",tunit,")"),paste0("AUCINF_obs (",aucunit,")"),paste0("AUCINF_D_obs (",aucunit,"/",dunit,")"),paste0("Vz_obs (",vlunit,")"),paste0("Cl_obs (",clunit,")"),paste0("AUCINF_pred (",aucunit,")"),paste0("AUCINF_D_pred (",aucunit,"/",dunit,")"),paste0("Vz_pred (",vlunit,")"),paste0("Cl_pred (",clunit,")"),paste0("AUMCINF_obs (",aumcunit,")"),paste0("AUMCINF_pred (",aumcunit,")"),paste0("MRTINF_obs (",tunit,")"),paste0("MRTINF_pred (",tunit,")"),paste0("Tau (",tunit,")"),paste0("Tmin (",tunit,")"),paste0("Cmin (",cunit,")"),paste0("Cavg (",cunit,")"),paste0("AUCtau (",aucunit,")"),paste0("AUMCtau (",aumcunit,")"),paste0("Clss (",clunit,")"),paste0("Vss_obs (",vlunit,")"),paste0("Vss_pred (",vlunit,")"))
    
    streamsEnv <- parent.frame()
    if(exists("outData")) assign("ncaOutput",  outData,   envir=streamsEnv)
    if(exists("grStat"))  assign("ObsStat",    grStat,    envir=streamsEnv)
    
    if(printOut==TRUE){
      write.table(outData, file=paste0(usrdir,"/ncaOutput-",outFileNm,".tsv"), sep="\t", row.names=F, col.names=T, quote=F)                          # write the output in a file
      fnOut <- list(arglist=match.call(), case=case, TXT=txt, pddf=pddf, prnTab1=prnTab1, prnTab2=prnTab2, spread=spread, conc=concplot, histobs=histobsplot)        # Function output list
    }
  }else{
    ########################################################################
    ############## Analyze the simulated data if exists ####################
    ########################################################################
    od <- paste(usrdir,"/SIMDATA",sep="")
    if (file.exists(od)){
      
      if(is.null(overwrite_SIMDATA)){
        dirTest <- readline("Directory \"SIMDATA\" already exists.\n
                            Overwrite it? (type 1)\n
                            Rename the existing folder and create a new one? (type 2)\n
                            Use data from the existing folder? (type 3)\n")
      } else if(overwrite_SIMDATA==T){
        dirTest <- "1"
      } else if(overwrite_SIMDATA==F) {
        dirTest <- "3"
      }
      
      if (dirTest == "1"){
        unlink(od, recursive=T)
        dir.create(od)
      }else if (dirTest == "2"){
        print("\nRenaming \"SIMDATA\" to \"SIMDATA_PREVIOUS\"\n")
        file.rename(from="SIMDATA", to="SIMDATA_PREVIOUS")
        dir.create(od)
      }else if (dirTest == "3"){
        file_list <- list.files(path="./SIMDATA/", pattern="sim_[0-9]*.csv", full.names=T)
      }else{setwd(usrdir);stop("Bad choice!!!\n")}
    }else{
      dir.create(od)
      dirTest <- "0"
    }
    
    if (dirTest != 3){
      # read NONMEM output into individual simulation data file
      IPSIM <- function(table.sim,MDV.rm=T){
        
        #browser()
        #table.sim <-  "nca_simulation.1.npctab.dta"
        # this is faster but doesn't read correctly with the extra lines between simulations
        #library(readr)
        #sim_2 <- read_table(table.sim, col_names = TRUE,skip = 1)
        
        sim <- read.table(table.sim,skip=1,header=T,fill=T,as.is=T)
        #sim <- as.data.frame(apply(sim,2,as.numeric))
        sim <- as.data.frame(apply(sim,2,function(x) suppressWarnings(as.numeric(x))))
        Nro <- min(which(is.na(sim[,1])))-1
        sim <- sim[!is.na(sim[,1]),]
        sim$NSUB <- rep(1:(nrow(sim)/Nro),each=Nro)
        #nsim <- max(sim$NSUB)
        if(MDV.rm==T){
          if(any(colnames(sim)=='MDV')){sim <- sim[sim[,'MDV']==0,]
          }else{cat('\nWarning MDV data item not listed in header,
                         Could not remove dose events!')}
        }
        assign("nmdf", sim)
        return(nmdf)
      }
      
      #if(new_data_method && (requireNamespace("readr", quietly = TRUE) && (packageVersion("readr") >= "0.2.2") && (requireNamespace("dplyr", quietly = TRUE)))){
      if(new_data_method){
        nmdf <- read_nm_table(simFile,sim_num = T,sim_name="NSUB")
        nmdf <- data.frame(nmdf)
      } else {
        nmdf <- IPSIM(simFile,MDV.rm=F)
      }
      
      #       # compare methods
      #       sim_1 <- read_nm_sim(simFile, only_obs = F) 
      #       sim_2 <- IPSIM(simFile,MDV.rm=F)
      #       library(dplyr)
      #       all(summarise_each(sim_1,funs(mean,sd)) == summarise_each(sim_2,funs(mean,sd))) 
      #       all(dim(sim_1)==dim(sim_2))
      
      simID <- unique(nmdf$NSUB)
      nsim <- length(simID)
      
      if (printOut==TRUE) write.table(nmdf, file=paste0(usrdir,"/ncaSimData-",outFileNm,".tsv"), row.names=F, quote=F, sep="\t")
      
      # exclude data based on specific values on filter column (optional)
      if (!is.null(filterNm)){
        if(filterNm%in%colnames(nmdf)==T & !is.null(filterExcl)){
          # filterExcl  == values to be excluded
          filterCol <- which(colnames(nmdf) == filterNm)[1]
          for (i in 1:length(filterExcl)){
            if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i])){
              nmdf <- nmdf[nmdf[,filterCol] != filterExcl[i],]
            }else if(!grepl("[<>!=]", filterExcl[i])){
              nmdf <- nmdf[nmdf[,filterCol] != filterExcl[i],]
            }else{
              nmdf <- eval(parse(text=paste0("subset(nmdf, !",filterNm,"%in% nmdf[nmdf[,",filterCol,"]",filterExcl[i],",filterCol])")))
            }
          }
        }else{
          print("Note: Incorrect filterNm or filterExcl specification. filterNm will not be used to process the simulated data.\n")
        }
      }
      
      srdf <- nmdf[nmdf$NSUB == 1,]  # copy simulated data before processing
      
      if (idNmSim%in%colnames(nmdf)==F | timeNmSim%in%colnames(nmdf)==F | concNmSim%in%colnames(nmdf)==F){
        setwd(usrdir);stop("Incorrect column names of ID, TIME and/or DV in simulation output\n")
      }else{
        idCol   <- which(colnames(nmdf) == idNmSim)[1]
        timeCol <- which(colnames(nmdf) == timeNmSim)[1]
        concCol <- which(colnames(nmdf) == concNmSim)[1]
      }
      
      # 1st level population stratification
      if (!is.null(str1Nm)){
        if (str1Nm%in%colnames(srdf)==F){setwd(usrdir);stop("Incorrect name for the 1st level stratification column in simulation output\n")}
        if (is.null(str1)){str1 <- unique(sort(srdf[,str1Nm]))}
        nstr1 <- length(str1)
        #for (i in 1:nstr1){
        #  if (nrow(srdf[srdf[,str1Nm]==str1[i],]) == 0){setwd(usrdir);stop("1st level stratification ID does not match the values within 1st level stratification column in simulation output\n")}
        #}
      }
      
      # 2nd level population stratification
      if (!is.null(str2Nm)){
        if (str2Nm%in%colnames(srdf)==F){setwd(usrdir);stop("Incorrect name for the 2nd level stratification column in simulation output\n")}
        if (is.null(str2)){str2 <- unique(sort(srdf[,str2Nm]))}
        nstr2 <- length(str2)
        #for (i in 1:nstr2){
        #  if (nrow(srdf[srdf[,str2Nm]==str2[i],]) == 0){setwd(usrdir);stop("2nd level stratification ID does not match the values within 2nd level stratification column in simulation output\n")}
        #}
      }
      
      # 3rd level population stratification
      if (!is.null(str3Nm)){
        if (str3Nm%in%colnames(srdf)==F){setwd(usrdir);stop("Incorrect name for the 3rd level stratification column in simulation output\n")}
        if (is.null(str3)){str3 <- unique(sort(srdf[,str3Nm]))}
        nstr3 <- length(str3)
        #for (i in 1:nstr3){
        #  if (nrow(srdf[srdf[,str3Nm]==str3[i],]) == 0){setwd(usrdir);stop("3rd level stratification ID does not match the values within 3rd level stratification column in simulation output\n")}
        #}
      }
      
      # Dose amount is extracted from doseAmtNm column
      if (!is.null(doseAmtNm)){
        if (doseAmtNm%in%colnames(nmdf)==T){
          doseAmtNm <- doseAmtNm
        }else if ("AMT"%in%colnames(nmdf)){
          doseAmtNm <- "AMT"
        }else{
          doseAmtNm <- NULL
          print("Note: Dose amount column name provided in doseAmtNm or AMT column does not exist in the simulated data file. Dose related NCA metrics will not be estimated for the simulated data.\n")
        }
      }else{
        if ("AMT"%in%colnames(nmdf)) doseAmtNm <- "AMT"
      }
      
      
      # ignore data with BLQ = 1 or user specified value (optional)
      if (!is.null(blqNm)){
        if(blqNm%in%colnames(nmdf) == T){
          blqCol <- which(colnames(nmdf) == blqNm)[1]
          for (i in 1:length(blqExcl)){
            if(grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i])){
              nmdf <- nmdf[nmdf[,blqNm] != blqExcl[i],]
            }else if(!grepl("[<>!=]", blqExcl[i])){
              nmdf <- nmdf[nmdf[,blqNm] != blqExcl[i],]
            }else{
              nmdf <- eval(parse(text=paste0("subset(nmdf, !",blqNm,"%in% nmdf[nmdf[,",blqCol,"]",blqExcl[i],",blqCol])")))
            }
          }
        }else{
          print("Note: Incorrect BLQ column name. BLQ will not be used to process the simulated data.\n")
        }
      }
      
      #if (!is.null(blqNm)){
      #  if(blqNm%in%colnames(nmdf)==T){
      #    blqCol <- which(colnames(nmdf) == blqNm)[1]
      #    for (i in 1:length(blqExcl)) {nmdf <- nmdf[nmdf[,blqNm] != blqExcl[i],]}
      #  }else{
      #    print("Note: Incorrect BLQ column name in simulation output. BLQ will not be used to process the data.\n")
      #  }
      #}
      
      # include data based on specific values on EVID column (optional) but keep rows with TIME == 0
      if (evid == TRUE){
        if("EVID"%in%colnames(nmdf) == T){
          # uevid == unique values in EVID column
          # evidIncl == EVID values to be included
          # ievid == EVID values to be ignored
          uevid <- unique(as.numeric(as.character(nmdf$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
          if (length(ievid) != 0){for (i in 1:length(ievid)){nmdf <- nmdf[nmdf$EVID != ievid[i],]}}
        }else{
          print("Note: EVID column is not present. EVID will not be used to process the simulated data.\n")
        }
        #if("EVID"%in%colnames(nmdf)==T){
        #  # uevid == unique values in EVID column
        #  # evidIncl == EVID values to be included
        #  # ievid == EVID values to be ignored
        #  uevid <- unique(as.numeric(as.character(nmdf$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
        #  if (length(ievid) != 0){
        #    for (i in 1:length(ievid)){
        #      if (ievid[i] != 1){
        #        nmdf <- nmdf[nmdf$EVID != ievid[i],]
        #      }else{
        #        if (length(which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$EVID)) == as.numeric(ievid[i]))) == 0) next
        #        nmdf <- nmdf[-which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$EVID)) == as.numeric(ievid[i])),]
        #        #if (length(which(nmdf[,timeCol] != 0 & nmdf$EVID == ievid[i])) == 0) next
        #        #nmdf <- nmdf[-which(nmdf[,timeCol] != 0 & nmdf$EVID == ievid[i]),]
        #      }
        #    }
        #  }
        #}else{
        #  print("Note: Incorrect EVID column name in simulation output. EVID will not be used to process the simulated data.\n")
        #}
      }
      
      # if MDV fiter is present, exclude data for MDV == 1 but keep rows with TIME == 0
      if (mdv == TRUE){
        if("MDV"%in%colnames(nmdf) == T){
          nmdf <- nmdf[nmdf$MDV == 0,]
        }else{
          print("Note: MDV column is not present. MDV will not be used to process the simulated data.\n")
        }
      }
      #if (mdv == TRUE){
      #  if("MDV"%in%colnames(nmdf) == T){
      #    if (length(which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$MDV)) == 1)) == 0) next
      #    nmdf <- nmdf[-which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$MDV)) == 1),]
      #  }else{
      #    print("Note: Incorrect MDV column name in simulation output. MDV will not be used to process the simulated data.\n")
      #  }
      #}
      
      
      # cpopStrNm = Combined stratifyting column names
      # npopStr   = Number of stratification levels
      # popStrNm1 = 1st level stratifying column name
      # popStr1   = 1st level stratification ID names
      # npopStr1  = Number of 1st level stratification ID names
      if(is.null(str1Nm)  & is.null(str2Nm)  & is.null(str3Nm)) {case<-1; npopStr<-0} # No stratification
      if(!is.null(str1Nm) & is.null(str2Nm)  & is.null(str3Nm)) {case<-2; cpopStrNm<-str1Nm; npopStr<-1; popStrNm1<-str1Nm; popStr1<-str1; npopStr1<-length(str1)} # Str1
      if(is.null(str1Nm)  & !is.null(str2Nm) & is.null(str3Nm)) {case<-2; cpopStrNm<-str2Nm; npopStr<-1; popStrNm1<-str2Nm; popStr1<-str2; npopStr1<-length(str2)} # Str2
      if(is.null(str1Nm)  & is.null(str2Nm)  & !is.null(str3Nm)){case<-2; cpopStrNm<-str3Nm; npopStr<-1; popStrNm1<-str3Nm; popStr1<-str3; npopStr1<-length(str3)} # Str3
      
      # Str1 & Str2
      if(!is.null(str1Nm) & !is.null(str2Nm) & is.null(str3Nm)){
        case<-3; cpopStrNm<-paste(str1Nm,str2Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStr1<-str1; popStr2<-str2; npopStr1<-length(str1); npopStr2<-length(str2)
      }
      # Str1 & Str3
      if(!is.null(str1Nm) & is.null(str2Nm) & !is.null(str3Nm)){
        case<-3; cpopStrNm<-paste(str1Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str3Nm; popStr1<-str1; popStr2<-str3; npopStr1<-length(str1); npopStr2<-length(str3)
      }
      # Str2 & Str3
      if(is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
        case<-3; cpopStrNm<-paste(str2Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str2Nm; popStrNm2<-str3Nm; popStr1<-str2; popStr2<-str3; npopStr1<-length(str2); npopStr2<-length(str3)
      }
      # Str1 & Str2 & Str3
      if(!is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
        case<-4; cpopStrNm<-paste(str1Nm,str2Nm,str3Nm,sep=", ")
        npopStr<-3
        popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStrNm3<-str3Nm
        popStr1<-str1; popStr2<-str2; popStr3<-str3
        npopStr1<-length(str1); npopStr2<-length(str2); npopStr3<-length(str3)
      }
      
      # Calculate AUC parameters for the simulation output
      dset = "sim"
      # Function to extract time and conc data for NCA metrics calculation for simulated data
      simNcaId <- function(ifdf,ID){
        if(adminType == "iv-infusion" & is.null(TI)){
          amt  <- ifdf[ifdf[,idCol]==ID & ifdf$AMT > 0,"AMT"][1]
          rate <- ifdf[ifdf[,idCol]==ID & ifdf$RATE > 0,"RATE"][1]
          if (is.na(amt) | is.na(rate) | rate==0){setwd(usrdir);stop(paste("Incorrect AMT and/or RATE value IN NONMEM output for ",ID,sep=""))}else{TI <- amt/rate}
        }else{TI <- "NaN"}
        conc <- as.numeric(as.character(ifdf[ifdf[,idCol]==ID,concCol]))
        if (simLog == TRUE){
          for (c in 1:length(conc)){conc[c] <- ifelse ((conc[c] != 0), exp(conc[c]), 0)}
        }
        time   <- as.numeric(as.character(ifdf[ifdf[,idCol]==ID,timeCol]))
        
        tc <- data.frame(time,conc)
        if(nrow(tc)>0){
          tc <- tc[order(tc$time),]
          tc$time <- tc$time + abs(min(tc$time))
        }
        return(tc)
      }
      
      for (s in 1:nsim){
        simData <- data.frame()
        smdf    <- nmdf[nmdf$NSUB == simID[s],]
        # Calculate NCA parameters
        if (case == 1){
          ifdf <- smdf
          if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
          if (nrow(ifdf) == 0){next}
          idd <- unique(as.character(ifdf[,idCol]))
          for (i in 1:length(idd)){
            if (!is.null(doseAmtNm)){
              idzAmt <- as.numeric(srdf[srdf[,idNmSim]==as.character(idd[i]) & srdf[,doseAmtNm]!="." & as.numeric(as.character(srdf[,doseAmtNm])) > 0, doseAmtNm][1])
            }else{
              idzAmt <- NA
            }
            stc  <- simNcaId(ifdf,as.character(idd[i]))
            if(nrow(stc)==0) next
            time <- as.numeric(stc$time)
            conc <- as.numeric(stc$conc)
            NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
            simData <- rbind(simData, data.frame(ID=as.character(idd[i]),Dose=idzAmt,t(NCAprm),NSIM=s))
          }
        }
        if (case == 2){
          for (s1 in 1:npopStr1){
            ifdf <- smdf[smdf[,popStrNm1]==as.character(popStr1[s1]),]
            if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
            if (nrow(ifdf) == 0){next}
            idd <- unique(as.character(ifdf[,idCol]))
            if(is.data.frame(idd)) idd <- idd[[1]]
            for (i in 1:length(idd)){
              if (!is.null(doseAmtNm)){
                idzAmt <- as.numeric(srdf[srdf[,popStrNm1]==as.character(popStr1[s1]) & srdf[,idNmSim]==as.character(idd[i]) & srdf[,doseAmtNm]!="." & as.numeric(as.character(srdf[,doseAmtNm])) > 0, doseAmtNm][1])
              }else{
                idzAmt <- NA
              }
              stc  <- simNcaId(ifdf,as.character(idd[i]))
              if(nrow(stc)==0) next
              time <- as.numeric(stc$time)
              conc <- as.numeric(stc$conc)
              NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
              simData <- rbind(simData, data.frame(ID=as.character(idd[i]),STRAT1=popStr1[s1],Dose=idzAmt,t(NCAprm),NSIM=s))
            }
          }
        }
        if (case == 3){
          for (s1 in 1:npopStr1){
            for (s2 in 1:npopStr2){
              ifdf <- smdf[smdf[,popStrNm1]==as.character(popStr1[s1]) & smdf[,popStrNm2]==as.character(popStr2[s2]),]
              if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
              if (nrow(ifdf) == 0){next}
              idd <- unique(as.character(ifdf[,idCol]))
              for (i in 1:length(idd)){
                if (!is.null(doseAmtNm)){
                  idzAmt <- as.numeric(srdf[srdf[,popStrNm1]==as.character(popStr1[s1]) & srdf[,popStrNm2]==as.character(popStr2[s2]) & srdf[,idNmSim]==as.character(idd[i]) & srdf[,doseAmtNm]!="." & as.numeric(as.character(srdf[,doseAmtNm])) > 0, doseAmtNm][1])
                }else{
                  idzAmt <- NA
                }
                stc  <- simNcaId(ifdf,as.character(idd[i]))
                if(nrow(stc)==0) next
                time <- as.numeric(stc$time)
                conc <- as.numeric(stc$conc)
                NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
                simData <- rbind(simData, data.frame(ID=as.character(idd[i]),STRAT1=popStr1[s1],STRAT2=popStr2[s2],Dose=idzAmt,t(NCAprm),NSIM=s))
              }
            }
          }
        }
        if (case == 4){
          for (s1 in 1:npopStr1){
            for (s2 in 1:npopStr2){
              for (s3 in 1:npopStr3){
                ifdf <- smdf[smdf[,popStrNm1]==as.character(popStr1[s1]) & smdf[,popStrNm2]==as.character(popStr2[s2]) & smdf[,popStrNm3]==as.character(popStr3[s3]),]
                if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
                if (nrow(ifdf) == 0){next}
                idd <- unique(as.character(ifdf[,idCol]))
                for (i in 1:length(idd)){
                  if (!is.null(doseAmtNm)){
                    idzAmt <- as.numeric(srdf[srdf[,popStrNm1]==as.character(popStr1[s1]) & srdf[,popStrNm2]==as.character(popStr2[s2]) & srdf[,popStrNm3]==as.character(popStr3[s3]) & srdf[,idNmSim]==as.character(idd[i]) & srdf[,doseAmtNm]!="." & as.numeric(as.character(srdf[,doseAmtNm])) > 0, doseAmtNm][1])
                  }else{
                    idzAmt <- NA
                  }
                  stc  <- simNcaId(ifdf,as.character(idd[i]))
                  if(nrow(stc)==0) next
                  time <- as.numeric(stc$time)
                  conc <- as.numeric(stc$conc)
                  NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
                  simData <- rbind(simData, data.frame(ID=as.character(idd[i]),STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Dose=idzAmt,t(NCAprm),NSIM=s))
                }
              }
            }
          }
        }
        
        if (case == 1) names(simData)[names(simData)%in%c("ID")] <- c(idNmSim)
        if (case == 2) names(simData)[names(simData)%in%c("ID","STRAT1")] <- c(idNmSim,popStrNm1)
        if (case == 3) names(simData)[names(simData)%in%c("ID","STRAT1","STRAT2")] <- c(idNmSim,popStrNm1,popStrNm2)
        if (case == 4) names(simData)[names(simData)%in%c("ID","STRAT1","STRAT2","STRAT3")] <- c(idNmSim,popStrNm1,popStrNm2,popStrNm3)
        write.csv(simData, file=paste(od,"/sim_",s,".csv",sep=""), row.names=F, quote=F)
      }
    }
    
    # Statistical analysis for each individual
    setwd(od)
    
    # read all simulated NCA parameters to a list
    lasdf <- lapply(list.files(pattern="sim_[0-9]*.csv",full.names=T),function(i){read.csv(i, header=T)})
    nsim <- length(lasdf)
    
    # Rename the ID and stratification column names to ID, STRAT1, STRAT2 and/STRAT3
    if(case==1) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim),names(x))] <- c("ID"); return(x)})
    if(case==2) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim,popStrNm1),names(x))] <- c("ID","STRAT1"); return(x)})
    if(case==3) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim,popStrNm1,popStrNm2),names(x))] <- c("ID","STRAT1","STRAT2"); return(x)})
    if(case==4) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim,popStrNm1,popStrNm2,popStrNm3),names(x))] <- c("ID","STRAT1","STRAT2","STRAT3"); return(x)})
    
    dasdf <- do.call(rbind, lapply(lasdf, as.data.frame))
    if (printOut==TRUE) write.table(dasdf, file=paste0(usrdir,"/ncaSimEst-",outFileNm,".tsv"), row.names=F, quote=F, sep="\t")
    
    # Population histogram
    if (case == 1){
      smeanData <- data.frame()
      for (i in 1:length(lasdf)){
        tmdf   <- subset(data.frame(lasdf[[i]]), select=param)
        tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
        smeanData <- rbind(smeanData, tmpPrm)
      }
      if(noPlot==FALSE){
        obsdata     <- subset(outData, select=param, ID!="")
        figlbl      <- NULL
        histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histpopgrob$gdr
        mylegend    <- histpopgrob$legend
        lheight     <- histpopgrob$lheight
        if(printOut==TRUE){
          fl <- paste0(usrdir,"/popMean")
          eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
          suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
        ggr <- grid.grab()
        popplot[[length(popplot)+1]] <- ggr
      }
    }
    
    # Group based histogram for no multiple flag
    if (case == 2){
      for (s1 in 1:npopStr1){
        if (nrow(dasdf[dasdf$STRAT1==as.character(popStr1[s1]),]) == 0) next
        smeanData <- data.frame()
        for (i in 1:length(lasdf)){
          tmdf   <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==as.character(popStr1[s1]))
          tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
          smeanData <- rbind(smeanData, tmpPrm)
        }
        if(noPlot==FALSE){
          obsdata     <- subset(outData, select=param, ID!="" & STRAT1==as.character(popStr1[s1]))
          figlbl      <- paste0(popStrNm1,"-",as.character(popStr1[s1]))
          histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histpopgrob$gdr
          mylegend    <- histpopgrob$legend
          lheight     <- histpopgrob$lheight
          if (printOut==TRUE){
            fl <- paste0(usrdir,"/popMean_",figlbl)
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
          ggr <- grid.grab()
          popplot[[length(popplot)+1]] <- ggr
        }
      }
    }
    
    # Flag based histogram
    if (case == 3){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          if (nrow(dasdf[dasdf$STRAT1==as.character(popStr1[s1]) & dasdf$STRAT2==as.character(popStr2[s2]),]) == 0) next
          smeanData <- data.frame()
          for (i in 1:length(lasdf)){
            tmdf   <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]))
            tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
            smeanData <- rbind(smeanData, tmpPrm)
          }
          if(noPlot==FALSE){
            obsdata     <- subset(outData, select=param, ID!="" & STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]))
            figlbl      <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]))
            histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histpopgrob$gdr
            mylegend    <- histpopgrob$legend
            lheight     <- histpopgrob$lheight
            if (printOut==TRUE){
              fl <- paste0(usrdir,"/popMean_",figlbl)
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            ggr <- grid.grab()
            popplot[[length(popplot)+1]] <- ggr
          }
        }
      }
    }
    
    # Group and flag based histogram
    if (case == 4){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          for (s3 in 1:npopStr3){
            if (nrow(dasdf[dasdf$STRAT1==as.character(popStr1[s1]) & dasdf$STRAT2==as.character(popStr2[s2]) & dasdf$STRAT3==as.character(popStr3[s3]),]) == 0) next
            smeanData <- data.frame()
            for (i in 1:length(lasdf)){
              tmdf   <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]) & STRAT3==as.character(popStr3[s3]))
              tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
              smeanData <- rbind(smeanData, tmpPrm)
            }
            if(noPlot==FALSE){
              obsdata     <- subset(outData, select=param, ID!="" & STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]) & STRAT3==as.character(popStr3[s3]))
              figlbl      <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]),"_",popStrNm3,"-",as.character(popStr3[s3]))
              histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
              gdr         <- histpopgrob$gdr
              mylegend    <- histpopgrob$legend
              lheight     <- histpopgrob$lheight
              if (printOut==TRUE){
                fl <- paste0(usrdir,"/popMean_",figlbl)
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                dev.off()
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              ggr <- grid.grab()
              popplot[[length(popplot)+1]] <- ggr
            }
          }
        }
      }
    }
    
    devcol  <- paste0("d",param)
    npdecol <- paste0("npde",param)
    # ggplot options for the forest plot
    ggOpt_forest <- list(scale_color_manual(name="",values=c("mean"="red","SD"="darkgreen")),
                         theme(plot.title = element_text(size=10,face="bold"),
                               axis.title.x = element_text(size=9,face="bold"),
                               axis.title.y = element_text(size=9,face="bold"),
                               axis.text.x  = element_text(size=7,face="bold",color="black",angle=45,vjust=1,hjust=1),
                               axis.text.y  = element_text(size=7,face="bold",color="black",hjust=0),
                               legend.text  = element_text(size=9,face="bold"),
                               legend.background = element_rect(),
                               legend.position = "bottom", legend.direction = "horizontal",
                               legend.key.size = unit(0.8, "cm"),
                               panel.margin = unit(0.5, "cm"),
                               plot.margin  = unit(c(0.5,0.5,0.5,0.5), "cm")),
                         facet_wrap(~type, scales="free", ncol=2),
                         theme(strip.text.x = element_text(size=9, face="bold")))
    
    OTL   <- data.frame(No_of_outliers=numeric(0),ID_metric=character(0))
    npde  <- data.frame()
    fpval <- data.frame(type=character(0),mean=numeric(0),mcil=numeric(0),mciu=numeric(0),sdu=numeric(0),sducil=numeric(0),sduciu=numeric(0),str=character(0))
    
    if (case == 1){
      tdasdf <- dasdf
      id     <- unique(tdasdf$ID)
      pde    <- data.frame()
      metric <- ""
      nout   <- 0
      for (i in 1:length(id)){
        obsdata <- subset(outData, ID==id[i])
        simdata <- subset(tdasdf, ID==id[i])
        figlbl  <- NULL
        pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit,noPlot=noPlot)
        pde     <- rbind(pde, cbind(data.frame(ID=id[i]), pdeout$pde))
        outData[outData$ID==id[i],] <- pdeout$obsdata
        if (pdeout$metric != ""){
          nout     <- nout + 1
          metric   <- paste(metric,pdeout$metric,sep=", ")
          if(noPlot==FALSE){
            gdr      <- pdeout$grob
            mylegend <- pdeout$legend
            lheight  <- pdeout$lheight
            if (printOut==TRUE){
              fl <- paste0(usrdir,"/Outlier_ID-",id[i])
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            ggr <- grid.grab()
            outlierplot[[length(outlierplot)+1]] <- ggr
          }
        }
      }
      if (metric != "") metric <- gsub("^, ", "", metric)
      OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
      npde <- rbind(npde,pde)
      
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste0("npde",alwprm)
      for (r in 1:nrow(outData)){
        if (nrow(npde[npde$ID==outData$ID[r],])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,npdeNm] <- npde[npde$ID==outData$ID[r],npdeNm]
        }
      }
      
      plotdata <- outData
      if (nrow(plotdata) == 0) next
      figlbl <- "All data"
      # Deviation plot
      if(noPlot==FALSE){
        ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
        if (!is.null(ggdev)){
          suppressMessages(suppressWarnings(print(ggdev)))
          if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
          devplot[[length(devplot)+1]] <- ggdev
        }
      }
      # NPDE plot
      npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
      if (is.null(npdeout$forestdata)) next
      forestdata <- npdeout$forestdata
      forestdata$str <- figlbl
      fpval <- rbind(fpval, forestdata)
      if(noPlot==FALSE){
        npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
        suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
        if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
        
        histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
        suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
        if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/histNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
      }
      
      if(noPlot==FALSE){
        # Forest plot for NPDE
        fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
          geom_point(aes(sd,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
          geom_text(aes(label=signif(mean,2),x=mean,y=str,color="mean",vjust=-1),size=2,show_guide=F) +
          geom_text(aes(label=signif(mcil,2),x=mcil,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(mciu,2),x=mciu,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sd,2),x=sd,y=str,color="SD",vjust=1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdcil,2),x=sdcil,y=str,color="SD",vjust=2),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdciu,2),x=sdciu,y=str,color="SD",vjust=2),size=2,show_guide=F)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/forestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=200)))
      }
    }else if (case==2){
      for (s1 in 1:npopStr1){
        if (nrow(dasdf[dasdf$STRAT1==as.character(popStr1[s1]),]) == 0) next
        tdasdf <- subset(dasdf, STRAT1==as.character(popStr1[s1]))
        id     <- unique(tdasdf$ID)
        pde    <- data.frame()
        metric <- ""
        nout   <- 0
        figlbl <- paste0(popStrNm1,"-",as.character(popStr1[s1]))
        for (i in 1:length(id)){
          obsdata <- subset(outData, ID==id[i] & STRAT1==as.character(popStr1[s1]))
          simdata <- subset(tdasdf, ID==id[i])
          pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit)
          pde     <- rbind(pde, cbind(data.frame(ID=id[i],STRAT1=as.character(popStr1[s1])), pdeout$pde))
          outData[(outData$ID==id[i] & outData$STRAT1==as.character(popStr1[s1])),] <- pdeout$obsdata
          if (pdeout$metric != ""){
            nout     <- nout + 1
            metric   <- paste(metric,pdeout$metric,sep=", ")
            if(noPlot==FALSE){
              gdr      <- pdeout$grob
              mylegend <- pdeout$legend
              lheight  <- pdeout$lheight
              if (printOut==TRUE){
                fl <- paste0(usrdir,"/Outlier_ID-",id[i],"_",figlbl)
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                dev.off()
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              ggr <- grid.grab()
              outlierplot[[length(outlierplot)+1]] <- ggr
            }
          }
        }
        if (metric != "") metric <- gsub("^, ", "", metric)
        OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
        npde <- rbind(npde,pde)
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste0("npde",alwprm)
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,paste0("npde",alwprm)] <- npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r]),paste0("npde",alwprm)]
        }
      }
      
      for (s1 in 1:npopStr1){
        plotdata <- subset(outData, STRAT1==as.character(popStr1[s1]))
        if (nrow(plotdata) == 0) next
        figlbl <- paste0(popStrNm1,"-",as.character(popStr1[s1]))
        # Deviation plot
        if(noPlot==FALSE){
          ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
          if (!is.null(ggdev)){
            suppressMessages(suppressWarnings(print(ggdev)))
            if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
            devplot[[length(devplot)+1]] <- ggdev
          }
        }
        # NPDE plot
        npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
        if (is.null(npdeout$forestdata)) next
        forestdata <- npdeout$forestdata
        forestdata$str <- figlbl
        fpval <- rbind(fpval, forestdata)
        if(noPlot==FALSE){
          npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
          suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
          if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
          
          histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
          suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
          if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/histNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
        }
      }
      
      if(noPlot==FALSE){
        # Forest plot for NPDE
        fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
          geom_point(aes(sd,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
          geom_text(aes(label=signif(mean,2),x=mean,y=str,color="mean",vjust=-1),size=2,show_guide=F) +
          geom_text(aes(label=signif(mcil,2),x=mcil,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(mciu,2),x=mciu,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sd,2),x=sd,y=str,color="SD",vjust=1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdcil,2),x=sdcil,y=str,color="SD",vjust=2),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdciu,2),x=sdciu,y=str,color="SD",vjust=2),size=2,show_guide=F)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/forestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=200)))
      }
    }else if (case == 3){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          if (nrow(dasdf[dasdf$STRAT1==as.character(popStr1[s1]) & dasdf$STRAT2==as.character(popStr2[s2]),]) == 0) next
          tdasdf <- subset(dasdf, STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]))
          id     <- unique(tdasdf$ID)
          pde    <- data.frame()
          metric <- ""
          nout   <- 0
          for (i in 1:length(id)){
            obsdata <- subset(outData, ID==id[i] & STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]))
            simdata <- subset(tdasdf, ID==id[i])
            figlbl  <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr1[s2]))
            pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit)
            pde     <- rbind(pde, cbind(data.frame(ID=id[i],STRAT1=as.character(popStr1[s1]),STRAT2=as.character(popStr2[s2])), pdeout$pde))
            outData[(outData$ID==id[i] & outData$STRAT1==as.character(popStr1[s1]) & outData$STRAT2==as.character(popStr2[s2])),] <- pdeout$obsdata
            if (pdeout$metric != ""){
              nout     <- nout + 1
              metric   <- paste(metric,pdeout$metric,sep=", ")
              if(noPlot==FALSE){
                gdr      <- pdeout$grob
                mylegend <- pdeout$legend
                lheight  <- pdeout$lheight
                if (printOut==TRUE){
                  fl <- paste0(usrdir,"/Outlier_ID-",id[i],"_",figlbl)
                  eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  dev.off()
                }
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                ggr <- grid.grab()
                outlierplot[[length(outlierplot)+1]] <- ggr
              }
            }
          }
          if (metric != "") metric <- gsub("^, ", "", metric)
          OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
          npde <- rbind(npde,pde)
        }
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste0("npde",alwprm)
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,paste0("npde",alwprm)] <- npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r]),paste0("npde",alwprm)]
        }
      }
      
      for (s1 in 1:npopStr1){  
        for (s2 in 1:npopStr2){
          plotdata <- subset(outData, STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]))
          if (nrow(plotdata) == 0) next
          figlbl <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]))
          # Deviation plot
          if(noPlot==FALSE){
            ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
            if (!is.null(ggdev)){
              suppressMessages(suppressWarnings(print(ggdev)))
              if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
              devplot[[length(devplot)+1]] <- ggdev
            }
          }
          # NPDE plot
          npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
          if (is.null(npdeout$forestdata)) next
          forestdata <- npdeout$forestdata
          forestdata$str <- figlbl
          fpval <- rbind(fpval, forestdata)
          if(noPlot==FALSE){
            npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
            suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
            if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
            
            histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
            suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
            if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/histNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
          }
        }
      }
      
      if(noPlot==FALSE){
        # Forest plot for NPDE
        fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
          geom_point(aes(sd,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
          geom_text(aes(label=signif(mean,2),x=mean,y=str,color="mean",vjust=-1),size=2,show_guide=F) +
          geom_text(aes(label=signif(mcil,2),x=mcil,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(mciu,2),x=mciu,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sd,2),x=sd,y=str,color="SD",vjust=1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdcil,2),x=sdcil,y=str,color="SD",vjust=2),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdciu,2),x=sdciu,y=str,color="SD",vjust=2),size=2,show_guide=F)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/forestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=200)))
      }
    }else if (case == 4){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          for (s3 in 1:npopStr3){
            if (nrow(dasdf[dasdf$STRAT1==as.character(popStr1[s1]) & dasdf$STRAT2==as.character(popStr2[s2]) & dasdf$STRAT3==as.character(popStr3[s3]),]) == 0) next
            tdasdf <- subset(dasdf, STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]) & STRAT3==as.character(popStr3[s3]))
            id     <- unique(tdasdf$ID)
            pde    <- data.frame()
            metric <- ""
            nout   <- 0
            for (i in 1:length(id)){
              obsdata <- subset(outData, ID==id[i] & STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]) & STRAT3==as.character(popStr3[s3]))
              simdata <- subset(tdasdf, ID==id[i])
              if(nrow(obsdata)==0 | nrow(simdata)==0) next
              figlbl  <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]),"_",popStrNm3,"-",as.character(popStr3[s3]))
              pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit)
              pde     <- rbind(pde, cbind(data.frame(ID=id[i],STRAT1=as.character(popStr1[s1]),STRAT2=as.character(popStr2[s2]),STRAT3=as.character(popStr3[s3])), pdeout$pde))
              outData[(outData$ID==id[i] & outData$STRAT1==as.character(popStr1[s1]) & outData$STRAT2==as.character(popStr2[s2]) & outData$STRAT3==as.character(popStr3[s3])),] <- pdeout$obsdata
              if (pdeout$metric != ""){
                nout     <- nout + 1
                metric   <- paste(metric,pdeout$metric,sep=", ")
                if(noPlot==FALSE){
                  gdr      <- pdeout$grob
                  mylegend <- pdeout$legend
                  lheight  <- pdeout$lheight
                  if (printOut==TRUE){
                    fl <- paste0(usrdir,"/Outlier_ID-",id[i],"_",figlbl)
                    eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)")))
                    suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                    dev.off()
                  }
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  ggr <- grid.grab()
                  outlierplot[[length(outlierplot)+1]] <- ggr
                }
              }
            }
            if (metric != "") metric <- gsub("^, ", "", metric)
            OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
            npde <- rbind(npde,pde)
          }
        }
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste0("npde",alwprm)
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r] & npde$STRAT3==outData$STRAT3[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,paste0("npde",alwprm)] <- npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r] & npde$STRAT3==outData$STRAT3[r]),paste0("npde",alwprm)]
        }
      }
      
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          for (s3 in 1:npopStr3){
            plotdata <- subset(outData, STRAT1==as.character(popStr1[s1]) & STRAT2==as.character(popStr2[s2]) & STRAT3==as.character(popStr3[s3]))
            if (nrow(plotdata) == 0) next
            figlbl  <- paste0(popStrNm1,"-",as.character(popStr1[s1]),"_",popStrNm2,"-",as.character(popStr2[s2]),"_",popStrNm3,"-",as.character(popStr3[s3]))
            # Deviation plot
            if(noPlot==FALSE){
              ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
              if (!is.null(ggdev)){
                suppressMessages(suppressWarnings(print(ggdev)))
                if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
                devplot[[length(devplot)+1]] <- ggdev
              }
            }
            # NPDE plot
            npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
            if (is.null(npdeout$forestdata)) next
            forestdata <- npdeout$forestdata
            forestdata$str <- figlbl
            fpval <- rbind(fpval, forestdata)
            if(noPlot==FALSE){
              npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
              suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
              if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
              
              histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
              suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
              if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/histNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=200)))
            }
          }
        }
      }
      
      if(noPlot==FALSE){
        # Forest plot for NPDE
        fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
          geom_point(aes(sd,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
          geom_text(aes(label=signif(mean,2),x=mean,y=str,color="mean",vjust=-1),size=2,show_guide=F) +
          geom_text(aes(label=signif(mcil,2),x=mcil,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(mciu,2),x=mciu,y=str,color="mean",vjust=-1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sd,2),x=sd,y=str,color="SD",vjust=1.5),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdcil,2),x=sdcil,y=str,color="SD",vjust=2),size=2,show_guide=F) +
          geom_text(aes(label=signif(sdciu,2),x=sdciu,y=str,color="SD",vjust=2),size=2,show_guide=F)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut==TRUE) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/forestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=200)))
      }
    }
    
    # Statistical analysis for each patient group
    stNm <- c("Tmax","Cmax","AUClast","AUClower_upper","AUCINF_obs","AUC_pExtrap_obs","AUCINF_pred","AUC_pExtrap_pred","AUMClast","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","HL_Lambda_z","Rsq","Rsq_adjusted","No_points_Lambda_z")
    simGrStat <- data.frame()
    if (case == 1){
      pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
      nm <- data.frame(character(0))
      counter <- 1
      for (i in 1:length(stNm)){
        Nm <- stNm[i]
        tdf <- as.numeric(as.character(dasdf[dasdf[,Nm]!="NaN",Nm]))
        if (length(tdf) < 2){
          nm <- rbind(nm, data.frame(Nm))
          pm[counter,] <- rep(NA,13)
        }else{
          nm <- rbind(nm, data.frame(Nm))
          stPrm <- calc.stat(x=tdf)        # calls calc.stat function
          stPrm <- unname(stPrm)
          pm[counter,] <- stPrm
        }
        counter <- counter + 1
      }
      if (nrow(pm) == 0) next
      pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
      tmpStat <- t(cbind(nm, pm))
      rownames(tmpStat)[1] <- "Name"
      tmpStat <- cbind(Stat=rownames(tmpStat),tmpStat)
      for (cnum in 2:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
      tmpStat <- tmpStat[-1,]
      simGrStat <- rbind(simGrStat,tmpStat)
    }
    if (case == 2){
      simGrStat <- data.frame()
      for (s1 in 1:npopStr1){
        pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
        nm <- data.frame(character(0))
        counter <- 1
        for (i in 1:length(stNm)){
          Nm <- stNm[i]
          tdf <- as.numeric(as.character(dasdf[(dasdf$STRAT1==as.character(popStr1[s1]) & dasdf[,Nm]!="NaN"),Nm]))
          if (length(tdf) < 2){
            nm <- rbind(nm, data.frame(Nm))
            pm[counter,] <- rep(NA,13)
          }else{
            nm <- rbind(nm, data.frame(Nm))
            stPrm <- calc.stat(x=tdf)        # calls calc.stat function
            stPrm <- unname(stPrm)
            pm[counter,] <- stPrm
          }
          counter <- counter + 1
        }
        if (nrow(pm) == 0) next
        pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
        tmpStat <- t(cbind(nm, pm))
        rownames(tmpStat)[1] <- "Name"
        tmpStat <- cbind(STRAT1=as.character(popStr1[s1]),Stat=rownames(tmpStat),tmpStat)
        for (cnum in 3:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
        tmpStat <- tmpStat[-1,]
        simGrStat <- rbind(simGrStat,tmpStat)
      }
      names(simGrStat)[1] <- popStrNm1
    }
    if (case == 3){
      simGrStat <- data.frame()
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
          nm <- data.frame(character(0))
          counter <- 1
          for (i in 1:length(stNm)){
            Nm <- stNm[i]
            tdf <- as.numeric(as.character(dasdf[(dasdf$STRAT1==as.character(popStr1[s1]) & dasdf$STRAT2==as.character(popStr2[s2]) & dasdf[,Nm]!="NaN"),Nm]))
            if (length(tdf) < 2){
              nm <- rbind(nm, data.frame(Nm))
              pm[counter,] <- rep(NA,13)
            }else{
              nm <- rbind(nm, data.frame(Nm))
              stPrm <- calc.stat(x=tdf)        # calls calc.stat function
              stPrm <- unname(stPrm)
              pm[counter,] <- stPrm
            }
            counter <- counter + 1
          }
          if (nrow(pm) == 0) next
          pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
          tmpStat <- t(cbind(nm, pm))
          rownames(tmpStat)[1] <- "Name"
          tmpStat <- cbind(STRAT1=as.character(popStr1[s1]),STRAT2=as.character(popStr2[s2]),Stat=rownames(tmpStat),tmpStat)
          for (cnum in 4:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
          tmpStat <- tmpStat[-1,]
          simGrStat <- rbind(simGrStat,tmpStat)
        }
      }
      names(simGrStat)[c(1,2)] <- c(popStrNm1,popStrNm2)
    }
    if (case == 4){
      simGrStat <- data.frame()
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          for (s3 in 1:npopStr3){
            pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
            nm <- data.frame(character(0))
            counter <- 1
            for (i in 1:length(stNm)){
              Nm <- stNm[i]
              tdf <- as.numeric(as.character(dasdf[(dasdf$STRAT1==as.character(popStr1[s1]) & dasdf$STRAT2==as.character(popStr2[s2]) & dasdf$STRAT3==as.character(popStr3[s3]) & dasdf[,Nm]!="NaN"),Nm]))
              if (length(tdf) < 2){
                nm <- rbind(nm, data.frame(Nm))
                pm[counter,] <- rep(NA,13)
              }else{
                nm <- rbind(nm, data.frame(Nm))
                stPrm <- calc.stat(x=tdf)        # calls calc.stat function
                stPrm <- unname(stPrm)
                pm[counter,] <- stPrm
              }
              counter <- counter + 1
            }
            if (nrow(pm) == 0) next
            pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){signif(x,digits=4)}else{x}}))
            tmpStat <- t(cbind(nm, pm))
            rownames(tmpStat)[1] <- "Name"
            tmpStat <- cbind(STRAT1=as.character(popStr1[s1]),STRAT2=as.character(popStr2[s2]),STRAT3=as.character(popStr3[s3]),Stat=rownames(tmpStat),tmpStat)
            for (cnum in 5:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
            tmpStat <- tmpStat[-1,]
            simGrStat <- rbind(simGrStat,tmpStat)
          }
        }
      }
      names(simGrStat)[c(1:3)] <- c(popStrNm1,popStrNm2,popStrNm3)
    }
    
    if (printOut==TRUE) write.table(simGrStat, file=paste0(usrdir,"/SimStat-",outFileNm,".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
    
    # Print output table
    # Raname ID and stratifier columns and format output table sigfig
    if(case == 1){
      names(outData)[names(outData)%in%c("ID")] <- c(idNmObs)
      outData[,c(2:ncol(outData))] <- as.data.frame(lapply(outData[,c(2:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    if(case == 2){
      names(outData)[names(outData)%in%c("ID","STRAT1")] <- c(idNmObs,popStrNm1)
      outData[,c(3:ncol(outData))] <- as.data.frame(lapply(outData[,c(3:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    if(case == 3){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2")] <- c(idNmObs,popStrNm1,popStrNm2)
      outData[,c(4:ncol(outData))] <- as.data.frame(lapply(outData[,c(4:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    if(case == 4){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2","STRAT3")] <- c(idNmObs,popStrNm1,popStrNm2,popStrNm3)
      outData[,c(5:ncol(outData))] <- as.data.frame(lapply(outData[,c(5:ncol(outData))], FUN=function(x) format(round(as.numeric(x),4),nsmall=4)))
    }
    
    
    # Subset table to print in the report
    if(case == 1) prnTab1 <- head(cbind(outData[,1:2], subset(outData, select = tabCol)),100)
    if(case == 2) prnTab1 <- head(cbind(outData[,1:3], subset(outData, select = tabCol)),100)
    if(case == 3) prnTab1 <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),100)
    if(case == 4) prnTab1 <- head(cbind(outData[,1:5], subset(outData, select = tabCol)),100)
    
    # Add unit to output table header
    names(outData)[names(outData)%in%c("Dose","C0","Tmax","simTmax","dTmax","Cmax","simCmax","dCmax","Cmax_D","Tlast","Clast","AUClast","simAUClast","dAUClast","AUMClast","simAUMClast","dAUMClast","MRTlast","AUClower_upper","simAUClower_upper","dAUClower_upper","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","simHL_Lambda_z","dHL_Lambda_z","AUCINF_obs","simAUCINF_obs","dAUCINF_obs","AUCINF_D_obs","Vz_obs","Cl_obs","AUCINF_pred","simAUCINF_pred","dAUCINF_pred","AUCINF_D_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMCINF_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred")] <- c(paste0(names(outData)[names(outData)%in%"Dose"]," (",dunit,")"),paste0("C0 (",cunit,")"),paste0("Tmax (",tunit,")"),paste0("simTmax (",tunit,")"),paste0("dTmax (",tunit,")"),paste0("Cmax (",cunit,")"),paste0("simCmax (",cunit,")"),paste0("dCmax (",cunit,")"),paste0("Cmax_D (",cunit,"/",dunit,")"),paste0("Tlast (",tunit,")"),paste0("Clast (",cunit,")"),paste0("AUClast (",aucunit,")"),paste0("simAUClast (",aucunit,")"),paste0("dAUClast (",aucunit,")"),paste0("AUMClast (",aumcunit,")"),paste0("simAUMClast (",aumcunit,")"),paste0("dAUMClast (",aumcunit,")"),paste0("MRTlast (",tunit,")"),paste0("AUClower_upper (",aucunit,")"),paste0("simAUClower_upper (",aucunit,")"),paste0("dAUClower_upper (",aucunit,")"),paste0("Lambda_z (/",tunit,")"),paste0("Lambda_z_lower (",tunit,")"),paste0("Lambda_z_upper (",tunit,")"),paste0("HL_Lambda_z (",tunit,")"),paste0("simHL_Lambda_z (",tunit,")"),paste0("dHL_Lambda_z (",tunit,")"),paste0("AUCINF_obs (",aucunit,")"),paste0("simAUCINF_obs (",aucunit,")"),paste0("dAUCINF_obs (",aucunit,")"),paste0("AUCINF_D_obs (",aucunit,"/",dunit,")"),paste0("Vz_obs (",vlunit,")"),paste0("Cl_obs (",clunit,")"),paste0("AUCINF_pred (",aucunit,")"),paste0("simAUCINF_pred (",aucunit,")"),paste0("dAUCINF_pred (",aucunit,")"),paste0("AUCINF_D_pred (",aucunit,"/",dunit,")"),paste0("Vz_pred (",vlunit,")"),paste0("Cl_pred (",clunit,")"),paste0("AUMCINF_obs (",aumcunit,")"),paste0("AUMCINF_pred (",aumcunit,")"),paste0("MRTINF_obs (",tunit,")"),paste0("MRTINF_pred (",tunit,")"),paste0("Tau (",tunit,")"),paste0("Tmin (",tunit,")"),paste0("Cmin (",cunit,")"),paste0("Cavg (",cunit,")"),paste0("AUCtau (",aucunit,")"),paste0("AUMCtau (",aumcunit,")"),paste0("Clss (",clunit,")"),paste0("Vss_obs (",vlunit,")"),paste0("Vss_pred (",vlunit,")"))
    
    # Add unit to report table header
    tabUnit <- data.frame(NAME=c(c("Dose","C0","Tmax","simTmax","dTmax","Cmax","simCmax","dCmax","Cmax_D","Tlast","Clast","AUClast","simAUClast","dAUClast","AUMClast","simAUMClast","dAUMClast","MRTlast","AUClower_upper","simAUClower_upper","dAUClower_upper","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","simHL_Lambda_z","dHL_Lambda_z","AUCINF_obs","simAUCINF_obs","dAUCINF_obs","AUCINF_D_obs","Vz_obs","Cl_obs","AUCINF_pred","simAUCINF_pred","dAUCINF_pred","AUCINF_D_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMCINF_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred")),
                          UNIT1=c(paste0(names(outData)[names(outData)%in%"Dose"],"\n(",dunit,")"),paste0("C0\n(",cunit,")"),paste0("Tmax\n(",tunit,")"),paste0("simTmax\n(",tunit,")"),paste0("dTmax\n(",tunit,")"),paste0("Cmax\n(",cunit,")"),paste0("simCmax\n(",cunit,")"),paste0("dCmax\n(",cunit,")"),paste0("Cmax_D\n(",cunit,"/",dunit,")"),paste0("Tlast\n(",tunit,")"),paste0("Clast\n(",cunit,")"),paste0("AUClast\n(",aucunit,")"),paste0("simAUClast\n(",aucunit,")"),paste0("dAUClast\n(",aucunit,")"),paste0("AUMClast\n(",aumcunit,")"),paste0("simAUMClast\n(",aumcunit,")"),paste0("dAUMClast\n(",aumcunit,")"),paste0("MRTlast\n(",tunit,")"),paste0("AUClower_upper\n(",aucunit,")"),paste0("simAUClower_upper\n(",aucunit,")"),paste0("dAUClower_upper\n(",aucunit,")"),paste0("Lambda_z (/",tunit,")"),paste0("Lambda_z_lower\n(",tunit,")"),paste0("Lambda_z_upper\n(",tunit,")"),paste0("HL_Lambda_z\n(",tunit,")"),paste0("simHL_Lambda_z\n(",tunit,")"),paste0("dHL_Lambda_z\n(",tunit,")"),paste0("AUCINF_obs\n(",aucunit,")"),paste0("simAUCINF_obs\n(",aucunit,")"),paste0("dAUCINF_obs\n(",aucunit,")"),paste0("AUCINF_D_obs\n(",aucunit,"/",dunit,")"),paste0("Vz_obs\n(",vlunit,")"),paste0("Cl_obs\n(",clunit,")"),paste0("AUCINF_pred\n(",aucunit,")"),paste0("simAUCINF_pred\n(",aucunit,")"),paste0("dAUCINF_pred\n(",aucunit,")"),paste0("AUCINF_D_pred\n(",aucunit,"/",dunit,")"),paste0("Vz_pred\n(",vlunit,")"),paste0("Cl_pred\n(",clunit,")"),paste0("AUMCINF_obs\n(",aumcunit,")"),paste0("AUMCINF_pred\n(",aumcunit,")"),paste0("MRTINF_obs\n(",tunit,")"),paste0("MRTINF_pred\n(",tunit,")"),paste0("Tau\n(",tunit,")"),paste0("Tmin\n(",tunit,")"),paste0("Cmin\n(",cunit,")"),paste0("Cavg\n(",cunit,")"),paste0("AUCtau\n(",aucunit,")"),paste0("AUMCtau\n(",aumcunit,")"),paste0("Clss\n(",clunit,")"),paste0("Vss_obs\n(",vlunit,")"),paste0("Vss_pred\n(",vlunit,")")),
                          UNIT2=c(paste0(names(outData)[names(outData)%in%"Dose"]," (",dunit,")"),paste0("C0 (",cunit,")"),paste0("Tmax (",tunit,")"),paste0("simTmax (",tunit,")"),paste0("dTmax (",tunit,")"),paste0("Cmax (",cunit,")"),paste0("simCmax (",cunit,")"),paste0("dCmax (",cunit,")"),paste0("Cmax_D (",cunit,"/",dunit,")"),paste0("Tlast (",tunit,")"),paste0("Clast (",cunit,")"),paste0("AUClast (",aucunit,")"),paste0("simAUClast (",aucunit,")"),paste0("dAUClast (",aucunit,")"),paste0("AUMClast (",aumcunit,")"),paste0("simAUMClast (",aumcunit,")"),paste0("dAUMClast (",aumcunit,")"),paste0("MRTlast (",tunit,")"),paste0("AUClower_upper (",aucunit,")"),paste0("simAUClower_upper (",aucunit,")"),paste0("dAUClower_upper (",aucunit,")"),paste0("Lambda_z (/",tunit,")"),paste0("Lambda_z_lower (",tunit,")"),paste0("Lambda_z_upper (",tunit,")"),paste0("HL_Lambda_z (",tunit,")"),paste0("simHL_Lambda_z (",tunit,")"),paste0("dHL_Lambda_z (",tunit,")"),paste0("AUCINF_obs (",aucunit,")"),paste0("simAUCINF_obs (",aucunit,")"),paste0("dAUCINF_obs (",aucunit,")"),paste0("AUCINF_D_obs (",aucunit,"/",dunit,")"),paste0("Vz_obs (",vlunit,")"),paste0("Cl_obs (",clunit,")"),paste0("AUCINF_pred (",aucunit,")"),paste0("simAUCINF_pred (",aucunit,")"),paste0("dAUCINF_pred (",aucunit,")"),paste0("AUCINF_D_pred (",aucunit,"/",dunit,")"),paste0("Vz_pred (",vlunit,")"),paste0("Cl_pred (",clunit,")"),paste0("AUMCINF_obs (",aumcunit,")"),paste0("AUMCINF_pred (",aumcunit,")"),paste0("MRTINF_obs (",tunit,")"),paste0("MRTINF_pred (",tunit,")"),paste0("Tau (",tunit,")"),paste0("Tmin (",tunit,")"),paste0("Cmin (",cunit,")"),paste0("Cavg (",cunit,")"),paste0("AUCtau (",aucunit,")"),paste0("AUMCtau (",aumcunit,")"),paste0("Clss (",clunit,")"),paste0("Vss_obs (",vlunit,")"),paste0("Vss_pred (",vlunit,")")))
    
    prnTab2 <- prnTab1
    names(prnTab1) <- unlist(lapply(names(prnTab1), FUN=function(x){if(x%in%tabUnit$NAME){x <- as.character(tabUnit[tabUnit$NAME==x,"UNIT1"])}else{x}}))
    names(prnTab2) <- unlist(lapply(names(prnTab2), FUN=function(x){if(x%in%tabUnit$NAME){x <- as.character(tabUnit[tabUnit$NAME==x,"UNIT2"])}else{x}}))
    
    streamsEnv <- parent.frame()
    if(exists("outData"))   assign("ncaOutput",  outData,   envir=streamsEnv)
    if(exists("grStat"))    assign("ObsStat",    grStat,    envir=streamsEnv)
    if(exists("simGrStat")) assign("SimStat",    simGrStat, envir=streamsEnv)
    if(exists("nmdf"))      assign("ncaSimData", nmdf,      envir=streamsEnv)
    if(exists("dasdf"))     assign("ncaSimEst",  dasdf,     envir=streamsEnv)
    
    if (printOut==TRUE){
      # Create HTML output
      simFileNm <- ifelse(is.data.frame(simFile), deparse(substitute(simFile)), simFile)
      txt <- paste(txt,paste0("Name of the NONMEM simulation output file: \"",simFileNm,"\""),sep="\n")
      txt <- paste(txt,paste0("Number of simulations performed: ",nsim),sep="\n")
      pddf <- cbind(pddf,OTL); names(pddf)[c((ncol(pddf)-1),ncol(pddf))] <- c("No. of outlier","Selected outliers ID and NCA metrics")
      
      write.table(outData, file=paste0(usrdir,"/ncaOutput-",outFileNm,".tsv"), sep="\t", row.names=F, col.names=T, quote=F)   # write output table
      fnOut <- list(arglist=match.call(), case=case, TXT=txt, pddf=pddf, prnTab1=prnTab1, prnTab2=prnTab2, NSIM=nsim, spread=spread, conc=concplot, histobs=histobsplot, pop=popplot, dev=devplot, outlier=outlierplot, forest=forestplot, npde=npdeplot, histnpde=histnpdeplot, phth=phth, pwth=pwth)
    }
  }
  setwd(usrdir)
  
  if (printOut==TRUE){
    misc <- system.file("misc", package = "ncappc")
    
    if(is.null(outFileNm)) outFileNm <- obsFileNm
    if (is.null(simFile)){
      mdFile <- paste(misc,"ncappcReport-NCA.Rmd",sep="/")
      nwFile <- paste(misc,"ncappcReport-NCA.Rnw",sep="/")
      outNm  <- paste0("ncappcReport-NCA-",outFileNm,".tex")
    }else{
      mdFile <- paste(misc,"ncappcReport-NCA-PPC.Rmd",sep="/")
      nwFile <- paste(misc,"ncappcReport-NCA-PPC.Rnw",sep="/")
      outNm  <- paste0("ncappcReport-NCA-PPC-",outFileNm,".tex")
    }
    
    knit2html(input=mdFile, output=outNm, style=paste(misc,"custom.css",sep="/"))
    knit(input=nwFile, output=outNm)
    if (.Platform$OS.type == "unix"){
      texcomp <- system('which texi2pdf')
      if (texcomp == 0){
        knit2pdf(input=nwFile, output=outNm)
      }else{
        print("Please install \"texi2pdf\" to compile the produced tex file into a PDF report")
      }
    }else if (.Platform$OS.type == "windows"){
      texcomp <- system('kpsewhich pdftex --version')
      if (texcomp == 0){
        knit2pdf(input=nwFile, output=outNm)
      }else{
        print("Please install \"pdftex\" to compile the produced tex file into a PDF report")
      }
    }
  }
  #unlink(list.files(pattern = "ncappcReport.[a,t,l,m,o]"))
  #unlink(list.files(pattern = "sum.tex"))
  #unlink(list.files(pattern = "tab.tex"))
}
