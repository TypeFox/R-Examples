#######################################################################################################
#               Written by Andrew R Sommerlot <andrewrs@vt.edu>, January 2016                         #
#######################################################################################################
#' Formats SWAT subbasin files to user specifications. Defaults support use of getSWATcfsr.
#' @param wd - String. location of the subbasin files to be formatted. Generally this is the TxtInOut folder generated when building a SWAT project.
#' @param outDir - String. The location to write the formatted subbasin files. By default, this is the same as the wd, and will overwrite exising subbasin files.
#' @param numPars - Integer. Number of measured weather parameters used in the model. Ranges from 1 to 5 with 1 as precipitation only, 2 as precipitation and temperature, 3 as precipitation, temperature, and solar radiation, 4 as precipitation, temperature, solar radiation, and relative humidity, and 5 as precipitaion, temperature, solar radiation, relative humididty, and wind speed. Default is 5, or all parameters are formatted as measured inputs.
#' @param basinCentroid - Logical. If TRUE, then all gage location flags are set to "1", meaning there is a single time series for each of the specified measured meterological inputs. Only use if a subbain center approximation of meteorological data is being used. Defaults to FALSE, or each subbasin has a corresponding time series in each measured meteorological input file.
#' @return returns formatted subbasin files in the outDir location. If outDir is not specified files are saved to the wd location and any exising files are overwritten.
#' @examples
#' \dontrun{
#' SWATsubGaged(wd = '~/SWAT_PROJECT_FOLDER/TxtInOut')
#' }
#' @export

SWATsubGage <- function(wd, outDir = '', numPars = 5, basinCentroid = FALSE){
  if(outDir == ''){
    outDir = wd
  }
  ##set the working directory to the TxtInOut foler (may be an example for now)
  setwd(wd)

  #list all files with .sub extention
  files = list.files(pattern = "sub")

  #remove the output.sub file from the list
  files = files[! files %in% "output.sub"]

  files = files[! files %in% "TxtInOutsub.dat"]

  files = files[! files %in% "sub.dat"]

  #get total number of files to modify
  subnum = length(files)

  ### test for FOR loop. If works, then set to loop through entire length of files
  # get file name
  for(i in 1:subnum){
    subfile = files[i]

    #read in file
    sub = readLines(subfile)

    #grab the first line
    headline = sub[1]

    #character gymnastics to get the subbasin number alone
    subvec = strsplit(headline,":")

    subvec2 = subvec[[1]][2]

    subvec3 = strsplit(subvec2," ")

    subnum.file = subvec3[[1]][2]

    #format lead space for gage name change
    subnum.format = formatC(subnum.file, digits = 0, width = 16, flag = ' ', format = 'f')
    if(basinCentroid == FALSE){
      if(numPars == 1){
        ##changes flags for just precip and temp to their respective subbasin number
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
      } else if(numPars == 2){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
      } else if(numPars == 3){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
        sub[9] =  paste(subnum.format,"    | ISGAGE: solar radiation gage data used in subbasin", sep = "")
      } else if(numPars == 4){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
        sub[9] =  paste(subnum.format,"    | ISGAGE: solar radiation gage data used in subbasin", sep = "")
        sub[10] =  paste(subnum.format,"    | IHGAGE: relative humidity gage data used in subbasin", sep = "")
      } else if(numPars == 5){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
        sub[9] =  paste(subnum.format,"    | ISGAGE: solar radiation gage data used in subbasin", sep = "")
        sub[10] =  paste(subnum.format,"    | IHGAGE: relative humidity gage data used in subbasin", sep = "")
        sub[11] =  paste(subnum.format,"    | IWGAGE: wind speed gage data used in subbasin", sep = "")
      } else(stop(''))
    }

    if(basinCentroid == TRUE){
      subnum.format = 1
      if(numPars == 1){
        ##changes flags for just precip and temp to their respective subbasin number
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
      } else if(numPars == 2){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
      } else if(numPars == 3){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
        sub[9] =  paste(subnum.format,"    | ISGAGE: solar radiation gage data used in subbasin", sep = "")
      } else if(numPars == 4){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
        sub[9] =  paste(subnum.format,"    | ISGAGE: solar radiation gage data used in subbasin", sep = "")
        sub[10] =  paste(subnum.format,"    | IHGAGE: relative humidity gage data used in subbasin", sep = "")
      } else if(numPars == 5){
        sub[7] = paste(subnum.format,"    | IRGAGE: precip gage data used in subbasin", sep = "")
        sub[8] = paste(subnum.format,"    | ITGAGE: temp gage data used in subbasin", sep = "")
        sub[9] =  paste(subnum.format,"    | ISGAGE: solar radiation gage data used in subbasin", sep = "")
        sub[10] =  paste(subnum.format,"    | IHGAGE: relative humidity gage data used in subbasin", sep = "")
        sub[11] =  paste(subnum.format,"    | IWGAGE: wind speed gage data used in subbasin", sep = "")
      }
    }

    #writes out the changes
    writeLines(sub, paste(outDir, subfile, sep = "/"))
  }
  cat(paste('Files written to', outDir, sep = ' '))
}
