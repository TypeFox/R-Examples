#' Import XSYG files to R
#'
#' Imports XSYG files produced by a Freiberg Instrument lexsyg reader into R.
#'
#' \bold{How does the import function work?}\cr\cr The function uses the
#' \code{\link{xml}} package to parse the file structure. Each sequence is
#' subsequently translated into an \code{\linkS4class{RLum.Analysis}}
#' object.\cr\cr
#'
#' \bold{General structure XSYG format}\cr\cr \code{<?xml?}\cr \code{
#' <Sample>}\cr \code{ <Sequence>}\cr \code{ <Record>}\cr \code{ <Curve
#' name="first curve" />}\cr \code{ <Curve name="curve with data">}\cr \code{
#' x0 , y0 ; x1 , y1 ; x2 , y2 ; x3 , y3}\cr \code{ </Curve>}\cr \code{
#' </Record>}\cr \code{ </Sequence>}\cr \code{ </Sample>}\cr\cr So far, each
#' XSYG file can only contain one \code{<Sample></Sample>}, but multiple
#' sequences. \cr\cr Each record may comprise several curves.\cr\cr
#'
#' \bold{TL curve recalculation}\cr
#'
#' On the FI lexsyg device TL curves are recorded as time against count values.
#' Temperature values are monitored on the heating plate and stored in a
#' separate curve (time vs. temperature). If the option
#' \code{recalculate.TL.curves = TRUE} is chosen, the time values for each TL
#' curve are replaced by temperature values.\cr
#'
#' Practically, this means combining two matrices (Time vs. Counts and Time vs.
#' Temperature) with different row numbers by their time values. Three cases
#' are considered:
#'
#' HE: Heating element\cr PMT: Photomultiplier tube\cr Interpolation is done
#' using the function \code{\link{approx}}\cr
#'
#' CASE (1): \code{nrow(matrix(PMT))} > \code{nrow(matrix(HE))} \cr
#'
#' Missing temperature values from the heating element are calculated using
#' time values from the PMT measurement.\cr
#'
#' CASE (2): \code{nrow(matrix(PMT))} < \code{nrow(matrix(HE))} \cr
#'
#' Missing count values from the PMT are calculated using time values from the
#' heating element measurement.\cr
#'
#' CASE (3): \code{nrow(matrix(PMT))} == \code{nrow(matrix(HE))} \cr
#'
#' A new matrix is produced using temperature values from the heating element
#' and count values from the PMT. \cr
#'
#' \emph{Note: Please note that due to the recalculation of the temperature
#' values based on values delivered by the heating element, it may happen that
#' mutiple count values exists for each temperature value and temperature
#' values may also decrease during heating, not only increase. }\cr
#'
#' \bold{Advanced file import}\cr
#'
#' To allow for a more efficient usage of the function, instead of single path to a file just
#' a directory can be passed as input. In this particular case the function tries to extract
#' all XSYG-files found in the directory and import them all. Using this option internally the function
#' constructs as list of the XSYG-files found in the directory. Please note no recursive detection
#' is supported as this may lead to endless loops.
#'
#' @param file \code{\link{character}} or \code{\link{list}} (\bold{required}): path and file name of the
#' XSYG file. If input is a \code{list} it should comprise only \code{character}s representing each valid
#' path and xsyg-file names. Alternatively the input character can be just a directory (path), in this case the
#' the function tries to detect and import all xsyg files found in the directory.
#'
#' @param recalculate.TL.curves \link{logical} (with default): if set to
#' \code{TRUE}, TL curves are returned as temperature against count values (see
#' details for more information) Note: The option overwrites the time vs. count
#' TL curve. Select \code{FALSE} to import the raw data delivered by the
#' lexsyg. Works for TL curves and spectra.
#'
#' @param fastForward \code{\link{logical}} (with default): if \code{TRUE} for a
#' more efficient data processing only a list of \code{RLum.Analysis} objects is returned.
#'
#' @param import \code{\link{logical}} (with default): if set to \code{FALSE}, only
#' the XSYG file structure is shown.
#'
#' @param pattern \code{\link{regex}} (with default): optional regular expression if \code{file} is
#' a link to a folder, to select just specific XSYG-files
#'
#' @param txtProgressBar \link{logical} (with default): enables \code{TRUE} or
#' disables \code{FALSE} the progression bar during import
#'
#' @return \bold{Using the option \code{import = FALSE}}\cr\cr A list
#' consisting of two elements is shown: \item{Sample}{\link{data.frame} with
#' information on file.} \item{Sequences}{\link{data.frame} with information on
#' the sequences stored in the XSYG file}.\cr\cr \bold{Using the option
#' \code{import = TRUE} (default)} \cr\cr A list is provided, the list elements
#' contain: \item{Sequence.Header}{\link{data.frame} with information on the
#' sequence.} \item{Sequence.Object}{\code{\linkS4class{RLum.Analysis}}
#' containing the curves.}
#'
#' @note This function is a beta version as the XSYG file format is not yet
#' fully specified. Thus, further file operations (merge, export, write) should
#' be done using the functions provided with the package \code{\link{xml}}.\cr
#'
#' \bold{So far, no image data import is provided!}\cr Corresponding values in
#' the XSXG file are skipped.
#'
#'
#' @section Function version: 0.5.4
#'
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#'
#' @seealso \code{\link{xml}}, \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Data.Curve}}, \code{\link{approx}}
#'
#'
#' @references Grehl, S., Kreutzer, S., Hoehne, M., 2013. Documentation of the
#' XSYG file format. Unpublished Technical Note. Freiberg, Germany \cr\cr
#'
#' \bold{Further reading} \cr\cr XML: \url{http://en.wikipedia.org/wiki/XML}
#'
#'
#' @keywords IO
#'
#' @aliases readXSYG2R
#'
#' @examples
#'
#'
#' ##(1) import XSYG file to R (uncomment for usage)
#'
#' #FILE <- file.choose()
#' #temp <- read_XSYG2R(FILE)
#'
#' ##(2) additional examples for pure XML import using the package XML
#' ##    (uncomment for usage)
#'
#'   ##import entire XML file
#'   #FILE <- file.choose()
#'   #temp <- XML::xmlRoot(XML::xmlTreeParse(FILE))
#'
#'   ##search for specific subnodes with curves containing 'OSL'
#'   #getNodeSet(temp, "//Sample/Sequence/Record[@@recordType = 'OSL']/Curve")
#'
#' ##(2) How to extract single curves ... after import
#' data(ExampleData.XSYG, envir = environment())
#'
#' ##grep one OSL curves and plot the first curve
#' OSLcurve <- get_RLum(OSL.SARMeasurement$Sequence.Object, recordType="OSL")[[1]]
#'
#' ##(3) How to see the structure of an object?
#' structure_RLum(OSL.SARMeasurement$Sequence.Object)
#'
#'
#' @export
read_XSYG2R <- function(
  file,
  recalculate.TL.curves = TRUE,
  fastForward = FALSE,
  import = TRUE,
  pattern = ".xsyg",
  txtProgressBar = TRUE
){

  # Self Call -----------------------------------------------------------------------------------
  # Option (a): Input is a list, every element in the list will be treated as file connection
  # with that many file can be read in at the same time
  # Option (b): The input is just a path, the function tries to grep ALL xsyg/XSYG files in the
  # directory and import them, if this is detected, we proceed as list

  if(is(file, "character")) {

    ##If this is not really a path we skip this here
    if (dir.exists(file) & length(dir(file)) > 0) {
      message("[read_XSYG2R()] Directory detected, trying to extract '*.xsyg' files ...\n")
      file <-
        as.list(paste0(file,dir(
          file, recursive = TRUE, pattern = pattern
        )))

    }

  }

  if (is(file, "list")) {
    temp.return <- lapply(1:length(file), function(x) {
      read_XSYG2R(
        file = file[[x]],
        recalculate.TL.curves = recalculate.TL.curves,
        fastForward = fastForward,
        import = import,
        txtProgressBar = txtProgressBar
      )
    })

    ##return
    if (fastForward) {
      return(unlist(temp.return, recursive = FALSE))

    }else{
      return(temp.return)

    }

  }


  # Consistency check -------------------------------------------------------


  ##check if file exists
  if(!file.exists(file)){

    warning("[read_XSYG2R()] Wrong file name or file does not exist, nothing imported!")
    return(NULL)
  }

  ##check if file is XML file
  if(tail(unlist(strsplit(file, split = "\\.")), 1) != "xsyg" &
     tail(unlist(strsplit(file, split = "\\.")), 1) != "XSYG" ){

    warning("[read_XSYG2R()] File is not of type 'XSYG', nothing imported!")
    return(NULL)

  }

  # (0) config --------------------------------------------------------------
  #version.supported <- c("1.0")

  #additional functions
  ##get curve values
  get_XSYG.curve.values <- function(curve.node){

    ##Four steps
    ##(1) split string to paris of xy-values
    ##(2) split string to xy-values itself
    ##(3) convert to numeric
    ##(4) transpose matrix
    curve.node <- t(
      sapply(
        strsplit(
          strsplit(
            XML::xmlValue(curve.node), split = ";", fixed = TRUE)[[1]],
          split = ",", fixed = TRUE),
        as.numeric))

  }

  get_XSYG.spectrum.values <- function(curve.node){

    ##1st grep wavelength table
    wavelength <- XML::xmlAttrs(curve.node)["wavelengthTable"]

    ##string split
    wavelength <- as.numeric(unlist(strsplit(wavelength, split = ";", fixed = TRUE)))

    ##2nd grep time values
    curve.node <- unlist(strsplit(XML::xmlValue(curve.node), split = ";", fixed = TRUE))
    curve.node <- unlist(strsplit(curve.node, split = ",", fixed = TRUE), recursive = FALSE)

    curve.node.time <- as.numeric(curve.node[seq(1,length(curve.node),2)])

    ##3rd grep count values
    curve.node.count <- as.character(curve.node[seq(2,length(curve.node),2)])

    ##remove from pattern...
    curve.node.count <- do.call("gsub", list(pattern="[[]|[]]", replacement=" ",
                                             x=curve.node.count))

    ##4th combine to spectrum matrix
    spectrum.matrix <- matrix(0,length(wavelength),length(curve.node.time))
    spectrum.matrix <- sapply(1:length(curve.node.time), function(x){

      as.numeric(unlist(strsplit(curve.node.count[x], "[|]")))

    })


    ##change row names (rows are wavelength)
    rownames(spectrum.matrix) <- round(wavelength, digits=3)

    ##change column names (columns are time/temp values)
    colnames(spectrum.matrix) <- round(curve.node.time, digits=3)


    return(spectrum.matrix)
  }

  # (1) Integrity tests -----------------------------------------------------

  ##parse XML tree using the package XML
  temp <- try(XML::xmlRoot(XML::xmlTreeParse(file, useInternalNodes = TRUE)), silent = TRUE)



  ##show error
  if(is(temp, "try-error") == TRUE){

    warning("[read_XSYG2R()] XML file not readable, nothing imported!)")
    return(NULL)

  }

  # (2) Further file processing ---------------------------------------------

  ##==========================================================================##
  ##SHOW STRUCTURE
  if(import == FALSE){

    ##sample information
    temp.sample <- as.data.frame(XML::xmlAttrs(temp), stringsAsFactors = FALSE)

    ##grep sequences files

    ##set data.frame
    temp.sequence.header <- data.frame(t(1:length(names(XML::xmlAttrs(temp[[1]])))),
                                       stringsAsFactors = FALSE)

    colnames(temp.sequence.header) <- names(XML::xmlAttrs(temp[[1]]))

    ##fill information in data.frame
    for(i in 1:XML::xmlSize(temp)){

      temp.sequence.header[i,] <- t(XML::xmlAttrs(temp[[i]]))

    }

    output <-  list(Sample = temp.sample, Sequences = temp.sequence.header)
    return(output)

  }else{


    ##==========================================================================##
    ##IMPORT XSYG FILE

    ##Display output
    message(paste0("[read_XSYG2R()]\n  Importing: ",file))

    ##PROGRESS BAR
    if(txtProgressBar == TRUE){
      pb <- txtProgressBar(min=0,max=XML::xmlSize(temp), char = "=", style=3)
    }

    ##loop over the entire sequence by sequence
    output <- lapply(1:XML::xmlSize(temp), function(x){

      ##read sequence header
      temp.sequence.header <- as.data.frame(XML::xmlAttrs(temp[[x]]), stringsAsFactors = FALSE)
      colnames(temp.sequence.header) <- ""

      ###-----------------------------------------------------------------------
      ##LOOP
      ##read records >> records are combined to one RLum.Analysis object
      temp.sequence.object <- unlist(lapply(1:XML::xmlSize(temp[[x]]), function(i){

        ##get recordType
        temp.sequence.object.recordType <- try(XML::xmlAttrs(temp[[x]][[i]])["recordType"],
                                               silent = TRUE)

        ##the XSYG file might be broken due to a machine error during the measurement, this
        ##control flow helps; if a try-error is observed NULL is returned
        if(!inherits(temp.sequence.object.recordType, "try-error")){

        ##correct record type in depending on the stimulator
        if(temp.sequence.object.recordType == "OSL"){

          if(XML::xmlAttrs(temp[[x]][[i]][[
            XML::xmlSize(temp[[x]][[i]])]])["stimulator"] == "ir_LED_850" |
            XML::xmlAttrs(temp[[x]][[i]][[
              XML::xmlSize(temp[[x]][[i]])]])["stimulator"] == "ir_LD_850"){

            temp.sequence.object.recordType  <- "IRSL"

          }

        }


        ##loop 3rd level
        lapply(1:XML::xmlSize(temp[[x]][[i]]), function(j){

          ##get values
          temp.sequence.object.curveValue <- temp[[x]][[i]][[j]]

          ##get curveType
          temp.sequence.object.curveType <- as.character(
            XML::xmlAttrs(temp[[x]][[i]][[j]])["curveType"])

          ##get detector
          temp.sequence.object.detector <- as.character(
            XML::xmlAttrs(temp[[x]][[i]][[j]])["detector"])


          ##get stimulator
          temp.sequence.object.stimulator <- as.character(
            XML::xmlAttrs(temp[[x]][[i]][[j]])["stimulator"])


          ##get parentID
          temp.sequence.object.parentID <- as.numeric(
            XML::xmlAttrs(temp[[x]][[i]][[j]])["partentID"])

          ##get additional information
          temp.sequence.object.info <- as.list(XML::xmlAttrs(temp.sequence.object.curveValue))

          ##add stimulator and detector and so on
          temp.sequence.object.info <- c(temp.sequence.object.info,
                                         partentID = temp.sequence.object.parentID,
                                         position = as.integer(as.character(temp.sequence.header["position",])),
                                         name = as.character(temp.sequence.header["name",]))




          ## TL curve recalculation ============================================
          if(recalculate.TL.curves == TRUE){

            ##TL curve heating values is stored in the 3rd curve of every set
            if(temp.sequence.object.recordType == "TL" && j == 1){

              #grep values from PMT measurement or spectrometer
              if("Spectrometer" %in% temp.sequence.object.detector == FALSE){

                temp.sequence.object.curveValue.PMT <- get_XSYG.curve.values(
                  temp[[x]][[i]][[j]])


                ##round values (1 digit is technical resolution of the heating element)
                temp.sequence.object.curveValue.PMT[,1] <- round(
                  temp.sequence.object.curveValue.PMT[,1], digits = 1)

                #grep values from heating element
                temp.sequence.object.curveValue.heating.element <- get_XSYG.curve.values(
                  temp[[x]][[i]][[3]])



              }else{

                temp.sequence.object.curveValue.spectrum <- get_XSYG.spectrum.values(
                  temp.sequence.object.curveValue)

                ##get time values which are stored in the row labels
                temp.sequence.object.curveValue.spectrum.time <- as.numeric(
                  colnames(temp.sequence.object.curveValue.spectrum))

                ##round values (1 digit is technical resolution of the heating element)
                temp.sequence.object.curveValue.spectrum.time <- round(
                  temp.sequence.object.curveValue.spectrum.time, digits = 1)

              }

              #grep values from heating element
              temp.sequence.object.curveValue.heating.element <- get_XSYG.curve.values(
                temp[[x]][[i]][[3]])


              if("Spectrometer" %in% temp.sequence.object.detector == FALSE){

                #reduce matrix values to values of the detection
                temp.sequence.object.curveValue.heating.element <-
                  temp.sequence.object.curveValue.heating.element[
                    temp.sequence.object.curveValue.heating.element[,1] >=
                      min(temp.sequence.object.curveValue.PMT[,1]) &
                      temp.sequence.object.curveValue.heating.element[,1] <=
                      max(temp.sequence.object.curveValue.PMT[,1]),]

              }else{

                #reduce matrix values to values of the detection
                temp.sequence.object.curveValue.heating.element <-
                  temp.sequence.object.curveValue.heating.element[
                    temp.sequence.object.curveValue.heating.element[,1] >=
                      min(temp.sequence.object.curveValue.spectrum.time) &
                      temp.sequence.object.curveValue.heating.element[,1] <=
                      max(temp.sequence.object.curveValue.spectrum.time),]

              }




              ## calculate corresponding heating rate, this makes only sense
              ## for linear heating, therefore is has to be the maximum value

              ##remove 0 values (not measured) and limit to peak
              heating.rate.values <- temp.sequence.object.curveValue.heating.element[
                temp.sequence.object.curveValue.heating.element[,2] > 0 &
                  temp.sequence.object.curveValue.heating.element[,2] <=
                  max(temp.sequence.object.curveValue.heating.element[,2]),]

              heating.rate <- (heating.rate.values[length(heating.rate.values[,2]), 2] -
                                 heating.rate.values[1,2])/
                (heating.rate.values[length(heating.rate.values[,1]), 1] -
                   heating.rate.values[1,1])


              ##round values
              heating.rate <- round(heating.rate, digits=1)

              ##add to info element
              temp.sequence.object.info <- c(temp.sequence.object.info,
                                             RATE = heating.rate)


              ##PERFORM RECALCULATION
              ##check which object contains more data

              if("Spectrometer" %in% temp.sequence.object.detector == FALSE){

                ##CASE (1)
                if(nrow(temp.sequence.object.curveValue.PMT) >
                   nrow(temp.sequence.object.curveValue.heating.element)){

                  temp.sequence.object.curveValue.heating.element.i <- approx(
                    x = temp.sequence.object.curveValue.heating.element[,1],
                    y = temp.sequence.object.curveValue.heating.element[,2],
                    xout = temp.sequence.object.curveValue.PMT[,1],
                    rule = 2)

                  temperature.values <-
                    temp.sequence.object.curveValue.heating.element.i$y

                  count.values <-
                    temp.sequence.object.curveValue.PMT[,2]

                  ##CASE (2)
                }else if((nrow(temp.sequence.object.curveValue.PMT) <
                          nrow(temp.sequence.object.curveValue.heating.element))){

                  temp.sequence.object.curveValue.PMT.i <- approx(
                    x = temp.sequence.object.curveValue.PMT[,1],
                    y = temp.sequence.object.curveValue.PMT[,2],
                    xout = temp.sequence.object.curveValue.heating.element[,1],
                    rule = 2)

                  temperature.values <-
                    temp.sequence.object.curveValue.heating.element[,2]

                  count.values <- temp.sequence.object.curveValue.PMT.i$y

                  ##CASE (3)
                }else{

                  temperature.values <-
                    temp.sequence.object.curveValue.heating.element[,2]

                  count.values <- temp.sequence.object.curveValue.PMT[,2]

                }

                ##combine as matrix
                temp.sequence.object.curveValue <- as.matrix(cbind(
                  temperature.values,
                  count.values))

                ##set curve identifier
                temp.sequence.object.info$curveDescripter <- "Temperature [\u00B0C]; Counts [a.u.]"


              }else{

                ##CASE (1) here different approach. in contrast to the PMT measurements, as
                ##         usually the resolution should be much, much lower for such measurements
                ##         Otherwise we would introduce some pseudo signals, as we have to
                ##         take care of noise later one

                if(length(temp.sequence.object.curveValue.spectrum.time) !=
                   nrow(temp.sequence.object.curveValue.heating.element)){

                  temp.sequence.object.curveValue.heating.element.i <- approx(
                    x = temp.sequence.object.curveValue.heating.element[,1],
                    y = temp.sequence.object.curveValue.heating.element[,2],
                    xout = temp.sequence.object.curveValue.spectrum.time,
                    rule = 2,
                    ties = -2)

                  temperature.values <-
                    temp.sequence.object.curveValue.heating.element.i$y

                  ##check for duplicated values and if so, increase this values
                  if(anyDuplicated(temperature.values)>0){

                    temperature.values[which(duplicated(temperature.values))] <-
                      temperature.values[which(duplicated(temperature.values))]+1

                    warning("read_XSYG2R()] Temperatures values are found to be duplicated and increased by 1 K")

                  }




                  ##CASE (2)  (equal)
                }else{

                  temperature.values <-
                    temp.sequence.object.curveValue.heating.element[,2]

                }

                ##reset values of the matrix
                colnames(temp.sequence.object.curveValue.spectrum) <- temperature.values
                temp.sequence.object.curveValue <- temp.sequence.object.curveValue.spectrum

                ##change curve descriptor
                temp.sequence.object.info$curveDescripter <- "Temperature [\u00B0C]; Wavelength [nm]; Counts [1/ch]"

              }




            }##endif
          }##endif recalculate.TL.curves == TRUE


          ##Set RLum.Data objects
          if("Spectrometer" %in% temp.sequence.object.detector == FALSE){

            if(is(temp.sequence.object.curveValue, "matrix") == FALSE){

              temp.sequence.object.curveValue <-
                get_XSYG.curve.values(temp.sequence.object.curveValue)

            }


            set_RLum(
              class = "RLum.Data.Curve",
              originator = "read_XSYG2R",
              recordType = paste(temp.sequence.object.recordType,
                                 " (", temp.sequence.object.detector,")",
                                 sep = ""),
              curveType = temp.sequence.object.curveType,
              data = temp.sequence.object.curveValue,
              info = temp.sequence.object.info)

          }else if("Spectrometer" %in% temp.sequence.object.detector == TRUE) {


            if(is(temp.sequence.object.curveValue, "matrix") == FALSE){

              temp.sequence.object.curveValue <-
                get_XSYG.spectrum.values(temp.sequence.object.curveValue)

            }


            set_RLum(
              class = "RLum.Data.Spectrum",
              originator = "read_XSYG2R",
              recordType = paste(temp.sequence.object.recordType,
                                 " (",temp.sequence.object.detector,")",
                                 sep = ""),
              curveType = temp.sequence.object.curveType,
              data = temp.sequence.object.curveValue,
              info = temp.sequence.object.info)

          }

        })

        }else{

         return(NULL)

        }##if-try condition

      }))


      ##if the XSYG file is broken we get NULL as list element
      if (!is.null(temp.sequence.object)) {
        ##set RLum.Analysis object
        temp.sequence.object <-  set_RLum(
          originator = "read_XSYG2R",
          class = "RLum.Analysis",
          records = temp.sequence.object,
          protocol = as.character(temp.sequence.header["protocol",1])
        )


        ##update progress bar
        if (txtProgressBar == TRUE) {
          setTxtProgressBar(pb, x)
        }


        ##merge output and return values
        if(fastForward){
          return(temp.sequence.object)

        }else{
          return(list(Sequence.Header = temp.sequence.header, Sequence.Object = temp.sequence.object))

        }

      }else{
        return(temp.sequence.object)

      }

    })##end loop for sequence list

    ##close ProgressBar
    if(txtProgressBar == TRUE){close(pb)}

    ##show output informatioj
    if(length(output[sapply(output, is.null)]) == 0){

      message(paste("\t >>",XML::xmlSize(temp), " sequence(s) loaded successfully.\n"), sep = "")

    }else{

      message(paste("\t >>",XML::xmlSize(temp), " sequence(s) in file.",
                XML::xmlSize(temp)-length(output[sapply(output, is.null)]), "sequence(s) loaded successfully. \n"), sep = "")

      warning(paste0(length(output[sapply(output, is.null)])), " incomplete sequence(s) removed.")

    }

    ##output
    invisible(output)

  }#end if

  ##get rid of the NULL elements (as stated before ... invalid files)
  return(output[!sapply(output,is.null)])

}
## ---- DEPRECATED GENERICS
# .Deprecated in package version 0.5.0
# .Defunct in 0.5.X
# Removed in 0.6.0
#' @noRd
#' @export
readXSYG2R <- function(...) {
  .Deprecated("read_XSYG2R")
  read_XSYG2R(...)
}
