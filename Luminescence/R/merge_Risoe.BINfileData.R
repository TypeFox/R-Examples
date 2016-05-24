  #' Merge Risoe.BINfileData objects or Risoe BIN-files
#'
#' Function allows merging Risoe BIN/BINX files or Risoe.BINfileData objects.
#'
#' The function allows merging different measurements to one file or one
#' object.\cr The record IDs are recalculated for the new object. Other values
#' are kept for each object. The number of input objects is not limited. \cr
#'
#' \code{position.number.append.gap} option \cr
#'
#' If the option \code{keep.position.number = FALSE} is used, the position
#' numbers of the new data set are recalculated by adding the highest position
#' number of the previous data set to the each position number of the next data
#' set. For example: The highest position number is 48, then this number will
#' be added to all other position numbers of the next data set (e.g. 1 + 48 =
#' 49)\cr
#'
#' However, there might be cases where an additional addend (summand) is needed
#' before the next position starts. Example: \cr
#'
#' Position number set (A): \code{1,3,5,7}\cr Position number set (B):
#' \code{1,3,5,7} \cr
#'
#' With no additional summand the new position numbers would be:
#' \code{1,3,5,7,8,9,10,11}. That might be unwanted. Using the argument
#' \code{position.number.append.gap = 1} it will become:
#' \code{1,3,5,7,9,11,13,15,17}.
#'
#' @param input.objects \code{\link{character}} or
#' \code{\linkS4class{Risoe.BINfileData}} (\bold{required}): Character vector
#' with path and files names (e.g. \code{input.objects = c("path/file1.bin",
#' "path/file2.bin")} or \code{\linkS4class{Risoe.BINfileData}} objects (e.g.
#' \code{input.objects = c(object1, object2)})
#'
#'
#' @param output.file \code{\link{character}} (optional): File output path and
#' name. \cr If no value is given, a \code{\linkS4class{Risoe.BINfileData}} is
#' returned instead of a file.
#'
#'
#' @param keep.position.number \code{\link{logical}} (with default): Allows
#' keeping the original position numbers of the input objects. Otherwise the
#' position numbers are recalculated.
#'
#'
#' @param position.number.append.gap \code{\link{integer}} (with default): Set
#' the position number gap between merged BIN-file sets, if the option
#' \code{keep.position.number = FALSE} is used. See details for further
#' information.
#'
#'
#' @return Returns a \code{file} or a \code{\linkS4class{Risoe.BINfileData}}
#' object.
#'
#'
#' @note The validity of the output objects is not further checked.
#'
#'
#' @section Function version: 0.2.4
#'
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#'
#' @seealso \code{\linkS4class{Risoe.BINfileData}}, \code{\link{read_BIN2R}},
#' \code{\link{write_R2BIN}}
#'
#'
#' @references Duller, G., 2007. Analyst.
#'
#'
#' @keywords IO manip
#'
#'
#' @examples
#'
#'
#' ##merge two objects
#' data(ExampleData.BINfileData, envir = environment())
#'
#' object1 <- CWOSL.SAR.Data
#' object2 <- CWOSL.SAR.Data
#'
#' object.new <- merge_Risoe.BINfileData(c(object1, object2))
#'
#'
#' @export
merge_Risoe.BINfileData <- function(
  input.objects,
  output.file,
  keep.position.number = FALSE,
  position.number.append.gap = 0
){


  # Integrity Checks --------------------------------------------------------

  if(length(input.objects) < 2){

    stop("[merge_Risoe.BINfileData()] At least two input objects are needed!")

  }

  if(is(input.objects, "character") == TRUE){

    for(i in 1:length(input.objects)){

      if(file.exists(input.objects[i])==FALSE){

        stop("[merge_Risoe.BINfileData()] File",input.objects[i],"does not exists!")

      }

    }

  }else{

    if(is(input.objects, "list") == TRUE){

      for(i in 1:length(input.objects)){

        if(is(input.objects[[i]], "Risoe.BINfileData") == FALSE){

          stop("[merge_Risoe.BINfileData()] Input list does not contain Risoe.BINfileData objects!")

        }

      }

    }else{

      stop("[merge_Risoe.BINfileData()]
                Input object is not a 'character' nor a 'list'!")

    }

  }


  # Import Files ------------------------------------------------------------



  ##set temp object
  temp <- list()


  ##loop over all files to store the results in a list
  ##or the input is already a list

  if(is(input.objects, "character") == TRUE){
    for(i in 1:length(input.objects)){

      temp[i] <- read_BIN2R(input.objects[i])

    }

  }else{

    temp <- input.objects

  }

  # Get POSITION values -------------------------------------------------------

  ##grep maximum position value from the first file
  temp.position.max <- max(temp[[1]]@METADATA[, "POSITION"])

  ##grep all position values except from the first file
  temp.position.values <- unlist(sapply(2:length(temp), function(x){

    temp <- temp[[x]]@METADATA[, "POSITION"] +
      temp.position.max +
      position.number.append.gap

    temp.position.max <<- max(temp)

    return(temp)
  }))

  temp.position.values <- c(temp[[1]]@METADATA[, "POSITION"], temp.position.values)


  # Get overall record length -----------------------------------------------
  temp.record.length <- sum(sapply(1:length(temp), function(x){

    length(temp[[x]]@METADATA[,"ID"])

  }))


  # Merge Files -------------------------------------------------------------

  ##loop for similar input objects
  for(i in 1:length(input.objects)){

    if(exists("temp.new.METADATA") == FALSE){

      temp.new.METADATA <- temp[[i]]@METADATA
      temp.new.DATA <- temp[[i]]@DATA


      if(inherits(try(temp[[i]]@.RESERVED, silent = TRUE), "try-error")){

        temp.new.RESERVED <- list()

      }else{

        temp.new.RESERVED <- temp[[i]]@.RESERVED

      }

    }else{

      temp.new.METADATA <- rbind(temp.new.METADATA, temp[[i]]@METADATA)
      temp.new.DATA <- c(temp.new.DATA, temp[[i]]@DATA)

      if(inherits(try(temp[[i]]@.RESERVED, silent = TRUE), "try-error")){

        temp.new.RESERVED <- c(temp.new.RESERVED, list())

      }else{

        temp.new.RESERVED <- c(temp.new.RESERVED, temp[[i]]@.RESERVED)

      }

    }
  }


  ##SET RECORD ID in METADATA
  temp.new.METADATA$ID <- 1:temp.record.length

  ##SET POSITION VALUES
  if(keep.position.number == FALSE){

    temp.new.METADATA$POSITION <- temp.position.values

  }

  ##TODO version number?
  # Produce BIN file object -------------------------------------------------

  temp.new <- set_Risoe.BINfileData(
    METADATA = temp.new.METADATA,
    DATA = temp.new.DATA,
    .RESERVED = temp.new.RESERVED

  )



  # OUTPUT ------------------------------------------------------------------

  if(missing(output.file) == FALSE){

    write_R2BIN(temp.new, output.file)

  }else{

    return(temp.new)

  }

}
