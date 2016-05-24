#' Find the block information 
#'
#' Takes the transect data and works out how many blocks of a given
#' size (in segment terms) fit into each.
#'
#' @param block.size number of segments per block
#' @param data data used to build the model
#' @param name.su names of the sampling units (ie. transects)
#'
#' @return a \code{data.frame} with the following columns
#'    \tabular{ll}{name        \tab the sample unit name (e.g. transect 
#'                                  label) \cr
#'                 num.seg     \tab number of segments in that transect \cr
#'                 num.block   \tab number of blocks available\cr
#'                 start.block \tab block # for first block\cr
#'                 end.block   \tab block # for last block\cr
#'                 num.req     \tab number of blocks needed for the unit\cr
#'                }
block.info.per.su <- function(block.size,data,name.su){

  unit <- data.frame(name=name.su)
  num.su <- length(name.su)

  for (i in 1:num.su) {
    unit$num.seg[i] <- length(
                     data$sampling.unit[data$sampling.unit == unit$name[i]])

    # if the number of segments in the unit is bigger than the block size,
    if(unit$num.seg[i]>block.size){
      unit$num.block[i] <- unit$num.seg[i] - block.size + 1
    # if the unit is smaller than the block size...
    }else{
      unit$num.block[i] <- unit$num.seg[i]
    }

    # the starting block is the last start block, plus the number of segments
    # in that last block, unless it's 1, of course
    if(i == 1){
      unit$start.block[i] <- 1
    }else{
      unit$start.block[i] <- unit$start.block[i - 1] + unit$num.block[i - 1]
    }
  }
  unit$end.block <- cumsum(unit$num.block)

  # Get number of blocks required for each sampling unit  --
  #  samp units different sizes so need to get more observations than required 
  unit$num.req <- ceiling(unit$num.seg/block.size)
  unit$num.req[unit$num.req==0] <- 1

  return(unit)
}
