#' Moving block bootstrap sampler
#'
#' Not usually used on its own, called from within \code{\link{dsm.var.movblk}}.
#'
#' @param  num.blocks.required number of blocks that we need.
#' @param  block.size number of segments per block.
#' @param  which.blocks which blocks should be sampled.
#' @param  dsm.data the \code{$data} element of the result of a call to \code{dsm}.
#' @param  unit.info result of calling \code{\link{block.info.per.su}}.
#' @param  n.units number of sampling units.
#'
#' @return vector of log-residuals
#'
#' @export
generate.mb.sample <- function(num.blocks.required, block.size, which.blocks,
                                dsm.data, unit.info, n.units){

  #bs <- NULL
  #bs$block <- which.blocks
  bs <- data.frame(block=which.blocks)
  # storage
  bs.response <- c()
  bs.data <- c()

  for(i in 1:num.blocks.required){
    ## find the sampling unit that the block is in
    # is the block a start or end point for sampling unit?
    j <- which(bs$block[i] == unit.info$start.block)
    if(length(j)==0){
      j <- which(bs$block[i] == unit.info$end.block)
    }
    # if not a start or end, then where is it?
    if(length(j)==0){
      find.block <- c(unit.info$start.block,bs$block[i])
      j <- which(sort(find.block)==bs$block[i])-1
    }

    bs$unit.name[i] <- as.character(unit.info$name[j])
    if(j == 1){
      bs$unit.block[i] <- bs$block[i]
    }else{
      bs$unit.block[i] <- bs$block[i] - unit.info$end.block[j - 1]
    }

    # pull out the data for this sampling unit
    x.unit <- data.frame(dsm.data[dsm.data$sampling.unit == bs$unit.name[i], ])

    # pull out the rows corresponding to this block
    # start.row is the block number and the end row is
    #  (block length) segments after that
    start.row <- bs$unit.block[i]
    end.row <- start.row + block.size - 1
    x.block <- data.frame(x.unit[start.row:end.row, ])

    # append this to the data
    bs.data <- rbind(bs.data, x.block)
  }

  # Now need to map this onto data vector of the same length
  # (ie chopping off unwanted bits of blocks)

  temp <- bs.data$log.resids

  # loop over the sample units
  for (j in 1:n.units) {
    # Get number of segments in the blocks
    tb <- unit.info$num.req[j] * block.size
    # grab enough residuals for this sample unit
    tran.response <- temp[1:unit.info$num.seg[j]]
    # remove all of the ones we sampled (i.e. if we over sampled
    # make sure that we get rid of them too)
    temp <- temp[-(1:tb)]

    # store the result
    bs.response <- c(bs.response,tran.response)
  }
  return(bs.response)
}
