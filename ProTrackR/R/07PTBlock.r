.validity.PTBlock <- function(object){
  if (class(object) != "matrix") stop("Not a valid PTBlock")
  if (typeof(object) != "list") stop("Not a valid PTBlock")
  if (any(unlist(lapply(object, class)) != "PTCell")) stop("Not a valid PTBlock")
  if (any(unlist(lapply(object, length)) != 1)) stop("Not a valid PTBlock")
  if (any(!unlist(lapply(object, validObject)))) stop("Not a valid PTBlock")
  return(T)
}

setGeneric("PTBlock", function(pattern, row, track) standardGeneric("PTBlock"))
#' Select and copy a range of PTCells into a PTBlock
#'
#' Select and copy a range of \code{\link{PTCell}}s from a
#' \code{\link{PTPattern}} into a \code{PTBlock}. This
#' allows a more flexible approach to select and modify
#' \code{\link{PTCell}}s and paste the modified cells back into
#' a \code{\link{PTPattern}}.
#'
#' Most objects in this \link{ProTrackR} package are very strict in the operations
#' that are allowed, in order to guarantee validity and compatibility with
#' the original ProTracker. This makes those objects not very flexible.
#'
#' This \code{\link{PTBlock}} is not a formal S4 object, in fact you
#' can hardly call it an object at all. It is just a \code{matrix}, where each
#' element holds a \code{list} with a single \code{\link{PTCell}}.
#'
#' This \code{matrix} is very flexible and makes it easier to select and modify
#' the cells. This flexibility comes at a cost as validity is only checked
#' at the level of the \code{\link{PTCell}}s. The \code{PTBlock}
#' can be pasted back into a \code{\link{PTPattern}} with the
#' \code{\link{pasteBlock}} method. At which point validity will be checked again. If your modifications
#' resulted in violation of ProTracker standards, you should not be able to
#' paste the block into a pattern.
#'
#' @rdname PTBlock
#' @name PTBlock
#' @aliases PTBlock,PTPattern,numeric,numeric-method
#' @param pattern A \code{\link{PTPattern}} object from which the
#' \code{PTBlock} needs to be selected.
#' @param row A \code{numeric} index or indices of rows that needs to be
#' copied from the \code{pattern} into the PTBlock.
#' @param track A \code{numeric} index or indices of tracks that needs to be
#' copied from the \code{pattern} into the PTBlock.
#' @return Returns a \code{matrix} from the selected \code{row}s and \code{track}s
#' from the \code{pattern}. Each element in the \code{matrix} is a \code{list} holding
#' a single \code{\link{PTCell}}.
#' @examples
#' data("mod.intro")
#'
#' ## in most ProTrackR methods you can only select a single row or track.
#' ## with a PTBlock your selection is more flexible.
#'
#' ## select rows 4 up to 8 and tracks 2 up to 4, from the first
#' ## pattern table in mod.intro:
#'
#' block <- PTBlock(PTPattern(mod.intro, 1), 4:8, 2:4)
#'
#' ## 'block' is now a matrix with in each a list with a PTCell.
#' ## These can now easily be accessed and modified:
#'
#' cell1 <- block[1, 1][[1]]
#'
#' print(cell1)
#' @family block.operations
#' @author Pepijn de Vries
#' @export
setMethod("PTBlock", c("PTPattern", "numeric", "numeric"), function(pattern, row, track){
  cells <- apply(pattern@data, 1, function(x){
    index <- as.list(((1:maximumTrackCount)- 1)*4 + 1)
    lapply(index, function(y) PTCell(x[y:(y+3)]))
  })
  cells <- matrix(unlist(cells), 64, byrow = T)
  return(cells[row, track, drop = F])
})

setGeneric("pasteBlock", function(pattern, block, row.start, track.start) standardGeneric("pasteBlock"))

#' Paste a block of PTCell data into a PTPattern
#'
#' Paste a block of \code{\link{PTCell}} data into a \code{\link{PTPattern}} at
#' a specified location.
#'
#' A \code{\link{PTBlock}} is not a formal S4 class. It is a \code{matrix} where
#' each element holds a \code{list} of a single \code{\link{PTCell}} object. As
#' explained at the \code{\link{PTBlock}} method documentation, this allows for
#' a flexible approach of manipulating \code{\link{PTCell}} objects. The
#' \code{pasteBlock} method allows you to paste a \code{\link{PTBlock}} back into
#' a \code{\link{PTPattern}}.
#'
#' The \code{\link{PTBlock}} will be pasted at the specified location and will
#' span the number of tracks and rows that are included in the \code{\link{PTBlock}}.
#' The \code{\link{PTCell}}s in the \code{pattern} will be replaced by those
#' of the \code{block}. Elements of the \code{bock} that are out of the range
#' of the \code{pattern} are not included in the \code{pattern}.
#' @rdname pasteBlock
#' @name pasteBlock
#' @aliases pasteBlock,PTPattern,matrix,numeric,numeric-method
#' @param pattern A \code{\link{PTPattern}} object into which the \code{block}
#' needs to be pasted.
#' @param block A \code{\link{PTBlock}} holding the \code{\link{PTCell}} data
#' that needs to be pasted into the \code{pattern}.
#' @param row.start A positive \code{integer} value (ranging from 1 up to 64)
#' indicating the starting position (row) in the \code{pattern} to paste the
#' \code{block} into.
#' @param track.start A positive \code{integer} value (ranging from 1 up to 4)
#' indicating the starting position (track) in the \code{pattern} to paste the
#' \code{block} into.
#' @return Returns a copy of \code{pattern} into which \code{block} is pasted.
#' @examples
#' data("mod.intro")
#'
#' block <- PTBlock(PTPattern(mod.intro, 1), 1:16, 1)
#'
#' ## Do some operations using lapply (the effect
#' ## code is set to "C10"):
#' block <- matrix(lapply(block, function(x) {(effect(x) <- "C10"); x}),
#'                 nrow(block), ncol(block), byrow = TRUE)
#'
#' ## Paste block back on the same position:
#' PTPattern(mod.intro, 1) <-
#'   pasteBlock(PTPattern(mod.intro, 1), block, 1, 1)
#'
#' ## You can also paste the block anywhere you like:
#' PTPattern(mod.intro, 1) <-
#'   pasteBlock(PTPattern(mod.intro, 1), block, 49, 2)
#'
#' @family block.operations
#' @family pattern.operations
#' @author Pepijn de Vries
#' @export
setMethod("pasteBlock", c("PTPattern", "matrix", "numeric", "numeric"),
          function(pattern, block, row.start, track.start){
  row.start     <- abs(as.integer(row.start[[1]]))
  if (!(row.start %in% 1:maximumPatternTableRowCount)) stop("Invalid row starting position")
  track.start <- abs(as.integer(track.start[[1]]))
  if (!(track.start %in% 1:maximumTrackCount)) stop("Invalid row starting position")
  block         <- .PTBlock.as.raw(block)
  nrow.end      <- nrow(block)
  ntrack.end  <- ncol(block)
  if ((nrow.end + row.start) > maximumPatternTableRowCount)
    nrow.end <- maximumPatternTableRowCount - row.start + 1
  if ((ntrack.end/4 + track.start) > maximumTrackCount)
    ntrack.end <- 4*(maximumTrackCount - track.start + 1)

  print(nrow.end)
  print(ntrack.end)

  pattern@data[row.start:(row.start + nrow.end - 1),
               (track.start*4 - 3):(track.start*4 + ntrack.end - 4)] <-
    block[1:nrow.end, 1:ntrack.end]
  return(pattern)
})

.PTBlock.as.raw <- function(block)
{
  # test if the block is valid:
  .validity.PTBlock(block)
  return(matrix(unlist(lapply(t(block), as.raw)),
                nrow(block), 4*ncol(block), byrow = T))
}
