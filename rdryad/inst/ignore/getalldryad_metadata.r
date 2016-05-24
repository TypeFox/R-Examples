#' Download metadata for all Dryad oai's for defined time period.
#' 
#' @import RCurl XML stringr gdata plyr
#' @param transform (logical) transform metadata to list, TRUE or FLALSE
#' @param progress print progress bar (built in to the call to llply, plyr package)
#' @param write (logical) write metadata to local file, TRUE or FALSE
#' @param dir FALSE (default) or give directory as e.g. "/Mac/dryad/" only if
#'    write argument == TRUE
#' @return A Dryad dataset in a data.frame.
#' @export
#' @examples \dontrun{
#' mymetdata <- getalldryad_metadata(T, progress = "text", T,
#'    "/path/to/dir/")
#' }
getalldryad_metadata <- function(transform, progress = 'text', write, dir = FALSE)
{
  myoailist <- dr_listidentifiers('r') # get all oai's
  myoailist <- llply(myoailist[[1]], function(x) x$identifier) # list of file identifers only
        # NOTE: myoailist[[1]] would give the data packages instead of files
  allmetadat <- llply(  # download metadata for all oai's
    myoailist,
    download_dryadmetadata,
    transform = T,
    .progress = progress
  )
  allmetadat_ <- allmetadat[unlist(lapply(allmetadat, function(x) length(x$metadata)) != 0)]
  dryadtodf <- function(x) { # fxn to transform list to data.frame
    temp <- xmlToList(x$metadata)
    temp2 <- llply(temp, function(x) ifelse(is.null(x) == TRUE, "no entry", x))
    data.frame(temp2)
  }
  df <- ldply(allmetadat_, dryadtodf) # transform the list to data.frame
  if (write == "TRUE") {
    if (!dir == "FALSE") {
      file <- paste(dir, "dryadmetadata.csv", sep = "")
      write.csv(df, file = file)
      return(df)
    } else
    if (dir == "FALSE") {
      file <- "dryadmetadata.csv"
      write.csv(df, file = file)
      return(df)
    }
  } else
  if (write == "FALSE") {
    return(df)
  }
}