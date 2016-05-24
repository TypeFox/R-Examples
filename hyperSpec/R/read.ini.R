##' Read INI files
##' 
##' \code{read.ini} reads ini files of the form
##'
##' [section]
##' key = value
##'
##' into a list.
##'
##' \code{read.ini} sanitizes the element names and tries to convert scalars and comma separated
##' numeric vectors to numeric.
##' @export
##' @rdname read-ini
##' @param con connection or file name
##' @author C. Beleites
##' @return a list with one element per section in the .ini file, each containing a list with elements
##' for the key-value-pairs.
##' @keywords IO file

read.ini <- function (con = stop ("Connection con needed.")){
  Lines  <- readLines(con)
  
  sections <- grep ("[[].*[]]", Lines)
  
  content <- Lines [- sections]
  ini <- as.list (gsub ("^.*=", "", content))
  names (ini) <- .sanitize.name (gsub ("=.*$", "", content))

  # try converting to numeric
  tmp <- lapply (ini, function (x) strsplit (x, ",") [[1]])
  tmp <- suppressWarnings (lapply (tmp, as.numeric))
  numbers <- ! sapply (tmp, function (x) any (is.na (x)))
  ini [numbers] <- tmp [numbers]
    
  tmp <- rep.int (seq_along (sections), diff (c (sections, length (Lines) + 1)) - 1)
  ini <- split (ini, tmp)
  
  sections <- Lines [sections]
  sections <- .sanitize.name (gsub ("^.(.*).$", "\\1", sections))
  names (ini) <- sections
  
  ini
}

.sanitize.name <- function (name){
  gsub ("[^a-zA-Z0-9._]", ".", name)
}
