###-----------------------------------------------------------------------------
###
### write.txt.wide
###
###
##' @param header.lines Toggle one or two line header (wavelengths in the
##'   second header line) for \code{write.txt.wide}
##' @aliases write.txt.wide
##' @rdname textio
##' @export

write.txt.wide <- function (object,
                            file = "",
                            cols = NULL,
                            quote = FALSE, sep = "\t",
                            row.names = FALSE,
                            col.names = TRUE,
                            header.lines = 1,   # 1 or 2 line header?
                                        # use labels instead of column names?
                            col.labels = if (header.lines == 1) FALSE else TRUE, 
                            append = FALSE,
                            ...){
  validObject (object)

  if (! is.null (cols))
    object <- object [, cols]

  if (col.names){
    col.spc <- match ("spc", colnames (object@data))

    if (col.labels){
      cln <- match (colnames (object@data), names (object@label))
      cln[!is.na (cln)] <- object@label [cln[!is.na(cln)]]
      cln[is.na (cln)] <- colnames (object@data) [is.na(cln)]
      cln <- sapply (cln, as.character)
                                        #cln [-col.spc] <- object@label []
    } else {
      cln <- colnames (object@data)
    }

    i <- seq_along (cln)

    if (header.lines == 1){
      write.table (matrix (c(if (row.names) "" else NULL,
                             cln [i < col.spc],
                             object@wavelength,
                             cln [i > col.spc]
                             ), nrow = 1),
                   file = file, append = append, quote = quote, sep = sep,
                   row.names = FALSE, col.names = FALSE)
      append = TRUE
    } else if (header.lines == 2) {
      ## 1st line
      write.table (matrix (c (
                              if (row.names) "" else NULL,
                              cln [i < col.spc],
                              if (col.labels) cln [col.spc] else "",
                              rep ("", length (object@wavelength) - 1),
                              cln [i > col.spc]), nrow = 1),
                   file = file, append = append, quote = quote, sep = sep,
                   row.names = FALSE, col.names = FALSE)
      append = TRUE 
      ## 2nd line
      write.table (matrix (c (if (row.names) (if (col.labels) as.character (object@label$.wavelength)
      else "wavelength")
      else NULL,
                              rep ("", sum (i < col.spc)),
                              object@wavelength,
                              rep ("", sum (i > col.spc))
                              ), nrow = 1),
                   file = file, append = append, quote, sep,
                   row.names = FALSE, col.names = FALSE)

    } else {
      stop ("Only 1 or 2 line headers supported.")
    }

  }

  # no AsIs columns!
  for (c in which (sapply (object@data, class) == "AsIs"))
    class (object@data [[c]]) <- NULL
   
  write.table (object@data, file = file, append = append, quote = quote, sep = sep,
               row.names = row.names, col.names = FALSE, ...)
}
