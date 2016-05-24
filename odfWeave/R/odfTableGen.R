"odfTableGen" <-
function(x, dataType, header = NULL, tableName, styles,
         cgroup = NULL, rgroup = NULL)
{
   # Sanity check dataType
   legaltypes <- c('float', 'percentage', 'currency', 'date', 'time',
                   'boolean', 'string')
   if (! all(dataType %in% legaltypes))
   {
      stop('illegal value of dataType argument')
   }

   # Sanity check cgroup and rgroup
   if (!is.null(cgroup))
   {
      if (!is.data.frame(cgroup) || ncol(cgroup) < 2 || ncol(cgroup) > 4)
         stop('cgroup must be a data frame with 2 to 4 columns')
      if (nrow(cgroup) == 0)
         cgroup <- NULL
   }

   if (!is.null(rgroup))
   {
      if (!is.data.frame(rgroup) || ncol(rgroup) < 2 || ncol(rgroup) > 4)
         stop('rgroup must be a data frame with 2 to 4 columns')
      if (nrow(rgroup) == 0)
         rgroup <- NULL
   }

   if (!is.null(cgroup))
   {
      t.cgroup <- as.character(cgroup[[1]])
      n.cgroup <- cgroup[[2]]
      ts.cgroup <- if (ncol(cgroup) > 2) as.character(cgroup[[3]]) else styles$cgroupText
      cs.cgroup <- if (ncol(cgroup) > 3) as.character(cgroup[[4]]) else styles$cgroupCell

      if (!is.numeric(n.cgroup) || sum(n.cgroup) != ncol(x))
         stop('the second column of cgroup must be a numeric vector summing to ', ncol(x))
   } else {
      t.cgroup <- NULL
      n.cgroup <- NULL
      ts.cgroup <- NULL
      cs.cgroup <- NULL
   }

   if (!is.null(rgroup))
   {
      t.rgroup <- as.character(rgroup[[1]])
      n.rgroup <- rgroup[[2]]
      ts.rgroup <- if (ncol(rgroup) > 2) as.character(rgroup[[3]]) else styles$rgroupText
      cs.rgroup <- if (ncol(rgroup) > 3) as.character(rgroup[[4]]) else styles$rgroupCell

      if (!is.numeric(n.rgroup) || sum(n.rgroup) != nrow(x))
         stop('the second column of rgroup must be a numeric vector summing to ', nrow(x))

      if (!is.null(header))
         header <- c("", header)

      # Adjust t.cgroup/n.cgroup due to the new column if they were specified
      if(!is.null(cgroup))
      {
         t.cgroup <- c("", t.cgroup)
         n.cgroup <- c(1, n.cgroup)
         ts.cgroup <- c(ts.cgroup[1], ts.cgroup)
         cs.cgroup <- c(cs.cgroup[1], cs.cgroup)
      }
   } else {
      t.rgroup <- NULL
      n.rgroup <- NULL
      ts.rgroup <- NULL
      cs.rgroup <- NULL
   }

   # Generate the cell matrix, which contains the "table:table-cell" elements
   cellMatrix <- genCellMatrix(x, dataType, styles, t.rgroup, n.rgroup,
                               ts.rgroup, cs.rgroup)

   # Generate the "header" element, which contains a
   # "table:table-header-rows" element with children
   headLine <- genHeadLine(header, styles, t.cgroup, n.cgroup, ts.cgroup, cs.cgroup)

   # Generate the "start" element, which contains the "table:table"
   # start tag and a complete "table:table-column" tag
   has <- function(x) !is.null(x) && x != ""
   tableStyle <- if(has(styles$table))
      paste(" table:style-name=\"", styles$table, "\" ", sep = "")
   else
      ""
   startText <- paste(
      "\n<table:table table:name=\"",  tableName, "\" ", tableStyle, ">",
      "\n  <table:table-column ",
      "table:number-columns-repeated=\"", ncol(cellMatrix), "\"/>",
      sep = "")

   # Gather together all the pieces to return as a list
   list(
      start = startText,
      header = headLine,
      cells = cellMatrix,
      end = "\n</table:table>\n")
}

# Generate a character matrix of the same dimensions as "x".
# Each element of this matrix is a partial XML document, with
# a "table:table-cell" element at the root.
#
# To test it by itself, you can call it as follows:
#
#   library(odfWeave)
#   x <- matrix(1:12, nrow=3)
#   xChar <- format(x, digits=3)
#   dataType <- rep("double", ncol(x))
#   styles <- tableStyles(xChar, useRowNames=FALSE, colnames(x))
#   odfWeave:::genCellMatrix(xChar, dataType, styles, NULL, NULL, NULL, NULL)
#
"genCellMatrix" <-
function(x, dataType, styles, t.rgroup, n.rgroup, ts.rgroup, cs.rgroup)
{
   # Function to generate a matrix containing a single repeated value
   # with the specified dimensions
   makeMatrix <- function(text, dims) matrix(rep(text,prod(dims)),nrow=dims[1])

   endQuote <- makeMatrix("\" ", dim(styles$cell))

   # Wrap each element of "x" in "text:p" tags with appropriate
   # "text:style-name" attributes
   textNameStart <- makeMatrix(" text:style-name=\"", dim(styles$text))
   textNameStart <- ifelse(styles$text == "", "", textNameStart)
   textNameEnd <- ifelse(styles$text == "", "", endQuote)
   textName <- matrixPaste(textNameStart, styles$text, textNameEnd,
                           sep = c("", ""))
   textStart <- makeMatrix("      <text:p", dim(styles$text))
   textEnd <- makeMatrix(">", dim(styles$text))
   tagEnd <- makeMatrix("</text:p>\n", dim(styles$text))
   textMatrix <- matrixPaste(textStart, textName, textEnd, x, tagEnd,
                             sep = rep("", 5))  # XXX repeat 4 times?

   # Wrap each element of "textMatrix" in "table:table-cell" tags
   # with appropriate "table:style-name" and "office:value-type" attributes
   cellNameStart <- makeMatrix(" table:style-name=\"", dim(styles$cell))
   cellNameStart <- ifelse(styles$cell == "", "", cellNameStart)
   cellNameEnd <- ifelse(styles$cell == "", "", endQuote)
   cellName <- matrixPaste(cellNameStart, styles$cell, cellNameEnd,
                           sep = c("", ""))
   valueType <- matrix(rep(dataType, each = nrow(x)), nrow = nrow(x))
   valueTypeStart <- makeMatrix(" <table:table-cell office:value-type=\"",
                                dim(styles$cell))
   cellStart <- matrixPaste(valueTypeStart, valueType, endQuote,
                            sep = c("", ""))
   cellEnd <- makeMatrix(">\n", dim(styles$cell))
   tagEnd <- makeMatrix(" </table:table-cell>\n", dim(styles$cell))
   cellMatrix <- matrixPaste(cellStart, cellName, cellEnd, textMatrix, tagEnd)

   # Process the "rgroup" argument if specified
   if(!is.null(t.rgroup))
   {
      # Compute "firstCol" which is used to modify the first column of the
      # "cells" element of the return value
      # XXX Should the cell style be hardcoded like this?
      firstCol0 <- paste(
         "<table:table-cell",
         " table:style-name=\"", cs.rgroup, "\"",
         " table:number-rows-spanned=\"", n.rgroup, "\"",
         " office:value-type=\"string\">",
         "<text:p text:style-name=\"", ts.rgroup, "\">", t.rgroup, "</text:p>",
         "</table:table-cell>",
         sep = "")
      firstCol <- rep("<table:covered-table-cell/>", nrow(cellMatrix))
      idx <- c(1, 1 + cumsum(n.rgroup)[-length(n.rgroup)])
      firstCol[idx] <- firstCol0

      # Concatenate firstCol to cellMatrix
      cellMatrix <- cbind(firstCol, cellMatrix)
   }

   # Complete the generation of the "cells" element by wrapping
   # "cellMatrix" in row tags
   leftRowTags <- matrix(
      c(
         rep("<table:table-row>", nrow(cellMatrix)),
         rep("", (ncol(cellMatrix) - 1) * nrow(cellMatrix))),
      nrow = nrow(cellMatrix))

   rightRowTags <- matrix(
      c(
         rep("", (ncol(cellMatrix) - 1) * nrow(cellMatrix)),
         rep("</table:table-row>\n", nrow(cellMatrix))),
      nrow = nrow(cellMatrix))

   matrixPaste(leftRowTags, cellMatrix, rightRowTags, sep = rep("\n", 2))
}

"genHeadLine" <-
function(header, styles, t.cgroup, n.cgroup, ts.cgroup, cs.cgroup)
{
   # Compute the style attributes for the "table:table-cell" and
   # "text:p" elements
   if (!is.null(styles$headerCell))
   {
      cellHeaderStyle <- paste(" table:style-name=\"", styles$headerCell,
                               "\" ", sep = "")
      cellHeaderStyle <- ifelse(styles$headerCell == "", "", cellHeaderStyle)
   }

   if (!is.null(styles$header))
   {
      textHeaderStyle <- paste(" text:style-name=\"", styles$header, "\" ",
                               sep = "")
      textHeaderStyle <- ifelse(styles$header == "", "", textHeaderStyle)
   }

   if (!is.null(cs.cgroup))
   {
      cgroupCellStyle <- paste(" table:style-name=\"", cs.cgroup,
                               "\" ", sep = "")
      cgroupCellStyle <- ifelse(cs.cgroup == "", "", cgroupCellStyle)
   }

   if (!is.null(ts.cgroup))
   {
      cgroupTextStyle <- paste(" text:style-name=\"", ts.cgroup, "\" ",
                               sep = "")
      cgroupTextStyle <- ifelse(ts.cgroup == "", "", cgroupTextStyle)
   }

   # Compute a "table:table-row" element based on "t.cgroup" and "n.cgroup"
   # if specified
   if(!is.null(t.cgroup))
   {
      idx <- c(1, 1 + cumsum(n.cgroup)[-length(n.cgroup)])

      firstRow0 <- paste(
         "\n      <text:p ",
         cgroupTextStyle,
         ">",
         t.cgroup,
         "</text:p>",
         sep = "")

      firstRow01 <- paste(
         "<table:table-cell ",
         cgroupCellStyle,
         " table:number-columns-spanned=\"", n.cgroup, "\"",
         " office:value-type=\"string\">",
         firstRow0,
         "</table:table-cell>",
         sep = "")

      nc <- sum(n.cgroup)
      firstRow02 <- rep("<table:covered-table-cell/>", nc)
      firstRow02[idx] <- firstRow01

      firstRow <- paste(
         "\n    <table:table-row>",
         paste(firstRow02, collapse = ""),
         "\n    </table:table-row>",
         sep = "")
   } else {
      firstRow <- NULL
   }

   if (!is.null(header))
   {
      headLine01 <- paste(
         "\n      <text:p ",
         textHeaderStyle,
         ">",
         header,
         "</text:p>",
         sep = "")
      headLine02 <- paste(
         "\n    <table:table-cell ",
         cellHeaderStyle,
         "office:value-type=\"string\">",
         headLine01,
         "\n    </table:table-cell>",
         sep = "")
      headLine03 <- paste(
         "\n    <table:table-row>\n",
         paste(headLine02, collapse = ""),
         "\n    </table:table-row>",
         sep = "")
   } else {
      headLine03 <- NULL
   }

   # Generate and return a "table:table-header-rows" element if we
   # generated either a normal "header" row, or a "cgroup" row.
   if (!is.null(firstRow) || !is.null(headLine03))
   {
      paste(
         "\n   <table:table-header-rows>\n",
         firstRow,
         headLine03,
         "\n   </table:table-header-rows>\n",
         sep = "")
   } else {
      NULL
   }
}
