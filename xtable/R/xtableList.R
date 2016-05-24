### Function to create lists of tables
xtableList <- function(x, caption = NULL, label = NULL, align = NULL,
                       digits = NULL, display = NULL, ...) {
  if (is.null(digits)){
    digitsList <- vector("list", length(x))
  } else {
    if (!is.list(digits)){
      digitsList <- vector("list", length(x))
      for (i in 1:length(x)) digitsList[[i]] <- digits
    }
  }
  if (is.null(display)){
    displayList <- vector("list", length(x))
  } else {
    if (!is.list(display)){
      displayList <- vector("list", length(x))
      for (i in 1:length(x)) displayList[[i]] <- display
    }
  }
  xList <- vector("list", length(x))
  for (i in 1:length(x)){
    xList[[i]] <- xtable(x[[i]], caption = caption, label = label,
                         align = align, digits = digitsList[[i]],
                         display = displayList[[i]], ...)
    attr(xList[[i]], 'subheading') <- attr(x, 'subheadings')[[i]]
  }
  attr(xList, "message") <- attr(x, "message")
  attr(xList, "caption") <- caption
  attr(xList, "label") <- label
  class(xList) <- c("xtableList")
  return(xList)
}

print.xtableList <- function(x,
  type = getOption("xtable.type", "latex"),
  file = getOption("xtable.file", ""),
  append = getOption("xtable.append", FALSE),
  floating = getOption("xtable.floating", TRUE),
  floating.environment = getOption("xtable.floating.environment", "table"),
  table.placement = getOption("xtable.table.placement", "ht"),
  caption.placement = getOption("xtable.caption.placement", "bottom"),
  caption.width = getOption("xtable.caption.width", NULL),
  latex.environments = getOption("xtable.latex.environments", c("center")),
  tabular.environment = getOption("xtable.tabular.environment", "tabular"),
  size = getOption("xtable.size", NULL),
  hline.after = NULL,
  NA.string = getOption("xtable.NA.string", ""),
  include.rownames = getOption("xtable.include.rownames", TRUE),
  colnames.format = "single",
  only.contents = getOption("xtable.only.contents", FALSE),
  add.to.row = NULL,
  sanitize.text.function = getOption("xtable.sanitize.text.function", NULL),
  sanitize.rownames.function = getOption("xtable.sanitize.rownames.function",
                                         sanitize.text.function),
  sanitize.colnames.function = getOption("xtable.sanitize.colnames.function",
                                         sanitize.text.function),
  sanitize.subheadings.function =
    getOption("xtable.sanitize.subheadings.function",
              sanitize.text.function),
  sanitize.message.function =
    getOption("xtable.sanitize.message.function",
              sanitize.text.function),
  math.style.negative = getOption("xtable.math.style.negative", FALSE),
  math.style.exponents = getOption("xtable.math.style.exponents", FALSE),
  html.table.attributes = getOption("xtable.html.table.attributes", "border=1"),
  print.results = getOption("xtable.print.results", TRUE),
  format.args = getOption("xtable.format.args", NULL),
  rotate.rownames = getOption("xtable.rotate.rownames", FALSE),
  rotate.colnames = getOption("xtable.rotate.colnames", FALSE),
  booktabs = getOption("xtable.booktabs", FALSE),
  scalebox = getOption("xtable.scalebox", NULL),
  width = getOption("xtable.width", NULL),
  comment = getOption("xtable.comment", TRUE),
  timestamp = getOption("xtable.timestamp", date()),
  ...)
{
  ## Get number of rows for each table in list of tables
  nCols <- dim(x[[1]])[2]
  rowNums <- sapply(x, dim)[1,]
  combinedRowNums <- cumsum(rowNums)
  combined <- do.call(rbind, x)
  if (type == "latex"){
    ## Special treatment if using booktabs
    if (booktabs){
      tRule <- "\\toprule"
      mRule <- "\\midrule"
      bRule <- "\\bottomrule"
    } else {
      tRule <- "\\hline"
      mRule <- "\\hline"
      bRule <- "\\hline"
    }
    ## Sanitize subheadings if required
    if (!is.null(sanitize.subheadings.function)) {
      for (i in 1:length(x)){
        attr(x[[i]], 'subheading') <-
          sanitize.subheadings.function(attr(x[[i]], 'subheading'))
      }
    }
    ## Sanitize message if required
    if (!is.null(sanitize.message.function)) {
      xMessage <- attr(x, 'message')
      xMessage <- sapply(xMessage, sanitize.message.function)
      attr(x, 'message') <- xMessage
    }
    if (colnames.format == "single"){

      add.to.row <- list(pos = NULL, command = NULL)
      add.to.row$pos <- as.list(c(0, combinedRowNums[-length(x)],
                                  dim(combined)[1]))
      command <- sapply(x, attr, "subheading")

      for (i in 1:length(x)){
        if( !is.null(command[[i]]) ){
          add.to.row$command[i] <-
            paste0(mRule,"\n\\multicolumn{", nCols, "}{l}{",
                   command[[i]],
                   "}\\\\\n")
        } else {
          add.to.row$command[i] <- paste0(mRule, "\n")
        }
      }
      ## Changed at request of Russ Lenth
      ## add.to.row$command[1:length(x)] <-
      ##   paste0(mRule,"\n\\multicolumn{", nCols, "}{l}{", command, "}\\\\\n")
      
      if ( (booktabs) & length(attr(x, "message") > 0) ){
        attr(x, "message")[1] <-
          paste0("\\rule{0em}{2.5ex}", attr(x, "message")[1])
      }
      add.to.row$command[length(x) + 1] <-
        paste0("\n\\multicolumn{", nCols, "}{l}{",
               attr(x, "message"), "}\\\\\n",
               collapse = "")
      add.to.row$command[length(x) + 1] <-
        paste0(bRule, add.to.row$command[length(x) + 1])

      class(combined) <- c("xtableList", "data.frame")
      hline.after <- c(-1)
      include.colnames <- TRUE
    }

    ## Create headings for columns if multiple headings are needed
    if (colnames.format == "multiple"){
      if (is.null(sanitize.colnames.function)) {
        colHead <- names(x[[1]])
      } else {
        colHead <- sanitize.colnames.function(names(x[[1]]))
      }
      if (rotate.colnames) {
        colHead <- paste("\\begin{sideways}", colHead, "\\end{sideways}")
      }
      colHead <- paste0(colHead, collapse = " & ")
      if (include.rownames) {
        colHead <- paste0(" & ", colHead)
      }
      colHead <- paste0(tRule, "\n", colHead, " \\\\", mRule, "\n")
      add.to.row <- list(pos = NULL, command = NULL)
      add.to.row$pos <- as.list(c(0, c(combinedRowNums[1:length(x)])))
      command <- sapply(x, attr, "subheading")


      add.to.row$command[1] <-
        if( !is.null(command[[1]]) ){
          add.to.row$command[1] <-
            paste0("\n\\multicolumn{", nCols, "}{l}{",
                   command[[1]],
                   "}\\\\ \n", colHead, "\n")
        } else {
          add.to.row$command[1] <- colHead
        }

      for (i in 2:length(x)) {
        add.to.row$command[i] <-
          if( !is.null(command[[i]]) ) {
            paste0(bRule,
                   "\\\\ \n\\multicolumn{", nCols, "}{l}{",
                   command[[i]], "}",
                   "\\\\ \n",
                   colHead)
          } else {
            add.to.row$command[i] <- paste0("\n", colHead)
          }
      }
      
      ## Changed at request of Russ Lenth
      ## add.to.row$command[1] <-
      ##   paste0("\n\\multicolumn{", nCols, "}{l}{", command[1],
      ##          "}", " \\\\ \n",
      ##          colHead)
      ## add.to.row$command[2:length(x)] <-
      ##   paste0(bRule,
      ##          "\\\\ \n\\multicolumn{", nCols, "}{l}{",
      ##          command[2:length(x)], "}",
      ##          "\\\\ \n",
      ##          colHead)
      if ( (booktabs) & length(attr(x, "message") > 0) ){
        attr(x, "message")[1] <-
          paste0("\\rule{0em}{2.5ex}", attr(x, "message")[1])
      }
      add.to.row$command[length(x) + 1] <-
        paste0("\n\\multicolumn{", nCols, "}{l}{",
               attr(x, "message"), "}\\\\\n",
               collapse = "")
      add.to.row$command[length(x) + 1] <-
        paste0(bRule, add.to.row$command[length(x) + 1])

      class(combined) <- c("xtableList", "data.frame")
      hline.after <- NULL

      include.colnames <- FALSE
    }

    print.xtable(combined,
                 type = type,
                 floating = floating,
                 floating.environment = floating.environment,
                 table.placement = table.placement,
                 caption.placement = caption.placement,
                 caption.width = caption.width,
                 latex.environments = latex.environments,
                 tabular.environment = tabular.environment,
                 size = size,
                 hline.after = hline.after,
                 NA.string = NA.string,
                 include.rownames = include.rownames,
                 include.colnames = include.colnames,
                 only.contents = only.contents,
                 add.to.row = add.to.row,
                 sanitize.text.function = sanitize.text.function,
                 sanitize.rownames.function = sanitize.rownames.function,
                 sanitize.colnames.function = sanitize.colnames.function,
                 math.style.negative = math.style.negative,
                 math.style.exponents = math.style.exponents,
                 html.table.attributes = html.table.attributes,
                 print.results = print.results,
                 format.args = format.args,
                 rotate.rownames = rotate.rownames,
                 rotate.colnames = rotate.colnames,
                 booktabs = booktabs,
                 scalebox = scalebox,
                 width = width,
                 comment = comment,
                 timestamp = timestamp,
                 ...)
  } else {
    stop("print.xtableList not yet implemented for this type")
  }
}


### Uses xtableList
xtableLSMeans <- function(x, caption = NULL, label = NULL,
                          align = NULL, digits = NULL,
                          display = NULL, auto = FALSE,
                          ...){
  if (attr(x, "estName") == "lsmean"){
    xList <- split(x, f = x[, 2])
    for (i in 1:length(xList)){
      xList[[i]] <- as.data.frame(xList[[i]][, -2])
    }
    attr(xList, "subheadings") <-
      paste0(dimnames(x)[[2]][2], " = ", levels(x[[2]]))
    attr(xList, "message") <- c("", attr(x, "mesg"))
    xList <- xtableList(xList, caption = caption, label = label,
                        align = align, digits = digits,
                        display = display, auto = auto, ...)
  } else {
    xList <- x
    xList <- xtable.data.frame(xList, caption = caption, label = label,
                               align = align, digits = digits,
                               display = display, auto = auto, ...)
  }
  return(xList)
}
