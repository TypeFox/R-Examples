writeListArraySDML <- function(x,
                               file,
                               textdata,
                               sep,
                               na.string,
                               null.string,
                               posinf.string,
                               neginf.string,
                               nan.string,
                               true,
                               false)
{
  if (is.null(x))
    catSDML("<empty/>\n", file = file)
  else if (is.list(x)) {
    catSDML("<list>\n", file = file)
    writeDimensionSDML(x, file = file)
    writePropertiesSDML(attributes(x), file = file, textdata = textdata,
           sep = sep, na.string = na.string, null.string = null.string,
           posinf.string = posinf.string, neginf.string = neginf.string,
           nan.string = nan.string, true = true, false = false)

    catSDML("<listdata>\n", file = file)
    lapply(x, writeListArraySDML, file = file, textdata = textdata,
           sep = sep, na.string = na.string, null.string = null.string,
           posinf.string = posinf.string, neginf.string = neginf.string,
           nan.string = nan.string, true = true, false = false)
    catSDML("</listdata>\n", file = file)
    catSDML("</list>\n", file = file)
  } else {
    catSDML("<array>\n", file = file)

    ## dimension tag
    writeDimensionSDML(x, file = file)

    ## prepare properties tag:
    ## check the datatype
    type <- if (inherits(x, "POSIXt"))
      "datetime"
    else if (is.character(x))
      "character"
    else if (is.factor(x))
      "categorical"
    else if (is.logical(x))
      "logical"
    else
      "numeric"

    ## check the mode
    mode <- NULL
    if (type == "numeric")
      mode <- if (is.integer(x))
        "integer"
      else if (is.double(x))
        "real"
      else
        "complex"
    if (type == "categorical")
      mode <- if (is.ordered(x))
        "ordered"
      else
        "unordered"

    ## write type tag
    writeTypeSDML(x, type, mode, file = file)

    ## write properties tag
    ## remove factor/POSIXxx attributes first
    xtmp <- x
    if (is.factor(x)) {
      attr(xtmp, "levels") <- NULL
      newclass <- class(xtmp)[!class(xtmp) %in% c("ordered","factor")]
      class(xtmp) <- if (!length(newclass)) NULL else newclass
    }
    if (inherits(x, "POSIXt")) {
      newclass <- class(xtmp)[!class(xtmp) %in% c("POSIXt","POSIXct","POSIXlt")]
      class(xtmp) <- if(!length(newclass)) NULL else newclass
    }

    writePropertiesSDML(attributes(xtmp), file = file, textdata = textdata,
                        sep = sep, na.string = na.string, null.string = null.string,
                        posinf.string = posinf.string, neginf.string = neginf.string,
                        nan.string = nan.string, true = true, false = false)

    if (is.null(textdata))
      textdata <- type != "character" && type != "datetime"

    x <- switch (type,
                 logical = if (textdata) ifelse(x, true, false) else x,
                 categorical = as.numeric(x),
                 numeric = x,
                 datetime = format(x, format="%Y-%m-%dT%H:%M:%S"),
                 character = x
                 )

    if (textdata) {
      ## textdata tag
      catSDML("<textdata sep=\"", sep, "\"",
              if (any(is.na(x[!is.nan(x)])))
                paste(" na.string=\"", markup(na.string), "\"", sep = ""),
              if (type == "character")
                paste(" null.string=\"", markup(null.string), "\"", sep = ""),
              if (type == "numeric") paste(
                    if (any(is.nan(x)))
                      paste(" nan.string=\"", markup(nan.string), "\"", sep = ""),
                    if (any(x == Inf, na.rm = TRUE))
                      paste(" posinf.string=\"", markup(posinf.string), "\"", sep = ""),
                    if (any(x == -Inf, na.rm = TRUE))
                      paste(" neginf.string=\"", markup(neginf.string), "\"", sep = ""),
                    sep = ""),
              if (type == "logical") paste(
                    if (any(x == true, na.rm = TRUE))
                      paste(" true=\"", markup(true), "\"", sep=""),
                    if (any(x == false, na.rm = TRUE))
                      paste(" false=\"", markup(false), "\"", sep=""),
                    sep = ""),
              ">",
              file = file
              )
      if (length(x)) {
        catSDML("\n", file = file)

        ## replacements
        x <- as.character(x)

        x[is.na(x)] <- na.string
        x[x == ""] <- null.string

        if (type == "numeric") {
          x <- gsub("NaN", nan.string, x)
          x <- gsub("^Inf", posinf.string, x)
          x <- gsub("-Inf", neginf.string, x)
        }

        if (type == "logical") {
          x[x == "TRUE"] <- true
          x[x == "FALSE"] <- false
        }

        ## write data
        cat(markup(x), sep=c(rep(substr(sep, 1, 1), 9), "\n"), file = file, append = TRUE)
      }

      catSDML("</textdata>\n", file = file)
    } else {
      ## data tag
      catSDML("<data>", file = file)

      if (length(x)) {
        catSDML("\n", file = file)

        if(!is.null(a <- attr(xtmp, "info")))
          attr(x, "info") <- a

        ## write data
        if (type == "numeric" && !is.null(mode) && mode == "complex")
          cetagsSDML(x, file = file)
        else
          etagsSDML(x, file = file)

        catSDML("\n", file = file)
      }
      catSDML("</data>\n", file = file)
    }
    catSDML("</array>\n", file=file)
  }
}

writeDimensionSDML <- function(x, file = "")
{
  ## make dataframes behave themselves
  x <- unclass(x)

  catSDML("<dimension>", file = file)
  if (length(x)) catSDML("\n", file = file)

  if (!is.null(dim(x))) {
    for (i in 1:length(dim(x))) {
      catSDML("<dim size=\"", dim(x)[i], "\"",
              if (!is.null(names(dimnames(x))[i]))
              paste(" name=\"", markup(names(dimnames(x))[i]), "\"", sep = ""),
              ">", file = file)

      if (!is.null(dimnames(x)[[i]]))
        etagsSDML(dimnames(x)[[i]], file = file)

      catSDML("</dim>\n", file = file)
    }
  } else if (length(x)) {
    catSDML("<dim size=\"", length(x), "\">", file = file)

    if (!is.null(names(x)))
      etagsSDML(names(x), file = file)

    catSDML("</dim>\n", file = file)
  }

  catSDML("</dimension>\n", file = file)
}

writePropertiesSDML <- function(attrib, file, textdata, ...)
{
  attrib[["dim"]] <- NULL
  attrib[["names"]] <- NULL
  attrib[["dimnames"]] <- NULL
  attrib[["length"]] <- NULL
  if (!is.null(textdata) && !textdata) attrib[["info"]] <- NULL

  if (!is.null(attrib) && length(attrib) > 0) {
    catSDML("<properties>\n", file = file)

    writeListArraySDML(attrib, file = file, textdata = textdata, ...)

    catSDML("</properties>\n", file = file)
  }
}

writeTypeSDML <- function(x, type, mode, file)
{
  catSDML("<type>", file = file)

  if (type %in% c("logical", "character", "datetime"))
    catSDML(" <", type, "/> ", file = file)
  else if (type == "numeric") {
    catSDML("<numeric>", file = file)
    if (mode == "complex") {
      catSDML("\n<complex/>\n", file = file)
    } else {
      catSDML("<",mode,">\n", file = file)
      tags(min(x, na.rm = TRUE), "min", file = file)
      tags(max(x, na.rm = TRUE), "max", file = file)
      catSDML("\n</",mode,">", file = file)
    }

    catSDML("</numeric>", file = file)
  } else {## categorical
    catSDML("\n<categorical mode=\"", mode, "\">\n", file = file)
    for (i in 1:length(levels(x)))
      catSDML("<label code=\"", i, "\">", markup(levels(x)[i]),
              "</label>\n", file = file)
    catSDML("</categorical>\n", file = file)
  }

  catSDML("</type>\n", file = file)
}


