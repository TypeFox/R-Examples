prettyPlotMain <- function(innervalues, format="expression") { ##
  values <- prettyUserValues(innervalues) ## lose their names
  charnames <- names(innervalues)
  listnoms <- lapply(charnames, formatName, format=format) ## returns a list of strings
  if (length(listnoms)>1) {
    if (format %in% c("ASCII", "charstring")) {
      noms <- paste(listnoms, collapse=", ")
      noms <- paste("(", noms, ")")
    } else { ## for expression
      mergenoms <- paste("(list(", paste(listnoms, collapse=", "), "))", sep="") ## string
    }
    values <- paste(values, collapse=", ~~")
    values <- paste("(list(", values, "))")
  } else {mergenoms <- listnoms[[1]]} ## 'values' is then given by first line
  if (("1D" %in% blackbox.getOption("DemographicModel")) && length(charnames)==1 && charnames=="latt2Ns2") {
    fullstring <- paste(" ", blackbox.getOption("GeoUnit"), sep="")
  } else fullstring <- ""
  if (format=="expression") fullstring <- gsub(pattern = " ", replacement = "~", x = fullstring)
  if(format %in% c("ASCII", "charstring")) {
    return(paste(mergenoms, "=", values, fullstring, sep=""))
  } else { ## expression
    total <- paste(mergenoms, "==", values, fullstring, sep="")
    return(parse(text=total))
  }
}
