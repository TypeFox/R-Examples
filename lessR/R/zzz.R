if (getRversion() >= "2.15.1") 
  globalVariables(c("mydata", "mycor"))


.onAttach <-
function(...) {

  packageStartupMessage("\n",
      "lessR 3.4.8      feedback: gerbing@pdx.edu      web: lessRstats.com\n",
      "-------------------------------------------------------------------\n",
      "1. Read a text, Excel, SPSS, SAS or R data file from your computer\n",
      "   into the mydata data table:  mydata <- Read()  \n",
      "2. For a list of help topics and functions, enter:  Help()\n",
      "3. Use theme function for global settings, Ex: theme(colors=\"gray\")\n")

  options(colors="dodgerblue")
  options(trans.fill.bar=0.25)
  options(trans.fill.pt=0.50)
  options(color.fill.bar = "#1874CCBF")  # .maketrans of dodgerblue3"  trans=.25
  options(color.fill.pt = "#1874CC80")  # .maketrans of "dodgerblue3"  trans=.50
  options(color.stroke.bar="steelblue4")
  options(color.stroke.pt="steelblue4")
  options(color.bg="#F2F4F5")  # rgb(242, 244, 245)
  options(color.grid="snow3")
  options(color.box="black")
  options(color.ghost=FALSE)
  options(color.heat="dodgerblue4")

  options(col.fill.bar = "#1874CCBF")  # .maketrans of dodgerblue3"  trans=.25
  options(col.fill.pt = "#1874CC80")  # .maketrans of "dodgerblue3"  trans=.50
  options(col.stroke.bar="steelblue4")
  options(col.stroke.pt="steelblue4")
  options(col.bg="#F2F4F5")  # rgb(242, 244, 245)
  options(col.grid="snow3")
  options(col.box="black")
  options(col.ghost=FALSE)
  options(col.heat="dodgerblue4")

  options(n.cat=8)
  options(lab.size=0.82)
  options(suggest=TRUE)
  options(quiet=FALSE)
  options(brief=FALSE)
  options(ghost=FALSE)

  options(explain=TRUE)
  options(interpret=TRUE)
  options(results=TRUE)
  options(document=TRUE)
  options(code=TRUE)

  options(show.signif.stars=FALSE)
  options(scipen=30)

}


.max.dd <- function(x) {

 n.dec <-function(xn) {
    xc <- format(xn)  # as.character(51.45-48.98) does not work
    nchar(xc)
    ipos <- 0
    for (i in 1:nchar(xc)) if (substr(xc,i,i)==".") ipos <- i
    n.dec <- ifelse (ipos > 0, nchar(xc)-ipos, 0)
    return(n.dec)
  }

  max.dd <- 0
  for (i in 1:length(x))
    if (!is.na(x[i])) if (n.dec(x[i]) > max.dd) max.dd <- n.dec(x[i])

  return(max.dd)
}


.getdigits <- function(x, min.digits) {
  digits.d <- .max.dd(x) + 1
  if (digits.d < min.digits) digits.d <- min.digits
  return(digits.d)
}


.dash <- function(ndash, cc, newline=TRUE) {
  if (missing(cc)) cc <- "-" 
  for (i in 1:(ndash)) cat(cc)
  if (newline) cat("\n") 
}


.dash2 <- function(ndash, cc="-") {
  tx <- ""
  if (!is.null(cc)) for (i in 1:(ndash)) tx <- paste(tx, cc, sep="")
  return(tx)
}


.plotList <- function(plot.i, plot.title) {
  mxttl <- 0
  for (i in 1:plot.i)
    if (nchar(plot.title[i]) > mxttl) mxttl <- nchar(plot.title[i])
  mxttl <- mxttl + 8
  cat("\n")
  .dash(mxttl, newline=FALSE)
  for (i in 1:plot.i) {
    cat("\n", "Plot ", i,": ", plot.title[i], sep="")
  }
  cat("\n")
  .dash(mxttl, newline=FALSE)
  cat("\n\n")

}


.plotList2 <- function(plot.i, plot.title) {
  tx <- character(length = 0)

  mxttl <- 0
  for (i in 1:plot.i)
    if (nchar(plot.title[i]) > mxttl) mxttl <- nchar(plot.title[i])
  mxttl <- mxttl + 8

  tx[length(tx)+1] <- .dash2(mxttl)
  for (i in 1:plot.i)
    tx[length(tx)+1] <- paste("Plot ", i,": ", plot.title[i], sep="")
  tx[length(tx)+1] <- .dash2(mxttl)

  return(tx)
}


.fmt <- function(k, d=getOption("digits.d"), w=0) {
  format(sprintf("%.*f", d, k), width=w, justify="right", scientific=FALSE)
}


.fmti <- function(k, w=0) {
  format(sprintf("%i", k), width=w, justify="right")
}


.fmtc <- function(k, w=0, j="right") {
  format(sprintf("%s", k), width=w, justify=j)
}


.fmtNS <- function(k) {
  format(k, scientific=FALSE )
}


.is.integer <- function(x, tol = .Machine$double.eps^0.5) {

  if (is.numeric(x)) {
    x <- na.omit(x)
    int.flg <- ifelse (abs(x - round(x)) < tol, TRUE, FALSE)  # each i of vector
    result.flg <- ifelse (all(int.flg), TRUE, FALSE)
  }
  else
    result.flg <- FALSE

  return(result.flg)
}


.xstatus <- function(var.name, dname, quiet=FALSE) {

  # see if analysis from data is based on a formula
  is.frml <- ifelse (grepl("~", var.name), TRUE, FALSE)

  # see if analysis is from descriptive stats or from data 
  from.data <- ifelse (var.name == "NULL", FALSE, TRUE)

  # see if the variable exists in the Global Environment
  in.global <- FALSE
  if (nchar(var.name)>0) if (exists(var.name, where=.GlobalEnv)) {
    if (!is.function(var.name)) { # a global "var" could be a function call 
      in.global <- TRUE
      #if (!quiet)
        #cat(">>> Note: ", var.name, "exists in the workspace, outside of",
            #"a data frame (table)\n")
    }
  }

  #see if "variable" is really an expression
  if (grepl("(", var.name, fixed=TRUE) ||  
      grepl("[", var.name, fixed=TRUE) ||  
      grepl("$", var.name, fixed=TRUE))  {
    txtA <- paste("A referenced variable in a lessR function can only be\n",
            "a variable name.\n\n", sep="")
    txtB <- "For example, this does not work:\n  > Histogram(rnorm(50))\n\n"
    txtC <- "Instead do this:\n  > Y <- rnorm(50)\n  > Histogram(Y)"
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        txtA, txtB, txtC, "\n")
  }

  if (!in.global && from.data) .nodf(dname)
 
  return(list(ifr=is.frml, fd=from.data, ig=in.global))
}


.getdfs <- function() {  # get list of data frames in global env

  objs <- function(x) class(get(x))

  dfs <- ls(.GlobalEnv)[sapply(ls(.GlobalEnv), objs) == 'data.frame']

  return(dfs)
}


.nodf <- function(dname) {

  # see if the data frame exists (mydata default), if x from data, not in Global Env
  if (!exists(dname, where=.GlobalEnv)) {
    dfs <- .getdfs()  # list of data frames in global env
    txtA <- ifelse(dname == "mydata", ", the default data table name, ", " ") 

    if ("Mydata" %in% dfs)
      txtM <- paste("Because you have a data table called Mydata,\n",
        " perhaps you meant to call it mydata, so just re-read \n",
        " into mydata instead of Mydata")
    else
      txtM <- paste(
        "If a data table is not named the default mydata, then to\n",
        "  analyze the variables in that data table, in the function call\n",
        "  for the analysis specify the actual data table name with\n",
        "  the data option\n",
        "For the data table called ", dfs[1], "\n",
        "  add the following to your function call:  , data=", dfs[1], "\n\n",
        "Or, just re-read the data into the mydata data table\n\n", sep="")

    if (length(dfs) == 0) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The analysis is of data values for one or more variables found\n",
        "  in a rectangular data table, with the data values for a \n",
        "  variable located in a column\n\n",
        "You have not yet read data into a data table for analysis,\n", 
        "  so the data table ", dname, txtA, "is not\n",
        "  available for analysis\n\n",
        "Create the data table using the Read function, usually\n",
        "  reading the data into an R data table called mydata\n\n",
        "To read a data file on your computer system into the mydata data\n", 
        "  table, in which you browse your file folders to locate the\n",
        "  desired date file, enter:\n",
        "     mydata <- Read()\n\n",
        "To read a data table from the web, enter:\n",
        "     mydata <- Read(\"web address\") \n",
        "In the web address include the http:// at the beginning\n",
        "  and also include the quotes around the web address\n\n")
    }
    else if (length(dfs) == 1) {
      nm <- parse(text=paste("names(", dfs[1],")"))
      nm <- eval(nm)
      for (i in 1:length(nm)) nm[i] <- paste(nm[i], " ")
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Data table ", dname, txtA, "does not exist\n\n",
        "You have read data into one data table, ", dfs[1], ", but that\n", 
        "  is not the data table ", dname, " that was to be analyzed\n\n", 
        "Following are the names of the variables that are available\n",
        "  for analysis in the ", dfs[1], " data table\n\n",
        "  ", nm, "\n\n",
        txtM, "\n\n")
    }
    else if (length(dfs) > 1) {
      dts <- ""
      for (i in 1:length(dfs)) dts <- paste(dts, dfs[i])
      if (dname == "mydata") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Data table ", dname, txtA, "does not exist\n\n",
          "Data tables you read and/or created: ", dts, "\n\n",
          "Perhaps you have a data table that contains the variables\n",
          "  of interest to be analyzed, but it is not named ", dname, "\n",
          "Can specify the actual name with the data option\n",
          "For example, for a data table named ", dfs[1], "\n",
          "  add the following to your function call:  , data=", dfs[1], "\n\n",
          "Or, just re-read the data into the mydata data table\n\n")
        }
      else {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Data table ", dname, txtA, "does not exist\n\n",
          "Perhaps you have a data table that contains the variables\n",
          "  of interest to be analyzed, but it is not named ", dname, "\n\n",
          "Data tables you have already read and/or created: ", dts, "\n\n",
          "Use an available data table, or read data into a new table\n\n")
      }
    }
  }
}


.xcheck <- function(var.name, dname, data) {

  if ( (!grepl(":", var.name, fixed=TRUE) &&
        !grepl("c(", var.name, fixed=TRUE)) ) { # x not var list

    if (grepl("\"", var.name)) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "No quotes with the variable name.\n\n")
    }

  # see if "variable" is really an expression
  if (grepl("(", var.name, fixed=TRUE) ||  grepl("[", var.name, fixed=TRUE))  {
    txtA <- paste("A referenced variable in a lessR function can only be\n",
            "a variable name\n\n", sep="")
    txtB <- paste("For example, for the Histogram function, this does not work:\n",
            "  > Histogram(rnorm(50))\n\n", sep="")
    txtC <- "Instead do this:\n  > Y <- rnorm(50)\n  > Histogram(Y)"
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        txtA, txtB, txtC, "\n")
  }

  # see if "variable" includes a $ 
  if ( grepl("$", var.name, fixed=TRUE))  {
    txtA <- paste("A referenced variable in a lessR function just includes\n",
            "the variable name\n\n", sep="")
    txtB <- paste("For example, for the Histogram function, this does not work:\n",
            "  > Histogram(mydata$Y)\n\n", sep="")
    txtC <- "Instead do this:\n  > Histogram(Y)"
    txtD <- "If you wish to specify a data table, use option: data"
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        txtA, txtB, txtC, "\n", txtD, "\n\n")
  }

    # see if variable exists in the data frame
    if (!exists(var.name, where=data)) { 
      dfs <- .getdfs()  # data frames in global

      txt1 <- ", the default name \n\n"
      txt2 <- "Either make sure to use the correct variable name, or\n"
      txt3 <- "specify the data table that contains the variable with option:  data\n"
      txt <- ifelse (dname == "mydata",  paste(txt1, txt2, txt3, sep=""), "\n")

      nm <- parse(text=paste("names(", dname,")"))
      nm <- eval(nm)
      for (i in 1:length(nm)) nm[i] <- paste(nm[i], " ")

      if (dname == "mydata")
         txtDef <- ", which is the default data table name\n"
      else
         txtDef <- ""

      if (length(dfs) == 1) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "You are attempting to analyze the variable ", var.name, " in the\n", 
          "  data table called ", dname, txtDef, "\n",
          "Unfortunately, variable ", var.name, " does not exist in ", dname, "\n\n",
          "The following variables are currently in the ", dname, " data table,\n",
          "  available for analysis:\n\n",
          "  ", nm,  "\n\n")
      }

      else if (length(dfs) > 1) {
        nm2 <- parse(text=paste("names(", dfs[1],")"))
        nm2 <- eval(nm2)
        for (i in 1:length(nm2)) nm2[i] <- paste(nm2[i], " ")
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "You are attempting to analyze the variable ", var.name, " in the\n", 
          "  data table called ", dname, txtDef, "\n",
          "Unfortunately, variable ", var.name, " does not exist in ", dname, "\n\n",
          "The following variables are currently in the ", dname, " data table,\n",
          "  available for analysis:\n\n",
          "  ", nm,  "\n\n",
          "You do have another data table, but it is named ", dfs[1], "\n",
          "The following variables are currently in the ", dfs[1], " data table,\n",
          "  available for analysis:\n\n",
          "  ", nm2,  "\n\n",
          "If a data table is not named the default mydata, then to\n",
          "  analyze the variables in that data table, in the function call\n",
          "  for the analysis specify the actual data table name with\n",
          "  the data option\n",
          "For the data table called ", dfs[1], "\n",
          "  add the following to your function call:  , data=", dfs[1],
          "\n\n", sep="") 
      }
    }
  }
}


# see if cor matrix exists as stand-alone or embedded in list structure
.cor.exists <- function(cor.nm) {

  if (!grepl("$cors", cor.nm, fixed=TRUE))  # no $cors in name
    is.there <- exists(cor.nm, where=.GlobalEnv)

  else {
    nm <- sub("$cors", "", cor.nm, fixed=TRUE)  # remove $cors from name
    if (!exists(nm, where=.GlobalEnv))  # root list exists?
      is.there <- FALSE
    else
      is.there  <- exists("cors", where=eval(parse(text=nm)))  #  cors inside?
  }
  if (!is.there) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "No correlation matrix entered.\n\n",
      "No object called ", cor.nm, " exists.\n\n",
      "Either enter the correct name, or calculate with: Correlation\n",
      "Or read the correlation matrix with: corRead\n\n", sep="")
  }

}


# get variable labels if they exist
.getlabels <- function(xlab=NULL, ylab=NULL, main=NULL, sub=NULL,
                       cex.lab=NULL, ...) {

  cutoff <- 0.76

  x.name <- getOption("xname")
  y.name <- getOption("yname")

  dname <- getOption("dname")  # not set for dependent option on tt
  if (!is.null(dname)) {
    if (exists(dname, where=.GlobalEnv))
      mylabels <- attr(get(dname, pos=.GlobalEnv), which="variable.labels")
    else
      mylabels <- NULL
  }
  else
    mylabels <- NULL

  if (!is.null(mylabels)) {
    x.lbl <- mylabels[which(names(mylabels) == x.name)]
    if (length(x.lbl) == 0) x.lbl <- NULL
    y.lbl <- mylabels[which(names(mylabels) == y.name)]
    if (length(y.lbl) == 0) y.lbl <- NULL
  }
  else {
    x.lbl <- NULL
    y.lbl <- NULL
  }

  # x-axis and legend labels
  if (is.null(x.lbl) && is.null(xlab))
    x.lab <- x.name

  else {  # process label
    if (!is.null(xlab))
      x.lab <- xlab  # xlab specified
    else if (!is.null(x.lbl)) 
      x.lab <- x.lbl

    if (length(x.lab) == 1  &&  !is.null(cex.lab)) {   # power.ttest: len > 1
      if (strwidth(x.lab, units="figure", cex=cex.lab) > cutoff) {
        brk <- nchar(x.lab)
        while (strwidth(substr(x.lab,1,brk), units="figure", cex=cex.lab) > cutoff)
          brk <- brk-1 
        while (substr(x.lab,brk,brk) != " ") brk <- brk-1
        x.lab <- paste(substr(x.lab,1,brk), "\n",
                       substr(x.lab,brk+1,nchar(x.lab)))
        while (strwidth(x.lab, units="figure", cex=cex.lab) > cutoff)
          cex.lab <- cex.lab-0.02
      }
    }

    var.nm <- ifelse(is.null(xlab), TRUE, FALSE)  # add var name to label?
    if (is.null(xlab))
      var.nm <- ifelse(is.null(x.name), FALSE, TRUE)
    if (grepl("Count", x.lab, fixed=TRUE)) var.nm <- FALSE 
    if (grepl("Proportion", x.lab, fixed=TRUE)) var.nm <- FALSE 
    if (grepl("Alternative", x.lab, fixed=TRUE)) var.nm <- FALSE 
    if (var.nm) {
      if (!grepl("\n", x.lab, fixed=TRUE))  # bquote removes the \n
        x.lab <- bquote(paste(italic(.(x.name)), ": ", .(x.lab))) 
      else
        x.lab <- paste(x.name, ":", x.lab)
    }
  }

  # y-axis and legend labels
  if (is.null(y.lbl) && is.null(ylab))
    y.lab <- y.name

  else {  # process label
    if (!is.null(ylab))
      y.lab <- ylab  # ylab specified
    else if (!is.null(y.lbl)) 
      y.lab <- y.lbl

    if (length(y.lab) == 1  &&  !is.null(cex.lab)) {   # power.ttest: len > 1
      if (strwidth(y.lab, units="figure", cex=cex.lab) > cutoff) {
        brk <- nchar(y.lab)
        while (strwidth(substr(y.lab,1,brk), units="figure", cex=cex.lab) > cutoff)
          brk <- brk-1 
        while (substr(y.lab,brk,brk) != " ") brk <- brk-1
        y.lab <- paste(substr(y.lab,1,brk), "\n",
                       substr(y.lab,brk+1,nchar(y.lab)))
        while (strwidth(y.lab, units="figure", cex=cex.lab) > cutoff)
          cex.lab <- cex.lab-0.02
      }
    }

    var.nm <- ifelse(is.null(ylab), TRUE, FALSE)  # add var name?
    if (is.null(ylab))
      var.nm <- ifelse(is.null(y.name), FALSE, TRUE)
    if (grepl("Frequency", y.lab, fixed=TRUE)) var.nm <- FALSE 
    if (grepl("Proportion", y.lab, fixed=TRUE)) var.nm <- FALSE 
    if (grepl("Power", y.lab, fixed=TRUE)) var.nm <- FALSE 
    if (var.nm) {
      if (!grepl("\n", y.lab, fixed=TRUE))  # bquote removes the \n
        y.lab <- bquote(paste(italic(.(y.name)), ": ", .(y.lab))) 
      else
        y.lab <- paste(y.name, ":", y.lab)
    }
  }


  if (!missing(main)) {
    if (!is.null(main)) main.lab <- main else main.lab <- ""
  }
  else
    main.lab <- NULL

  if (!missing(sub)) {
    if (!is.null(sub)) sub.lab <- sub else sub.lab <- NULL
  }
  else
    sub.lab <- NULL

  return(list(xn=x.name, xl=x.lbl, xb=x.lab, yn=y.name, yl=y.lbl, yb=y.lab,
              mb=main.lab, sb=sub.lab, cex.lab=cex.lab))
}


# get values for passed par parameters because ... cannot be passed to title 
#   in calling routine as each title will invoke sub if active,
#   and setting the line forces everything on the same line
.getdots <- function(...) {

  col.main <- NULL
  col.lab <- NULL
  sub.lab <- NULL
  col.sub <- NULL
  cex.main <- NULL

  dots <- list(...)  
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      if (names(dots)[i] == "col.main")  col.main <- dots[[i]] 
      if (names(dots)[i] == "col.lab")  col.lab <- dots[[i]] 
      if (names(dots)[i] == "sub.lab")  sub.lab <- dots[[i]] 
      if (names(dots)[i] == "col.sub")  col.sub <- dots[[i]] 
      if (names(dots)[i] == "cex.main")  cex.main <- dots[[i]] 
    }
  }

  return(list(col.main=col.main, col.lab=col.lab, sub.lab=sub.lab,
              col.sub=col.sub, cex.main=cex.main))

}


.varlist <- function(n.pred, i, var.name, pred.lbl, n.obs, n.keep, lvls=NULL) {

  if (i == 1)
    txt <- "Response Variable:  "
  else
      txt <- paste(pred.lbl, " Variable ", toString(i-1), ": ", sep="")
  cat(txt, var.name)

  dname <- getOption("dname")
  if (exists(dname, where=.GlobalEnv))
    mylabels <- attr(get(dname, pos=.GlobalEnv), which="variable.labels")
  else
    mylabels <- NULL

  if (!is.null(mylabels)) {
    lbl <- mylabels[which(names(mylabels) == var.name)]
    if (!is.null(lbl)) cat(", ", as.character(lbl))
  }

  if (!is.null(lvls)) if (i > 1) cat("\n  Levels:", lvls)
  cat("\n")

  if (i == n.pred+1) {
    cat("\n")
    cat("Number of cases (rows) of data: ", n.obs, "\n")
    cat("Number of cases retained for analysis: ", n.keep, "\n")
  }
}


.varlist2 <- function(n.pred, ind, var.name, pred.lbl, n.obs, n.keep, lvls=NULL) {
  tx <- character(length = 0)

  if (ind == 1)
    txt <- "Response Variable: "
  else {
    if (n.pred > 1)
      txt <- paste(pred.lbl, " Variable ", toString(ind-1), ": ", sep="")
    else
      txt <- paste(pred.lbl, " Variable: ", sep="")
  }
  if (pred.lbl == "Factor"  &&  ind > 1) tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(txt, var.name, sep="")

  dname <- getOption("dname")

  if (exists(dname, where=.GlobalEnv))
    mylabels <- attr(get(dname, pos=.GlobalEnv), which="variable.labels")
  else
    mylabels <- NULL
  if (exists(dname, where=.GlobalEnv))
    myunits <- attr(get(dname, pos=.GlobalEnv), which="variable.units")
  else
    myunits <- NULL

  if (!is.null(mylabels)) {
    lbl <- mylabels[which(names(mylabels) == var.name)]
    unt <- myunits[which(names(myunits) == var.name)] 
    if (!is.null(unt)) if (nzchar(unt))  if(!is.na(unt))
      lbl <- paste(lbl, " (", unt, ")", sep="")
    if (!is.null(lbl))
      tx[length(tx)] <- paste(tx[length(tx)], ", ", as.character(lbl), sep="")
  }

  if (!is.null(lvls)) {
    tx2 <- "  Levels:"
    for (i in 1:length(lvls)) tx2 <- paste(tx2, lvls[i])
    tx[length(tx)+1] <- tx2 
  }

  if (ind == n.pred+1) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("Number of cases (rows) of data: ", n.obs)
    tx[length(tx)+1] <- paste("Number of cases retained for analysis: ", n.keep)
  }

  return(tx)
}


.title <- function(x.name, y.name, x.lbl, y.lbl, isnullby) {

  txt1 <- x.name
  if (!is.null(x.lbl)) txt1 <- paste(txt1, ": ", x.lbl, sep="")

  if (isnullby) txt1 <- paste("---", txt1, "---")
  else {
    txt2 <- paste(y.name, sep="")
    if (!is.null(y.lbl)) txt2 <- paste(txt2, ": ", y.lbl, sep="") 
  }

  cat("\n")
  cat(txt1, "\n")
  if (!isnullby) {
    cat("  - by levels of - \n")
    cat(txt2, "\n")
    ndash <- max(nchar(txt1),nchar(txt2))
    .dash(ndash)
  }
  cat("\n")

}


.title2 <- function(x.name, y.name, x.lbl, y.lbl, isnullby, new.ln=TRUE) {

  txt1 <- x.name
  if (!is.null(x.lbl)) txt1 <- paste(txt1, ": ", x.lbl, sep="")

  if (isnullby) {
    txt1 <- paste("---", txt1, "---")
    if (new.ln) txt1 <- paste(txt1, "\n", sep="")
  }
  else {
    txt2 <- paste(y.name, sep="")
    if (!is.null(y.lbl)) txt2 <- paste(txt2, ": ", y.lbl, sep="") 
  }

  tx <- character(length = 0)

  tx[length(tx)+1] <- txt1
  if (!isnullby) {
    if (is.null(y.lbl))
      tx[length(tx)+1] <- "  - by levels of - "
    else
      tx[length(tx)+1] <- "\n  - by levels of - \n"
    tx[length(tx)+1] <- txt2
    if (is.null(y.lbl))
      tx[length(tx)+1] <- .dash2(max(nchar(txt1),nchar(txt2)))
  }

  return(tx)

}


.axes <- function(x.lvl, y.lvl, axT1, axT2, par1, par3,
         cex.axis, col.axis, rotate.values=0, offset=0.5, ...) {

  if (is.null(x.lvl)  &&  !is.null(axT1)) {  # numeric, uses axT1)
    axis(1, at=axT1, labels=FALSE, tck=-.01)
    dec.d <- .getdigits(round(axT1,6),1) - 1
    text(x=axT1, y=par3, labels=.fmt(axT1,dec.d),
         pos=1, xpd=TRUE, cex=cex.axis, col=col.axis, srt=rotate.values,
         offset=offset, ...)
  }
  else if (!is.null(x.lvl)) {  # categorical, uses x.lvl
    axis(1, at=axT1, labels=FALSE, tck=-.01)
    text(x=axT1, y=par3, labels=x.lvl,
         pos=1, xpd=TRUE, cex=cex.axis, col=col.axis, srt=rotate.values,
         offset=offset, ...)
  }

  if (is.null(y.lvl)  &&  !is.null(axT2)) {
    axis(2, at=axT2, labels=FALSE, tck=-.01)
    dec.d <- .getdigits(round(axT2,6),1) - 1
    text(x=par1, y=axT2, labels=.fmt(axT2,dec.d),
         pos=2, xpd=TRUE, cex=cex.axis, col=col.axis, ...)
  }
  else if (!is.null(y.lvl)) {
    axis(2, at=axT2, labels=FALSE, tck=-.01)
    text(x=par1, y=axT2, labels=y.lvl,
         pos=2, xpd=TRUE, cex=cex.axis, col=col.axis, ...)
  }
}


# axis labels
.axlabs <- function(x.lab, y.lab, main.lab, sub.lab, max.lbl.y,
                    xy.ticks=TRUE, offset=0.5, ...) {

  lbl.lns <- ifelse(xy.ticks, 3, 1)

  # xlab positioning
  lblx.lns <- ifelse (grepl("\n", x.lab, fixed=TRUE), lbl.lns + 0.4, lbl.lns)
  lblx.lns <- ifelse (!is.null(sub.lab), lblx.lns - 1.4, lblx.lns - 0.7)
  if (!xy.ticks) lblx.lns <- lblx.lns + .75
  if (offset > 0.5) lblx.lns <- lblx.lns + 0.5

  # ylab positioning
  multi <- FALSE
  for (i in 1:length(y.lab))
    if (!is.null(y.lab))
      if (grepl("\n", y.lab[i], fixed=TRUE)) multi <- TRUE  # multi-line
  lm <- par("mar")[2]  # get the current left margin
  lbly.lns <- ifelse(multi, lm - 2, lm - 1.2)

  title(xlab=x.lab, line=lblx.lns, ...)
  title(sub=sub.lab, line=lblx.lns+1, cex.sub=0.76, ...)
  title(ylab=y.lab, line=lbly.lns, ...)
  title(main=main.lab, ...)

}


.RSadj <- function(bubble.size=0.25, cex.axis) {

  # scale for regular R or RStudio
  if (options("device") == "RStudioGD") bubble.size <- bubble.size*1.5
  size.axis <- ifelse (options("device") != "RStudioGD", cex.axis, cex.axis*1.118)
  sz.lab <- getOption("lab.size")
  size.lab <- ifelse (options("device") != "RStudioGD", sz.lab, sz.lab*1.119)
  size.txt <- ifelse (options("device") != "RStudioGD", 0.7, 0.8)

  return(list(bubble.size=bubble.size, size.axis=size.axis, size.lab=size.lab,
              size.txt=size.txt))

}


.showfile <- function(fname, txt) {
  if (getwd() == "/")
    workdir <- "top level (root) of your file system"
  else
    workdir <- getwd()
  cat("The", txt, "written at the current working directory\n")
  cat("       ", fname, " in:  ", workdir, "\n")
  cat("\n")
}


.showfile2 <- function(fname, txt) {
  if (getwd() == "/")
    workdir <- "top level (root) of your file system"
  else
    workdir <- getwd()

  tx <- character(length = 0)

  tx[length(tx)+1] <- paste("The", txt, "written at the current working directory.")
  tx[length(tx)+1] <- paste("       ", fname, " in:  ", workdir)

  return(tx)

}


.pdfname <- function(analysis, x.name, go.pdf, pdf.nm, pdf.file) {
  if (go.pdf) {
    if (pdf.nm)
      if (!grepl(".pdf", pdf.file))
        pdf.fnm <- paste(pdf.file, ".pdf", sep="")
      else
        pdf.fnm <- pdf.file
    else
      pdf.fnm <- paste(analysis, "_", x.name, ".pdf", sep="")
  }
  else
    pdf.fnm <- NULL

  return(pdf.fnm)
}


# see if manage graphics or just sequentially plot
.graphman <- function() {

  in.RStudio <- ifelse (options("device") != "RStudioGD", FALSE, TRUE)

  in.knitr <- ifelse (is.null(options()$knitr.in.progress), FALSE, TRUE)

  manage.gr <- ifelse (!in.RStudio  &&  !in.knitr, TRUE, FALSE)

  return(manage.gr)
}


# manages the graphics system (not in RStudio or knitr)
.graphwin <- function(wnew=1, d.w=NULL, d.h=NULL) {

  dl <- dev.list()
  dl2 <- dl[which(dl==2)]  # device #2
  dl.more <- dl[which(dl>2)]  # devices larger than #2

  # remove all open windows past device 2
  if (length(dl.more) > 0) {
    min.dd <- dl.more[which(dl.more==min(dl.more))]
    max.dd <- dl.more[which(dl.more==max(dl.more))]
    for (i in min.dd:max.dd) dev.off(which=i)
  }

  off.two <- ifelse (length(dl2) == 0, TRUE, FALSE)

  # open graphics windows
  # if not already present, generate a null window for #2 and then remove
  if (off.two) wnew <- wnew + 1
  RS <- ifelse (options("device") == "RStudioGD", TRUE, FALSE)
  for (i in 1:wnew) 
    if (is.null(d.w) && is.null(d.h))
      dev.new(noRStudioGD=RS)
    else
      dev.new(width=d.w, height=d.h, noRStudioGD=RS)
  if (off.two) dev.off(which=2)

}


.opendev <- function(pdf.fnm, pdf.width, pdf.height) {

  if (is.null(pdf.fnm)) {
    if (options("device") != "RStudioGD"  &&  is.null(options()$knitr.in.progress)) {
      .graphwin(1)
      orig.params <- par(no.readonly=TRUE)
      on.exit(par(orig.params))
    }
  }
  else 
    pdf(file=pdf.fnm, width=pdf.width, height=pdf.height)

}


.is.num.cat <- function(x, n.cat) {

  x <- sort(unique(na.omit(x)))

  nu.x <- length(x)

  # num.cat var is integer with small number of unique values
  if (.is.integer(x)  &&  nu.x <= n.cat) {
    eq.int <- TRUE
    d.x <- diff(x)  # check for equal intervals
    if (nu.x > 2) {
      for (i in 2:(length(x)-1)) {
        if ((abs(d.x[i-1] - d.x[i]) > 0.0000000001)) eq.int <- FALSE
      }
      status <- eq.int  # num.cat var has equal intervals
    }
    else
      status <- TRUE

  }
  else
    status <- FALSE 

  return(status)

}


.ncat <- function(analysis, x.name, nu, n.cat, brief=FALSE) {

  cat("\n")
  cat(">>> ", x.name, " has only only ", nu, " equally spaced unique ",
      "integer values <= n.cat=", n.cat, "\n",
      "    so treat as categorical, and perhaps convert to an R factor\n", sep="")

  if (!brief)
    cat("    For numeric, set n.cat smaller than ", nu, 
        " with ", analysis, " or globally with  theme", sep="")

  cat("\n")

}


.corcolors <- function(R, NItems, main, bm=3, rm=3, diag=NULL,
                       pdf.file, pdf.width, pdf.height) {

    if (!is.null(diag)) {
      for (i in 1:NItems) R[i,i] <- diag
      cat("\nNote: To provide more color separation for off-diagonal\n",
          "      elements, the diagonal elements of the matrix for\n",
          "      computing the heat map are set to 0.\n", sep="")
    }

    max.color <- getOption("col.heat")
    hmcols <- colorRampPalette(c("white",max.color))(256)
    
    .opendev(pdf.file, pdf.width, pdf.height)  # set up graphics

    heatmap(R[1:NItems,1:NItems], Rowv=NA, Colv="Rowv", symm=TRUE,
      col=hmcols, margins=c(bm,rm), main=main)

    if (!is.null(pdf.file)) {  # terminate pdf graphics
      dev.off()
      .showfile(pdf.file, "plot")
    }
}


.maketrans <- function(col.name, trans.level) {
  r.tr <- col2rgb(col.name)[1]
  g.tr <- col2rgb(col.name)[2]
  b.tr <- col2rgb(col.name)[3]

  #trans.level <- (1-trans.level) * 256
  col.trans <- rgb(r.tr, g.tr, b.tr, alpha=trans.level, maxColorValue=256)

  return(col.trans)
}


# discrete color steps with no order
.col.discrete <- function(bright=FALSE) {

  # based on rainbow_hcl(24,c=40,l=70) from colorspace
  if (!bright)
    clr <- c("#92ADD6", "#D59E93", "#7FB88B", "#BDAA78", "#D09AC8",
             "#64B8C1", "#C19FD3", "#D899B8")
  # based on rainbow_hcl(8,c=80,l=65) from colorspace
  else
    clr <- c( "#57A2F0", "#E68369", "#1DB556", "#BD9B00", "#E475D6",
              "#00B7C7", "#CA80EA", "#F072B8")

  return(clr)

}


# for BarChart 1-var and PieChart for an ordered factor
.ordcolors <- function(colors, col.low, col.hi) {

    if (colors == "dodgerblue") { 
      if (is.null(col.low)) col.low <- rgb(.765,.824,.886)
      if (is.null(col.hi)) col.hi <- "dodgerblue4"
    }
    else if (colors == "sienna") { 
      if (is.null(col.low)) col.low <- "#F7E9E2"
      if (is.null(col.hi)) col.hi <- "sienna3"
    }
    else if (colors == "blue") { 
      if (is.null(col.low)) col.low <- "slategray2"
      if (is.null(col.hi)) col.hi <- "slategray4"
    }
    else if (colors == "gray") {
      if (is.null(col.low)) col.low <- "gray90"
      if (is.null(col.hi)) col.hi <- "gray25"
    }
    else if (colors == "gray.black") {
      if (is.null(col.low)) col.low <- "gray70"
      if (is.null(col.hi)) col.hi <- "gray25"
    }
    else if (colors == "green") {
      if (is.null(col.low)) col.low <- "darkseagreen1"
      if (is.null(col.hi)) col.hi <- "darkseagreen4"
    }
    else if (colors == "rose") {
      if (is.null(col.low)) col.low <- "mistyrose1"
      if (is.null(col.hi)) col.hi <- "mistyrose4"
    }
    else if (colors == "gold") {
      if (is.null(col.low)) col.low <- "goldenrod1"
      if (is.null(col.hi)) col.hi <- "goldenrod4"
    }
    else if (colors == "red") { 
      if (is.null(col.low)) col.low <- "coral1"
      if (is.null(col.hi)) col.hi <- "coral4"
    }
    else if (colors == "orange.black") { 
      if (is.null(col.low)) col.low <- rgb(255,173,91, maxColorValue=256)
      if (is.null(col.hi)) col.hi <- rgb(169,66,2, maxColorValue=256)
    }
    else if (colors == "purple") { 
      if (is.null(col.low)) col.low <- "purple1"
      if (is.null(col.hi)) col.hi <- "purple4"
    }
    else if (colors == "white") { 
      if (is.null(col.low)) col.low <- "gray30"
      if (is.null(col.hi)) col.hi <- "gray90"
    }

    return(list(col.low=col.low, col.hi=col.hi))
}


.to256 <- function(trans.level)
   trn <- (1-getOption(trans.level))*256


# change class call to class character
.fun.call.deparse <- function(fun.call) {

  fc.d <- deparse(fun.call)
  if (length(fc.d) > 1) {  # multiple lines
    fc <- fc.d[1]
    for (i in 2:length(fc.d)) fc <- paste(fc, fc.d[i], sep="")  
  }
  else
    fc <- fc.d

  fc <- sub("     ", " ", fc, fixed=TRUE)
  fc <- sub("    ", " ", fc, fixed=TRUE)
  fc <- sub("  ", " ", fc, fixed=TRUE)

  return(fc)

}


# get the value for a specified function argument
.get.arg <- function(argm, fc) {

  loc <- regexec(argm, fc)
  strt1 <- loc[[1]]  # beginning of argument
  if (strt1 > 0) {
    j <- strt1
    while(substr(fc, start=j, stop=j) != "\"") j <- j + 1 
    strt <- j
    j <- j + 1  # first " after ,
    while(substr(fc, start=j, stop=j) != "\"") j <- j + 1 
    stp <- j  # second " after ,
    value <- substr(fc, start=strt, stop=stp)
  }
  else
    value <- ""

  return(value)
}


# remove argument and character value from a function call
.rm.arg <-  function(argm, fc) {

  loc <- regexec(argm, fc)[[1]]  # beginning of argument

  if (loc > 0) {

    first.arg <- ifelse (substr(fc, loc-1, loc-1) == "(", TRUE, FALSE)

    j <- loc
    if (!first.arg)  # is not first argument, start at preceding comma
      while(substr(fc, start=j, stop=j) != ",") if (j > 0) j <- j - 1 
    strt <- j  #  closing parentheses or comma before argument

    while(substr(fc, start=j, stop=j) != "\"") if (j < 1000) j <- j + 1 
    j <- j + 1  # first " after ,
    while(substr(fc, start=j, stop=j) != "\"") if (j < 1000) j <- j + 1 
    stp <- j  # second " after ,

    if (first.arg) stp <- stp + 2  # remove trailing comma and space

    remv <- substr(fc, start=strt, stop=stp)
    fc.new <- sub(remv, "", fc, fixed=TRUE)

  }

  return(fc.new)
}


# remove argument and logical value from a function call
.rm.arg.l <-  function(argm, fc) {

  loc <- regexec(argm, fc)[[1]]  # beginning of argument

  if (loc > 0) {

    first.arg <- ifelse (substr(fc, loc-1, loc-1) == "(", TRUE, FALSE)

    j <- loc
    if (!first.arg)  # is not first argument, start at preceding comma
      while(substr(fc, start=j, stop=j) != "," &&  substr(fc, start=j, stop=j) != "")
         if (j < 1000) j <- j + 1 
    stp <- j  #  closing parentheses or comma before argument
    if (first.arg) stp <- stp + 2  # remove trailing comma and space
    strt <- loc - 1

    remv <- substr(fc, start=strt, stop=stp)
    fc.new <- sub(remv, "", fc, fixed=TRUE)
    fc.new <- sub(",,", "", fc.new, fixed=TRUE)  # hack

  }

  return(fc.new)
}


# remove argument and Non-String value from a function call
.rm.arg.ns <-  function(argm, fc) {

  loc <- regexec(argm, fc)[[1]]  # beginning of argument

  if (loc > 0) {

    first.arg <- ifelse (substr(fc, loc-1, loc-1) == "(", TRUE, FALSE)

    j <- loc
    if (!first.arg)  # is not first argument, start at preceding comma
      while(substr(fc, start=j, stop=j) != ",") if (j > 0) j <- j - 1 
    strt <- j  #  closing parentheses or comma before argument

    dlm <- c(",", ")")

    j <- j + 1
    while(!(substr(fc, start=j, stop=j) %in% dlm))
      if (j < 1000) j <- j + 1 

    stp <- j  # got a "," or a ")" 
    stp <- stp - 1  # retain the "," or ")"

    if (first.arg) stp <- stp + 2  # remove trailing comma and space

    remv <- substr(fc, start=strt, stop=stp)
    fc.new <- sub(remv, "", fc, fixed=TRUE)

  return(fc.new)
  }

}


.outliers <- function(x) {

  tx <- character(length = 0)

  outliers <- boxplot.stats(x)$out

  if (length(outliers>0) && length(unique(na.omit(x)>3))) {
    tx[length(tx)+1] <- paste("(Box plot) Outliers:", length(outliers))

    lo.whisker <- boxplot.stats(x)$stats[1]
    lo.out <- outliers[outliers < lo.whisker]
    lo.out <- sort(lo.out, decreasing=FALSE)
    lo.len <- length(lo.out)
    tx[length(tx)+1] <- "Small: "
    if (lo.len > 0) {
      if (lo.len <= 25) {
        for (i in 1:lo.len)
          tx[length(tx)] <- paste(tx[length(tx)], .fmtNS(lo.out[i]))
      }
      else {
        for (i in 1:16) 
          tx[length(tx)] <- paste(tx[length(tx)], .fmtNS(lo.out[i]))
        tx[length(tx)] <- paste(tx[length(tx)], "... ")
        for (i in (lo.len-5):lo.len)
          tx[length(tx)] <- paste(tx[length(tx)], .fmtNS(lo.out[i]))
      }
    }
    else
      tx[length(tx)] <- paste(tx[length(tx)], "none")

    hi.whisker <- boxplot.stats(x)$stats[5]
    hi.out <- outliers[outliers > hi.whisker]
    hi.out <- sort(hi.out, decreasing=FALSE)
    hi.len <- length(hi.out)
    tx[length(tx)+1] <- "Large:"
    if (hi.len > 0) {
      if (hi.len <= 25) {
        for (i in 1:hi.len) 
          tx[length(tx)] <- paste(tx[length(tx)], .fmtNS(hi.out[i]))
      }
      else {
        for (i in 1:16)
          tx[length(tx)] <- paste(tx[length(tx)], .fmtNS(hi.out[i]))
        tx[length(tx)] <- paste(tx[length(tx)], "... ")
        for (i in (hi.len-5):hi.len)
          tx[length(tx)] <- paste(tx[length(tx)], .fmtNS(hi.out[i]))
      }
    }
    else
      tx[length(tx)] <- paste(tx[length(tx)], "none")

  }

  else
    tx <- ""

  return(tx)
}


.prntbl <- function(x, digits.d=2, cut=0, cc="-", cors=FALSE,
                    brk=NULL, bnd=NULL, v1.nm=NULL, v2.nm=NULL) {

# brk: ... replaces rows not printed 
# bnd: boundary between groups

  tx <- character(length = 0)

  max.ch <- ifelse (cors, 3, 0)  # max char per column, 0 is not applicable

  # width of column 1
  max.c1 <- 0
  for (i in 1:nrow(x)) {
    c1 <- nchar(rownames(x)[i])
    if (c1 > max.c1) max.c1 <- c1
  }
  if (!is.null(v2.nm)) if (nchar(v2.nm) > max.c1) max.c1 <- nchar(v2.nm) 
  max.c1 <- max.c1 + 2

  # widths of variable names
  colnm.w <- integer(length=ncol(x))
  for (i in 1:ncol(x)) 
    colnm.w[i] <- nchar(colnames(x)[i])

  # width of columns
  max.ln <- integer(length=ncol(x))
  for (j in 1:ncol(x)) {
    if (is.numeric(x[,j])) {
      c.val <- 0        
      for (i in 1:nrow(x)) {
        i.val <- nchar(formatC(x[i,j], digits=digits.d, format="f"))
        if (i.val > c.val) c.val <- i.val
      }
    }
    else {
      c.val <- 0        
      for (i in 1:nrow(x)) {
        i.val <- nchar(as.character(x[i,j]))
        if (i.val > c.val) c.val <- i.val
      }
    }
      #c.val <- 4
    if (!cors)
      max.ln[j] <- max(colnm.w[j], c.val) + 1
    else {
      max.ln[j] <- max(colnm.w[j], 4)
      if (max.ch > 0) max.ln[j] <- max.ch
      if (max.ln[j] > 4) max.ln[j] <- max.ln[j] + 1
    }
    if (max.ln[j] < 4) max.ln[j] <- 4
  }

  if (!is.null(cc)) tx[length(tx)+1] <- .dash2(sum(max.ln)+max.c1, cc=cc)

  # matrix for potentially multi-row column names
  if (max.ch > 0) {
    nr.ind.lbl <- integer(length=ncol(x))
    for (i in 1:ncol(x))
      nr.ind.lbl[i] <- ((nchar(colnames(x)[i]) + (max.ch-1)) %/% max.ch)

    nr.lbl <- max(nr.ind.lbl)  # n row of labels
    col.nm <- matrix(nrow=nr.lbl, ncol=ncol(x))
    for (i in 1:nrow(col.nm)) {
      for (j in 1:ncol(col.nm)) { 
        srt <- ((i-1)*max.ch) + 1
        stp <- srt + (max.ch - 1) 
        col.nm[i,j] <- substr(colnames(x)[j], srt, stp) 
        #if (nchar(col.nm[i,j]) > 0)  # left adjust within column
          #while (nchar(col.nm[i,j]) <= (max.ch-1))
            #col.nm[i,j] <- paste(col.nm[i,j], " ", sep="")
      }
    }
  }
  else {
    nr.lbl <- 1
    col.nm <- matrix(nrow=1, ncol=ncol(x))
    for (j in 1:ncol(col.nm)) col.nm[1,j] <- colnames(x)[j] 
  }
  # for each row, shift down value if next row is "", repeat
  if (nr.lbl > 1) { 
    for (k in 2:nrow(col.nm)) {  # repeat for each row
      for (i in 2:nrow(col.nm)) {
        for (j in 1:ncol(col.nm)) { 
          if (nchar(col.nm[i,j]) == 0) {
            col.nm[i,j] <- col.nm[i-1,j]  
            col.nm[i-1,j] <- ""
          }
        }
      }
    }
  }

  blnk <- format("", width=max.c1-1)
  if (!is.null(v1.nm)) tx[length(tx)+1] <- paste(blnk, v1.nm)
  # write col labels
  for (i in 1:nr.lbl) {  # for each row of column labels
    if (is.null(v2.nm))
      tx[length(tx)+1] <- format("", width=max.c1)
    else
      tx[length(tx)+1] <- paste(" ", v2.nm, format("", width=max.c1-nchar(v2.nm)-2), sep="")
    for (j in 1:ncol(x)) {
      wd <- max.ln[j]
      tx[length(tx)] <- paste(tx[length(tx)], .fmtc(col.nm[i,j], w=wd), sep="")
      if (!is.null(bnd)) if (j %in% bnd)
        if (i == nr.lbl)
          tx[length(tx)] <- paste(tx[length(tx)], "|", sep="")
        else
          tx[length(tx)] <- paste(tx[length(tx)], " ", sep="")
    }
  }
  if (!is.null(bnd))
    tx[length(tx)+1] <- .dash2(sum(max.ln)+max.c1+length(bnd), cc="-")

  # factor vars to char vars
  if (is.data.frame(x)) { 
    i.col <- sapply(x, is.factor)
    x[i.col] <- lapply(x[i.col], as.character)
  }

  # write values
  for (i in 1:nrow(x)) {
    if (i %in% brk) tx[length(tx)+1] <- "..."
    rwnm <- paste(" ", rownames(x)[i])
    if (is.null(v2.nm))
      tx[length(tx)+1] <- format(rwnm, width=max.c1, justify="right")
    else
      tx[length(tx)+1] <- format(rwnm, width=max.c1-1, justify="right")

    for (j in 1:ncol(x)) {
      if (is.integer(x[i,j]))
        tx[length(tx)] <- paste(tx[length(tx)], .fmti(x[i,j], w=max.ln[j]), sep="")

      else if (is.numeric(x[i,j])) {
        wd <- max.ln[j] 
        if (cors) {
          if (max.ln[j] > 5) wd <- max(6, max.ln[j]+1) + 1 
          else wd <- max(6, max.ln[j]+1)
          cs <- .fmt(x[i,j], d=digits.d, w=wd)
          cs <- sub("0.", "", cs, fixed=TRUE)
          cs <- sub(" 1.00", "100", cs, fixed=TRUE)
        }
        else
          cs <- .fmt(x[i,j], d=digits.d, w=wd)
        if (abs(x[i,j]) < cut) cs <- paste(rep(" ", wd-2), collapse="")
        tx[length(tx)] <- paste(tx[length(tx)], cs, sep="")
        if (!is.null(bnd)) if (j %in% bnd)
          tx[length(tx)] <- paste(tx[length(tx)], "|", sep="") 
      }

      else if (is.character(x[i,j]))
        tx[length(tx)] <- paste(tx[length(tx)], .fmtc(x[i,j], w=max.ln[j]) , sep="") 
    }

    if (!is.null(bnd)) if (i %in% bnd)
      tx[length(tx)+1] <- .dash2(sum(max.ln)+max.c1+length(bnd), cc="-")
  }

  return(tx)

}  # end .prntbl
