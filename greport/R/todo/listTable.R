#' List Table
#'
#' Convert a data.frame object into a LaTeX table.
#'
#' @param fileName character. A description of the file connection.
#' @param longtable logical. Toggle \sQuote{longtable} or \sQuote{table} environment. Defaults to \sQuote{TRUE}.
#' @param landscape logical. Use \sQuote{landscape} environment. Defaults to \sQuote{FALSE}.
#' @param caption character. Main table caption.
#' @param fontSize character. Define font size from one of the following: \sQuote{tiny},
#' \sQuote{scriptsize}, \sQuote{footnotesize}, \sQuote{small}, \sQuote{normalsize}, 
#' \sQuote{large}, \sQuote{Large}, \sQuote{LARGE}, \sQuote{huge}, \sQuote{Huge}.  Default is \sQuote{small}.
#' @param dataframe data.frame. Provides content for list table.
#' @param zebraPattern character. Defaults to \sQuote{none}, other options are \sQuote{plain}, \sQuote{group} and \sQuote{plaingroup.}
#' \sQuote{plaingroup} is only recommended for large groups (more than four objects in a group).
#' @param by character. Column used to generate zebra pattern, defaulting to the first column in \sQuote{dataframe}.
#' @param orderGroups logical. If \sQuote{TRUE} order the data by \sQuote{by}.  This is recommened
#' when \sQuote{zebraPattern} is set to \sQuote{group} or \sQuote{plaingroup}.
#' @param colNames character vector. Define column name headers for \sQuote{dataframe}.
#' @param vars character vector. Column names to select from \sQuote{dataframe}.
#' @param fixedColVars character vector. Column variables.
#' @param fixedColWdths character vector. Fixed width for each column.
#' @param markVar character vector. Marker variable.
#' @param markVarVal character vector. Value for each marker variable.
#' @param toLatexChar logical. If \sQuote{TRUE} text will be checked and escaped should they contain special LaTeX characters.
#' @param appendix logical. If \sQuote{TRUE} the function will require the \code{subsection} and \code{marker}
#' arguments since it will have to reference the tables in the future.
#' @param subsection character. Name of document subsection that refers to this table.
#' @param marker character. Marker for document subsection that refers to this table.
#' @param append logical. If \sQuote{TRUE} output will be appended instead of overwritten.
#' @export
#' @examples
#' listTable(fileName='', caption="\\label{table:listtable}Table of groupings", zebraPattern='group', 
#' dataframe=data.frame(code=c('(a)', '(b)', '(c)'), def=c("apple, orange", "dog, cat, horse", "windows, linux, mac")),
#' appendix=FALSE) # print generated table to standard out

listTable <- function(fileName,
                       longtable=TRUE, landscape=FALSE,
                       caption = "", fontSize="small",
                       dataframe, zebraPattern="none", by=names(dataframe)[1], orderGroups=FALSE,
                       colNames = names(dataframe),
                       vars =names(dataframe), fixedColVars=c(), fixedColWdths=c(),
                       markVar="", markVarVal="",
                       toLatexChar=TRUE,
                       appendix=TRUE, subsection=NULL, marker=NULL, append=FALSE){
  ### appendix: if TRUE the function will require the subsection and marker arguments
  ### since it will have to reference the tables in the future.
  ### append: if TRUE all the out put will be appended to a file specified in fileName 

  #internal constants and functions definitions
  fontSizes <- c("tiny","scriptsize","footnotesize","small",
                 "normalsize","large","Large","LARGE","huge","Huge")
  zebraPatterns <- c("none", "plain", "group", "plaingroup")
  #NOTE: pattern "plaingroup" is recommended only for large groups (more than 4 objects in a group)

  #!!! not to remove: adjusted pallet white&gray
  pallet1 <- list(lightwhite=c(0,0,0,0), darkwhite=c(0,0,0,0.07),
                  lightgray =c(0,0,0,0.2), darkgray=c(0,0,0,0.27),
                  middlegray=c(0,0,0,0.18), red=c(0,0.9,0.3,0))
  #!!! not to remove: adjusted pallet: orange&blue1
  pallet2 <- list(lightwhite=c(0,0.05,0.15,0), darkwhite=c(0,0.1,0.3,0),
                  lightgray=c(0.3,0.15,0.075,0), darkgray=c(.4,0.2,0.1,0),
                  middlegray=c(0.3,0.15,0.075,0), red=c(0,0.7,1,0))
  #!!! not to remove: adjusted pallet: orange&blue2
  pallet3 <- list(lightwhite=c(0,0.04,0.12,0), darkwhite=c(0,0.08,0.24,0),
                  lightgray=c(0.22,0.11,0.055,0), darkgray=c(0.3,0.15,0.075,0),
                  middlegray=c(0.22,0.11,0.055,0), red=c(0,0.7,1,0))
  #!!! not to remove: adjusted pallet: purple&yellow
  pallet4 <- list(lightwhite=c(0.0025,0.05,0.2,0), darkwhite=c(0.05,0.1,0.4,0),
                  lightgray=c(.235,0.235,0.1125,0), darkgray=c(.3,0.3,0.15,0),
                  middlegray=c(.235,0.235,0.1125,0), red=c(0,0.9,0.3,0))
  pallet <- pallet1

  latexTextMode <- function(str) {
    #internal constants and functions definitions
    latexTextModeSpec <- function(char){
      if (char == "\\"){
        retStr <- "$\\backslash$"
      }else{
        if (char == "^" | char == "~"){
          retStr <- paste("\\verb*+",char,"+",sep="")
        }else{
          retStr <- paste("\\",char,sep="")
        }
      }
      retStr 
    }
    latexTextModeMath <- function(char){
      retStr <- paste("$",char,"$",sep="")
      retStr
    }
    charToLatexTextChar <- function(char){
      if (is.na(latexCharText[char])){
        char
      }else{
        latexCharText[char]
      }
    }

    latexTextModeText <- function(char){
      if (char %in% latexSpecialChar){
        retStr <- latexTextModeSpec(char)
      }else{
        if (char %in% latexMathChar){
          retStr <- latexTextModeMath(char)
        }else{
          retStr <- NULL
        }
      }
      retStr
    }

    latexSpecialChar <- c("#","$","%","^","_","{","}","~","&","\\")
    latexMathChar <- c("<", ">", "|")
    latexSpecAndMath <- c(latexSpecialChar,latexMathChar)
    latexCharText <- sapply(X = latexSpecAndMath, FUN=latexTextModeText)

    #beginning of the function latexTextMode
    spl <- unlist(strsplit(str,""))
    paste(sapply(X = spl, FUN=charToLatexTextChar),collapse="")    
  }

  processBeginCommand <- function(beginCom=c(), outFile,
                                  caption = "", fontSize="small", colNames = c(),
                                  dataframe, zebraPattern,
                                  fixedColVars=c(), fixedColWdths=c(),
                                  markVar="", markVarVal="") {
    #internal constants and functions definitions
    commandBegin <- function(command, outFile){
      if (command %in% fontSizes){
        cat("\n{\\",command,"\n",sep="", file=outFile)
      }else{
        cat("\n\\begin{",command,"}",sep="", file=outFile)
      }
    }
    commandEnd <- function(command, outFile){
      if (command %in% fontSizes){
        cat("}\n", file=outFile)
      }else{
        cat("\\end{",command,"}\n",sep="", file=outFile)
      }
    }
    latexCaption <- function(captionFill, longtable, outFile){
      cap <- paste("\n\\caption{",captionFill,"}",sep="")
      if (longtable){
        cat(cap,"\\\\\n", sep="", file=outFile)
      }else{
        cat(cap,"\n", sep="", file=outFile)
      }
    }
    processColsFormat <- function(data, fixedColVars=c(), fixedColWdths=c(), outFile){
      colFormat <- c()
      for (n in names(data)){
        if (n %in% fixedColVars){
          colFormat <- c(colFormat,paste("p{",fixedColWdths[fixedColVars==n],"pt}", sep=""))
        }else{
          colFormat <- c(colFormat,"l")
        }
      }
      cat(" {",paste(colFormat, collapse=""),"}", sep="", file=outFile)
    }
    hline <- function(outFile, number=1){
      cat(paste(rep("\\hline",number),collapse=""),"\n", file=outFile)
    }

    processColsHead <- function(colNames, longtable, caption="", outFile){
      processColNames <- function(colNames, style, outFile){
        hline(outFile,2)
        cat(style, colNames[1], file=outFile)
        for (i in 2:length(colNames)){
          cat("&", style, colNames[i], file=outFile)
        }
        cat("\\\\\n", file=outFile)
        hline(outFile)
      }
      headStyle <- "\\bfseries"
      otherStyle <- "\\bfseries\\em"
      colNum <- length(colNames)
      processColNames(colNames, headStyle, outFile)
      hline(outFile)
      if (longtable){
        cat("\\endfirsthead\n", file=outFile)
        cat("\\caption[]{",caption,"{",otherStyle," (continued)}} \\\\\n", file=outFile)
        processColNames(colNames, headStyle, outFile)
        cat("\\endhead\n", file=outFile)
        hline(outFile)
        cat("\\multicolumn{",colNum,"}{r}{",otherStyle," Continued on next page}\\\\\n",
            sep="", file=outFile)
        cat("\\endfoot\n", file=outFile)
        hline(outFile)
        cat("\\multicolumn{",colNum,"}{r}{",otherStyle," End}\\\\\n",
            sep="", file=outFile)
        cat("\\endlastfoot\n", file=outFile)
      }
    }

    processRows <- function(data, zebraPattern, outFile, markVar, markVarVal){
      #internal constants and functions definitions
      if(FALSE) {
        processRow <- function(row, color, outFile, markVarIndex){
          if (!is.na(color)){
            cat("\\rowcolor{",color,"}\n",sep="", file=outFile)
          }
          cat(row[[1]], file=outFile)
          for (i in c(2:length(row))){
            if (i==markVarIndex){
              cat(" &", "\\color{red}{\\bfseries\\em ", row[[i]],"}", file=outFile)
            }else{
              cat(" &", row[[i]], file=outFile)
            }
          }
          cat("\\\\\n", file=outFile)
        }
      }
      if(TRUE) {
        processRow <- function(row, color, outFile, markVarIndex){
          rowStr <- ""
          if (!is.na(color)){
            rowStr <- paste(rowStr,"\\rowcolor{",color,"}\n",sep="")
          }
          rowStr <- paste(rowStr,row[[1]],sep="")
          for (i in c(2:length(row))){
            if (i==markVarIndex){
              rowStr <- paste(rowStr," &", "\\color{red}{\\bfseries\\em ", row[[i]],"}",sep="")
            }else{
              rowStr <- paste(rowStr," &", row[[i]],sep="")
            }
          }
          rowStr <- paste(rowStr,"\\\\\n", sep="")
          cat(rowStr, file=outFile)
        }
      }
      #beginning of the function processRows
      markVarIndex <- match(markVar, names(data))
      if (markVar!="" & is.na(markVarIndex)) stop("Error: Variable to mark is not in dataframe names\n")
      if (length(data[[1]]) != 0){
        for (i in c(1:length(data[[1]]))){
          if (!is.na(markVarIndex)){
            if (data[i,markVarIndex]==markVarVal){
              processRow(data[i,], zebraPattern[i], outFile, markVarIndex)
            }else{
              processRow(data[i,], zebraPattern[i], outFile, -1)
            }
          }else{
            processRow(data[i,], zebraPattern[i], outFile, -1)
          }
        }
      }
      hline(outFile)
    }

    #beginning of the function processBeginCommand
    if (length(beginCom[!is.na(beginCom)])>0) {
      commandBegin(beginCom[1],outFile)
      if (beginCom[1] == "table"){
        latexCaption(caption, beginCom[1] == "longtable",outFile)
      }
      if (beginCom[1] == "tabular"){
        processColsFormat(data=dataframe, fixedColVars, fixedColWdths, outFile)
        processColsHead(colNames = colNames, beginCom[1] == "longtable", outFile=outFile)
        processRows(data=dataframe, zebraPattern=zebraPattern, outFile,
                    markVar, markVarVal)
      }
      if (beginCom[1] == "longtable"){
        processColsFormat(data=dataframe, fixedColVars, fixedColWdths, outFile)
        latexCaption(caption, beginCom[1] == "longtable",outFile)
        processColsHead(colNames = colNames, beginCom[1] == "longtable", caption, outFile)
        processRows(data=dataframe, zebraPattern=zebraPattern, outFile,
                    markVar, markVarVal)
      }
      processBeginCommand(beginCom[2:(length(beginCom)+1)], outFile,
                          caption, fontSize, colNames,
                          dataframe, zebraPattern,
                          fixedColVars, fixedColWdths,
                          markVar, markVarVal)
      commandEnd(beginCom[1], outFile)
    }
  }

  makePattern <- function(zebraPattern, by){
    #internal constants and functions definitions
    plain <- function(col1, col2, len){
      pattern <- rep(c(col1,col2), len)
      pattern <- pattern[1:len]  
      pattern
    }
    group <- function(by){
      col <- 1
      pattern <- rep(col, length(by))
      current <- by[1]
      for (i in 2:length(by)){
        nextg <- by[i]
        if (current==nextg) {
          pattern[i] <- col          
        }else{
          col <- col*(-1)
          pattern[i] <- col          
        }
        current <- nextg
      }
      pattern
    }

    #beginning of the function makePattern
    if (zebraPattern != "none") {
      if (!(zebraPattern %in% zebraPatterns)) stop("Error: Illigal Zebra Pattern\n")
      if (zebraPattern=="plain"){
        pattern <- plain("lightwhite","middlegray",length(by))
      }
      if (zebraPattern=="group"){
        pattern <- group(by)
        pattern <- ifelse(pattern>0,"middlegray","lightwhite")
      }
      if (zebraPattern=="plaingroup"){
        pattern <- group(by)
        light <- plain("darkwhite","lightwhite",length(by))
        dark <- plain("lightgray","darkgray",length(by))
        pattern <- ifelse(pattern>0,dark,light)
      }
      pattern
    } else {
      rep(NA, length(by))
    }
  }

  defineColors <- function(outFile){
    for (n in names(pallet)){
      cat("\\definecolor{",n,"}{cmyk}{",pallet[[n]][1],",",
                                        pallet[[n]][2],",",
                                        pallet[[n]][3],",",
                                        pallet[[n]][4],"}\n",sep="",file=outFile)
    }
  }

  #beginning of the function listTable
  if (appendix && is.null(marker)) {
    stop("Argument 'marker' has to be provided\n")
  }

  if(append) {
    openMode <- "at"
  } else {
    openMode <- "wt"
  }
  # allow file to be standard out
  if(fileName == "") {
    outFile <- ""
  } else {
    outFile <- file(fileName, open=openMode)
    on.exit(close(outFile))
  }
  if(appendix) {
    if (!is.null(subsection)){
      cat(paste("\\subsection{", subsection, "}\n", sep="") , file=outFile)
    }
    cat(paste("\\label{", marker, "}\n", sep=""), file=outFile)
  }

  dataframe <- dataframe[,vars]
  if (orderGroups){
    dataframe <- dataframe[order(dataframe[[by]]),]
  }
  if (zebraPattern %in% c("group","plaingroup")){
    ordered <- dataframe[[by]][order(dataframe[[by]])]
    if (!all(dataframe[[by]]==ordered) && !all(dataframe[[by]]==ordered[length(ordered):1])){
      cat("\nWARNING: It is recommended to order the data by", by, "\n",
          "         when argument 'zebraPattern' is set to 'group' or 'plaingroup'.\n",
          "         It can be done by setting argument 'orderGroups' to TRUE.\n")
    }
  }
  for (n in names(dataframe)){
    dataframe[[n]] <- as.character(as.character(dataframe[[n]]))
    if (toLatexChar)
      dataframe[[n]] <- sapply(X = dataframe[[n]], FUN=latexTextMode)
  }
  pattern <- makePattern(zebraPattern, dataframe[[by]])
  beginCommands <- c()
  if (landscape) beginCommands <- c(beginCommands,"landscape")
  beginCommands <- c(beginCommands, fontSize)
  if (!landscape) beginCommands <- c(beginCommands,"center")
  if (longtable) beginCommands <- c(beginCommands,"longtable")
  else beginCommands <- c(beginCommands,c("table","tabular"))
  defineColors(outFile)
  processBeginCommand(beginCommands, outFile,
                      caption, fontSize, colNames,
                      dataframe, pattern,
                      fixedColVars, fixedColWdths,
                      markVar, markVarVal)
  if (appendix){
    cat("\\clearpage\n", file=outFile)
  }
}
