HHscriptnames <- function(chapternumbers=NULL, edition=2) {
  if (edition != 1 && edition != 2) stop("edition must be either 1 or 2", call. = FALSE)
  if (is.null(chapternumbers))
    return(system.file(paste("scripts/hh", edition, "/", sep=""), package="HH"))
  HH.chapternames <- list(
    "1"=c(
      "1"   ="Ch01-intr",
      "2"   ="Ch02-data",
      "3"   ="Ch03-conc",
      "4"   ="Ch04-grap",
      "5"   ="Ch05-iinf",
      "6"   ="Ch06-oway",
      "7"   ="Ch07-mcomp",
      "8"   ="Ch08-rega",
      "9"   ="Ch09-regb",
      "10"  ="Ch10-regbb",
      "11"  ="Ch11-regc",
      "12"  ="Ch12-tway",
      "13"  ="Ch13-dsgn",
      "14"  ="Ch14-dsgntwo",
      "15"  ="Ch15-twtb",
      "16"  ="Ch16-npar",
      "17"  ="Ch17-logi",
      "18"  ="Ch18-tser",
      "intr"   ="Ch01-intr",
      "data"   ="Ch02-data",
      "conc"   ="Ch03-conc",
      "grap"   ="Ch04-grap",
      "iinf"   ="Ch05-iinf",
      "oway"   ="Ch06-oway",
      "mcomp"  ="Ch07-mcomp",
      "rega"   ="Ch08-rega",
      "regb"   ="Ch09-regb",
      "regbb"  ="Ch10-regbb",
      "regc"   ="Ch11-regc",
      "tway"   ="Ch12-tway",
      "dsgn"   ="Ch13-dsgn",
      "dsgntwo"="Ch14-dsgntwo",
      "twtb"   ="Ch15-twtb",
      "npar"   ="Ch16-npar",
      "logi"   ="Ch17-logi",
      "tser"   ="Ch18-tser"),
    "2"=c(
      "1"   ="intr",
      "2"   ="data",
      "3"   ="conc",
      "4"   ="grap",
      "5"   ="iinf",
      "6"   ="oway",
      "7"   ="mcomp",
      "8"   ="rega",
      "9"   ="regb",
      "10"  ="regbb",
      "11"  ="regc",
      "12"  ="tway",
      "13"  ="dsgn",
      "14"  ="dsgntwo",
      "15"  ="twtb",
      "16"  ="npar",
      "17"  ="logi",
      "18"  ="tser",
     ## "19"  ="likert",
     ## "20"  ="medphss",
      A = "RApx",
      B = "HHApx",
      C = "RcmdrApx",
      D = "RExcelApx",
      E = "ShinyApx",
      F = "Rpack",
      G = "PrcnApx",
      H = "otherApx",
      I = "mthp",
      J = "dstr",
      K = "edit",
      L = "typg",
      M = "emcs",
      N = "latexApx",
      O = "MSword")
    ##  = "grapb"
     )



  chapternumbers.char <- as.character(chapternumbers)
  names(chapternumbers.char) <- chapternumbers.char
  chapternames <- matrix(chapternumbers.char,
                      length(chapternumbers), 1,
                      dimnames=list(chapternumbers.char, "Absolute Pathname"))
  validname <- (chapternumbers.char == chapternumbers.char) ## initialize to all TRUE
  old.warn <- options(warn=-1)  ## suppress warning from as.numeric("abcd")
  for (chapternumber in chapternumbers) {
    if ((!is.na(as.numeric(chapternumber)) &&
         (as.numeric(chapternumber) >= 1) &&
         (as.numeric(chapternumber) <= length(HH.chapternames[[edition]]))) ||
        chapternumber %in% names(HH.chapternames[[edition]]))
      HHname <- HH.chapternames[[edition]][chapternumber]
    else {
      if (chapternumber %in% HH.chapternames[[edition]])
        HHname <- chapternumber
      else {
        HHname <- "Warning: Not a valid chapter"
        validname[as.character(chapternumber)] <- FALSE
      }
    }
    chapternames[as.character(chapternumber), 1] <- HHname
  }
  options(old.warn)
  pathnames <- chapternames
  for (chapseq in seq(along=chapternames)[validname]) {
    pathnames[chapseq] <- system.file(paste("scripts/hh", edition, "/",
                                              pathnames[chapseq],
                                              c(".r", ".R")[as.numeric(edition)], sep=""),
                                        package="HH")
  }
  pathnames[pathnames == ""] <- "Warning: No script file for chapter"
  pathnames
}


WindowsPath <- function(x, display=TRUE) {
  result <- gsub("/", "\\\\", x)
  if (!display) return(result)
  cat(paste(rownames(result), result, "\n"))
  invisible(result)
}

if (FALSE) {
HHscriptnames()
HHscriptnames(c(2, "RExcelApx"))
HHscriptnames(c("iinf","RApx"))
HHscriptnames(c(1:5))
HHscriptnames(c(1:32))
WindowsPath(HHscriptnames(c(1:32)))
WindowsPath(HHscriptnames(c(6:8)))
HHscriptnames(2)
HHscriptnames("2")
HHscriptnames("abcd")
HHscriptnames("mcomp")
HHscriptnames(42)
}

if (FALSE) { ## new
system.file("scripts/hh2", package="HH") ## Macintosh
## [1] "/Library/Frameworks/R.framework/Versions/3.1/Resources/library/HH/scripts/hh2"

system.file("scripts/hh2", package="HH") ## Windows
## [1] "C:/Program Files/R/R-devel/library/HH/scripts/hh2"
WindowsPath(system.file("scripts/hh2", package="HH"))
## C:\Program Files\R\R-devel\library\HH\scripts\hh2
WindowsPath(system.file("scripts/hh2", package="HH"), display=FALSE)
## [1] "C:\\Program Files\\R\\R-devel\\library\\HH\\scripts\\hh2"
tmp <- WindowsPath(system.file("scripts/hh2", package="HH"))
## C:\Program Files\R\R-devel\library\HH\scripts\hh2
tmp
## [1] "C:\\Program Files\\R\\R-devel\\library\\HH\\scripts\\hh2"
}
