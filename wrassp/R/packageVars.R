##' package variable to force the usage of the logger
##' set to FALSE by default
##' @author Raphael Winkelmann
##' @export
useWrasspLogger <- FALSE



##' list of default output extensions,
##' track names and output type 
##' for each signal processing function in wrassp
##' @author Raphael Winkelmann
##' @export
wrasspOutputInfos = list("acfana" = list("ext"= c("acf"), "tracks"=c("acf"), "outputType"="SSFF"),
  "afdiff" = list("ext"= c("dwav"), "tracks"=c(""), "outputType"="wav"),
  "affilter" = list("ext"= c("hpf", "lpf", "bpf", "bsf"), "tracks"=c(""), "outputType"="wav"),
  "cepstrum" = list("ext"= c("cep"), "tracks"=c("cep"), "outputType"="SSFF"),
  "cssSpectrum" = list("ext"= c("css"), "tracks"=c("css"), "outputType"="SSFF"),
  "dftSpectrum" = list("ext"= c("dft"), "tracks"=c("dft"), "outputType"="SSFF"),
  "ksvF0" = list("ext"= c("f0"), "tracks"=c("F0"), "outputType"="SSFF"),
  "mhsF0" = list("ext"= c("pit"), "tracks"=c("pitch"), "outputType"="SSFF"),
  "forest" = list("ext"= c("fms"), "tracks"=c("fm", "bw"), "outputType"="SSFF"),
  "lpsSpectrum" = list("ext"= c("lps"), "tracks"=c("lps"), "outputType"="SSFF"),
  "rfcana" = list("ext"= c("rfc", "arf", "lar", "lpc"), "tracks"=c("rms", "gain", "arf|lar|lpc|rfc"), "outputType"="SSFF"),
  "rmsana" = list("ext"= c("rms"), "tracks"=c("rms"), "outputType"="SSFF"),
  "zcrana" = list("ext"= c("zcr"), "tracks"=c("zcr"), "outputType"="SSFF")
  )


##' list of possibly useful file formats for AsspDataObj corresponding to the
##' first element of the fileInfo attribute
##' @author Lasse Bombien
##' @docType data
##' @seealso \code{\link{AsspFileFormat}}
##' @format 
##' \tabular{rll}{
##' Code Name \tab code number \tab description\cr
##' RAW     \tab  1\tab headerless or unsupported format \cr
##' ASP_A   \tab  2\tab ASP with ASCII data \cr
##' ASP_B   \tab  3\tab ASP with binary data \cr
##' XASSP   \tab  4\tab xassp ASCII \cr
##' IPDS_M  \tab  5\tab labels in IPdS `MIX' format \cr
##' IPDS_S  \tab  6\tab labels in IPdS `SAMPA' format \cr
##' AIFF    \tab  7\tab Apple Audio Interchange File Format \cr
##' AIFC    \tab  8\tab AIFF extended for compressed data \cr
##' CSL     \tab  9\tab Kay Elemetrics Computerized Speech Lab \cr
##' CSRE    \tab 10\tab Computerized Speech Research Environment \cr
##' ESPS    \tab 11\tab Entropic Signal Processing System \cr
##' ILS     \tab 12\tab \cr
##' KTH     \tab 13\tab Kungliga Tekniska Hoegskolan Stockholm \cr
##' SWELL   \tab 13\tab commercial version of KTH \cr
##' SNACK   \tab 13\tab as Tcl extension \cr
##' SFS     \tab 14\tab University College London Speech Filing System \cr
##' SND     \tab 15\tab NeXT version of `SND' format \cr
##' AU      \tab 15\tab Sun version of `SND' format \cr
##' NIST    \tab 16\tab National Institute of Standards and Technology \cr
##' SPHERE  \tab 16\tab SPeech HEader REsources \cr
##' PRAAT_S \tab 17\tab UvA praat 'short' text file \cr
##' PRAAT_L \tab 18\tab UvA praat 'long' text file \cr
##' PRAAT_B \tab 19\tab UvA praat binary file \cr
##' SSFF    \tab 20\tab Simple Signal File Format \cr
##' WAVE    \tab 21\tab IBM/Microsoft RIFF-WAVE \cr
##' WAVE_X  \tab 22\tab RIFF-WAVE extended format (Revision 3) \cr
##' XLABEL  \tab 24\tab ESPS xlabel \cr
##' YORK    \tab 25\tab University of York (Klatt'80 parameters) \cr
##' UWM     \tab 26\tab University of Wisconsin at Madison (microbeam data) )\cr
##' }
##' @export
AsspFileFormats <- c(  
  RAW     =  1, ## headerless or unsupported format 
  ASP_A   =  2, ## ASP with ASCII data 
  ASP_B   =  3, ## ASP with binary data 
  XASSP   =  4, ## xassp ASCII 
  IPDS_M  =  5, ## labels in IPdS `MIX' format 
  IPDS_S  =  6, ## labels in IPdS `SAMPA' format 
  AIFF    =  7, ## Apple Audio Interchange File Format 
  AIFC    =  8, ##   AIFF extended for compressed data 
  CSL     =  9, ## Kay Elemetrics Computerized Speech Lab 
  CSRE    = 10, ## Computerized Speech Research Environment 
  ESPS    = 11, ## Entropic Signal Processing System 
  ILS     = 12, ## 
  KTH     = 13, ## Kungliga Tekniska Hoegskolan Stockholm 
  SWELL   = 13, ##   commercial version of KTH 
  SNACK   = 13, ##   as Tcl extension 
  SFS     = 14, ## University College London Speech Filing System 
  SND     = 15, ## NeXT version of `SND' format 
  AU      = 15, ## Sun version of `SND' format 
  NIST    = 16, ## National Institute of Standards and Technology 
  SPHERE  = 16, ##   SPeech HEader REsources 
  PRAAT_S = 17, ## UvA praat 'short' text file 
  PRAAT_L = 18, ## UvA praat 'long' text file 
  PRAAT_B = 19, ## UvA praat binary file 
  SSFF    = 20, ## Simple Signal File Format 
  WAVE    = 21, ## IBM/Microsoft RIFF-WAVE 
  WAVE_X  = 22, ##   RIFF-WAVE extended format (Revision 3) 
  XLABEL  = 24, ## ESPS xlabel 
  YORK    = 25, ## University of York (Klatt'80 parameters) 
  UWM     = 26  ## University of Wisconsin at Madison (microbeam data) )
  )


