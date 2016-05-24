#' Meta-Field Tag Extraction
#'
#' It extracts some new useful field tags from standard ISI/SCOPUS fied tag codify.
#' @param M is a data frame obtained by the converting function \code{\link{convert2df}}.
#'        It is a data matrix with cases corresponding to articles and variables to Field Tag in the original ISI or SCOPUS file.
#' @param Field is a character object. New tag exctracted from aggregated data is specified by this string.
#' Field can be equal to one of this tags:
#' \tabular{lll}{
#' \code{"CR_AU"}\tab   \tab First Author of each cited reference\cr
#' \code{"CR_SO"}\tab   \tab Source of each cited reference\cr
#' \code{"AU_CO"}\tab   \tab Affiliation Country of each co-author}
#'
#' @param sep is the field separator character. This character separates strings in each column of the data frame. The default is \code{sep = ";"}.
#' @return the bibliometric data frame with a new column containing data about new field tag indicated in the argument \code{Field}.
#'
#'
#'
#' @examples
#' # Example 1: First Authors for each cited reference
#'
#' data(scientometrics)
#' scientometrics <- metaTagExtraction(scientometrics, Field = "CR_AU", sep = ";")
#' unlist(strsplit(scientometrics$CR_AU[1], ";"))
#'
#'
#' #Example 2: Source for each cited reference
#'
#' data(scientometrics)
#' scientometrics <- metaTagExtraction(scientometrics, Field = "CR_SO", sep = ";")
#' unlist(strsplit(scientometrics$CR_SO[1], ";"))
#'
#' #Example 3: Affiliation country for co-author
#'
#' data(scientometrics)
#' scientometrics <- metaTagExtraction(scientometrics, Field = "AU_CO", sep = ";")
#' scientometrics$AU_CO[1:10]
#'
#' @seealso \code{\link{scopus2df}} for converting ISO or SCPUS Export file into a data frame.
#' @seealso \code{\link{biblioAnalysis}} function for bibliometric analysis
#'

metaTagExtraction<-function(M, Field = "CR_AU", sep = ";"){


size=dim(M)
FCAU=list(NULL)
CCR=NULL
M$CR=str_replace_all(as.character(M$CR),"DOI;","DOI ")
CR=M$CR


if (Field=="CR_AU"){

listCAU=strsplit(as.character(CR),sep)
listCAU=lapply(listCAU,function(l) l=l[nchar(l)>10])  ## delete not congruent references

  # vector of cited authors
  for (i in 1:size[1]){
    FCAU[[i]]=str_replace_all(trim.leading(sub(",.*", "", listCAU[[i]])), "[[:punct:]]", "")
    CCR[i]=paste(FCAU[[i]],collapse=";")}

  M$CR_AU=CCR
  }

if (Field=="CR_SO"){
  listCAU=strsplit(as.character(CR),sep)
  FCAU=list(NULL)

  # vector of cited Journals
  if (M$DB[1]=="ISI"){
  for (i in 1:size[1]){

    elem=strsplit(as.character(listCAU[[i]]),",")
    ind=lengths(elem)
    if (max(ind)>2) {
    elem=elem[ind>2]
    FCAU[[i]]=trim.leading(unlist(lapply(elem,function(l) l[[3]])))
    CCR[i]=paste(FCAU[[i]],collapse=";")}
    else {CCR[[i]]=NA}}

    } else if (M$DB[1]=="SCOPUS") {


    for (i in 1:size[1]){

      listCAU[[i]]=gsub(".*?\\) ", "", listCAU[[i]])
      elem=strsplit(as.character(listCAU[[i]]),",")
      ind=lengths(elem)
      CCR[[i]]=NA
      if (length(ind)>0){
      if (max(ind)>2) {
        elem=elem[ind>2]
        FCAU[[i]]=trim.leading(unlist(lapply(elem,function(l) l[[1]])))
        CCR[i]=paste(FCAU[[i]],collapse=";")}}
      } 
  }

  M$CR_SO=CCR
}

if (Field=="AU_CO"){
  # Countries
  data("countries",envir=environment())
  countries=as.character(countries[[1]])
  if (M$DB[1]=="ISI"){
  countries=as.character(sapply(countries,function(s) paste0(s,".",collapse="")))}
  else if (M$DB[1]=="SCOPUS"){
    countries=as.character(sapply(countries,function(s) paste0(s,";",collapse="")))}

  M$AU_CO=NA
  C1=M$C1
  C1[which(is.na(C1))]=M$RP[which(is.na(C1))]
  C1=gsub("\\[.*?\\] ", "", C1)
  C1=paste(C1,";",sep="")
  for (i in 1:size[1]){
    if (!is.na(C1[i])){
    ind=unlist(sapply(countries, function (l) (gregexpr ( l , C1[i],fixed=TRUE))))
    if (sum(ind>-1)>0) {M$AU_CO[i]=paste(names(ind[ind>-1]),collapse=";")}
    }
  }
  M$AU_CO=gsub("[[:digit:]]","",M$AU_CO)
  M$AU_CO=gsub(".", "", M$AU_CO, fixed = TRUE)
  M$AU_CO=gsub(";;", ";", M$AU_CO, fixed = TRUE)


}

return(M)
}
