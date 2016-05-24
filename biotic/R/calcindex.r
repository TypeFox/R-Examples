#' Calculate biotic indices for invertebrate samples
#' @description Calculates a range of freshwater invertebrate biotic indices
#'  in use in the UK (based on family level identification). Currently
#'  calculates BMWP (N-taxa and ASPT), Whalley revised BMWP, Whalley habitat
#'   specific BWMP (riffle, pool and riffle/pool), LIFE, PSI, WHPT
#'   (presence-absence and abundance-weighted) and AWIC.  For details of these
#'    indices see the vignette.
#'
#' @param df A dataframe containing list of taxon names and their abundances
#' in samples, along with sample identifiers.  Default format is for taxon
#' names to be in the first column and sample abundances in subsequent
#' columns with identifers as column headers. See built-in \code{\link{almond}}
#'  dataset for an example. If data are in the transposed format i.e taxa as
#'  columns and samples as rows, the \code{\link{transposedata}} function can
#'  be used prior to calculation.
#'
#' @param index A choice of index to be calculated. Defaults to \code{"BMWP"}.
#' Options are: \code{"BMWP"}, \code{"Whalley"}, \code{"Riffle"}, \code{"Pool"},
#' \code{"RiffPool"}, \code{"LIFE"}, \code{"PSI"}, \code{"WHPT"},
#' \code{"WHPT_AB"} and \code{"AWIC"}.
#' @param type Indicates type of data being processed. Options are "num" for
#' numeric data, "log" for integer log abundance categories (1-5) or "alpha"
#' for alphabetic abundance categories (A-E). Default value is "num".
#'
#' @return A data frame consisting of columns of index values with samples in
#' rows. The number of columns returned depends on the index selected.
#' @export calcindex
#' @examples
#' # use the built-in River Almond dataset, numeric abundances
#'
#' # calculate the BMWP index for this dataset
#' # 'index' and 'type' do not have to specified as defaults are used
#' # ("BMWP" and "num")
#'
#' calcindex(almond)
#'
#' # calculate the PSI index for this dataset
#' # type does not have to specified as default is used ("num")
#'
#' calcindex(almond, "PSI")
#'
#' # calculate the WHPT abundance-weighted index for this dataset
#'
#' calcindex(almond, "WHPT_AB")
#'
#' # example of processing data in alphabetic log abundance categories
#' # using the 'type' argument
#'
#' # 'braidburn' dataset contains alphabetic log category data
#' # see ?braidburn for details
#'
#' # calculate the Whalley revised BMWP index (including N-taxa and ASPT)
#'
#' calcindex(braidburn, "Whalley", "alpha")
#'
#' # example of processing data in numeric log abundance categories
#' # using the 'type' argument
#'
#' # 'greenburn' dataset contains numeric log category data
#' # see ?greenburn for details
#'
#' # calculate the LIFE index for this dataset
#'
#' calcindex(greenburn, "LIFE", "log")
#'

calcindex<-function(df, index="BMWP", type="num"){

  # check that a correct method has been specified
  INDICES<-c("BMWP", "Whalley", "Riffle", "Pool", "RiffPool", "PSI", "LIFE", "WHPT", "WHPT_AB", "AWIC")
  indextype<-pmatch(index, INDICES)
  if (is.na(indextype))
    stop("Invalid index choice")

  # check that a correct type has been specified
  TYPES<-c("num", "log", "alpha")
  datatype<-pmatch(type, TYPES)
  if (is.na(datatype))
    stop("Invalid data type specified")


  # if abundances are recorded as alphabetic categories, convert them
  if (type=="alpha"){
    df<-convertalpha(df)
  }
  # if abundances are recorded as log categories, convert them
  if (type=="log"){
    # check maximum value is 5
    maxabund<-max(df[,2:ncol(df)], na.rm=TRUE)
    if (maxabund>5){
        stop("Maximum value is > 5; not log categories")
    }
    df<-convertlog(df)
  }

  # check for and combine oligochaete families to class level

  # set up vector of oligochaete taxa
  families<-c("Lumbricidae", "Lumbriculidae", "Enchytraeidae", "Haplotaxidae", "Naididae", "Tubificidae", "Oligochaeta")

  # create logical vector of rows with oligochaetes
  present<- df[,1] %in% families

  # extract non oligochaete rows
  rest<-df[!present,]

  # subset rows with worms present
  worms<-df[present,]

  if (nrow(worms)!=0){
    # convert taxon to character for replacement
    worms[,1]<-as.character(worms[,1])

    # if there is more than one row of worms
    if(nrow(worms)>1){
      # if there is more than one sample
      if (ncol(worms)>2){
        # sum abundance across all oligochaetes and add to first row
        worms[1,-1]<-colSums(worms[,-1],na.rm=TRUE)
      } else {
        worms[1,-1]<-sum(worms[,2], na.rm=FALSE)
      }
    }
    # add taxon string back in
    worms[1,1]<-"Oligochaeta"

    # convert back to factor
    worms[,1]<-as.factor(worms[,1])

    # just take first row (sum)
    worms<-worms[1,]

    # recombine with rest of taxa
    df<-rbind(rest, worms)
  }

  # separate out sample taxon list
  taxonlist<-df[,1]

  # and samples
  samples<-df[,2:ncol(df)]

  # check whether the data look like presence-absence and warn if calculating abundance-based index
  maxabund<-max(samples,na.rm=TRUE)
  if (maxabund==1){
    if (index=="PSI" || index=="LIFE" || index=="WHPT_AB"){
      warning("Maximum abundance in samples is 1. Abundance-weighted indices may not be meaningful")
    }
  }

  # if only one sample present, need to process differently
  if (ncol(df)==2){

    # calculate scores
    output<-calcscore(samples, taxonlist=taxonlist, index=index)

    # transpose output
    output<-t(output)

    # add on sample name to row
    row.names(output)<-names(df[2])
  }
  # if there is more than one sample, apply across columns
  else{
    output<-apply(samples, 2, calcscore, taxonlist=taxonlist, index=index)
  }

  # only need to bind rows for multiple samples
  if (ncol(df)>2){
    output<-rbind(output)
    output<-t(output)
  }

  # add column names depending on index
  if (index=="BMWP"){
    colnames(output)<-c("BMWP", "Ntaxa", "ASPT")
  }
  if (index=="Whalley"){
    colnames(output)<-c("Whalley_BMWP", "Whalley_Ntaxa", "Whalley_ASPT")
  }
  if (index=="Riffle"){
    colnames(output)<-c("Riffle_BMWP", "Riffle_Ntaxa", "Riffle_ASPT")
  }
  if (index=="Pool"){
    colnames(output)<-c("Pool_BMWP", "Pool_Ntaxa", "Pool_ASPT")
  }
  if (index=="RiffPool"){
    colnames(output)<-c("RiffPool_BMWP", "RiffPool_Ntaxa", "RiffPool_ASPT")
  }
  if (index=="LIFE"){
    colnames(output)<-c("LIFE")
  }
  if (index=="PSI"){
    colnames(output)<-c("PSI")
  }
  if (index=="AWIC"){
    colnames(output)<-c("AWIC")
  }
  if (index=="WHPT"){
    colnames(output)<-c("WHPT_ASPT", "WHPT_Ntaxa")
  }
  if (index=="WHPT_AB"){
    colnames(output)<-c("WHPT_AB_ASPT", "WHPT_AB_Ntaxa")
  }

  # add on sample identifier column and set row names to null
  output<-as.data.frame(cbind.data.frame(row.names(output), output))
  row.names(output)<-NULL
  colnames(output)[1]<-"Sample"

  return(output)
}
