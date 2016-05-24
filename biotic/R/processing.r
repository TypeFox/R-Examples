# Locates BMWP composite families in sample data
# @description Takes a listing of taxa and abundances and locates those
# taxa that are represented as composites in the BMWP index to avoid double
# counting of scores.
#
# @param composites A dataframe containing the list of composite taxa with
# each pair of composites in two columns.
# @param df A dataframe containing list of taxa in first column, followed by
#  columns of abundances with sample names in header row.
# @return An indexing vector, giving the rows of the sample data that need
# to be deleted to avoid double counting.

findcomposites<-function(composites, df){

  # check for match with each composite pair. Use column index to avoid
  # assumptions over title of column
  firstrow<-match(composites[1], df[,1])
  secondrow<-match(composites[2], df[,1])

  # is the second of the pair present?
  if (!is.na(secondrow)){

    # if so then if the first is also present, then need to delete the row
    if (!is.na(firstrow)){
      rowforremoval<-firstrow
    } else {
      rowforremoval<-NA
      }
  } else {
    rowforremoval<-NA
  }
  return(rowforremoval)
  # end of function
}


# Extract taxon scores based on presence.
# @description Extracts a subset of biotic index score table data based on
# taxon names in first row of dataframe passed in which represent taxa
# present in any of the samples being analysed.
#
# @param df A dataframe or atomic vector representing the sample. If in vector
#  form, should be a string vector of taxon names. If a dataframe, should
#  contain taxon names in strong format in first column, followed by columns
#  of numeric abundances in samples.
# @param indextable A dataframe of biotic index values, with strings of taxon
#  names in the first column and scores in subsequent columns. Usually
#  loaded from object in the sysdata.rda file.
# @return A dataframe consisting of the rows of the index table matching the
# taxa present in the sample.

extractrows<-function(df, indextable){

  # create logical vector of taxa present in sample
  present<-indextable$Taxon %in% df[,1]

  # subset of indices for present taxa
  scorerows<-indextable[present,]

  return(scorerows)
}


# Classify taxa into log class score sets based on their sample abundance.
# @description Takes a dataframe of taxon names and abundances and returns a dataframe of taxon index scores (through a call to extractrows), along with a 'class' index column indicating the log abundance class for each taxon in the sample.
#
# @param df A dataframe containing list of taxa in first column, followed by a column of abundances with sample names in header row.
# @return A data frame consisting of taxon names in the first column followed by all taxon index values and a log abundance class index column.

logslice<-function(df){
  # 1-9 abundance
  logclass1<-df[df[,2]>0 & df[,2]<10,]
  log1scores<-extractrows(logclass1, indextable)
  log1scores$class<-rep(1, times=nrow(log1scores))

  # 10-99 abundance
  logclass2<-df[df[,2]>9 & df[,2]<100,]
  log2scores<-extractrows(logclass2, indextable)
  log2scores$class<-rep(2, times=nrow(log2scores))

  # 100-999 abundance
  logclass3<-df[df[,2]>99 & df[,2]<1000,]
  log3scores<-extractrows(logclass3, indextable)
  log3scores$class<-rep(3, times=nrow(log3scores))

  # 1000+ abundance
  logclass4<-df[df[,2]>999,]
  log4scores<-extractrows(logclass4, indextable)
  log4scores$class<-rep(4, times=nrow(log4scores))
  sliced<-list(log1scores, log2scores, log3scores, log4scores)
  sliced<-do.call(rbind, sliced)

  return(sliced)
}

#' Transpose data layout
#' @description Transposes a dataset, correctly processing column and
#' row labels.
#' @param df A dataframe containing abundances of invertebrate taxa in
#' different samples.
#' @return A data frame transposing the input data, with row and column
#' labels processed correctly.
#' @export transposedata
#' @examples
#' # transpose the built-in River Almond dataset
#' # this would have to be transposed back to original format for calculation
#'
#' transposedata(almond)

transposedata<-function(df){

  rowlabs<-df[,1]

  collabs<-names(df)

  df_t<-as.data.frame(t(df[,-1]))

  names(df_t)<-rowlabs

  df_t<-cbind(as.factor(collabs[-1]), df_t)

  row.names(df_t)<-NULL

  names(df_t)[1]<-""

  return(df_t)
}

# Convert alphabetic log abundance classes.
# @description Takes data recorded as alphabetic log abundance classes
# and converts them to numeric format, with arbitary abundances within
# each class: A=1, B=10, C=100, D=1000, E=10000.
#
# @param df A dataframe containing abundances of invertebrate taxa in
# different samples recorded as alphabetic categories (character format).
# @return A data frame with alphabetic categories converted to
# appropriate numeric values for index calculation.

convertalpha<-function(df){
  # check the number of samples being processed
  numsamples<-ncol(df)-1
  # if just one, extract sample name
  if (numsamples==1){
    samplename<-names(df[2])
  }
  # separate first column (taxa/labels) from data
  firstcol<-df[,1]
  # rest of df is alphabetic abundance categories
  df<-df[,-1]
  # convert to character format
  if (numsamples>1){
    df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  } else {
    df<-data.frame(as.character(df), stringsAsFactors=FALSE)
    names(df)<-samplename
  }

  # convert categories to integers within classes
  df[df=="A"]<-1
  df[df=="B"]<-10
  df[df=="C"]<-100
  df[df=="D"]<-1000
  df[df=="E"]<-10000

  # convert to numeric values using data.matrix
  df<-as.data.frame(data.matrix(df))

    # add sample labels back on
  df<-cbind.data.frame(firstcol, df)

  # return converted values
  return(df)
}

# Convert numeric log abundance classes.
# @description Takes data recorded as integer log abundance classes
# and converts them to actual numeric values, with arbitary abundances within
# each class: 1=1, 2=10, 3=100, 4=1000, 5=10000.
#
# @param df A dataframe containing abundances of invertebrate taxa in
# different samples recorded as integer categories.
# @return A data frame with integer log categories converted to
# appropriate numeric values for index calculation.

convertlog<-function(df){
  # check the number of samples being processed
  numsamples<-ncol(df)-1
  # if just one, extract sample name
  if (numsamples==1){
    samplename<-names(df[2])
  }

  # separate first column (taxa/labels) from data
  firstcol<-df[,1]

  # rest of df is integer abundance categories
  df<-df[,-1]

  # convert log categories to integers within classes
  df[df==1]<-1
  df[df==2]<-10
  df[df==3]<-100
  df[df==4]<-1000
  df[df==5]<-10000

  # add sample labels back on
  df<-cbind.data.frame(firstcol, df)
  if (numsamples==1){
    names(df)[2]<-samplename
  }

  # return converted values
  return(df)
}

#' Check taxa against scoring list
#' @description Check the list of taxa present in the sample dataset against
#' the list of scoring taxa within package to identify any non-scoring taxa
#' in the samples (or spelling mistakes).
#' @param df A dataframe containing abundances of invertebrate taxa in
#' different samples.
#' @return A data frame containing the names of taxa that are not in the
#' list of scoring taxa, or \code{NA} if all taxa are scoring.
#' @export checktaxa
#' @examples
#' # check the taxa in the built-in Braid Burn dataset
#' # returns 'NA' if all taxa present have scores and are spelt correctly
#'
#' checktaxa(braidburn)

checktaxa<-function(df){

  # extract scoring taxa from indextable
  scoring<-as.character(indextable$Taxon)

  # create logical vector of rows with scoring taxa
  onthelist<- df[,1] %in% scoring

  # extract taxa not on the scoring list
  notonthelist<-df[!onthelist,]

  # return list of taxa not scoring, if there are any
  if (nrow(notonthelist)!=0){
    return(notonthelist[1])
  } else {
    notonthelist<-NA
    return(notonthelist)
  }
}
