enterratdata <- function(dataframe=NULL, resplabel=F){
  if (is.null(dataframe)) {
    if (resplabel==T) {
    dataframe=data.frame(resplabel=character(), schoolid=numeric(0), groupid=numeric(0), respid=numeric(0),
                      r01=numeric(0), r02=numeric(0), r03=numeric(0), stringsAsFactors=FALSE)}
    else {dataframe=data.frame(schoolid=numeric(0), groupid=numeric(0), respid=numeric(0),
                              r01=numeric(0), r02=numeric(0), r03=numeric(0), stringsAsFactors=FALSE)}
  } else {
    if (!("resplabel" %in% colnames(dataframe)) & resplabel){
      # add column resplabel
      resplabel=rep("", nrow(dataframe))
      dataframe=data.frame(resplabel, dataframe, stringsAsFactors=FALSE)
    }
  } 
  #suppressWarnings(edit(dataframe))
  edit(dataframe)
}

# #create a new data frame with rating data
# df=enterratdata()
# 
# # edit existing data frame
# df=enterratdata(df)
# 
# # add respondent names to the dataframe
# df=enterratdata(df, resplabel=T)
# 
# #create a new data frame with rating data and respondent names
# df=enterratdata(resplabel=T)
