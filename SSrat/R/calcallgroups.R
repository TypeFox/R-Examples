#dataframe=DF.ex6
calcallgroups <- function(dataframe, alpha = c(0.10, 0.05, 0.01), NBcriteria=F) {
  schoolgroups = unique(cbind(dataframe$schoolid, dataframe$groupid))  #sch=4;gr=20

  lastcol = ncol(dataframe)
  ratingnames = grepl("^r\\d", names(dataframe))
  # count columns with name rn
  maxrated = sum(ratingnames)
  srn = sort(names(dataframe)[ratingnames])
  # ratings may be padded with NA
  dataframe = data.frame(dataframe[!ratingnames],dataframe[srn])
  maxnrrated = maxrated - 1
  # get ratings
  grr = dataframe[, (lastcol - maxnrrated):lastcol]
  
  scale = max(grr, na.rm=T)
  switch(EXPR=scale,{scale=3},{scale=3},{scale=3},{scale=5},{scale=5},{scale=7},{scale=7},
         {scale=9},{scale=9})
  DF = data.frame()
  #i=2
  for (i in 1:nrow(schoolgroups)) {
    out = calcgroup(schoolid = schoolgroups[i, 1], groupid = schoolgroups[i,2], 
                    dataframe , scale, alpha, NBcriteria)
    DF = rbind(DF, out$dataframe)
  }
  DF
}
