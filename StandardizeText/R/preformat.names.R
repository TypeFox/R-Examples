preformat.names <-
function(x) {
  #Geenral formatting
  #Convert to lowercase
  x <- tolower(x)
  #Replace spaces with underscore
  x <- gsub("\\s+", '_', x)
  #Remove non alphanumeric chars
  x <- gsub("\\W+", '', x)
  #Replace and/of/the with underscore
  x <- gsub('(_|\\b)and_the(_|\\b)', '_', x)
  x <- gsub('(_|\\b)of_the(_|\\b)', '_', x)
  x <- gsub('(_|\\b)and(_|\\b)', '_', x)
  x <- gsub('(_|\\b)of(_|\\b)', '_', x)
  x <- gsub('(_|\\b)the(_|\\b)', '_', x)
  #Replace multiple underscores with one
  x <- gsub('_+', '_', x)
  #Remove underscores at the beginning or end of the name
  x <- gsub('(\\b_+)|(_+\\b)', '', x)
  #Country specific formatting
  #Standardize country abbreviations
  x <- gsub('democratic', 'dem', x)
  x <- gsub('republic', 'rep', x)
  x <- gsub('people_?s_dem_rep|dem_people_?s_rep', 'pdr', x)
  x <- gsub('(federated|federal)', 'fed', x)
  x <- gsub('yugoslav(ia)?', 'yug', x)
  x <- gsub('former_yug_rep', 'fyr', x)
  x <- gsub('special_admini?strative_region', 'sar', x)
  x <- gsub('^congo$','congo_rep',x)
  return(x)
}
