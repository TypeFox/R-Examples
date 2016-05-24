infile1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/10001/3/6e52deaa45c1695e7742c923ba04d16b" 
infile1 <- sub("^https","http",infile1) 
gltc_values <-read.csv(infile1, 
                       col.names=c(
                     "recordID",     
                     "variable",     
                     "year",     
                     "siteID",     
                     "value"    ), check.names=TRUE,
                   stringsAsFactors = FALSE)


infile2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/10001/3/6167b9938e8dc99e9ee75251c70776a9"
infile2 <- sub("^https","http",infile2) 
gltc_metadata <-read.csv(infile2, 
                     ,quot='"' 
                     , col.names=c(
                       "siteID",     
                       "Lake.name",     
                       "Other.names",     
                       "lake.or.reservoir",     
                       "location",     
                       "region",     
                       "latitude",     
                       "longitude", 
                       "geospatial.accuracy.km",
                       "elevation.m",     
                       "mean.depth.m",     
                       "max.depth.m",     
                       "surface.area.km2",     
                       "volume.km3",     
                       "source",     
                       "sampling.depth",  
                       "sampling.time",
                       "time.period",     
                       "contributor"), check.names=TRUE,
                     stringsAsFactors = FALSE)


make_names = function(df, col_name){
  
  for (i in 1:length(col_name)){
    nm <- df[[col_name[i]]]
    nm <- gsub(pattern = '_', replacement = '.', x = nm)
    nm <- gsub(pattern = '..', replacement = '.', x = nm, fixed = T)
    nm <- gsub(pattern = '..', replacement = '.', x = nm, fixed = T)
    nm <- make.names(nm)
    nm <- pop_dot(nm)

    df[[col_name[i]]] <- nm
  }
  
  return(df)
}

pop_dot <- function(text){
  
  
  for (i in 1:length(text)){
    if (substr(text[i],nchar(text[i]),nchar(text[i])) == '.'){
      text[i] <- substr(text[i],1,(nchar(text[i]) - 1))
    }
  }
  return(text)
}
gltc_metadata <- make_names(gltc_metadata, c("Lake.name"))
gltc_values <- make_names(gltc_values, c("variable"))
save(gltc_metadata, gltc_values, file = 'sysdata.rda', compress = "xz")



