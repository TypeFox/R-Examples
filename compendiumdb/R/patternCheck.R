## Function to create tags for improperly formatted sample characteristics (missing a ':')
patternCheck <- function(y){
  t <- unlist(strsplit(y,";;"))
  indx <-  setdiff(1:length(t),grep(":",t))
  for (i in (1:length(indx))) t[indx[i]] <- paste0("tag",i,":",t[indx[i]])
  paste0(t,collapse=";;")
}
