simplifyPath =
# Could use strsplit, etc.
#   simplifyPath2("XMLBasics/../WebTechData/RawData/KyphosisRpartExtract.xml")
#   simplifyPath2("../WebTechData/RawData/KyphosisRpartExtract.xml")
#   simplifyPath2("../../WebTechData/RawData/KyphosisRpartExtract.xml")  
#   simplifyPath2("a/b/../../WebTechData/RawData/KyphosisRpartExtract.xml")
#   simplifyPath2("top/a/b/../../WebTechData/RawData/KyphosisRpartExtract.xml")  
#   simplifyPath2("abc/../../WebTechData/RawData/KyphosisRpartExtract.xml")    
function(path)
{
 els = strsplit(path, "/")[[1]]
 GoOn = TRUE

 els = els[ els != "."]
 
 while(GoOn && length(i <- which(els == ".."))) {
   i = min(i)
   if(length(i) == 1 && i == 1)
     break

   if(all(els[ seq( 1, i) ] == ".."))
     break
   if(i == 2 && els[1] == "..")
     break

   els = els[ - c(i, i - 1L) ]
 }

 paste(els, collapse = "/")
}



if(FALSE) {
simplifyPath =
function(path)
{
  path = gsub("/\\./", "/", path)  
  path = gsub("^\\./", "", path)  
  # Doesn't  handle "../foo"
  while(grepl("[^./]/\\.\\.", path)) {
     path = gsub("/[^/.]+/\\.\\./?", "/", path)
 }

  
 path = gsub("^(\\./)+", "", path)
  
 path
}


simplifyPath1 =
# Could use strsplit, etc.  
function(path)
{
 els = strsplit(path, "/")[[1]]
 
 while(length(i <- which(els == ".."))) {
   i = max(i)
   if(length(i) == 1 && i == 1)
     break

   i = i[i != 1]
 }

 paste(els, sep = "/")
}
 
}
