thisFileName <-
function(){
  return(gsub("\\.Rmd$","",thisfile_knit()) )
}
