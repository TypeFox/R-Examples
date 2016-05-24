cleanFun <-
function(htmlString) {
  return(gsub("<.*?>", "", htmlString))
}
