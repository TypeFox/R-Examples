sws <-
function(w,only.possible=TRUE,check.valid=FALSE){
  #Scrabble words have no case
  w <- tolower(w)
  wv <- strsplit(w,"")
  init.score <- sapply(lapply(wv,sls),sum)
  if(only.possible){
    init.score <- init.score - impossible.points(wv)
  }
  if(check.valid){
    init.score[!is.twl06.word(w)] <- 0
  }
  return(init.score)
}
