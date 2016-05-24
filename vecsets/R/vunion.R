vunion <-
function (x, y, multiple=TRUE) {
trueun <- c(vintersect(x,y), vsetdiff(x,y), vsetdiff(y,x) )
if( !multiple) trueun <- unique(trueun)
return(trueun)
}
