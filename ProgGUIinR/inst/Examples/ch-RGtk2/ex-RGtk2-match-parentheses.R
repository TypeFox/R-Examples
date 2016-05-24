### R code from vignette source 'ex-RGtk2-match-parentheses.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-match-parentheses.Rnw:17-21
###################################################
window <- gtkWindowNew()
view <- gtkTextView()
window$add(view)
buffer <- view$getBuffer()


###################################################
### code chunk number 2: ex-RGtk2-match-parentheses.Rnw:24-32
###################################################
gtkTextIterGetTextAtIter <- function(iter) {
  siter <- iter$copy();  eiter <- iter$copy()
  val <- siter$backwardChar()
  if(val) 
    siter$getText(eiter)                # GetText method of an iter
  else 
    return(NA)
}


###################################################
### code chunk number 3: ex-RGtk2-match-parentheses.Rnw:43-73
###################################################
## match up tokens
tokenMatch <- list("("=")",")"="(",
                   "["="]","]"="[",
                   "{"="}","}"="{")

findLeftMatch <- function(iter) {
                    
  titer <- iter$copy()
  token = titer$getTextAtIter()

  if(!token %in% c(")","}","]"))
    return(list(retval=FALSE, iter= NA))

  lt <- rt <- 0                      # count left, right tokens
  while(!is.na(char <- titer$getTextAtIter())) {
    if(char == token)
      rt <- rt + 1
    else if(char == tokenMatch[[token]])
      lt <- lt + 1
    if(lt == rt) {
      titer$backwardChar()              # step back 1
      return(list(retval=TRUE, iter=titer))
    } else if(lt > rt) {
      return(list(retval=FALSE, iter=NA))
    }

    titer$backwardChar()
  }
  return(list(retval=FALSE, iter=NA))   # no match
}  


###################################################
### code chunk number 4: ex-RGtk2-match-parentheses.Rnw:78-80
###################################################
tag.highlight <-
  buffer$createTag(tag.name="highlight", background="yellow")


###################################################
### code chunk number 5: ex-RGtk2-match-parentheses.Rnw:90-104
###################################################
ID <- 
  gSignalConnect(buffer,"changed",
                 f = function(buffer, data) {
                   ## get iter at cursor
                   iter <- buffer$getIterAtMark(buffer$getInsert())$iter
                   char <- iter$getTextAtIter()
                   val <- findLeftMatch(iter)
                   if(val$retval) {
                     buffer$applyTagByName("highlight",val$iter,iter)
                   } else {
                     bnds <- buffer$getBounds() # remove all
                     buffer$removeTagByName("highlight",bnds$start, bnds$end)
                   }
                 })


