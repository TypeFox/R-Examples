
###################################################
### code chunk number 125: Controls.Rnw:1296-1307
###################################################
offspring <- function(path = character(0), lst, ...) {
  if(length(path))
    obj <- lst[[path]]
  else
      obj <- lst
  #
  f <- function(i) is.recursive(i) && !is.null(names(i))
  data.frame(comps = names(obj), 
             hasOffspring = sapply(obj, f),
             stringsAsFactors = FALSE)
}


###################################################
### code chunk number 126: Controls.Rnw:1317-1321
###################################################
lst <- list(a = "1", b =  list(a = "2", b = "3", 
                       c = list(a = "4")))
window <- gwindow("Tree test")
tree <- gtree(offspring, offspring.data = lst, cont = window)
