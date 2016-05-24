#######################################################################
##
## Function: replace.list()
## Author  : Jonathan Wand, http://wand.stanford.edu
## Created : 2002-10-01
##
## Function to replace selected items in lists (or lists of lists)
## with contents of another list.  Note: this function is recursive.
## 
## INPUT:
##   old = original list
##   new = contains only stuff that you want to have inserted/replaced
##         into 'old' list
## OUTPUT:
##   old list in updated form
##
## MODIFIED
##    2008-05-01 : JW
##    - added is.list(tmp.new) to recursive case
##      this ensures that the structure of the new list dominates
#######################################################################
replace.list <- function( old,new ) {

  for (n.new in names(new)) {
    tmp.old <- eval(parse(text=paste("old$",n.new,sep="")))
    tmp.new <- eval(parse(text=paste("new$",n.new,sep="")))
    if (is.list(tmp.old) && is.list(tmp.new) ) {
      tmp.new <- replace.list(tmp.old,tmp.new)
    } 
    eval(parse(text=paste("old$",n.new,"<- tmp.new",sep="")))
  }
  return( old )
}
