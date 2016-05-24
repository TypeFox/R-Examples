if(interactive()) {

  w <- gwindow("glayout example", visible=FALSE)
  lyt <- glayout(cont=w)

  ## character -> glabel
  lyt[1,1] <- "Label"

  ## put lyt on both sides here (parent container on right)
  lyt[1,2] <- gedit("", cont=lyt)

  ## alignment options
  lyt[2,1, expand=TRUE, anchor=c(1,0)] <- "Label 2"

  lyt[2,2] <- gslider(cont=lyt)

  ## can reference children via [
  f <- function(h,...) {
    lst <- sapply(lyt[1:2, 2], svalue)
    print(lst)
  }
  
  ## stretch
  lyt[3, 1:2] <- gbutton("click me", cont=lyt, handler=f)

  visible(w) <- TRUE


}
