if(require(RUnit)) {
test.gcombobox <- function() {
    w <- gwindow()
    g <- ggroup(cont = w, horiz = FALSE)
    
    items <- letters
    m <- data.frame(letters[1:4], rep("quit",4), toupper(letters[1:4]))
    
    ## basic
    widget <- gcombobox(items, selected = 1, cont = g)
    
    ## svalue
    checkEquals(svalue(widget), items[1])
  checkEquals(svalue(widget, index=TRUE), 1) # index=TRUE

    ## svalue<-
    svalue(widget) <- "b"
    checkEquals(svalue(widget), "b")
    
    ## index
    svalue(widget, index=TRUE) <- 3
  checkEquals(svalue(widget), items[3])
    
    
    ## [
    checkEquals(widget[], items)
    
    ## [<-
    widget[] <- m[,1, drop=TRUE]
    checkEquals(widget[], m[,1,drop=TRUE])
    
    
    ## data frame
    widget <- gcombobox(m, cont = g)
    
    svalue(widget, index=TRUE) <- 1
    checkEquals(svalue(widget, index=TRUE), 1)
    
  ## [<- with data frame
    widget[] <- m[1:2,]
    svalue(widget, index=TRUE) <- 1
    checkEquals(svalue(widget, index=TRUE), 1)
    
  }

}
