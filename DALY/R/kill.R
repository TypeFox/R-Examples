## Disable table cells

kill <-
function(thisTbl, thisVar, y, x){
  for (Y in y){
    for (X in x){
      thisVar[[Y, X]] <- NULL
      tcl(.Tk.ID(thisTbl), "tag", "celltag", "dead",
          paste(Y, ",", X, sep = ""))
      tcl(.Tk.ID(thisTbl), "tag", "configure", "dead",
          state = "disabled", bg = "#777")
    }
  }
}