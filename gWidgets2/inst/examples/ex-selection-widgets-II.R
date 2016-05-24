if(interactive()) {
  ## a tour of the selection widgets

  w <- gwindow("Selection widgets")
  lyt <- glayout(cont=w)

  ## checkbox
  lyt[1,1] <- "checkbox"
  lyt[1,2] <- gcheckbox("checkbox", checked=FALSE, cont=lyt)

  ## radio button group
  lyt[2,1] <- "radio buttons"
  lyt[2,2] <- gradio(state.name[1:3], selected=2, cont=lyt, horizontal=TRUE)
  
  ## checkbox group
  lyt[3,1] <- "checkbox group"
  lyt[3,2] <- gcheckboxgroup(state.name[1:3], checked=c(TRUE, FALSE, TRUE), cont=lyt, horizontal=TRUE)
  
  ## checkbox group, using atable
  lyt[4,1] <- "checkbox group"
  lyt[4,2] <- gcheckboxgroup(state.name[1:3], checked=c(TRUE, FALSE, TRUE), cont=lyt, use.table=TRUE)
  
  ## combobox (drop list)
  lyt[5,1] <- "combobox"
  lyt[5,2] <- gcombobox(state.name, selected=match("New York", state.name), cont=lyt)
  
  
  ## editable combobox
  lyt[6,1] <- "editable combobox"
  lyt[6,2] <- gcombobox(state.name, selected=match("New York", state.name), cont=lyt, editable=TRUE)

  lyt[7,1:2] <- gseparator(cont=lyt)

  lyt[8,2] <- gbutton("values", cont=lyt, handler=function(h,...) {
    print(sapply(lyt[1:6, 2], svalue))
  })
                      
  visible(w) <- TRUE
}
