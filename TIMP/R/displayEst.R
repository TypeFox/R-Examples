"displayEst" <- function (plotoptions) 
{
  if (length(plotoptions@makeps) != 0) 
    filename <- paste(plotoptions@makeps, "_paramEst.txt", sep = "")
  else
    filename <- ".currParamEst"	
       
  if(dev.interactive()) {
   writeLines("\n ### Parameter estimates ###")
    writeLines(readLines(filename))
# ## Got rid of tcltk dependency here:       
#    local({
#      tt <- tktoplevel()
#       tkwm.title(tt, "Parameter estimates")
#      txt <- tktext(tt, bg = "white", font = "courier 15")
#      scr <- tkscrollbar(tt, repeatinterval = 5, command = function(...) tkyview(txt, 
#                                                   ...))
#      tkconfigure(txt, yscrollcommand = function(...) tkset(scr, 
#                         ...))
#      tkpack(txt, side = "left", fill = "both", expand = TRUE)
#      tkpack(scr, side = "right", fill = "y")
#      chn <- tclopen(filename)
#      tkinsert(txt, "end", tclread(chn))
#      tclclose(chn)
#      tkconfigure(txt, state = "disabled")
#      tkmark.set(txt, "insert", "0.0")
#      tkfocus(txt)
#  })
  }
}
