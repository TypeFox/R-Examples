## Helper function for 'setData'

defineTable <-
function(frame, var){
  table <- tkwidget(frame, "table", resizeborders = "none",
                    selectmode = "extended", multiline = "0",
                    rowseparator = "\n", colseparator = "\t",
                    variable = var, rows = 6, cols = 7, colwidth = 11,
                    titlerows = 1, titlecols = 1, background = "white")
  return(table)
}