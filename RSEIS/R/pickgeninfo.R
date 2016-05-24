pickgeninfo<-function()
  {
doc = c("#############################", "swig is a generic Picking program",
"Use left Mouse click to select traces and windows.",
"Use the buttons to operate on the selected traces/windows",
"To zoom, click twice with the left mouse on the traces, then right mouse zooms",
"Right Mouse with no left mouse clicks is equivalent to \"Done\" ",
  "Some Buttons require selection of a trace and/or a window",
"#############################"
  )
cat(doc, sep="\n")
  }


