### R code from vignette source 'ex-tcltk-text.Rnw'

###################################################
### code chunk number 1: makeGUI
###################################################
w <- tktoplevel(); tkwm.title(w, "Text buffer example")
f <- ttkframe(w, padding = c(3,3,12,12))
tkpack(f, expand = TRUE, fill = "both")
txt <- tktext(f, width = 80, height = 24)   # default size
addScrollbars(f, txt)


###################################################
### code chunk number 2: ex-tcltk-text.Rnw:24-29
###################################################
tktag.configure(txt, "commandTag", foreground = "blue", 
                font = "courier 12 italic")
tktag.configure(txt, "outputTag", font = "courier 12")
tktag.configure(txt, "errorTag", foreground = "red", 
                font = "courier 12 bold")


###################################################
### code chunk number 3: ex-tcltk-text.Rnw:36-60
###################################################
eval_cmd_chunk <- function(txt, cmds) {
  
  cmd_chunks <- try(parse(text = cmds), silent = TRUE)
  if(inherits(cmd_chunks,"try-error")) {
    tkinsert(t, "end", "Error", "errorTag") # add markup tag
  }

  for(cmd in cmd_chunks) {
    cutoff <- 0.75 * getOption("width")
    dcmd <- deparse(cmd, width.cutoff = cutoff)
    command <- 
      paste(getOption("prompt"),
            paste(dcmd, collapse = paste("\n", 
                          getOption("continue"), sep = "")),
            sep = "", collapse = "")
    tkinsert(txt, "end", command, "commandTag")
    tkinsert(txt, "end","\n")
    ## output, should check for errors in eval!
    output <- capture.output(eval(cmd, envir = .GlobalEnv))
    output <- paste(output, collapse = "\n")
    tkinsert(txt, "end", output, "outputTag")
    tkinsert(txt, "end","\n")
  }
}


###################################################
### code chunk number 4: ex-tcltk-text.Rnw:63-80
###################################################
## function to add scrollbars to a widget
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}


###################################################
### code chunk number 5: TryIt
###################################################
eval_cmd_chunk(txt, "2 + 2; lm(mpg ~ wt, data = mtcars)")


