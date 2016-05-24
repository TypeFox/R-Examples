# A GUI for viewing or editing text.

EditText <- function(txt, read.only=FALSE, win.title="View Text",
                     is.fixed.width.font=FALSE, parent=NULL) {

  ## Additional functions

  # Close GUI and return edited text
  SaveText <- function() {
    txt <- as.character(tclvalue(tkget(frame1.txt.1.1, "1.0", "end-1c")))
    txt <- strsplit(txt, split="\n", fixed=TRUE)[[1]]
    rtn <<- txt
    tclvalue(tt.done.var) <- 1
  }

  # Open file in text console
  OpenFile <- function(is.appended=FALSE) {
    txt <- as.character(tclvalue(tkget(frame1.txt.1.1, "1.0", "end-1c")))
    if (txt != "" & !is.appended) {
      msg <- paste("This action will delete existing console text.",
                   "Would you like to continue?", sep="\n")
      ans <- as.character(tkmessageBox(icon="question", message=msg,
                                       title="Warning", type="yesno",
                                       parent=tt))
      if (ans == "yes")
        ClearConsole()
      else
        return()
    }
    f <- GetFile(cmd="Open", exts="txt", win.title="Open Text File",
                 parent=tt)
    if (is.null(f))
      return()
    txt <- paste(readLines(f), collapse="\n")
    tkinsert(frame1.txt.1.1, "end", txt)
  }

  # Save current text to file
  SaveAs <- function() {
    txt <- as.character(tclvalue(tkget(frame1.txt.1.1, "1.0", "end-1c")))
    f <- GetFile(cmd="Save As", exts="txt", win.title="Save Text As",
                 defaultextension="txt", parent=tt)
    if (is.null(f))
      return()
    cat(txt, file=f, sep="\n")
  }

  # Edit menu functions
  EditUndo <- function() {
    tkfocus(frame1.txt.1.1)
    try(tcl(frame1.txt.1.1, "edit", "undo"), silent=TRUE)
  }
  EditRedo <- function() {
    tkfocus(frame1.txt.1.1)
    try(tcl(frame1.txt.1.1, "edit", "redo"), silent=TRUE)
  }
  EditCut <- function() {
    tkfocus(frame1.txt.1.1)
    tcl("tk_textCut", frame1.txt.1.1)
  }
  EditCopy <- function() {
    tkfocus(frame1.txt.1.1)
    tcl("tk_textCopy", frame1.txt.1.1)
  }
  EditPaste <- function() {
    tkfocus(frame1.txt.1.1)
    tcl("tk_textPaste", frame1.txt.1.1)
  }
  EditSelectAll <- function() {
    tkfocus(frame1.txt.1.1)
    tktag.add(frame1.txt.1.1, "sel", "1.0", "end")
  }
  ClearConsole <- function() {
    tcl(frame1.txt.1.1, "delete", "1.0", "end")
    tkfocus(frame1.txt.1.1)
  }

  ## Main program

  # Assign missing values
  if (missing(txt) || is.null(txt) || length(txt) == 0)
    txt <- ""
  if (!is.character(txt))
    stop("input text argument is not of class character")

  # Add end-of-line for vector of character strings and
  # determine the maximum number of characters in a line
  if (length(txt) > 1) {
    txt <- paste(txt, collapse="\n")
    n <- max(vapply(strsplit(txt, split="\n", fixed=TRUE)[[1]], nchar, 0L))
  } else {
    n <- 0
  }

  # Determine the width of the text window
  txt.width <- 80

  # Determine font type
  if (is.fixed.width.font)
    font.type <- "TkFixedFont"
  else
    font.type <- "TkTextFont"

  # Assigin global variables
  rtn <- NULL

  # Assign variables linked to Tk widgets
  tt.done.var <- tclVar(0)

  # Open GUI
  tclServiceMode(FALSE)
  tt <- tktoplevel()
  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }
  tktitle(tt) <- win.title

  # Start top menu
  top.menu <- tkmenu(tt, tearoff=0)

  # Edit menu

  menu.file <- tkmenu(tt, tearoff=0, relief="flat")
  menu.edit <- tkmenu(tt, tearoff=0, relief="flat")

  tkadd(top.menu, "cascade", label="File", menu=menu.file, underline=0)
  tkadd(top.menu, "cascade", label="Edit", menu=menu.edit, underline=0)


  if (read.only) {
    tkadd(menu.file, "command", label="Save as\u2026", accelerator="Ctrl+s",
          command=SaveAs)
    tkadd(menu.edit, "command", label="Select all", accelerator="Ctrl+a",
          command=EditSelectAll)
    tkadd(menu.edit, "separator")
    tkadd(menu.edit, "command", label="Copy", accelerator="Ctrl+c",
          command=EditCopy)
  } else {
    tkadd(menu.file, "command", label="Open\u2026", accelerator="Ctrl+o",
          command=function() OpenFile())
    tkadd(menu.file, "command", label="Open and append\u2026",
          command=function() OpenFile(is.appended=TRUE))
    tkadd(menu.file, "command", label="Save as\u2026", accelerator="Ctrl+s",
          command=SaveAs)

    tkadd(menu.edit, "command", label="Undo", accelerator="Ctrl+z",
          command=EditUndo)
    tkadd(menu.edit, "command", label="Redo", accelerator="Ctrl+y",
          command=EditRedo)
    tkadd(menu.edit, "separator")
    tkadd(menu.edit, "command", label="Cut", accelerator="Ctrl+x",
          command=EditCut)
    tkadd(menu.edit, "command", label="Copy", accelerator="Ctrl+c",
          command=EditCopy)
    tkadd(menu.edit, "command", label="Paste", accelerator="Ctrl+v",
          command=EditPaste)
    tkadd(menu.edit, "separator")
    tkadd(menu.edit, "command", label="Select all", accelerator="Ctrl+a",
          command=EditSelectAll)
    tkadd(menu.edit, "command", label="Clear console", accelerator="Ctrl+l",
          command=ClearConsole)
  }

  # Finish top menu
  tkconfigure(tt, menu=top.menu)

  # Frame 0, ok and cancel buttons

  frame0 <- tkframe(tt, relief="flat")

  if (read.only) {
    frame0.but.1.2 <- "x"
    frame0.but.1.3 <- ttkbutton(frame0, width=12, text="Close",
                            command=function() tkdestroy(tt))
  } else {
    frame0.but.1.2 <- ttkbutton(frame0, width=12, text="OK",
                                command=SaveText)
    frame0.but.1.3 <- ttkbutton(frame0, width=12, text="Cancel",
                                command=function() tclvalue(tt.done.var) <- 1)
  }
  frame0.grp.1.4 <- ttksizegrip(frame0)

  tkgrid("x", frame0.but.1.2, frame0.but.1.3, frame0.grp.1.4)

  tkgrid.configure(frame0.but.1.3, columnspan=2, padx=c(4, 10), pady=10)
  tkgrid.configure(frame0.grp.1.4, sticky="se")

  tkraise(frame0.but.1.3, frame0.grp.1.4)

  tkgrid.columnconfigure(frame0, 0, weight=1)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1, text window

  frame1 <- tkframe(tt, relief="flat")

  frame1.txt.1.1 <- tktext(frame1, bg="white", font=font.type, padx=2,
                           pady=2, width=txt.width, height=20, undo=1,
                           autoseparators=1, wrap="none", foreground="black",
                           relief="flat",
                           yscrollcommand=function(...)
                                            tkset(frame1.ysc.1.2, ...),
                           xscrollcommand=function(...)
                                            tkset(frame1.xsc.2.1, ...))

  frame1.ysc.1.2 <- ttkscrollbar(frame1, orient="vertical")
  frame1.xsc.2.1 <- ttkscrollbar(frame1, orient="horizontal")

  tkconfigure(frame1.ysc.1.2, command=paste(.Tk.ID(frame1.txt.1.1), "yview"))

  tkgrid(frame1.txt.1.1, frame1.ysc.1.2)
  tkgrid.configure(frame1.txt.1.1, sticky="nswe")
  tkgrid.configure(frame1.ysc.1.2, sticky="ns")

  if (!read.only || (read.only & n > 120)) {
    tkconfigure(frame1.xsc.2.1, command=paste(.Tk.ID(frame1.txt.1.1), "xview"))
    tkgrid(frame1.xsc.2.1, "x")
    tkgrid.configure(frame1.xsc.2.1, sticky="we")
  }

  tkgrid.columnconfigure(frame1, 0, weight=1)
  tkgrid.rowconfigure(frame1, 0, weight=1)

  tkpack(frame1, fill="both", expand="yes")

  tkinsert(frame1.txt.1.1, "end", txt)
  tcl(frame1.txt.1.1, "edit", "reset")
  tcl(frame1.txt.1.1, "edit", "separator")
  if (read.only)
    tkconfigure(frame1.txt.1.1, state="disabled")

  # Bind events

  tclServiceMode(TRUE)

  if (!read.only) {
    tkbind("Text", "<Control-s>", SaveAs)
    tkbind("Text", "<Control-z>", EditUndo)
    tkbind("Text", "<Control-y>", EditRedo)
    tkbind("Text", "<Control-v>", EditPaste)
    tkbind("Text", "<Control-l>", ClearConsole)
  }
  tkbind("Text", "<Control-a>", EditSelectAll)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  # GUI control

  tkfocus(frame1.txt.1.1)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  invisible(rtn)
}
