# A GUI for managing and manipulating data.

ManageVariables <- function(cols, vars, query, changelog, parent=NULL) {

  ## Additional functions

  # Save changes and close GUI
  SaveChanges <- function(type) {
    SaveNb()
    if (!identical(cols, old.cols)) {
      rtn <<- list(cols=cols, vars=vars, query=query, changelog=changelog)
      old.cols <<- cols
    }
    if (type == "ok")
      tclvalue(tt.done.var) <- 1
  }

  # Set variable id and update functions to reflect this change
  SetVarId <- function(idx=NULL) {
    if (is.null(idx))
      idx <- as.integer(tkcurselection(frame1.lst)) + 1
    if (length(idx) == 0)
      return()

    # Save name
    nam <- tclvalue(name.var)
    cols[[idx]]$name <<- nam
    if (nam == "")
      nam <- "Unknown"

    # Account for duplicate ids
    new.id <- nam
    old.id <- cols[[idx]]$id
    old.ids <- vapply(cols, function(i) i$id, "")
    i <- 1L
    hold.new.id <- new.id
    while (new.id %in% old.ids[-idx]) {
      new.id <- paste0(hold.new.id, " (", i, ")")
      i <- i + 1L
    }
    cols[[idx]]$id <<- new.id
    tclvalue(list.var) <- tcl("lreplace", tclvalue(list.var),
                              idx - 1, idx - 1, new.id)

    # Update variable id's used in functions, query, and changelog
    if (!is.null(old.id)) {
      old.fun <- cols[[idx]]$fun
      str.1 <- paste0("\"", old.id, "\"")
      str.2 <- paste0("\"", new.id, "\"")
      funs <- sapply(cols, function(i) gsub(str.1, str.2, i$fun, fixed=TRUE))
      sapply(seq_along(cols), function(i) cols[[i]]$fun <<- funs[[i]])
      new.fun <- cols[[idx]]$fun
      if (!identical(old.fun, new.fun)) {
        tkconfigure(frame2.txt.4.2, state="normal")
        tcl(frame2.txt.4.2, "delete", "1.0", "end")
        tkinsert(frame2.txt.4.2, "end", new.fun)
        tkconfigure(frame2.txt.4.2, state="disabled")
      }
      if (!is.null(query))
        query <<- gsub(str.1, str.2, query, fixed=TRUE)
      if (!is.null(changelog))
        changelog[changelog[, "variable"] %in% old.id, "variable"] <<- new.id
    }
  }

  # Save notebook content
  SaveNb <- function() {
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()

    # Save format
    old.fmt <- cols[[idx]]$format
    new.fmt <- as.character(tclvalue(fmt.var))
    cols[[idx]]$format <<- new.fmt

    # Save function
    old.fun <- cols[[idx]]$fun
    new.fun <- as.character(tclvalue(tkget(frame2.txt.4.2, "1.0", "end-1c")))
    cols[[idx]]$fun <<- new.fun

    # Save summary string
    if (!identical(old.fun, new.fun))
      cols[[idx]]$summary <- summary(EvalFunction(new.fun, cols))

    # Save name
    SetVarId(idx)
  }

  # Update notebook content
  UpdateNb <- function() {
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()

    # Update name
    saved.name <- cols[[idx]]$name
    if (is.null(saved.name))
      saved.name <- ""
    tclvalue(name.var) <- saved.name

    # Update class
    tkconfigure(frame2.ent.3.2, state="normal")
    saved.class <- cols[[idx]]$class
    tclvalue(class.var) <- paste(saved.class)
    tkconfigure(frame2.ent.3.2, state="readonly")

    # Update format
    saved.fmt <- cols[[idx]]$format
    tkconfigure(frame2.ent.2.2, state="normal")
    tclvalue(fmt.var) <- saved.fmt
    tkconfigure(frame2.ent.2.2, state="readonly")

    # Update function
    tkconfigure(frame2.txt.4.2, state="normal")
    tcl(frame2.txt.4.2, "delete", "1.0", "end")
    tkinsert(frame2.txt.4.2, "end", cols[[idx]]$fun)
    tkconfigure(frame2.txt.4.2, state="disabled")
    s <- "disabled"
    if (is.na(cols[[idx]]$index))
      s <- "normal"
    tkconfigure(frame2.but.4.3, state=s)

    # Update summary
    tkconfigure(frame3.txt, state="normal")
    tcl(frame3.txt, "delete", "1.0", "end")
    if (!is.null(cols[[idx]]$summary)) {
      txt <- paste(c("", capture.output(cols[[idx]]$summary)), collapse="\n")
      tkinsert(frame3.txt, "end", txt)
    }
    tkconfigure(frame3.txt, state="disabled")
  }

  # Account for change in notebook tab
  ChangeTab <- function() {
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()
    SaveNb()
    UpdateNb()
    tabid <- tclvalue(tcl(nb, "select"))

    # Arrive at tab
    if (tabid == frame2$ID) {
      tkfocus(frame2)
    } else if (tabid == frame3$ID) {
      tkfocus(frame3)
    }
  }

  # Delete existing variable

  DeleteVar <- function() {
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()

    var.str <- paste0("\"", cols[[idx]]$id, "\"")

    if (!is.null(query) &&  grepl(var.str, query, fixed=TRUE))
      query <<- NULL
    if (!is.null(changelog) &&  grepl(var.str, changelog, fixed=TRUE))
      changelog <<- NULL

    funs.with.var <- grep(var.str, sapply(cols, function(i) i$fun), fixed=TRUE)
    dependent.vars <- funs.with.var[!funs.with.var %in% idx]

    if (length(dependent.vars) > 0) {
      ids <- vapply(cols, function(i) i$id, "")[dependent.vars]
      msg <- paste0("Variables dependent on variable \"", cols[[idx]]$id,
                    "\" include:\n\n  ", paste(ids, collapse=", "),
                    "\n\nThese variables must first be removed before this ",
                    "operation can be completed.")
      tkmessageBox(icon="error", message=msg, title="Deletion Prevented",
                   type="ok", parent=tt)
      return()
    }
    if (!is.na(cols[[idx]]$index)) {
      msg <- paste0("Variable \"", cols[[idx]]$id,
                    "\" corresponds with menu.view.unpr data.\n\n",
                    "Are you sure you want to remove it?")
      ans <- tkmessageBox(icon="question", message=msg, title="Question",
                          type="yesno", parent=tt)
      if (as.character(ans) == "no")
        return()
    }

    tclvalue(list.var) <- tcl("lreplace", tclvalue(list.var),
                              idx - 1L, idx - 1L)

    cols <<- cols[-idx]
    vars <<- vars[!vars %in% idx]

    if (length(cols) == 0)
      cols <<- NULL
    if (length(vars) == 0)
      vars <<- NULL

    for (i in seq_along(vars)) {
      if (vars[[i]] > idx)
        vars[[i]][1] <<- vars[[i]] - 1
    }

    tkselection.clear(frame1.lst, 0, "end")

    n <- length(cols)
    if (n > 0) {
      if (idx > n)
        tkselection.set(frame1.lst, idx - 2)
      else
        tkselection.set(frame1.lst, idx - 1)
      UpdateNb()
    } else {
      tclvalue(name.var) <- ""
      tclvalue(class.var) <- ""
      tclvalue(fmt.var) <- ""
      tkconfigure(frame2.txt.4.2, state="normal")
      tcl(frame2.txt.4.2, "delete", "1.0", "end")
      tkconfigure(frame2.txt.4.2, state="disabled")
      tkconfigure(frame3.txt, state="normal")
      tcl(frame3.txt, "delete", "1.0", "end")
      tkconfigure(frame3.txt, state="disabled")
    }
  }

  # Save new variable

  SaveNewVar <- function() {
    SaveNb()

    new.name <- "New Variable"
    idx <- length(cols) + 1L

    cols[[idx]] <- list(id="", class="")

    m <- Data("data.raw", which.attr="nrows")
    if (is.null(m) && length(cols) > 1)
      m <- length(EvalFunction(cols[[1]]$fun, cols))

    f <- EditFunction(cols, index=idx, value.length=m, win.title="New Variable",
                      parent=tt)
    if (is.null(f$fun) || f$fun == "")
      return()

    cols[[idx]] <<- list(id="", name="New Variable", format="", class=f$class,
                         index=NA, fun=f$fun, sample=f$sample,
                         summary=f$summary)

    tcl("lappend", list.var, new.name)
    tkselection.clear(frame1.lst, 0, "end")
    tkselection.set(frame1.lst, idx - 1L, idx - 1L)
    tkyview(frame1.lst, idx - 1L)

    UpdateNb()
    SetVarId(idx)
  }

  # Edit a variables function formula

  CallEditFunction <- function() {
    SaveNb()

    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()

    m <- Data("data.raw", which.attr="nrows")
    f <- EditFunction(cols, index=idx, value.length=m, parent=tt)

    if (is.null(f$fun))
      return()
    if (f$fun == "") {
      msg <- paste0("Nothing has been defined for this function; therefore,\n",
                    "the variable '", cols[[idx]]$name, "' will be removed.")
      ans <- as.character(tkmessageBox(icon="question", message=msg,
                                       title="Warning", type="okcancel",
                                       parent=tt))
      if (ans == "ok")
        DeleteVar()
      return()
    }

    if (!identical(f$class, cols[[idx]]$class))
        cols[[idx]]$format <<- ""

    cols[[idx]]$fun     <<- f$fun
    cols[[idx]]$class   <<- f$class
    cols[[idx]]$summary <<- f$summary
    cols[[idx]]$sample  <<- f$sample

    UpdateNb()
  }

  # Edit format

  CallFormat <- function() {
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()

    sample.value <- cols[[idx]]$sample

    old.fmt <- as.character(tclvalue(fmt.var))
    if (inherits(sample.value, c("POSIXt", "Date"))) {
      new.fmt <- FormatDateTime(sample=sample.value, fmt=old.fmt, parent=tt)
    } else {
      if (is.null(sample.value))
        sample.value <- NA
      new.fmt <- Format(sample=sample.value, fmt=old.fmt, parent=tt)
    }

    if (is.null(new.fmt))
      new.fmt <- ""

    tclvalue(fmt.var) <- new.fmt
  }

  # Arrange variables in listbox

  Arrange <- function(type) {
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()

    n <- length(cols)
    idxs <- seq_len(n)

    if (type == "back") {
      if (idx == 1)
        return()
      new.idxs <- c(idx, idxs[-idx])
      new.idx <- 1
    } else if (type == "front") {
      if (idx == n)
        return()
      new.idxs <- c(idxs[-idx], idx)
      new.idx <- n
    } else if (type == "backward") {
      if (idx == 1)
        return()
      new.idxs <- seq_len(n)
      new.idxs[c(idx - 1L, idx)] <- c(idx, idx - 1L)
      new.idx <- idx - 1L
    } else if (type == "forward") {
      if (idx == n)
        return()
      new.idxs <- seq_len(n)
      new.idxs[c(idx, idx + 1L)] <- c(idx + 1L, idx)
      new.idx <- idx + 1L
    }

    cols <<- cols[new.idxs]

    for (i in seq_along(vars)) {
      vars[[i]][1] <<- idxs[new.idxs %in% vars[[i]][1]]
    }

    ids <- vapply(cols, function(i) i$id, "")

    for (i in seq_len(n))
      tclvalue(list.var) <- tcl("lreplace", tclvalue(list.var),
                                i - 1, i - 1, ids[i])
    tkselection.clear(frame1.lst, 0, "end")
    tkselection.set(frame1.lst, new.idx - 1L)
    tkyview(frame1.lst, new.idx - 1L)
  }

  # View data for selected variable
  CallEditData <- function() {
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"))
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE), add=TRUE)
    idx <- as.integer(tkcurselection(frame1.lst)) + 1L
    if (length(idx) == 0)
      return()
    SaveNb()
    d <- list(EvalFunction(cols[[idx]]$fun, cols))
    rows <- Data("data.raw", which.attr="row.names")
    EditData(d, col.names=cols[[idx]]$id, row.names=rows,
             col.formats=cols[[idx]]$format, read.only=TRUE,
             win.title="View Raw Data", parent=tt)
    return()
  }


  ## Main program

  # Assign variables

  rtn <- NULL

  old.cols <- cols
  ids <- vapply(cols, function(i) i$id, "")

  w <- 300
  h <- 50

  # Assign the variables linked to Tk widgets

  list.var <- tclVar()
  for (i in seq_along(ids))
    tcl("lappend", list.var, ids[i])

  name.var  <- tclVar()
  fmt.var   <- tclVar()
  class.var <- tclVar()

  tt.done.var <- tclVar(0)

  # Open gui

  tclServiceMode(FALSE)
  tt <- tktoplevel()

  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }
  tktitle(tt) <- "Manage Variables"

  # Create menus

  top.menu <- tkmenu(tt, tearoff=0)

  menu.edit <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="Edit", menu=menu.edit, underline=0)
  tkadd(menu.edit, "command", label="New\u2026", accelerator="Ctrl+n",
        command=SaveNewVar)
  tkadd(menu.edit, "command", label="Delete", command=DeleteVar)

  menu.view <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="View", menu=menu.view, underline=0)
  tkadd(menu.view, "command", label="Raw data", command=CallEditData)

  menu.arrange <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Arrange", menu=menu.arrange, underline=0)
  tkadd(menu.arrange, "command", label="Send to top",
        accelerator="Shift+Ctrl+[", command=function() Arrange("back"))
  tkadd(menu.arrange, "command", label="Send upward",
        accelerator="Ctrl+[", command=function() Arrange("backward"))
  tkadd(menu.arrange, "command", label="Bring downward",
        accelerator="Ctrl+]", command=function() Arrange("forward"))
  tkadd(menu.arrange, "command", label="Bring to bottom",
        accelerator="Shift+Ctrl+]", command=function() Arrange("front"))

  tkconfigure(tt, menu=top.menu)

  # Frame 0, ok and cancel buttons, and size grip

  frame0 <- ttkframe(tt, relief="flat")
  frame0.but.1 <- ttkbutton(frame0, width=2, image=GetBitmapImage("top"),
                            command=function() Arrange("back"))
  frame0.but.2 <- ttkbutton(frame0, width=2, image=GetBitmapImage("up"),
                            command=function() Arrange("backward"))
  frame0.but.3 <- ttkbutton(frame0, width=2, image=GetBitmapImage("down"),
                            command=function() Arrange("forward"))
  frame0.but.4 <- ttkbutton(frame0, width=2, image=GetBitmapImage("bottom"),
                            command=function() Arrange("front"))
  frame0.but.5 <- ttkbutton(frame0, width=2, image=GetBitmapImage("view"),
                            command=CallEditData)
  frame0.but.6 <- ttkbutton(frame0, width=2, image=GetBitmapImage("plus"),
                            command=SaveNewVar)
  frame0.but.7 <- ttkbutton(frame0, width=2, image=GetBitmapImage("delete"),
                            command=DeleteVar)

  frame0.but.9 <- ttkbutton(frame0, width=12, text="OK",
                            command=function() SaveChanges("ok"))
  frame0.but.10 <- ttkbutton(frame0, width=12, text="Cancel",
                             command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.11 <- ttkbutton(frame0, width=12, text="Apply",
                             command=function() SaveChanges("apply"))
  frame0.but.12 <- ttkbutton(frame0, width=12, text="Help",
                             command=function() {
                               print(help("ManageVariables", package="RSurvey",
                                          verbose=FALSE))
                             })
  frame0.grp.12 <- ttksizegrip(frame0)

  tkgrid(frame0.but.1, frame0.but.2, frame0.but.3, frame0.but.4, frame0.but.5,
         frame0.but.6, frame0.but.7, "x", frame0.but.9, frame0.but.10,
         frame0.but.11, frame0.but.12, frame0.grp.12)

  tkgrid.columnconfigure(frame0, 7, weight=1)

  tkgrid.configure(frame0.but.1, frame0.but.2, frame0.but.3, frame0.but.4,
                   frame0.but.5, frame0.but.6, frame0.but.7, sticky="n",
                   padx=c(0, 2), pady=c(0, 0))
  tkgrid.configure(frame0.but.1, padx=c(10, 2))
  tkgrid.configure(frame0.but.9, frame0.but.10, frame0.but.11, frame0.but.12,
                   padx=c(0, 4), pady=c(15, 10))
  tkgrid.configure(frame0.but.12, columnspan=2, padx=c(0, 10))
  tkgrid.configure(frame0.grp.12, sticky="se")

  tkraise(frame0.but.12, frame0.grp.12)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Paned window
  pw <- ttkpanedwindow(tt, orient="horizontal")

  # Frame 1, listbox with variable names

  frame1 <- tkframe(pw, relief="flat")

  frame1.lst <- tklistbox(frame1, selectmode="browse", activestyle="none",
                          relief="flat", borderwidth=5, width=25,
                          exportselection=FALSE, listvariable=list.var,
                          highlightthickness=0)
  frame1.ysc <- ttkscrollbar(frame1, orient="vertical")
  tkconfigure(frame1.lst, background="white",
              yscrollcommand=paste(.Tk.ID(frame1.ysc), "set"))
  tkconfigure(frame1.ysc, command=paste(.Tk.ID(frame1.lst), "yview"))
  tkpack(frame1.lst, side="left",  fill="both", expand=TRUE, pady=c(2, 2))
  tkpack(frame1.ysc, side="right", fill="y", anchor="w",
         padx=c(0, 2), pady=c(2, 2))

  tkselection.set(frame1.lst, 0)

  tkadd(pw, frame1, weight=0)

  # Notebook with tabs

  nb <- ttknotebook(pw)

  # Frame 2, variable

  frame2 <- ttkframe(nb, relief="flat", padding=10, borderwidth=2)
  tkadd(nb, frame2, text="   Variable   ")

  frame2.lab.1.1 <- ttklabel(frame2, text="Name")
  frame2.lab.2.1 <- ttklabel(frame2, text="Format")
  frame2.lab.3.1 <- ttklabel(frame2, text="Class")
  frame2.lab.4.1 <- ttklabel(frame2, text="Function")

  frame2.ent.1.2 <- ttkentry(frame2, textvariable=name.var)
  frame2.ent.2.2 <- ttkentry(frame2, textvariable=fmt.var)
  frame2.ent.3.2 <- ttkentry(frame2, textvariable=class.var)

  frame2.txt.4.2 <- tktext(frame2, padx=2, pady=2, width=45, height=6,
                           undo=1, wrap="none", foreground="black",
                           background="#ebebe4", borderwidth=1,
                           font="TkFixedFont")

  frame2.but.2.3 <- ttkbutton(frame2, text="Edit", width=5,
                              command=CallFormat)
  frame2.but.4.3 <- ttkbutton(frame2, text="Edit", width=5,
                              command=CallEditFunction)

  tkgrid(frame2.lab.1.1, frame2.ent.1.2, "x")
  tkgrid(frame2.lab.2.1, frame2.ent.2.2, frame2.but.2.3)
  tkgrid(frame2.lab.3.1, frame2.ent.3.2, "x")
  tkgrid(frame2.lab.4.1, frame2.txt.4.2, frame2.but.4.3)

  tkgrid.configure(frame2.lab.1.1, frame2.lab.2.1, frame2.lab.3.1, sticky="w")

  tkgrid.configure(frame2.lab.4.1, sticky="ne")

  tkgrid.configure(frame2.ent.1.2, frame2.ent.2.2, frame2.ent.3.2,
                   sticky="we", padx=2, pady=2)

  tkgrid.configure(frame2.txt.4.2, padx=2, pady=2, sticky="nswe")

  tkgrid.configure(frame2.but.2.3, sticky="w")
  tkgrid.configure(frame2.lab.4.1, pady=c(4, 0))
  tkgrid.configure(frame2.but.4.3, sticky="nw", pady=c(1, 0))

  tkgrid.columnconfigure(frame2, 1, weight=1, minsize=25)
  tkgrid.rowconfigure(frame2, 3, weight=1, minsize=25)

  # Frame 3, summary

  frame3 <- ttkframe(nb, relief="flat", padding=0, borderwidth=0)
  tkadd(nb, frame3, text="   Summary   ")

  frame3.ysc <- ttkscrollbar(frame3, orient="vertical")

  frame3.txt <- tktext(frame3, bg="white", padx=2, pady=2, width=60, height=8,
                undo=1, wrap="none", foreground="black", relief="flat",
                font="TkFixedFont",
                yscrollcommand=function(...) tkset(frame3.ysc, ...))

  tkconfigure(frame3.ysc, command=paste(.Tk.ID(frame3.txt), "yview"))

  tkgrid(frame3.txt, frame3.ysc)

  tkgrid.configure(frame3.txt, sticky="news")
  tkgrid.configure(frame3.ysc, sticky="ns")

  tkgrid.columnconfigure(frame3, 0, weight=1, minsize=25)
  tkgrid.rowconfigure(frame3, 0, weight=1, minsize=25)

  # Insert notebook and paned window

  tkadd(pw, nb, weight=1)
  tkpack(pw, fill="both", expand="yes", padx=10, pady=c(10, 2))

  # Update Notebook
  UpdateNb()

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Control-n>", SaveNewVar)
  tkbind(tt, "<Control-]>", function() Arrange("forward"))
  tkbind(tt, "<Shift-Control-}>", function() Arrange("front"))
  tkbind(tt, "<Control-[>", function() Arrange("backward"))
  tkbind(tt, "<Shift-Control-{>", function() Arrange("back"))
  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(nb, "<<NotebookTabChanged>>", ChangeTab)

  tkbind(frame1.lst, "<ButtonPress-1>", SaveNb)
  tkbind(frame1.lst, "<Up>", SaveNb)
  tkbind(frame1.lst, "<Down>", SaveNb)
  tkbind(frame1.lst, "<<ListboxSelect>>", UpdateNb)

  tkbind(frame2.ent.1.2, "<Return>", function() SetVarId())

  # GUI control

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  return(rtn)
}
