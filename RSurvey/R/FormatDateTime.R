# Build calendar date and time string formats.

FormatDateTime <- function(sample=as.POSIXct("1991-08-25 20:57:08"), fmt="",
                           parent=NULL) {

  ## Additional functions

  # Save format
  SaveFormat <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    if (as.character(tclvalue(sample.var)) == "") {
      msg <- "Format results in empty character string, please try again."
      tkmessageBox(icon="error", message=msg, title="Error", type="ok",
                   parent=tt)
      return()
    }
    fmt <- as.character(tclvalue(fmt.var))
    new.fmt <<- fmt
    tclvalue(tt.done.var) <- 1
  }

  # Selection change in the viewtree
  SelectionChange <- function() {
    cur.sel <- tcl(frame1.tre, "selection")
    cur.val <- as.character(tcl(frame1.tre, "item", cur.sel, "-values"))[1]
    tkfocus(frame2a.ent)
    if (!is.na(cur.val))
      AddString(cur.val)
  }

  # Update sample date-time entry
  UpdateSample <- function() {
    fmt <- sub("%$", "", tclvalue(fmt.var))
    tclvalue(sample.var) <- format(sample, format=fmt)
  }

  # Add string to format entry
  AddString <- function(txt) {
    if (as.logical(tcl(frame2a.ent, "selection", "present")))
      tcl(frame2a.ent, "delete", "sel.first", "sel.last")
    tkinsert(frame2a.ent, "insert", txt)
    UpdateSample()
    tkfocus(frame2a.ent)
  }

  # Expand or collapse nodes in treeview
  ToggleTreeView <- function(open.nodes) {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    if (open.nodes)
      img <- img.minus
    else
      img <- img.plus
    tkconfigure(frame0.but.1, image=img,
                command=function() ToggleTreeView(!open.nodes))
    tcl(frame1.tre, "item", id.dt, "-open", open.nodes)

    tcl(frame1.tre, "item", id.yr, "-open", open.nodes)
    tcl(frame1.tre, "item", id.mo, "-open", open.nodes)
    tcl(frame1.tre, "item", id.wk, "-open", open.nodes)
    tcl(frame1.tre, "item", id.dy, "-open", open.nodes)
    if (is.posixt) {
      tcl(frame1.tre, "item", id.tm, "-open", open.nodes)
      tcl(frame1.tre, "item", id.hr, "-open", open.nodes)
      tcl(frame1.tre, "item", id.mn, "-open", open.nodes)
      tcl(frame1.tre, "item", id.sc, "-open", open.nodes)
    }
  }

  # Copy format to clipboard
  CopyFormat <- function() {
    txt <- as.character(tclvalue(fmt.var))
    cat(txt, file="clipboard")
  }

  # Paste format from clipboard
  PasteFormat <- function() {
    cb <- try(scan(file="clipboard", what="character", sep="\n", quiet=TRUE),
              silent=TRUE)
    if (inherits(cb, "try-error"))
      return()
    tclvalue(fmt.var) <- cb
    UpdateSample()
    tkfocus(frame2a.ent)
  }

 # Clear format from entry
  ClearFormat <- function() {
    tclvalue(fmt.var) <- ""
    UpdateSample()
    tkfocus(frame2a.ent)
  }

  ## Main program

  if (!inherits(sample, c("POSIXt", "Date")))
    stop("Sample object must be of class POSIXt or Date.")
  if (!is.character(fmt))
    stop("format argument must be of class character")
  is.posixt <- inherits(sample, "POSIXt")

  cur.val <- NA
  new.fmt <- NULL
  img.plus  <- GetBitmapImage("plus")
  img.minus <- GetBitmapImage("minus")

  # Assign variables linked to Tk widgets

  sample.var <- tclVar()
  fmt.var <- tclVar(fmt)
  tt.done.var <- tclVar(0)

  UpdateSample()

  # Open GUI

  tclServiceMode(FALSE)
  tt <- tktoplevel()

  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }
  tktitle(tt) <- ifelse(is.posixt, "Format Date and Time", "Format Date")

  # Frame 0, load and cancel buttons, and size grip

  frame0 <- ttkframe(tt, relief="flat")

  frame0.but.1 <- ttkbutton(frame0, width=2, image=GetBitmapImage("plus"),
                            command=function() ToggleTreeView(TRUE))
  frame0.but.3 <- ttkbutton(frame0, width=12, text="OK", command=SaveFormat)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.5 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("FormatDateTime", package="RSurvey"))
                            })
  frame0.grp.6 <- ttksizegrip(frame0)

  tkgrid(frame0.but.1, "x", frame0.but.3, frame0.but.4, frame0.but.5,
         frame0.grp.6)

  tkgrid.columnconfigure(frame0, 1, weight=1)

  tkgrid.configure(frame0.but.1, padx=c(10, 0), pady=4, sticky="n")
  tkgrid.configure(frame0.but.3, frame0.but.4, frame0.but.5,
                   padx=c(0, 4), pady=c(15, 10))
  tkgrid.configure(frame0.but.5, columnspan=2, padx=c(0, 10))
  tkgrid.configure(frame0.grp.6, sticky="se")

  tkraise(frame0.but.5, frame0.grp.6)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Paned window
  pw <- ttkpanedwindow(tt, orient="horizontal")

  # Frame 1, treeview for conversion specifications

  frame1 <- ttkframe(pw, relief="flat", padding=0, borderwidth=0)

  frame1.tre <- ttktreeview(frame1, selectmode="browse",
                            columns=c("spec", "exam"))

  frame1.ysc <- ttkscrollbar(frame1, orient="vertical")
  tkconfigure(frame1.tre, yscrollcommand=paste(.Tk.ID(frame1.ysc), "set"))
  tkconfigure(frame1.ysc, command=paste(.Tk.ID(frame1.tre), "yview"))

  tcl(frame1.tre, "column", "#0",   width=230, anchor="center")
  tcl(frame1.tre, "column", "spec", width=80,  minwidth=80, anchor="center")
  tcl(frame1.tre, "column", "exam", width=120, minwidth=80, anchor="center")

  tcl(frame1.tre, "heading", "#0", text="Select", anchor="w")
  tcl(frame1.tre, "heading", "spec", text="Specification")
  tcl(frame1.tre, "heading", "exam", text="Example")

  id.dt <- tkinsert(frame1.tre, "", "end", tags="bg", text="date")
  id.yr <- tkinsert(frame1.tre, "", "end", tags="bg", text="year")
  id.mo <- tkinsert(frame1.tre, "", "end", tags="bg", text="month")
  id.wk <- tkinsert(frame1.tre, "", "end", tags="bg", text="week")
  id.dy <- tkinsert(frame1.tre, "", "end", tags="bg", text="day")
  tkinsert(frame1.tre, id.dt, "end", tags="bg", text="month day year",
           values=c("%m/%d/%Y", format(sample, format="%m/%d/%Y")))
  tkinsert(frame1.tre, id.dt, "end", tags="bg", text="month day year",
           values=c("%m/%d/%y", format(sample, format="%m/%d/%y")))
  tkinsert(frame1.tre, id.dt, "end", tags="bg", text="year month day",
           values=c("%Y-%m-%d", format(sample, format="%Y-%m-%d")))
  tkinsert(frame1.tre, id.dt, "end", tags="bg", text="day month year",
           values=c("%d%b%Y", format(sample, format="%d%b%Y")))
  tkinsert(frame1.tre, id.dt, "end", tags="bg", text="weekday month day",
           values=c("%A %B %d", format(sample, format="%A %B %d")))
  tkinsert(frame1.tre, id.dt, "end", tags="bg", text="weekday month day",
           values=c("%a %b %d", format(sample, format="%a %b %d")))
  tkinsert(frame1.tre, id.yr, "end", tags="bg",
           text="year without century (00-99)",
           values=c("%y", format(sample, format="%y")))
  tkinsert(frame1.tre, id.yr, "end", tags="bg", text="year with century",
           values=c("%Y", format(sample, format="%Y")))
  tkinsert(frame1.tre, id.mo, "end", tags="bg", text="month (01-12)",
           values=c("%m", format(sample, format="%m")))
  tkinsert(frame1.tre, id.mo, "end", tags="bg", text="abbreviated month name",
           values=c("%b", format(sample, format="%b")))
  tkinsert(frame1.tre, id.mo, "end", tags="bg", text="full month name",
           values=c("%B", format(sample, format="%B")))
  tkinsert(frame1.tre, id.wk, "end", tags="bg",
           text="week of the year (00-53), US",
           values=c("%U", format(sample, format="%U")))
  tkinsert(frame1.tre, id.wk, "end", tags="bg",
           text="week of the year (00-53), UK",
           values=c("%W", format(sample, format="%W")))
  tkinsert(frame1.tre, id.dy, "end", tags="bg", text="day of the month (01-31)",
           values=c("%d", format(sample, format="%d")))
  tkinsert(frame1.tre, id.dy, "end", tags="bg",
           text="day of the year (001-366)",
           values=c("%j", format(sample, format="%j")))
  tkinsert(frame1.tre, id.dy, "end",  tags="bg",
           text="weekday (0-6, Sunday is 0)",
           values=c("%w", format(sample, format="%w")))
  tkinsert(frame1.tre, id.dy, "end", tags="bg", text="abbreviated weekday name",
           values=c("%a", format(sample, format="%a")))
  tkinsert(frame1.tre, id.dy, "end", tags="bg", text="full weekday name",
           values=c("%A", format(sample, format="%A")))

  if (is.posixt) {
    id.tm <- tkinsert(frame1.tre, "", "end", tags="bg", text="time")
    id.hr <- tkinsert(frame1.tre, "", "end", tags="bg", text="hour")
    id.mn <- tkinsert(frame1.tre, "", "end", tags="bg", text="minute")
    id.sc <- tkinsert(frame1.tre, "", "end", tags="bg", text="second")
    id.tz <- tkinsert(frame1.tre, "", "end", tags="bg", text="timezone")
    tkinsert(frame1.tre, id.tm, "end", tags="bg", text="hour minute second",
             values=c("%H:%M:%S", format(sample, format="%H:%M:%S")))
    tkinsert(frame1.tre, id.tm, "end", tags="bg", text="hour minute fractional-second",
             values=c("%H:%M:%OS3", format(sample, format="%H:%M:%OS3")))
    tkinsert(frame1.tre, id.tm, "end", tags="bg", text="hour minute",
             values=c("%I:%M %p", format(sample, format="%I:%M %p")))
    tkinsert(frame1.tre, id.hr, "end", tags="bg", text="hours (00-23)",
             values=c("%H", format(sample, format="%H")))
    tkinsert(frame1.tre, id.hr, "end", tags="bg", text="hours (01-12)",
             values=c("%I", format(sample, format="%I")))
    tkinsert(frame1.tre, id.hr, "end", tags="bg", text="AM/PM indicator",
             values=c("%p", format(sample, format="%p")))
    tkinsert(frame1.tre, id.mn, "end", tags="bg", text="minute (00-59)",
             values=c("%M", format(sample, format="%M")))
    tkinsert(frame1.tre, id.sc, "end", tags="bg", text="second (00-61)",
             values=c("%S", format(sample, format="%S")))
    tkinsert(frame1.tre, id.sc, "end",  tags="bg",
             text="decisecond precision",
             values=c("%OS1", format(sample, format="%OS1")))
    tkinsert(frame1.tre, id.sc, "end",  tags="bg",
             text="centisecond precision",
             values=c("%OS2", format(sample, format="%OS2")))
    tkinsert(frame1.tre, id.sc, "end",  tags="bg",
             text="millisecond precision",
             values=c("%OS3", format(sample, format="%OS3")))
     tkinsert(frame1.tre, id.tz, "end",  tags="bg",
             text="time zone (output only)",
             values=c("%Z", format(sample, format="%Z")))
  }

  tktag.configure(frame1.tre, "bg", background="white")

  # Frame 2

  frame2 <- ttkframe(pw, relief="flat")

  frame2a <- ttklabelframe(frame2, relief="flat", borderwidth=5, padding=5,
                           text="Conversion specification format")

  frame2a.ent <- ttkentry(frame2a, textvariable=fmt.var, width=30)
  tkicursor(frame2a.ent, "end")

  frame2a.but.1 <- ttkbutton(frame2a, width=2, text="/",
                             command=function() AddString("/"))
  frame2a.but.2 <- ttkbutton(frame2a, width=2, text="-",
                             command=function() AddString("-"))
  frame2a.but.3 <- ttkbutton(frame2a, width=2, text=",",
                             command=function() AddString(","))
  frame2a.but.4 <- ttkbutton(frame2a, width=2, text=":",
                             command=function() AddString(":"))
  frame2a.but.5 <- ttkbutton(frame2a, width=2, text=" ",
                             command=function() AddString(" "))

  frame2a.but.6 <- ttkbutton(frame2a, width=2, image=GetBitmapImage("copy"),
                             command=CopyFormat)
  frame2a.but.7 <- ttkbutton(frame2a, width=2, image=GetBitmapImage("paste"),
                             command=PasteFormat)
  frame2a.but.8 <- ttkbutton(frame2a, width=2, image=GetBitmapImage("delete"),
                             command=ClearFormat)

  tkgrid(frame2a.ent, pady=c(5, 0))
  tkgrid(frame2a.but.1, frame2a.but.2, frame2a.but.3, frame2a.but.4,
         frame2a.but.5, frame2a.but.6, frame2a.but.7, frame2a.but.8, "x",
         padx=c(0, 2), pady=c(8, 5), sticky="w")
  tkgrid.configure(frame2a.ent, sticky="we", columnspan=9)

  tkgrid.configure(frame2a.but.6, padx=c(5, 2))
  tkgrid.configure(frame2a.but.8)

  tcl("grid", "anchor", frame2a, "w")
  tkgrid.columnconfigure(frame2a, 8, weight=1, minsize=0)
  tkpack(frame2a, fill="x", padx=c(5, 0), pady=c(5, 2))

  frame2b <- ttklabelframe(frame2, relief="flat", borderwidth=5, padding=5,
                           text="Sample")
  frame2b.ent <- ttkentry(frame2b, textvariable=sample.var, width=30,
                          state="readonly", takefocus=FALSE)
  tkgrid(frame2b.ent, pady=5)
  tkgrid.configure(frame2b.ent, sticky="we")
  tcl("grid", "anchor", frame2b, "w")
  tkgrid.columnconfigure(frame2b, 0, weight=1, minsize=10)
  tkpack(frame2b, fill="x", padx=c(5, 0), pady=c(0, 2))

  frame2c <- ttklabelframe(frame2, relief="flat", borderwidth=5, padding=5,
                           text="Example")
  fg <- "#414042"
  fmt <- ifelse(is.posixt, "%Y-%m-%d %H:%M:%S", "%Y-%m-%d")
  frame2c.lab.1 <- ttklabel(frame2c, foreground=fg, text=fmt)
  frame2c.lab.2 <- ttklabel(frame2c, foreground=fg,
                            text=format(sample, format=fmt))
  tkgrid(frame2c.lab.1, padx=0, pady=c(5, 1))
  tkgrid(frame2c.lab.2, padx=0, pady=c(1, 5))
  tcl("grid", "anchor", frame2c, "w")
  tkpack(frame2c, fill="x", padx=c(5, 0), pady=c(5, 0))

  # Layout paned window

  tkgrid(frame1.tre, frame1.ysc, frame2)

  tkgrid.configure(frame1.ysc, pady=c(20, 0), sticky="nws")
  tkgrid.configure(frame1.tre, sticky="news")

  tkgrid.rowconfigure(frame1, frame1.tre, weight=1)
  tkgrid.columnconfigure(frame1, frame1.tre, weight=1)

  tkadd(pw, frame1, weight=2)
  tkadd(pw, frame2, weight=0)

  tkpack(pw, fill="both", expand=TRUE, padx=10, pady=c(10, 0))

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)
  tkbind(frame1.tre, "<<TreeviewSelect>>", SelectionChange)
  tkbind(frame2a.ent, "<KeyRelease>", UpdateSample)

  # GUI control

  tkfocus(frame2a.ent)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  return(new.fmt)
}
