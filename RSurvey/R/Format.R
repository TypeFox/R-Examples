# Build C-style string formats.

Format <- function(sample=pi, fmt="", parent=NULL) {

  ## Additional functions

  # Save conversion specification format
  SaveFormat <- function() {
    opt <- as.integer(tclvalue(opt.var))
    fmt <- as.character(tclvalue(fmt.var))
    sam <- as.character(tclvalue(sample.var))
    if (opt == 3L && sam == "") {
      msg <- paste0("Invalid format '", fmt, "'; please try again.")
      tkmessageBox(icon="error", message=msg, title="Error", type="ok",
                   parent=tt)
      tkfocus(frame3.ent.1)
    } else {
      new.fmt <<- ifelse(opt == 1L, "", fmt)
      tclvalue(tt.done.var) <- 1
    }
  }

  # Translate format string

  TranslateFormat <- function() {
    tclvalue(opt.var) <- 1L
    if (nchar(fmt) == 0L)
      return()
    tclvalue(opt.var) <- 3L

    if (nchar(fmt) < 2L || substr(fmt, 1L, 1L) != "%")
      return()
    code <- substr(fmt, nchar(fmt), nchar(fmt))
    if (inherits(sample, "integer")) {
      if (code != "d")
        return()
    } else if (inherits(sample, "numeric")) {
      if (!code %in% c("f", "e"))
        return()
    } else {
      if(code != "s")
        return()
    }

    width <- ""
    precision <- ""
    is.scientific <- ifelse(code == "e", TRUE, FALSE)

    len <- attr(regexpr("\\-", fmt), "match.length")
    is.left <- ifelse(len > 0L, TRUE, FALSE)
    if (len > 1L)
      return()
    len <- attr(regexpr("\\.", fmt), "match.length")
    is.period <- ifelse(len > 0L, TRUE, FALSE)
    if (len > 1L | (is.period & code %in% c("d", "s")))
      return()
    len <- attr(regexpr("[[:space:]]", fmt), "match.length")
    is.space <- ifelse(len > 0L, TRUE, FALSE)
    if (len > 1L | (is.space & code == "s"))
      return()
    len <- attr(regexpr("\\+", fmt), "match.length")
    is.sign <- ifelse(len > 0L, TRUE, FALSE)
    if (len > 1L | (is.sign & code == "s"))
      return()

    is.pad <- FALSE
    if (nchar(fmt) > 2L) {
      fmt <- substr(fmt, 2L, nchar(fmt) - 1L)
      if (is.period) {
        fmt.split <- strsplit(fmt, "\\.")[[1L]]
        if (length(fmt.split) > 1L) {
          precision <- suppressWarnings(as.integer(fmt.split[2L]))
          if (is.na(precision))
            return()
        }
        fmt <- fmt.split[1L]
      }
      fmt <- sub("[[:space:]]", "", fmt)
      fmt <- sub("\\+", "", fmt)
      fmt <- sub("\\-", "", fmt)
      if (nchar(fmt) > 1L) {
        if (substr(fmt, 1L, 1L) == "0") {
          is.pad <- TRUE
          fmt <- substr(fmt, 2L, nchar(fmt))
        }
      }
      if (nchar(fmt) > 0L) {
        width <- suppressWarnings(as.integer(fmt))
        if (is.na(width))
          return()
      }
    }

    tclvalue(opt.var) <- 2L
    tclvalue(width.var) <- width
    tclvalue(precision.var) <- precision
    tclvalue(scientific.var) <- is.scientific
    tclvalue(left.var) <- is.left
    tclvalue(sign.var) <- is.sign
    tclvalue(space.var) <- is.space
    tclvalue(pad.var) <- is.pad
    BuildFormat()
    return()
  }

  # Add string to conversion format entry
  AddString <- function(txt) {
    if (as.logical(tcl(frame3.ent.1, "selection", "present")))
      tcl(frame3.ent.1, "delete", "sel.first", "sel.last")
    tkinsert(frame3.ent.1, "insert", txt)
    UpdateSample()
    tkfocus(frame3.ent.1)
  }

  # Toggle GUI state based on custom check box

  ToggleState <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))

    opt <- as.integer(tclvalue(opt.var))

    s <- ifelse(opt == 2L, "normal", "readonly")
    tkconfigure(frame2.ent.1.2, state=s)
    tkconfigure(frame2.ent.1.4, state=s)

    s <- ifelse(opt == 2L, "normal",  "disabled")
    tkconfigure(frame2.lab.1.1, state=s)
    tkconfigure(frame2.lab.1.3, state=s)
    tkconfigure(frame2.chk.2.1, state=s)
    tkconfigure(frame2.chk.3.1, state=s)
    tkconfigure(frame2.chk.4.1, state=s)
    tkconfigure(frame2.chk.5.1, state=s)
    tkconfigure(frame2.chk.6.1, state=s)

    s <- ifelse(opt == 3L, "normal", "readonly")
    tkconfigure(frame3.ent.1, state=s)

    s <- ifelse(opt == 3L, "normal", "disabled")
    tkconfigure(frame3a.but.01, state=s)
    tkconfigure(frame3a.but.02, state=s)
    tkconfigure(frame3a.but.03, state=s)
    tkconfigure(frame3a.but.04, state=s)
    tkconfigure(frame3a.but.05, state=s)
    tkconfigure(frame3a.but.06, state=s)
    tkconfigure(frame3a.but.07, state=s)
    tkconfigure(frame3a.but.08, state=s)
    tkconfigure(frame3a.but.09, state=s)
    tkconfigure(frame3a.but.10, state=s)

    if (opt == 1L) {
      UpdateSample()
    } else if (opt == 2L) {
      BuildFormat()
      if (inherits(sample, "numeric"))
        tkfocus(frame2.ent.1.4)
      else
        tkfocus(frame2.ent.1.2)
    } else if (opt == 3L) {
      UpdateSample()
      tkfocus(frame3.ent.1)
    }
    return()
  }

  # Update sample value
  UpdateSample <- function() {
    opt <- as.integer(tclvalue(opt.var))
    fmt <- as.character(tclvalue(fmt.var))
    if (fmt == "" || opt == 1L) {
      tclvalue(sample.var) <- format(sample)
    } else {
      ans <- try(sprintf(fmt, sample), silent=TRUE)
      tclvalue(sample.var) <- ifelse(inherits(ans, "try-error"), "", ans)
    }
  }

  # Build format

  BuildFormat <- function() {
    tclvalue(width.var) <- CheckEntry("integer", tclvalue(width.var))
    width <- as.character(tclvalue(width.var))

    is.left <- as.logical(as.integer(tclvalue(left.var)))
    left <- ifelse(is.left, "-", "")

    if (is.numeric(sample)) {
      is.space <- as.logical(as.integer(tclvalue(space.var)))
      is.pad   <- as.logical(as.integer(tclvalue(pad.var)))
      is.sign  <- as.logical(as.integer(tclvalue(sign.var)))

      space <- ifelse(is.space, " ", "")
      pad <- ifelse((is.pad && width != ""), "0", "")
      sign <- ifelse(is.sign, "+", "")
      period <- ""
      precision <- ""

      if (is.integer(sample)) {
        letter <- "d"
      } else {
        tclvalue(precision.var) <- CheckEntry("integer",
                                              tclvalue(precision.var))
        precision <- as.character(tclvalue(precision.var))
        if (!(width == "" & precision == ""))
          period <- "."

        is.scientific <- as.logical(as.integer(tclvalue(scientific.var)))
        letter <- ifelse(is.scientific, "e", "f")
      }
      fmt <- paste0("%", space, pad, left, sign, width, period, precision,
                    letter)
    } else {
      fmt <- paste0("%", left, width, "s")
    }
    tclvalue(fmt.var) <- fmt
    UpdateSample()
  }

  # Copy format to clipboard
  CopyFormat <- function() {
    txt <- as.character(tclvalue(fmt.var))
    cat(txt, file="clipboard")
  }

  # Paste format from clipboard
  PasteFormat <- function() {
    opt <- as.integer(tclvalue(opt.var))
    if (opt != 3L)
      return()
    cb <- try(scan(file="clipboard", what="character", sep="\n", quiet=TRUE),
              silent=TRUE)
    if (inherits(cb, "try-error"))
      return()
    tclvalue(fmt.var) <- cb
    UpdateSample()
    tkfocus(frame3.ent.1)
  }


  ## Main program

  ans <- try(as.character(sample), silent=TRUE)
  if (inherits(ans, "try-error") || length(ans) != 1)
    stop("class of sample object is not valid")
  if (!is.character(fmt))
    stop("format argument must be of class character")

  new.fmt <- NULL

  # Assign variables linked to Tk widgets

  opt.var        <- tclVar()
  sample.var     <- tclVar()
  fmt.var        <- tclVar(fmt)
  width.var      <- tclVar()
  precision.var  <- tclVar()
  scientific.var <- tclVar(0)
  left.var       <- tclVar(0)
  sign.var       <- tclVar(0)
  space.var      <- tclVar(0)
  pad.var        <- tclVar(0)
  tt.done.var    <- tclVar(0)

  TranslateFormat()
  UpdateSample()

  # Open GUI

  tclServiceMode(FALSE)
  tt <- tktoplevel()
  tkwm.resizable(tt, 1, 0)

  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }
  tktitle(tt) <- "Format"

  # Frame 0, ok and cancel buttons

  frame0 <- ttkframe(tt, relief="flat")
  frame0.but.2 <- ttkbutton(frame0, width=12, text="OK", command=SaveFormat)
  frame0.but.3 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("Format", package="RSurvey"))
                            })
  tkgrid("x", frame0.but.2, frame0.but.3, frame0.but.4,
         sticky="se", pady=10, padx=c(4, 0))
  tkgrid.columnconfigure(frame0, 0, weight=1)
  tkgrid.configure(frame0.but.4, padx=c(4, 10))
  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1, format using

  frame1 <- ttklabelframe(tt, relief="flat", borderwidth=5, padding=5,
                          text="Format using")

  frame1.rbt.1.1 <- ttkradiobutton(frame1, variable=opt.var, value=1L,
                                   text="Common format",
                                   command=ToggleState)
  frame1.rbt.1.2 <- ttkradiobutton(frame1, variable=opt.var, value=2L,
                                   text="C-style options",
                                   command=ToggleState)
  frame1.rbt.1.3 <- ttkradiobutton(frame1, variable=opt.var, value=3L,
                                   text="C-style string",
                                   command=ToggleState)

  tkgrid(frame1.rbt.1.1, frame1.rbt.1.2, frame1.rbt.1.3)

  tkgrid.configure(frame1.rbt.1.2, padx=10)

  tkpack(frame1, fill="x", padx=10, pady=10)

# Frame 2

  frame2 <- ttkframe(tt, relief="flat", borderwidth=0, padding=0)

  frame2.lab.1.1 <- ttklabel(frame2, text="Field width")
  frame2.lab.1.3 <- ttklabel(frame2, text="Precision")

  frame2.ent.1.2 <- ttkentry(frame2, textvariable=width.var, width=15)
  frame2.ent.1.4 <- ttkentry(frame2, textvariable=precision.var, width=15)

  txt <- "Use scientific notation"
  frame2.chk.2.1 <- ttkcheckbutton(frame2, text=txt, variable=scientific.var,
                                   command=BuildFormat)
  txt <- "Left adjustment of converted argument in its field"
  frame2.chk.3.1 <- ttkcheckbutton(frame2, text=txt, variable=left.var,
                                   command=BuildFormat)
  txt <- "Always print number with sign (\u002b/\u2212)"
  frame2.chk.4.1 <- ttkcheckbutton(frame2, text=txt, variable=sign.var,
                                   command=BuildFormat)
  txt <- "Prefix a space if the first character is not a sign"
  frame2.chk.5.1 <- ttkcheckbutton(frame2, text=txt, variable=space.var,
                                   command=BuildFormat)
  txt <- "Pad to the field width with leading zeros"
  frame2.chk.6.1 <- ttkcheckbutton(frame2, text=txt, variable=pad.var,
                                   command=BuildFormat)

  if (is.numeric(sample) && !is.integer(sample)) {
    tkgrid(frame2.lab.1.1, frame2.ent.1.2, frame2.lab.1.3, frame2.ent.1.4)
    tkgrid(frame2.chk.2.1, columnspan=4, sticky="w", padx=c(10, 0))
    tkgrid.configure(frame2.lab.1.3, padx=c(10, 2))
  } else {
    tkgrid(frame2.lab.1.1, frame2.ent.1.2, "x", "x", sticky="w")
    tkgrid.columnconfigure(frame2, 2, weight=1)
    tkgrid(frame2.chk.3.1, columnspan=3, sticky="w", padx=c(10, 0))
  }
  tkgrid.configure(frame2.lab.1.1, padx=c(10, 2))

  if (is.numeric(sample)) {
    tkgrid(frame2.chk.3.1, columnspan=4, sticky="w", padx=c(10, 0))
    tkgrid(frame2.chk.4.1, columnspan=4, sticky="w", padx=c(10, 0))
    tkgrid(frame2.chk.5.1, columnspan=4, sticky="w", padx=c(10, 0))
    tkgrid(frame2.chk.6.1, columnspan=4, sticky="w", padx=c(10, 0))
  }

  tkpack(frame2, padx=10, pady=0, anchor="w")

  # Frame 3, conversion specification format

  frame3 <- ttklabelframe(tt, relief="flat", borderwidth=5, padding=5,
                          text="Conversion specification string")

  frame3.ent.1 <- ttkentry(frame3, textvariable=fmt.var, width=30)

  frame3a <- ttkframe(frame3, relief="flat", borderwidth=0, padding=0)

  frame3a.but.01 <- ttkbutton(frame3a, width=2, text="%",
                              command=function() AddString("%"))
  frame3a.but.02 <- ttkbutton(frame3a, width=2, text="\u002b",
                              command=function() AddString("+"))
  frame3a.but.03 <- ttkbutton(frame3a, width=2, text="\u2212",
                              command=function() AddString("-"))
  frame3a.but.04 <- ttkbutton(frame3a, width=2, text=" ",
                              command=function() AddString(" "))
  frame3a.but.05 <- ttkbutton(frame3a, width=2, text="0",
                              command=function() AddString("0"))
  frame3a.but.06 <- ttkbutton(frame3a, width=2, text=".",
                              command=function() AddString("."))
  frame3a.but.07 <- ttkbutton(frame3a, width=2, text="f",
                              command=function() AddString("f"))
  frame3a.but.08 <- ttkbutton(frame3a, width=2, text="e",
                              command=function() AddString("e"))
  frame3a.but.09 <- ttkbutton(frame3a, width=2, text="d",
                              command=function() AddString("d"))
  frame3a.but.10 <- ttkbutton(frame3a, width=2, text="s",
                              command=function() AddString("s"))

  frame3a.but.11 <- ttkbutton(frame3a, width=2, image=GetBitmapImage("copy"),
                             command=CopyFormat)
  frame3a.but.12 <- ttkbutton(frame3a, width=2, image=GetBitmapImage("paste"),
                             command=PasteFormat)

  if (is.numeric(sample)) {
    if (is.integer(sample))
      tkgrid(frame3a.but.01, frame3a.but.02, frame3a.but.03, frame3a.but.04,
             frame3a.but.05, frame3a.but.06, frame3a.but.07, frame3a.but.08,
             frame3a.but.09, frame3a.but.10, frame3a.but.11, frame3a.but.12,
             pady=c(2, 0), padx=c(0, 2))
    else
      tkgrid(frame3a.but.01, frame3a.but.02, frame3a.but.03, frame3a.but.04,
             frame3a.but.05, frame3a.but.06, frame3a.but.07, frame3a.but.08,
             frame3a.but.10, frame3a.but.11, frame3a.but.12,
             pady=c(2, 0), padx=c(0, 2))
  } else {
    if (is.logical(sample))
      tkgrid(frame3a.but.01, frame3a.but.02, frame3a.but.03, frame3a.but.04,
             frame3a.but.05, frame3a.but.06, frame3a.but.09, frame3a.but.10,
             frame3a.but.11, frame3a.but.12, pady=c(2, 0), padx=c(0, 2))
    else
      tkgrid(frame3a.but.01, frame3a.but.02, frame3a.but.03, frame3a.but.04,
             frame3a.but.05, frame3a.but.06, frame3a.but.10, frame3a.but.11,
             frame3a.but.12, pady=c(2, 0), padx=c(0, 2))
  }

  tkgrid(frame3.ent.1)
  tkgrid(frame3a, "x", pady=c(2, 0), sticky="w")
  tkgrid.configure(frame3a.but.10, padx=c(0, 10))

  tkgrid.configure(frame3.ent.1, sticky="we", columnspan=2, padx=c(0, 2))

  tkgrid.columnconfigure(frame3, 1, weight=1)

  tkpack(frame3, fill="x", padx=10, pady=10)

  # Frame 3, sample

  frame3 <- ttklabelframe(tt, relief="flat", borderwidth=5, padding=5,
                          text="Sample")
  frame3.ent <- ttkentry(frame3, textvariable=sample.var, width=30,
                         state="readonly", takefocus=FALSE)
  tkgrid(frame3.ent)
  tkgrid.configure(frame3.ent, sticky="we")
  tcl("grid", "anchor", frame3, "w")
  tkgrid.columnconfigure(frame3, 0, weight=1, minsize=13)
  tkpack(frame3, fill="x", expand=TRUE, padx=10)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame3.ent.1, "<KeyRelease>", UpdateSample)

  tkbind(frame2.ent.1.2, "<KeyRelease>", BuildFormat)
  tkbind(frame2.ent.1.4, "<KeyRelease>", BuildFormat)

  # GUI control

  ToggleState()

  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  return(new.fmt)
}
