# A GUI for defining a function in the R language with a focus on table data.

EditFunction <- function(cols, index=NULL, fun=NULL, value.length=NULL,
                         value.class=NULL, win.title="Edit Function",
                         parent=NULL) {

  ## Additional functions

  # Save function
  SaveFunction <- function() {
    txt <- as.character(tclvalue(tkget(frame2.txt.2.1, "1.0", "end-1c")))
    if (txt == "") {
      rtn <<- list(fun="")
    } else {
      fun <- txt
      pattern <- paste0("\"", ids, "\"")
      replacement <- paste0("DATA[[", seq_along(ids), "]]")
      for (i in seq_along(ids))
        fun <- gsub(pattern[i], replacement[i], fun, fixed=TRUE)
      fun <- paste0("function(DATA) {", fun, "}")

      fun <- try(parse(text=fun), silent=TRUE)
      if (inherits(fun, "try-error")) {
        msg <- "There's a problem with function syntax, try revising."
        tkmessageBox(icon="error", message=msg, detail=fun, title="Error",
                     type="ok", parent=tt)
        return()
      }
      obj <- EvalFunction(txt, cols)
      if (inherits(obj, "try-error")) {
        msg <- "Function results in error during evaluation, try revising."
        tkmessageBox(icon="error", message=msg, detail=obj, title="Error",
                     type="ok", parent=tt)
        return()
      }

      if (!is.null(value.length) && length(obj) != value.length) {
        msg <- paste0("Evaluated function must be of length ", value.length,
                      ", try revising.")
        dtl <- paste0("Resulting vector is currently of length ", length(obj),
                      ".")
        tkmessageBox(icon="error", message=msg, detail=dtl, title="Error",
                     type="ok", parent=tt)
        return()
      }

      if (!is.null(value.class) && !inherits(obj, value.class)) {
        msg <- paste0("Filter must result in a vector of class \"",
                      value.class,
                      "\". The evaluated function is a vector of class \"",
                      class(obj), "\", please revise.")
        tkmessageBox(icon="error", message=msg, title="Error", type="ok",
                     parent=tt)
        return()
      }

      rtn <<- list(fun=txt, class=class(obj)[1], sample=na.omit(obj)[1],
                   summary=summary(obj))
    }
    tclvalue(tt.done.var) <- 1
  }

  # Rebuild list box based on selected class type to show
  RebuildList <- function() {
    idx <- as.integer(tcl(frame1.box.3.1, "current"))
    if (idx > 0)
      show.ids <- ids[vapply(cols, function(i) classes[idx] %in% i$class, TRUE)]
    else
      show.ids <- ids

    tclvalue(variable.var) <- ""
    for (i in seq_along(show.ids))
      tcl("lappend", variable.var, show.ids[i])

    tkselection.clear(frame1.lst.2.1, 0, "end")

    tclvalue(value.var) <- ""

    tkconfigure(frame1.but.5.1, state="disabled")
    tkfocus(frame2.txt.2.1)
  }

  # Insert character string into text box
  InsertString <- function(txt, sel="<variable>") {
    tcl(frame2.txt.2.1, "edit", "separator")
    seltxt <- as.character(tktag.ranges(frame2.txt.2.1, "sel"))
    if (length(seltxt) > 1)
      tcl(frame2.txt.2.1, "delete", seltxt[1], seltxt[2])

    cur <- as.character(tkindex(frame2.txt.2.1, "insert"))
    cur <- as.integer(strsplit(cur, ".", fixed=TRUE)[[1]])
    cur.line <- cur[1]
    cur.char <- cur[2]

    tkinsert(frame2.txt.2.1, "insert", txt)
    tkfocus(frame2.txt.2.1)

    if (txt %in% c("()", "[]")) {
      cursor.insert <- paste(cur.line, cur.char + 1, sep=".")
      tkmark.set(frame2.txt.2.1, "insert", cursor.insert)
    } else {
      search.txt <- gregexpr(pattern=sel, txt)[[1]]
      if (search.txt[1] > 0) {
        match.idx <- search.txt[1]
        match.len <- attr(search.txt, "match.length")[1]
        tkfocus(frame2.txt.2.1)

        char <- c(match.idx, match.idx + match.len) + cur.char - 1
        sel0 <- paste(cur.line, char[1], sep=".")
        sel1 <- paste(cur.line, char[2], sep=".")
        tktag.add(frame2.txt.2.1, "sel", sel0, sel1)
        tkmark.set(frame2.txt.2.1, "insert", sel1)
      }
    }
  }

  # Insert variable into text box
  InsertVar <- function() {
    idx <- as.integer(tkcurselection(frame1.lst.2.1))
    if (length(idx) == 0)
      return()
    id <- as.character(tkget(frame1.lst.2.1, idx, idx))
    txt <- paste0("\"", id, "\"")
    InsertString(txt)
  }

  # Call date and time format editor
  CallFormatDateTime <- function(sample) {
    fmt <- FormatDateTime(sample=sample, parent=tt)
    tkfocus(frame2.txt.2.1)
    if(!is.null(fmt))
      InsertString(gsub("%OS[[:digit:]]+", "%OS", fmt))
  }

  # Text edit functions
  EditUndo <- function() {
    tkfocus(frame2.txt.2.1)
    try(tcl(frame2.txt.2.1, "edit", "undo"), silent=TRUE)
  }
  EditRedo <- function() {
    tkfocus(frame2.txt.2.1)
    try(tcl(frame2.txt.2.1, "edit", "redo"), silent=TRUE)
  }
  EditCut <- function() {
    tkfocus(frame2.txt.2.1)
    tcl("tk_textCut", frame2.txt.2.1)
  }
  EditCopy <- function() {
    tkfocus(frame2.txt.2.1)
    tcl("tk_textCopy", frame2.txt.2.1)
  }
  EditPaste <- function() {
    tkfocus(frame2.txt.2.1)
    tcl("tk_textPaste", frame2.txt.2.1)
  }
  EditSelectAll <- function() {
    tkfocus(frame2.txt.2.1)
    tktag.add(frame2.txt.2.1, "sel", "1.0", "end")
  }

  # Clear console
  ClearConsole <- function() {
    tcl(frame2.txt.2.1, "delete", "1.0", "end")
    tkfocus(frame2.txt.2.1)
  }

  # Show unique values
  ShowUniqueValues <- function() {
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"))
    idx <- as.integer(tkcurselection(frame1.lst.2.1))
    if (length(idx) == 0)
      return()
    id <- as.character(tkget(frame1.lst.2.1, idx, idx))
    idx <- which(vapply(cols, function(i) i$id, "") == id)

    var.fmt <- cols[[idx]]$format
    if (is.null(var.fmt))
      var.fmt <- ""

    var.vals <- unique(EvalFunction(cols[[idx]]$fun, cols))
    var.class <- cols[[idx]]$class

    n <- length(var.vals)
    if (n > 50000) {
      msg <- paste("There are", n, "unique values; this operation can be",
                   "computationally expensive. Would you like to continue?")
      ans <- as.character(tkmessageBox(icon="question", message=msg,
                                       title="Warning", type="yesno",
                                       parent=tt))
      if (ans == "no")
        return()
    }

    var.vals <- sort(var.vals, na.last=TRUE)
    if (var.fmt == "") {
      var.vals.txt <- format(var.vals)
    } else if ("POSIXt" %in% var.class) {
      var.vals.txt <- format(var.vals, format=var.fmt)
    } else {
      var.vals.txt <- try(sprintf(var.fmt, var.vals), silent=TRUE)
      if (inherits(var.vals.txt, "try-error"))
        var.vals.txt <- format(var.vals)
    }
    var.vals.txt <- gsub("^\\s+|\\s+$", "", var.vals.txt)

    tclvalue(value.var) <- ""
    for (i in seq_along(var.vals.txt))
      tcl("lappend", value.var, var.vals.txt[i])
    tkselection.clear(frame1.lst.4.1, 0, "end")
    tkconfigure(frame1.but.5.1, state="disabled")
    tkfocus(frame2.txt.2.1)
  }

  # Change variable selection
  ChangeVar <- function() {
    tclvalue(value.var) <- ""
    idx <- as.integer(tkcurselection(frame1.lst.2.1))
    if (length(idx) == 0)
      return()
    tkconfigure(frame1.but.5.1, state="normal")
  }

  # Insert value into text box
  InsertValue <- function() {
    idx <- as.integer(tkcurselection(frame1.lst.2.1))
    if (length(idx) == 0)
      return()
    id <- as.character(tkget(frame1.lst.2.1, idx, idx))
    idx <- which(vapply(cols, function(i) i$id, "") == id)
    var.fmt <- cols[[idx]]$format
    var.class <- cols[[idx]]$class
    idx <- as.integer(tkcurselection(frame1.lst.4.1))
    if (length(idx) == 0)
      return()
    val <- as.character(tkget(frame1.lst.4.1, idx, idx))
    if ("factor" %in% var.class && is.na(suppressWarnings(as.numeric(val))))
      var.class <- "character"
    if ("POSIXt" %in% var.class) {
      if (var.fmt == "")
        var.fmt <- "%Y-%m-%d %H:%M:%S"
      txt <- paste0("as.POSIXct(\"", val, "\", format = \"", var.fmt,
                    "\", tz = \"GMT\")")
    } else if ("Date" %in% var.class) {
      if (var.fmt == "")
        var.fmt <- "%Y-%m-%d"
      txt <- paste0("as.Date(\"", val, "\", format = \"", var.fmt, "\")")
    } else if ("integer" %in% var.class &&
               !val %in% c("NA", "NaN", "Inf", "-Inf")) {
      txt <- paste0(val, "L")
    } else if ("character" %in% var.class &&
               !val %in% c("NA", "NaN", "Inf", "-Inf")) {
      txt <- paste0("\"", val, "\"")
    } else {
      txt <- val
    }
    InsertString(txt)
  }

  # Insert trigonometric function
  InsertTrigFunction <- function(fun) {
    is.inverse <- as.logical(as.integer(tclvalue(inverse.var)))
    is.hyperbolic <- as.logical(as.integer(tclvalue(hyperbolic.var)))
    angles <- as.character(tclvalue(angles.var))
    if (is.inverse)
      fun <- paste0("a", fun)
    if (is.hyperbolic)
      fun <- paste0(fun, "h")
    if (angles == "deg")
      InsertString(paste0(fun, "(<variable> * pi / 180)"))
    else
      InsertString(paste0(fun, "(<variable>)"))
  }

  ## Main program

  if (is.null(index)) {
    old.fun <- fun
  } else {
    old.fun <- cols[[as.integer(index)]]$fun
  }
  rtn <- NULL
  ids <- vapply(cols, function(i) i$id, "")

  # Remove variable being defined
  if (!is.null(index)) {
    edit.fun.id <- ids[index]
    ids <- ids[-index]
  }

  # Class types
  classes <- NULL
  for (i in seq_along(cols)) {
    if (!i %in% index)
      classes <- c(classes, cols[[i]]$class)
  }
  classes <- suppressWarnings(sort(unique(classes)))

  # Required vector length
  if (!is.null(value.length))
    value.length <- as.integer(value.length)

  # Assign variables linked to Tk widgets
  variable.var <- tclVar()
  for (i in seq_along(ids))
    tcl("lappend", variable.var, ids[i])  # must be unique
  value.var <- tclVar()
  inverse.var <- tclVar(0)
  hyperbolic.var <- tclVar(0)
  angles.var <- tclVar("rad")
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

  # Top menu
  top.menu <- tkmenu(tt, tearoff=0)

  # Project menu

  menu.edit <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Edit", menu=menu.edit, underline=0)
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

  menu.convert <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Convert", menu=menu.convert, underline=0)

  menu.convert.char <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.char, "command", label="Factor",
        command=function() InsertString("as.factor(<variable>)"))
  tkadd(menu.convert.char, "command", label="Numeric",
        command=function() InsertString("as.numeric(<variable>)"))
  tkadd(menu.convert.char, "command", label="Integer",
        command=function() InsertString("as.integer(<variable>)"))
  tkadd(menu.convert.char, "command", label="Logical",
        command=function() InsertString("as.logical(<variable>)"))
  tkadd(menu.convert.char, "command", label="POSIXct",
        command=function() InsertString("as.POSIXct(strptime(<variable>, format = \"<format>\", tz = \"GMT\"))"))
  tkadd(menu.convert.char, "command", label="Date",
        command=function() InsertString("as.Date(<variable>, format = \"<format>\")"))
  tkadd(menu.convert, "cascade", label="Character to", menu=menu.convert.char)

  menu.convert.factor <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.factor, "command", label="Character",
        command=function() InsertString("as.character(<variable>)"))
  tkadd(menu.convert.factor, "command", label="Integer",
        command=function() InsertString("as.integer(<variable>)"))
  tkadd(menu.convert, "cascade", label="Factor to", menu=menu.convert.factor)

  menu.convert.num <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.num, "command", label="POSIXct",
        command=function() InsertString("as.POSIXct(<variable>, origin = \"1970-01-01 00:00:00.00\", tz = \"GMT\")"))
  tkadd(menu.convert, "cascade", label="Numeric to", menu=menu.convert.num)

  menu.convert.int <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.int, "command", label="POSIXct",
        command=function() InsertString("as.POSIXct(<variable>, origin = \"1970-01-01 00:00:00\", tz = \"GMT\")"))
  tkadd(menu.convert.int, "command", label="Date",
        command=function() InsertString("as.Date(<variable>, origin = \"1899-12-30\")"))
  tkadd(menu.convert, "cascade", label="Integer to", menu=menu.convert.int)

  menu.convert.log <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.log, "command", label="Integer",
        command=function() InsertString("as.integer(<variable>)"))
  tkadd(menu.convert, "cascade", label="Logical to", menu=menu.convert.log)

  menu.convert.posix <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.posix, "command", label="Numeric",
        command=function() InsertString("as.numeric(<variable>)"))
  tkadd(menu.convert.posix, "command", label="Date",
        command=function() InsertString("as.Date(<variable>, tz = \"UTC\")"))
  tkadd(menu.convert, "cascade", label="POSIXct to", menu=menu.convert.posix)

  menu.convert.date <- tkmenu(tt, tearoff=0)
  tkadd(menu.convert.date, "command", label="Integer",
        command=function() InsertString("as.integer(<variable>)"))
  tkadd(menu.convert.date, "command", label="POSIXct",
        command=function() InsertString("as.POSIXct(<variable>, tz = \"GMT\")"))
  tkadd(menu.convert, "cascade", label="Date to", menu=menu.convert.date)

  menu.math <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Math", menu=menu.math, underline=0)
  tkadd(menu.math, "command", label="Square root",
        command=function() InsertString("sqrt(<variable>)"))
  tkadd(menu.math, "command", label="Absolute value",
        command=function() InsertString("abs(<variable>)"))
  tkadd(menu.math, "separator")
  tkadd(menu.math, "command", label="Floor",
        command=function() InsertString("floor(<variable>)"))
  tkadd(menu.math, "command", label="Ceiling",
        command=function() InsertString("ceiling(<variable>)"))
  tkadd(menu.math, "command", label="Truncation",
        command=function() InsertString("trunc(<variable>)"))
  menu.math.round <- tkmenu(tt, tearoff=0)
  tkadd(menu.math.round, "command", label="Decimal places",
        command=function() InsertString("round(<variable>, digits = 0)"))
  tkadd(menu.math.round, "command", label="Significant digits",
        command=function() InsertString("signif(<variable>, digits = 6)"))
  tkadd(menu.math, "cascade", label="Round to", menu=menu.math.round)
  tkadd(menu.math, "separator")
  tkadd(menu.math, "command", label="Exponential",
        command=function() InsertString("exp(<variable>)"))
  menu.math.log <- tkmenu(tt, tearoff=0)
  tkadd(menu.math.log, "command", label="Common (base 10)",
        command=function() InsertString("log10(<variable>)"))
  tkadd(menu.math.log, "command", label="Natural (base \u0065)",
        command=function() InsertString("log(<variable>, base = exp(1))"))
  tkadd(menu.math.log, "command", label="Binary (base 2)",
        command=function() InsertString("log2(<variable>)"))
  tkadd(menu.math, "cascade", label="Logarithm", menu=menu.math.log)
  tkadd(menu.math, "separator")
  tkadd(menu.math, "command", label="Sine",
        command=function() InsertTrigFunction("sin"))
  tkadd(menu.math, "command", label="Cosine",
        command=function() InsertTrigFunction("cos"))
  tkadd(menu.math, "command", label="Tangent",
        command=function() InsertTrigFunction("tan"))
  tkadd(menu.math, "separator")
  menu.math.cum <- tkmenu(tt, tearoff=0)
  tkadd(menu.math.cum, "command", label="Sum",
        command=function() InsertString("cumsum(<variable>)"))
  tkadd(menu.math.cum, "command", label="Product",
        command=function() InsertString("cumprod(<variable>)"))
  tkadd(menu.math, "cascade", label="Cumulative", menu=menu.math.cum)

  tkadd(menu.math, "separator")
  tkadd(menu.math, "command", label="Sum",
        command=function() InsertString("sum(<variable>, na.rm = TRUE)"))
  tkadd(menu.math, "command", label="Product",
        command=function() InsertString("prod(<variable>, na.rm = TRUE)"))

  menu.stats <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Stats", menu=menu.stats, underline=0)
  tkadd(menu.stats, "command", label="Minimum",
        command=function() InsertString("min(<variable>, na.rm = TRUE)"))
  tkadd(menu.stats, "command", label="Maximum",
        command=function() InsertString("max(<variable>, na.rm = TRUE)"))
  tkadd(menu.stats, "separator")
  tkadd(menu.stats, "command", label="Mean",
        command=function() InsertString("mean(<variable>, na.rm = TRUE)"))
  tkadd(menu.stats, "command", label="Median",
        command=function() InsertString("median(<variable>, na.rm = TRUE)"))
  tkadd(menu.stats, "command", label="Standard deviation",
        command=function() InsertString("sd(<variable>, na.rm = TRUE)"))

  tkadd(menu.stats, "separator")
  tkadd(menu.stats, "command", label="Set seed",
        command=function() InsertString("set.seed(124)"))
  menu.stats.ran <- tkmenu(tt, tearoff=0)
  nobs <- ifelse(is.null(value.length), "<integer>", value.length)
  tkadd(menu.stats.ran, "command", label="Normal distribution",
        command=function() InsertString(paste0("rnorm(n = ", nobs, ", mean = 0, sd = 1)")))
  tkadd(menu.stats.ran, "command", label="Uniform distribution",
        command=function() InsertString(paste0("runif(n = ", nobs, ", min = 0, max = 1)")))
  tkadd(menu.stats, "cascade", label="Random samples from a",
        menu=menu.stats.ran)

  menu.operator <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Operator", menu=menu.operator, underline=0)
  tkadd(menu.operator, "command", label="And",
        command=function() InsertString(" & "))
  tkadd(menu.operator, "command", label="Or",
        command=function() InsertString(" | "))
  tkadd(menu.operator, "command", label="Not",
        command=function() InsertString("!"))
  tkadd(menu.operator, "separator")
  tkadd(menu.operator, "command", label="In",
        command=function() InsertString(" %in% "))
  tkadd(menu.operator, "command", label="Match",
        command=function() {
          InsertString("match(<variable>, <values>, nomatch = NA)")
        })
  tkadd(menu.operator, "separator")
  tkadd(menu.operator, "command", label="Exponentiation",
        command=function() InsertString("<variable>^<power>"))

  menu.const <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Constant", menu=menu.const, underline=0)
  tkadd(menu.const, "command", label="\u03C0",
        command=function() InsertString("pi"))
  tkadd(menu.const, "command", label="\u0065",
        command=function() InsertString("exp(1)"))
  tkadd(menu.const, "separator")
  tkadd(menu.const, "command", label="True",
        command=function() InsertString("TRUE"))
  tkadd(menu.const, "command", label="False",
        command=function() InsertString("FALSE"))
  tkadd(menu.const, "separator")
  tkadd(menu.const, "command", label="Not available",
        command=function() InsertString("NA"))
  tkadd(menu.const, "command", label="Not a number",
        command=function() InsertString("NaN"))
  tkadd(menu.const, "separator")
  tkadd(menu.const, "command", label="Positive infinity",
        command=function() InsertString("Inf"))
  tkadd(menu.const, "command", label="Negative infinity",
        command=function() InsertString("-Inf"))
  tkadd(menu.const, "separator")
  menu.const.is <- tkmenu(tt, tearoff=0)
  tkadd(menu.const.is, "command", label="Not available",
        command=function() InsertString("is.na(<variable>)"))
  tkadd(menu.const.is, "command", label="Not a number",
        command=function() InsertString("is.nan(<variable>)"))
  tkadd(menu.const.is, "separator")
  tkadd(menu.const.is, "command", label="Finite (not infinite and not missing)",
        command=function() InsertString("is.finite(<variable>)"))
  tkadd(menu.const.is, "command", label="Infinite",
        command=function() InsertString("is.infinite(<variable>)"))
  tkadd(menu.const, "cascade", label="Which elements are ", menu=menu.const.is)
  tkadd(menu.const, "separator")
  tkadd(menu.const, "command", label="Are all values true",
        command=function() InsertString("all(<variable>, na.rm = FALSE)"))
  tkadd(menu.const, "command", label="Are any values true",
        command=function() InsertString("any(<variable>, na.rm = FALSE)"))

  menu.string <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="String", menu=menu.string, underline=0)
  tkadd(menu.string, "command", label="Concatenate",
        command=function() {
          InsertString("paste(<variable>, <variable>, sep = \" \")")
        })
  tkadd(menu.string, "command", label="Extract substring",
        command=function() {
          InsertString(paste("substr(<variable>, start = 1, stop = 2)"))
        })
  tkadd(menu.string, "command", label="Number of characters",
        command=function() {
          InsertString("nchar(<variable>)")
        })
  tkadd(menu.string, "command", label="Which elements are non-empty strings",
        command=function() {
          InsertString("nzchar(<variable>)")
        })

  menu.tools <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Tools", menu=menu.tools, underline=0)
  menu.tools.time <- tkmenu(tt, tearoff=0)
  tkadd(menu.tools.time, "command", label="POSIXct\u2026",
        command=function() CallFormatDateTime(Sys.time()))
  tkadd(menu.tools.time, "command", label="Date\u2026",
        command=function() CallFormatDateTime(Sys.Date()))
  tkadd(menu.tools, "cascade", label="Build format for", menu=menu.tools.time)

  # Finalize top menu
  tkconfigure(tt, menu=top.menu)

  # Frame 0, ok and cancel buttons, and size grip

  frame0 <- tkframe(tt, relief="flat")

  frame0a <- tkframe(frame0, relief="flat")
  frame0a.chk.1 <- ttkcheckbutton(frame0a, text="Inverse", variable=inverse.var)
  frame0a.chk.2 <- ttkcheckbutton(frame0a, text="Hyperbolic",
                                  variable=hyperbolic.var)
  frame0a.rad.3 <- ttkradiobutton(frame0a, variable=angles.var, value="rad",
                                  text="Radians")
  frame0a.rad.4 <- ttkradiobutton(frame0a, variable=angles.var, value="deg",
                                  text="Degrees")
  tkgrid(frame0a.chk.1, frame0a.chk.2, frame0a.rad.3, frame0a.rad.4)
  tkgrid.configure(frame0a.chk.2, padx=c(2, 12))
  tkgrid.configure(frame0a.rad.3, padx=c(0, 2))

  frame0.but.3 <- ttkbutton(frame0, width=12, text="OK",
                            command=SaveFunction)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.5 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("EditFunction", package="RSurvey"))
                            })
  frame0.grp.6 <- ttksizegrip(frame0)

  tkgrid(frame0a, "x", frame0.but.3, frame0.but.4, frame0.but.5,
         frame0.grp.6)

  tkgrid.columnconfigure(frame0, 1, weight=1)

  tkgrid.configure(frame0a, padx=10, pady=c(0, 10), sticky="sw")
  tkgrid.configure(frame0.but.3, frame0.but.4, frame0.but.5,
                   padx=c(0, 4), pady=c(15, 10))
  tkgrid.configure(frame0.but.5, columnspan=2, padx=c(0, 10))
  tkgrid.configure(frame0.grp.6, sticky="se")

  tkraise(frame0.but.5, frame0.grp.6)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Paned window
  pw <- ttkpanedwindow(tt, orient="horizontal")

  # Frame 1

  frame1 <- tkframe(pw, relief="flat")

  txt <- "Double click to insert variable"
  frame1.lab.1.1 <- ttklabel(frame1, text=txt, foreground="#141414")
  frame1.lst.2.1 <- tklistbox(frame1, selectmode="browse", activestyle="none",
                              relief="flat", borderwidth=5, width=25, height=8,
                              exportselection=FALSE, listvariable=variable.var,
                              highlightthickness=0)
  frame1.ysc.2.2 <- ttkscrollbar(frame1, orient="vertical")
  box.vals <- "{Show all classes}"
  if (length(classes) > 0)
    box.vals <- c("Show all classes", classes)
  frame1.box.3.1 <- ttkcombobox(frame1, state="readonly", value=box.vals)
  tkconfigure(frame1.lst.2.1, background="white",
              yscrollcommand=paste(.Tk.ID(frame1.ysc.2.2), "set"))
  tkconfigure(frame1.ysc.2.2, command=paste(.Tk.ID(frame1.lst.2.1), "yview"))
  tcl(frame1.box.3.1, "current", 0)

  frame1.lst.4.1 <- tklistbox(frame1, selectmode="browse", activestyle="none",
                              relief="flat", borderwidth=5, width=25, height=5,
                              exportselection=FALSE, listvariable=value.var,
                              highlightthickness=0)
  frame1.ysc.4.2 <- ttkscrollbar(frame1, orient="vertical")
  frame1.but.5.1 <- ttkbutton(frame1, width=15, text="Unique Values",
                              command=ShowUniqueValues)
  tkconfigure(frame1.lst.4.1, background="white",
              yscrollcommand=paste(.Tk.ID(frame1.ysc.4.2), "set"))
  tkconfigure(frame1.ysc.4.2, command=paste(.Tk.ID(frame1.lst.4.1), "yview"))
  tkconfigure(frame1.but.5.1, state="disabled")

  tkgrid(frame1.lab.1.1, "x")
  tkgrid(frame1.lst.2.1, frame1.ysc.2.2)
  tkgrid(frame1.box.3.1, "x")
  tkgrid(frame1.lst.4.1, frame1.ysc.4.2)
  tkgrid(frame1.but.5.1, "x")

  tkgrid.configure(frame1.lab.1.1, padx=c(10, 0), pady=c(10, 0), sticky="w")
  tkgrid.configure(frame1.lst.2.1, padx=c(10, 0), pady=c(2, 0),  sticky="nsew")
  tkgrid.configure(frame1.ysc.2.2, padx=c(0, 0),  pady=c(2, 0),  sticky="ns")
  tkgrid.configure(frame1.box.3.1, padx=c(10, 0), pady=c(4, 0),  sticky="we")
  tkgrid.configure(frame1.lst.4.1, padx=c(10, 0), pady=c(15, 0), sticky="nsew")
  tkgrid.configure(frame1.ysc.4.2, padx=c(0, 0),  pady=c(15, 0), sticky="ns")
  tkgrid.configure(frame1.but.5.1, padx=c(10, 0), pady=c(4, 0))

  tkgrid.rowconfigure(frame1, 1, weight=1)
  tkgrid.rowconfigure(frame1, 3, weight=1)
  tkgrid.columnconfigure(frame1, 0, weight=1, minsize=20)

  # Frame 2

  frame2 <- tkframe(pw, relief="flat")

  txt <- "Define function"
  if (!is.null(index) && edit.fun.id != "")
    txt <- paste0(txt, " for \"", edit.fun.id, "\"")
  frame2.lab.1.1 <- ttklabel(frame2, text=txt, foreground="#141414")

  if (is.null(value.length)) {
    txt <- ""
  } else {
    txt <- paste("(resulting vector must be of length", format(value.length,
                                                               big.mark=","))
    if (!is.null(value.class))
      txt <- paste(txt, "and class", value.class)
    txt <- paste0(txt, ")")
  }
  frame2.lab.1.2 <- ttklabel(frame2, text=txt, foreground="#A40802")

  frame2.txt.2.1 <- tktext(frame2, bg="white", font="TkFixedFont",
                           padx=2, pady=2, width=80, height=12, undo=1,
                           autoseparators=1, wrap="none", foreground="black",
                           relief="flat",
                           yscrollcommand=function(...)
                                            tkset(frame2.ysc.2.2, ...))

  frame2.ysc.2.2 <- ttkscrollbar(frame2, orient="vertical")
  tkconfigure(frame2.ysc.2.2, command=paste(.Tk.ID(frame2.txt.2.1), "yview"))

  frame2a <- tkframe(frame2, relief="flat")
  frame2a.but.01 <- ttkbutton(frame2a, width=3, text="\u002b",
                              command=function() InsertString(" + "))
  frame2a.but.02 <- ttkbutton(frame2a, width=3, text="\u2212",
                              command=function() InsertString(" - "))
  frame2a.but.03 <- ttkbutton(frame2a, width=3, text="\u00d7",
                              command=function() InsertString(" * "))
  frame2a.but.04 <- ttkbutton(frame2a, width=3, text="\u00f7",
                              command=function() InsertString(" / "))
  frame2a.but.05 <- ttkbutton(frame2a, width=3, text=">",
                              command=function() InsertString(" > "))
  frame2a.but.06 <- ttkbutton(frame2a, width=3, text="<",
                              command=function() InsertString(" < "))
  frame2a.but.07 <- ttkbutton(frame2a, width=3, text="\u2265",
                              command=function() InsertString(" >= "))
  frame2a.but.08 <- ttkbutton(frame2a, width=3, text="\u2264",
                              command=function() InsertString(" <= "))
  frame2a.but.09 <- ttkbutton(frame2a, width=3, text="=",
                              command=function() InsertString(" == "))
  frame2a.but.10 <- ttkbutton(frame2a, width=3, text="\u2260",
                              command=function() InsertString(" != "))
  frame2a.but.11 <- ttkbutton(frame2a, width=3, text="in",
                              command=function() InsertString(" %in% "))
  frame2a.but.12 <- ttkbutton(frame2a, width=3, text="( )",
                              command=function() InsertString("()"))
  frame2a.but.13 <- ttkbutton(frame2a, width=3, text="[ ]",
                              command=function() InsertString("[]"))
  frame2a.but.14 <- ttkbutton(frame2a, width=3, text="\u25C0\u2212",
                              command=function() InsertString(" <- "))
  frame2a.but.15 <- ttkbutton(frame2a, width=3, text="\u0078\u207F",
                              command=function() InsertString("^"))

  tkgrid(frame2a.but.01, frame2a.but.02, frame2a.but.03, frame2a.but.04,
         frame2a.but.05, frame2a.but.06, frame2a.but.07, frame2a.but.08,
         frame2a.but.09, frame2a.but.10, frame2a.but.11, frame2a.but.12,
         frame2a.but.13, frame2a.but.14, frame2a.but.15, pady=c(4, 0))

  tkgrid.configure(frame2a.but.01, frame2a.but.02, frame2a.but.03,
                   frame2a.but.04, frame2a.but.06, frame2a.but.07,
                   frame2a.but.08, frame2a.but.09, frame2a.but.10,
                   frame2a.but.11, frame2a.but.13, frame2a.but.14,
                   frame2a.but.15, padx=c(2, 0))
  tkgrid.configure(frame2a.but.05, frame2a.but.12, padx=c(12, 0))

  tkgrid(frame2.lab.1.1, frame2.lab.1.2, "x")
  tkgrid(frame2.txt.2.1, "x", frame2.ysc.2.2)
  tkgrid(frame2a, "x", "x")

  tkgrid.configure(frame2.lab.1.1, frame2.lab.1.2, padx=c(2, 0), pady=c(10, 0),
                   sticky="w")
  tkgrid.configure(frame2.txt.2.1, padx=c(2, 0),  pady=c(2, 0), columnspan=2,
                   sticky="nsew")
  tkgrid.configure(frame2.ysc.2.2, padx=c(0, 10), pady=c(2, 0), columnspan=2,
                   sticky="ns")
  tkgrid.configure(frame2a, columnspan=2, sticky="we")

  tkgrid.rowconfigure(frame2, 1, weight=1)
  tkgrid.columnconfigure(frame2, 1, weight=1, minsize=20)

  if (!is.null(old.fun) && old.fun != "")
    tkinsert(frame2.txt.2.1, "end", old.fun)

  tcl(frame2.txt.2.1, "edit", "reset")

  tkmark.set(frame2.txt.2.1, "insert", "end")

  # Pack frames into paned window
  tkadd(pw, frame1, weight=0)
  tkadd(pw, frame2, weight=1)
  tkpack(pw, fill="both", expand="yes")

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame1.lst.2.1, "<<ListboxSelect>>", ChangeVar)
  tkbind(frame1.lst.2.1, "<Double-ButtonRelease-1>", InsertVar)
  tkbind(frame1.lst.4.1, "<Double-ButtonRelease-1>", InsertValue)
  tkbind(frame1.box.3.1, "<<ComboboxSelected>>", RebuildList)

  tkbind("Text", "<Control-z>", EditUndo)
  tkbind("Text", "<Control-y>", EditRedo)
  tkbind("Text", "<Control-v>", EditPaste)
  tkbind("Text", "<Control-a>", EditSelectAll)
  tkbind("Text", "<Control-l>", ClearConsole)

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
