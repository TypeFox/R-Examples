showData <- function (dataframe, colname.bgcolor = "grey50", rowname.bgcolor = "grey50", 
          body.bgcolor = "white", colname.textcolor = "white", rowname.textcolor = "white", 
          body.textcolor = "black", font = "Courier 12", maxheight = 30, 
          maxwidth = 80, title = NULL, rowname.bar = "left", colname.bar = "top", 
          rownumbers = FALSE, placement = "-20-40", suppress.X11.warnings = TRUE) 
{
  object.name <- deparse(substitute(dataframe))
  if (!is.data.frame(dataframe)) {
    temp <- try(dataframe <- as.data.frame(dataframe), silent = FALSE)
    if (inherits(temp, "try-error")) {
      stop(paste(object.name, "cannot be coerced to a data frame"))
    }
    object.name <- paste("as.data.frame(", object.name, ")", 
                         sep = "")
  }
  if (is.numeric(rownumbers) && length(rownumbers) != nrow(dataframe)) {
    stop("rownumbers argument must be TRUE, FALSE or have length nrow(dataframe)")
  }
  if (!capabilities("tcltk")) 
    stop("tcltk capability missing")
  requireNamespace("tcltk", quietly = TRUE)
  oldwidth <- unlist(options("width"))
  options(width = 10000)
  conn <- file()
  sink(conn)
  print(dataframe)
  sink()
  zz <- scan(conn, sep = "\n", what = character(0), quiet = TRUE)
  close(conn)
  if (length(zz) > 1 + nrow(dataframe)) 
    stop("data frame too wide")
  options(width = oldwidth)
  if (suppress.X11.warnings) {
    messages.connection <- textConnection(".messages", open = "w", 
                                          local = TRUE)
    sink(messages.connection, type = "message")
    on.exit({
      sink(type = "message")
      close(messages.connection)
    })
  }
  base <- tcltk::tktoplevel()
  tcltk::tkwm.geometry(base, placement)
  tcltk::tkwm.title(base, {
    if (is.null(title)) 
      object.name
    else title
  })
  nrows <- length(zz) - 1
  if (is.numeric(rownumbers)) 
    rowname.text <- paste(rownumbers, row.names(dataframe))
  else if (rownumbers) 
    rowname.text <- paste(1:nrows, row.names(dataframe))
  else rowname.text <- row.names(dataframe)
  namewidth = max(nchar(rowname.text))
  yy <- substring(zz, 2 + max(nchar(row.names(dataframe))))
  datawidth <- max(nchar(yy))
  winwidth <- min(1 + datawidth, maxwidth)
  hdr <- tcltk::tktext(base, bg = colname.bgcolor, fg = colname.textcolor, 
                       font = font, height = 1, width = winwidth, takefocus = TRUE)
  ftr <- tcltk::tktext(base, bg = colname.bgcolor, fg = colname.textcolor, 
                       font = font, height = 1, width = winwidth, takefocus = TRUE)
  textheight <- min(maxheight, nrows)
  txt <- tcltk::tktext(base, bg = body.bgcolor, fg = body.textcolor, 
                       font = font, height = textheight, width = winwidth, setgrid = 1, 
                       takefocus = TRUE)
  lnames <- tcltk::tktext(base, bg = rowname.bgcolor, fg = rowname.textcolor, 
                          font = font, height = textheight, width = namewidth, 
                          takefocus = TRUE)
  rnames <- tcltk::tktext(base, bg = rowname.bgcolor, fg = rowname.textcolor, 
                          font = font, height = textheight, width = namewidth, 
                          takefocus = TRUE)
  xscroll <- tcltk::tkscrollbar(base, orient = "horizontal", 
                                repeatinterval = 1, command = function(...) {
                                  tcltk::tkxview(txt, ...)
                                  tcltk::tkxview(hdr, ...)
                                  tcltk::tkxview(ftr, ...)
                                })
  string.to.vector <- function(string.of.indices) {
    string.of.indices <- tcltk::tclvalue(string.of.indices)
    as.numeric(strsplit(string.of.indices, split = " ")[[1]])
  }
  tcltk::tkconfigure(txt, xscrollcommand = function(...) {
    tcltk::tkset(xscroll, ...)
    xy <- string.to.vector(tcltk::tkget(xscroll))
    tcltk::tkxview.moveto(hdr, xy[1])
    tcltk::tkxview.moveto(ftr, xy[1])
  })
  tcltk::tkconfigure(hdr, xscrollcommand = function(...) {
    tcltk::tkset(xscroll, ...)
    xy <- string.to.vector(tcltk::tkget(xscroll))
    tcltk::tkxview.moveto(txt, xy[1])
    tcltk::tkxview.moveto(ftr, xy[1])
  })
  tcltk::tkconfigure(ftr, xscrollcommand = function(...) {
    tcltk::tkset(xscroll, ...)
    xy <- string.to.vector(tcltk::tkget(xscroll))
    tcltk::tkxview.moveto(hdr, xy[1])
    tcltk::tkxview.moveto(txt, xy[1])
  })
  yscroll <- tcltk::tkscrollbar(base, orient = "vertical", 
                                repeatinterval = 1, command = function(...) {
                                  tcltk::tkyview(txt, ...)
                                  tcltk::tkyview(lnames, ...)
                                  tcltk::tkyview(rnames, ...)
                                })
  tcltk::tkconfigure(txt, yscrollcommand = function(...) {
    tcltk::tkset(yscroll, ...)
    xy <- string.to.vector(tcltk::tkget(yscroll))
    tcltk::tkyview.moveto(lnames, xy[1])
    tcltk::tkyview.moveto(rnames, xy[1])
  })
  tcltk::tkconfigure(lnames, yscrollcommand = function(...) {
    tcltk::tkset(yscroll, ...)
    xy <- string.to.vector(tcltk::tkget(yscroll))
    tcltk::tkyview.moveto(txt, xy[1])
    tcltk::tkyview.moveto(rnames, xy[1])
  })
  tcltk::tkconfigure(rnames, yscrollcommand = function(...) {
    tcltk::tkset(yscroll, ...)
    xy <- string.to.vector(tcltk::tkget(yscroll))
    tcltk::tkyview.moveto(txt, xy[1])
    tcltk::tkyview.moveto(lnames, xy[1])
  })
  tcltk::tkbind(txt, "<B2-Motion>", function(x, y) {
    tcltk::tkscan.dragto(txt, x, y)
  })
  {
    copyText.hdr <- function() {
      tcltk::tcl("event", "generate", tcltk::.Tk.ID(hdr), 
                 "<<Copy>>")
    }
    tcltk::tkbind(hdr, "<Button-1>", function() tcltk::tkfocus(hdr))
    editPopupMenu.hdr <- tcltk::tkmenu(hdr, tearoff = FALSE)
    tcltk::tkadd(editPopupMenu.hdr, "command", label = "Copy <Ctrl-C>", 
                 command = copyText.hdr)
    RightClick.hdr <- function(x, y) {
      rootx <- as.integer(tcltk::tkwinfo("rootx", hdr))
      rooty <- as.integer(tcltk::tkwinfo("rooty", hdr))
      xTxt <- as.integer(x) + rootx
      yTxt <- as.integer(y) + rooty
      tcltk::tcl("tk_popup", editPopupMenu.hdr, xTxt, yTxt)
    }
    tcltk::tkbind(hdr, "<Button-3>", RightClick.hdr)
    tcltk::tkbind(hdr, "<Control-KeyPress-c>", copyText.hdr)
    copyText.ftr <- function() {
      tcltk::tcl("event", "generate", tcltk::.Tk.ID(ftr), 
                 "<<Copy>>")
    }
    tcltk::tkbind(ftr, "<Button-1>", function() tcltk::tkfocus(ftr))
    editPopupMenu.ftr <- tcltk::tkmenu(ftr, tearoff = FALSE)
    tcltk::tkadd(editPopupMenu.ftr, "command", label = "Copy <Ctrl-C>", 
                 command = copyText.ftr)
    RightClick.ftr <- function(x, y) {
      rootx <- as.integer(tcltk::tkwinfo("rootx", ftr))
      rooty <- as.integer(tcltk::tkwinfo("rooty", ftr))
      xTxt <- as.integer(x) + rootx
      yTxt <- as.integer(y) + rooty
      tcltk::tcl("tk_popup", editPopupMenu.ftr, xTxt, yTxt)
    }
    tcltk::tkbind(ftr, "<Button-3>", RightClick.ftr)
    tcltk::tkbind(ftr, "<Control-KeyPress-c>", copyText.ftr)
    copyText.txt <- function() {
      tcltk::tcl("event", "generate", tcltk::.Tk.ID(txt), 
                 "<<Copy>>")
    }
    tcltk::tkbind(txt, "<Button-1>", function() tcltk::tkfocus(txt))
    editPopupMenu.txt <- tcltk::tkmenu(txt, tearoff = FALSE)
    tcltk::tkadd(editPopupMenu.txt, "command", label = "Copy <Ctrl-C>", 
                 command = copyText.txt)
    RightClick.txt <- function(x, y) {
      rootx <- as.integer(tcltk::tkwinfo("rootx", txt))
      rooty <- as.integer(tcltk::tkwinfo("rooty", txt))
      xTxt <- as.integer(x) + rootx
      yTxt <- as.integer(y) + rooty
      tcltk::tcl("tk_popup", editPopupMenu.txt, xTxt, yTxt)
    }
    tcltk::tkbind(txt, "<Button-3>", RightClick.txt)
    tcltk::tkbind(txt, "<Control-KeyPress-c>", copyText.txt)
    copyText.lnames <- function() {
      tcltk::tcl("event", "generate", tcltk::.Tk.ID(lnames), 
                 "<<Copy>>")
    }
    tcltk::tkbind(lnames, "<Button-1>", function() tcltk::tkfocus(lnames))
    editPopupMenu.lnames <- tcltk::tkmenu(lnames, tearoff = FALSE)
    tcltk::tkadd(editPopupMenu.lnames, "command", label = "Copy <Ctrl-C>", 
                 command = copyText.lnames)
    RightClick.lnames <- function(x, y) {
      rootx <- as.integer(tcltk::tkwinfo("rootx", lnames))
      rooty <- as.integer(tcltk::tkwinfo("rooty", lnames))
      xTxt <- as.integer(x) + rootx
      yTxt <- as.integer(y) + rooty
      tcltk::tcl("tk_popup", editPopupMenu.lnames, xTxt, 
                 yTxt)
    }
    tcltk::tkbind(lnames, "<Button-3>", RightClick.lnames)
    tcltk::tkbind(lnames, "<Control-KeyPress-c>", copyText.lnames)
    copyText.rnames <- function() {
      tcltk::tcl("event", "generate", tcltk::.Tk.ID(rnames), 
                 "<<Copy>>")
    }
    tcltk::tkbind(rnames, "<Button-1>", function() tcltk::tkfocus(rnames))
    editPopupMenu.rnames <- tcltk::tkmenu(rnames, tearoff = FALSE)
    tcltk::tkadd(editPopupMenu.rnames, "command", label = "Copy <Ctrl-C>", 
                 command = copyText.rnames)
    RightClick.rnames <- function(x, y) {
      rootx <- as.integer(tcltk::tkwinfo("rootx", rnames))
      rooty <- as.integer(tcltk::tkwinfo("rooty", rnames))
      xTxt <- as.integer(x) + rootx
      yTxt <- as.integer(y) + rooty
      tcltk::tcl("tk_popup", editPopupMenu.rnames, xTxt, 
                 yTxt)
    }
    tcltk::tkbind(rnames, "<Button-3>", RightClick.rnames)
    tcltk::tkbind(rnames, "<Control-KeyPress-c>", copyText.rnames)
  }
  tcltk::tktag.configure(hdr, "notwrapped", wrap = "none")
  tcltk::tktag.configure(ftr, "notwrapped", wrap = "none")
  tcltk::tktag.configure(txt, "notwrapped", wrap = "none")
  tcltk::tktag.configure(lnames, "notwrapped", wrap = "none")
  tcltk::tktag.configure(rnames, "notwrapped", wrap = "none")
  tcltk::tkinsert(txt, "end", paste(paste(yy[-1], collapse = "\n"), 
                                    sep = ""), "notwrapped")
  tcltk::tkgrid(txt, row = 1, column = 1, sticky = "nsew")
  if ("top" %in% colname.bar) {
    tcltk::tkinsert(hdr, "end", paste(yy[1], sep = ""), "notwrapped")
    tcltk::tkgrid(hdr, row = 0, column = 1, sticky = "ew")
  }
  if ("bottom" %in% colname.bar) {
    tcltk::tkinsert(ftr, "end", paste(yy[1], sep = ""), "notwrapped")
    tcltk::tkgrid(ftr, row = 2, column = 1, sticky = "ew")
  }
  if ("left" %in% rowname.bar) {
    tcltk::tkinsert(lnames, "end", paste(rowname.text, collapse = "\n"), 
                    "notwrapped")
    tcltk::tkgrid(lnames, row = 1, column = 0, sticky = "ns")
  }
  if ("right" %in% rowname.bar) {
    tcltk::tkinsert(rnames, "end", paste(rowname.text, collapse = "\n"), 
                    "notwrapped")
    tcltk::tkgrid(rnames, row = 1, column = 2, sticky = "ns")
  }
  tcltk::tkconfigure(txt, state = "disabled")
  tcltk::tkconfigure(lnames, state = "disabled")
  tcltk::tkconfigure(rnames, state = "disabled")
  if (maxheight < nrows) {
    tcltk::tkgrid(yscroll, row = 1, column = 3, sticky = "ns")
  }
  if (maxwidth < datawidth) {
    tcltk::tkgrid(xscroll, row = 3, column = 1, sticky = "ew")
  }
  tcltk::tkgrid.rowconfigure(base, 1, weight = 1)
  tcltk::tkgrid.columnconfigure(base, 1, weight = 1)
  tcltk::tkwm.maxsize(base, 1 + datawidth, nrows)
  tcltk::tkwm.minsize(base, 1 + nchar(names(dataframe)[1]), 
                      1)
  
  onClose <- function(){
      if ("Rcmdr" %in% loadedNamespaces()){
          open.showData.windows <- Rcmdr::getRcmdr("open.showData.windows", fail=FALSE)
          if (!is.null(open.showData.windows)){
              open.showData.windows[[object.name]] <- NULL
              Rcmdr::putRcmdr("open.showData.windows", open.showData.windows)
          }
      }
      tcltk::tkdestroy(base)
  }
  tcltk::tkwm.protocol(base, "WM_DELETE_WINDOW", onClose)
  
  invisible(base)
}
