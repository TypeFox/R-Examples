slider <- function (sl.functions, sl.names, sl.mins, sl.maxs, sl.deltas,
                    sl.defaults, but.functions, but.names, no, set.no.value,
                    obj.name, obj.value, reset.function, title) {
  # slightly modified by J. Fox from the TeachingDemos package
  requireNamespace('tcltk', quietly = TRUE)
  
  if (!missing(no))
      return(as.numeric(tcltk::tclvalue(get(paste("slider", no, sep = ""),
                                     envir = slider.env))))
  if (!missing(set.no.value)) {
      try(eval(parse(text = paste("tcltk::tclvalue(slider", set.no.value[1],
                     ")<-", set.no.value[2], sep = "")), envir = slider.env))
      return(set.no.value[2])
  }
  if (!exists("slider.env"))
      slider.env <<- new.env()
  if (!missing(obj.name)) {
      if (!missing(obj.value))
          assign(obj.name, obj.value, envir = slider.env)
      else obj.value <- get(obj.name, envir = slider.env)
      return(obj.value)
  }
  if (missing(title))
      title <- "slider control widget"
  
  nt <- tcltk::tktoplevel()
  tcltk::tkwm.title(nt, title)
  tcltk::tkwm.geometry(nt, "+0+0")
  if (missing(sl.names))
      sl.names <- NULL
  if (missing(sl.functions))
      sl.functions <- function(...) {
      }
  for (i in seq(sl.names)) {
      eval(parse(text = paste("assign('slider", i,
                 "',tcltk::tclVar(sl.defaults[i]),envir=slider.env)",
                 sep = "")))
      tcltk::tkpack(fr <- tcltk::tkframe(nt))
      lab <- tcltk::tklabel(fr, text = sl.names[i], width = "25")
      sc <- tcltk::tkscale(fr, from = sl.mins[i], to = sl.maxs[i],
                    showvalue = T, resolution = sl.deltas[i], orient = "horiz")
      tcltk::tkpack(lab, sc, side = "right")
      assign("sc", sc, envir = slider.env)
      eval(parse(text = paste("tcltk::tkconfigure(sc,variable=slider",
                 i, ")", sep = "")), envir = slider.env)
      sl.fun <- if (length(sl.functions) > 1)
          sl.functions[[i]]
      else sl.functions
      if (!is.function(sl.fun))
          sl.fun <- eval(parse(text = paste("function(...){",
                               sl.fun, "}")))
      tcltk::tkconfigure(sc, command = sl.fun)
  }
  assign("slider.values.old", sl.defaults, envir = slider.env)
  tcltk::tkpack(f.but <- tcltk::tkframe(nt), fill = "x")
  tcltk::tkpack(tcltk::tkbutton(f.but, text = "Exit", command = function()
                  tcltk::tkdestroy(nt)),
         side = "right")
  if (!missing(reset.function)){
      if (!is.function(reset.function))
          reset.function <- eval(parse(text = paste("function(...){",
                                       reset.function, "}")))
      tcltk::tkpack(tcltk::tkbutton(f.but, text = "Reset", command = function() {
          for (i in seq(sl.names)) eval(parse(text =
                                              paste("tcltk::tclvalue(slider",
                                                    i, ")<-", sl.defaults[i], sep = "")), envir = slider.env)
          reset.function()
      }), side = "right")
  }
  if (missing(but.names))
      but.names <- NULL
  for (i in seq(but.names)) {
      but.fun <- if (length(but.functions) > 1)
          but.functions[[i]]
      else but.functions
      if (!is.function(but.fun))
          but.fun <- eval(parse(text = paste("function(...){",
                                but.fun, "}")))
      tcltk::tkpack(tcltk::tkbutton(f.but, text = but.names[i], command = but.fun),
             side = "left")
      cat("button", i, "eingerichtet")
  }
  invisible(nt)
}
