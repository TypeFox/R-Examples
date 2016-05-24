### R code from vignette source 'ex-tcltk-t-test.Rnw'

###################################################
### code chunk number 1: ex-tcltk-t-test.Rnw:16-19
###################################################
## t.test dialog
## using basic widgets -- no entry widgets yet
library(tcltk)


###################################################
### code chunk number 2: ex-tcltk-t-test.Rnw:22-40
###################################################
## helper functions
## not shown

get_numeric_vars <- function(DF) {
  if(missing(DF))
    return(c(""))
  ProgGUIinR:::find_vars(DF, is.numeric)
}
get_two_level_factor <- function(DF) {
  if(missing(DF))
    return(c(""))
  nms <- names(DF)
  ind <- sapply(DF, function(i) length(levels(as.factor(i))) == 2)
  if(length(ind) > 0)
    nms[ind]
  else
    c("")
}


###################################################
### code chunk number 3: dataModel
###################################################
e <- new.env()
e$x <- tclVar(""); e$f <- tclVar(""); e$data <- tclVar("")
e$mu <- tclVar(0); e$alternative <- tclVar("two.sided")
e$conf.level <- tclVar(95); e$var.equal <- tclVar(FALSE)


###################################################
### code chunk number 4: ex-tcltk-t-test.Rnw:58-88
###################################################
## We don't show the function runTTest.
## It is a bit long, as care must be taken as it isn't clear if a formula should be used.  
runTTest <- function() {
  l <- lapply(e, tclvalue)
  
  ## ugly function to run t test
  if(l$data == "" || l$x == "")
    return()

  l$data <- get(l$data, envir = .GlobalEnv)

  if(l$f != "") {
    l$formula <- formula(paste(l$x,l$f, sep = "~"))
    l$x <- l$f <- NULL
    l$mu <- NULL
    l$var.equal <- as.logical(as.numeric(l$var.equal))

    TTest <- stats:::t.test.formula
  } else {
    l$x <- l$data[, l$x]
    l$f <- NULL
    l$mu = as.numeric(l$mu)
    l$var.equal <- NULL
    
    TTest <- stats:::t.test.default
  }
  l$conf.level <- as.numeric(l$conf.level)/100
  ret <- capture.output(do.call("TTest", l))
  cat(paste(ret, collapse = "\n"))
}


###################################################
### code chunk number 5: notShown
###################################################
### GUI Our standard setup
window <- tktoplevel()
tkwm.title(window, "t-test Dialog")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 6: layout
###################################################
label_frame <- ttklabelframe(frame, text = "t.test()", 
                             padding = 10)
tkpack(label_frame, expand = TRUE, fill = "both", 
       padx = 5, pady = 5)


###################################################
### code chunk number 7: ex-tcltk-t-test.Rnw:115-119
###################################################
tkgrid.columnconfigure(label_frame, 0, weight = 1)
tkgrid.columnconfigure(label_frame, 1, weight = 10)
tkgrid.columnconfigure(label_frame, 2, weight = 1)
tkgrid.columnconfigure(label_frame, 1, weight = 10)


###################################################
### code chunk number 8: tkgridHelper
###################################################
put_label <- function(parent, text, row, column) {
  label <- ttklabel(parent, text = text)
  tkgrid(label, row = row, column = column, sticky = "e")
}


###################################################
### code chunk number 9: readonly
###################################################
put_label(label_frame, "data:",0,0)
data_combo <- ttkcombobox(label_frame, state = "readonly", 
                         values = ProgGUIinR:::avail_dfs(), 
                         textvariable = e$data)
tkgrid(data_combo, row = 0, column = 1, sticky="ew", padx = 2)
tkfocus(data_combo)                      # give focus


###################################################
### code chunk number 10: notShown
###################################################
## not shown
put_label(label_frame, "x:",1,0)
x_combo <-  ttkcombobox(label_frame, 
                       values = get_numeric_vars(), textvariable = e$x)
tkgrid(x_combo, row = 1, column = 1, sticky = "ew", padx = 2)


###################################################
### code chunk number 11: notShown
###################################################
## not shown
put_label(label_frame, "~ f:",1,2)
factor_combo <-  ttkcombobox(label_frame, values = get_two_level_factor(), textvariable = e$f)
tkgrid(factor_combo, row = 1, column = 3, sticky = "ew", padx = 2)


###################################################
### code chunk number 12: mu
###################################################
put_label(label_frame, "mu:", 2, 0)
mu_combo <-  ttkentry(label_frame,  textvariable = e$mu)
tkgrid(mu_combo, row = 2, column = 1, sticky = "ew", padx = 2)


###################################################
### code chunk number 13: ex-tcltk-t-test.Rnw:172-174
###################################################
ttkscale <- function(parent, ...) tkwidget(parent, "ttk::scale", ...)
tkspinbox <- function(parent, ...) tkwidget(parent, "tk::spinbox", ...)


###################################################
### code chunk number 14: ex-tcltk-t-test.Rnw:176-185
###################################################
put_label(label_frame, "alternative:", 3, 0)
rb_frame <- ttkframe(label_frame)
sapply(c("two.sided","less","greater"), function(i) {
  radio_button <- 
    ttkradiobutton(rb_frame, variable = e$alternative, 
                   text = i, value = i)
  tkpack(radio_button, side = "left")
})
tkgrid(rb_frame, row = 3, column = 1, sticky = "ew", padx = 2)


###################################################
### code chunk number 15: ex-tcltk-t-test.Rnw:191-206
###################################################
put_label(label_frame, "conf.level:", 3, 2)
conf_level_frame <- ttkframe(label_frame)
tkgrid(conf_level_frame, row = 3, column = 3, columnspan = 2, 
       sticky = "ew", padx = 2)
##
conf_level_scale <- ttkscale(conf_level_frame, 
                     from = 75, to = 100,  
                     variable = e$conf.level)
conf_level_spin <- tkspinbox(conf_level_frame, 
                     from = 75, to = 100, increment = 1, 
                     textvariable = e$conf.level, width = 5)
##
tkpack(conf_level_scale, expand = TRUE, fill = "y",
       side = "left")
tkpack(conf_level_spin, side = "left")


###################################################
### code chunk number 16: ex-tcltk-t-test.Rnw:210-215
###################################################
put_label(label_frame, "var.equal:", 4, 0)
var_equal_check <- 
  ttkcheckbutton(label_frame, variable = e$var.equal)
tkgrid(var_equal_check, row = 4, column = 1, stick = "w", 
       padx = 2)


###################################################
### code chunk number 17: ok-cancel
###################################################
button_frame <- ttkframe(frame)
cancel_button <- ttkbutton(button_frame, text = "cancel")
ok_button <- ttkbutton(button_frame, text = "ok")
#
tkpack(button_frame, fill = "x", padx = 5, pady = 5)
tkpack(ttklabel(button_frame, text = " "), expand = TRUE,
       fill = "y", side = "left")               # add a spring
sapply(list(cancel_button, ok_button), tkpack, 
       side = "left", padx = 6)


###################################################
### code chunk number 18: ex-tcltk-t-test.Rnw:237-241
###################################################
tkconfigure(ok_button, command = runTTest)
tkconfigure(cancel_button, 
            command = function() tkdestroy(window))
tkbind("TButton", "<Return>", function(W) tcl(W, "invoke"))


###################################################
### code chunk number 19: ex-tcltk-t-test.Rnw:251-297
###################################################
update_ui <- function() {
  dfName <- tclvalue(e$data)
  curDfs <- ProgGUIinR:::avail_dfs()
  tkconfigure(data_combo, values = curDfs)
  if(!dfName %in% curDfs) {
    dfName <- ""
    tclvalue(e$data) <- ""
  }

  if(dfName == "") {
    ## 3 ways to disable needed!!
    x <- list(x_combo, factor_combo, mu_combo,  
              conf_level_scale, var_equal_check, ok_button)
    sapply(x, function(W) tcl(W, "state", "disabled"))
    sapply(as.character(tkwinfo("children", rb_frame)), 
           function(W) tcl(W, "state", "disabled"))
    tkconfigure(conf_level_spin, state = "disabled")
  } else {
    ## enable univariate, ok
    sapply(list(x_combo,mu_combo,conf_level_scale,ok_button),
           function(W) tcl(W, "state", "!disabled"))
    sapply(as.character(tkwinfo("children", rb_frame)), 
           function(W) tcl(W, "state", "!disabled"))
    tkconfigure(conf_level_spin, state = "normal")
    
    DF <- get(dfName, envir = .GlobalEnv)
    numVars <- get_numeric_vars(DF)
    tkconfigure(x_combo, values = numVars)
    if(! tclvalue(e$x) %in% numVars)
      tclvalue(e$x) <- ""

    ## bivariate
    avail_factors <- get_two_level_factor(DF)
    sapply(list(factor_combo, var_equal_check),
           function(W) {
             val <- if(length(avail_factors)) "!" else ""
             tcl(W, "state", sprintf("%sdisabled", val))
           })
    tkconfigure(factor_combo, values = avail_factors)
    if(!tclvalue(e$f) %in% avail_factors)
      tclvalue(e$f) <- ""
      
         }
}
update_ui()
tkbind(data_combo, "<<ComboboxSelected>>", update_ui)


###################################################
### code chunk number 20: digest
###################################################
require(digest)
create_function <- function() {
  .DFs <- digest(ProgGUIinR:::avail_dfs())
  f <- function(...) {
    if((val <- digest(ProgGUIinR:::avail_dfs())) != .DFs) {
      .DFs <<- val
      update_ui()
    }
    return(TRUE)
  }
}


###################################################
### code chunk number 21: taskcallback (eval = FALSE)
###################################################
## id <- addTaskCallback(create_function())


