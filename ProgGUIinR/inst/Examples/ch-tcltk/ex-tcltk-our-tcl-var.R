
###################################################
### code chunk number 129: ourTclVar
###################################################
our_tcl_var <- function(...) {
  var <- tclVar(...)
  .TkRoot$env[[as.character(var)]] <- var
  var
}
## lookup function
get_tcl_var_by_id <- function(id) {
  .TkRoot$env[[as.character(id)]]
}


###################################################
### code chunk number 130: our_tcl_varExample
###################################################
window <- tktoplevel(); tkwm.title(window, "Check button example")
frame <- ttkframe(window); tkpack(frame, expand = TRUE, fill = "both")
value_var <- our_tcl_var(TRUE)


###################################################
### code chunk number 131: Widgets.Rnw:142-146
###################################################
callback <- function(W) {
  id <- tkcget(W, "variable" = NULL)
  print(get_tcl_var_by_id(id))
}


###################################################
### code chunk number 132: Widgets.Rnw:149-154
###################################################
## command does not pass back in widget, we use tkbind instead
label_var <- tclVar("check button label")
check_button<- 
  ttkcheckbutton(frame, variable = value_var, 
                 textvariable = label_var)
tkpack(check_button)

tkbind(check_button, "<Button-1>", callback)
