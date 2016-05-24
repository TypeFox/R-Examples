### R code from vignette source 'ex-tcltk-validation.Rnw'

###################################################
### code chunk number 1: ex-tcltk-validation.Rnw:6-40
###################################################
## Example of using validation to adjust the date
## in case a user doesn't use desired format

## Docs on validation
## VALIDATION
## The -validate, -validatecommand, and -invalidcommand options are used to enable entry widget validation.
## VALIDATION MODES
## There are two main validation modes: prevalidation, in which the -validatecommand is evaluated prior to each edit and the return value is used to determine whether to accept or reject the change; and revalidation, in which the -validatecommand is evaluated to determine whether the current value is valid.

## The -validate option determines when validation occurs; it may be set to any of the following values:

## none
##     Default. This means validation will only occur when specifically requested by the validate widget command.

## key
##     The entry will be prevalidated prior to each edit (specifically, whenever the insert or delete widget commands are called). If prevalidation fails, the edit is rejected.

## focus
##     The entry is revalidated when the entry receives or loses focus.

## focusin
##     The entry is revalidated when the entry receives focus.

## focusout
##     The entry is revalidated when the entry loses focus.

## all
##     Validation is performed for all above conditions.

## The -invalidcommand is evaluated whenever the -validatecommand returns a false value.

## The -validatecommand and -invalidcommand may modify the entry widget's value via the widget insert or delete commands, or by setting the linked -textvariable. If either does so during prevalidation, then the edit is rejected regardless of the value returned by the -validatecommand.

## If -validatecommand is empty (the default), validation always succeeds.


###################################################
### code chunk number 2: ex-tcltk-validation.Rnw:44-49
###################################################
## test of validation command
## tricky bit is that validation commands must return TRUE or FALSE
## we can do this with tcl("eval","FALSE") or .Tcl("expr false")

library(tcltk)


###################################################
### code chunk number 3: ex-tcltk-validation.Rnw:85-92
###################################################
date_patterns <- c()
for(i in list(c("%m","%d","%Y"),        # U.S. style
              c("%m","%d","%y"))) {
  for(j in c("/","-"," ") )
    date_patterns[length(date_patterns)+1] <- 
      paste(i,sep = "", collapse = j)
}


###################################################
### code chunk number 4: setValidDateCallback
###################################################
is_valid_date <- function(W, P) { # P is the current value
  for(i in date_patterns) {
    date <- try( as.Date(P, format = i), silent = TRUE)
    if(!inherits(date, "try-error") && !is.na(date)) {
      tkconfigure(W, foreground = "black")  # or use style
      tkdelete(W, 0,"end")
      tkinsert(W, 0, format(date, format = "%m/%d/%y"))
      return(tcl("expr","TRUE"))        
    } 
  }
  return(tcl("expr","FALSE"))
}


###################################################
### code chunk number 5: setInvalidCallback
###################################################
indicate_invalid_date <- function(W) {
  tkconfigure(W, foreground = "red")
  tcl("expr", "TRUE")
}


###################################################
### code chunk number 6: notShown
###################################################
## A simple GUI to show the entry widget.
window <- tktoplevel(); tkwm.title(window, "Validation of date example")
frame <- ttkframe(window, padding = c(3,3,12,12)); tkpack(frame, expand = TRUE, fill = "both")
tkpack(ttklabel(frame, text = "Enter a date (mm/dd/yy):"), side = "left", padx = 2)


###################################################
### code chunk number 7: entryWidgetWithValidation
###################################################
entry <- ttkentry(frame, validate = "focusout",
                  validatecommand = is_valid_date,
                  invalidcommand = indicate_invalid_date)
button <- ttkbutton(frame, text = "click")  # focus target
sapply(list(entry, button), tkpack, side = "left", padx = 2)


