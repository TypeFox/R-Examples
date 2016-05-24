### R code from vignette source 'ex-tcltk-initial-message.Rnw'

###################################################
### code chunk number 1: InitialMsg
###################################################
## "R5" class for ttk entry with initial message
library(tcltk)


###################################################
### code chunk number 2: ex-tcltk-initial-message.Rnw:16-23
###################################################
setOldClass(c("tkwin", "tclVar"))
TtkEntry <- setRefClass("TtkEntry",
                        fields = list(
                          entry = "tkwin",     # entry
                          tcl_var  = "tclVar", # text variable
                          init_msg = "character"
                          ))


###################################################
### code chunk number 3: init-msg-style
###################################################
.Tcl("ttk::style configure Gray.TEntry -foreground gray") 


###################################################
### code chunk number 4: init-msg-methods
###################################################
TtkEntry$methods(
                 is_init_msg = function() {
                   "Is the init text showing?"
                   as.character(tclvalue(tcl_var)) == init_msg
                 },
                 hide_init_msg = function() {
                   "Hide the initial text"
                   if(is_init_msg()) {
                     tkconfigure(entry, style = "TEntry")
                     set_text("", hide = FALSE)
                   }
                 },
                 show_init_msg = function() {
                   "Show the initial text"
                   tkconfigure(entry, style = "Gray.TEntry")
                   set_text(init_msg, hide = FALSE)
                 })


###################################################
### code chunk number 5: get-set-text
###################################################
TtkEntry$methods(
                 set_text = function(text, hide = TRUE) {
                   "Set text into widget"
                   if(hide) hide_init_msg()
                   tcl_var_local <- tcl_var   # avoid warning
                   tclvalue(tcl_var_local) <- text
                 },
                 get_text = function() {
                   "Get the text value"
                   if(!is_init_msg())
                     as.character(tclvalue(tcl_var))
                   else
                     ""
                 })


###################################################
### code chunk number 6: add-bindings
###################################################
TtkEntry$methods(
                 add_bindings = function() {
                   "Add focus bindings to make this work"
                   tkbind(entry, "<FocusIn>", hide_init_msg)
                   tkbind(entry, "<FocusOut>", function() {
                     if(nchar(get_text()) == 0)
                       show_init_msg()
                   })
                 })


###################################################
### code chunk number 7: ex-tcltk-initial-message.Rnw:96-109
###################################################
TtkEntry$methods(
   initialize = function(parent, text, init_msg = "", ...) {
     tcl_var <<- tclVar()
     entry <<- ttkentry(parent, textvariable = tcl_var)
     init_msg <<- init_msg
     ##
     if(missing(text))
       show_init_msg()
     else
       set_text(text)
     add_bindings()
     callSuper(...)
   })


###################################################
### code chunk number 8: ex-tcltk-initial-message.Rnw:116-126
###################################################
window <- tktoplevel()
widget <- TtkEntry$new(parent = window, 
                       init_msg = "type value here")
tkpack(widget$entry)
#
button <- ttkbutton(window, text = "focus out onto this", 
               command = function() {
                 print(widget$get_text())
               })
tkpack(button)


