
###################################################
### code chunk number 167: ttkscale
###################################################
ttkscale <- function(parent, ...) 
  tkwidget(parent, "ttk::scale", ...)


###################################################
### code chunk number 168: ttksliderclass
###################################################
Slider <-
  setRefClass("TtkSlider",
     fields = c("frame", "widget", "var", "x", "FUN"),
     methods = list(
       initialize = function(parent, x, ...) {
         initFields(x = x, var = tclVar(1),
                    FUN = NULL, frame = ttkframe(parent))
         widget <<- ttkscale(frame, from = 1, to = length(x),
                       variable = var, orient = "horizontal")
         ## For this widget, the callback is passed a value 
         ## which we ignore here
         tkconfigure(widget, command = function(...) {
           if(is.function(FUN)) FUN(.self)
         })
         layout_gui()
         callSuper(...)
       },
       layout_gui = function() {         
         tkgrid(widget, row = 0, column = 0, columnspan = 3, 
                sticky = "we")
         tkgrid(ttklabel(frame, text = x[1]), 
                row = 1, column = 0)
         tkgrid(ttklabel(frame, text = x[length(x)]), 
                row = 1, column = 2)
         tkgrid.columnconfigure(frame, 1, weight = 1)
       },
       add_callback = function(FUN) FUN <<- FUN,
       get_value = function() x[as.numeric(tclvalue(var))],
       set_value = function(value) {
         "Set value. Value must be in x"
         ind <- match(value, x)
         if(!is.na(ind)) {
           var_local <- var
           tclvalue(var_local) <- ind
         }
       }
       ))


###################################################
### code chunk number 169: Widgets.Rnw:495-505
###################################################
window <- tktoplevel()
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
x <- seq(0,1,by = 0.05)
##
slider <- Slider$new(parent = window, x = x)
tkpack(slider$frame, expand = TRUE, fill = "x", anchor = "n")
##
slider$set_value(0.5)
print(slider$get_value())


###################################################
### code chunk number 170: use-slider-command
###################################################
slider$add_callback(function(obj) print(obj$get_value()))

