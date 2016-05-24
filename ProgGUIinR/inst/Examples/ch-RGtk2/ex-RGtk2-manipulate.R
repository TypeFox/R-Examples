### R code from vignette source 'ex-RGtk2-manipulate.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-manipulate.Rnw:1-27
###################################################
## manipulate for RGtk2
#
# Original license for manipulate package
#
# Copyright (C) 2009-11 by RStudio, Inc.
#
# This program is licensed to you under the terms of version 3 of the
# GNU Affero General Public License. This program is distributed WITHOUT
# ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
# AGPL (http://www.gnu.org/licenses/agpl-3.0.txt) for more details.
#
#

## THe main point of AGPL:
##   The GNU Affero General Public License is designed specifically to
## ensure that, in such cases, the modified source code becomes available
## to the community.  It requires the operator of a network server to
## provide the source code of the modified version running there to the
## users of that server.  Therefore, public use of a modified version, on
## a publicly accessible server, gives the public access to the source
## code of the modified version.
## This is satisfied by the ProgGUIinR package, which will contain this entire example.

require(RGtk2)
require(cairoDevice)


###################################################
### code chunk number 2: resolveVariableArguments
###################################################
## Not shown
resolveVariableArguments <- function(args) {
  # if the first argument is an unnamed list then just use this list
  if ( (length(args) == 1L) &&
       is.list(args[[1L]])  &&
       (is.null(names(args)) || (names(args)[[1L]] == "")) )  {
    return (args[[1L]])
  } else {
    return (args)
  }
}


###################################################
### code chunk number 3: manipulate_example (eval = FALSE)
###################################################
## manipulate(## expression
##            plot(cars, xlim = c(x.min, x.max), type = type, 
##                 axes = axes, ann = label),
##            ## controls
##            x.min = slider(0, 15),
##            x.max = slider(15, 30, initial = 25),
##            type = picker("p", "l", "b", "c", "o", "h", "s"),
##            axes = checkbox(TRUE, label = "Draw Axes"),
##            label = checkbox(FALSE, label = "Draw Labels")
##            )


###################################################
### code chunk number 4: ManipulateClass
###################################################
Manipulate <- setRefClass("Manipulate",
                          fields=list(
                            .code="ANY",
                            .controls="list"
                            ))


###################################################
### code chunk number 5: manipulate_validate_controls
###################################################
Manipulate$methods(validate_controls = function() {
  "Validate that controls are specified properly"
  ## validate that all controls have unique names
  controlNames <- names(.controls)
  duplicatedIndex <- anyDuplicated(controlNames)
  if (duplicatedIndex > 0)
    stop(paste("duplicated control name:", controlNames[[duplicatedIndex]]))
  
  ## iterate over the names and controls, adding the default values to the env
  for (name in names(.controls)) {
    ## check the name
    if (name == "")
      stop("all controls passed to manipulate must be named")
    ## confirm that this is in fact a control
    if(!is(.controls[[name]], "ManipulateControls"))
      stop(paste("argument", name, "is not a control"))
    ## default label is control name
    if(length(.controls[[name]]$label) == 0) 
      .controls[[name]]$label <<- name
  }
})


###################################################
### code chunk number 6: Manipulate_change_handler
###################################################
Manipulate$methods(
           get_values = function() {
             sapply(.controls, 
                    function(control) control$get_value(), 
                    simplify=FALSE)     # return a list
           },
           change_handler = function(...) {
             "Evaluate code with current values"
             values <- get_values()
             result <- withVisible(eval(.code, envir=values))
             if (result$visible) {
               eval(print(result$value))
             }
           })


###################################################
### code chunk number 7: Manipulate_execute
###################################################
Manipulate$methods(  
           execute=function() {
             "Make the GUI"
             window <- gtkWindow(show=FALSE)
             window$setTitle("ManipulateR")
             ## Set up graphic device
             hpaned <- gtkHPaned()
             window$add(hpaned)
             device <- gtkDrawingArea()
             device$setSizeRequest(480, 480)
             asCairoDevice(device)
             hpaned$add(device)
             ## Controls frame
             frame <- gtkFrame("Controls")
             control_table <- gtkTableNew()
             control_table$setHomogeneous(FALSE)
             control_table['column-spacing'] <- 10
             ## insert horizontal strut
             control_table$attach(strut <- gtkHBox(), 1,2,0,1,
                           xoptions="", yoptions="shrink")
             strut$setSizeRequest(75, -1)
             frame$add(control_table)
             hpaned$add(frame)
             ## add each control
             sapply(.controls, function(control) {
               control$make_gui(cont=control_table, 
                                handler=.self$change_handler)
             })
             window$show()
             change_handler()                    # initial
           })


###################################################
### code chunk number 8: Manipulate_Initialize
###################################################
Manipulate$methods(  
           initialize = function(code, ...) {
             controls <- resolveVariableArguments(list(...))
             initFields(.code = code,
                        .controls = controls)
             validate_controls()
             callSuper()
           })


###################################################
### code chunk number 9: manipulate_constructor
###################################################
manipulate <- function(`_expr`,...) {
  manip <- Manipulate$new(substitute(`_expr`),...)
  manip$execute()
}


###################################################
### code chunk number 10: ManipulateControls
###################################################
ManipulateControls <- setRefClass("ManipulateControls",
                        fields=list(
                          l="list",
                          widget = "ANY",
                          label="ANY",
                          initial="ANY"
                          ))


###################################################
### code chunk number 11: MC_Interface
###################################################
ManipulateControls$methods(
            validate_inputs = function(...) {
              "Validate input code"
            },
            get_value = function(...) {
              "Get value of widget"
            })


###################################################
### code chunk number 12: MC_make_gui
###################################################
ManipulateControls$methods(make_gui = function(cont) {
            "Create widget, then add to table"
            ## cont a GtkTable instance
            nrows <- cont['n-rows']
            label_widget <- gtkLabel(label)
            label_widget['xalign'] <- 1
            cont$attach(label_widget, 0, 1, nrows, nrows + 1,
                        xoptions = "shrink", yoptions="shrink"
                        )
            cont$attach(widget, 1, 2, nrows, nrows + 1,
                        xoptions = c("expand", "fill"),
                        yoptions = "")
          })


###################################################
### code chunk number 13: Slider_constructor
###################################################
slider <- function(min, max, initial = min, label=NULL, 
                   step = -1, ticks = TRUE) {
  Slider$new(min, max, initial = initial, label = label, 
             step = step, ticks = ticks)
}


###################################################
### code chunk number 14: Slider
###################################################
Slider <- setRefClass("Slider",
                      contains = "ManipulateControls")


###################################################
### code chunk number 15: Slider_validate
###################################################
Slider$methods(validate_inputs = function(min, max, initial, step, ticks, label) {
                            ## validate inputs
                          if (!is.numeric(initial) || !is.numeric(min) || !is.numeric(max))
                            stop("min, max, amd initial must all be numeric values")
                          else if (initial < min)
                            stop(paste("slider initial value", initial, "is less than the specified minimum"))
                          else if (initial > max)
                            stop(paste("slider initial value", initial, "is greater than the specified maximum"))
                          else if (min > max)
                            stop(paste("slider maximum is greater than minimum"))
                          else if ( !is.null(step) ) {
                            if ( !is.numeric(step) )
                              stop("step is not a numeric value")
                            if ( step > (max - min) )
                              stop("step is greater than range")
                          } else if ( !is.logical(ticks) )
                            stop("ticks is not a logical value")
                        })


###################################################
### code chunk number 16: Slider_initialize
###################################################
Slider$methods(
       initialize = function(min, max, initial = min, 
         label = NULL, step = -1, ticks = TRUE, ...) {
           validate_inputs(min, max, initial, step, ticks)
           ## create slider and return it
           slider <- list(type = 0,
                          min = min,
                          max = max,
                          step = step,
                          ticks = ticks)
           initFields(l = slider, label = label, 
                      initial = initial)
           callSuper()
         })


###################################################
### code chunk number 17: Slider_make_gui
###################################################
Slider$methods(
       make_gui = function(cont, handler, ...) {
         widget <<- gtkHScale(min = l$min, max = l$max, 
                              step = l$step)
         widget$setValue(initial)
         gSignalConnect(widget, "value-changed", handler)
         callSuper(cont)
       },
       get_value = function() {
         as.numeric(widget$getValue())
       })


###################################################
### code chunk number 18: Picker
###################################################
## Not shown -- too long
Picker <- setRefClass("Picker",
                      contains="ManipulateControls",
                      methods=list(
                        initialize=function(..., initial = NULL, label = NULL) {
                          
                          ## get values
                          values <- resolveVariableArguments(list(...))
                          
                          ## get value names
                          valueNames <- names(values)
                          if (is.null(valueNames))
                            valueNames <- character(length(values))
                          
                          ## default missing names to choice values
                          missingNames <- valueNames == ""
                          valueNames[missingNames] <- paste(values)[missingNames]
                          names(values) <- valueNames
                          validate_inputs(values, valueNames, initial,label)
                          
                          if(is.null(initial)) 
                            initial <<- valueNames[1]
                          else
                            initial <<- initial
                          ## create picker
                          picker <- list(type = 1,
                                         choices = valueNames,
                                         values = values
                                         )
                          initFields(l=picker, label=label)
                          callSuper()
                        },
                        make_gui=function(cont, handler, ...) {
                          widget <<- gtkComboBoxNewText()
                          sapply(l$choices, widget$appendText) # visible ones
                          ## initialize
                          ind <- match(initial, l$choices)
                          if(is.na(ind)) ind <- 1
                          widget$setActive(ind - 1L)
                          ## add signal
                          gSignalConnect(widget, "changed", handler)
                          callSuper(cont)
                        },
                        get_value=function() {
                          ind <- widget$getActive()
                          l$values[[ind + 1L]]
                        },
                        validate_inputs=function(values, valueNames, initial,label) {
                          if ( length(values) < 1 ) {
                            stop("picker choices must contain at least one value")
                          } else if ( length(valueNames) != length(unique(valueNames)) ) {
                            stop("picker choices must have unique names (duplicate detected)")
                          } else if ( !is.null(initial) ) {
                            if (length(initial) != 1)
                              stop("initial must be a single object")
                            else if ( !(as.character(initial) %in% valueNames) )
                              stop("initial doesn't match one of the supplied choices") 
                          }
                        }
                        
                        ))
picker <- function(..., initial = NULL, label = NULL) 
  Picker$new(..., initial=initial, label=label)


###################################################
### code chunk number 19: Checkbox
###################################################
Checkbox <- setRefClass("Checkbox", contains="ManipulateControls")
Checkbox$methods(validate_inputs=function(initial, label) {
                   if ( !is.logical(initial) )
                     stop("initial must be a logical")
                 })


###################################################
### code chunk number 20: ex-RGtk2-manipulate.Rnw:437-453
###################################################
Checkbox$methods(
         initialize = function(initial=FALSE, label=  NULL) {
           validate_inputs(initial, label)
           checkbox <- list(type = 2)
           initFields(l = checkbox, label = label, 
                      initial = initial)
           .self
         },
         make_gui = function(cont, handler, ...) {
           widget <<- gtkCheckButton() # no label
           widget$setActive(initial)
           gSignalConnect(widget, "toggled", handler)
           callSuper(cont)
         },
         get_value = function() widget['active']
         )


###################################################
### code chunk number 21: Checkbox_constructor
###################################################
checkbox <- function(initial = FALSE, label = NULL) Checkbox$new(initial, label)                            


###################################################
### code chunk number 22: ex-RGtk2-manipulate.Rnw:464-465
###################################################
manipulate(## expression
           plot(cars, xlim = c(x.min, x.max), type = type, 
                axes = axes, ann = label),
           ## controls
           x.min = slider(0, 15),
           x.max = slider(15, 30, initial = 25),
           type = picker("p", "l", "b", "c", "o", "h", "s"),
           axes = checkbox(TRUE, label = "Draw Axes"),
           label = checkbox(FALSE, label = "Draw Labels")
           )


