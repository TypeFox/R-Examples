### R code from vignette source 'ex-tcltk-subset-filter.Rnw'

###################################################
### code chunk number 1: ex-tcltk-subset-filter.Rnw:34-35
###################################################
library(tcltk)


###################################################
### code chunk number 2: FilterList
###################################################
setOldClass("tkwin")
setOldClass("tclVar")
FilterList <- setRefClass("FilterList",
                          fields = list(
                            DF = "data.frame",
                            l = "list",
                            id = "ANY",
                            frame = "tkwin" 
                            ))


###################################################
### code chunk number 3: FilterListMethods
###################################################
FilterList$methods(
          setup_gui = function(parent) {
            enc_frame <- ttkframe(parent, padding = 5)
            tkpack(enc_frame, expand = TRUE, fill = "both")
            frame <<- ttkframe(enc_frame)
            button_frame <- ttkframe(enc_frame)
            ## use grid to manage these
            tkgrid(frame, sticky = "news")
            tkgrid(button_frame, sticky = "new")
            tkgrid.rowconfigure(enc_frame, 1, weight = 1)
            tkgrid.columnconfigure(enc_frame, 0, weight = 1)
            ##
            add_button <- 
              ttkbutton(button_frame, text = "Add", 
                        command = function() .self$add())
            preview_button <- 
              ttkbutton(button_frame, text = "Preview", 
                        command = function() .self$preview())
            ##
            sapply(list(add_button, preview_button), tkpack, 
                   side = "left", padx = 5)
          })


###################################################
### code chunk number 4: FilterListInitialize
###################################################
FilterList$methods(
           initialize = function(DF, parent, ...) {
            initFields(DF = DF, l = list(), id = 0L)
            setup_gui(parent)
             callSuper(...)
           })


###################################################
### code chunk number 5: FilterListSelectVariable
###################################################
FilterList$methods(
           select_variable = function() {
             "Return a variable name from the data frame"
             x <- sapply(DF, function(i) class(i)[1])
             m <- cbind(Variables = names(x), Type = x)
             window <- tktoplevel()
             fr <- ttkframe(window, padding = c(3,3,3,12))
             tkpack(fr, expand = TRUE, fill = "both")
             ##
             a <- populate_rectangular_treeview(fr, m)
             tkconfigure(a$frame, width = 300, height = 200)
             tkpack(a$frame, expand = TRUE, fill = "both")
             ## select a value, store in out
             out <- NA
             tkbind(a$tr, "<<TreeviewSelect>>", function(W) {
               sel <- tcl(W, "selection")
               val <- tcl(W, "item", sel, "-values")
               assign("out", as.character(val)[1], 
                      inherits = TRUE)
               tkdestroy(window)
             })
             tkwait.window(window)
             return(out)
           })


###################################################
### code chunk number 6: FilterListAdd
###################################################
FilterList$methods(
           add = function(variable_name, ...) {
             if(missing(variable_name)) 
               variable_name <- select_variable()
             x <- get(variable_name, DF)
             ## new item
             id <<- id + 1
             item <- newFilterItem(x,variable_name, id, .self)
             ## make frame
             enc_frame <- ttkframe(frame)
             tkpack(enc_frame, 
                    expand = TRUE, fill = "both", pady = 2)
             l[[as.character(id)]] <<- list(frame = enc_frame, 
                                            item = item)
             item$make_gui(enc_frame)
           })


###################################################
### code chunk number 7: FilterListRemove
###################################################
FilterList$methods(
           remove=function(id_obj, ...) {
             "Remove. id is character or item object"
             if(!is.character(id_obj))
               id_obj <- id_obj$id
             tkpack.forget(l[[id_obj]]$frame)
             l[[id_obj]] <<- NULL
           })


###################################################
### code chunk number 8: FilterListGetValue
###################################################
FilterList$methods(
           get_value = function() {
             "Return logical value for all filter items"
             if(length(l) == 0)
               return(rep(TRUE, length=nrow(DF)))
             ##
             out <- sapply(l, function(i) i$item$get_value())
             out[is.na(out)] <- FALSE   ## coerce NA to FALSE
             apply(out, 1, all)
           })


###################################################
### code chunk number 9: ex-tcltk-subset-filter.Rnw:227-252
###################################################
FilterList$methods(
           preview = function() {
             "Preview data frame"
             ind <- get_value()
             if(!any(ind)) {
               message("No matches")
               return()
             }
             ## coerce to character
             m <- DF[ind,]
             for(i in seq_along(m)) 
               m[,i] <- as.character(m[,i])
             ##
             window <- tktoplevel()
             fr <- ttkframe(window, padding = c(3,3,3,12))
             tkpack(fr, expand = TRUE, fill = "both")
             a <- populate_rectangular_treeview(fr, m)
             tkconfigure(a$frame, width = 400, height = 300)
             tkpack(a$frame, expand = TRUE, fill = "both")
             ##
             button <- ttkbutton(fr, text = "dismiss", 
                         command=function() tkdestroy(window))
             tkpack(button, anchor = "sw")
             tkwait.window(window)
           })


###################################################
### code chunk number 10: runIt (eval = FALSE)
###################################################
## window <- tktoplevel()
## require(MASS)
## filter_list <- FilterList$new(DF = Cars93, parent = window)


###################################################
### code chunk number 11: newFilterItem
###################################################
newFilterItem <- function(x, nm = deparse(substitute(x)), id, 
                          list_ref) UseMethod("newFilterItem")
newFilterItem.default <- function(x,nm=deparse(substitute(x)), 
                                  id, list_ref) {
  FilterItemNumeric$new(x = x, nm = nm, id = id, 
                        list_ref = list_ref)
}


###################################################
### code chunk number 12: ex-tcltk-subset-filter.Rnw:280-288
###################################################
## not shown
newFilterItem.character <- function(x,  nm = deparse(substitute(x)), id, list_ref) {
  FilterItemCharacter$new(x = x, nm = nm, id = id, list_ref = list_ref)
}

newFilterItem.factor <- function(x,  nm = deparse(substitute(x)), id, list_ref) {
  newFilterItem(as.character(x), nm, id, list_ref)
}


###################################################
### code chunk number 13: FilterItem
###################################################
FilterItem <- setRefClass("FilterItem",
                          fields = list(
                            x = "ANY",
                            nm = "character",
                            id = "character",
                            list_ref = "ANY"
                            ))


###################################################
### code chunk number 14: FilterItemInitialize
###################################################
FilterItem$methods(
           initialize = function(...) {
             initFields(...)
             .self
           },
           get_value = function() {
             "Return logical value of length x"
             stop("Must be subclassed")
           },
           remove = function() list_ref$remove(.self),
           make_gui = function(parent, ...) {
             "Set up GUI, including defining widgets"
             remove_button <- ttkbutton(parent, text="remove",
                                    command = function() {
                                      .self$remove()
                                    })
             tkpack(remove_button, side = "right")
           })


###################################################
### code chunk number 15: FilterItemNumeric
###################################################
FilterItemNumeric <- setRefClass("FilterItemNumeric",
                                 contains = "FilterItem",
                                 fields = list(
                                   ineq_variable = "tclVar",
                                   value_variable = "tclVar"
                                   ))


###################################################
### code chunk number 16: FilterItemNumericGetValue
###################################################
FilterItemNumeric$methods(
      get_value = function() {
        xpr <- paste(nm, tclvalue(ineq_variable), 
                     tclvalue(value_variable))
        eval(parse(text = xpr), 
             envir = list_ref$DF, parent.frame())
      })



###################################################
### code chunk number 17: FilterItemNumericMakeGui
###################################################
FilterItemNumeric$methods(
      make_gui = function(parent) {
        ## standard width for label
        label_width <- max(sapply(names(list_ref$DF), nchar))
        label <- ttklabel(parent, text=nm, width=label_width)
        ## ineq combo
        vals <- c(">=", ">", "==", "!=", "<", "<=")
        ineq_variable <<- tclVar("<=")
        ineq <- ttkcombobox(parent, values = vals, 
                   textvariable = ineq_variable, width = 4)
        ## entry
        value_variable <<- tclVar(max(x, na.rm = TRUE))
        val <- ttkentry(parent, textvariable = value_variable)
        ##
        sapply(list(label, ineq, val), tkpack, side = "left",
               padx = 5)
        callSuper(parent)
      })



###################################################
### code chunk number 18: FilterItemCharacter
###################################################
FilterItemCharacter <- 
  setRefClass("FilterItemCharacter",
              contains = "FilterItem",
              fields = list(
                tr = "tkwin",
                button = "tkwin",
                poss_vals = "character",
                cur_vals = "character"
                ))


###################################################
### code chunk number 19: ex-tcltk-subset-filter.Rnw:443-447
###################################################
FilterItemCharacter$methods(
          get_value = function() {
            x %in% cur_vals
          })


###################################################
### code chunk number 20: sel_by_name
###################################################
sel_by_name <- function(tr, nms) {
  all_ind <- as.character(tcl(tr, "children", ""))
  vals <- sapply(all_ind, function(i) {
    as.character(tcl(tr, "item", i, "-values"))
  })
  ind <- names(vals[vals %in% nms])
  sapply(ind, function(i) tcl(tr, "selection", "add", i))
  sapply(setdiff(all_ind, ind), 
         function(i) tcl(tr, "selection", "remove", i))
}


###################################################
### code chunk number 21: FilterItemShortenPoss
###################################################
FilterItemCharacter$methods(ellipsize = function() {
            tmp <- paste(cur_vals, collapse = ", ")
            if((N <- nchar(tmp)) > 50)
              tmp <- sprintf("%s...%s", substr(tmp, 0, 15),
                             substr(tmp, N-12, N))
            sprintf("%50s", tmp)
          })


###################################################
### code chunk number 22: FilterItemCharacterSelectValuesDialog
###################################################
FilterItemCharacter$methods(
          select_values_dialog = function() {
            window <- tktoplevel()
            fr <- ttkframe(window, padding = c(3,3,12,12))
            tkpack(fr, expand = TRUE, fill = "both")
            tkpack(ttklabel(fr, 
              text = "Select values by extending selection"))
            ## selection
            m <- matrix(poss_vals)
            colnames(m) <- "Values"
            a <- populate_rectangular_treeview(fr, m)
            tkconfigure(a$tr, selectmode = "extended")
            tkconfigure(a$frame, width = 200, height = 300)
            tkpack(a$frame, expand = TRUE, fill = "both")
            
            sel_by_name(a$tr, cur_vals)         # see above
            
            tkbind(a$tr, "<<TreeviewSelect>>", function() {
              ind <- as.character(tcl(a$tr, "selection"))
              cur <- sapply(ind, function(i) {
                as.character(tcl(a$tr, "item", i, "-values"))
              })
              if(length(cur) == 0)
                cur <- character(0)
              cur_vals <<- cur
            })
            ## buttons
            frame_1 <- ttkframe(fr)
            tkpack(frame_1)
            toggle_button <- ttkbutton(frame_1, text="toggle", 
                         command=function() toggle_sel(a$tr))
            set_button <- ttkbutton(frame_1, text = "set", 
                         command=function() tkdestroy(window))
            sapply(list(toggle_button, set_button), tkpack, 
                   side = "left", padx = 5)
            ## make modal
            tkwait.window(window)
            tkconfigure(button, text = ellipsize())
          })


###################################################
### code chunk number 23: FilterItemCharacterMakeGUI
###################################################
FilterItemCharacter$methods(make_gui = function(parent) {
            poss_vals <<- sort(unique(as.character(x)))
            cur_vals <<- poss_vals
            ## label, ineq, val
            l_width <- max(sapply(names(list_ref$DF), nchar))
            label <- ttklabel(parent, text = nm, 
                              width = l_width)
            ##
            in_label <- ttklabel(parent, text = "%in%")
            ##
            button <<- ttkbutton(parent, text = ellipsize(), 
                       command = .self$select_values_dialog)
            ##
            sapply(list(label, in_label), tkpack,
                   side = "left", padx = 5)
            tkpack(button, 
                   expand = TRUE, fill = "x", side = "left")
            callSuper(parent)
          })


###################################################
### code chunk number 24: Helpers
###################################################
## helpers
## We don't show these
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}

populate_rectangular_treeview <- function(parent, m) {
  enc_frame <- ttkframe(parent)
  frame <- ttkframe(enc_frame)
  tkpack(frame, expand = TRUE, fill = "both")
  tr <- ttktreeview(frame,
                    columns = seq_len(ncol(m)),
                    show = "headings",
                    selectmode = "browse"
                    )
  addScrollbars(frame, tr)
  tkpack.propagate(enc_frame, FALSE)


  ## headings,widths
  charWidth <- as.integer(tclvalue(tcl("font","measure","TkTextFont","0")))
  sapply(seq_len(ncol(m)), function(i) {
    tcl(tr, "heading", i, text = colnames(m)[i])
    tcl(tr, "column", i, width = 10 + charWidth*max(apply(m, 2, nchar)))
  })
  tcl(tr, "column", ncol(m), stretch = TRUE)
  
  ## values
  if(ncol(m) == 1)  m <- as.matrix(paste("{",m ,"}", sep=""))
  apply(m, 1, function(vals) 
    tcl(tr, "insert", "", "end", values = vals)
        )
  return(list(tr = tr, frame = enc_frame))
}

   

cur_sel <- function(tr) {
  ind <- as.character(tcl(tr, "selection"))
  sapply(ind, function(i) {
    as.character(tcl(tr, "item", i, "-values"))
  })
}

toggle_sel <- function(tr) {
  children <- as.character(tcl(tr, "children", ""))
  tcl(tr, "selection", "toggle", children) 
}


###################################################
### code chunk number 25: ex-tcltk-subset-filter.Rnw:618-622
###################################################
## Call when all is said and done
window <- tktoplevel()
require(MASS)
filter_list <- FilterList$new(DF = Cars93, parent = window)


