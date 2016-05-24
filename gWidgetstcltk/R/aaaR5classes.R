## A widget framework for tcltk. Mostly just for fun, but useful for 3 widgets so far.
## This provides some R reference classes for working with tcltk widgets. A few new ones are used within gWidgetstcltk
## Namely: gradio and gcheckboxgroup, as we can now extend and lengthen the number of options
##         and gedit for the autocompletion code


##################################################
## Class structures
setRefClass("TcltkWidget",
            fields=list(
              widget="ANY",                  # main widget
              value = "ANY",                 # main value of widget
              ## handler stuff
              handlers = "list",             # list of handlers
              handler_args = "character",   # list of substitution values
              handler_id = "numeric",        # count of ids
              block_all_handlers= "logical", # flag to block all handlers
              blocked_handlers = "character" # all blocked handlers
              ),
            methods=list(
              init_widget = function(parent, ...) {
                "Initialize widget"
                NULL
              },
              initialize=function(parent, ...) {
                "Set up widget any other necessary values. Should be subclassed."
                handlers <<- list()
                handler_args <<-  character()
                handler_id <<- 0
                block_all_handlers <<- FALSE
                ## now initialize widget
                init_widget(parent, ...)
                .self
              },
              get_widget=function() {
                "Return widget"
                widget
              },
              set_widget=function(widget) {
                "Set widget"
                widget <<- widget
              },
              ##################################################
              ## handler methods
              new_handler_id = function(signal) {
                "Create a new handler id. An id is a list with a unique id and the signal"
                handler_id <<- handler_id + 1
                list(signal=signal, id=as.character(handler_id))
              },
              ##
              ##
              run_handlers = function(h) {
                ##' run handlers
                ##' @param h list passing in signal, user.data percent subs. values
                ##' @note much hackery to environment of handler to get percent substitutions passed in
                ##' also pass in reference to \code{.self} and \code{user.data}
                "Run handlers for a given signal"
                if(block_all_handlers || !is_enabled())
                  return()

                lapply(handlers[[h$signal]], function(i) {
                  if(!i$id %in% blocked_handlers) {
                    FUN <- i$handler
                    e <- environment(FUN)
                    for(j in names(formals(FUN))) e[[j]] <- h[[j]]
                    if(!exists(".self", e))
                      e[[".self"]] <- .self
                    e[["user.data"]] <- i$user.data
                    environment(FUN) <- e
                    formals(FUN) <- alist()
                    FUN()
                  }
                })
              },
              ##
              ##
              bind_handler = function(signal) {
                ##' @param signal signal o find handlers for
                ##' @param return id
                ##' @note hackery to get substitution variables passed along
                "Bind run_handlers for give signal."
                if(signal == "command") 
                  tkconfigure(widget, command=make_f(signal))
                else
                  tkbind(widget, signal, make_f(signal))
              },
              ##
              ##
              make_f = function(signal) {
                ##' Make a function for the handler. Involves hackery to get signature correct
                f <- function() {
                  h <- list()
                  for(i in handler_args) h[[i]] <- if(exists(i)) get(i) else NULL
                  h[['signal']] <- signal
                  run_handlers(h)
                }
                if(length(handler_args)) {
                  ## do an eval/parse hack. Not sure how else to work with alist.
                  txt <- paste("alist","(",
                               paste(c(handler_args), "=", sep="", collapse=","),
                               ")", sep="")
                  formals(f) <- eval(parse(text=txt))
                }
                f
              },
              add_handler = function(signal, handler, user.data=NULL) {
                ##' add a handler for a given signal
                ##' @param signal signal for tkbind
                ##' @param handler function. Percent substitution values used. Also user.data passed in
                ##' if given as an argument
                ##' @param user.data will be passed into handler if it includes \code{user.data} argument
                ##' @note much hackery to get percent substitutions passed to handlers. Handler has
                ##' \code{.self} reference and \code{user.data} if argument present
                ##' @return returns an id
                "Add a handler to the list"
                l <- handlers
                if(is.null(l[[signal]]))
                  l[[signal]] <- list()
                id <- new_handler_id(signal)
                l[[signal]][[id$id]] <- list(handler=handler, user.data=user.data, id=id$id)
                handlers <<- l
                handler_args <<- unique(c(handler_args, setdiff(names(formals(handler)), "user.data")))
                bind_handler(signal)
                id
              },
              ##
              ##
              remove_handler = function(id) {
                ##' remove a handler by id
                ##' @param id an id returned by add_handler method
                "Remove a handler by its id"
                l <- handlers[[id$signal]]
                l[[id$id]] <- NULL
                handlers[[id$signal]] <<- l
                if(length(l) == 0)
                  tkbind(widget, id$signal, "")
              },
              ##
              ##
              block_handler = function(id) {
                ##' block a handler by id (or all handlers)
                ##' @param id an id to identify which handler to block. If missing all handlers will be blocked
                "block a handler by its id or all handlers"
                if(missing(id))
                  block_all_handlers <<- TRUE
                else
                  blocked_handlers <<- unique(c(blocked_handlers, id$id))
              },
              ##
              ##
              unblock_handler = function(id) {
                ##' unblock an handler by id (or all handlers)
                ##' @param id an id to identify which handler to unblock. If missing all handlers are unblocked
                "Unblock handler by id, or all handlers"
                if(missing(id)) {
                  block_all_handlers <<- FALSE
                  blocked_handlers <<- character(0)
                } else {
                  blocked_handlers <<- setdiff(blocked_handlers, id$id)
                }
              },
              ##################################################
              ## API for widgets. These are state properties
              get_value = function() {
                ##' Return primary value of widget. Possibly coerced in subclasses
                ##' @return value
                "Return value of widget"
                value
              },
              set_value = function(value, index=TRUE) {
                ##' set the main value for the widget
                ##' @param value value to be set
                ##' @note subclass should update the widget, this just updates value field
                "Set value for widget. Updates value field. Subclass should update widget's value"
                value <<- value
              },
              is_enabled = function() {
                ##' Is the widget enabled (or disabled)
                "Is widget enabled (as opposed to disabled)"
                as.logical(tcl(widget, "instate", "!disabled"))
              },
              set_enabled = function(value) {
                ##' set the enabled status of widget
                ##' @param value if TRUE enable widget, if FALSE disable (not sensitive to user input)
                "Set enabled or disabled as per value (logical)"
                tcl(widget, "state", ifelse(value, "!disabled", "disabled"))
                invisible()
              },
              is_focus = function() {
                ##' does the widget have the focus
                ##' @return logical
                "Does widget have focus?"
                as.logical(tcl(widget, "instate", "focus"))
              },
              set_focus = function(value=TRUE) {
                ##' set focus
                ##' @param value if TRUE set focus on widget
                if(value)
                  tcl(widget, "state", "focus")
                invisible()
              },
              is_readonly = function() {
                ##' is widget readonly?
                "Is widget read only?"
                as.logical(tcl(widget, "instate", "readonly"))
              },
              set_readonly = function(value=TRUE) {
                ##' set widget as readonly (not editable)
                ##' @note readonly mya vary for some widgets
                "Set widget as readonly. That is not editable."
                if(value)
                  tcl(widget, "state", "readonly")
                invisible()
              }
              
              )
            )

setRefClass("Button",
            contains=c("TcltkWidget"),
            methods= list(
              ## subclass overrides
              init_widget = function(parent, text="", ...) {
                widget <<- ttkbutton(parent, ...)
                set_value(text)
              },
              set_value = function(value) {
                tkconfigure(widget, text=value)
                callSuper(value)
              },
              ## new API
              set_image = function(image, compound="left") {
                if(!tclObj_exists(image)) {
                  ## a stock icon
                  image <- "XXX get icon name from stock icons"
                }
                tkconfigure(widget, "image"=image, compound=compound)
              }
                
              
              )
            )

setRefClass("TcltkWidgetWithTclvariable",
            contains=c("TcltkWidget"),
            fields=list(
              v = "tclVar",
              coerce_with="function"
              ),
            methods=list(
              initialize=function(...) {
                coerce_with <<- as.character
                callSuper(...)
              },
              set_coerce_with = function(f) {
                "Set function to coerce value with function(value) {...}"
                if(is.character(f))
                  f <- get(f, inherits=TRUE)
                coerce_with <<- f
              },
              get_value=function() {
                coerce_with(tclvalue(v))
              },
              set_value=function(value) {
                "Set value"
                a <- v                    # avoid local assignment warning
                tclvalue(a) <- value
              }
              )
            )

## simple label
setRefClass("Label",
            contains=c("TcltkWidgetWithTclvariable"),
            methods=list(
              init_widget = function(parent, text="") {
                ##' @param parent parent widget
                ##' @param text text for label or tclVariable 
                if(is.character(text))
                  v <<- tclVar(text)
                else
                  v <<- text
                widget <<- ttklabel(parent, textvariable=v)
              },
              set_value = function(value) {
                .value <- paste(value, collapse="\n")
                callSuper(.value)
              }
              )
            )

## A check button
setRefClass("CheckButton",
            contains=c("TcltkWidgetWithTclvariable"),
            methods=list(
              init_widget=function(parent, text, checked=FALSE, image, compound="none") {
                v <<- tclVar(as.numeric(checked))
                widget <<- ttkcheckbutton(parent, variable=v)
                set_label(text)
                if(!missing(image))
                  set_image(image, compound)
                set_coerce_with(function(x) as.logical(as.numeric(x)))
              },
              set_label=function(value) {
                tkconfigure(widget, text=value)
              },
              set_image=function(image, compound) {
#                if(tclObj_exists(image))
                  tkconfigure(widget, image=image, compound=compound)
              } ##,
              ## make_f = function(signal) {
              ##   ##' Make a function for the handler. Involves hackery to get signature correct
              ##   f <- function() {
              ##     h <- list()
              ##     for(i in handler_args) h[[i]] <- get(i)
              ##     h[['signal']] <- signal
              ##     ##tcl("after", 150, function(...) {
              ##       run_handlers(h)
              ##     ##})
              ##   }
              ##   if(length(handler_args)) {
              ##     ## do an eval/parse hack. Not sure how else to work with alist.
              ##     txt <- paste("alist","(",
              ##                  paste(c(handler_args), "=", sep="", collapse=","),
              ##                  ")", sep="")
              ##     formals(f) <- eval(parse(text=txt))
              ##   }
              ##   f
              ## }

              ),
            )


## entry with type ahead
##' configuration property. (tkconfigure(widget, foreground="gray")??
setRefClass("Entry",
            contains=c("TcltkWidgetWithTclvariable"),
            fields=list(
              m="tkwin",
              l="tkwin",
              lindex = "numeric",  # index of selection widget
              no.wds = "numeric",  # track number of possible wds to choose from
              words = "character",
              max.words = "numeric", # maximum words in a display
              init_msg = "character" # an initial message
              ),
              
            methods=list(
              init_widget = function(parent, text="", coerce.with, max.words=20, words) {
                ##' @param parent parent widget
                ##' @param text text for label or tclVariable 
                if(is.character(text))
                  v <<- tclVar(text)
                else
                  v <<- text
                widget <<- ttkentry(parent, textvariable=v)
                if(!missing(coerce.with))
                  set_coerce_with(coerce.with)
                ## popup stuff
                tclServiceMode(FALSE)                                      
                m <<- tktoplevel()
                tkwm.transient(m, parent)
                tkwm.overrideredirect(m, TRUE)
                tkwm.withdraw(m)
                tclServiceMode(TRUE)                      
                l <<- tktext(m); tkpack(l)
                lindex <<- 0      # index of selected
                max.words <<- max.words
                if(!missing(words))
                        set_words(words)
                ##
                set_init_msg("")
                addBindings()
              },
              set_words = function(words) {
                words <<- unique(as.character(words))
              },
              set_value = function(value) {
                old_value <- tclvalue(v)
                if(old_value == init_msg)
                    tkconfigure(widget, foreground="black")
                v_local <- v
                tclvalue(v_local) <- value
                lindex <<- 0
                tcl(widget, "icursor", "end")
                if(old_value != tclvalue(v))
                  tcl("event","generate", widget, "<<Changed>>")
                callSuper(value)
              },
              ## find match in word list
              findMatch = function(x) {
                ind <- grepl(sprintf("^%s", tolower(x)), tolower(words))
                words[ind]
              },
              showWordList = function(str) {
                ##' show the word list
                ##' @param str a string. If
                ##' missing do nothing, otherwise match against
                ##' string to generate word list. Popup menu
                ##' depending on length
                
                char.height <- 16 ## or compute from font metrics
                wds <- findMatch(str)
                if(length(wds) == 0) {
                  no.wds <<- 0
                  hideWordList()
                  return()
                }
                
                ## compute max.height -- number of words that can be shown
                screenheight <- as.numeric(tkwinfo("screenheight", widget))
                y <- as.numeric(tclvalue(tkwinfo("rooty",widget)))
                max_words <- min(max.words, floor((screenheight - y)/char.height))
                if(length(wds) > max_words)
                  wds <- c(wds[1:max_words], "...")
                tkdelete(l, "0.0", "end")
                tkinsert(l, "end", paste(wds, collapse="\n"))
                lindex <<- 1; no.wds <<- length(wds)
                
                ## set geometry
                x <- as.numeric(tclvalue(tkwinfo("rootx", widget)))
                y <- as.numeric(tclvalue(tkwinfo("rooty",widget)))
                geo <- as.character(tkwinfo("geometry",widget))
                geo <- as.numeric(strsplit(geo, "[x+]")[[1]])
                tkwm.geometry(m, sprintf("%sx%s+%s+%s", geo[1], 10 + char.height*length(wds), x, y + geo[2]))
                ## popup
                tcl("wm","attributes", m, "topmost"=TRUE) 
                tcl("wm","attributes", m, "alpha"=0.8)
                tkwm.deiconify(m)
                tcl("raise", m)
                highlightWordList()
              },
              ## hide the word list
              hideWordList = function() {
                tcl("wm","attributes", m, "topmost"=FALSE) # not working!
                tkwm.withdraw(m)
              },
              ## highlight word on lindex
              highlightWordList = function() {
                if(lindex > 0) {
                  tktag.remove(l, "selectedWord", "0.0", "end")
                         tktag.add(l,"selectedWord",sprintf("%s.0", lindex), sprintf("%s.end", lindex))
                  tktag.configure(l, "selectedWord", font="bold")
                }
              },
              no_items = function() {
                length(no.wds)
              },
              ## get current word. From lineindex if applicable, or from entry widget itself
              getCurrentWord = function() {
                if(no.wds > 0)
                  if(lindex > 0) {
                    tclvalue(tkget(l, sprintf("%s.0", lindex), sprintf("%s.end", lindex)))
                  } else {
                    ""
                  }
                else
                  tclvalue(v)
              },
              ##' initial message code
              get_value = function() {
                "Get the text value"
                if(!is_init_msg())
                  as.character(tclvalue(v))
                else
                  ""
              },
              set_text = function(text, hide=TRUE) {
                "Set text into widget"
                if(hide)
                  hide_init_msg()
                set_value(text)
              },
              set_init_msg=function(txt) {
                "Set the initial message"
                init_msg <<- txt
              },
              is_init_msg=function() {
                "Is the init text showing?"
                if(nchar(init_msg) == 0)
                  FALSE
                else
                  as.character(tclvalue(v)) == init_msg
              },
              hide_init_msg= function() {
                "Hide the initial text"
                if(is_init_msg()) {
                  tkconfigure(widget, foreground="black")
                  set_text("", hide=FALSE)
                }
              },
              show_init_msg=function() {
                "Show the intial text"
                tkconfigure(widget, foreground="gray")
                set_text(init_msg, hide=FALSE)
              },
              ##' Add bindings to entry box
              addBindings = function() {
                add_handler("<KeyRelease>", function(W, K) {
                  ## set out virtual event, as otherwise we can;t have addHandlerKeystroke
                  tcl("event","generate", .self$widget, "<<KeyRelease>>", "data"=K) 
                  ## Main bindings
                  if(nchar(K) == 1 || K == "BackSpace") {
                    ## single letter, popup menu
                    val <- tclvalue(tcl(W, "get"))
                    showWordList(val)
                  } else if(K == "Down") {
                    ## down arrow. Open if empty, but also scroll down list
                    if(nchar(val <- getCurrentWord()) == 0) {
                      showWordList(".")
                      lindex <<- 0
                    }
                    lindex <<- min(lindex + 1, no.wds)
                    highlightWordList()
                  } else if(K == "Up") {
                    ## move up list
                    lindex <<- max(lindex - 1, 1)
                    highlightWordList()
                  } else if(K == "Return") {
                    ## get value and put into e widget
                    hideWordList()
                    if(lindex > 0) {
                      set_value(getCurrentWord())
                    } else {
                      tcl("event","generate", .self$widget, "<<Changed>>")
                    }
                  } else if(K == "Escape") {
                    ## close the word list
                    hideWordList()
                    lindex <<- 0
                  }
                })
                ## show or hide, depending
                add_handler("<Map>", showWordList)
                tkbind(tcl("winfo", "toplevel", widget), "<Configure>", hideWordList)
                tkbind(widget,"<Destroy>", hideWordList)
                add_handler("<FocusIn>", hide_init_msg)
                add_handler("<FocusOut>", hideWordList)
                add_handler("<FocusOut>", function() {
                  if(nchar(get_value()) == 0 && nchar(init_msg) > 0)
                    show_init_msg()
                })
                add_handler("<Unmap>", hideWordList)

                tkbind(l, "<Motion>", function(x, y) {
                  tmp <- as.character(tcl(l, "index", sprintf("@%s,%s", x, y)))
                  lindex <<- as.numeric(strsplit(tmp, "\\.")[[1]][1])
                  highlightWordList()
                })
                      
                ## bind to text widget
                tkbind(l, "<Button-1>", function(x,y) {
                  wd <- getCurrentWord()
                  hideWordList()
                  if(wd != "...") {
                    set_value(getCurrentWord())
                  }
                })
                ## we don't want focus on l
                tkbind(l, "<FocusIn>", function() {
                  tkfocus(widget)
                })
              }
              
              )
            )

## A class for widget with items (radio buttons, checkbuttons)
setRefClass("TcltkWidgetWithItems",
            contains=c("TcltkWidget"),
            fields=list(
              items = "ANY",            # hold items
              button_items = "list"     # item instances
              ),
            methods=list(
              ### Handlers are important
              bind_handler=function(signal) {
                ##' override. We bind to button items
              },
              run_handlers=function() {
                ##' override. Handlers run by button item
              },
              add_handler=function(signal, handler, user.data=NULL) {
                ##' add handler to each item
                ##' @return id of handler. A list.
                ## stash handler into list
                callSuper(signal, handler, user.data)
                ## add handler to each button item
                id <- lapply(button_items, function(i) {
                  list(widget=i, id=i$add_handler(signal, handler, user.data))
                })
                return(id)
              },
              remove_handler=function(id) {
                ##' @return NULL
                if(missing(id)) {
                  lapply(button_items, function(i) i$remove_handler())
                } else {
                  lapply(id, function(i) {
                    (i$widget)$remove_handler(i$id)
                  })
                }
                invisible()
              },
              block_handler=function(id) {
                ##' @return NULL
                if(missing(id)) {
                  lapply(button_items, function(i) i$block_handler())
                } else {
                  lapply(id, function(i) {
                    (i$widget)$block_handler(i$id)
                  })
                }
                invisible()
              },
              unblock_handler=function(id) {
                ##' @return NULL
                if(missing(id)) {
                  lapply(button_items, function(i) i$unblock_handler())
                } else {
                  lapply(id, function(i) {
                    (i$widget)$unblock_handler(i$id)
                  })
                }
                invisible()                  
              },
              transfer_handlers=function() {
                ##' @return NULL
                "Copy handlers onto child items"
                lapply(button_items, function(i) {
                  i$handlers <- handlers
                  i$handler_args <- handler_args
                  for(signal in names(handlers)) # signal
                    i$bind_handler(signal) 
                })
                invisible()
              },
              ##################################################
              get_items = function(drop=TRUE) {
                ##' @param drop if TRUE just items, else items and images as a data frame (if present)
                "Get items to select from. Drops images by default."
                if(drop)
                  items[,1, drop=TRUE]
                else
                  items                 # a data frame
              },
              no_items = function() {
                "Number of items"
                length(get_items())
              },
              make_new_item = function(i, has_image, compound) {
                ##' @param i index from item list
                ##' @param has_image logical indicating if an image is specified
                ##' @param compound if an image, how is it shown
                "Make a new item. Called from set_items"
                ## OVerride in subclasses
              },
              set_items = function(new_items, compound) {
                ##' @param new_items a vector of data frame (character, images)
                ##' @param compound if images, then specifies how.
                "Set new items for object"
                if(!is.data.frame(new_items)) {
                  new_items <- data.frame(new_items, stringsAsFactors=FALSE)
                  items <<- new_items
                  has_image <- FALSE
                } else {
                  items <<- new_items
                  has_image <- TRUE
                }
                ## clear any children
                tclServiceMode(FALSE)
                lapply(as.character(tkwinfo("children", widget)), function(window_id) {
                  tkpack.forget(window_id)
                })

                ## add in button items
                button_items <<-
                  lapply(seq_len(no_items()), function(i) {
                    new_item <- make_new_item(i, has_image, compound)
                    tkpack(new_item$get_widget(), side=orientation, anchor="nw")
                    new_item
                  })
                transfer_handlers()
                tclServiceMode(TRUE)
              },
              ##
              set_enabled=function(value) {
                callSuper(value)
                lapply(button_items, function(i) i$set_enabled(value))
                invisible()
              }
                     
              

              )
            )
## Radio button group
setRefClass("RadioButton",
            contains=c("TcltkWidgetWithItems"),
            fields=list(
              orientation="character",  # which orientation for packing
              state_variable="tclVar"   # holds state
              ),
            methods=list(
              init_widget = function(parent, items, selected=1, horizontal=TRUE, compound="none") {
                ##' @param items a vector of items or a data frame with columns items and images (names)
                ##' @param compound if images specified, how to configure
                
                widget <<- ttkframe(parent)
                selected <- max(1, min(as.integer(selected), length(items)))
                orientation <<- ifelse(horizontal, "left", "top") # pack arguments for side
                state_variable <<- tclVar(ifelse(is.null(dim(items)), items[selected], items[selected,1]))
                set_items(items, compound)
              },
              get_index = function() {
                "Get selected value by index"
                as.integer(which(get_items() %in% as.character(get_value())))
              },
              set_index = function(i) {
                "Set selected value by index"
                i <- as.integer(i)
                tmp <- get_items()
                if(i < 1 || i > length(tmp))
                  i <- 1                # default
                set_value(tmp[i])
              },
              get_value = function(index=FALSE) {
                "Get selected value"
                if(index) return(get_index())
                as.character(tclvalue(state_variable))
              },
              set_value = function(value, index=FALSE) {
                "Set selected value by label"
                if(index) return(set_index(value))
                
                if(value %in% get_items()) {
                  state_variable_local <- state_variable
                  tclvalue(state_variable_local) <- value
                }
                if(is_enabled())
                  tcl(button_items[[get_index()]]$get_widget(), "invoke")
              },
              make_new_item=function(i, has_image, compound) {
                "Make a new item. Pass in image information"
                if(!has_image)
                  new_item <- getRefClass("RadioButtonItem")$new(parent=widget,
                                                                 state_variable=state_variable,
                                                                 text=items[i,1])
                else
                  new_item <- getRefClass("RadioButtonItem")$new(parent=widget,
                                                                 state_variable=state_variable,
                                                                 text=items[i,1],
                                                                 image=items[i,2],
                                                                 compound)
                new_item
              },
              set_items = function(new_items, compound="none") {
                ##' @param new_items character vector or data frame
                ##' @param compound
                "Set items for radio button group, configure handlers, set state, and add to frame"
                ##
                
                ## store selected in case we are replacing
                if(!is(items, "uninitializedField"))  {
                  selected <- get_index()
                } else {
                  selected <- NULL
                }

                callSuper(new_items, compound)
                
                if(!is.null(selected))
                  set_index(selected)

              }
              ))

setRefClass("RadioButtonItem",
            contains=c("TcltkWidget"),
            fields=list(
              state_variable="tclVar"
              ),
            methods=list(
              init_widget = function(parent, state_variable, text, image, compound="none") {
                ##' @param parent parent container
                ##' @param state_variable tclvariable holding the state
                ##' @param text value for the widget
                ##' @param image optional image
                ##' @param compound how to configure image if present
                widget <<- ttkradiobutton(parent, variable=state_variable)
                state_variable <<- state_variable
                set_label(text)
                if(!missing(image) && tclObj_exists(image)) 
                  set_image(image, compound)
              },
              is_checked=function() {
                "Is this item the checked one?"
                as.character(tclvalue(state_variable)) == get_label()
              },
              set_checked=function() {
                "Set this item as the checked one"
                state_variable_local <- state_variable
                tclvalue(state_variable_local) <- get_label()
              },
              get_label=function() {
                "Get label text"
                as.character(tkcget(widget, "-text"))
              },
              set_label=function(text) {
                "Set label text"
                tkconfigure(widget, value=text, text=text)
              },
              set_image=function(image, compound="none") {
                "Configure image for radio button"
                tkconfigure(widget, image=image, compound=compound)
              } ##,
              ## make_f = function(signal) {
              ##   ##' Make a function for the handler. Involves hackery to get signature correct
              ##   f <- function() {
              ##     h <- list()
              ##     for(i in handler_args) h[[i]] <- get(i)
              ##     h[['signal']] <- signal
              ##     ##tcl("after", 150, function(...) {
              ##       if(.self$is_checked())
              ##         run_handlers(h)
              ##     ##})
              ##   }
              ##   if(length(handler_args)) {
              ##     ## do an eval/parse hack. Not sure how else to work with alist.
              ##     txt <- paste("alist","(",
              ##                  paste(c(handler_args), "=", sep="", collapse=","),
              ##                  ")", sep="")
              ##     formals(f) <- eval(parse(text=txt))
              ##   }
              ##   f
              ## }

              )
            )


## checkbuttongroup
setRefClass("CheckButtonGroup",
            contains=c("TcltkWidgetWithItems"),
            fields=list(
              orientation="character"
              ),
            methods=list(
              init_widget = function(parent, items, selected=FALSE, horizontal=TRUE, compound="none") {
                ##' @param items a vector of items or a data frame with columns items and images (names)
                ##' @param compound if images specified, how to configure
                widget <<- ttkframe(parent)
                orientation <<- ifelse(horizontal, "left", "top") # pack arguments for side
                set_items(items, compound)
                set_value(selected)
              },
              make_new_item = function(i, has_image, compound) {
                "Make a new item. Pass in image information"
                if(!has_image)
                  new_item <- getRefClass("CheckButton")$new(parent=widget,
                                                                      checked=FALSE,
                                                                      text=items[i,1])
                else
                  new_item <- getRefClass("CheckButton")$new(parent=widget,
                                                             checked=FALSE,                                                                      text=items[i,1],
                                                                   image=items[i,2],
                                                                   compound)
                new_item
              },
              get_value=function(index=FALSE) {
                ind <- get_index()
                if(index)
                  return(ind)
                get_items()[ind]
              },
              set_value=function(value, index=FALSE) {
                ##' @param value vector of values from items
                if(index) return(set_index(value))
                
                
                if(is.logical(value)) {
                  .value <- rep(value, length.out=no_items())
                  ind <- which(.value)
                } else {
                  ind <- which(get_items() %in% value)
                }
                set_index(ind)
              },
              get_index=function() {
                ## return indices of logical
                which(sapply(button_items, function(i) i$get_value()))
              },
              set_index=function(ind) {
                lapply(seq_len(no_items()), function(i) {
                  if(i %in% ind)
                    ## invoke calls command **and** changes state so we set_value after
                    tcl(button_items[[i]]$get_widget(), "invoke")
                  button_items[[i]]$set_value(i %in% ind)
                })
                invisible()
              }
              )
            )


##' spinbutton class. Spinbutton is not a themed widget!
setRefClass("SpinButton",
            contains=c("TcltkWidget"),
            methods= list(
              ## subclass overrides
              init_widget = function(parent, from=0, to=100, by=1, ...) {
                ## ttk spinbox new as of 8.5.9 
                out <- try(tkwidget(parent, "ttk::spinbox", from=from, to=to, increment=by), silent=TRUE)
                if(inherits(out, "try-error"))
                  out <- tkwidget(parent, "spinbox", from=from, to=to, increment=by)
                widget <<- out
              },
              get_value = function() {
                "Return specified value"
                as.numeric(tcl(widget,"get"))
              },
              set_value = function(value, index=FALSE, notify=FALSE) {
                "set spinner values"
                tcl(widget,"set", as.numeric(value))
                if(notify) {
                  message("notify handler")
                }
              },
              set_items = function(items) {
                "Set items to select from. A regular sequence"
                ## check that value is a regular sequence
                if(length(items) <=1) {
                  warning("Can only assign a vector with equal steps, as produced by seq")
                  return(obj)
                }
                if(length(items) > 2 &&
                   !all.equal(diff(diff(items)), rep(0, length(items) - 2))) {
                  warning("Can only assign a vector with equal steps, as produced by seq")
                  return(obj)
                }
                
                ## get current value, increment
                curValue <- get_value()
                inc <- head(diff(items), n=1)

                tkconfigure(widget, from=min(items), to=max(items), increment=inc)
                tcl(widget, "set", curValue)

              },
              ## override the default for this, spinbox is old widget
              is_enabled = function() {
                as.character(tkcget(widget, "-state")) == "normal"
              }
              
              )
            )

## This gets folded into the gComponentR5tcltk class, so we cna define generic methods for these:


### methods for gComponentR5tcltk. Most are shared, but we have them hard coded.
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gComponentR5tcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            
            r5_widget <- obj@R5widget
            index <- getWithDefault(index, FALSE)
            if(index) {
              return(r5_widget$get_value(index=TRUE))
            } else {
              val <- r5_widget$get_value()
              if(.hasSlot(obj, "coercewith") && !is.null(obj@coercewith))
                return(obj@coercewith(val))
              else
                return(val)
            }
          })
          

## toggles state to be T or F
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gComponentR5tcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   r5_widget <- obj@R5widget
                   index <- getWithDefault(index, FALSE)

                   r5_widget$set_value(value, index=index)
                   return(obj)
                 })

setMethod("[",
          signature(x="gComponentR5tcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })

setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gComponentR5tcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            r5_widget <- x@R5widget
            items <- r5_widget$get_items()
            if(missing(i))
              items
            else
              items[i]
          })

setReplaceMethod("[",
                 signature(x="gComponentR5tcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })
setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gComponentR5tcltk"),
          function(x, toolkit, i, j, ..., value) {

            r5_widget <- x@R5widget
            if(!missing(i)) {
              items <- r5_widget$get_items()
              items[i] <- value
              value <- items
            }
            r5_widget$set_items(value)

             return(x)
          })

setMethod(".length",
          signature(toolkit="guiWidgetsToolkittcltk",x="gComponentR5tcltk"),
          function(x,toolkit) {
            r5_widget <- x@R5widget
            r5_widget$no_items()
          })


## inherited enabled isn't workgin                
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gComponentR5tcltk"),
                 function(obj, toolkit, ..., value) {

                   r5_widget <- obj@R5widget                   
                   r5_widget$set_enabled(value)
                   return(obj)
                  
                 })


## ##################################################
## w <- tktoplevel()
## b <- getRefClass("Button")$new(parent=w, text="boom chica boom")
## tkpack(b$get_widget())

## l <- getRefClass("Label")$new(parent=w, text="huh")
## tkpack(l$get_widget())


## e <- getRefClass("Entry")$new(parent=w, text="", coerce.with=as.character)
## e$set_words(state.name)
## tkpack(e$get_widget())

## id <- e$add_handler("<<Changed>>", handler=function(W, user.data) {
##   print(user.data)
##   print(.self$get_value())
## }, user.data="x")

## rb <- getRefClass("RadioButton")$new(parent=w, items=state.name[1:3], horizontal=TRUE)
## tkpack(rb$get_widget())
## rb$add_handler("command", handler=function(user.data) {
##   print(user.data$get_value())
## }, user.data=rb)

## cbg <- getRefClass("CheckButtonGroup")$new(parent=w, horizontal=FALSE, items=state.name[1:3])
## tkpack(cbg$get_widget())

## id <- cbg$add_handler(command", handler=function(user.data) {
##   print(user.data$get_value())
## }, user.data=cbg)



## f <- "~/Downloads/dumb.gif"
## nm <- make_tcl_image("dump",f)
## cb <- getRefClass("CheckButton")$new(parent=w, text="test", image=nm, compound="left")
## tkpack(cb$get_widget())


## items <- data.frame(a=state.name[1:3], b=c(nm,nm,nm), stringsAsFactors=FALSE)

## rb1 <- getRefClass("CheckButtonGroup")$new(parent=w, items=items, horizontal=TRUE, compound="left")
## tkpack(rb1$get_widget())

