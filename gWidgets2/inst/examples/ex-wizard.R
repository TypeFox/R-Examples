## an example of a simple wizard. Here we implement classes for the wizard and a wizard page interface. 

library(gWidgets2)

##' Sample wizard page interface. See below for two subclasses used to make a given page
WizardPage <- setRefClass("WizardPage",
                          fields=list(
                            wizard="ANY",
                            prev_button="ANY",
                            next_button="ANY",
                            widgets="list"
                            ),
                          methods=list(
                            initialize=function(wizard=NULL, ...) {
                              initFields(wizard=wizard)
                              callSuper(...)
                            },
                            can_next=function() {
                              "Return logical if we can go to next page"
                              TRUE
                            },
                            can_prev=function() {
                              "Return logical if we can go to previous"
                              TRUE
                            },
                            make_page=function(content_area) {
                              "Add content here"
                            },
                            update_page=function() {
                              "Called by change handlers in widgets"
                              wizard$update_page()
                            },
                            get_values=function() {
                              "Return values in a named list"
                              sapply(widgets, svalue, simplify=FALSE)
                            }
                          ))


Wizard <- setRefClass("Wizard",
                      fields=list(
                        pages="list",
                        cur_page="ANY",
                        main_window="ANY",
                        nb="ANY"
                        ),
                      methods=list(
                        initialize=function(title="", ...) {
                          initFields(pages=list(),
                                     main_window=gwindow(title, visible=FALSE)
                                     )
                          g <- ggroup(cont=main_window)
                          g$set_borderwidth(10)
                          nb <<- gstackwidget(cont=g)
                          
                          callSuper(...)
                        },
                        no_pages=function() length(pages),
                        add_page=function(page, title="") {
                          page$wizard <- .self
                          pages <<- c(pages, page)
                          box <- ggroup(cont=nb, label=title, horizontal=FALSE)
                          content_area <- ggroup(container=box, expand=TRUE, fill=TRUE, horizontal=FALSE)
                          page$make_page(content_area)
                          
                          button_area <- ggroup(cont=box, horizontal=TRUE)
                          addSpring(button_area)
                          page$prev_button <- gbutton("previous", cont=button_area, handler=function(h,...) {
                            h$action$prev_page()
                          }, action=.self)
                          page$next_button <- gbutton("next", cont=button_area, handler=function(h,...) {
                            h$action$next_page()
                          }, action=.self)
                        },
                        
                        prev_page=function() {
                          cur_page_no <- Filter(function(i) identical(pages[[i]], cur_page), seq_along(pages))
                          if(cur_page_no > 1) {
                            .nb <- nb
                            svalue(.nb) <- cur_page_no - 1
                            cur_page <<- pages[[cur_page_no - 1]]
                            update_page()
                          }
                        },
                        next_page = function() {
                          cur_page_no <- Filter(function(i) identical(pages[[i]], cur_page), seq_along(pages))
                          if(cur_page_no < length(pages)) {
                            .nb <- nb
                            svalue(.nb) <- cur_page_no + 1
                            cur_page <<- pages[[cur_page_no + 1]]
                            update_page()
                          } else {
                            call_finalizer()
                          }
                        },
                        call_finalizer=function() {
                          "replace me in subclass, or make configurable..."
                          print("Tada, all done")
                          print(get_values())
                          close_ui()
                        },
                        update_page=function() {
                          "Update buttons"
                          p_b <- cur_page$prev_button; n_b <- cur_page$next_button
                          enabled(p_b) <- cur_page$can_prev()
                          enabled(n_b) <- cur_page$can_next()
                        },
                        make_ui=function() {
                          .nb <- nb
                          svalue(.nb) <- 1
                          cur_page <<- pages[[1]]
                          update_page()
                          .main_window <- main_window
                          visible(.main_window) <- TRUE
                        },
                        close_ui=function() {
                          dispose(main_window)
                        },
                        get_values=function() {
                          "Return values from page"
                          out <- sapply(pages, function(page) page$get_values(), simplify=FALSE)
                          unlist(out, recursive=FALSE)
                        }
                        ))
                      

## Example implementation
wizard <- Wizard$new(title="test")

page1 <- setRefClass("Page1",
                     contains="WizardPage",
                     methods=list(
                       can_next = function() {
                         with(widgets,
                              nchar(svalue(ed1)) > 0 && nchar(svalue(ed2)) > 0
                              )
                       },
                       can_prev=function() FALSE,
                       make_page=function(content_area) {
                         lyt <- gformlayout(cont=content_area)
                         widgets$ed1 <<- gedit("", initial.msg="Enter some text", cont=lyt, label="Area 1")
                         widgets$ed2 <<- gedit("", initial.msg="and some more...", cont=lyt, label="And 2")

                         sapply(widgets, function(obj)
                                addHandlerKeystroke(obj, handler=function(h,...) h$action$update_page(), action=.self))
                       }
                       ))$new()
                       



page2 <- setRefClass("Page2",
                     contains="WizardPage",
                     methods=list(
                       can_next = function() {
                         with(widgets,
                              svalue(cb)
                              )
                       },
                       can_prev=function() TRUE,
                       make_page=function(content_area) {
                         lyt <- gformlayout(cont=content_area)
                         widgets$rb <<- gradio(state.name[1:4], cont=lyt, label="radio group")
                         widgets$cb <<- gcheckbox("all ready?", checked=FALSE, cont=lyt, label="checkbox")

                         sapply(widgets, function(obj)
                                addHandlerChanged(obj, handler=function(h,...) h$action$update_page(), action=.self))
                       }
                       ))$new()



wizard$add_page(page1, "Page 1")
wizard$add_page(page2, "Page 2")

wizard$make_ui()
