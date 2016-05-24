VVMain <-
  setRefClass("RzVVMain",
              fields = c("main"),
              contains = "RzVVCore",
              methods = list(
                initialize            = function(...) {
                  initFields(...)
                  callSuper(...)
                  
                  main <<- gtkVBox()
                  main$packStart(widget, expand=TRUE, fill=TRUE)
                }
              )
  )
VVMain$accessors("main")
