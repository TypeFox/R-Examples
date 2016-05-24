duplicateData <-
  setRefClass("RzVVDuplicateData",
              fields = c("main", "vvcore", "data"),
              methods = list(
                initialize = function(...) {
                  initFields(...)
                  
                  radio.button1 <- gtkRadioButtonNewWithLabel(label = gettext("All variables"))
                  radio.button2 <- gtkRadioButtonNewWithLabelFromWidget(radio.button1, gettext("Only selected variables"))
                  radio.button3 <- gtkRadioButtonNewWithLabelFromWidget(radio.button1, gettext("Only not selected variables"))
                  check.button1 <- gtkCheckButton(gettext("Only selected cases"))
                  check.button2 <- gtkCheckButton(gettext("Exclude missing values"))
                  button.execute <- gtkButtonNewFromStock(GTK_STOCK_EXECUTE)
                  hbox <- gtkHBoxNew()
                  hbox$packStart(button.execute, expand=FALSE)
                  
                  vbox1 <- gtkVBoxNew(spacing=4)
                  vbox1$setBorderWidth(2)
                  vbox1$packStart(radio.button1, expand=FALSE)
                  vbox1$packStart(radio.button2, expand=FALSE)
                  vbox1$packStart(radio.button3, expand=FALSE)
                  
                  frame.var <- gtkFrameNew(gettext("Variables"))
                  frame.var$add(vbox1)
                  
                  vbox2 <- gtkVBoxNew(spacing=4)
                  vbox2$setBorderWidth(2)
                  vbox2$packStart(check.button1, expand=FALSE)
                  vbox2$packStart(check.button2, expand=FALSE)
                  
                  frame.cases <- gtkFrameNew(gettext("Cases"))
                  frame.cases$add(vbox2)
                  
                  vbox.main <- gtkVBoxNew(spacing=4)
                  vbox.main$setBorderWidth(2)
                  vbox.main$packStart(frame.var, expand=FALSE)
                  vbox.main$packStart(frame.cases, expand=FALSE)
                  vbox.main$packStart(hbox, expand=FALSE)
                  
                  main <<- gtkScrolledWindowWithViewportNew()
                  main$add(vbox.main)
                  
                  gSignalConnect(button.execute, "clicked", function(...){
                    data.set <- data$getData.set()
                    
                    if (check.button1$getActive()) {
                      if (data$getSubset.on()) {
                        data.set <- data$getData.set.subset()
                      } else {
                        rzTools$runDialog(gettext("\"Select Cases\" isn't enabled."), type = "error")
                        return()
                      }
                    }
                    
                    if (radio.button2$getActive()) {
                      ind <- vvcore$getSelectedRows()
                      if (length(ind) == 0) {
                        rzTools$runDialog(gettext("No variables are selected."), type = "error")
                        return()
                      }
                      data.set <- data.set[ind]
                      
                    } else if (radio.button3$getActive()) {
                      ind <- vvcore$getSelectedRows()
                      if (length(ind) == 0) {
                        rzTools$runDialog(gettext("No variables are selected."), type = "error")
                        return()
                      }
                      data.set <- data.set[-ind]                      
                    }
                    
                    if (check.button2$getActive()) {
                      nas <- apply(as.data.frame(data.set), 1, function(x) any(is.na(x)))
                      nas <- which(!nas)
                      data.set <- data.set[nas, ]
                    }
                    
                    rzTools$addData(data.set,
                              name = gettextf("[Duplicated from %s]", data$getData.set.name()))
                  })
                  
                }                
              ))
duplicateData$accessors("main")
