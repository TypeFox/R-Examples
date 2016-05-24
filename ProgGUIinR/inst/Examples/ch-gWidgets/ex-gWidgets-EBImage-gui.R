### R code from vignette source 'ex-gWidgets-EBImage-gui.Rnw'

## How to install EBImage
# source("http://bioconductor.org/biocLite.R")
#    biocLite("EBImage")
#    biocLite("EBImage", type="source")

###################################################
### code chunk number 1: EBImageGUI
###################################################
## A simple GUI for EBImage as an example for basic controls and gWidgets
## Could be extended quite a bit

require(gWidgets)
options(guiToolkit="RGtk2")
require(EBImage)

setRefClass("ImageData",
            fields = list(                
              image = "Image",            # Original
              copy = "Image",             # copy
              tmp_filename = "character", # some filename
              
              ## widget instances
              toplevel = "gWindow",
              widgets = "list",
              img = "gImage",                # gimage instance
              actions = "list",            # gaction list
              
              ## parameter values
              colormode = "character",
              brightness = "numeric",
              contrast = "numeric",
              gamma = "numeric",
              flipflop = "character",
              thresh = "logical",
              filter = "character"
              ),
            methods = list(
              #
              initialize = function(f, tmp_filename = "/tmp/file.gif") {
                "Initialize GUI"
                ##' @param f filename
                copy <<- image <<- readImage(f)
                tmp_filename <<- tmp_filename

                init_gui()

                .self
              },
              #
              set = function(key, value, update = TRUE) {
                "Set key as value, update image if requested"
                assign(key, value, inherits = TRUE)
                if(update)
                  update_graphic()
              },
              #
              update_graphic = function() {
                "Update graphic for give field values"
                copy <- image
                colorMode(copy) <- colormode
                copy <- copy + brightness
                copy <- copy * contrast
                copy <- copy^gamma
                if("Flip" %in% flipflop) copy <- EBImage:::flip(copy)
                if("Flop" %in% flipflop) copy <- EBImage:::flop(copy)
                if(thresh) copy <- thresh(copy)
                if(filter == "Low") {
                  f <- makeBrush(21, shape = 'disc', step = FALSE)
                  f <- f/sum(f)
                  copy <- filter2(copy, f)
                } else if(filter == "High") {
                  f <- matrix(1, nc = 3, nr = 3); f[2,2] <- -8
                  copy <- filter2(copy, f)
                }
                
                copy <<- copy
                save()
                svalue(img) <- tmp_filename
              },
              #
              init_gui = function() {
                "Initialize GUI"
                toplevel <<- gwindow("EBImage GUI", visible = FALSE)
                g <- gpanedgroup(cont = toplevel, horizontal = TRUE)
                cg <- ggroup(cont = g, horizontal = FALSE)
                img <<- gimage(cont = g, width = 500, height = 480)

                actions <<- list(quit = gaction("Quit", icon = "quit", key.accel = "Control-q", parent = toplevel,
                                   handler = function(h,...) dispose(toplevel)),
                                 sep = gseparator(),
                                 open = gaction("Open", icon = "open", key.accel = "Control-o",parent = toplevel,
                                   handler = function(h,...) h$action$open(), action = .self),
                                 save = gaction("Save", icon = "save", key.accel = "Control-s", parent = toplevel,
                                   handler = function(h,...) h$action$save(), action = .self))
                gtoolbar(actions, cont = toplevel)

                make_controls(cg)
                init_controls()
                update_graphic()
                
                visible(toplevel) <- TRUE
              },
              #
              make_controls = function(g) {
                "Make the controls, save in widgets field"
                tbl <- glayout(cont = g)
                l <- list(); i <- 1
                def_handler <- function(h,...)  .self$set(h$action, svalue(h$obj))

                tbl[i,1] <- "ColorMode"
                tbl[i,2] <- (l$colormode <- gradio(c("Color","Grayscale"), 
                                                   selected = 1, cont = tbl, 
                                                   handler = def_handler, action = "colormode"))
                i <- i + 1
                

                
                tbl[i,1] <- "Brightness"
                tbl[i,2] <- (l$brightness <- gslider(from = -1, to = 1, by = .05, value = 0, 
                                                     handler = def_handler, action = "brightness"))
                i <- i+1
                
                tbl[i,1] <- "Contrast"
                tbl[i,2] <- (l$contrast <- gspinbutton(from = 0, to = 10, by = .05, value = 1, 
                                                       handler = def_handler, action = "contrast"))
                i <- i + 1

                tbl[i,1] <- "Gamma"
                tbl[i,2] <- (l$gamma <- gslider(from = 0, to = 10, by = .05, value = 1, 
                                                handler = def_handler, action = "gamma"))
                i <- i + 1

                tbl[i,1] <- "Flip/Flop"
                tbl[i,2] <- (l$flipflop <- gcheckboxgroup(c("Flip","Flop"), 
                                                          handler = def_handler, action = "flipflop"))
                i <- i + 1

                tbl[i,2] <- (l$thresh <- gcheckbox("Thresh", use.togglebutton = TRUE,
                                      selected = FALSE, handler = def_handler, action = "thresh"))
                i <- i + 1

                tbl[i,1] <- "Filter"
                tbl[i,2] <- (l$filter <- gcombobox(c("None", "Low", "High"), cont = tbl, 
                                                   handler = def_handler, action = "filter"))
                i <- i + 1
                
                widgets <<- l
              },
              #
              init_controls = function(l) {
                "Initialize the controls and the fields. Does not update graphic."
                if(missing(l)) {
                  l <- list(colormode = "Grayscale",
                            brightness = 0,
                            contrast = 1,
                            gamma = 1,
                            flipflop  =  character(0),
                            thresh = FALSE,
                            filter = "None")
                }
                sapply(names(l), function(i) {
                  set(i, l[[i]], update = FALSE)
                  blockHandler(widgets[[i]]) # don't update graphic
                  svalue(widgets[[i]]) <- l[[i]]
                  unblockHandler(widgets[[i]])
                })
              },
              #
              open = function() {
                "Open a new file"
                f <- gfile("Open an image file",
                           filter = list("Image file" = list(patterns = c("*.gif", "*.jpeg", "*.png")),
                             "All files" = list(patterns = c("*"))
                             ))
                if(is.na(f)) return()
                copy <<- image <<- readImage(f)
                update_graphic()
              },
              #
              save = function(quality = 100) {
                "Write copy to tmp_filename"
                writeImage(copy, tmp_filename, quality = quality)
              }
              ))


## Test it out
f <- system.file("images","lena.gif", package = "EBImage")
img <- getRefClass("ImageData")$new(f)


