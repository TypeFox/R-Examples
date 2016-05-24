### R code from vignette source 'ch-RGtk2.Rnw'

###################################################
### code chunk number 1: ch-RGtk2.Rnw:13-38
###################################################
options(prompt=" ")
options(continue=" ")
options(width=60)

findmethod <- function (obj, name, where=.GlobalEnv) 
{
    classes <- c(attr(obj, "interfaces"), class(obj))
    sym <- paste(tolower(substring(classes, 1, 1)), substring(classes, 
        2), toupper(substring(name, 1, 1)), substring(name, 2), 
        sep = "")
    which <- sapply(sym, exists, where)
    if (!any(which)) 
        stop(paste("No such method", name, "for classes", paste(class(obj), 
            collapse = ", ")))
    return(sym[which][1])
}

## override NULL in output
args <- function(name) {
  body(name) <- NULL
  environment(name) <- .GlobalEnv
  name
  ## out <- capture.output(base::args(name))
  ## invisible(cat(out[-length(out)], "\n"))
}


###################################################
### code chunk number 2: Introduction.Rnw:2-3
###################################################
require(RGtk2)


###################################################
### code chunk number 3: gtk-overview-initial-example
###################################################
button <- gtkButton("Click Me")
button['image'] <- gtkImage(stock = "gtk-apply", 
                            size = "button")
gSignalConnect(button, "clicked", function(button) {
  message("Hello World!")
})
##
window <- gtkWindow(show = FALSE)
window$add(button)
window$showAll()


###################################################
### code chunk number 4: gtk-intro-classes-ancestors
###################################################
gTypeGetAncestors("GtkWidget")


###################################################
### code chunk number 5: gtk-intro-class-interfaces
###################################################
gTypeGetInterfaces("GtkWidget")


###################################################
### code chunk number 6: intro-constructor-gtkWindow
###################################################
window <- gtkWindow("toplevel", show = FALSE)


###################################################
### code chunk number 7: gtk-overview-construct-image
###################################################
gtkImage(stock = "gtk-apply", size = "button")


###################################################
### code chunk number 8: gtk-overview-construct-image-args (eval = FALSE)
###################################################
## args(gtkImage)


###################################################
### code chunk number 9: gtk-overview-objects-value
###################################################
a <- -1
abs(a)
a


###################################################
### code chunk number 10: gtk-overview-objects-ref
###################################################
gtkButtonSetLabel(button, "New text")
gtkButtonGetLabel(button)


###################################################
### code chunk number 11: intro-constructor-classes (eval = FALSE)
###################################################
## class(window)


###################################################
### code chunk number 12: intro-constructor-interfaces
###################################################
interface(window)


###################################################
### code chunk number 13: intro-methods-button
###################################################
button <- gtkButton("Hello World")
window$add(button)
window$setDefaultSize(200, 200)


###################################################
### code chunk number 14: gtkButtonAddApi
###################################################
gtkButtonSayHello <- function(obj, target) 
  obj$setLabel(paste("Hello", target))
button$sayHello("World")
button$getLabel()


###################################################
### code chunk number 15: showProperties
###################################################
head(names(button), n = 8)             # or b$getPropInfo()


###################################################
### code chunk number 16: intro-props-get-set
###################################################
image <- gdkPixbuf(filename = imagefile("rgtk-logo.gif"))
window$set(icon = image[[1]], title = "Hello World 1.0")


###################################################
### code chunk number 17: Introduction.Rnw:354-356
###################################################
window$setTitle("Hello World 1.0")
window$getTitle()


###################################################
### code chunk number 18: intro-props-visible
###################################################
window["visible"]


###################################################
### code chunk number 19: intro-props-show
###################################################
window["visible"] <- TRUE 
window$show() # same effect


###################################################
### code chunk number 20: Introduction.Rnw:388-389
###################################################
names(gTypeGetSignals("GtkButton"))


###################################################
### code chunk number 21: Introduction.Rnw:399-400 (eval = FALSE)
###################################################
## args(gSignalConnect)


###################################################
### code chunk number 22: intro-signals-hello-world
###################################################
gSignalConnect(button, "clicked", 
               function(button) message("Hello World!"))


###################################################
### code chunk number 23: Introduction.Rnw:446-455
###################################################
window <- gtkWindow(); window['title'] <- "test signals"
x <- 1; 
button <- gtkButton("click me"); window$add(button)
gSignalConnect(button, signal = "clicked", 
               f = function(button) {
                 button$setData("x", 2)
                 x <- 2
                 return(TRUE)
               })


###################################################
### code chunk number 24: Introduction.Rnw:458-459
###################################################
button$setData("x", 2)                        # fix non-interactivity


###################################################
### code chunk number 25: Introduction.Rnw:462-463
###################################################
cat(x, button$getData("x"), "\n") # 1 and 2


###################################################
### code chunk number 26: Introduction.Rnw:474-489
###################################################
button <- gtkButton("click")
window <- gtkWindow()
window$add(button)
gSignalConnect(button, "button-press-event", 
               function(button, event, data) {
                 message("hi"); return(FALSE)
               })
gSignalConnect(button, "button-press-event", 
               function(button, event, data) {
                 message("and"); return(TRUE)
               })
gSignalConnect(button, "button-press-event", 
               function(button, event, data) {
                 message("bye"); return(TRUE)
               })


###################################################
### code chunk number 27: intro-enum-window (eval = FALSE)
###################################################
## window <- gtkWindow("toplevel", show = FALSE)


###################################################
### code chunk number 28: intro-enum-GtkWindowType
###################################################
GtkWindowType


###################################################
### code chunk number 29: intro-enum-GtkWidgetFlags
###################################################
GtkWidgetFlags


###################################################
### code chunk number 30: intro-enum-gtkWidgetFlags
###################################################
window$flags()


###################################################
### code chunk number 31: intro-enum-toplevel
###################################################
(window$flags() & GtkWidgetFlags["toplevel"]) > 0


###################################################
### code chunk number 32: Introduction.Rnw:576-578
###################################################
while(gtkEventsPending()) 
  gtkMainIteration()


###################################################
### code chunk number 33: Glade.Rnw:16-18
###################################################
builder <- gtkBuilder()
builder$addFromFile("buildable.xml")


###################################################
### code chunk number 34: Glade.Rnw:25-27
###################################################
dialog1 <- builder$getObject("dialog1")
dialog1$showAll()


###################################################
### code chunk number 35: Glade.Rnw:39-43
###################################################
ok_button_clicked <- function(button, userData) {
  message("hello world")
}
builder$connectSignals()


###################################################
### code chunk number 36: Containers.Rnw:9-10
###################################################
library(RGtk2)


###################################################
### code chunk number 37: Containers.Rnw:47-52
###################################################
window <- gtkWindow(show=FALSE)          # use default type
window$setTitle("Window title")          # set window title
window['title']                          # or  use getTitle
window$setDefaultSize(250,300)           # 250 wide, 300 high
window$show()                            # show window


###################################################
### code chunk number 38: basic-window-label
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Hello World")
label <- gtkLabel("Hello World")
window$add(label)


###################################################
### code chunk number 39: gtk-container-window-delete
###################################################
gSignalConnect(window, "delete-event", function(event,...) {
  dialog <- gtkMessageDialog(parent = window, flags = 0, 
                             type = "question", 
                             buttons = "yes-no",
                             "Are you sure you want to quit?")
  out <- dialog$run(); dialog$destroy()
  out != GtkResponseType["yes"]
})


###################################################
### code chunk number 40: gtk-container-window-destroy
###################################################
window$destroy()


###################################################
### code chunk number 41: Containers.Rnw:131-142
###################################################
## create a window and a dialog window
window <- gtkWindow(show = FALSE)
window$setTitle("Top level window")
##
dialog <- gtkWindow(show = FALSE)
dialog$setTitle("dialog window")
dialog$setTransientFor(window)
dialog$setPosition("center-on-parent")
dialog$setDestroyWithParent(TRUE)
window$show()
dialog$show()


###################################################
### code chunk number 42: Containers.Rnw:169-173
###################################################
window <- gtkWindow(show=FALSE)
window$setTitle("Hello World")
label <- gtkLabel("Hello World")
window$add(label)


###################################################
### code chunk number 43: Containers.Rnw:181-182
###################################################
window$getChild()['label']


###################################################
### code chunk number 44: gtk-container-brackets
###################################################
window[[1]]['label']


###################################################
### code chunk number 45: Containers.Rnw:195-197 (eval = FALSE)
###################################################
## ## leave out?
## l$getParent()


###################################################
### code chunk number 46: remove-child-widget-3
###################################################
window$remove(label)
window$add(label)


###################################################
### code chunk number 47: layout-window-show-first
###################################################
window <- gtkWindow()
window$setTitle("Hello World")
label <- gtkLabel("Hello World")
window$add(label)


###################################################
### code chunk number 48: layout-window-show-first-alloc
###################################################
label$getAllocation()$allocation


###################################################
### code chunk number 49: layout-window-show-later
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Hello World")
label <- gtkLabel("Hello World")
window$add(label)
window$show()
label$getAllocation()$allocation


###################################################
### code chunk number 50: Containers.Rnw:268-269
###################################################
sapply(label$getAllocation()$allocation, function(i) i)


###################################################
### code chunk number 51: basic-box-homo
###################################################
box <- gtkHBox(TRUE, 5)


###################################################
### code chunk number 52: basic-box-homo-nofill
###################################################
button_a <- gtkButton("Button A")
button_b <- gtkButton("Button B")
box$packStart(button_a, fill = FALSE)
box$packStart(button_b, fill = FALSE)


###################################################
### code chunk number 53: basic-box-hetero
###################################################
box <- gtkHBox(FALSE, 5)


###################################################
### code chunk number 54: Containers.Rnw:414-417
###################################################
## re create buttons
button_a <- gtkButton("Button A")
button_b <- gtkButton("Button B")


###################################################
### code chunk number 55: basic-box-expand
###################################################
box$packStart(button_a, expand = TRUE, fill = FALSE)
box$packStart(button_b, expand = FALSE, fill = FALSE)


###################################################
### code chunk number 56: Containers.Rnw:463-465
###################################################
hbox <- gtkHBox()
sapply(1:3, function(i) hbox$packStart(gtkLabel(i)))


###################################################
### code chunk number 57: Containers.Rnw:467-469
###################################################
b3 <- hbox[[3]]
hbox$reorderChild(b3, 2 - 1)               # second is 2 - 1


###################################################
### code chunk number 58: basic-layout-align-window
###################################################
window <- gtkWindow(); window$setTitle("Hello World")
label <- gtkLabel("Hello World")
window$add(label)


###################################################
### code chunk number 59: basic-layout-align-left
###################################################
label["xalign"] <- 0


###################################################
### code chunk number 60: basic-layout-align-GtkAlignment
###################################################
window <- gtkWindow(); window$setTitle("Hello World")
alignment <- gtkAlignment()
alignment$set(xalign = 0, yalign = 0.5, xscale = 0, yscale=1)
window$add(alignment)
label <- gtkLabel("Hello World")
alignment$add(label)


###################################################
### code chunk number 61: Pre-defined-dialogs.Rnw:21-29
###################################################
window <- gtkWindow(); window['title'] <- "Parent window"
#
dialog <- gtkMessageDialog(parent=window, 
                           flags="destroy-with-parent",
                           type="question", 
                           buttons="ok",
                           "My message")
dialog['secondary-text'] <- "A secondary message"


###################################################
### code chunk number 62: Pre-defined-dialogs.Rnw:49-58
###################################################
response <- dialog$run()
if(response == GtkResponseType["cancel"] ||
   response == GtkResponseType["close"] ||
   response == GtkResponseType["delete-event"]) {
  ## pass
} else if(response == GtkResponseType["ok"]) {
  message("Ok")
}
dialog$destroy()


###################################################
### code chunk number 63: Pre-defined-dialogs.Rnw:86-91
###################################################
dialog <- gtkDialogNewWithButtons(title = "Enter a value", 
                   parent = NULL, flags = 0,
                   "gtk-ok", GtkResponseType["ok"],
                   "gtk-cancel", GtkResponseType["cancel"],
                   show = FALSE)


###################################################
### code chunk number 64: OurDialogsLayout
###################################################
hbox <- gtkHBox()
hbox['spacing'] <- 10
#
hbox$packStart(gtkLabel("Enter a value:"))
entry <- gtkEntry()
hbox$packStart(entry)
#
vbox <- dialog$getContentArea()
vbox$packStart(hbox)


###################################################
### code chunk number 65: connectResponse
###################################################
gSignalConnect(dialog, "response", 
               f=function(dialog, response, user.data) {
                 if(response == GtkResponseType["ok"])
                   print(entry$getText()) # Replace this
                 dialog$Destroy()
               })
dialog$showAll()
dialog$setModal(TRUE)


###################################################
### code chunk number 66: openFileDialog
###################################################
dialog <- gtkFileChooserDialog(title = "Open a file", 
                     parent = NULL, action = "open",
                     "gtk-ok", GtkResponseType["ok"],
                     "gtk-cancel", GtkResponseType["cancel"],
                     show = FALSE)


###################################################
### code chunk number 67: Pre-defined-dialogs.Rnw:165-173
###################################################
gSignalConnect(dialog, "response", 
               f = function(dialog, response, data) {
                 if(response == GtkResponseType["ok"]) {
                   filename <- dialog$getFilename()
                   print(filename)
                 }
                 dialog$destroy()
               })


###################################################
### code chunk number 68: Pre-defined-dialogs.Rnw:183-188
###################################################
fileFilter <- gtkFileFilter()
fileFilter$setName("R files")
fileFilter$addPattern("*.R")
fileFilter$addPattern("*.Rdata")
dialog$addFilter(fileFilter)


###################################################
### code chunk number 69: gtk-container-frame
###################################################
frame <- gtkFrame("Options")
vbox <- gtkVBox()
vbox$packStart(gtkCheckButton("Option 1"), FALSE)
vbox$packStart(gtkCheckButton("Option 2"), FALSE)
frame$add(vbox)


###################################################
### code chunk number 70: gtk-container-expander
###################################################
expander <- gtkExpander("Advanced")
expander$add(frame)


###################################################
### code chunk number 71: qt-layout-notebook
###################################################
notebook <- gtkNotebook()
notebook$appendPage(gtkLabel("Page 1"), gtkLabel("Tab 1"))
notebook$appendPage(gtkLabel("Page 2"), gtkLabel("Tab 2"))


###################################################
### code chunk number 72: qt-layout-notebook-pos
###################################################
notebook['tab-pos'] <- "bottom"


###################################################
### code chunk number 73: qt-layout-notebook-current
###################################################
notebook['page'] <- 1
notebook['page']


###################################################
### code chunk number 74: Containers.Rnw:656-674
###################################################
gtkNotebookInsertPageWithCloseButton <- 
  function(object, child, label.text="", position=-1) {
    icon <- gtkImage(pixbuf = 
      object$renderIcon("gtk-close", "button", size = "menu"))
    closeButton <- gtkButton()
    closeButton$setImage(icon)
    closeButton$setRelief("none")
    ##
    label <- gtkHBox()
    label$packStart(gtkLabel(label.text))
    label$packEnd(closeButton)
    ##
    gSignalConnect(closeButton, "clicked", function(button) {
      index <- object$pageNum(child)
      object$removePage(index)
    })
    object$insertPage(child, label, position)
  }


###################################################
### code chunk number 75: Containers.Rnw:679-685
###################################################
window <- gtkWindow()
notebook <- gtkNotebook(); window$add(notebook)
notebook$insertPageWithCloseButton(gtkButton("hello"), 
                                   label.text = "page 1")
notebook$insertPageWithCloseButton(gtkButton("world"), 
                                   label.text = "page 2")


###################################################
### code chunk number 76: gtk-container-scrolled-device
###################################################
library(cairoDevice)
device <- gtkDrawingArea()
device$setSizeRequest(600, 400)
asCairoDevice(device)


###################################################
### code chunk number 77: gtk-container-scrolled-construct
###################################################
scrolled <- gtkScrolledWindow()
scrolled$addWithViewport(device)


###################################################
### code chunk number 78: gtk-container-scrolled-zoom
###################################################
zoomPlot <- function(x = 2.0) {
  allocation <- device$getAllocation()$allocation
  device$setSizeRequest(allocation$width * x, 
                        allocation$height * x)
  updateAdjustment <- function(adjustment) {
    adjustment$setValue(x * adjustment$getValue() + 
                        (x - 1) * adjustment$getPageSize()/2)
  }
  updateAdjustment(scrolled$getHadjustment())
  updateAdjustment(scrolled$getVadjustment())
}


###################################################
### code chunk number 79: gtk-container-scrolled-key-press
###################################################
gSignalConnect(scrolled, "key-press-event", 
               function(scrolled, event) {
                 key <- event[["keyval"]]
                 if (key == GDK_plus)
                   zoomPlot(2.0)
                 else if (key == GDK_minus)
                   zoomPlot(0.5)
                 TRUE
               })


###################################################
### code chunk number 80: gtk-container-scrolled-window
###################################################
win <- gtkWindow(show = FALSE)
win$add(scrolled)
win$showAll()


###################################################
### code chunk number 81: gtk-container-scrolled-plot
###################################################
plot(mpg ~ hp, data = mtcars)


###################################################
### code chunk number 82: gtk-container-paned-construct
###################################################
paned <- gtkHPaned()


###################################################
### code chunk number 83: gtk-container-paned-add
###################################################
paned$add1(gtkLabel("Left (1)"))
paned$add2(gtkLabel("Right (2)"))


###################################################
### code chunk number 84: gtk-container-paned-pack
###################################################
paned$pack1(gtkLabel("Left (1)"), resize = TRUE, shrink=TRUE)
paned$pack2(gtkLabel("Right (2)"), resize = TRUE, shrink=TRUE)


###################################################
### code chunk number 85: ex-RGtk2-dialog-layout.Rnw:4-6
###################################################
## layout a basic dialog -- center align
library(RGtk2)


###################################################
### code chunk number 86: gtk-container-table-construct
###################################################
table <- gtkTable(rows = 3, columns = 2, homogeneous = FALSE)


###################################################
### code chunk number 87: ex-RGtk2-dialog-layout.Rnw:26-40
###################################################
size_label <- gtkLabel("Sample size:")
size_combo <- gtkComboBoxNewText()
sapply(c(5, 10, 15, 30), size_combo$appendText)
##
diag_label <- gtkLabel("Diagnostic:")
diag_radio <- gtkVBox()
radiogp <- list()
radiogp$t <- gtkRadioButton(label = "t-statistic")
radiogp$mean <- gtkRadioButton(radiogp, label = "mean")
radiogp$median <- gtkRadioButton(radiogp, label = "median")
sapply(radiogp, diag_radio$packStart)
##
submit_vbox <- gtkVBox()
submit_vbox$packEnd(gtkButton("Run simulation"), expand=FALSE)


###################################################
### code chunk number 88: gtk-container-layout-align
###################################################
size_label['xalign'] <- 1
diag_label['xalign'] <- 1; diag_label['yalign'] <- 0
diag_align <- gtkAlignment(xalign = 0)
diag_align$add(diag_radio)


###################################################
### code chunk number 89: ex-RGtk2-dialog-layout.Rnw:80-94
###################################################
table$attach(size_label, left.attach = 0,1, top.attach = 0,1, 
             xoptions = c("expand", "fill"), yoptions = "")
table$attach(size_combo, left.attach = 1,2, top.attach = 0,1, 
             xoptions = "fill", yoptions = "")
##
table$attach(diag_label, left.attach = 0,1, top.attach = 1,2, 
             xoptions = c("expand", "fill"), 
             yoptions = c("expand", "fill"))
##
table$attach(diag_align, left.attach = 1,2, top.attach = 1,2, 
             xoptions = c("expand", "fill"), yoptions = "")
##
table$attach(submit_vbox, left.attach = 1,2, top.attach = 2,3, 
             xoptions = "", yoptions = c("expand", "fill"))


###################################################
### code chunk number 90: gtk-container-table-spacing
###################################################
table$setColSpacing(0, 5)


###################################################
### code chunk number 91: ex-RGtk2-dialog-layout.Rnw:111-115
###################################################
window <- gtkWindow(show=FALSE)
window['border-width'] <- 14
window$setTitle("GtkTable Example")
window$add(table)


###################################################
### code chunk number 92: ex-RGtk2-dialog-layout.Rnw:118-119
###################################################
window$show()


###################################################
### code chunk number 93: ButtonConstructors
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Various buttons")
window$setDefaultSize(400, 25)
hbox <- gtkHBox(homogeneous = FALSE, spacing = 5)
window$add(hbox)
button <- gtkButtonNew() 
button$setLabel("long way")
hbox$packStart(button)
hbox$packStart(gtkButton(label = "label only") )
hbox$packStart(gtkButton(stock.id = "gtk-ok") )
hbox$packStart(gtkButtonNewWithMnemonic("_Mnemonic") )
window$show()


###################################################
### code chunk number 94: CallbackExampleForButton
###################################################
window <- gtkWindow(); button <- gtkButton("click me");
window$add(button)

gSignalConnect(button, "button-press-event", # just mouse
               f = function(widget, event, data) {
                 print(event$getButton())    # which button
                 return(FALSE)               # propagate
               })
gSignalConnect(button, "clicked",            # keyboard too
               f = function(widget, ...) {
                 print("clicked")
               })


###################################################
### code chunk number 95: gtk-widget-button-sensitive
###################################################
button$setSensitive(FALSE)


###################################################
### code chunk number 96: MacOSXstyleButton
###################################################
## not shown
window <- gtkWindow(show=FALSE)
window$setTitle("MAC OS X style buttons")
fg <- gtkVBox()
fg$setSizeRequest(width=800, height=-1)
window$add(fg)

hbox <- gtkHBox()
fg$packStart(hbox, padding=15)              # for size grip


###################################################
### code chunk number 97: StockButtons
###################################################
ok <- gtkButton(stock.id="gtk-ok")
cancel <- gtkButton(stock.id="gtk-cancel")
delete <- gtkButton(stock.id="gtk-delete")


###################################################
### code chunk number 98: macButtonPack
###################################################
hbox$packEnd(ok, padding = 0)
hbox$packEnd(cancel, padding = 12)
hbox$packEnd(delete, padding = 12)
hbox$packEnd(gtkLabel(""), expand = TRUE, fill = TRUE)
##
ok$grabFocus()


###################################################
### code chunk number 99: ex-RGtk2-mac-buttons.Rnw:60-62
###################################################
## not shown
window$showAll()


###################################################
### code chunk number 100: gtkHButtonBoxExample
###################################################
## not shown
## Had we only wanted to use a button box
button_box <- gtkHButtonBox()
fg$packStart(button_box, padding=15)              # for size grip

button_box$add(gtkButton(stock.id="gtk-delete"))
button_box$add(gtkButton(stock.id="gtk-cancel"))
button_box$add(gtkButton(stock.id="gtk-ok"))


###################################################
### code chunk number 101: gtk-widget-label-window
###################################################
window <- gtkWindow(show=FALSE)
window$setTitle("Label formatting")
window$setSizeRequest(250,300)               # narrow
vbox <- gtkVBox(spacing=2); vbox$setBorderWidth(5); window$add(vbox)


###################################################
### code chunk number 102: LabelFormatting
###################################################
string <- "the quick brown fox jumped over the lazy dog"
## wrap by setting number of characters
basicLabel <- gtkLabel(string)
basicLabel$setLineWrap(TRUE)
basicLabel$setWidthChars(35)            # no. characters

## Set ellipsis to shorten long text
ellipsized <- gtkLabel(string)
ellipsized$setEllipsize("middle")

## Right justify text lines
## use xalign property for aligning entire block
rightJustified <- gtkLabel("right justify")
rightJustified$setJustify("right")
rightJustified['xalign'] <- 1

## PANGO markup
pangoLabel <- gtkLabel()
tmpl <- "<span foreground='blue' size='x-small'>%s</span>"
pangoLabel$setMarkup(sprintf(tmpl, string))
#
sapply(list(basicLabel,ellipsized,rightJustified, pangoLabel), 
       vbox$packStart, expand = TRUE, fill = TRUE)
window$showAll()


###################################################
### code chunk number 103: ex-RGtk2-ImageForGraphics.Rnw:1-2
###################################################
library(RGtk2)


###################################################
### code chunk number 104: ex-RGtk2-ImageForGraphics.Rnw:12-17
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Graphic window")
window$setSizeRequest(400, 400)
hbox <- gtkHBox(); window$add(hbox)
window$showAll()


###################################################
### code chunk number 105: ex-RGtk2-ImageForGraphics.Rnw:25-27
###################################################
theSize <- hbox$getAllocation()$allocation
width <- theSize$width; height <- theSize$height


###################################################
### code chunk number 106: ex-RGtk2-ImageForGraphics.Rnw:33-38
###################################################
require(cairoDevice)
pixmap <- gdkPixmap(drawable = NULL, 
                    width = width, height = height, depth=24)
asCairoDevice(pixmap)
hist(rnorm(100))


###################################################
### code chunk number 107: ex-RGtk2-ImageForGraphics.Rnw:43-45
###################################################
image <- gtkImage(pixmap = pixmap)
hbox$packStart(image, expand = TRUE, fill = TRUE)


###################################################
### code chunk number 108: notShown
###################################################
## Work this into an example ###
makeIconRGtk2 <- function(widget, giffile) {
  if(checkPtrType(w, "GtkWindow")) {
    img <- gdkPixbufNewFromFile(giffile)
    if(!is.null(img$retval))
      widget$setIcon(img$retval)
  }
}


###################################################
### code chunk number 109: gtkStockListIds
###################################################
head(unlist(gtkStockListIds()), n=3)   


###################################################
### code chunk number 110: gtk-widget-entry
###################################################
entry <- gtkEntry()


###################################################
### code chunk number 111: gtk-widget-entry-activate
###################################################
gSignalConnect(entry, "activate", function() {
  message("Text entered: ", entry$getText())
})


###################################################
### code chunk number 112: gtk-widget-entry-validate
###################################################
validatedEntry <- gtkEntry()
gSignalConnect(validatedEntry, "changed", function(entry) {
  text <- entry$getText()
  if (nzchar(gsub("[a-zA-Z]", "", text))) {
    entry$setIconFromStock("primary", "gtk-no")
    entry$setIconTooltipText("primary", 
                                 "Only letters are allowed")
  } else { 
    entry$setIconFromStock("primary", "gtk-yes")
    entry$setIconTooltipText("primary", NULL)
  }
})
validatedEntry$setIconFromStock("primary", "gtk-yes")


###################################################
### code chunk number 113: BasicComponents.Rnw:430-433
###################################################
w <- gtkWindow(show=FALSE)
w$add(validatedEntry)
w$showAll()


###################################################
### code chunk number 114: gtk-widget-check-button
###################################################
checkButton <- gtkCheckButton("Option")


###################################################
### code chunk number 115: gtk-widget-check-button-active
###################################################
checkButton['active']
checkButton['active'] <- TRUE


###################################################
### code chunk number 116: gtk-widget-check-button-toggle
###################################################
gSignalConnect(checkButton, "toggled", function(button) {
  state <- ifelse(button$active, "active","inactive")
  message("Button is ", state)
})


###################################################
### code chunk number 117: RadioGroupExample
###################################################
labels <- c("two.sided", "less", "greater")
radiogp <- list()                           # list for group
radiogp[[labels[1]]] <- gtkRadioButton(label=labels[1]) 
for(label in labels[-1]) 
  radiogp[[label]] <- gtkRadioButton(radiogp, label=label)


###################################################
### code chunk number 118: BasicComponents.Rnw:509-512
###################################################
window <- gtkWindow(); window$setTitle("Radio group example")
vbox <- gtkVBox(FALSE, 5); window$add(vbox)
sapply(radiogp, gtkBoxPackStart, object = vbox)


###################################################
### code chunk number 119: BasicComponents.Rnw:516-518
###################################################
vbox[[3]]$setActive(TRUE)           
sapply(radiogp, `[`, "active") 


###################################################
### code chunk number 120: BasicComponents.Rnw:523-528
###################################################
sapply(radiogp, gSignalConnect, "toggled",     # connect each
       f = function(button, data) {
         if(button['active']) # set before callback
           message("clicked", button$getLabel(),"\n")
       })


###################################################
### code chunk number 121: BasicComponents.Rnw:537-544
###################################################
radiogp <- gtkRadioButton(label=labels[1])
btns <- sapply(labels[-1], gtkRadioButtonNewWithLabelFromWidget, 
               group = radiogp)
window <- gtkWindow()
window['title'] <- "Radio group example"
vbox <- gtkVBox(); window$add(vbox)
sapply(rev(radiogp$getGroup()), gtkBoxPackStart, object = vbox)


###################################################
### code chunk number 122: gtk-widget-combo
###################################################
combo <- gtkComboBoxNewText()
sapply(c("two.sided", "less", "greater"), combo$appendText)


###################################################
### code chunk number 123: gtk-widget-combo-active
###################################################
combo['active']


###################################################
### code chunk number 124: gtk-widget-combo-changed
###################################################
gSignalConnect(combo, "changed",
               f = function(button, ...) {
                 if(button$getActive() < 0) 
                   message("No value selected")
                 else
                   message("Value is", button$getActiveText())
               })


###################################################
### code chunk number 125: ComboBoxExample
###################################################
## An example of two comboboxes where 1 updates the other
require(RGtk2)
data(mtcars); library(MASS); data(Cars93) # need some data frames


###################################################
### code chunk number 126: ex-RGtk2-comboboxes.Rnw:11-12
###################################################
library(ProgGUIinR)                     # for avail_dfs, find_vars


###################################################
### code chunk number 127: Widgets
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("gtkComboBox example")

df_combo <- gtkComboBoxNewText()
var_combo <- gtkComboBoxNewText()


###################################################
### code chunk number 128: Layout
###################################################
vbox <- gtkVBox(); window$add(vbox)
#
vbox1 <- gtkHBox(); vbox$packStart(vbox1)
vbox1$packStart(gtkLabel("Data frames:"))
vbox1$packStart(df_combo)
#
vbox2 <- gtkHBox(); vbox$packStart(vbox2)
vbox2$packStart(gtkLabel("Variable:"))
vbox2$packStart(var_combo)
vbox2$hide()


###################################################
### code chunk number 129: configureComboBoxes
###################################################
sapply(avail_dfs(), df_combo$appendText)
df_combo$setActive(-1)
#
gSignalConnect(df_combo, "changed", function(df_combo, ...) {
  var_combo$getModel()$clear()
  sapply(find_vars(df_combo$getActiveText()),  
         var_combo$appendText)
  vbox2$show()
})


###################################################
### code chunk number 130: ex-RGtk2-comboboxes.Rnw:56-58
###################################################
## show window
window$show()


###################################################
### code chunk number 131: ex-RGtk2-range-widget.Rnw:15-17
###################################################
## make a range widget combining both a slider and spinbutton to choose a number
library(RGtk2)


###################################################
### code chunk number 132: ex-RGtk2-range-widget.Rnw:22-23
###################################################
from <- 0; to <- 100; by <- 1


###################################################
### code chunk number 133: ex-RGtk2-range-widget.Rnw:30-34
###################################################
slider <- gtkHScale(min = from, max = to, step = by)
slider['draw-value'] <- FALSE
adjustment <- slider$getAdjustment()
spinbutton <- gtkSpinButton(adjustment = adjustment)


###################################################
### code chunk number 134: ex-RGtk2-range-widget.Rnw:41-44
###################################################
hbox <- gtkHBox()
hbox$packStart(slider, expand=TRUE, fill = TRUE, padding = 5)
hbox$packStart(spinbutton, expand = FALSE, padding = 5)


###################################################
### code chunk number 135: ex-RGtk2-range-widget.Rnw:48-53
###################################################
w <- gtkWindow(show=FALSE)
w['title'] <- "Example of a range widget"
w$setSizeRequest(width=200, height=-1)
w$add(hbox)
w$show()


###################################################
### code chunk number 136: BasicComponents.Rnw:769-779
###################################################
window <- gtkWindow(); window$setTitle("Progress bar example")
progress_bar <- gtkProgressBar()
window$add(progress_bar)
#
progress_bar$setText("Please be patient...")
for(i in 1:100) {
  progress_bar$setFraction(i/100)
  Sys.sleep(0.05) ## replace with a step in the process
}
progress_bar$setText("All done.")


###################################################
### code chunk number 137: gtk-widget-progress-pulse
###################################################
progress_bar$pulse()


###################################################
### code chunk number 138: gtk-widget-spinner (eval = FALSE)
###################################################
## spinner <- gtkSpinner()
## spinner$start()
## spinner$stop()


###################################################
### code chunk number 139: installPackagesWizard
###################################################
## gtk Assistant example
require(RGtk2)


###################################################
### code chunk number 140: defineAssistant
###################################################
assistant <- gtkAssistant(show=FALSE)
assistant$setSizeRequest(500, 500)
gSignalConnect(assistant, "cancel", 
               function(assistant) assistant$destroy())


###################################################
### code chunk number 141: makePages
###################################################
pages <- lapply(1:5, gtkVBox, spacing=5, homogeneous = FALSE)
page_types <- c("intro", rep("confirm", 3), "summary")
sapply(pages, gtkAssistantAppendPage, object = assistant)
sapply(pages, gtkAssistantSetPageType, object = assistant, 
       type=page_types)


###################################################
### code chunk number 142: sideLogo1
###################################################
image <- gdkPixbuf(filename = imagefile("rgtk-logo.gif"))[[1]]
sapply(pages, gtkAssistantSetPageSideImage, object=assistant, 
       pixbuf = image)


###################################################
### code chunk number 143: ex-RGtk2-install-wizard.Rnw:52-59
###################################################
populate_page <- list()                
gSignalConnect(assistant, "prepare", 
       function(assistant, page, data) {
         page_no <- which(sapply(pages, identical, page))
         if(!length(page$getChildren()))
           populate_page[[page_no]]()
       })


###################################################
### code chunk number 144: ex-RGtk2-install-wizard.Rnw:68-74
###################################################
assistant$setForwardPageFunc(function(page_index, data) {
  if(page_index == 0 && have_CRAN()) 
    2L 
  else 
    as.integer(page_index + 1)
}, data=NULL)


###################################################
### code chunk number 145: ex-RGtk2-install-wizard.Rnw:78-80
###################################################
CRAN_package <- NA
install_options <- list() #type, dependencies, lib


###################################################
### code chunk number 146: HelperFunctions
###################################################
## Helper functions
##' return value or NA
##'
gtkTreeViewGetSelectedValue <- function(object, column) {
  cur <- object$getSelection()$getSelected()
  if(cur$retval)
    with(cur, object$getModel()$getValue(iter, column -1 )$value)
  else
    NA
}


have_CRAN <- function() getOption("repos")["CRAN"] != "@CRAN@"

##' from getCRANmirrors
set_CRAN <- function(url) {
  if(is.null(url)) return()
  repos <- getOption("repos")
  repos["CRAN"] <- gsub("/$", "", url)
  options(repos=repos)
}


###################################################
### code chunk number 147: page1
###################################################
populate_page[[1]] <- function() {
  assistant$setPageTitle(pages[[1]], "Install a CRAN package")
  pages[[1]]$packStart(label <- gtkLabel())
  pages[[1]]$packStart(gtkLabel(), expand=TRUE) # a spring
  
  label$setMarkup(paste(
       "<span font='x-large'>Install a CRAN package</span>",
       "This wizard will help install a package from",
       "<b>CRAN</b>. If you have not already specified a",
       "CRAN repository, you will be prompted to do so.",
       sep="\n"))
  assistant$setPageComplete(pages[[1]], TRUE)
}


###################################################
### code chunk number 148: CRANMirror
###################################################
## Not shown
populate_page[[2]] <- function() {
  assistant$setPageTitle(pages[[2]], "Select a CRAN mirror")

  CRAN_mirrors <- getCRANmirrors(all = FALSE, local.only = FALSE)[, c(1,2,4)]
  nms <- names(CRAN_mirrors)
  d <- rGtkDataFrame(CRAN_mirrors)
  #
  view <- gtkTreeView()
  mapply(view$insertColumnWithAttributes, -1, nms[1:2], 
         list(gtkCellRendererText()), text = 0:1)
  view$setModel(d)
  view$getSelection()$unselectAll()     # no selection
  gSignalConnect(view$getSelection(), "changed", function(view, ...) {
    CRAN_repos <- view$getSelectedValue(3)
    set_CRAN(CRAN_repos)
    assistant$setPageComplete(pages[[2]], TRUE)
  }, data=view, user.data.first=TRUE)
  
  
  sw <- gtkScrolledWindow(); sw$add(view)
  sw$setPolicy("automatic", "automatic")
  
  pages[[2]]$packStart(gtkLabel("Select a CRAN mirror"), expand=FALSE)
  pages[[2]]$packStart(sw, expand=TRUE, fill=TRUE)

}


###################################################
### code chunk number 149: SelectPacakge
###################################################
## Not shown
populate_page[[3]] <- function() {
  assistant$setPageTitle(pages[[3]], "Select a CRAN package")
  #
  avail_packages <- available.packages()[, c(1,2)]
  nms <- colnames(avail_packages)
  avail_packages_store <- rGtkDataFrame(avail_packages)
  #
  view <- gtkTreeView()
  mapply(view$insertColumnWithAttributes, -1, nms, 
         list(gtkCellRendererText()), text = 0:1)
  view$setModel(avail_packages_store)
  view$getSelection()$unselectAll()     # no selection
  gSignalConnect(view$getSelection(), "changed", function(view, ...) {
    CRAN_package <<- view$getSelectedValue(1)
    assistant$setPageComplete(pages[[3]], TRUE)
  }, data=view, user.data.first=TRUE) 
  #
  sw <- gtkScrolledWindow(); sw$add(view)
  sw$setPolicy("automatic", "automatic")
  #
  pages[[3]]$packStart(gtkLabel("Select a package to install"), expand=FALSE)
  pages[[3]]$packStart(sw, expand=TRUE, fill=TRUE)
}


###################################################
### code chunk number 150: ex-RGtk2-install-wizard.Rnw:193-262
###################################################
populate_page[[4]] <- function() {
  assistant$setPageTitle(pages[[4]], "Install a CRAN package")
  ##
  get_desc <- function(pkgname) {
    o <- "http://cran.r-project.org/web/packages/%s/%s"
    x <- readLines(sprintf(o, pkgname, "DESCRIPTION"))
    f <- tempfile(); cat(paste(x, collapse="\n"), file=f)
    read.dcf(f)
  }
  desc <- get_desc(CRAN_package)
  #
  label <- gtkLabel()
  label$setLineWrap(TRUE)
  label$setWidthChars(40)
  label$setMarkup(paste(
    sprintf("Install package: <b>%s</b>", desc[1,'Package']),
    "\n",
    sprintf("%s", gsub("\\n", " ", desc[1,'Description'])),
    sep="\n"))
  
  pages[[4]]$packStart(label)
  ##
  table <- gtkTable()
  pages[[4]]$packStart(table, expand=FALSE)
  pages[[4]]$packStart(gtkLabel(), expand=TRUE)
  
  ##
  combo <- gtkComboBoxNewText()
  pkg_types <- c("source", "mac.binary", "mac.binary.leopard",
                 "win.binary", "win64.binary")
  sapply(pkg_types, combo$appendText)
  combo$setActive(which(getOption("pkgType") == pkg_types)-1)
  gSignalConnect(combo, "changed", function(combo, ...) {
    cur <- 1L + combo$getActive()
    install_options[['type']] <<- pkg_types[cur]
  })
  table$attachDefaults(gtkLabel("Package type:"), 0, 1, 0, 1)
  table$attachDefaults(combo, 1, 2, 0, 1)

  ##
  checkButton <- gtkCheckButton()
  checkButton$setActive(TRUE)
  gSignalConnect(checkButton, "toggled", function(ck_btn) {
    install_options$dependencies <<- ck_btn$getActive()
  })
  table$attachDefaults(gtkLabel("Install dependencies"),
                       0, 1, 1, 2)
  table$attachDefaults(checkButton, 1, 2, 1, 2)

  ##
  file_chooser <- gtkFileChooserButton("Select directory...", 
                                      "select-folder")
  file_chooser$setFilename(.libPaths()[1])
  gSignalConnect(file_chooser, "selection-changed", 
                 function(file_chooser) {
                   dir <- file_chooser$getFilename()
                   install_options[['lib']] <<- dir
                 })
  table$attachDefaults(gtkLabel("Where"), 0, 1, 2, 3)
  table$attachDefaults(file_chooser, 1, 2, 2, 3)
  ## align labels to right and set spacing
  sapply(table$getChildren(), function(child) {
    widget <- child$getWidget()
    if(is(widget, "GtkLabel"))  widget['xalign'] <- 1
  })
  table$setColSpacing(0L, 5L)
  ##
  assistant$setPageComplete(pages[[4]], TRUE)
}


###################################################
### code chunk number 151: ex-RGtk2-install-wizard.Rnw:271-290
###################################################
populate_page[[5]] <- function() {
  assistant$setPageTitle(pages[[5]], "Done")
  install_options$pkgs <- CRAN_package
  out <- try(do.call("install.packages", install_options), 
             silent=TRUE)

  label <- gtkLabel(); pages[[5]]$packStart(label)
  if(!inherits(out, "try-error")) {
    label$setMarkup(sprintf("Package %s was installed.", 
                            CRAN_package))
  } else {
    label$setMarkup(paste(sprintf("Package %s, failed install", 
                                  CRAN_package),
                          paste(out, collapse="\n"),
                          sep="\n"))
  }

  assistant$setPageComplete(pages[[5]], FALSE)
}


###################################################
### code chunk number 152: showAssistant
###################################################
populate_page[[1]]()
assistant$show()


###################################################
### code chunk number 153: gtk-cairo-device
###################################################
library(cairoDevice)
device <- gtkDrawingArea()
asCairoDevice(device)
##
window <- gtkWindow(show=FALSE)
window$add(device)
window$showAll()
plot(mpg ~ hp, data = mtcars)


###################################################
### code chunk number 154: gtk-cairo-print-operation
###################################################
print_op <- gtkPrintOperation()


###################################################
### code chunk number 155: gtk-cairo-draw-page
###################################################
gSignalConnect(print_op, "draw-page", 
               function(print_op, context, page_nr) {
                 asCairoDevice(context)
                 plot(mpg ~ wt, data = mtcars)
               })


###################################################
### code chunk number 156: gtk-cairo-run-dialog
###################################################
print_op$run(action = "print-dialog", parent = NULL)


###################################################
### code chunk number 157: ex-RGtk2-manipulate.Rnw:1-27
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
### code chunk number 158: resolveVariableArguments
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
### code chunk number 159: manipulate_example (eval = FALSE)
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
### code chunk number 160: ManipulateClass
###################################################
Manipulate <- setRefClass("Manipulate",
                          fields=list(
                            .code="ANY",
                            .controls="list"
                            ))


###################################################
### code chunk number 161: manipulate_validate_controls
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
### code chunk number 162: Manipulate_change_handler
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
### code chunk number 163: Manipulate_execute
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
### code chunk number 164: Manipulate_Initialize
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
### code chunk number 165: manipulate_constructor
###################################################
manipulate <- function(`_expr`,...) {
  manip <- Manipulate$new(substitute(`_expr`),...)
  manip$execute()
}


###################################################
### code chunk number 166: ManipulateControls
###################################################
ManipulateControls <- setRefClass("ManipulateControls",
                        fields=list(
                          l="list",
                          widget = "ANY",
                          label="ANY",
                          initial="ANY"
                          ))


###################################################
### code chunk number 167: MC_Interface
###################################################
ManipulateControls$methods(
            validate_inputs = function(...) {
              "Validate input code"
            },
            get_value = function(...) {
              "Get value of widget"
            })


###################################################
### code chunk number 168: MC_make_gui
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
### code chunk number 169: Slider_constructor
###################################################
slider <- function(min, max, initial = min, label=NULL, 
                   step = -1, ticks = TRUE) {
  Slider$new(min, max, initial = initial, label = label, 
             step = step, ticks = ticks)
}


###################################################
### code chunk number 170: Slider
###################################################
Slider <- setRefClass("Slider",
                      contains = "ManipulateControls")


###################################################
### code chunk number 171: Slider_validate
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
### code chunk number 172: Slider_initialize
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
### code chunk number 173: Slider_make_gui
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
### code chunk number 174: Picker
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
### code chunk number 175: Checkbox
###################################################
Checkbox <- setRefClass("Checkbox", contains="ManipulateControls")
Checkbox$methods(validate_inputs=function(initial, label) {
                   if ( !is.logical(initial) )
                     stop("initial must be a logical")
                 })


###################################################
### code chunk number 176: ex-RGtk2-manipulate.Rnw:437-453
###################################################
Checkbox$methods(
         initialize = function(initial=FALSE, label=  NULL) {
           validate_inputs(initial, label)
           checkbox <- list(type = 2)
           initFields(l = checkbox, label = label, 
                      initial = initial)
           callSuper()
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
### code chunk number 177: Checkbox_constructor
###################################################
checkbox <- function(initial = FALSE, label = NULL) Checkbox$new(initial, label)                            


###################################################
### code chunk number 178: ex-RGtk2-manipulate.Rnw:464-465
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


###################################################
### code chunk number 179: BasicComponents.Rnw:861-868
###################################################
TARGET.TYPE.TEXT   <- 80                # our enumeration
TARGET.TYPE.PIXMAP <- 81                  
widgetTargetTypes <- 
  list(text = gtkTargetEntry("text/plain", 0, 
         TARGET.TYPE.TEXT),
       pixmap = gtkTargetEntry("image/x-pixmap", 0, 
         TARGET.TYPE.PIXMAP))


###################################################
### code chunk number 180: BasicComponents.Rnw:878-886
###################################################
window <- gtkWindow(); window['title'] <- "Drag Source"
drag_source_widget <-  gtkButton("Text to drag")
window$add(drag_source_widget)

gtkDragSourceSet(drag_source_widget,
       start.button.mask=c("button1-mask", "button3-mask"),
       targets=widgetTargetTypes[["text"]],
       actions="copy")


###################################################
### code chunk number 181: BasicComponents.Rnw:899-903
###################################################
gSignalConnect(drag_source_widget, "drag-data-get", 
               function(widget, context, sel, tType, eTime) {
                 sel$setText(widget$getLabel()) 
               })


###################################################
### code chunk number 182: BasicComponents.Rnw:914-922
###################################################
window <- gtkWindow(); window['title'] <- "Drop Target"
drop_target_widget <- gtkButton("Drop here")
window$add(drop_target_widget)

gtkDragDestSet(drop_target_widget,
               flags="all", 
               targets=widgetTargetTypes[["text"]],
               actions="copy")


###################################################
### code chunk number 183: BasicComponents.Rnw:940-945
###################################################
gSignalConnect(drop_target_widget, "drag-data-received", 
       function(widget, context, x, y, sel, tType, eTime) {
         dropdata <- sel$getText()
         widget$setLabel(rawToChar(dropdata))
       })


###################################################
### code chunk number 184: WidgetsWithModels.Rnw:1-2
###################################################
library(RGtk2)


###################################################
### code chunk number 185: WidgetsWithModels.Rnw:84-88
###################################################
data(Cars93, package="MASS")             # mix of classes
model <- rGtkDataFrame(Cars93)
model[1, 4] <- 12
model[1, 4]                              # get value


###################################################
### code chunk number 186: WidgetsWithModels.Rnw:105-106
###################################################
model$setFrame(Cars93[1:5, 1:5])


###################################################
### code chunk number 187: rgtk2-mvc-treeview-construc
###################################################
view <- gtkTreeView(model)


###################################################
### code chunk number 188: rgtk2-mvc-insert-column-hardway
###################################################
column <- gtkTreeViewColumn()
column$setTitle("Manufacturer")
cell_renderer <- gtkCellRendererText()
column$packStart(cell_renderer)
column$addAttribute(cell_renderer, "text", 0)
view$insertColumn(column, 0)


###################################################
### code chunk number 189: rgtk2-mvc-insert-column-easyway
###################################################
view$insertColumnWithAttributes(position = -1, 
                                title = "Model", 
                                cell = gtkCellRendererText(), 
                                text = 2 - 1) # second column


###################################################
### code chunk number 190: rgtk2-mvc-insert-all-columns
###################################################
view <- gtkTreeView(model)
mapply(view$insertColumnWithAttributes,  
       position = -1, 
       title = colnames(model), 
       cell = list(gtkCellRendererText()), 
       text = seq_len(ncol(model)) - 1
       )


###################################################
### code chunk number 191: scrollView
###################################################
window <- gtkWindow()
window$setTitle("Tabular view of data frame")
scrolled_window <- gtkScrolledWindow()
window$add(scrolled_window)
scrolled_window$add(view)


###################################################
### code chunk number 192: rgtk2-mvc-path-constructor-list
###################################################
second_row <- gtkTreePathNewFromIndices(1)


###################################################
### code chunk number 193: rgtk2-mvc-path-constructor-tree
###################################################
abc_path <- gtkTreePathNewFromIndices(c(0, 2, 1))
abc_path <- gtkTreePathNewFromString("0:2:1")


###################################################
### code chunk number 194: rgtk2-mvc-iter-traverse
###################################################
iter <- model$getIterFirst()
manufacturer <- character()
while(iter$retval) {
  manufacturer <- c(manufacturer, model$get(iter$iter,0)[[1]])
  iter$retval <- model$iterNext(iter$iter)
}


###################################################
### code chunk number 195: rgtk2-mvc-iter-apply
###################################################
nrows <- model$iterNChildren(NULL)
manufacturer <- sapply(seq(nrows) - 1L, function(i) {
  iter <- model$iterNthChild(NULL, i)
  model$get(iter$iter, 0)[[1]]
})


###################################################
### code chunk number 196: WidgetsWithModels.Rnw:403-407
###################################################
model <- rGtkDataFrame(mtcars)
view <- gtkTreeView(model)
selection <- view$getSelection()
selection$setMode("single")


###################################################
### code chunk number 197: WidgetsWithModels.Rnw:415-428
###################################################
column <- gtkTreeViewColumn()
view$insertColumnWithAttributes(0, "title", gtkCellRendererText(), text = 0)
## pack in GUI
scrolled_window <- gtkScrolledWindow()
scrolled_window$add(view)
##
window <- gtkWindow(show=FALSE)
window['title'] <- "Multiple selection example"
window$add(scrolled_window)
window$show()
## some selection
selection$selectPath(gtkTreePathNewFromIndices(3)) # set 
# 


###################################################
### code chunk number 198: WidgetsWithModels.Rnw:433-435
###################################################
selected <- selection$getSelected()
with(selected, model$getValue(iter, 0)$value)


###################################################
### code chunk number 199: WidgetsWithModels.Rnw:447-455
###################################################
gSignalConnect(selection, "changed", function(selection) {
  selected_rows <- selection$getSelectedRows()
  if(length(selected_rows$retval)) {
    rows <- sapply(selected_rows$retval, 
                   gtkTreePathGetIndices) + 1L
    selected_rows$model[rows, 1]
  }
})


###################################################
### code chunk number 200: WidgetsWithModels.Rnw:467-468 (eval = FALSE)
###################################################
## sapply(view$getColumns(), function(i) i == column)


###################################################
### code chunk number 201: rgtk2-mvc-sorting-clickable
###################################################
column <- view$getColumn(0)
column$setSortColumnId(0)


###################################################
### code chunk number 202: rgtk2-mvc-sorting-sortable
###################################################
model$setSortColumnId(0, "ascending")


###################################################
### code chunk number 203: WidgetsWithModels.Rnw:504-505
###################################################
require(MASS)


###################################################
### code chunk number 204: basicSort
###################################################
model <- rGtkDataFrame(Cars93)
sorted_model <- gtkTreeModelSortNewWithModel(model)
view <- gtkTreeView(sorted_model)
mapply(view$insertColumnWithAttributes,
       position = -1,
       title = colnames(model),
       cell = list(gtkCellRendererText()),
       text = seq_len(ncol(model)) - 1)
sapply(seq_len(ncol(model)), function(i)
       view$getColumn(i - 1)$setSortColumnId(i - 1))


###################################################
### code chunk number 205: sort-example
###################################################
f <- function(model, iter1, iter2, user.data) {
  types <- c("Compact", "Small", "Sporty", "Midsize", 
             "Large", "Van")
  column <- user.data
  val1 <- model$getValue(iter1, column)$value
  val2 <- model$getValue(iter2, column)$value
  as.integer(match(val1, types) - match(val2, types))
}
sorted_model$setSortFunc(sort.column.id = 3 - 1, sort.func=f, 
                         user.data = 3 - 1)


###################################################
### code chunk number 206: notShown
###################################################
## basic GUI
sw <- gtkScrolledWindow()
sw$add(view)
w <- gtkWindow(show=FALSE)
w['title'] <- "Example of sortable treeview"
w$add(sw)
w$show()


###################################################
### code chunk number 207: WidgetsWithModels.Rnw:571-578
###################################################
DF <- Cars93
model <- rGtkDataFrame(cbind(DF, .vis=rep(TRUE, nrow(DF))))
filtered_model <- model$filter()
filtered_model$setVisibleColumn(length(DF))          # 0-based
view <- gtkTreeView(filtered_model)
## Adjust filter
model[,".vis"] <- DF$MPG.highway >= 30


###################################################
### code chunk number 208: notShown
###################################################
mapply(view$insertColumnWithAttributes,  
       position=-1, 
       title=colnames(DF), 
       cell=list(gtkCellRendererText()), 
       text = seq_len(length(DF)) - 1
       )    
##
sw <- gtkScrolledWindow()
sw$add(view)
w <- gtkWindow(show=FALSE)
w$add(sw)
w$show()


###################################################
### code chunk number 209: ex-RGtk2-filtered.Rnw:5-6
###################################################
library(RGtk2)


###################################################
### code chunk number 210: ex-RGtk2-filtered.Rnw:20-23
###################################################
DF <- data.frame(state.name)
DF$visible <- rep(TRUE, nrow(DF))
model <- rGtkDataFrame(DF)


###################################################
### code chunk number 211: ex-RGtk2-filtered.Rnw:28-31
###################################################
filtered_model <- model$filter()
filtered_model$setVisibleColumn(ncol(DF) - 1)      # offset
view <- gtkTreeView(filtered_model)


###################################################
### code chunk number 212: ex-RGtk2-filtered.Rnw:35-37
###################################################
view$insertColumnWithAttributes(0, "Col", 
                 gtkCellRendererText(), text = 0)


###################################################
### code chunk number 213: ex-RGtk2-filtered.Rnw:45-52
###################################################
entry <- gtkEntry()
gSignalConnect(entry, "changed", function(entry, user.data) {
  pattern <- entry$getText()
  DF <- user.data$getModel()
  values <- DF[, "state.name"]
  DF[, "visible"] <- grepl(pattern, values)
}, data=filtered_model)


###################################################
### code chunk number 214: ex-RGtk2-filtered.Rnw:58-74
###################################################
## not shown, but this places widgets into a simple GUI
window <- gtkWindow(show=FALSE)
window['title'] <- "A filtered data model"
window$setSizeRequest(width=300, height=400)

box <- gtkVBox()
window$add(box)
box$packStart(entry, expand=FALSE)

## add scroll window
sw <- gtkScrolledWindow()
sw$setPolicy("automatic", "automatic")
sw$add(view)
box$packStart(sw, expand=TRUE, fill=TRUE)

window$show()


###################################################
### code chunk number 215: gtk-mvc-cell-explicit
###################################################
cell_renderer <- gtkCellRendererText()
cell_renderer['cell-background'] <- "gray"


###################################################
### code chunk number 216: cr-right-aligned
###################################################
cell_renderer <- gtkCellRendererText()
cell_renderer['xalign'] <- 1          # default 0.5 = centered
cell_renderer['family'] <- "Helvetica"  


###################################################
### code chunk number 217: WidgetsWithModels.Rnw:696-698
###################################################
cell_renderer <- gtkCellRendererText()
store <- model


###################################################
### code chunk number 218: editedSignal
###################################################
cell_renderer['editable'] <- TRUE
gSignalConnect(cell_renderer, "edited", 
       f=function(cell_renderer, path, newtext, user.data) {
         i <- as.numeric(path) + 1
         j <- user.data$column
         model <- user.data$model
         model[i, j] <- newtext
       }, data=list(model=store, column=1))


###################################################
### code chunk number 219: editableTableForCollectingOptions
###################################################
## GUI for configuring options -- in a table
library(RGtk2)


###################################################
### code chunk number 220: ex-RGtk2-options-in-table.Rnw:21-28
###################################################
opts <- c("main", "sub", "xlab", "ylab", "line", "outer")
DF <- data.frame(option = opts,
           value = c("", "", "", "", "0", "FALSE"),
           class = c(rep("character",4),"integer", "logical"),
           edit_color = rep("gray95", 6),
           dirty = rep(FALSE, 6),
           stringsAsFactors = FALSE)


###################################################
### code chunk number 221: model
###################################################
model <- rGtkDataFrame(DF)
view <- gtkTreeView(model)
##
cell_renderer <- gtkCellRendererText()
cell_renderer['background'] <- 'gray80'
view$insertColumnWithAttributes(position = -1,
                                title = "Option",
                                cell = cell_renderer,
                                text = 1 - 1)


###################################################
### code chunk number 222: secondColumn
###################################################
cell_renderer <- gtkCellRendererText()
cell_renderer['editable'] <- TRUE
view$insertColumnWithAttributes(position = -1,
                                title = "Value",
                                cell = cell_renderer,
                                text = 2 - 1,
                                background = 4 - 1
                                )


###################################################
### code chunk number 223: editConnect
###################################################
gSignalConnect(cell_renderer, "edited", 
    function(cell_renderer, path, new.text, user.data) {
      model <- user.data$model
      i <- as.numeric(path) + 1; j <- user.data$column
      val <- as(new.text, model[i, 'class'])
      model[i,j] <- as(val, "character")   
      model[i, 'dirty'] <- TRUE                 # mark dirty
      model[i, 'edit_color'] <- 'gray70'        # change color
    }, data=list(model=model, column=2))


###################################################
### code chunk number 224: ex-RGtk2-options-in-table.Rnw:85-92
###################################################
window <- gtkWindow(show=FALSE)
window['title'] <- "Option editor"
window$setSizeRequest(300,500)
scrolled_window <- gtkScrolledWindow()
window$add(scrolled_window)
scrolled_window$add(view)
window$show()


###################################################
### code chunk number 225: ex-RGtk2-options-in-table.Rnw:114-120
###################################################
require(helpr, quietly=TRUE)
package <- "graphics"; topic <- "title"
rd <- helpr:::parse_help(helpr:::pkg_topic(package, topic), 
                         package = package)
descs <- rd$params$args
names(descs) <- sapply(descs, function(i) i$param)


###################################################
### code chunk number 226: ex-RGtk2-options-in-table.Rnw:129-143
###################################################
view["has-tooltip"] <- TRUE
gSignalConnect(view, "query-tooltip", 
       function(view, x, y, key_mode, tooltip, user.data) {
         out <- view$getTooltipContext(x, y, key_mode)
         if(out$retval) {
           model <- view$getModel()
           i <- as.numeric(out$path$toString()) + 1
           val <- model[i, "option"]
           txt <- descs[[val]]$desc
           txt <- gsub("code>","b>", txt)  # no code in Pango
           tooltip$setMarkup(txt)
         }
         out$retval
       })


###################################################
### code chunk number 227: WidgetsWithModels.Rnw:742-747
###################################################
cell_renderer <- gtkCellRendererCombo()
model <- rGtkDataFrame(state.name)
cell_renderer['model'] <- model
cell_renderer['text-column'] <- 0
cell_renderer['editable'] <- TRUE                  # needed


###################################################
### code chunk number 228: VariableSelectionExample
###################################################
## Example showing implementation of variable selection widget where two tables show possible selections
## and selection. Similar to SPSS widget
## Illustrates filtered models, icons in view column
library(RGtk2)


###################################################
### code chunk number 229: ex-RGtk2-select-variables.Rnw:24-25
###################################################
DF <- get(data(Cars93, package="MASS"))


###################################################
### code chunk number 230: ex-RGtk2-select-variables.Rnw:41-43
###################################################
library(ProgGUIinR)                     # for make_icon
#source("../ProgGUIInR/R/misc.R")     # for make_icon


###################################################
### code chunk number 231: make_icon
###################################################
make_icon_pixmap <- function(x, ...) {
  require(grid); require(cairoDevice)
  pixmap <- gdkPixmap(drawable = NULL, width = 16, height=16, 
                      depth = 24)
  asCairoDevice(pixmap)
  grid.newpage()
  grid.draw(make_icon(x))
  dev.off()
  gdkPixbufGetFromDrawable(NULL, pixmap, NULL, 0,0,0,0,-1,-1)
}


###################################################
### code chunk number 232: model
###################################################
model_df <- data.frame(Variables = I(sort(names(DF))),
                       icon = I(sapply(DF, make_icon_pixmap)),
                       selected = rep(FALSE, ncol(DF)))
model <- rGtkDataFrame(model_df)


###################################################
### code chunk number 233: filterModels
###################################################
selected_filter <- model$filter()
selected_filter$setVisibleColumn(2)
unselected_filter <- model$filter()
unselected_filter$setVisibleFunc(function(model, iter) {
  !model$get(iter, 2)[[1]]
})


###################################################
### code chunk number 234: views
###################################################
views <- list()
views$unselected_view <- gtkTreeView(unselected_filter)
views$selected_view <- gtkTreeView(selected_filter)
##
sapply(views, function(view) {
  selection <- view$getSelection()
  selection$setMode('multiple')
})


###################################################
### code chunk number 235: viewColumns
###################################################
make_view_column <- function() {
  column <- gtkTreeViewColumn()
  column$setTitle("Variable")
  column$packStart(cell_renderer <- gtkCellRendererPixbuf())
  column$addAttribute(cell_renderer, "pixbuf", 1L)
  column$packStart(cell_renderer <- gtkCellRendererText())
  column$addAttribute(cell_renderer, "text", 0L)
  column
}
sapply(views, function(view) 
       view$insertColumn(make_view_column(), 0))


###################################################
### code chunk number 236: extendAPI
###################################################
## add to the gtkTreeView API for convenience
gtkTreeViewSelectedIndices <- function(object) {
  model <- object$getModel()          # Filtered!
  paths <- object$getSelection()$getSelectedRows()$retval
  path_strings <- sapply(paths, function(i) {
    model$convertPathToChildPath(i)$toString()
  })
  if(length(path_strings) == 0)
    integer(0)
  else
    as.numeric(path_strings) + 1 # 1-based
}
## does object have selection?
gtkTreeViewHasSelection <-
  function(obj) length(obj$selectedIndices()) > 0


###################################################
### code chunk number 237: buttons
###################################################
buttons <- list()
buttons$unselect_button <- gtkButton("<")
buttons$select_button <- gtkButton(">")
toggleSelectionOnClick <- function(button, view) {
  gSignalConnect(button, "clicked", function(button) {
    message("clicked")
    ind <- view$selectedIndices()
    model[ind, "selected"] <- !model[ind, "selected"]
  })
}
sapply(1:2, function(i) toggleSelectionOnClick(buttons[[i]], 
                                               views[[3-i]]))


###################################################
### code chunk number 238: sensitiveButtons
###################################################
sapply(buttons, gtkWidgetSetSensitive, FALSE)
trackSelection <- function(button, view) {
  gSignalConnect(view$getSelection(), "changed", 
     function(x) button['sensitive'] <- view$hasSelection())
}
sapply(1:2, function(i) trackSelection(buttons[[i]], 
                                       views[[3-i]]))


###################################################
### code chunk number 239: guiLayout
###################################################
window <- gtkWindow(show=FALSE)
window$setTitle("Select variables example")
window$setDefaultSize(600, 400)
hbox <- gtkHBox()
window$add(hbox)
## scrollwindows
scrolls <- list()
scrolls$unselected_scroll <- gtkScrolledWindow()
scrolls$selected_scroll <- gtkScrolledWindow()
mapply(gtkContainerAdd, object = scrolls, widget = views)
mapply(gtkScrolledWindowSetPolicy, scrolls, 
       "automatic", "automatic")
## buttons
button_box <- gtkVBox()
centered_box <- gtkVBox()
button_box$packStart(centered_box, expand=TRUE, fill = FALSE)
centered_box$setSpacing(12)
sapply(buttons, centered_box$packStart, expand = FALSE)
##
hbox$packStart(scrolls$unselected_scroll, expand = TRUE)
hbox$packStart(button_box, expand = FALSE)
hbox$packStart(scrolls$selected_scroll, expand = TRUE)


###################################################
### code chunk number 240: packButtons
###################################################
window$show()


###################################################
### code chunk number 241: cellRendererToggle
###################################################
cell_renderer <- gtkCellRendererToggle()
cell_renderer['activatable'] <- TRUE   # cell can be activated
cell_renderer['active'] <- TRUE
gSignalConnect(cell_renderer, "toggled", function(w, path) {
  model$active[as.numeric(path) + 1] <- w['active']
})


###################################################
### code chunk number 242: ex-RGtk2-add-toggle-to-df.Rnw:11-13
###################################################
## example showing how to add a toggle button on left of data display
library(RGtk2)


###################################################
### code chunk number 243: FixACRANforSweave
###################################################
repos <- getOption("repos")
repos["CRAN"] <- "http://streaming.stat.iastate.edu/CRAN"
options(repos = repos)


###################################################
### code chunk number 244: getUpgradablePackages
###################################################
old_packages <- 
  old.packages()[,c("Package", "Installed", "ReposVer")]
DF <- as.data.frame(old_packages)


###################################################
### code chunk number 245: ex-RGtk2-add-toggle-to-df.Rnw:33-35
###################################################
doUpdate <- function(old_packages) 
  install.packages(old_packages$Package)


###################################################
### code chunk number 246: ex-RGtk2-add-toggle-to-df.Rnw:47-48
###################################################
model <- rGtkDataFrame(cbind(DF, .toggle=rep(FALSE, nrow(DF))))


###################################################
### code chunk number 247: ex-RGtk2-add-toggle-to-df.Rnw:53-66
###################################################
view <- gtkTreeView()
cell_renderer <- gtkCellRendererToggle() # add toggle
view$insertColumnWithAttributes(0, "", cell_renderer, 
                                active = ncol(DF))
cell_renderer['activatable'] <- TRUE
gSignalConnect(cell_renderer, "toggled", 
               function(cell_renderer, path, user.data) {
                 view <- user.data
                 row <- as.numeric(path) + 1
                 model <- view$getModel()
                 n <- dim(model)[2]
                 model[row, n] <- !model[row, n]
               }, data=view)


###################################################
### code chunk number 248: ex-RGtk2-add-toggle-to-df.Rnw:70-72 (eval = FALSE)
###################################################
## mapply(view$insertColumnWithAttributes, -1, colnames(DF), 
##        list(gtkCellRendererText()), text = seq_along(DF) -1L)


###################################################
### code chunk number 249: ex-RGtk2-add-toggle-to-df.Rnw:76-77
###################################################
view$setModel(model)


###################################################
### code chunk number 250: ex-RGtk2-add-toggle-to-df.Rnw:86-94
###################################################
button <- gtkButton("Update packages")
gSignalConnect(button, "clicked", function(button, data) {
  view <- data
  model <- view$getModel()
  old_packages <- 
    model[model[, ncol(model)], -ncol(model), drop = FALSE]
  doUpdate(old_packages)
}, data=view)


###################################################
### code chunk number 251: ex-RGtk2-add-toggle-to-df.Rnw:100-112
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Installed packages that need upgrading")
window$setSizeRequest(300, 300)

vbox <- gtkVBox(); window$add(vbox)
scrolled_window <- gtkScrolledWindow()
vbox$packStart(scrolled_window, expand = TRUE, fill = TRUE)

scrolled_window$add(view)
scrolled_window$setPolicy("automatic", "automatic")
vbox$packStart(button, expand = FALSE)
window$show()


###################################################
### code chunk number 252: comboEditor
###################################################
cell_renderer <- gtkCellRendererProgress()
cell_renderer["value"] <- 50


###################################################
### code chunk number 253: WidgetsWithModels.Rnw:867-873
###################################################
func <- function(column, cell_renderer, model, iter, data) {
  val <- model$getValue(iter, 0)$value
  f_val <- sprintf("%.3f", val)
  cell_renderer['text'] <- f_val
  cell_renderer['xalign'] <- 1
}


###################################################
### code chunk number 254: WidgetsWithModels.Rnw:879-885
###################################################
view <- gtkTreeView(rGtkDataFrame(data.frame(rnorm(100))))
cell_renderer <- gtkCellRendererText()
view$insertColumnWithAttributes(0, "numbers", cell_renderer, 
                                text = 0)
column <- view$getColumn(0)
column$setCellDataFunc(cell_renderer, func)


###################################################
### code chunk number 255: WidgetsWithModels.Rnw:940-952
###################################################
model <- gtkTreeStore("gchararray")
by(Cars93, Cars93$Manufacturer, function(DF) {
  parent_iter <- model$append()
  model$setValue(parent_iter$iter, column = 0, value = 
                 DF$Manufacturer[1])
  sapply(DF$Model, function(car_model) {
    child_iter <- model$append(parent = parent_iter$iter)
    if (is.null(child_iter$retval)) 
      model$setValue(child_iter$iter, column = 0, 
                     value = car_model)
  })
})


###################################################
### code chunk number 256: WidgetsWithModels.Rnw:957-959
###################################################
iter <- model$getIterFromString("0:0")
model$getValue(iter$iter, column = 0)$value


###################################################
### code chunk number 257: rgtk2-mvc-tree-traverse
###################################################
iter <- model$getIterFirst()
values <- NULL
while(iter$retval) {
  child_iter <- model$iterChildren(iter$iter)
  while(child_iter$retval) {
    values <- c(values, model$get(child_iter$iter, 0)[[1]])
    child_iter$retval <- model$iterNext(child_iter$iter)
  }
  iter$retval <- model$iterNext(iter$iter)
}


###################################################
### code chunk number 258: notShown
###################################################
## define tstore, but aslo in earlier example so not shown
data(Cars93, package="MASS")
model <- gtkTreeStore("gchararray")
Manufacturers <- Cars93$Manufacturer
Makes <- split(Cars93[,"Model"], Manufacturers)
for(i in unique(Manufacturers)) {
  piter <- model$append()              # parent
  model$setValue(piter$iter, column=0, value=i)
  for(j in Makes[[i]]) { 
    sibiter <- model$append(parent=piter$iter) # child
    if(is.null(sibiter$retval)) 
      model$setValue(sibiter$iter,column=0, value=j)
  }
}


###################################################
### code chunk number 259: makeView
###################################################
view <- gtkTreeView()
view$insertColumnWithAttributes(0, "Make", 
           gtkCellRendererText(), text = 0)
view$setModel(model)


###################################################
### code chunk number 260: makeGUI
###################################################
w <- gtkWindow(show=FALSE)
w['title'] <- "Example of changing models"
sw <- gtkScrolledWindow()
sw$add(view)
w$add(sw)
w$show()


###################################################
### code chunk number 261: ex-RGtk2-simple-tree.Rnw:45-47
###################################################
model <- rGtkDataFrame(Cars93[,"Model", drop=FALSE])
view$setModel(model)


###################################################
### code chunk number 262: ex-RGtk2-combobox-entry.Rnw:1-4
###################################################
## a combobox that learns as you go.
## no tooltip per item, but here we add as detail
library(RGtk2)


###################################################
### code chunk number 263: ex-RGtk2-combobox-entry.Rnw:14-18
###################################################
model <- rGtkDataFrame(data.frame(filename = character(0), 
                                  visits = character(0), 
                                  nvisits = integer(0), 
                                  stringsAsFactors = FALSE))


###################################################
### code chunk number 264: ex-RGtk2-combobox-entry.Rnw:32-34
###################################################
combo_box <- gtkComboBoxEntryNewWithModel(model, 
                                          text.column = 0)


###################################################
### code chunk number 265: ConfigureCellRenderers
###################################################
cell_renderer <- gtkCellRendererText()
combo_box$packStart(cell_renderer)
combo_box$addAttribute(cell_renderer, "text", 1)
cell_renderer['foreground'] <- "gray50"
cell_renderer['ellipsize'] <- "end"
cell_renderer['style'] <- "italic"
cell_renderer['alignment'] <- "right"


###################################################
### code chunk number 266: helperFunction2
###################################################
callHelpFunction <- function(combo_box, value) {
  model <- combo_box$getModel()
  ind <- match(value, model[,1,drop=TRUE])
  nvisits <- model[ind, "nvisits"] <- model[ind, "nvisits"]+1
  model[ind, "visits"] <- 
    sprintf(ngettext(nvisits,"%s visit","%s visits"), nvisits)
  ## select for easier editing
  combo_box$getChild()$selectRegion(start = 0, end = -1)
  help(value)
}
gSignalConnect(combo_box, "changed", 
               f = function(combo_box, ...) {
                 if(combo_box$getActive() >= 0) {
                   value <- combo_box$getActiveText()
                   callHelpFunction(combo_box, value)
                 }
               })


###################################################
### code chunk number 267: ex-RGtk2-combobox-entry.Rnw:98-110
###################################################
gSignalConnect(combo_box$getChild(), "activate", 
       f = function(combo_box, entry, ...) {
         value <- entry$getText()
         if(!any(value == combo_box$getModel()[,1])) {
           model <- combo_box$getModel()
           tmp <- data.frame(filename = value, visits = "", 
                             nvisits = 0, 
                             stringsAsFactors = FALSE)
           model$appendRows(tmp)
         }
         callHelpFunction(combo_box, value)
       }, data = combo_box, user.data.first = TRUE)


###################################################
### code chunk number 268: Layout
###################################################
window <- gtkWindow(show = FALSE)
window['border-width'] <- 15
hbox <- gtkHBox(); window$add(hbox)
hbox$packStart(gtkLabel("Help on:"))
hbox$packStart(combo_box, expand = TRUE, fill = TRUE)
#
window$show()


###################################################
### code chunk number 269: ex-RGtk2-entry-completion.Rnw:2-3
###################################################
require(RGtk2)


###################################################
### code chunk number 270: AppendWords
###################################################
entry <- gtkEntry(); completion <- gtkEntryCompletion()
entry$setCompletion(completion)


###################################################
### code chunk number 271: SetCompletion
###################################################
model <- rGtkDataFrame(state.name)
completion$setModel(model)
completion$setTextColumn(0)
completion['inline-completion'] <- TRUE
completion['popup-single-match'] <- FALSE


###################################################
### code chunk number 272: SetMatchFunc
###################################################
matchAnywhere <- function(completion, key, iter, user.data) {
  model <- completion$getModel()
  row_value <- model$getValue(iter, 0)$value
  key <- completion$getEntry()$getText() # case sensitivity
  grepl(key, row_value)
}
completion$setMatchFunc(matchAnywhere)


###################################################
### code chunk number 273: notShown
###################################################
## Our basic GUI is basic:
w <- gtkWindow(show=FALSE)
w$setTitle("Test of entry with completion")
w$add(entry)
w$showAll()


###################################################
### code chunk number 274: gtk-mvc-entry-buffer
###################################################
buffer <- gtkEntryBuffer()        
entry1 <- gtkEntry(buffer = buffer)
entry2 <- gtkEntry(buffer = buffer)
entry1$setText("echo")
entry2$getText()


###################################################
### code chunk number 275: WidgetsWithModels.Rnw:1190-1198
###################################################
view <- gtkTextView()
scrolled_window <- gtkScrolledWindow()
scrolled_window$add(view)
scrolled_window$setPolicy("automatic", "automatic")
##
window <- gtkWindow()
window['border-width'] <- 15
window$add(scrolled_window)


###################################################
### code chunk number 276: setText
###################################################
buffer <- view$getBuffer()
buffer$setText("Lorem ipsum dolor sit amet ...")


###################################################
### code chunk number 277: bufferGetText
###################################################
start <- buffer$getStartIter()$iter    
end <- buffer$getEndIter()$iter
buffer$getText(start, end)


###################################################
### code chunk number 278: gtk-mvc-text-noneditable
###################################################
view['editable'] <- FALSE
view['cursor-visible'] <- FALSE


###################################################
### code chunk number 279: gtk-mvc-buffer-iter-bounds
###################################################
bounds <- buffer$getBounds()
bounds


###################################################
### code chunk number 280: gtk-mvc-buffer-iter-atLineOffset
###################################################
iter <- buffer$getIterAtLineOffset(0, 6)
iter$iter$getChar()                     # unicode, not text


###################################################
### code chunk number 281: gtk-mvc-buffer-iter-getChar
###################################################
bounds$start$getChar()                  # unicode


###################################################
### code chunk number 282: gtk-mvc-buffer-iter-getText
###################################################
bounds$start$getText(bounds$end)


###################################################
### code chunk number 283: gtk-mvc-buffer-iter-insert
###################################################
buffer$insert(bounds$start, "prefix")


###################################################
### code chunk number 284: WidgetsWithModels.Rnw:1392-1396
###################################################
## setup example, not shown
w <- gtkWindow()
view <- gtkTextView()
w$add(view)


###################################################
### code chunk number 285: FindWordAtMouseClick
###################################################
gSignalConnect(view, "button-press-event", 
       f=function(view, event, ...) {
         start <- view$getIterAtLocation(event$getX(), 
                                         event$getY())$iter
         end <- start$copy()
         start$backwardWordStart()
         end$forwardWordEnd()
         val <- start$getText(end)
         print(val)
         return(FALSE) # call next handler
       })


###################################################
### code chunk number 286: gtk-mvc-text-mark-insert
###################################################
insert <- buffer$getMark("insert")


###################################################
### code chunk number 287: gtk-mvc-text-mark-getIter
###################################################
insert_iter <- buffer$getIterAtMark(insert)$iter
bounds$start$getText(insert_iter)


###################################################
### code chunk number 288: gtk-mvc-text-mark-gravity
###################################################
insert_iter$getOffset()
buffer$insert(insert_iter, "at insertion point")
buffer$getIterAtMark(insert)$iter$getOffset()


###################################################
### code chunk number 289: gtk-mvc-text-mark-construct
###################################################
mark <- buffer$createMark(mark.name = "start", 
                          where = buffer$getStartIter()$iter, 
                          left.gravity = TRUE)


###################################################
### code chunk number 290: gtk-mvc-text-tags-create
###################################################
tag_bold <- buffer$createTag(tag.name="bold", 
                             weight=PangoWeight["bold"])
tag_emph <- buffer$createTag(tag.name="emph", 
                             style=PangoStyle["italic"])
tag_large <- buffer$createTag(tag.name="large", 
                              font="Serif normal 18")


###################################################
### code chunk number 291: gtk-mvc-text-tags-apply
###################################################
iter <- buffer$getBounds()
buffer$applyTag(tag_bold, iter$start, iter$end) # iters update
buffer$applyTagByName("emph", iter$start, iter$end)


###################################################
### code chunk number 292: gtk-mvc-text-selectRange
###################################################
start_iter <- buffer$getStartIter()$iter
end_iter <- start_iter$copy(); end_iter$forwardWordEnd()
buffer$selectRange(start_iter, end_iter)


###################################################
### code chunk number 293: gtk-mvc-text-clipboard-get
###################################################
clipboard <- gtkClipboardGet()


###################################################
### code chunk number 294: gtk-mvc-text-clipboard-copy-paste
###################################################
buffer$copyClipboard(clipboard)
buffer$pasteClipboard(clipboard, 
            override.location = buffer$getEndIter()$iter, 
            default.editable = TRUE)


###################################################
### code chunk number 295: gtk-mvc-text-anchor
###################################################
anchor <- buffer$createChildAnchor(buffer$getEndIter()$iter)


###################################################
### code chunk number 296: gtk-mvc-text-addChild
###################################################
button <- gtkButton("click me")
view$addChildAtAnchor(button, anchor)


###################################################
### code chunk number 297: ex-RGtk2-terminal.Rnw:13-15
###################################################
## make a *basic* terminal in RGtk2
library(RGtk2)


###################################################
### code chunk number 298: TextViewWidget
###################################################
view <- gtkTextView()
buffer <- view$getBuffer()
font <- pangoFontDescriptionFromString("Monospace")
view$modifyFont(font)                     # widget wide


###################################################
### code chunk number 299: ex-RGtk2-terminal.Rnw:30-36
###################################################
buffer$createTag(tag.name = "cmdInput")
buffer$createTag(tag.name = "cmdOutput", 
                 weight = PangoWeight["bold"])
buffer$createTag(tag.name = "cmdError", 
       weight = PangoStyle["italic"], foreground = "red")
buffer$createTag(tag.name = "uneditable", editable = FALSE)


###################################################
### code chunk number 300: ex-RGtk2-terminal.Rnw:41-46
###################################################
start_cmd <- buffer$createMark("start_cmd", 
                              buffer$getStartIter()$iter, 
                              left.gravity = TRUE)
bufferEnd <- buffer$createMark("bufferEnd", 
                               buffer$getEndIter()$iter)


###################################################
### code chunk number 301: ex-RGtk2-terminal.Rnw:53-67
###################################################
add_prompt <- function(obj, prompt = c("prompt", "continue"),
                      set_mark = TRUE) 
{
  prompt <- match.arg(prompt)
  prompt <- getOption(prompt)
  
  end_iter <- obj$getEndIter()
  obj$insert(end_iter$iter, prompt)
  if(set_mark)
    obj$moveMarkByName("start_cmd", end_iter$iter)
  obj$applyTagByName("uneditable", obj$getStartIter()$iter, 
                     end_iter$iter)
}
add_prompt(buffer) ## place an initial prompt


###################################################
### code chunk number 302: add_ouput
###################################################
add_ouput <- function(obj, output, tag_name = "cmdOutput") {
  end_iter <- obj$getEndIter()
  if(length(output) > 0)  
    sapply(output, function(i)  {
      obj$insertWithTagsByName(end_iter$iter, i, tag_name)
      obj$insert(end_iter$iter, "\n", len=-1)
    })
}


###################################################
### code chunk number 303: ex-RGtk2-terminal.Rnw:90-98
###################################################
find_cmd <- function(obj) {
  end_iter <- obj$getEndIter()
  start_iter <- obj$getIterAtMark(start_cmd)
  cmd <- obj$getText(start_iter$iter, end_iter$iter, TRUE)
  regex <- paste("\n[", getOption("continue"), "] ", sep = "")
  cmd <- unlist(strsplit(cmd, regex))
  cmd
}


###################################################
### code chunk number 304: evalCmd
###################################################
require(evaluate)
eval_cmd <- function(view, cmd) {
  buffer <- view$getBuffer()
  out <- try(evaluate:::evaluate(cmd, .GlobalEnv), 
             silent = TRUE)

  if(inherits(out, "try-error")) {
    ## parse error
    add_ouput(buffer, out, "cmdError")
  } else if(inherits(out[[2]], "error")) {
    if(grepl("end", out[[2]])) {        # a hack here
      add_prompt(buffer, "continue", set_mark = FALSE)
      return()
    } else {
      add_ouput(buffer, out[[2]]$message, "cmdError")
    }
  } else {
    add_ouput(buffer, out[[2]], "cmdOutput")
  }
  add_prompt(buffer, "prompt", set_mark = TRUE)
}


###################################################
### code chunk number 305: connectBinding
###################################################
gSignalConnect(view, "key-release-event", 
               f=function(view, event) {
                 buffer <- view$getBuffer()
                 keyval <- event$getKeyval()
                 if(keyval == GDK_Return) {
                   cmd <- find_cmd(buffer)
                   if(length(cmd) && nchar(cmd) > 0)
                     eval_cmd(view, cmd)
                 }
               })


###################################################
### code chunk number 306: ex-RGtk2-terminal.Rnw:154-160
###################################################
scroll_viewport <- function(view, ...) {
  view$scrollToMark(bufferEnd, within.margin = 0)
  return(FALSE)
}
gSignalConnect(buffer, "changed", scroll_viewport, data=view, 
               after = TRUE, user.data.first = TRUE)


###################################################
### code chunk number 307: makeGUI
###################################################
## scroll window
sw <- gtkScrolledWindow()
sw$setPolicy("automatic", "automatic")
sw$add(view)

## top-level window
w <- gtkWindow(show=FALSE)
w$setTitle("A terminal")
w$add(sw)
w$setSizeRequest(400,200)
w$showAll()


###################################################
### code chunk number 308: ex-RGtk2-terminal.Rnw:179-251
###################################################
## History features
## This is not illustrated in text, but is added here to illustrate how this might be implemented
## The major issue with this example is we can't trap the return or arrow keys before they move 
## the cursor so any thing ends up looking jerky

## store the stack and a pointer to the current command with the text buffer
buffer$setData("history", list())
buffer$setData("ptr", 0)


## replace cmd with that in str.
replace_cmd <- function(obj, str) {
  end_iter <- obj$getEndIter()
  start_iter <- obj$getIterAtMark(start_cmd)
  obj$delete(start_iter$iter, end_iter$iter)
  end_iter <- obj$getEndIter()
  obj$insertWithTagsByName(end_iter$iter, str[1], "cmdInput")
  if(length(str) > 1) {
    for(i in str[-1]) {
      obj$insert(end_iter$iter, "\n")
      obj$insertWithTagsByName(end_iter$iter, getOption("continue"), "cmdInput")
      obj$insertWithTagsByName(end_iter$iter, i, "cmdInput")
    }
  }
  moveViewport(obj)
}

## This adds the command to the history stack and moves the pointer.
add_history <- function(obj, cmd) {
  history <- obj$GetData("history"); ptr <- obj$GetData("ptr")
  history <- c(history, cmd)
  ptr <- length(history)
  obj$SetData("ptr", ptr)
  obj$SetData("history", history)
}

## these next two functions scroll through the history
scroll_history_up <- function(obj) {
  ## move through history
  ptr <- obj$GetData("ptr") - 1
  if(ptr > 0)
    replace_cmd(obj, obj$GetData("history")[[ptr]])
  obj$SetData("ptr", max(ptr,0))
  obj$PlaceCursor(obj$GetEndIter()$iter)
}

scroll_history_down <- function(obj) {
  ## move through history
  ptr <- obj$GetData("ptr") + 1
  history <- obj$GetData("history")
  if(ptr <= length(history)) 
    replace_cmd(obj, history[[ptr]])
  obj$SetData("ptr", min(ptr,length(history)))
  obj$PlaceCursor(obj$GetEndIter()$iter)
}

## History bindings
## this uses Control-p and Control-n to move
ID <- gSignalConnect(view, "key-release-event", f=function(w, e, data) {
  if(e$GetState() != GdkModifierType['control-mask'])
    return(TRUE)

  obj <- w$GetBuffer()
  keyval <- e$GetKeyval()

  if(keyval == GDK_p) {
    scroll_history_up(obj)
  } else if(keyval == GDK_n) {
    scroll_history_down(obj)
  }
  return(TRUE)
})


###################################################
### code chunk number 309: menus.Rnw:3-4
###################################################
require(RGtk2)


###################################################
### code chunk number 310: rgtk2-menus-actions-constructor
###################################################
action <- gtkAction(name = "ok", label = "_Ok", 
             tooltip = "An OK button", stock.id = "gtk-ok")


###################################################
### code chunk number 311: rgtk2-menus-actions-activate
###################################################
gSignalConnect(action, "activate", 
               f = function(action, data) {
                 print(action$getName())
               })


###################################################
### code chunk number 312: ConnectAction
###################################################
button <- gtkButton()
button$setRelatedAction(action)


###################################################
### code chunk number 313: rgtk2-menus-action-group
###################################################
group <- gtkActionGroup()
group$addActionWithAccel(action, "<control>O")


###################################################
### code chunk number 314: rgtk2-menus-toggle-action
###################################################
full_screen_act <- 
  gtkToggleAction("fullscreen", "Full screen", 
                  "Toggle full screen",
                  stock.id = "gtk-fullscreen")
gSignalConnect(full_screen_act, "toggled", function(action) {
  if(full_screen_action['active'])
    window$fullscreen()
  else
    window$unfullscreen()
})


###################################################
### code chunk number 315: showGUI
###################################################
window <- gtkWindow(show=FALSE)
window['title'] <- "Action with button example"
window$add(button)
window$showAll()


###################################################
### code chunk number 316: rgtk2-menus-menu- bar
###################################################
menubar <- gtkMenuBar()


###################################################
### code chunk number 317: rgtk2-menus-menu
###################################################
file_menu <- gtkMenu()


###################################################
### code chunk number 318: rgtk2-menus-menuitem
###################################################
file_item <- gtkMenuItemNewWithMnemonic(label = "_File")
file_item$setSubmenu(file_menu)


###################################################
### code chunk number 319: rgtk2-menus-append
###################################################
menubar$append(file_item)


###################################################
### code chunk number 320: rgtk2-menus-open
###################################################
open_item <- gtkMenuItemNewWithMnemonic("_Open")


###################################################
### code chunk number 321: rgtk2-menus-open-activate
###################################################
gSignalConnect(open_item, "activate", function(item) {
  file.show(file.choose())
})


###################################################
### code chunk number 322: rgtk2-menus-append-item
###################################################
file_menu$append(open_item)


###################################################
### code chunk number 323: rgtk2-menus-save-action
###################################################
save_action <- 
  gtkAction("save", "Save", "Save object", "gtk-save")


###################################################
### code chunk number 324: rgtk2-menus-save-item
###################################################
save_item <- save_action$createMenuItem()
file_menu$append(save_item)


###################################################
### code chunk number 325: rgtk2-menus-separator
###################################################
file_menu$append(gtkSeparatorMenuItem())


###################################################
### code chunk number 326: rgtk2-menus-toggle-item
###################################################
auto_save_action <- gtkToggleAction("autosave", "Autosave", 
                                    "Enable autosave")
auto_save_item <- auto_save_action$createMenuItem()
file_menu$append(auto_save_item)


###################################################
### code chunk number 327: rgtk2-menus-window
###################################################
main_mindow <- gtkWindow()
vbox <- gtkVBox()
main_mindow$add(vbox)
vbox$packStart(menubar, FALSE, FALSE)


###################################################
### code chunk number 328: "menubar-ex"
###################################################
popup <- gtkMenu()                       # top level
popup$append(gtkMenuItem("cut"))
popup$append(gtkMenuItem("copy"))
popup$append(gtkSeparatorMenuItem())
popup$append(gtkMenuItem("paste"))


###################################################
### code chunk number 329: rgtk2-menus-popup-button
###################################################
button <- gtkButton("Click me with right mouse button")
window <- gtkWindow(); window$setTitle("Popup menu example")
window$add(button)


###################################################
### code chunk number 330: ex-RGtk2-menu-popup.Rnw:22-32
###################################################
gSignalConnect(button, "button-press-event",
  f = function(button, event, menu) {
    if(event$getButton() == 3 ||
       (event$getButton() == 1 && # a mac
        event$getState() == GdkModifierType['control-mask'])) 
      gtkMenuPopup(menu, 
                   button = event$getButton(),
                   activate.time = event$getTime())
    return(FALSE)
  }, data = popup)


###################################################
### code chunk number 331: ex-RGtk2-menu-popup.Rnw:43-48
###################################################
sapply(popup$getChildren(), function(child) {
  if(!inherits(child, "GtkSeparatorMenuItem")) # skip these
    gSignalConnect(child, "activate",
           f = function(child, data) message("replace me"))
})


###################################################
### code chunk number 332: rgtk2-menus-toolbar-construct
###################################################
toolbar <- gtkToolbar()


###################################################
### code chunk number 333: rgtk2-menus-toolbar-open-item
###################################################
open_button <- gtkToolButton(stock.id = "gtk-open") 


###################################################
### code chunk number 334: rgtk2-menus-toolbar-add
###################################################
toolbar$add(open_button)


###################################################
### code chunk number 335: rgtk2-menus-toolbar-save-item
###################################################
save_button <- save_action$createToolItem()
toolbar$add(save_button)


###################################################
### code chunk number 336: rgtk2-menus-toolbar-separator
###################################################
toolbar$add(gtkSeparatorToolItem())


###################################################
### code chunk number 337: rgtk2-menus-toolbar-toggle
###################################################
full_screen_button <- full_screen_act$createToolItem()
toolbar$add(full_screen_button)


###################################################
### code chunk number 338: rgtk2-menus-toolbar-style (eval = FALSE)
###################################################
## toolbar$setStyle("icon")


###################################################
### code chunk number 339: rgtk2-menus-toolbar-is-important
###################################################
full_screen_act["is-important"] <- TRUE


###################################################
### code chunk number 340: rgtk2-menus-toolbar-expand (eval = FALSE)
###################################################
## expander <- gtkSeparatorToolItem()
## expander["draw"] <- FALSE
## toolbar$add(expander)
## toolbar$childSet(expander, expand = TRUE)


###################################################
### code chunk number 341: rgtk2-menus-toolbar-help
###################################################
help_action <- gtkAction("help","Help","Get help","gtk-help")
toolbar$add(help_action$createToolItem())


###################################################
### code chunk number 342: rgtk2-menus-toolbar-place
###################################################
vbox$packStart(toolbar, FALSE, FALSE)


###################################################
### code chunk number 343: rgtk2-mennus-toolbar-color-button
###################################################
gdk_color <- gdkColorParse(palette()[1])$color
color_button <- gtkColorButton(gdk_color)


###################################################
### code chunk number 344: rgtk2-menus-toolbar-color-menu
###################################################
colorMenuItem <- function(color) {
  drawing_area <- gtkDrawingArea()
  drawing_area$setSizeRequest(20, 20)
  drawing_area$modifyBg("normal", color)
  image_item <- gtkImageMenuItem(color)
  image_item$setImage(drawing_area)
  image_item
}
color_items <- sapply(palette(), colorMenuItem)
color_menu <- gtkMenu()
for (item in color_items)
  color_menu$append(item)


###################################################
### code chunk number 345: rgtk2-menus-toolbar-color-cb
###################################################
colorMenuItemActivated <- function(item) {
  color <- gdkColorParse(item$getLabel())$color
  color_button$setColor(color)
}
sapply(color_items, gSignalConnect, "activate", 
       colorMenuItemActivated)


###################################################
### code chunk number 346: ex-RGtk2-color-tool-button.Rnw:54-55
###################################################
toolbar <- gtkToolbar()


###################################################
### code chunk number 347: rgtk2-menus-toolbar-menu
###################################################
menu_button <- gtkMenuToolButton(color_button, "Color")
menu_button$setMenu(color_menu)
toolbar$add(menu_button)


###################################################
### code chunk number 348: rgtk2-menus-tool-item-group (eval = FALSE)
###################################################
## file_group <- gtkToolItemGroup("File")
## file_group$add(gtkToolButton(stock.id = "gtk-open"))
## file_group$add(save_action$createToolItem())
## help_group <- gtkToolItemGroup("Help")
## help_group$add(help_action$createToolItem())


###################################################
### code chunk number 349: rgtk2-menus-tool-palette (eval = FALSE)
###################################################
## palette <- gtkToolPalette()
## palette$add(file_group)
## palette$add(help_group)


###################################################
### code chunk number 350: rgtk2-menus-tool-palette-collapse (eval = FALSE)
###################################################
## help_group$setCollapsed(TRUE)


###################################################
### code chunk number 351: gtk-app-status-bar
###################################################
statusbar <- gtkStatusbar()
io_id <- statusbar$getContextId("I/O")
statusbar$push(io_id, "Incomplete final line")
## ...
statusbar$pop(io_id)


###################################################
### code chunk number 352: menus.Rnw:445-447
###################################################
info_bar <- gtkInfoBar(show=FALSE)
info_bar$setNoShowAll(TRUE)


###################################################
### code chunk number 353: menus.Rnw:456-459
###################################################
label <- gtkLabel("Warning, Warning ....")
info_bar$setMessageType("warning")            
info_bar$getContentArea()$add(label)


###################################################
### code chunk number 354: menus.Rnw:463-465
###################################################
info_bar$addButton(button.text = "gtk-ok",
                   response.id = GtkResponseType['ok'])


###################################################
### code chunk number 355: menus.Rnw:473-475
###################################################
gSignalConnect(info_bar, "response", 
               function(info_bar, resp.id) info_bar$hide())


###################################################
### code chunk number 356: addToWinodw
###################################################
vbox$packStart(info_bar, expand = FALSE)
info_bar$show()


###################################################
### code chunk number 357: helperFUnction
###################################################
require(RGtk2)

##' helper function to bypass lack of cached value in method call
##'
##' @param meth method name
##' @param obj method of object's class
##' @return the method
call_meth <- function(meth, obj) {
  if(exists(meth, obj, inherits=FALSE))
    get(meth, obj)
  else
    methods:::envRefInferField(obj, meth, getClass(class(obj)), obj)
}


###################################################
### code chunk number 358: ex-RGtk2-UImanager-II.Rnw:34-36
###################################################
## Stub UI Manager instance for use with examples
uimanager <- gtkUIManager()


###################################################
### code chunk number 359: ui-xml
###################################################
ui.xml <- readLines(out <- textConnection('
<ui>
  <menubar name="menubar">
    <menu name="FileMenu" action="File">
      <menuitem action="Save"/>
      <menuitem action="SaveAs" />
      <menu name="Export" action="Export">
        <menuitem action="ExportToCSV" />
        <menuitem action="ExportToSaveFile" />
      </menu>
      <separator />
      <menuitem name="FileQuit" action="CloseWindow" />
    </menu>
    <menu action="Edit">
      <menuitem name="EditUndo" action="Undo" />
      <menuitem name="EditRedo" action="Redo" />
      <menuitem action="ChangeColumnName" />
    </menu>
    <menu action="Tools">
      <menuitem action="Filter" />
      <menuitem action="Sort" />
    </menu>
  </menubar>
  <toolbar name="toolbar">
    <toolitem action="Save"/>
    <toolitem action="SaveAs"/>
    <separator />
    <toolitem action="CloseWindow"/>
  </toolbar>
</ui>'), warn=FALSE)
close(out)


###################################################
### code chunk number 360: loadUIFromString
###################################################
id <- uimanager$addUiFromString(ui.xml)


###################################################
### code chunk number 361: ex-RGtk2-UImanager-II.Rnw:93-94
###################################################
fun <- function(...) {}


###################################################
### code chunk number 362: ex-RGtk2-UImanager-II.Rnw:102-118
###################################################
file_list <- 
  list(## name, ID, label, accelerator, tooltip, callback
       list("File",NULL,"_File",NULL,NULL,NULL),
       list("Save", "gtk-save", "Save", "<ctrl>S", 
            "Save data to variable", fun),
       list("SaveAs", "gtk-save", "Save as...", NULL, 
            "Save data to variable", fun),
       list("Export", NULL, "Export", NULL, NULL, NULL),
       list("ExportToCSV", "gtk-export", "Export to CSV", 
            NULL, "Save data to CSV file", fun),
       list("ExportToSaveFile", "gtk-export", 
            "Export to save file", NULL, 
            "Save data to save() file", fun),
       list("CloseWindow", "gtk-close", "Close window", 
            "<ctrl>W", "Close current window", fun)
       )


###################################################
### code chunk number 363: addActionGroup
###################################################
action_group <- gtkActionGroup("FileGroup")
action_group$addActions(file_list)


###################################################
### code chunk number 364: ex-RGtk2-UImanager-II.Rnw:129-130
###################################################
uimanager$insertActionGroup(action_group, 0)


###################################################
### code chunk number 365: GUILayout (eval = FALSE)
###################################################
## window <- gtkWindow(show = FALSE)
## ##
## vbox <- gtkVBox()
## window$add(vbox)
## ##
## menubar <- uimanager$getWidget("/menubar")
## vbox$packStart(menubar, FALSE)
## toolbar <- uimanager$getWidget("/toolbar")
## vbox$packStart(toolbar, FALSE)
## ## ...


###################################################
### code chunk number 366: ex-RGtk2-UImanager-II.Rnw:171-172 (eval = FALSE)
###################################################
## window$addAccelGroup(uimanager$getAccelGroup())


###################################################
### code chunk number 367: ex-RGtk2-UImanager-II.Rnw:193-200
###################################################
Command <- setRefClass("Command",
                       fields = list(
                         receiver="ANY",
                         meth="character",
                         params="list",
                         old_params="list"
                         ))


###################################################
### code chunk number 368: ex-RGtk2-UImanager-II.Rnw:208-218
###################################################
Command$methods(
        initialize = function(receiver, meth, ...) {
          .params <- list(...)
          initFields(receiver = receiver, meth = meth, 
                     params = .params, old_params = .params)
          callSuper()
        },
        execute = function(params) {
          do.call(call_meth(meth, receiver), params)
        })


###################################################
### code chunk number 369: ex-RGtk2-UImanager-II.Rnw:225-233
###################################################
Command$methods(
                do = function() {
                  out <- execute(params)
                  old_params$value <<- out
                },
                undo = function() execute(old_params)
                )



###################################################
### code chunk number 370: illustrateCommand
###################################################
x <- 1
set_x <- function(value) {
  old <- x
  x <<- value
  old
}
cmd <- Command$new(.GlobalEnv, "set_x", value = 2)
cmd$do(); x


###################################################
### code chunk number 371: ex-RGtk2-UImanager-II.Rnw:252-253
###################################################
cmd$undo();


###################################################
### code chunk number 372: ex-RGtk2-UImanager-II.Rnw:255-256
###################################################
x


###################################################
### code chunk number 373: typicalAction (eval = FALSE)
###################################################
## cmd <- Command$new(df_model, "set_col_name", j=j, value=value)
## command_stack$add(cmd)


###################################################
### code chunk number 374: col_name_methods (eval = FALSE)
###################################################
## DfModel$methods(
##                 get_col_name = function(j) varnames[j,1],
##                 get_col_names = function() varnames[ ,1],
##                 set_col_name = function(j, value) {
##                   "Set name, return old"
##                   old_col_name <- get_col_name(j)
##                   varnames[j,1] <<- value
##                   old_col_name
##                 })


###################################################
### code chunk number 375: ensure_type
###################################################
##' S3 generic to ensure we don't change data type when assigning into column
##'
##' @param x column values
##' @param value new value
##' @return coerced new value
ensure_type <- function(x, value) UseMethod("ensure_type")
ensure_type.default <- function(x, value) value
ensure_type.character <- function(x, value) as.character(value)
ensure_type.factor <- function(x, value) {x[length(x) + 1] <- value; tail(x, n=1)}
ensure_type.numeric <- function(x, value) as.numeric(value)
ensure_type.integer <- function(x, value) as.integer(value)
ensure_type.logical <- function(x, value) as.logical(value)


###################################################
### code chunk number 376: DfModel
###################################################
## Define a model to hold the model for an editable data frame
sapply(c("RGtkDataFrame"), setOldClass)
DfModel <- setRefClass("DfModel",
                       fields=list(
                         store="RGtkDataFrame",
                         filtered="ANY",
                         name="character",
                         varnames="RGtkDataFrame",
                         casenames="RGtkDataFrame"
                         ))

## Initialize along with a column for filtering
DfModel$methods(
                initialize=function(DF, nm, ...) {
                  store <<- rGtkDataFrame(cbind(DF, `_visible`=rep(TRUE, nrow(DF))))
                  varnames <<- rGtkDataFrame(data.frame(names(DF), stringsAsFactors=FALSE))
                  casenames <<- rGtkDataFrame(data.frame(rownames(DF), stringsAsFactors=FALSE))
                  if(missing(nm))
                    name <<- deparse(substitute(DF))
                  else
                    name <<- nm
                  filtered <<- store$filter()
                  filtered$setVisibleColumn(length(DF))
                  callSuper()
                })

## Methods to work with the underlying data frame (Get, save, ...)
DfModel$methods(
                get_dataframe=function() {
                  DF <- store[,seq_len(ncol(store)-1)]
                  dimnames(DF) <- list(casenames[,1], varnames[,1])
                  DF
                },
                save=function(nm) {
                  "Save to global workspace"
                  if(!missing(nm))
                    name <<- nm
                  assign(name, get_dataframe(), envir=.GlobalEnv)
                  
                },
                export_to_csv=function(f)  {
                  "Export to csv file"
                  write.csv(get_dataframe(), file=f)
                },
                export_to_save=function(f) {
                  "Export using save()"
                  assign(name, get_dataframe())
                  save(list=name, file=f)
                },
                no_rows=function() dim(store)[1],
                no_cols=function() dim(store)[2] - 1L
                )
## Methods to get and set a cell value. 
DfModel$methods(
                get_cell=function(i,j) {
                  "Return cell value"
                  store[i,j]
                },
                set_cell=function(i, j, value) {
                  "Set cell, return old_value"
                  old <- get_cell(i,j)
                  store[i,j] <<- ensure_type(store[1,j], value)
                  old
                })
## Methods for column names. Similar one for rownames could be implemented, but we
## don't show these in our view. So leave to the reader/
DfModel$methods(
                get_col_name=function(j) varnames[j,1],
                get_col_names=function() varnames[,1],
                set_col_name=function(j, value) {
                  "Set name, return old"
                  old <- get_col_name(j)
                  varnames[j,1] <<- value
                  old
                })
## Code for filtering the display.
DfModel$methods(
                get_filter=function() {
                  "Return logical indicating filter"
                  store[,ncol(store)]
                },
                set_filter=function(value) {
                  "Filter by value. Return old filter value"
                  if(!is.logical(value)) stop("Filter requires a logical variable")
                  ind <- rep(value, length.out=no_rows())
                  old <- get_filter()
                  store[,ncol(store)] <<- value
                  old
                })

## In RGtk2, one can't both sort and filter by proxy. Since R makes sorting easy, 
## we let Gtk handle the filtering and implement sorting below. The "old" value 
## returned by this is what is needed to reverse a sort.
DfModel$methods(
                reorder=function(value) {
                  "Reorder data frame. Return order(value)"
                  perm <- as.integer(value)
                  if(length(perm) != nrow(store)) stop("reorder requires a permutation")
                  if(length(perm) != length(unique(perm))) stop("value has repeated values")
                  if(min(perm) != 1 || max(perm) != nrow(store)) stop("value is not permutation of row indices")

                  store[,] <<- store[perm,]
                  order(perm)  # will revers a[ind][order(ind)] is involution
                })


###################################################
### code chunk number 377: CommandStack
###################################################
## Command Stack
## A list with ptr. delegates call of do or undo to appropriate command
CommandStack <- setRefClass("CommandStack",
                            fields=list(
                              l="list",
                              ptr="integer"
                              ))
## initialize method just sets the list and pointer to a default.
CommandStack$methods(
                     initialize=function() {
                       initFields(l=list(), ptr=0L)
                       callSuper()
                     })
## do method finds the right command then delegates to the commands do method
## undo is similar
## The can_do and can_undo commands are used to check if the command stack allows for
## these operations
CommandStack$methods(
                     do=function() {
                       if(!can_do()) return()
                       cmd <- l[[ptr]]
                       ptr <<- ptr + 1L
                       cmd$do()
                     },
                     undo=function() {
                       if(!can_undo()) return()
                       cmd <- l[[ptr-1]]
                       ptr <<- ptr - 1L
                       cmd$undo()
                     },
                     can_do=function() ptr > 0 && ptr <= length(l),
                     can_undo=function() ptr > 1
                     )
## Methods to add to and clear the command stack
CommandStack$methods(
                     add=function(cmd, call=TRUE) {
                       if(ptr <= 1) {
                         l <<- list(cmd)
                         ptr <<- 1L
                       } else {
                         l <<- l[1:(ptr-1)]
                         l[[length(l) + 1]] <<- cmd
                       }
                       if(call)
                         do()
                     },
                     clear=function(cmd) {
                       l <<- list(); ptr <<- 0L
                     })


###################################################
### code chunk number 378: addCellRenderer
###################################################
## We create our cellrenderers using an S3 generic to dispatch based on the class of the column. This
## works out well, as the view is column based as well. The editable commands have 
## to find a row, a column and a value before make a command to add to the command stack.
## The row comes from the path, but must be "unfiltered" to point to the original data store. 
## The column is passed into the function by the caller.


##' Create an appropriate cell renderer
##'
##' @param x vector to display in column
##' @param nm name of vector for title
##' @param obj a DfModel instance
##' @param view GtkTreeView instance we add the cellrenderer to
##' @param command_stack a CommandStack instance needed for the callback
##' @return NULL
add_cellrenderer_by_class <- function(x, nm, obj, view, j, command_stack) UseMethod("add_cellrenderer_by_class")
add_cellrenderer_by_class.default <- function(x, nm, obj, view, j, command_stack) {
  cr <- gtkCellRendererText()
  cr['editable'] <- TRUE
  gSignalConnect(cr, "edited", f=function(cr, path, newtext) {
    i <- as.numeric(path) + 1
    i <- which(obj$get_filter())[i]     # in regular    
    value <- newtext
    cmd <- Command$new(obj, "set_cell", i=i, j=j, value=value)
    command_stack$add(cmd)
  })
  view$insertColumnWithAttributes(position=-1, 
                                  title=nm,
                                  cell=cr,
                                  text=j-1)
}

add_cellrenderer_by_class.logical <- function(x, nm, obj, view, j, command_stack) {
  cr <- gtkCellRendererToggle()
  cr['activatable'] <- TRUE
  gSignalConnect(cr, "toggled", function(w, path) {
    i <- as.numeric(path) + 1           # in filtered
    i <- which(obj$get_filter())[i]     # in regular
    value <- !obj$get_cell(i,j)
    cmd <- Command$new(obj, "set_cell", i=i, j=j, value=value)
    command_stack$add(cmd)
  })
  view$insertColumnWithAttributes(position=-1, 
                                  title=nm,
                                  cell=cr,
                                  active=j-1)
}

add_cellrenderer_by_class.factor <- function(x, nm, obj, view, j, command_stack) {
  cr <- gtkCellRendererCombo()
  cr_store <- rGtkDataFrame(sort(levels(x)))
  cr['model'] <- cr_store
  cr['text-column'] <- 0
  cr['has-entry'] <- FALSE
  cr['editable'] <- TRUE
  gSignalConnect(cr, "changed", function(w, path, iter, user.data) {
    i <- as.numeric(path) + 1
    i <- which(obj$get_filter())[i]     # in regular    
    value <- cr_store$getValue(iter, 0)$value
    cmd <- Command$new(obj, "set_cell", i=i, j=j, value=value)
    command_stack$add(cmd)
  })
  view$insertColumnWithAttributes(position=-1, 
                                  title=nm,
                                  cell=cr,
                                  text=j-1L)
}


###################################################
### code chunk number 379: EditDataFrame
###################################################
## Main reference class to edit a data frame within a GUI
## The view relies on a DataFrameModel and CommandStack instance, each of which is 
## defined within the initialize method.
EditDataFrame <- setRefClass("EditDataFrame",
                             fields=list(
                               df_model="ANY",
                               command_stack="ANY",
                               actions="list",
                               ## layout
                               mainwindow="ANY",
                               statusbar="ANY",
                               uimanager="ANY", 
                               view="ANY"
                               ))
## The initialize method makes several different calls. Here we initialize the actions into action group.
EditDataFrame$methods(
                      initialize_actions=function(box) {
                        ## our callback. Calls an appropriately named method of this class.
                        fun=function(action) {
                          meth <- action$getName()
                          out <- try(do.call(call_meth(meth, .self), list()), silent=TRUE)
                        }

                        ## Define action groups in a list
                        fileL <- list(## name, ID, label, accelerator, tooltip, callback
                                      list("File",NULL,"_File",NULL,NULL,NULL),
                                      list("Save", "gtk-save", "Save", "<ctrl>S", "Save data to variable", fun),
                                      list("SaveAs", "gtk-save", "Save as...", NULL, "Save data to variable", fun),
                                      list("Export", NULL, "Export", NULL, NULL, NULL),
                                      list("ExportToCSV", "gtk-export", "Export to CSV", NULL, "Save data to CSV file", fun),
                                      list("ExportToSaveFile", "gtk-export", "Export to save() file", NULL, "Save data to save() file", fun),
                                      list("CloseWindow", "gtk-close", "Close window", "<ctrl>W", "Close current window", fun)
                                      )

                        editL <- list(## name, ID, label, accelerator, tooltip, callback
                                      list("Edit", NULL, "_Edit", NULL, NULL, NULL),
                                      list("Undo", "gtk-undo", "Undo", "<ctrl>Z",  "Undo last command", fun),
                                      list("Redo", "gtk-redo", "Redo", "<ctrl>U", "Redo undo command", fun),
                                      list("ChangeColumnName", "gtk-change", "Change column name",
                                                          NULL, "Change a column name", fun)
                                      )
                        
                        toolL <- list(
                                      list("Tools", NULL, "_Tools", NULL, NULL, NULL),
                                      list("Filter", "gtk-filter", "Filter", NULL, "Filter data frame", fun),
                                      list("Sort", "gtk-sort", "Sort", NULL, "Sort data frame by column name", fun)                                      )
                        l <- list(fileL, editL, toolL)

                        ## create UI manager, insert action groups
                        uimanager <<- gtkUIManager()
                        for(i in seq_along(l)) {
                          ag <- gtkActionGroup(sprintf("Group%s",i))
                          ag$addActions(l[[i]])
                          uimanager$insertActionGroup(ag, i-1)
                        }
                      })

## Here we initialize the UI 
EditDataFrame$methods(
                      initialize_ui=function() {
                        ## define xml specifying menu bars and toolbars
                        ui.xml <- readLines(out <- textConnection('
<ui>
  <menubar name="menubar">
    <menu name="FileMenu" action="File">
      <menuitem action="Save"/>
      <menuitem action="SaveAs" />
      <menu name="Export" action="Export">
        <menuitem action="ExportToCSV" />
        <menuitem action="ExportToSaveFile" />
      </menu>
      <separator />
      <menuitem name="FileQuit" action="CloseWindow" />
    </menu>
    <menu action="Edit">
      <menuitem name="EditUndo" action="Undo" />
      <menuitem name="EditRedo" action="Redo" />
      <menuitem action="ChangeColumnName" />
    </menu>
    <menu action="Tools">
      <menuitem action="Filter" />
      <menuitem action="Sort" />
    </menu>
  </menubar>
  <toolbar name="toolbar">
    <toolitem action="Save"/>
    <toolitem action="SaveAs"/>
    <separator />
    <toolitem action="CloseWindow"/>
  </toolbar>
</ui>'), warn=FALSE)
                        close(out)
                        ## specify the UI using XML specification
                        id <- uimanager$addUiFromString(paste(ui.xml, collapse="\n")) 
                      })

## Here we layout the GUI
EditDataFrame$methods(
                      make_gui=function() {
                        DF <- df_model$get_dataframe()
                        nms <- names(DF)
                        view <<- gtkTreeView(df_model$filtered)
                        sapply(seq_len(length(DF)), function(j) {
                          add_cellrenderer_by_class(DF[[j]], nms[j], df_model, view, j, command_stack)
                        })
                        ##
                        ## place into GUI
                        mainwindow <<- w <- gtkWindow(show=FALSE)
                        #
                        vbox <- gtkVBox()
                        w$add(vbox)
                        #
                        menubar <- uimanager$getWidget("/menubar")
                        vbox$packStart(menubar, FALSE)
                        toolbar <- uimanager$getWidget("/toolbar")
                        vbox$packStart(toolbar, FALSE)
                        w$addAccelGroup(uimanager$getAccelGroup())
                        ##
                        sw <- gtkScrolledWindow()
                        sw$add(view)
                        vbox$PackStart(sw, TRUE, TRUE)
                        ##
                        statusbar <<- gtkStatusbar()
                        statusbar$getChildren()[[1]]$setSizeRequest(-1, 25)
                        vbox$PackStart(statusbar, FALSE)
                        w$show()
                      })

## This method call updates the GUI: sets the redo/undo buttons and the status bar.
EditDataFrame$methods(
                      update_UI=function(event="") {
                        ## update actions
                        ## Could save undo/redo actions as we look them up inefficiently below
                        undo <- redo <- NULL
                        for(i in uimanager$getActionGroups()) {
                          tmp <- i$getAction("Redo")
                          if(!is.null(tmp)) redo <- tmp
                          tmp <- i$getAction("Undo")
                          if(!is.null(tmp)) undo <- tmp
                        }
                        undo$setSensitive(command_stack$can_undo())
                        redo$setSensitive(command_stack$can_do())

                        ## update status bar
                        tpl <- "Editing %s. Showing %s lines of %s."
                        statusbar$push(statusbar$getContextId("message"),
                                       sprintf(tpl, df_model$name,
                                               sum(df_model$get_filter()),
                                               df_model$no_rows()))
                        
                      })
## This sets up a callback when parts of the DataFrameModel change. Here is how
## we synchronize column names and why we used a RGtkDataFrame class to hold them in the definition
## of the DataFrameModel class
EditDataFrame$methods(
                      synchronize_view=function() {
                        gSignalConnect(df_model$store, "row-changed", function(model, path, iter) {
                          update_UI()
                        })
                        gSignalConnect(df_model$varnames, "row-changed", function(model, path, iter) {
                          j <- as.numeric(path$toString()) + 1
                          value <- df_model$varnames[j,1]
                          col <- view$getColumn(j-1)
                          col['title'] <- value
                          update_UI()
                        })
                      })
## Finally an initialization method
EditDataFrame$methods(
                      initialize=function(DF) {
                        if(!is.data.frame(DF))
                          stop("Requires a data frame")
                        initFields(df_model=DfModel$new(DF),
                                   command_stack=CommandStack$new())

                        initialize_actions()
                        initialize_ui()
                        make_gui()
                        synchronize_view()
                        update_UI()
                        callSuper()
                      })


###################################################
### code chunk number 380: Actions
###################################################
## Actions are defined here
## Basically we delegate down to data frame model
## We are lazy about some dialogs, so use the gWidgets package
require(gWidgets); options(guiToolkits="RGtk2")
EditDataFrame$methods(
                      Save=function() {
                        df_model$save()
                        command_stack$clear()
                      })

EditDataFrame$methods(
                      Undo=function() {command_stack$undo()},
                      Redo=function() {command_stack$do()}
                      )

EditDataFrame$methods(
                      SaveAs=function() {
                        current_vars <- ls(envir=.GlobalEnv)
                        dlg <- gbasicdialog("Select a variable name...", parent=mainwindow, handler=function(h,...) {
                          var <- svalue(e)
                          if(nchar(var)) {
                            if(exists(var, .GlobalEnv)) {
                              if(!gconfirm(c("Variable exists", "Really overwrite?"), parent=dlg))
                                return()
                            }
                            df_model$save(var)
                            update_UI()
                            command_stack$clear()
                          }
                        })
                        g <- ggroup(cont=dlg, horizontal=FALSE)
                        glabel("Variable name to save as:", cont=g)
                        e <- gedit("", cont=g)
                        e[] <- current_vars
                        visible(dlg, set=TRUE)
                      })

EditDataFrame$methods(
                      ExportToCSV=function() {
                        f <- gfile("Select a filename", type="save")
                        if(!is.na(f))
                          df_model$export_to_csv(f)
                      })

EditDataFrame$methods(
                      ExportToSaveFile=function() {
                        f <- gfile("Select a filename", type="save")
                        if(!is.na(f))
                          df_model$export_to_save(f)
                      })

EditDataFrame$methods(
                      CloseWindow=function() {
                        if(command_stack$can_undo() || command_stack$can_do()) {
                          if(!gconfirm(c("Really quit", "There are pending changes"),
                                       parent=mainwindow))
                            return()
                        }
                        mainwindow$destroy()
                      })

EditDataFrame$methods(
                      ChangeColumnName=function() {
                        j <- NA; value <- character(0)
                        ## get column and new names
                        dlg <- gbasicdialog("Rename a column", parent=mainwindow)
                        g <- ggroup(horizontal=FALSE, cont=dlg)
                        varnames <- df_model$get_col_names()
                        tbl <- gtable(data.frame(Variables=varnames), cont=g, expand=TRUE)
                        size(tbl) <- c(300, 250)
                        l <- glabel("Select a variable", cont=g)
                        e <- gedit("", cont=g); enabled(e) <- FALSE
                        addHandlerClicked(tbl, handler=function(h,...) {
                          val <- svalue(h$obj)
                          assign("j", match(val, varnames), inherits=TRUE)
                          if(!is.na(j)) {
                            svalue(l) <- sprintf("Change %s to:", val)
                            enabled(e) <- TRUE
                          } else {
                            svalue(l) <- "Select a variable"
                            svalue(e) <- ""
                            enabled(e) <- FALSE
                          }
                        })
                        addHandlerKeystroke(e, handler=function(h,...) {
                          assign("value", svalue(h$obj), inherits=TRUE)
                        })
                        ret <- visible(dlg, set=TRUE)
                        if(ret && !is.na(j)) {
                          cmd <- Command$new(df_model, "set_col_name", j=j, value=value)
                          command_stack$add(cmd)
                        }
                      })
                      
EditDataFrame$methods(
                      Filter=function() {
                        ind <- NULL
                        dlg <- gbasicdialog("Enter an expression", parent=mainwindow, handler=function(h,...) {
                          val <- svalue(e)
                          DF <- df_model$get_dataframe()
                          out <- try(eval(parse(text=val), DF), silent=FALSE)
                          if(!inherits(out, "try-error"))
                            assign("ind", out, inherits=TRUE)
                        })
                        g <- ggroup(cont=dlg, horizontal=FALSE)
                        glabel("Enter an expression to filter by:", cont=g)
                        e <- gedit("", cont=g)
                        ret <- visible(dlg, set=TRUE)
                        if(ret && is.logical(ind)) {
                          cmd <- Command$new(df_model, "set_filter", value=ind)
                          command_stack$add(cmd)
                        }
                      })

EditDataFrame$methods(
                      Sort=function() {
                        perm <- integer(0)
                        DF <- df_model$get_dataframe()
                        varnames <- df_model$get_col_names()
                        dlg <- gbasicdialog("Sort by:", parent=mainwindow, handler=function(h,...) {
                          var <- svalue(tbl, index=TRUE)
                          if(length(var) == 0) return()
                          x <- DF[,var]
                          assign("perm", order(x, decreasing=svalue(decreasing)), inherits=TRUE)
                        })
                        g <- ggroup(horizontal=FALSE, cont=dlg)
                        tbl <- gtable(data.frame(Variables=varnames), cont=g)
                        size(tbl) <- c(300, 250)
                        decreasing <- gcheckbox("Decreasing?", checked=FALSE, cont=g)
                        ret <- visible(dlg, set=TRUE)
                        ##
                        if(ret && length(perm)) {
                          cmd <- Command$new(df_model, "reorder",  value=perm)
                          command_stack$add(cmd)
                        }
                      })


###################################################
### code chunk number 381: testItOut
###################################################
## Test it out....
require(MASS)
DF <- Cars93[sample(1:93, 20),c(1, 5, 26)]; DF$American <- DF$Origin == "USA"
a = EditDataFrame$new(DF)


###################################################
### code chunk number 382: oldWay
###################################################
## Old way to add actions, menu bar, For comparison
## not called by the initialize method
EditDataFrame$methods(
                      initialize_actions_old=function() {
                           ## actions. Must have a matching method
                        al <- list()
                        al$save <- gtkAction("Save", "Save", "Save data to variable", "gtk-save")
                        al$saveas <- gtkAction("SaveAs", "Save as...", "Save data to variable", "gtk-save")
                        al$exportAsCSV <- gtkAction("ExportToCSV", "Export to CSV", "Save data to CSV file", "gtk-export")
                        al$exportAsSaveFile <- gtkAction("ExportToSaveFile", "Export to save() file", "Save data to save() file", "gtk-export")
                        al$close <- gtkAction("CloseWindow", "Close window", "Close current window", "gtk-close")
                        ## Edit menu
                        al$undo <- gtkAction("Undo", "Undo", "Undo last command", "gtk-undo")
                        al$redo <- gtkAction("Redo", "Redo", "Redo undo command", "gtk-redo")
                        al$change_column_name <- gtkAction("ChangeColumnName", "Change column name",
                                                          "Change a column name", "gtk-change")
                        ## Tools
                        al$filter <- gtkAction("Filter", "Filter", "Filter data frame", "gtk-filter")
                        al$sort <- gtkAction("Sort", "Sort", "Sort data frame by column name", "gtk-sort")
                        
                        ## stub handler
                        sapply(al, gSignalConnect, "activate", function(action) {
                          meth <- action$getName()
                          out <- try(do.call(get(meth, .self), list()), silent=TRUE)
                          print(out)
                        })
                        
                        actions <<- al
                      },
                      make_menu=function(box) {
                        mb <- gtkMenuBar()

                        fileMenu <- gtkMenu()
                        fileItem <- gtkMenuItem("File")
                        fileItem$setSubmenu(fileMenu)
                        sapply(c("save","saveas", "exportAsCSV","exportAsSaveFile","close"),
                               function(act)
                               fileMenu$append(actions[[act]]$createMenuItem()))

                        editMenu <- gtkMenu()
                        editItem <- gtkMenuItem("Edit")
                        editItem$setSubmenu(editMenu)
                        sapply(c("undo","redo", "change_column_name"),
                               function(act)
                               editMenu$append(actions[[act]]$createMenuItem()))

                        toolsMenu <- gtkMenu()
                        toolsItem <- gtkMenuItem("Tools")
                        toolsItem$setSubmenu(toolsMenu)
                        sapply(c("filter", "sort"),
                               function(act)
                               toolsMenu$append(actions[[act]]$createMenuItem()))

                        sapply(list(fileItem, editItem, toolsItem), mb$append)
                        box$packStart(mb, FALSE)

                      }
)


###################################################
### code chunk number 383: gtk-class-def
###################################################
tform_scale_type <- 
  gClass("RTransformedHScale", "GtkHScale",
         .props = list(
           gParamSpec(type = "R", name = "expr", nick = "e", 
                      blurb = "Transformation of scale value",                 
                      default.value = expression(x))
           ),
         GtkScale = list(
           format_value = function(self, x) 
             as.character(self$transformValue(x))
           ),
         .public = list(
           getExpr = function(self) self["expr"],
           getTransformedValue = function(self) 
             self$transformValue(self$value)
           ),
         .private = list(
           transformValue = function(self, x) 
             eval(self$expr, list(x = x))
          )
         )


###################################################
### code chunk number 384: gtk-class-madata
###################################################
n <- 5000
backbone <- rnorm(n)
ma_data <- cbind(backbone + c(rnorm(3 * (n / 4), sd = 0.1), 
                              rt(n/4, 80)), 
                 backbone + c(rnorm(3 * (n / 4), , 0.1), 
                              rt(n / 4, 80)))
ma_data <- apply(ma_data, 2, function(col) col - min(col))


###################################################
### code chunk number 385: gtk-class-instance
###################################################
adj <- gtkAdjustment(0.5, 0.15, 1.00, 0.05, 0.5, 0)
s <- gObject(tform_scale_type, adjustment = adj, 
             expr = expression(x^3))
gSignalConnect(s, "value_changed", function(scale) {
  plot(ma_data, col = rgb(0,0,0, scale$getTransformedValue()),
       xlab = "Replicate 1", ylab = "Replicate 2", 
       main = "Expression levels of WT at time 0", pch = 19)
})


###################################################
### code chunk number 386: gtk-class-window (eval = FALSE)
###################################################
## win <- gtkWindow(show = FALSE)
## da <- gtkDrawingArea()
## vbox <- gtkVBox()
## vbox$packStart(da)
## vbox$packStart(s, FALSE)
## win$add(vbox)
## win$setDefaultSize(400, 400)
## #
## require(cairoDevice)
## asCairoDevice(da)
## #
## win$showAll()
## par(pty = "s")
## s$setValue(0.7)


###################################################
### code chunk number 387: ch-RGtk2.Rnw:162-165
###################################################
options(prompt="> ")
options(continue="+ ")
options(width=80)


