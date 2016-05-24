
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
