### R code from vignette source 'ex-RGtk2-button-new-stock-hardway.Rnw'

###################################################
### code chunk number 1: ButtonNewFromStockHardWay
###################################################
b <- gtkButton()
g <- gtkHBox()
pbuf <- b$renderIcon("gtk-ok", size=GtkIconSize["button"]) 
i <- gtkImageNewFromPixbuf(pbuf)
i['xalign'] <- 1; i['xpad'] <- 5        # right align with padding
g$packStart(i, expand=FALSE)
l <- gtkLabel(gettext("ok")); 
l['xalign'] <- 0 # left align
g$packStart(l, expand=TRUE, fill=TRUE)
b$add(g)
## show it
w <- gtkWindow(); w$add(b)


