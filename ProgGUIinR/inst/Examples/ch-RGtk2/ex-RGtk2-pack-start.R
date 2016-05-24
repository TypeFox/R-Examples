### R code from vignette source 'ex-RGtk2-pack-start.Rnw'

###################################################
### code chunk number 1: PackStartExamples
###################################################
## example of differences between spacing -- argument to gtk*Box and padding -- argument to PackStart
## example of expand=FALSE; expand=TRUE, fill=FALSE; and expand=TRUE, fill=TRUE
## Modified from tutorial http://www.gtk.org/tutorial1.2/gtk_tut-4.html

library(RGtk2)
w <- gtkWindow()
w$setTitle("Example of packStart options")
w$setSizeRequest(700, 300)
w$modifyBg(state="normal", color="gray30")
g <- gtkVBox(homogeneous=FALSE); w$add(g)

g1 <- gtkHBox()
l <- list(gtkButton("button"),
          gtkButton("expand=FALSE"),
          gtkButton("padding=0")
          )
sapply(l, function(i) g1$PackStart(i, expand=FALSE, padding=0))

g1 <- gtkHBox(); g$packStart(g1)
l <- list(gtkButton("button"),
          gtkButton("expand=FALSE"),
          gtkButton("fill=FALSE"),
          gtkButton("padding=0")
          )
sapply(l, function(i) g1$PackStart(i, expand=FALSE, padding=0))

g1 <- gtkHBox(); g$packStart(g1)
l <- list(gtkButton("button"),
          gtkButton("expand=FALSE"),
          gtkButton("fill=FALSE"),
          gtkButton("padding=10")
          )
sapply(l, function(i) g1$PackStart(i, expand=FALSE, padding=10))

g1 <- gtkHBox(spacing=10); g$packStart(g1)
l <- list(gtkButton("button"),
          gtkButton("expand=FALSE"),
          gtkButton("fill=FALSE"),
          gtkButton("spacing=10")
          )
sapply(l, function(i) g1$PackStart(i, expand=FALSE, padding=0))

g1 <- gtkHBox();  g$packStart(g1)
l <- list(gtkButton("button"),
          gtkButton("expand=TRUE"),
          gtkButton("fill=FALSE"),
          gtkButton("padding=0")
          )
sapply(l, function(i) g1$PackStart(i, expand=TRUE, fill=FALSE, padding=0))

g1 <- gtkHBox();  g$packStart(g1)
l <- list(gtkButton("button"),
          gtkButton("expand=TRUE"),
          gtkButton("fill=TRUE"),
          gtkButton("padding=0")
          )
sapply(l, function(i) g1$PackStart(i, expand=TRUE, fill=TRUE, padding=0))


w$showAll()


