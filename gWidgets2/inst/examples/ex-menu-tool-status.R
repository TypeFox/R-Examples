require(gWidgets2)

##' A main window (gwindow) can be "decorated" with
##' * a toolbar. This shows in the window at the top or below the menu. In tcltk, the style
##'   option is not heeded
##' * a menubar. This shows at the top of the window in most cases, but may appear in the
##'   system menubar (Qt and Mac, Aqua tcltk, ..)
##' * a statusbar. This shows at the bottom of the page, usually set off by appearing in a gutter,
##'   though this is style dependent
##' * an infobar. In RGtk2 this is an actual widget that shows transient messages below the toolbar.
##'   The other toolkits may have something similar hacked in. The infobar -- if present -- is called
##'   by gaction when the parent is a main window.

w <- gwindow("Example of various 'bars'", visible=FALSE)

## some list of actions
acts <- list(open=gaction("open", icon="open", handler=function(...) svalue(sb) <- "open", parent=w),
             save=gaction("save", icon="save", handler=function(...) svalue(sb) <- "save", parent=w),
             undo=gaction("undo", icon="undo", handler=function(...) svalue(sb) <- "undo", parent=w),
             redo=gaction("redo", icon="redo", handler=function(...) svalue(sb) <- "redo", parent=w),
             quit=gaction("quit", icon="quit", handler=function(...) dispose(w), parent=w)
             )

## menu bar list. Use nested lists to get heirarchy (not shown)
mb_list <- list(File=list(
                  acts[[1]],
                  acts[[2]],
                  gseparator(parent=w),
                  acts[[5]]
                  ),
                Edit=list(
                  acts[[3]],
                  acts[[4]]
                  )
                )

## toolbar lists are "flat"
tb_list <- list(acts[[1]], acts[[2]], acts[[5]])

mb <- gmenu(mb_list, cont=w)
tb <- gtoolbar(tb_list, cont=w)
sb <- gstatusbar("some status text", cont=w)

g <- gvbox(cont=w)
glabel("
Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Suspendisse eu magna mi. Praesent nisl dolor, consectetur
fermentum eleifend nec, dignissim non sapien. Fusce egestas
fringilla congue. Praesent consequat pretium velit, a
eleifend erat fermentum a. Cras id risus in tellus sagittis
imperdiet.
", cont=g)

button_gp <- ggroup(cont=g)
addSpring(button_gp)
gbutton("info bar", cont=button_gp, handler=function(h,...) {
  galert("Lorem Ipsum is simply dummy text of the printing and typesetting industry.", parent=w)
})


visible(w) <- TRUE
