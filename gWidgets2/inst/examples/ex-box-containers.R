## Expand, fill, anchor,  arguments
library(gWidgets2)


##' The box containers (`ggroup`, its alias `gvbox`, `gframe`,
##' `gexpandgroup`) have children added through their `add` method, which
##' really is the `add_child` reference class method.  This is usually
##' called by a widget constructor, as is the case when a parent container
##' is passed to the `container` argument of the constructor.  The `add`
##' method has arguments `expand`, `fill`, and `anchor` to control how the
##' child widgets are placed when added.
##' 
##' * When a parent has more space than is requested by its children, the
##'   space is shared amongst children that request it by having
##'   expand=TRUE. For some of the toolkits, the expansion is
##'   weighted. One can pass a logical or a non-negative number to this
##'   argument. If non-negative *and* the toolkit uses weights, then this
##'   will be used.
##' 
##' * When a child widget\'s space is expanded, the widget can fill all
##'   the space, or fill in one direction (or neither). The `fill`
##'   argument takes care of this specification with values `both`, `x`,
##'   `y` or NULL. This does not always work. Some toolkits will always
##'   set the fill value so that expansion happens orthogonally to the
##'   direction of packing (for a vertical box, the x direction).
##' 
##' * When no fill but there is space, then the `anchor` argument can be
##'   used to anchor the child within the allocated region. This uses points
##'   in {-1, 0, 1} cross {-1, 0, 1} to specify one of the 9 compass points.




w <- gwindow("Box containers", visible=FALSE)


nb <- gnotebook(cont=w)

g <- ggroup(cont=nb, horizontal=FALSE, label="fill")

gbutton("e=F, f=NULL", expand=FALSE, fill=NULL, cont=g)
gbutton("e=T, f=NULL", expand=TRUE, fill=NULL, cont=g)
gbutton("e=F, f=x", expand=FALSE, fill="x", cont=g)
gbutton("e=F, f=y", expand=FALSE, fill="y", cont=g)
gbutton("e=F, f=both", expand=FALSE, fill="both", cont=g)
gbutton("e=T, f=x", expand=TRUE, fill="x", cont=g)
gbutton("e=T, f=y", expand=TRUE, fill="y", cont=g)
gbutton("e=T, f=both", expand=TRUE, fill="both", cont=g)


## * all toolkits exand the button horizontally unless fill="y" is given. In RGtk2, this does not override
## * expand for tcltk and RGtk2 will give same allocatin of vertical space to the bottom three buttons. Not
##   so with Qt -- that depends also on the fill value


f <- gframe("frame", cont=nb, horizontal=FALSE, label="anchor")
glabel("-1,-1", expand=TRUE, anchor=c(-1,-1), cont=f)
glabel("1,1", expand=TRUE, anchor=c(1,1), cont=f)
glabel("0,0", expand=TRUE, anchor=c(0,0), cont=f)

## * In Qt, works as expected. As expand is given, all share same vertical space. The anchor then works within that
## * In RGtk2 the x component works, but not the y, which is centered in each case
## * In tcltk the y component works, but not the x. Each is left aligned

vb <- gvbox(cont=nb, label="spacing")
svalue(vb) <- 20                        # between child space
vb$set_borderwidth(10)                  # margin border. No S3 method.
gbutton("button 1", cont=vb)
gbutton("button 2", cont=vb)
addSpring(vb)

## * tcltk: spacing around child object. svalue<- doesn't work after children added. 
## * RGtk2 spacing in vertical pixels between children. 
## * Qt like tcltk, vertical spacing between children.

exp <- gexpandgroup("expand group", cont=nb, horizontal=FALSE, label="exp")
glabel("Lorem ipsum ...", cont=exp)
glabel("sic dolor ...", cont=exp)
addSpring(exp)                          # tightens up Qt display

## * In tcltk labels are left aligned
## * In RGtk2 labels are centered
## * In Qt  use addSpring to adjust weight given to children. LEft aligned by default.

visible(w) <- TRUE
