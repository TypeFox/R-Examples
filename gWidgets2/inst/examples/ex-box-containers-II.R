if(interactive()) {

  ## basic container:
  w <- gwindow("Box container examples", visible=FALSE)
  g <- ggroup(cont=w, horizontal=FALSE)

  ## simple horizontal box
  g1 <- ggroup(cont=g)
  sapply(c("one", "two", "three"), function(i) gbutton(i, cont=g1))

  ## framed box
  g2 <- gframe("Framed box", cont=g)
  sapply(c("one", "two", "three"), function(i) gbutton(i, cont=g2))

  ## expanding box
  g3 <- gexpandgroup("expanding box", cont=g)
  sapply(c("one", "two", "three"), function(i) gbutton(i, cont=g3))
  visible(g3) <- TRUE


  ## visible
  g5 <- ggroup(cont=g)
  l <- glabel("click button to hide/show label", cont=g5)
  gbutton("hide/show label", cont=g5, handler=function(h,...) {
    l <- g5[1]
    visible(l) <- !visible(l)
  })

  ## delete
  g6 <- ggroup(cont=g)
  l <- glabel("click button to delete label", cont=g6)
  gbutton("delete label", cont=g6, handler=function(h,...) {
    delete(g6, l)
    enabled(h$obj) <- FALSE
  })

  ## parent, child
  identical(l$parent[1], l)             # true
  

  visible(w) <- TRUE


  ## anchor, expand, fill
  w <- gwindow("expand, anchor, fill", visible=FALSE)
  g <- ggroup(cont=w, horizontal=FALSE)

  ## no expand
  g1 <- ggroup(cont=g, expand=TRUE)
  glabel("no expand", cont=g1)

  ## no expand
  g1 <- ggroup(cont=g, expand=TRUE)

  ## expand, anchor
  g1 <- ggroup(cont=g, expand=TRUE)
  b1 <- glabel("expand, c(-1,1)", cont=g1, expand=TRUE, anchor=c(-1,1), fill=FALSE)
  
  ## expand, fill=TRUE
  g1 <- ggroup(cont=g, expand=TRUE)
  b1 <- glabel("expand=fill=TRUE", cont=g1, expand=TRUE, fill=TRUE)

  ## expand, fill=x
  g1 <- ggroup(cont=g, expand=TRUE)
  b1 <- glabel("expand=TRUE, fill='x'", cont=g1, expand=TRUE, fill="x")


  ## expand, fill=y
  g1 <- ggroup(cont=g, expand=TRUE)
  b1 <- glabel("expand=TRUE, fill='y'", cont=g1, expand=TRUE, fill="y")

  visible(w) <- TRUE
  
}
