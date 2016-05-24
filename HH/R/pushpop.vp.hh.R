"push.vp.hh" <-
function(scale=100) {
  c.v <- current.viewport()
  ##  print(current.viewport()[c("clip","xscale","yscale","width","height")])
  c.v$clip <- TRUE
  c.v$height <- unit(scale,"npc")
  c.v$width <- unit(scale,"npc")         
  pushViewport(c.v)
  c.v$clip <- FALSE
  c.v$height <- unit(1/scale,"npc")
  c.v$width <- unit(1/scale,"npc")         
  pushViewport(c.v)
}

"pop.vp.hh" <-
function() {
  popViewport()
  popViewport()
}
