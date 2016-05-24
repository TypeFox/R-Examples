
## extracted and simplified from lattice:::plot.trellis

panel.key <-
    function(text, ..., corner = c(0, 1), x = corner[1], y = corner[2])
{
    key <- simpleKey(text, ...)
    key.gf <- draw.key(key, draw = FALSE)
    vp <- viewport(x = unit(x, "npc") + unit(0.5 - corner[1], "grobwidth", list(key.gf)),
                   y = unit(y, "npc") + unit(0.5 - corner[2], "grobheight", list(key.gf)))
    pushViewport(vp)
    grid.draw(key.gf)
    upViewport()
}
