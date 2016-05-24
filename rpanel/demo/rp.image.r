click.capture <- function(panel, x, y) 
{
print(x)
print(y)
print(panel)
  if (!is.null(panel$x))
  { 
    rp.line(panel, "gulls.image", panel$x, panel$y, as.numeric(x), as.numeric(y), color="red", width=3, id = "current")
  }
  panel$x <- x
  panel$y <- y
  panel
}

gulls.panel <- rp.control()
image.file <- file.path(system.file(package = "rpanel"), "images", "gulllmks.gif")
gulls.image <- rp.image(gulls.panel, image.file, name="gulls.image", action=click.capture, pos=list(row=0, column=0))
rp.button(gulls.panel, function(panel) { rp.clearlines(gulls.panel, "gulls.image"); panel }, "Remove lines", pos=list(row=1, column=0))