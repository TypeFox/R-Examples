options(guiToolkit = "tcltk"); require(tkrplot)
window <- gwindow("How to embed tkrplot", visible = FALSE)
group <- ggroup(cont = window, horizontal = FALSE)
bb <- 1
img <- tkrplot(getToolkitWidget(group), 
               fun = function() plot(1:20,(1:20)^bb))
add(group, img)
f <- function(...) {
  b <- svalue(slider)
  print(b)
  if (b != bb) {
    bb <<- b
    tkrreplot(img)
  }
}
slider <- gslider(from = 0.05, to = 2, by = 0.05, cont = group, 
                  handler = f, expand = TRUE)
visible(window) <- TRUE
