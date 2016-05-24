library(tcltk)
if (exists("Z") && (class(Z) == "DynamicGraph")) DynamicGraph(frameModels = Z, redraw = TRUE) else warning("Run first, e.g., 'demo(Circle)', quit R with saving, and restart.")
