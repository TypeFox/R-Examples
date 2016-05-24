library(tcltk)
if (exists("U") && (class(U) == "DynamicGraph")) DynamicGraph(frameModels = U, redraw = TRUE) else warning("Run first 'demo(addModelsAndViews)', quit R with saving, and restart.")
