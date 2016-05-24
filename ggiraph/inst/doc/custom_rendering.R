## ----message=FALSE-------------------------------------------------------
library(ggiraph)
dataset <- mtcars
dataset$carname <- row.names(dataset)
gg_point_1 <- ggplot(dataset, aes(x = disp, y = qsec, tooltip = carname, data_id = carname, color= wt) ) + 
	geom_point_interactive(size=3)

# htmlwidget call
ggiraph(code = {print(gg_point_1)}, tooltip_offx = 20, tooltip_offy = -10 )

## ------------------------------------------------------------------------
tooltip_css <- "background-color:transparent;font-style:italic;"

## ----message=FALSE-------------------------------------------------------
ggiraph(code = {print(gg_point_1)}, tooltip_extra_css = tooltip_css )

## ------------------------------------------------------------------------
tooltip_css <- "background-color:white;font-style:italic;padding:10px;border-radius:10px 20px 10px 20px;"

ggiraph(code = {print(gg_point_1)}, tooltip_extra_css = tooltip_css, tooltip_opacity = .75 )

## ------------------------------------------------------------------------
ggiraph(code = {print(gg_point_1)}, hover_css = "fill:red;" )

## ------------------------------------------------------------------------
ggiraph(code = {print(gg_point_1)}, hover_css = "fill:red;r:6pt" )

