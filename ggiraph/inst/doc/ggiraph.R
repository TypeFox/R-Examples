## ----echo=FALSE, warning=FALSE, message=FALSE----------------------------
library(ggiraph)

# gg theme
gg_th <- theme(axis.ticks = element_line(colour = "gray90"), 
    panel.grid.major = element_line(colour = "gray90"), 
    panel.grid.minor = element_line(colour = "transparent"), 
    panel.background = element_rect(fill = "transparent"))

# geom_point_interactive example
gg_point_1 <- ggplot(mtcars, aes(x = disp, y = qsec, 
		color = wt, tooltip = row.names(mtcars), data_id = row.names(mtcars) ) ) + 
	geom_point_interactive(size=3) + 
  scale_color_gradient(low = "#F3C899", high = "#8C120A") 

gg_point_1 <- gg_point_1 + gg_th
# htmlwidget call
ggiraph(code = {print(gg_point_1)}, 
        tooltip_extra_css = "padding:2px;background:rgba(70,70,70,0.1);color:black;border-radius:2px 2px 2px 2px;",
        hover_css = "fill:#1279BF;stroke:#999999;stroke-width:1pt;r:6pt;cursor:pointer;")


## ----message=FALSE-------------------------------------------------------
library(ggiraph)

dataset <- mtcars
head(dataset)
dataset$tooltip <- row.names(dataset)

# geom_point_interactive example
gg_point_1 <- ggplot(dataset, aes(x = disp, y = qsec, 
		color = wt, tooltip = tooltip ) ) + 
	geom_point_interactive(size=3) + gg_th

# htmlwidget call
ggiraph(code = {print(gg_point_1)}, width = "400px", height = "400px")

## ------------------------------------------------------------------------
dataset$data_id <- tolower(row.names(dataset))

# geom_point_interactive example
gg_point_2 <- ggplot(dataset, aes(x = disp, y = qsec, 
		color = wt, tooltip = tooltip, data_id = data_id ) ) + 
	geom_point_interactive(size=4) + gg_th

# htmlwidget call
ggiraph(code = {print(gg_point_2)}, 
        hover_css = "fill:orange;r:6px;cursor:pointer;")

## ----message=FALSE, warning=FALSE----------------------------------------
crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
crimes$onclick <- sprintf(
  "window.open(\"%s%s\")",
  "http://en.wikipedia.org/wiki/",
  as.character(crimes$state)
)

gg_point_3 <- ggplot(crimes, aes(x = Murder, y = Assault, size = UrbanPop, colour = Rape )) + 
  geom_point_interactive(
    aes( data_id = state, tooltip = state, onclick = onclick ) ) + 
  scale_colour_gradient(low = "#999999", high = "#FF3333") + gg_th

ggiraph(code = print(gg_point_3), 
        hover_css = "fill-opacity:.3;cursor:pointer;")

