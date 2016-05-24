# Internal - Prepare canvas for charts

.ss.prepCanvas <- function(main = "Six Sigma graph", sub = "My Six Sigma Project",
                           ss.col = c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")){
  #Plot
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(col = ss.col[2], 
                            lwd = 2, 
                            fill = ss.col[5]))
  vp.canvas <- grid::viewport(name = "canvas",
                              width = unit(1,"npc") - unit(6, "mm"),
                              height = unit(1, "npc") - unit(6, "mm"),
                              layout = grid::grid.layout(3, 1,
                                                   heights = unit(c(3, 1, 2),
                                                                  c("lines", "null", "lines"))
                              ))
  grid::pushViewport(vp.canvas)
  grid::grid.rect(gp = grid::gpar(col = "#FFFFFF", 
                            lwd = 0, 
                            fill = "#FFFFFF"))
  
  #Title
  vp.title <- grid::viewport(layout.pos.col = 1, 
                             layout.pos.row = 1, 
                             name = "title")
  grid::pushViewport(vp.title)
  grid::grid.text (main, 
                   gp = grid::gpar(fontsize = 20))
  grid::popViewport()
  
  #Subtitle
  vp.subtitle <- grid::viewport(layout.pos.col = 1,
                                layout.pos.row = 3,
                                name="subtitle")
  grid::pushViewport(vp.subtitle)
  grid::grid.text(sub, 
                  gp = grid::gpar(col = ss.col[1]))
  grid::popViewport()
  
  #Container
  vp.container <- grid::viewport(layout.pos.col = 1,
                           layout.pos.row = 2,
                           name = "container")
  grid::pushViewport(vp.container)
}


