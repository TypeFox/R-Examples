RAM.border <- function(base_size = 12, base_family = "") { 
  # Starts with theme_grey and then modify some parts 
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme( 
    axis.line =         element_blank(),
    axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "grey50", angle = 45, vjust = 1, hjust=1),
    axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "grey50", hjust = 1),
    axis.ticks =        element_line(size=1, colour = "grey50"),
    axis.title.x =      element_text(size = base_size, vjust = 0.5),
    axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.5, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
 
    legend.background = element_rect(colour="white"), 
    legend.key =        element_rect(fill = "grey95", colour = "white"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       element_text(size = base_size * 0.8),
    legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
 
    panel.background =  element_rect(fill = "grey90", colour = NA), 
    panel.border =      element_rect(fill=NA, color="black", size=0.25, linetype="solid"), 
    #panel.grid.major =  element_line(colour = "white"),
     panel.grid.major =  element_blank(),
    #panel.grid.minor =  element_line(colour = "grey95", size = 0.25),
     panel.grid.minor =  element_blank(),
    panel.margin =      unit(0.25, "lines"),
 
    strip.background =  element_rect(fill = "grey80", colour = "black", size=0.25),
    strip.text.x =      element_text(size = base_size * 0.8, face="bold"),
    strip.text.y =      element_text(size = base_size * 0.8, angle = -90, face="bold"),
 
    #aspect.ratio =      1,
    plot.background =   element_rect(colour = NA, fill = "white",
    linetype = "longdash"),
    plot.title =        element_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")

  )   
  
}


RAM.color <- function(base_size = 12, base_family = "") { 
  # Starts with theme_grey and then modify some parts 
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme( 
    axis.line =         element_blank(),
    axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "grey50", angle = 45, vjust = 1, hjust=1),
    axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "grey50", hjust = 1),
    axis.ticks =        element_line(size=1, colour = "grey50"),
    axis.title.x =      element_text(size = base_size, vjust = 0.5),
    axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.5, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
 
    legend.background = element_rect(colour="white"), 
    legend.key =        element_rect(fill = "grey95", colour = "white"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       element_text(size = base_size * 0.8),
    legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
 
    panel.background =  element_rect(fill = "aliceblue", colour = NA), 
    panel.border =      element_rect(fill=NA, color="black", size=0.25, linetype="solid"), 
    #panel.grid.major =  element_line(colour = "white"),
     panel.grid.major =  element_blank(),
    #panel.grid.minor =  element_line(colour = "grey95", size = 0.25),
     panel.grid.minor =  element_blank(),
    panel.margin =      unit(0.25, "lines"),
 
    strip.background =  element_rect(fill = "blue", colour = "black", size=0.25),
    strip.text.x =      element_text(colour="white", size = base_size * 0.8, face="bold"),
    strip.text.y =      element_text(colour="white", size = base_size * 0.8, angle = -90, face="bold"),
 
    #aspect.ratio =      1,
    plot.background =   element_rect(colour = NA, fill = "white",
    linetype = "longdash"),
    plot.title =        element_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")

  )   
  
}

