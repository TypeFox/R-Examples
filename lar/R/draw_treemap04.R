# File: draw_treemap04.R
# Date: March 25, 2014
# Author: richard.zijdeman@iisg.nl
# Last change: 
  # 2014-04-29, replaced subsetting by indexing to prevent NOTE in CRAN checking
  # 2014-04-28, cleaned file syntax for acceptance on CRAN
  # 2014-04-03, added as.character to title value, to avodi issues with 
  #             with special characters in names of title value

draw.lar <- function(data.frame, sex="total") {
  title.value <- as.character(unique(data.frame$ctry.time))
  vp.tree     <- viewport(x = 0.5, y = 0.55, width = 1, height = .9)
  pushViewport(vp.tree)
  
  # first make sure, all 6 elements (Commodified, EitherOr, Non working,
  # Reciprocal, Tributary and Unknown are available in the selected data,
  # and adding them if this was not the case)
  
  temp <- as.data.frame(c("Non working",
                          "Reciprocal",
                          "Tributary",
                          "Commodified",
                          "EitherOr",
                          "Unknown"))
  names(temp) <- "alltxt1.1"
  data.frame <- merge(data.frame, temp,
                by.x = "txt1.1",
                by.y = "alltxt1.1",
                all  = TRUE)
  #rm(temp)
  
  ## Create proper sample from options
  if(sex == "total"){
    data.frame <- data.frame
  }
  else if(sex == "male"){
    data.frame <- data.frame[data.frame$gender == "M", ]
  }
  else if(sex == "female"){
    data.frame <- data.frame[data.frame$gender == "F", ]
  }
  else {
    stop("Invalid sex specfication.")
  }
  
  
  ## Now, draw the treemap
  treemap(data.frame, 
          title           = title.value ,
          index           = c(paste0("txt1.",c(1:2)), 
                             "txt1.3.ext",       # first labrel: 3 levels deep
                             "txt2.1","txt3.1"), # 2nd, 3rd labrel, 1 deep
          vSize           = "total",   
          palette         = c("#00B888", "#4F81BD", "#99E412", "#EF703D", 
                               "#FF89B0","#FFD92F"), 
          vColor          = "txt1.1",
          type            = "categorical",
          sortID          = "sortID2",
          position.legend = "none",
          vp              = vp.tree,
          border.lwds     = c(3,2,1,1,1),
          border.col      = c("black", "black", "black", "white", "white" ), 
          align.labels    = list(c("left", "top"), 
                              c("right", "bottom"),
                              c("left", "top"),
                              c("left", "center"),
                              c("left", "bottom")),
          fontsize.labels = c(0, 11, 11, 11, 11),
          fontface.labels = c("bold", "bold", "italic", "plain", "plain"),
          bg.labels       = 0,
          overlap.labels  = 0, # value 0-1 (0 is less printed labels)
          # inflate.labels = TRUE, # DO NOT USE! (makes labels too huge)
          drop.unused.levels = FALSE)
  popViewport(1)
  
  ## Drawing the legend ####
  vp.legend <- viewport(x = 0.5, y = 0.05, width = 1, height = .1)
  pushViewport(vp.legend)
  
  # Non working # Lime "#99E412"
  vp.col1 <- viewport(x = 0.02, y = 0.5, width = .02, height = 0.4)
  pushViewport(vp.col1)
  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = "#99E412"))
  popViewport(1)
  vp.lab1 <- viewport(x = 0.14, y = 0.5, width = .125, height = 0.3)
  pushViewport(vp.lab1)
  grid.text("Non working", x = 0.3, y = 0.5, gp = gpar(fontsize=16))
  popViewport(1)
  
  # Reciprocal # Orange "#EF703D"
  vp.col2 <- viewport(x = 0.22, y = 0.5, width = .02, height = 0.4)
  pushViewport(vp.col2)
  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = "#EF703D"))
  popViewport(1)
  vp.lab2 <- viewport(x = 0.34, y = 0.5, width = .125, height = 0.3)
  pushViewport(vp.lab2)
  grid.text("Reciprocal", x = 0.2, y = 0.5, gp = gpar(fontsize = 16))
  popViewport(1)
  
  # Tributary # Pink "#FF89B0"
  vp.col3 <- viewport(x = 0.39, y = 0.5, width = .02, height = 0.4)
  pushViewport(vp.col3)
  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = "#FF89B0"))
  popViewport(1)
  vp.lab3 <- viewport(x = 0.50, y = 0.5, width = .125, height = 0.3)
  pushViewport(vp.lab3)
  grid.text("Tributary", x = 0.2,y = 0.5, gp = gpar(fontsize = 16))
  popViewport(1)
  
  # Commodified # Green "#00B888"
  vp.col4 <- viewport(x = 0.54, y = 0.5, width = .02, height = 0.4)
  pushViewport(vp.col4)
  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = "#00B888"))
  popViewport(1)
  vp.lab4 <- viewport(x = 0.67, y = 0.5, width = .125, height = 0.3)
  pushViewport(vp.lab4)
  grid.text("Commodified", x = 0.2,y= 0.5, gp = gpar(fontsize = 16))
  popViewport(1)
  
  # EitherOr # Blue "#4F818BD"
  vp.col5 <- viewport(x = 0.73, y = 0.5, width = .02, height = 0.4)
  pushViewport(vp.col5)
  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = "#4F81BD"))
  popViewport(1)
  vp.lab5 <- viewport(x = 0.83, y = 0.5, width = .125, height = 0.3)
  pushViewport(vp.lab5)
  grid.text("EitherOr", x = 0.2, y = 0.5, gp = gpar(fontsize = 16))
  popViewport(1)
  
  # Unknown # Yellow
  vp.col6 <- viewport(x = 0.86, y = 0.5, width = .02, height = 0.4)
  pushViewport(vp.col6)
  grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = "#FFD92F"))
  popViewport(1)
  vp.lab6 <- viewport(x = 0.97, y = 0.5, width = .125, height = 0.3)
  pushViewport(vp.lab6)
  grid.text("Unknown", x = 0.2, y = 0.5, gp = gpar(fontsize = 16))
  popViewport(1)
  
  popViewport(1)
}

# EOF
