flemming.plot <-
  function(x = x, lang = lang, axis.labels = axis.labels, show.names = show.names, 
           show.lines = show.lines, show.legend=show.legend, show.grid = show.grid, 
           cex.axis = cex.axis, cex.names = cex.names, col.labels= col.labels,
           col.axis = col.axis, col.names = col.names, col.lines = col.lines, col.grid = col.grid, 
           lty.grid = par("lty"))
  {
    at <- seq(0.05,0.95,by=0.05)
    tick.labels <- list(l = c("95","","","","75","","","","","50","","","","","25","","","","5"),
                        r = c("","90","","","75","","","","","50","","","","","25","","","10",""),
                        b = c("95","","","","75","","","","","50","","","","","25","","","","5"))
    if (is.null(axis.labels)) axis.labels <- colnames(x)
    
    sin60 <- sin(pi/3)
    bx1 <- at
    bx2 <- at
    by1 <- rep(0, 9)
    by2 <- rep(-0.02 * sin60, 9)
    
    ly1 <- at * sin60
    lx1 <- bx1 * 0.5
    lx2 <- lx1 - 0.015  
    ly2 <- ly1 + 0.008
    
    rx1 <- at * 0.5 + 0.5
    rx2 <- rx1 + 0.015
    ry1 <- rev(ly1)
    ry2 <- rev(ly2) 
    
    if (show.grid) 
    {
      segments(bx1, by1, lx1, ly1, lty = lty.grid, col = col.grid)
      segments(lx1, ly1, rev(rx1), rev(ry1), lty = lty.grid, col = col.grid)
      segments(rx1, ry1, bx1, by1, lty = lty.grid, col = col.grid)
    }
    
    par(xpd = TRUE)
    
    par(srt = 57) 
    xoffset <- 0.03 
    yoffset <- 0.02    
    text(lx1 - xoffset, ly1 + yoffset, tick.labels$l, cex = cex.axis, adj = 0.5, col= col.axis)
    
    par(srt = 303)
    xoffset <- 0.022
    yoffset <- 0.017
    text(rx2 + xoffset, ry1 + yoffset, tick.labels$r, cex = cex.axis, adj=0.5, col= col.axis)
    
    par(srt = 0)
    xoffset <- 0.008
    yoffset <- 0.035
    text(bx1 + xoffset, by1 - yoffset, rev(tick.labels$b), cex = cex.axis, adj=0.5, col= col.axis)
    
    x1 <- c(0, 0, 0.5)
    x2 <- c(1, 0.5, 1)
    y1 <- c(0, 0, sin60)
    y2 <- c(0, sin60, 0)
    segments(x1, y1, x2, y2, col = col.lines)
    
    if (show.lines)
    {
      triax.segments <- function(h1, h3, t1, t3, col)
        segments(1 - h1 - h3/2, h3 * sin(pi/3), 1 - t1 - t3/2, t3 * sin(pi/3), col = col)
      h1 <- c(5, 25, 50, 75, 95, 95, 50 , 50, 75, 75)/100
      h3 <- c(0, 0, 0, 0, 0, 2.5, 45, 5, 18.75, 6.25)/100
      t1 <- c(5, 25, 50, 75, 95, 0, 0, 0, 0, 0)/100
      t3 <- c(95, 75, 50, 25, 5, 50, 90, 10, 75, 25)/100
      triax.segments(h1, h3, t1, t3, col.lines)
    }
    if (show.names)
    {
      if(show.legend==FALSE)
        warning ("if is TRUE show.names show.legend should also be TRUE")
      par(srt = 0)
      xpos <- c(0.027, 0.110, 0.135, 0.225, 0.275, 0.325, 0.375) 
      ypos <- c(0.021, 0.110, 0.040, 0.350, 0.250, 0.150, 0.050) * sin(pi/3)
      snames <- c("S", "A-II", "A-I", "B-IV", "B-III", "B-II", "B-I")
      text(xpos, ypos, snames, col = col.names, cex = cex.names)
      
      xpos <- c(0.345, 0.385, 0.450, 0.525, 0.590, 0.640) 
      ypos <- c(0.625, 0.540, 0.400, 0.250, 0.110, 0.025) * sin(pi/3)
      snames <- c("C-VI", "C-V","C-IV", "C-III", "C-II", "C-I") 
      text(xpos, ypos, snames, col = col.names, cex = cex.names)
      
      xpos <- c(0.450, 0.500, 0.580, 0.680, 0.775, 0.825) 
      ypos <- c(0.800, 0.700, 0.525, 0.325, 0.150, 0.050) * sin(pi/3)
      snames <- c("D-VI", "D-V","D-IV", "D-III", "D-II", "D-I") 
      text(xpos, ypos, snames, col = col.names, cex = cex.names)
      
      par(srt = 303)
      xpos <- c(0.510, 0.572, 0.670, 0.790, 0.890, 0.950) 
      ypos <- c(0.930, 0.800, 0.610, 0.372, 0.172, 0.050) * sin(pi/3)
      snames <- c("E-VI", "E-V","E-IV", "E-III", "E-II", "E-I") 
      text(xpos, ypos, snames, col = col.names, cex = cex.names)
      par(srt = 0)
    }
    
    xpos <- c(-0.04, 1.04, 0.5)
    ypos <- c(-0.03, -0.03, 1.04) * sin(pi/3) 
    snames <- axis.labels 
    text(xpos, ypos, snames, col = col.axis, cex = cex.axis)
    }
