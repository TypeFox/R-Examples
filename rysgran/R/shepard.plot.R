shepard.plot <-
function(x = x, lang = lang, axis.labels = axis.labels, show.names = show.names, 
show.lines = show.lines, show.legend = show.legend,show.grid = show.grid, 
cex.axis = cex.axis, cex.names = cex.names, col.labels= col.labels,
col.axis = col.axis, col.names = col.names, col.lines = col.lines, col.grid = col.grid, 
lty.grid = par("lty"))
{
 at <- seq(0.25, 0.75, by = 0.25)
 tick.labels <- list(l = c("25%","50%","75%"), r = c("25%","50%","75%"), b = c("25%","50%","75%"))
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

 segments(bx1, by1, bx2, by2, col = col.lines)
 segments(lx1, ly1, lx2, ly2, col = col.lines)
 segments(rx1, ry1, rx2, ry2, col = col.lines)

 if (show.lines)
 {
  triax.segments <- function(h1, h3, t1, t3, col)
   segments(1 - h1 - h3/2, h3 * sin(pi/3), 1 - t1 - t3/2, t3 * sin(pi/3), col = col)
  h1 <- c(75, 25, 0, 50, 0, 50, 50, 50, 25, 12.5, 75, 12.5)/100
  h3 <- c(0, 75, 25, 0, 50, 50, 25, 25, 25, 12.5, 12.5, 75)/100
  t1 <- c(75, 0, 25, 33.33, 33.33, 33.33, 25, 25, 25, 25, 50, 25)/100
  t3 <- c(25, 75, 0, 33.33, 33.33, 33.33, 25, 50, 50, 25, 25, 50)/100
  triax.segments(h1, h3, t1, t3, col.lines)
 }

 if (show.names)
 {
  if(show.legend==FALSE)
   warning ("if is TRUE show.names show.legend should also be TRUE")
  xpos <- c(0.5, 0.4, 0.6, 0.5, 0.25, 0.45, 0.55, 0.75, 0.12, 0.35, 0.65, 0.88) 
  ypos <- c(0.85, 0.55, 0.55, 0.4, 0.31, 0.3, 0.3, 0.31, 0.1, 0.12, 0.12, 0.1) * sin(pi/3) 
  snames <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12") 
  text(xpos, ypos, snames, col = col.names, cex = cex.names)
 }

 xpos <- c(-0.04, 1.04, 0.5)
 ypos <- c(-0.03, -0.03, 1.04) * sin(pi/3) 
 snames <- axis.labels 
 text(xpos, ypos, snames, col = col.axis, cex = cex.axis)

 if (show.legend)
 {
  if(show.names==FALSE)
   warning ("if is TRUE show.legend show.names should also be TRUE")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")
   legend("topleft",c("1 - Clay","2 - Sandy clay","3 - Silty clay",
"4 - Sandy silty clay","5 - Clayey sand","6 - Silty clayey sand",
"7 - Sandy clayey silt","8 - Clayey silt","9 - Sand",
"10 - Silty sand","11 - Sandy silt","12 - Silt"),  bty="n",
col = "black" ,horiz=FALSE, cex=cex.axis)
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
   legend("topleft",c("1 - Argila","2 - Argila arenosa","3 - Argila s\u00EDltica",
"4 - Argila s\u00EDltico-arenosa","5 - Areia argilosa","6 - Areia s\u00EDltico-argilosa",
"7 - Silte argilo-arenoso","8 - Silte argiloso","9 - Areia",
"10 - Areia s\u00EDltica","11 - Silte arenoso","12 - Silte"),bty="n",
          col = "black" ,horiz=FALSE, cex=cex.axis)
 }
}

