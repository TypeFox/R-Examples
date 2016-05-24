pejrup.plot <-
function (x = x, lang = lang, axis.labels = axis.labels, show.names = show.names, 
show.lines = show.lines, show.legend = show.legend,show.grid = show.grid, 
cex.axis = cex.axis, cex.names = cex.names, col.labels= col.labels,
col.axis = col.axis, col.names = col.names, col.lines = col.lines, col.grid = col.grid, 
lty.grid = par("lty"))
{
 at <- seq(0.1,0.9,by=0.1)
 tick.labels <- list(l = c("90%","","","","50%","","","","10%"),
r = c("","80%","","","50%","","","20%",""),
b = c("10%","","","","50%","","","","90%"))
 if (is.null(axis.labels)) axis.labels <- colnames(x)[1:3]

 sin60 <- sin(pi/3)
 bx <- at
 by <- rep(0, 9)

 ly <- at * sin60
 lx <- at * 0.5

 rx <- at * 0.5 + 0.5
 ry <- rev(at * sin60)

 if (show.grid) 
 {
  segments(bx, by, lx, ly, lty = lty.grid, col = col.grid)
  segments(lx, ly, rev(rx), rev(ry), lty = lty.grid, col = col.grid)
  segments(rx, ry, bx, by, lty = lty.grid, col = col.grid)
 }

 par(xpd = TRUE)
 
 par(srt = 303) 
 xoffset <- 0.022 
 yoffset <- 0.017 
 text(rx + xoffset, ry + yoffset, tick.labels$r, cex= cex.axis,adj=0.5, col= col.axis)

 par(srt = 0)
 xoffset <- 0.008
 text(bx + xoffset, by - 0.03, rev(tick.labels$b), cex= cex.axis,adj=0.5, col= col.axis)

 x1 <- c(0, 0, 0.5)
 x2 <- c(1, 0.5, 1)
 y1 <- c(0, 0, sin60)
 y2 <- c(0, sin60, 0)
 segments(x1, y1, x2, y2, col = col.lines)

 if (show.lines)
 {
  triax.segments <- function(h1, h3, t1, t3, col)
   segments(1 - h1 - h3/2, h3 * sin(pi/3), 1 - t1 - t3/2, t3 * sin(pi/3), col = col)
  h1 <- c(0, 0, 0, 10, 50, 90)/100
  h3 <- c(50, 80, 20, 0, 0, 0)/100
  t1 <- c(100, 100, 100, 10, 50, 90)/100
  t3 <- c(0, 0, 0, 90, 50, 10)/100
  triax.segments(h1, h3, t1, t3, col.lines)
 }

 if (show.names)
 {
  if(show.legend==FALSE)
   warning ("if is TRUE show.names show.legend should also be TRUE")
  xpos <- c(0.57, 0.69, 0.85, 0.98, 0.018, 0.135, 0.32, 0.46) 
  ypos <- c(0.92, 0.68, 0.37, 0.12, 0.08, 0.32, 0.70, 0.975) * sin(pi/3) 
  snames <- c("I", "II", "III", "IV", "A", "B", "C", "D")
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
   legend("topleft",c("Hydrodynamics","I - Low","II - Moderate",
"III - High", "IV - Very high","","Sand","A - 90-100%","B - 50-90%",
"C - 10-50%","D - 00-10%"), bty="n", col = "black",horiz=FALSE, cex=cex.axis)
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
   legend("topleft",c("Hidrodin\u00E2mica","I - Baixa","II - Moderada",
"III - Alta", "IV - Muito alta","","Areia","A - 90-100%","B - 50-90%",
"C - 10-50%","D - 00-10%"),  bty="n", col = "black",horiz=FALSE, cex=cex.axis)
 }
}

