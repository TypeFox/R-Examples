draw.rect <- function (col = "grey", lty = 2, ...) 
{
    axis(3, -49.5:69.5, paste(rep(LETTERS[1:12],each=10),0:9,sep=''), 
         tick = F, line = -0.75)
    axis(4, seq(36.25, 85.25, by = 0.5), sprintf("%02d",c(1:99)), 
         tick = F,las = 1, line = -0.75)
    abline(v = -49:70, col = col, lty = lty, ...)
    abline(h = seq(36, 90, by = 0.5), col = col, lty = lty, 
        ...)
}