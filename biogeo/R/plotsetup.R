plotsetup <-
function (xi, yi) 
{
    graphics.off()
    dev.new(width=xi, height=yi, xpos = 10,noRStudioGD = TRUE)
    d <- dev.size(units = "px")
    d10 <- (d[1] + 10)
    d2 <- ceiling(d[2]/2)
    dev.new(width=xi, height=yi, xpos = (d[1] + 154),noRStudioGD = TRUE)
    dev.new(width=1.5, height=1.7, xpos = d10, ypos = 0,noRStudioGD = TRUE)
}
