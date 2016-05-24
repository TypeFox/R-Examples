cooks20x<-
function (lmfit) 
## Cooks Distance plot
{
    y <- cooks.distance(lmfit)
    show.r <- order(-y)[1:3]
    plot(1:length(y), y, type = "h", main = "Cook's Distance plot", 
        xlab = "observation number", ylab = "cook's distance")
    text(show.r, y[show.r] + 0.4 * 0.75 * strheight(" "), show.r)
}

