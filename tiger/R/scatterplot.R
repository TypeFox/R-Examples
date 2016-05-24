scatterplot <- function(measures, show.measures=1:num.measures){

#    require(hexbin)
    num.measures <- NCOL(measures)
    panel.hist <- function(x, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y,  ...)
    }

    panel.hex <- function(x,y,...){
        panel.hexbinplot(x,y,...)
        #browser()
        #rahmen mit lrect. Aber weis nicht, wie 
        #die richtigen Ränder für trellis
    }
    no.plot <- function(x,...){}

#Scatterplot und statistik
    splom(~measures[,show.measures], data=measures[,show.measures], lower.panel=panel.hex, upper.panel=no.plot)
}
