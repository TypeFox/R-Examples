chalandise.DEP0<-function (CP, inc = 0.25,arrondi=FALSE) 
{
    if(arrondi){CP <- floor(CP/1000)}
    tCP.round <- table(CP)
    CP.coord.x <- XY.DEP$x[as.numeric(names(tCP.round))]
    CP.coord.y <- XY.DEP$y[as.numeric(names(tCP.round))]
    map("france")
    symbols(x = CP.coord.x, y = CP.coord.y, circles = tCP.round, 
        inches = inc, add = TRUE, bg =rev(brewer.pal(3,name="PuRd"))[1])
}

