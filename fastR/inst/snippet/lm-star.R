require(faraway); require(grid)
star.plot1 <- xyplot(light~temp,star)
hotstar <- star[star$temp > 3.7,]      # select all but 4 coolest stars
star.model1 <- lm(light~temp,star)
star.model2 <- lm(light~temp,hotstar)
star.plot2 <- xyplot(light~temp,star,
    panel = function(x,y,...){
        panel.abline(reg=star.model1, lwd=2, lty=1,
            col=trellis.par.get('superpose.line')$col[2])
        panel.abline(reg=star.model2, lwd=2, lty=1,
            col=trellis.par.get('superpose.line')$col[1])
        panel.xyplot(x,y,...)
        ids <- which(star$temp < 4.0)
        grid.text(x=x[ids] + 0.04, y=y[ids],
            as.character(ids),
            default.units="native",gp=gpar(cex=0.7))
    })
