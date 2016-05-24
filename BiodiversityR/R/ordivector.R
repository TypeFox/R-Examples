`ordivector` <-
function(ordiplot,spec,lty=2,...) {
    speciescoord <- scores(ordiplot, display="species")
    speciesselect <- speciescoord[rownames(speciescoord)==spec]
    sitescoord <- scores(ordiplot, display="sites")
    b1 <- speciesselect[2]/speciesselect[1]
    b2 <- -1/b1
    calc <- array(dim=c(nrow(sitescoord),3))
    calc[,3] <- sitescoord[,2]-b2*sitescoord[,1]
    calc[,1] <- calc[,3]/(b1-b2)
    calc[,2] <- b1*calc[,1]
    for (i in 1:nrow(sitescoord)) {
        graphics::segments(sitescoord[,1], sitescoord[,2], calc[,1], calc[,2], lty=lty)
    }
    graphics::abline(0, b1, lty=lty)
    graphics::arrows(0, 0, speciesselect[1], speciesselect[2],lty=1,...)
}

