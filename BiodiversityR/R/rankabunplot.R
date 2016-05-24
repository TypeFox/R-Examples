`rankabunplot` <-
function(xr,addit=F,labels="",scale="abundance",scaledx=F,type="o",xlim=c(min(xpos),max(xpos)),ylim=c(0,max(x[,scale])),specnames=c(1:5),...) {
    x <- xr
    xpos <- 1:nrow(x)
    if (scaledx==T) {xpos <- xpos/nrow(x)*100}
    if (scale=="accumfreq") {type <- "o"}
    if (addit==F) {
        if (scale=="logabun") {
            graphics::plot(xpos, x[,"abundance"], xlab="species rank", ylab="abundance", type=type, bty="l", log="y", xlim=xlim,...)
        }else{
            graphics::plot(xpos, x[,scale], xlab="species rank", ylab=scale, type=type, bty="l", ylim=ylim, xlim=xlim,...)
        }
    }else{
        if (scale=="logabun") {
            graphics::points(xpos,x[,"abundance"], type=type,...)
        }else{
            graphics::points(xpos,x[,scale], type=type,...)
        }
    }
    if (length(specnames) > 0) {
        for (i in specnames) {
            if (scale=="logabun") {
                graphics::text(i+0.5, x[i,"abundance"], rownames(x)[i], pos=4)
            }else{
                graphics::text(i+0.5, x[i,scale], rownames(x)[i], pos=4)
            }
        }
    }
    if (labels!="") {
            if (scale=="logabun") {
                graphics::text(1, x[1,"abundance"], labels=labels, pos=2)
            }else{
                graphics::text(1, x[1,scale], labels=labels, pos=2)
            }
    }
}

