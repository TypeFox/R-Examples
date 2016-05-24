genome.plot = function(mydata, style=1, type="h", sig.line=c(4, -4), sig.color=c("red", "red"), ...)
{
	mydata = na.omit(mydata)
    t.chr = as.character(mydata$chr)
    t.chr = ifelse(t.chr %in% c("1","2","3","4","5","6","7","8","9"), paste("0", t.chr, sep=""), t.chr)
    t.chr = factor(t.chr)
    chr.name = levels(t.chr)
    chr.name = ifelse(substr(chr.name, 1,1) == "0", substring(chr.name, 2), chr.name)
    levels(t.chr) = chr.name
    n.chr = length(chr.name)
    mins = as.vector(tapply(mydata$pos, t.chr, FUN="min"))
    maxs = as.vector(tapply(mydata$pos, t.chr, FUN="max"))

    if (style == 1){
    	xyplot(y ~ pos|t.chr, data=mydata, xlab = "Chromosome", type=type, ..., 
               layout=c(n.chr,1), scales=list(x=list(relation="free", draw=FALSE)), 
               par.settings = list(layout.widths=list(panel=maxs-mins), axis.line=list(lwd=0.1),                
                                   strip.border=list(lwd=0.1)),
               strip = function(..., bg, par.strip.text) 
                               strip.default(..., bg="pink",  par.strip.text=list(cex=0.75)),
               abline=list(h=sig.line, col=sig.color))
    }
    else if(style == 2){
    	xyplot(y ~ pos|t.chr, data=mydata, xlab = "Chromosome", type=type, ..., 
               layout=c(n.chr,1), strip=FALSE,
               scales=list(x=list(relation="free", tck=c(0,0),
                           at=as.vector((maxs+mins)/2, mode="list"), 
                           labels=as.vector(chr.name, mode="list"))), 
               par.settings = list(layout.widths=list(panel=maxs-mins), axis.line=list(lwd=0.1)), 
               abline=list(h=sig.line, col=sig.color))
    }
    else stop("The value of shape should be either 1 or 2")
}



