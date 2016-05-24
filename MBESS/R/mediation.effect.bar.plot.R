mediation.effect.bar.plot <- function(x, mediator, dv, main="Mediation Effect Bar Plot", width=1, left.text.adj=0, right.text.adj=0, rounding=3, file="", save.pdf = FALSE, save.eps = FALSE, save.jpg = FALSE, ...)
{	
Mediation.Results <- mediation(x=x, mediator=mediator, dv=dv, conf.level=.95)

observed.c <- Mediation.Results$Y.on.X$Regression.Table[2,1]

observed.c.prime <- Mediation.Results$Y.on.X.and.M$Regression.Table[2,1]

max.possible.c <- sqrt(var(dv))/sqrt(var(x))

if(observed.c < 0) max.possible.c <- -max.possible.c

if(width < 1)
{
width <- .5*(1-width)
}
if(width > 1)
{
width <- .5*(1+width)
}

if(save.pdf==TRUE)
{
if(save.eps==TRUE) stop("Only one file format for saving figure may be used at a time (you have both PDF and EPS specified).")
if(save.jpg==TRUE) stop("Only one file format for saving figure may be used at a time (you have both PDF and JPG specified).")
}
if(save.eps==TRUE)
{
if(save.jpg==TRUE) stop("Only one file format for saving figure may be used at a time (you have both EPS and JPG specified).")
}


if(save.pdf==TRUE | save.eps==TRUE | save.jpg==TRUE)
{
no.file.name <- FALSE
if (file == "") 
{
file <- "mediation.effect.bar.plot"
no.file.name <- TRUE
}
}

if(save.pdf==TRUE) pdf(file = paste(file, ".pdf", sep = ""), ...)
if(save.eps == TRUE) jpeg(filename = paste(file, ".eps", sep = ""), ...)
if(save.jpg == TRUE) jpeg(filename = paste(file, ".jpg", sep = ""), ...)

plot(c(-2, 2), seq(0, 1), ylab="", xlab="", xaxt="n", yaxt="n", bty="n", type="n", main=main, ...)

segments(x0=-.5*width, y0=0, x1=-.5*width, y1=1)
segments(x0=.5*width, y0=0, x1=.5*width, y1=1)

segments(x0=.5*width, y0=0, x1=-.5*width, y1=0)
segments(x0=.5*width, y0=1, x1=-.5*width, y1=1)

segments(x0=.5*width, y0=observed.c/max.possible.c, x1=-.5*width, y1=observed.c/max.possible.c)
segments(x0=.5*width, y0=observed.c.prime/max.possible.c, x1=-.5*width, y1=observed.c.prime/max.possible.c) 

rect(xleft=-.5*width, ybottom=0, xright=.5*width, ytop=observed.c.prime/max.possible.c, density = 10, angle = 45, border=NA)
rect(xleft=-.5*width, ybottom=observed.c.prime/max.possible.c, xright=.5*width, ytop=observed.c/max.possible.c, density = 10, angle = 135, border=NA)


if(left.text.adj==0) 
{
left.text.adj <- -.5*width - (.5*width/3)
}
if(left.text.adj != 0)
{
left.text.adj <- -.5*width - (.5*width/3) + left.text.adj
}

if(right.text.adj==0) 
{
right.text.adj <- .5*width + (.5*width/20)
}
if(right.text.adj != 0)
{
right.text.adj <- .5*width + (.5*width/20) + right.text.adj
}

use.this <- round(max.possible.c, rounding)
text(x=right.text.adj*1.3, y=1, bquote(paste(plain("max possible"), phantom(x), italic(c)==.(use.this))))

use.this <- round(observed.c, rounding)
text(x=left.text.adj, y=observed.c/max.possible.c, bquote(paste(plain(observed), phantom(x), italic(c)==.(use.this))))

use.this <- round(observed.c.prime, rounding)
text(x=left.text.adj, y=observed.c.prime/max.possible.c, bquote(paste(plain(observed), phantom(x), italic(c), phantom(x), plain(prime)==.(use.this))))

use.this <- round(observed.c - observed.c.prime, rounding)
text(x=right.text.adj, y=observed.c/max.possible.c - observed.c.prime/max.possible.c, bquote(italic(ab)==.(use.this)))

segments(x0=right.text.adj*.6, y0=observed.c/max.possible.c, x1=right.text.adj*.6, y1=observed.c.prime/max.possible.c)

segments(x0=right.text.adj*.6, y0=observed.c/max.possible.c, x1=right.text.adj*.55, y1=observed.c/max.possible.c)
segments(x0=right.text.adj*.6, y0=observed.c.prime/max.possible.c, x1=right.text.adj*.55, y1=observed.c.prime/max.possible.c)


text(x=right.text.adj*.8, y=0, "zero")


if (save.pdf == TRUE) 
{
dev.off()
if (no.file.name == TRUE) 
print(paste("'mediation.effect.bar.plot.pdf' file saved at the directory", 
getwd()))
}

if (save.eps == TRUE) 
{
dev.off()
if (no.file.name == TRUE) 
print(paste("'mediation.effect.bar.plot.eps' file saved at the directory", 
getwd()))
}

if (save.jpg == TRUE) 
{
dev.off()
if (no.file.name == TRUE) 
print(paste("'mediation.effect.bar.plot.jpg' file saved at the directory", 
getwd()))
}

}

