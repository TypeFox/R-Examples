mediation.effect.plot <- function(x, mediator, dv, ylab="Dependent Variable", xlab="Mediator", main="Mediation Effect Plot", pct.from.top.a=.05, pct.from.left.c=.05, arrow.length.a=.05, arrow.length.c=.05, legend.loc="topleft", file="", pch=20,
xlim=NULL, ylim=NULL,
save.pdf = FALSE, save.eps = FALSE, save.jpg = FALSE, ...)
{	
if(getRversion() >= "3.1.0") utils::suppressForeignCheck(package="stats")
if(getRversion() >= "3.1.0") utils::suppressForeignCheck(package="grDevices")
if(getRversion() >= "3.1.0") utils::suppressForeignCheck(package="graphics")
	
# file name and file path for the graph(s) to save, if file="" a file would be saved in the current working directory

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
file <- "mediation.effect.plot"
no.file.name <- TRUE
}
}


Mediation.Results <- mediation(x=x, mediator=mediator, dv=dv, conf.level=.95)

Int_M.X <- Mediation.Results$M.on.X$Regression.Table[1,1]
Slope_M.X <-  Mediation.Results$M.on.X$Regression.Table[2,1]
##########
Int_Y.X <- Mediation.Results$Y.on.X$Regression.Table[1,1]
Slope_Y.X <-  Mediation.Results$Y.on.X$Regression.Table[2,1]
##########
Int_Y.XM <- Mediation.Results$Y.on.X.and.M$Regression.Table[1,1]
Slope_Y.XM.X <-  Mediation.Results$Y.on.X.and.M$Regression.Table[2,1]
Slope_Y.XM.M <-  Mediation.Results$Y.on.X.and.M$Regression.Table[3,1]


M.x <- mean(x)
s.x <- var(x)^.5

s.mediator <- var(mediator)^.5


if(save.pdf==TRUE) pdf(file = paste(file, ".pdf", sep = ""), ...)
if(save.eps == TRUE) jpeg(filename = paste(file, ".eps", sep = ""), ...)
if(save.jpg == TRUE) jpeg(filename = paste(file, ".jpg", sep = ""), ...)
	
plot(dv~mediator, ylab=ylab, xlab=xlab, ylim=ylim, xlim=xlim, pch=pch, ...)

# Horizontal Lines
abline(h=(Int_Y.X + Slope_Y.X*M.x))
abline(h=(Int_Y.X + Slope_Y.X*(M.x + 1)))

# Vertical Lines
abline(v=(Int_M.X + Slope_M.X*M.x))
abline(v=(Int_M.X + Slope_M.X*(M.x + 1)))

if(!is.null(ylim))
{
here.value.a <- min(max(dv), ylim[2])
print(here.value.a)
print(here.value.a*(1-pct.from.top.a))
}
if(is.null(ylim))
{
here.value.a <-max(dv)
}

text(x=((Int_M.X + Slope_M.X*M.x) + .5*Slope_M.X), y=here.value.a*(1-pct.from.top.a), labels="a", pos=3)

arrows(x0=(Int_M.X + Slope_M.X*M.x), y0=here.value.a*(1-pct.from.top.a), x1=(Int_M.X + Slope_M.X*(M.x + 1)), y1=here.value.a*(1-pct.from.top.a), 
code=3, length=arrow.length.a)

here.c1 <- (Int_Y.X + Slope_Y.X*M.x)
here.c2 <- (Int_Y.X + Slope_Y.X*(M.x + 1))

if(!is.null(xlim))
{
here.value.1 <- max(min(mediator), xlim[1])
}
if(is.null(xlim))
{
here.value.1 <-min(mediator)
}

text(x=here.value.1*(1+pct.from.left.c), y=(Int_Y.X + Slope_Y.X*M.x)+.5*Slope_Y.X, labels="c", pos=2)

arrows(x0=here.value.1*(1+pct.from.left.c), y0=min(here.c1, here.c2), x1=here.value.1*(1+pct.from.left.c), y1=max(here.c1, here.c2), code=3, length=arrow.length.c)

# Prediction lines for Y on X and 

These.m.values <- seq(min(mediator), max(mediator), s.mediator/1000)

lines(These.m.values, Int_Y.XM + These.m.values*Slope_Y.XM.M + Slope_Y.XM.X*M.x, lty=1)
lines(These.m.values, Int_Y.XM + These.m.values*Slope_Y.XM.M + Slope_Y.XM.X*(M.x - s.x), lty=3)
lines(These.m.values, Int_Y.XM + These.m.values*Slope_Y.XM.M + Slope_Y.XM.X*(M.x + s.x), lty=2)

legend(legend.loc, legend=c("1 SD Unit Above Mean on IV", "At Mean of IV", "1 SD Unit Below Mean on IV"), lty=c(2, 1, 3), bty="n")

if (save.pdf == TRUE) 
{
dev.off()
if (no.file.name == TRUE) 
print(paste("'mediation.effect.plot.pdf' file saved at the directory", 
getwd()))
}

if (save.eps == TRUE) 
{
dev.off()
if (no.file.name == TRUE) 
print(paste("'mediation.effect.plot.eps' file saved at the directory", 
getwd()))
}

if (save.jpg == TRUE) 
{
dev.off()
if (no.file.name == TRUE) 
print(paste("'mediation.effect.plot.jpg' file saved at the directory", 
getwd()))
}
}

