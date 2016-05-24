a.form.ade <-
function(x = 0, y = 0, h=0, b = 1, b2 = 1, side=1, col='Gray50', nslices=100, border=FALSE, density = NULL, angle = 45, density2 = NULL, angle2 = 45 , gradient=FALSE, bcol=1, form='r', fill=NULL, lwd=1, blwd=1, lty=1, blty=1, a.tab){

cubus<-F
zylinder<-F
arect=FALSE
x<-as.numeric(x)
y<-as.numeric(y)
h<-as.numeric(h)
b<-as.numeric(b)
if(form=='z' | form=='Zylinder' | form=='zylinder' | form=='Zyl' | form=='zyl') zylinder<-T
if(form=='c' | form=='cubus' | form=='Cubus' | form=='Cub' | form=='cub'| form=='Cube' | form=='box'| form=='Box') cubus<-T
if(form=='r' | form=='rect' | form=='rectangle' | form=='rec' | form=='R') arect<-T

if(side==2 & zylinder){
cubus=TRUE
zylinder=FALSE
}

################################################################################
################################################################################

# Rect
if(arect){
tograd<-'x'
if(is.character(gradient)){
tograd<-gradient
gradient<-T
}


if(side==1){
if(!border)  bcol<-F
if(!is.null(fill) & !gradient) rect(x-(b/2), y, x+(b/2), y+h, col=fill, border=FALSE      )
if(!is.null(angle) & !is.null(angle2))   rect(x-(b/2), y, x+(b/2), y+h, col=col, border=FALSE,    density=density2,   angle =angle2,  lwd=lwd, lty=lty)
if(!gradient)                  rect(x-(b/2), y, x+(b/2), y+h, col=col, border=bcol, density=density,    angle = angle , lwd=lwd, lty=lty)
if(border)                     rect(x-(b/2), y, x+(b/2), y+h, col=NULL, border=bcol,                lwd=blwd, lty=blty)


if(gradient & length(dim(a.tab))==1){
for(kl in 1:length(col)){
cylindrect(x[kl]-(b/2), y, x[kl]+(b/2), y+h[kl], col=col[kl], gradient=tograd, nslices=nslices, border=bcol)
}
}

if(gradient & length(dim(a.tab))>1) cylindrect(x-(b/2), y, x+(b/2), y+h, col=col, gradient=tograd, nslices=nslices, border=bcol)

if(gradient & !is.null(angle2) )   rect(x-(b/2), y, x+(b/2), y+h, col=col, border=FALSE,    density=density2, angle =angle2,  lwd=lwd, lty=lty)
if(gradient & !is.null(angle)  )   rect(x-(b/2), y, x+(b/2), y+h, col=col, border=bcol,     density=density,  angle = angle , lwd=lwd, lty=lty)

}


if(side==2){
if(!border)  bcol<-F
if(!is.null(fill) & !gradient) rect(y, x-(b/2), y+h, x+(b/2),  col=fill, border=FALSE      )
if(!gradient & !is.null(angle) & !is.null(angle2))                  rect(y, x-(b/2), y+h, x+(b/2),  col=col, border=FALSE,   density=density2,   angle =angle2, lwd=lwd, lty=lty)
if(!gradient)                  rect(y, x-(b/2), y+h, x+(b/2),  col=col, border=bcol, density=density,   angle =angle , lwd=lwd, lty=lty)
if(border)                     rect(y, x-(b/2), y+h, x+(b/2),  col=NULL, border=bcol,                lwd=blwd, lty=blty)
if(gradient)             cylindrect(y, x-(b/2), y+h, x+(b/2),  col=col, gradient=tograd, nslices=nslices, border=bcol)
}

}



################################################################################
################################################################################



if(cubus){


#########################################################
# if verctor of independet numbers
if(length(x)>1){
if(length(y)==1)  y<-rep(y, length(x))
if(length(b2)==1) b2<-rep(b2, length(x))
if(length(h)==1)  h<-rep(h, length(x))
if(length(b)==1)  b<-rep(b, length(x))
if(length(blwd)==1)  blwd<-rep(blwd, length(x))
if(length(lwd)==1)  lwd<-rep(lwd, length(x))
if(length(blty)==1)  blty<-rep(blty, length(x))
if(length(lty)==1)  lty<-rep(lty, length(x))
if(length(density)==1)  density<-rep(density, length(x))
if(length(density2)==1)  density2<-rep(density2, length(x))
if(length(angle)==1)  angle<-rep(angle, length(x))
if(length(angle2)==1)  angle2<-rep(angle2, length(x))
if(length(fill)==1)  fill<-rep(fill, length(x))
if(length(bcol)==1)  bcol<-rep(bcol, length(x))
if(length(col)==1)  col<-rep(col, length(x))
for(i in 1:length(x)){
a.form.ade(x = x[i], y = y[i], h=h[i], b = b[i],  b2 = b2[i], side=side, col=col[i], nslices=nslices, border=border, density = density[i], angle = angle[i], density2 = density2[i], angle2 = angle2[i], gradient=gradient, bcol=bcol[i], form=form, fill=fill[i], lwd=lwd[i], blwd=blwd[i], lty=lty[i], blty=blty[i], a.tab=a.tab)

}
}
#########################################################

if(length(x)==1){
if(is.character(gradient)){
tograd<-gradient
gradient<-T
}


if(side==1){
xrange<- diff(range(par('usr')[1:2]))
yrange<- diff(range(par('usr')[3:4]))
hx<-(xrange/40)*(b2)
hy<-(yrange/40)*(b2)

################
# Colors
col2<- a.coladd.ade(col   , -50)
col3<- a.coladd.ade(col2  , -50)
##############

wd<-(b/2)
y<-y-(hy/2)
x<-x-(hx/2)

x5<-c(x-wd,x-wd+hx,x+hx+wd,x+wd)
y5<-c(y,y+hy,y+hy,y)

x4<-c(x-wd,x-wd,x+hx-wd,x+hx-wd)
y4<-c(y,y+h,y+hy+h,y+hy)

x3<-c(x+wd,x+wd,x+hx+wd,x+hx+wd)
y3<-c(y,y+h,y+h+hy,y+hy)

x2<-c(x-wd,x-wd+hx,x+wd+hx,x+wd)
y2<-c(y+h, (y+h+hy), (y+h+hy), y+h)

x1<-c((x-wd),(x-wd),(x+wd),(x+wd))
y1<-c(y,y+h,y+h,y)

if(!is.null(fill)){
fill2<- a.coladd.ade(fill   , -50)
fill3<- a.coladd.ade(fill2  , -50)
polygon(x5,y5, col=fill2, lty=blty,  border=FALSE, lwd=1)
polygon(x4,y4, col=fill3, lty=blty,  border=FALSE, lwd=1)
polygon(x3,y3, col=fill3, lty=blty,  border=FALSE, lwd=1)
polygon(x2,y2, col=fill2, lty=blty,  border=FALSE, lwd=1)
polygon(x1,y1, col=fill , lty=blty,  border=FALSE, lwd=1)
}


polygon(x5,y5, col=col2, lty=blty,  border=border, lwd=blwd)
polygon(x4,y4, col=col3, lty=blty,  border=border, lwd=blwd)
polygon(x3,y3, col=col3, lty=blty,  border=border, lwd=blwd)
polygon(x2,y2, col=col2, lty=blty,  border=border, lwd=blwd)
polygon(x1,y1, col=col , lty=blty,  border=border, lwd=blwd)


if(!is.null(density)){
polygon(x5,y5, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x4,y4, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x3,y3, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x2,y2, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x1,y1, col=col3 , lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
}

if(!is.null(density2)){
polygon(x5,y5, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x4,y4, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x3,y3, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x2,y2, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x1,y1, col=col3 , lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
}
}



if(side==2){
yrange<- diff(range(par('usr')[1:2]))
xrange<- diff(range(par('usr')[3:4]))
hx<-(xrange/40)*(b2)
hy<-(yrange/40)*(b2)

################
# Colors
col2<- a.coladd.ade(col   , -50)
col3<- a.coladd.ade(col2  , -50)
##############

wd<-(b/2)
y<-y-(hy/2)
x<-x-(hx/2)

y5<-c(x+wd, x-wd,x-wd+hx,x+hx+wd)
x5<-c(y, y, y+hy,y+hy)

y4<-c(x+hx-wd, x-wd,x-wd,x+hx-wd)
x4<-c(y+hy, y,y+h,y+hy+h)

y3<-c(x+hx+wd, x+wd, x+wd, x+hx+wd)
x3<-c(y+hy, y,y+h,y+h+hy)

y2<-c(x+wd,x-wd ,x-wd+hx,x+wd+hx)
x2<-c(y+h, y+h , (y+h+hy), (y+h+hy))

y1<-c((x-wd),(x+wd),(x+wd),(x-wd))
x1<-c(y, y, y+h,y+h)

if(!is.null(fill)){
fill2<- a.coladd.ade(fill   , -50)
fill3<- a.coladd.ade(fill2  , -50)
polygon(x5,y5, col=fill2, lty=blty,  border=FALSE, lwd=1)
polygon(x4,y4, col=fill3, lty=blty,  border=FALSE, lwd=1)
polygon(x3,y3, col=fill3, lty=blty,  border=FALSE, lwd=1)
polygon(x2,y2, col=fill2, lty=blty,  border=FALSE, lwd=1)
polygon(x1,y1, col=fill , lty=blty,  border=FALSE, lwd=1)
}


polygon(x5,y5, col=col2, lty=blty,  border=border, lwd=blwd)
polygon(x4,y4, col=col3, lty=blty,  border=border, lwd=blwd)
polygon(x3,y3, col=col3, lty=blty,  border=border, lwd=blwd)
polygon(x2,y2, col=col2, lty=blty,  border=border, lwd=blwd)
polygon(x1,y1, col=col , lty=blty,  border=border, lwd=blwd)


if(!is.null(density)){
polygon(x5,y5, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x4,y4, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x3,y3, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x2,y2, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
polygon(x1,y1, col=col3 , lty=lty,  border=FALSE, lwd=lwd,  density=density, angle=angle )
}

if(!is.null(density2)){
polygon(x5,y5, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x4,y4, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x3,y3, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x2,y2, col=col3, lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
polygon(x1,y1, col=col3 , lty=lty,  border=FALSE, lwd=lwd,  density=density2, angle=angle2 )
}
}

#######################################################################################################
#######################################################################################################
}
}



#######################################################################################################
#######################################################################################################

if(zylinder){



#########################################################
# if verctor of independet numbers
if(length(x)>1){
if(length(y)==1)  y<-rep(y, length(x))
if(length(b2)==1) b2<-rep(b2, length(x))
if(length(h)==1)  h<-rep(h, length(x))
if(length(b)==1)  b<-rep(b, length(x))
if(length(blwd)==1)  blwd<-rep(blwd, length(x))
if(length(lwd)==1)  lwd<-rep(lwd, length(x))
if(length(blty)==1)  blty<-rep(blty, length(x))
if(length(lty)==1)  lty<-rep(lty, length(x))
if(length(density)==1)  density<-rep(density, length(x))
if(length(density2)==1)  density2<-rep(density2, length(x))
if(length(angle)==1)  angle<-rep(angle, length(x))
if(length(angle2)==1)  angle2<-rep(angle2, length(x))
if(length(fill)==1)  fill<-rep(fill, length(x))
if(length(bcol)==1)  bcol<-rep(bcol, length(x))
if(length(col)==1)  col<-rep(col, length(x))

for(i in 1:length(x)){
a.form.ade(x = x[i], y = y[i], h=h[i], b = b[i],  b2 = b2[i], col=col[i], nslices=nslices, border=border, density = density[i], angle = angle[i], density2 = density2[i], angle2 = angle2[i], gradient=gradient, bcol=bcol[i], form=form, fill=fill[i], lwd=lwd[i], blwd=blwd[i], lty=lty[i], blty=blty[i], a.tab=a.tab)
}


}
#########################################################


if(length(x)==1){


if(is.character(gradient)){
tograd<-gradient
gradient<-T
}
xrange<- diff(range(par('usr')[1:2]))
yrange<- diff(range(par('usr')[3:4]))
#hx<-(xrange/9)*  b
hy<-(yrange/20)* (b2)
hx<-(b/2)



################################################################################
##                                                                            ##
sub.ellipse.ade<-function(x = 0, y = 0, h=0, hx = 1, hy = hx, col='Gray20', nslices=100, border=border, bcol=1, gradient=gradient){
library(plotrix)

if(gradient) col<- a.coladd.ade(col, -50)
theta <- 0
hlaxa <- hx
hlaxb <- hy
if(gradient) col2<- a.coladd.ade(col, -50)
xc <- x
yc <- y




################################################################################
##  Nur Ellipse funktion                                                      ##
##----------------------------------------------------------------------------##
ellipse <-function(hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, npoints = 100, col='Gray80', bcol=1, toplot=TRUE, col2=NULL, border=TRUE)
{
        xp <- NULL
        yp <- NULL
        for(i in 0:npoints) {
                a <- (2 * pi * i)/npoints
                x <- hlaxa * cos(a)
                y <- hlaxb * sin(a)
                if(theta != 0) {
                        alpha <- angle(x, y)
                        rad <- sqrt(x^2 + y^2)
                        x <- rad * cos(alpha + theta)
                        y <- rad * sin(alpha + theta)
                }
                xp <- c(xp, x + xc)
                yp <- c(yp, y + yc)
        }
        if(border) lines(xp, yp, type = "l", col=col2, lwd=3)
        if(!is.null(col) & toplot) {
          polygon(xp, yp, border = NA, col = col)
        }
        k<-1
        if(length(col)>1){
        if(!is.null(col2))
        k<-1
        for(i in 1:length(col)){
        k<-k+0.2
        xpw<- xp-mean(xp)
        ypw<- yp-mean(yp)
        polygon((xpw/(k))+mean(xp), (ypw/(k))+mean(yp), border = NA, col = col[i])

        }
        }
        invisible()

        return(list(xp, yp))
}
##____________________________________________________________________________##
################################################################################

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
##----------------------------------------------------------------------------##
angle <- function(x, y){
if(x > 0) { atan(y/x) }
 else { if(x < 0 & y != 0) { atan(y/x) + sign(y) * pi }
  else { if(x < 0 & y == 0) { pi }
   else { if(y != 0) { (sign(y) * pi)/2 }
    else { NA } } } }
}
##____________________________________________________________________________##
################################################################################

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
## cylinderect  # Hier nur Color vorbereitun                                  ##
##----------------------------------------------------------------------------##
cylinder.ade <-function(x1, y1, x2, y2, col, border = NA, gradient = "x", nslices = 50)
{
    # Nur Color
    rgbval <- col2rgb(col)/255
    maxrgb <- max(rgbval)
    cols <- matrix(0, nrow = 3, ncol = 6)


    for (i in 1:3) {
        if (rgbval[i] == maxrgb) delta <- 1 - rgbval[i]
        else delta <- (1 - 0.2 * (maxrgb - rgbval[i])/sum(maxrgb -
            rgbval)) - rgbval[i]
        cols[i, ] <- c(rgbval[i] + 0.3 * delta, rgbval[i] + 0.6 *
            delta, rgbval[i] + 0.9 * delta, rgbval[i] + 0.6 *
            delta, rgbval[i] + 0.3 * delta, rgbval[i])
    }


    col<-gradient.rect.ade(x1, y1, x2, y2, cols[1, ], cols[2, ], cols[3, ], gradient = gradient, nslices = nslices, border = border)

    invisible(cols)
    return(col)
}
##____________________________________________________________________________##
################################################################################

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
## Gradient Rects                                                             ##
##----------------------------------------------------------------------------##
gradient.rect.ade<- function(x1, y1, x2, y2, reds, greens, blues, col = NULL, nslices = 50, gradient = "x", border = par("fg"))
{

if (is.null(col)) col <- color.gradient(reds, greens, blues, (nslices/2)) else nslices <- length(col)

x1<- x1[round(length(x1)/2):length(x1)]
x2<- x2[round(length(x2)/2):length(x2)]
x1<-x1[c(-1, -length(x1))]
x2<-x2[c(-1, -2)]

y1<- y1[round(length(y1)/2):length(y1)]
y2<- y2[round(length(y2)/2):length(y2)]
y1<-y1[c(-1, -length(y1))]
y2<-y2[c(-1, -2)]

rect(x1[1], y2[1], x2[length(x2)], y1[length(y1)], col = col, lty = 0, border=NA)

rect(x1, y1, x2, y2, col = col, lty = 0, border=NA)  #<--------  HIER!!!

    invisible(col)
    return(col)
}
##____________________________________________________________________________##
################################################################################



if(is.null(h)) ellipse(hlaxa = hlaxa, hlaxb = hlaxb, theta = theta, xc = xc, yc = yc, npoints = nslices, col=col, col2=bcol, border=border)

# Zylinder darstellung
if(!is.null(h)) {

p1<-ellipse(hlaxa = hlaxa, hlaxb = hlaxb, theta = theta, xc = xc, yc = yc, npoints = nslices, col=col, toplot=FALSE, col2=bcol, bcol=bcol,border=border)

p2<-ellipse(hlaxa = hlaxa, hlaxb = hlaxb, theta = theta, xc = xc, yc = yc+h, npoints = nslices, col=col, toplot=FALSE, col2=bcol, bcol=bcol,border=border)



if(border) {

segments(p1[[1]][1], p1[[2]][1], p2[[1]][1], p2[[2]][1], lwd=3, col=bcol)
segments(p1[[1]][length(p1[[2]])/2], p1[[2]][length(p1[[2]])/2], p2[[1]][length(p1[[2]])/2], p2[[2]][length(p1[[2]])/2], lwd=5, col=bcol)
}


if(gradient){
myfrad<- cylinder.ade(p1[[1]], p1[[2]], p2[[1]], p2[[2]], col=col, border = border, nslices = nslices, gradient = "x")
g<-ellipse(hlaxa = hlaxa, hlaxb = hlaxb, theta = theta, xc = xc, yc = yc+h, npoints = nslices, col=myfrad, toplot=TRUE, col2=bcol, border=border)
}


if(!gradient){
ellipse(hlaxa = hlaxa, hlaxb = hlaxb, theta = theta, xc = xc, yc = yc, npoints = nslices, col=col, toplot=TRUE, col2=bcol, border=border)
rect(p1[[1]][length(p1[[2]])/2], p1[[2]][1], p1[[1]][1], p2[[2]][length(p2[[2]])], col=col, border=NA)
ellipse(hlaxa = hlaxa, hlaxb = hlaxb, theta = theta, xc = xc, yc = yc+h, npoints = nslices, col=col, toplot=TRUE, col2=bcol, border=border)
}


}

}
################################################################################
################################################################################
sub.ellipse.ade(x = x, y = y, h=h, hx = hx, hy = hy, col=col, nslices=nslices, bcol=bcol, border=border, gradient=gradient)
}
}
}
