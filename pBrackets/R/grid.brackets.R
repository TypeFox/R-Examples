grid.brackets <-
function(x1, y1, x2, y2, h=NULL, ticks=0.5, curvature=0.5, type=1, col=1, lwd=1, lty='solid')
{

if(is.null(h)) h <- 0.05
if(!is.unit(x1) | !is.unit(y1) | !is.unit(x2) | !is.unit(y2) | !is.unit(h)){
x1 <- unit(x1, 'native')
x2 <- unit(x2, 'native')
y1 <- unit(y1, 'native')
y2 <- unit(y2, 'native')
h  <- unit(h,  'npc')
}

if(!is.numeric(curvature)) stop('curvature must be numeric')
if(!is.numeric(type))      stop('type must be a integer, 1 to 5')
if(length(ticks)==1)       if(is.na(ticks)) ticks<- NULL
if(!is.numeric(ticks) & !is.null(ticks))     stop('ticks must be numeric or NULL')
if(length(ticks)>1){
if(any(duplicated(abs(ticks)))) stop('duplicated ticks')
}
if(curvature<0) curvature<- 0
if(curvature>1) curvature<- 1


myangle <- 180*atan2((convertY(y2, unitTo='inches', valueOnly = TRUE)-convertY(y1, unitTo='inches', valueOnly = TRUE)),(convertX(x2, unitTo='inches', valueOnly = TRUE)-convertX(x1, unitTo='inches', valueOnly = TRUE)))/pi
mywidth <- sqrt((convertX(x2, unitTo='inches', valueOnly = TRUE)-convertX(x1, unitTo='inches', valueOnly = TRUE))^2+(convertY(y2, unitTo='inches', valueOnly = TRUE)-convertY(y1, unitTo='inches', valueOnly = TRUE))^2)
mywidth <- convertWidth(unit(mywidth, units='inches'),  unitTo='npc', valueOnly = TRUE)
myheight <- convertHeight(h,  unitTo='npc', valueOnly = TRUE)


################################################################################
x1<-convertX(x1, unitTo='npc', valueOnly = TRUE)
x2<-convertX(x2, unitTo='npc', valueOnly = TRUE)
y1<-convertY(y1, unitTo='npc', valueOnly = TRUE)
y2<-convertY(y2, unitTo='npc', valueOnly = TRUE)

xd <- (x2-x1)
yd <- (y2-y1)
v1 <- viewport(x=x1, y=y1, width=mywidth, height=myheight, angle=myangle, just=c("left", "bottom"))
pushViewport(v1)
################################################################################

brackets<- a_cb_brackets(phi=curvature, ticks=ticks, type=type)

grid.lines(brackets[1,], brackets[2,], gp=gpar(col=col, lwd=lwd, lty=lty))

popViewport()
}
