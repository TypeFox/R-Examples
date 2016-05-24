see.mixedII <- function(){

if(dev.capabilities("locator")$locator == FALSE) stop("Device cannot implement function")

#--------------------- SETUP -------------------#

setup.mixed <- function(){
if(as.numeric(dev.cur())!=1)dev.off()
dev.new(width=9, height = 4)

make.table <- function(nr, nc) {
     savepar <- par(mar=rep(0, 4),pty = "m")
     plot(c(0, nc + 6.9), c(1, -(nr + 1)),
     type="n", xlab="", ylab="", axes=FALSE)
     savepar
      }
 
draw.cell<-function(text, r, c, cex = .9){
     rect(c, -r, c+1, -r+1) 
     text(c+.5, -r+.5, text, cex = cex)
      }

make.table(4,19)

cn <- c("",paste("R",1:18,sep=""))
c1 <- c(7,6,5,7,6,5,4,4,4,1,2,3,4,4,4,1,2,3)
c2 <- c(4,4,4,1,2,3,7,6,5,7,6,5,1,2,3,4,4,4)
c3 <- c(1,2,3,4,4,4,1,2,3,4,4,4,7,6,5,7,6,5)
data <- t(cbind(c1,c2,c3))

for(i in 1:19){
    draw.cell(cn[i], 1, i)
    }
    draw.cell("F1",2,1)

for(i in 2:19){
    draw.cell(c1[i-1], 2, i)
    }
    draw.cell("F2",3,1)

for(i in 2:19){
    draw.cell(c2[i-1], 3, i)
    }
    draw.cell("F3",4,1)

for(i in 2:19){
    draw.cell(c3[i-1], 4, i)
    }

text(11, 0.6, "Random", cex=1.2, font = 2)
text(-0.2, -2, "Fixed", cex=1.2, font = 2, srt = 90)
text(23.25, 0, "Mean", cex=1.1, font = 2)
text(21.5, -0.5, expression(underline("All levels")))
text(25, -0.5, expression(underline("Selected levels")))

mn <- apply(data, 1, mean)
for(i in 2:4)text(21.5,-i+.5,round(mn[i-1],0))

rect(8.2, -5.1, 11.8, -4.2) 
text(10, -4.65, "SAMPLE", font = 2, col = 1)

fl <- function(){
ans <- locator(1)
    yv <- ans$y
    xv <- ans$x
    if((yv>=-5.1)&(yv<=-4.2)&(xv >= 8.2)&(xv <= 11.8))
        {rect(8.2, -5.1, 11.8, -4.2, col = rgb(red=0.5,blue=0.5,green=0.5,alpha=.6));sample.mixed()}
    else fl1()
        }

fl1 <- function(){
ans <- locator(1)
    yv <- ans$y
    xv <- ans$x
    if((yv>=-5.1)&(yv<=-4.2)&(xv >= 8.2)&(xv <= 11.8))
        {rect(8.2, -5.1, 11.8, -4.2, col = rgb(red=0.5,blue=0.5,green=0.5,alpha=.6));sample.mixed()}
    else fl() 
    }   # infinite loop    

fl()
}

#--------------- SAMPLE -----------------#

sample.mixed <- function(r = 4, c = 19, n = 3){
col.row<- function(col=rgb(red=0.5, blue=0.5, green=0.5, alpha=.6), r, c) 
    {
        rect(c, -r-3, c+1, -r+1, col = col)
    }

sn <- sample(seq(2, 19), size = n)

for(i in 1:length(sn)) col.row(c = sn[i], r=1)
   
rect(8.2, -5.1, 11.8, -4.2) 
text(10, -4.65, "SAMPLE", font = 2, col = 1)

c1 <- c(7,6,5,7,6,5,4,4,4,1,2,3,4,4,4,1,2,3)
c2 <- c(4,4,4,1,2,3,7,6,5,7,6,5,1,2,3,4,4,4)
c3 <- c(1,2,3,4,4,4,1,2,3,4,4,4,7,6,5,7,6,5)
data <- t(cbind(c1,c2,c3))

mn <- apply(data[,sn-1], 1, mean)

for(i in 2:4)text(25,-i+.5,round(mn[i-1],1))

fl <- function(){
ans <- locator(1)
    yv <- ans$y
    xv <- ans$x
        if((yv>=-5.1)&(yv<=-4.2)&(xv >= 8.2)&(xv <= 11.8)) 
        {dev.off(); setup.mixed()}
        else fl1()
        }

fl1 <- function(){
ans <- locator(1)
    yv <- ans$y
    xv <- ans$x
        if((yv>=-5.1)&(yv<=-4.2)&(xv >= 8.2)&(xv <= 11.8)) 
        {dev.off(); setup.mixed()}
        else fl()
        }
fl()
}
setup.mixed()
}






  

