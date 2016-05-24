## ----fig.width=7, fig.retina=2, fig.align='left', fig.cap="Rectangular Cartogram of the U.S. election 2004; The area corresponds to the number of electors (color indicates the party red: democrates / blue: republican; the color intensity ~ outcome of the vote.). The graphic was computed by using the original implementation of the construction heuristic RecMap MP2 introduced in [@recmap].", echo=FALSE, warning=FALSE, comment="ccc", error=FALSE, message=FALSE----

library(recmap)
op<-par(mar=c(0,0,0,0), bg='black')
recmap:::.draw_recmap_us_state_ev()
par(op)

## ------------------------------------------------------------------------
usa <- data.frame(x=state.center$x, 
    y = state.center$y, 
    # make the rectangles overlapping by correcting lines of longitude distance
    dx = sqrt(state.area) / 2 / (0.8 * 60 * cos(state.center$y*pi/180)), 
    dy = sqrt(state.area) / 2 / (0.8 * 60) , 
    z = sqrt(state.area),
    name = state.name)

## ----fig.width=7, fig.height=3-------------------------------------------
op<-par(mfrow=c(1,1), mar=c(0,0,0,0))
library(recmap)
plot_recmap(M <- usa[!usa$name %in% c("Hawaii", "Alaska"), ],  col.text = 'black', lwd=2)

## ----echo=FALSE, fig.width=2, fig.height=2, fig.align='center'-----------

draw_and_place_rectangle <- function(alpha=0.0, a.x=2, a.y=-5, a.dx = 20, a.dy = 5, b.dx=1.5, b.dy=5, ...){
   
   rect(a.x - a.dx - b.dx, a.y - a.dy - b.dy, 
        a.x + a.dx + b.dx, a.y + a.dy + b.dy, 
        lwd=3, border='grey')
  
   rect(a.x - a.dx, a.y - a.dy, 
        a.x + a.dx, a.y + a.dy)
   
  
  c <- recmap:::place_rectanle(a.x, a.y, a.dx, a.dy, b.dx, b.dy, alpha)

   rect(c$x - c$dx, c$y - c$dy, c$x + c$dx, c$y + c$dy, ...)
}

op <- par(mar=c(0,0,0,0))
plot(0,0 , xlim=c(-25,25), ylim=c(-25,25), asp=1, xlab='', ylab='', axes=FALSE); abline(v=0,h=0)
r<-lapply(seq(0, pi/2, length=90), function(alpha){draw_and_place_rectangle(alpha, col='#DD001199', border='#88888877')})
r<-lapply(seq(pi/2, pi, length=90), function(alpha){draw_and_place_rectangle(alpha, col='#00DD1199', border='#88888877')})
r<-lapply(seq(-pi, -pi/2, length=90), function(alpha){draw_and_place_rectangle(alpha, col='#11DDDD99', border='#88888877')})
r<-lapply(seq(-pi/2, 0, length=90), function(alpha){draw_and_place_rectangle(alpha, col='#0011DD99', border='#88888877')})

par(op)

## ------------------------------------------------------------------------
head(Cartogram <- recmap(Map <- usa[!usa$name %in% c("Hawaii", "Alaska"), ]))

## ----fig.width=8, fig.height=4, fig.align='left', fig.cap="Area ~ population estimate as of July 1, 1975;"----

op<-par(mfrow=c(1,1), mar=c(0,0,0,0))
usa$z <- state.x77[, 'Population']
M <- usa[!usa$name %in% c("Hawaii", "Alaska"), ]
plot_recmap(Cartogram.Population <- recmap(M[order(M$x),]), col.text = 'black', lwd=2)

## ----fig.width=8, fig.height=4, fig.align='left', fig.cap="Area ~ population estimate as of July 1, 1975; a better index order has been choosen to minimize the relative position error."----
op<-par(mfrow=c(1,1), mar=c(0,0,0,0))
# index order
smp <- c(20,47,4,40,9,6,32,33,3,10,34,22,2,28,15,12,39,7,42,45,19,13,43,30,24,
         25,11,17,37,41,26,29,21,35,8,36,14,16,31,48,46,38,23,18,1,5,44,27)
plot_recmap(Cartogram.Population <- recmap(M[smp,]), col.text = 'black', lwd=2)

## ----fig.width=8, fig.height=4, fig.align='left', fig.cap="Area ~ capita income (1974);"----
op<-par(mfrow=c(1,1), mar=c(0,0,0,0))
usa$z <- state.x77[, 'Income']
M <- usa[!usa$name %in% c("Hawaii", "Alaska"), ]
plot_recmap(Cartogram.Income <- recmap(M[order(M$x),]), col.text = 'black', lwd=2)

## ----fig.width=8, fig.height=4, fig.align='left', fig.cap="Area ~ mean number of days with minimum temperature below freezing (1931–1960) in capital or large city;"----
op<-par(mfrow=c(1,1), mar=c(0,0,0,0))
usa$z <- state.x77[, 'Frost'] 
M <- usa[!usa$name %in% c("Hawaii", "Alaska"), ]
plot_recmap(Cartogram.Income <- recmap(M[order(M$x),]), col.text = 'black', lwd=2)

## ----fig.width=7, fig.height=3.5, fig.align='center', fig.retina=2, fig.cap="chess board fun"----
op<-par(mar=c(0,0,0,0), mfrow=c(1, 2), bg='white')

plot_recmap(ChessBoard <- recmap:::.checker_board(8),
            col=c('white','white','white','black')[ChessBoard$z])

# found by a Non-deterministic Turing machine
smp <- c(8,56,18,5,13,57,3,37,62,58,7,16,40,59,17,34,29,41,46,27,54,43,2,21,
         38,52,31,20,28,48,1,22,55,11,25,19,50,10,24,53,47,30,45,44,32,35,51,
         15,64,12,14,39,26,6,42,33,4,36,63,49,60,61,9,23)

plot_recmap(Cartogram.ChessBoard <- recmap(ChessBoard[smp,]), 
            col=c('white','white','white','black')[Cartogram.ChessBoard$z])

