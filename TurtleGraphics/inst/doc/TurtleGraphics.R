## ----echo=FALSE,results='hide',warning=FALSE,message=FALSE,cache=FALSE----
options(digits=7)
options(width=73)
require('knitr')
# require('tikzDevice')
#
# options(tikzDefaultEngine = 'pdftex')
#
# options(tikzLatexPackages = c( # dolaczanie uzywanych pakietow TeX-a
#    '\\usepackage{amsmath,amssymb,amsfonts}', # pakiety AMS
#    '\\usepackage{tikz}',
# #   '\\usepackage[MeX,T1,plmath]{polski}', # obsluga m.in. polskich ogonkow
#    '\\usepackage[utf8]{inputenc}',
#    '\\usepackage[T1]{fontenc}',
#    '\\usetikzlibrary{calc}',
#    '\\usepackage[english]{babel}',
#    '\\selectlanguage{english}',
#    '\\usepackage{standalone}'
# ))
#
# options(tikzMetricsDictionary='~/R/tikzMetrics')
#
# options(tikzDocumentDeclaration = '\\documentclass[11pt]{standalone}\n')
#
# options(tikzMetricPackages = c(
#    '\\usepackage[utf8]{inputenc}',
#    '\\usepackage[T1]{fontenc}',
#    '\\usepackage{amsmath,amssymb,amsfonts}',
#    '\\usetikzlibrary{calc}',
#    '\\usepackage[english]{babel}',
#    '\\selectlanguage{english}'
# ))



# opts_knit$set(progress = TRUE, verbose = TRUE)

opts_chunk$set(
   keep.source=TRUE,
   out.width='4.5in',
   fig.width=6,
   fig.height=6/sqrt(3),
#    fig.path='figures-knitr/',
#    cache.path='cache-knitr/',
   cache=TRUE,
   tidy=FALSE,
#    dev='cairo_pdf',
#    dev.args=list(pointsize=11),
#    dev='tikz',
#    external=TRUE,
   fig.align='center',
   size='small'
)

# knit_theme$set(knit_theme$get('solarized-light'))

## ----results='hide', eval=FALSE----------------------------------------
#  install.packages("TurtleGraphics")

## ----results='hide',message=FALSE--------------------------------------
library("TurtleGraphics")

## ----results='hide', eval=FALSE----------------------------------------
#  turtle_init()

## ----fig.keep='last',fig1----------------------------------------------
turtle_init()
turtle_forward(dist=30)

## ----fig.keep='last',echo=3,fig2---------------------------------------
turtle_init()
turtle_forward(dist=30)
turtle_backward(dist=10)

## ----fig.keep='last',echo=-(1:3),fig3----------------------------------
turtle_init()
turtle_forward(dist=30)
turtle_backward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=10)
turtle_left(angle=135)
turtle_forward(dist=14)
turtle_left(angle=90)
turtle_forward(dist=14)
turtle_left(angle=135)
turtle_forward(dist=10)

## ----fig.keep='last',echo=-(1:11),fig4---------------------------------
turtle_init()
turtle_forward(dist=30)
turtle_backward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=10)
turtle_left(angle=135)
turtle_forward(dist=14)
turtle_left(angle=90)
turtle_forward(dist=14)
turtle_left(angle=135)
turtle_forward(dist=10)
turtle_up()
turtle_right(90)
turtle_forward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=17)
turtle_down()
turtle_left(angle=180)
turtle_forward(dist=34)

## ----fig.keep='last',fig5,echo=-(1:19)---------------------------------
turtle_init()
turtle_forward(dist=30)
turtle_backward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=10)
turtle_left(angle=135)
turtle_forward(dist=14)
turtle_left(angle=90)
turtle_forward(dist=14)
turtle_left(angle=135)
turtle_forward(dist=10)
turtle_up()
turtle_right(90)
turtle_forward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=17)
turtle_down()
turtle_left(angle=180)
turtle_forward(dist=34)
turtle_hide()
turtle_col(col="green")
turtle_left(angle=150)
turtle_forward(dist=20)
turtle_left(angle=60)
turtle_forward(dist=20)
turtle_show()

## ----fig.keep='last',fig6,echo=-(1:26)---------------------------------
turtle_init()
turtle_forward(dist=30)
turtle_backward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=10)
turtle_left(angle=135)
turtle_forward(dist=14)
turtle_left(angle=90)
turtle_forward(dist=14)
turtle_left(angle=135)
turtle_forward(dist=10)
turtle_up()
turtle_right(90)
turtle_forward(dist=10)
turtle_right(angle=90)
turtle_forward(dist=17)
turtle_down()
turtle_left(angle=180)
turtle_forward(dist=34)
turtle_hide()
turtle_col(col="green")
turtle_left(angle=150)
turtle_forward(dist=20)
turtle_left(angle=60)
turtle_forward(dist=20)
turtle_show()
turtle_left(angle=150)
turtle_lty(lty=4)
turtle_forward(dist=17)
turtle_lwd(lwd=3)
turtle_forward(dist=15)

## ----fig.keep='none'---------------------------------------------------
turtle_init()
turtle_status()

## ----fig.keep='none'---------------------------------------------------
turtle_init()
turtle_getpos()
turtle_getangle()

## ----fig.keep='last'---------------------------------------------------
turtle_init()
turtle_do({
   turtle_move(10)
   turtle_turn(45)
   turtle_move(15)
})

## ----fig.keep='last'---------------------------------------------------
turtle_init()
turtle_setpos(x=30, y=50)
turtle_do({
   for(i in 1:180) {
      turtle_forward(dist=1)
      turtle_right(angle=2)
   }
})

## ----fig.keep='last',echo=-1-------------------------------------------
set.seed(5671)
turtle_init()
turtle_do({
   for (i in 1:5) {
      x <- runif(1) # this function returns a random value between 0 and 1, see ?runif
      if (x>0.5) {
         turtle_right(angle=45)
         turtle_lwd(lwd=1)
         turtle_col(col="red")
      } else {
         turtle_left(angle=45)
         turtle_lwd(lwd=3)
         turtle_col(col="purple")
      }
      turtle_forward(dist=10)
   }
})

## ----results='hide',fun1-----------------------------------------------
turtle_square <- function(r) {
   for (i in 1:4) {
      turtle_forward(r)
      turtle_right(90)
   }
}

## ----fig.keep='last',dependson='fun1'----------------------------------
turtle_init()
turtle_square(10)
turtle_left(90)
turtle_forward(30)
turtle_square(5)

## ----fig.keep='last'---------------------------------------------------
set.seed(124) # assure reproducibility
turtle_init(100, 100, mode = "cycle")
turtle_do({
   for (i in 1:10) {
      turtle_left(runif(1, 0, 360))
      turtle_forward(runif(1, 0, 1000))
   }
})

## ----fig.keep='last'---------------------------------------------------
drawSpiral <- function(lineLen) {
   if (lineLen > 0) {
      turtle_forward(lineLen)
      turtle_right(50)
      drawSpiral(lineLen-5)
   }
   invisible(NULL) # return value: nothing interesting
}

turtle_init(500, 500, mode="clip")
turtle_do({
   turtle_setpos(x=0, y=0)
   turtle_col("blue")
   drawSpiral(500)
   turtle_setpos(x=250, y=0)
   turtle_left(45)
   turtle_col("green")
   drawSpiral(354)
   turtle_setangle(0)
})

## ----fig.keep='last'---------------------------------------------------
turtle_star <- function(intensity=1) {
   y <- sample(1:657, 360*intensity, replace=TRUE)
   for (i in 1:(360*intensity)) {
      turtle_right(90)
      turtle_col(colors()[y[i]])
      x <- sample(1:100,1)
      turtle_forward(x)
      turtle_up()
      turtle_backward(x)
      turtle_down()
      turtle_left(90)
      turtle_forward(1/intensity)
      turtle_left(1/intensity)
   }
}

set.seed(124)
turtle_init(500, 500)
turtle_do({
   turtle_left(90)
   turtle_star(5)
})

## ----fig.keep='last'---------------------------------------------------
turtle_brownian <- function(steps=100, length=10) {
   turtle_lwd(2)
   angles <- sample(c(90,270,180,0), steps,replace=TRUE)
   coll <- sample(1:657, steps, replace=TRUE)
   for (i in 1:steps){
      turtle_left(angles[i])
      turtle_col(colors()[coll[i]])
      turtle_forward(length)
   }
}

set.seed(124)
turtle_init(800, 800, mode="clip")
turtle_do(turtle_brownian(1000, length=25))

## ----fig.keep='last'---------------------------------------------------
fractal_tree <- function(s=100, n=2) {
   if (n <= 1) {
      turtle_forward(s)
      turtle_up()
      turtle_backward(s)
      turtle_down()
   }
   else {
      turtle_forward(s)

      a1 <- runif(1, 10, 60)
      turtle_left(a1)
      fractal_tree(s*runif(1, 0.25, 1), n-1)
      turtle_right(a1)

      a2 <- runif(1, 10, 60)
      turtle_right(a2)
      fractal_tree(s*runif(1, 0.25, 1), n-1)
      turtle_left(a2)

      turtle_up()
      turtle_backward(s)
      turtle_down()
   }
}

set.seed(123)
turtle_init(500, 600, "clip")
turtle_do({
   turtle_up()
   turtle_backward(250)
   turtle_down()
   turtle_col("darkgreen")
   fractal_tree(100, 12)
})

## ----fig.keep='last'---------------------------------------------------
koch <- function(s=50, n=6) {
   if (n <= 1)
      turtle_forward(s)
   else {
      koch(s/3, n-1)
      turtle_left(60)
      koch(s/3, n-1)
      turtle_right(120)
      koch(s/3, n-1)
      turtle_left(60)
      koch(s/3, n-1)
   }
}


turtle_init(600, 400, "error")
turtle_do({
   turtle_up()
   turtle_left(90)
   turtle_forward(250)
   turtle_right(180)
   turtle_down()
   koch(500, 6)
})

## ----fig.keep='last'---------------------------------------------------
drawTriangle <- function(points) {
   turtle_setpos(points[1,1], points[1,2])
   turtle_goto(points[2,1], points[2,2])
   turtle_goto(points[3,1], points[3,2])
   turtle_goto(points[1,1], points[1,2])
}

getMid <- function(p1, p2)
   (p1+p2)*0.5

sierpinski <- function(points, degree){
   drawTriangle(points)
   if (degree > 0) {
      p1 <- matrix(c(points[1,], getMid(points[1,], points[2,]),
                     getMid(points[1,], points[3,])), nrow=3, byrow=TRUE)

      sierpinski(p1, degree-1)
      p2 <- matrix(c(points[2,], getMid(points[1,], points[2,]),
                     getMid(points[2,], points[3,])), nrow=3, byrow=TRUE)

      sierpinski(p2, degree-1)
      p3 <- matrix(c(points[3,], getMid(points[3,], points[2,]),
                     getMid(points[1,], points[3,])), nrow=3, byrow=TRUE)
      sierpinski(p3, degree-1)
   }
   invisible(NULL)
}

turtle_init(520, 500, "clip")
turtle_do({
   p <- matrix(c(10, 10, 510, 10, 250, 448), nrow=3, byrow=TRUE)
   turtle_col("red")
   sierpinski(p, 6)
   turtle_setpos(250, 448)
})

