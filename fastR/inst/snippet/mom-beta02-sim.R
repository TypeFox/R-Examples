plot3 <- xyplot(shape2~shape1, results, panel=function(x,y,...){
            panel.abline(a=0,b=5/2)
            panel.xyplot(x,y,...)
            })
plot4 <- xhistogram(~shape2/shape1, results, type='density', v=2.5)
