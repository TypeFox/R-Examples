##
## Computing and ploting a bpca object with 'rgl' package - 3d
##

open3d()
plot(pca <- bpca(iris[-5],
                 d=1:3),
     rgl.use=TRUE,
     var.col='brown',
     var.factor=.4,
     var.cex=1.2,
     obj.names=FALSE,
     obj.cex=.8,
     obj.col=c('blue', 'green', 'red')[unclass(iris$Species)],
     simple.axes=TRUE)

scores <- pca$coord$objects

ell <- ellipse3d(cov(scores),
                 center=mean(scores),
                 level=0.68)

plot3d(ell,
       col='gray',
       alpha=0.2,
       add=TRUE)

play3d(rgl:::spin3d(axis=c(1,2,3)),
       duration=12)

# This graphic was suggested by Michael Friendly (York University). 
# Suggestion: Interact with the graphic with the mouse
# left button: press, maintain and movement it to interactive rotation;
# right button: press, maintain and movement it to interactive zoom.
# Enjoy it!

