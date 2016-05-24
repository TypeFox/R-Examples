fig_PharmPow_crossover <- function(data,
                              power=80,
                              colorabove="black",
                              colorbelow="red",
                              designAlab="Number of patients sampling design A",
                              designBlab="Number of patients sampling design B",
                              zaxeslab="power (%)",
                              axessize=1.0,
                              labsize=1.0){

### Prepare data and create figure ###

data <- read.csv(data,header=FALSE)

coll <- data[2:ncol(data),1]
rows <- data[1,2:ncol(data)]

d <- data[2:ncol(data),2:ncol(data)]

colnames(d) <- coll
rownames(d) <- rows
z <- as.matrix(d)

### Number of ID's 3D plot
xx <- rep(0:(ncol(z)-1),each=nrow(z))
yy <- rep(0:(nrow(z)-1),ncol(z))
zz <- as.vector(z)
clc <- rep(NA,length(zz))
clc[zz<power] <- colorbelow
clc[zz>=power] <- colorabove
cl <- rep(1,length(zz))

s3d <- scatterplot3d(xx,yy,z,
                     type="p",
                     color=rep(clc,cl),
                     highlight.3d=FALSE,
                     angle=55,
                     cex.axis=axessize,
                     cex.lab=labsize,
                     zlab=zaxeslab,
                     xlab=designAlab,
                     ylab=designBlab,
                     zlim=c(0,100)
                     )

fit <- c(80,0,0)
s3d$plane3d(fit,col="black",lwd=2)

}
