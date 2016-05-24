opar <- par(ask = dev.interactive(orNone = TRUE))

###############################################
###############################################
###
### loads the data

library(adehabitat)
data(puechabon)
xy<-puechabon$locs[,c("X","Y")]
id<-puechabon$locs$Name

## The data are:
xy[1:4,]     ## relocations coordinates
id[1:4]      ## ID


###############################################
###############################################
###
### Home ranges

## MCP
hr<-mcp(xy, id)         ## home range estimation
plot(hr)                ## displays the MCP
(jj<-mcp.area(xy, id))  ## home range size
plot(jj)                ## plots home-range size


## Kernel home range
(hr<-kernelUD(xy, id, h="LSCV", grid=100)) ## UD estimation
plotLSCV(hr)                               ## LSCV criterion
image(hr)                                  ## displays the UD
(jj<-kernel.area(xy, id))                  ## home range size
plot(jj)                                   ## Plots home range size
ver <- getverticeshr(hr)                   ## home-range contours
plot(ver)                                  ## Plots contours


## Ellipse home range
if (require(ellipse)) {

    foo<-function(x) {
        u<-ellipse(cov(x), centre=apply(x,2,mean))
        plot(u, ty="n", axes=FALSE, asp=1)
        points(x, pch=16)
        lines(u)
        box()
    }
    opar <- par(mfrow=c(2,2), mar=c(0,0,0,0))
    lapply(split(xy,id), foo)
    par(opar)
}


## Cluster home range
(res <- clusthr(xy, id))       ## Computes the home range
plot(res)                      ## Displays the home range
clusthr.area(res)              ## Computes the home-range size
ver <- getverticesclusthr(res) ## Computes the home-range contours
plot(ver)                      ## Plots the contours



## Nearest neighbour convex hull
nn <- NNCH(xy,id, k=c(5:7))      ## Home range estimation
summary(nn)                      ## Summary of the analysis
plot(nn, k=5, border=1)          ## Plot the home range
uu <- NNCH.select(nn, id="Chou") ## Select one animal
## Rasterize the home range
asc <- ascgen(xy,nrcol=100)
asc <- NNCH.asciigrid(uu, k=7, asc=asc)
image(asc)
points(xy[id=="Chou",], pch=16, col="red")



cat("******************************************************\n",
    "The deeply commented source for this demo can be found in the file:\n",
    file.path(system.file(package = "adehabitat"), "demo", "homerange.r"),
    "\n******************************************************\n")
