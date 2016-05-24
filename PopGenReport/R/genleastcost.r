#leatcost function
genleastcost <- function(cats, fric.raster, gen.distance, NN=NULL, pathtype="leastcost", plotpath=TRUE, theta=1)
{
### helper mmod function until mmod is updated

pairwise_D2 <- function (x, linearized = FALSE) 
{
      pops <- seppop(x)
      n.pops <- length(pops)
      allP <- utils::combn(1:n.pops, 2)
      pair <- function(index.a, index.b) {
        a <- pops[[index.a]]
        b <- pops[[index.b]]
        temp <- repool(a, b)
        return(D_Jost(temp)$global.het)
      }
      res <- sapply(1:dim(allP)[2], function(i) pair(allP[, i][1], 
                                                     allP[, i][2]))
      attributes(res) <- list(class = "dist", Diag = FALSE, Upper = FALSE, 
                              Labels = popNames(x), Size = n.pops)
      if (linearized) {
        return(res/(1 - res))
      }
      else (return(res))
}
  
  pairwise_Gst_Hedrick2 <-
    function (x, linearized = FALSE) 
    {
      pops <- seppop(x)
      n.pops <- length(pops)
      allP <- utils::combn(1:n.pops, 2)
      pair <- function(index.a, index.b) {
        a <- pops[[index.a]]
        b <- pops[[index.b]]
        temp <- repool(a, b)
        return(Gst_Hedrick(temp)$global)
      }
      res <- sapply(1:dim(allP)[2], function(i) pair(allP[, i][1], 
                                                     allP[, i][2]))
      attributes(res) <- list(class = "dist", Diag = FALSE, Upper = FALSE, 
                              Labels = popNames(x), Size = n.pops)
      if (linearized) {
        return(res/(1 - res))
      }
      return(res)
    }
  
  pairwise_Gst_Nei2 <-
    function (x, linearized = FALSE) 
    {
      pops <- seppop(x)
      n.pops <- length(pops)
      allP <- utils::combn(1:n.pops, 2)
      pair <- function(index.a, index.b) {
        a <- pops[[index.a]]
        b <- pops[[index.b]]
        temp <- repool(a, b)
        return(Gst_Nei(temp)$global)
      }
      res <- sapply(1:dim(allP)[2], function(i) pair(allP[, i][1], 
                                                     allP[, i][2]))
      attributes(res) <- list(class = "dist", Diag = FALSE, Upper = FALSE, 
                              Labels = popNames(x), Size = n.pops)
      if (linearized) {
        return(res/(1 - res))
      }
      return(res)
    }
  

  
  
  
  if (is.null(NN) & pathtype=="leastcost") 
  {
    cat("NN is not specified!\nPlease specify the number of nearest neighbour to use for the least-cost path calculations (NN=4 or NN=8). If linear features are tested you may want to consider NN=4 otherwise NN=8 is the most commonly used and prefered option. In any case check the actual least-cost paths for artefacts by inspecting the plot on least-cost paths.\n")
    return()
  }
  
dist.type<-NA
if (gen.distance=="D" || gen.distance=="Gst.Hedrick" || gen.distance=="Gst.Nei") dist.type<- "pop" 

if (gen.distance=="Kosman" || gen.distance=="Smouse" || gen.distance=="propShared") dist.type<- "ind" 

if (is.na(dist.type)) 
  {cat("No valid genetic distance type was provided. Please check ?landgenreport for valid options\n")
   return(-1)
}

if (dist.type=="pop")
{
#calculate the centers if population meassurment is wanted
c.x <- tapply(cats@other$xy[,1],cats@pop, mean)
c.y <- tapply(cats@other$xy[,2],cats@pop, mean)
cp<-cbind(c.x, c.y)
} else 
{
cp <- cbind(cats@other$xy[,1], cats@other$xy[,2])
}


eucl.mat <- round(as.matrix(dist(cp)),3)

if (dist.type=="pop") 
{

  dimnames(eucl.mat) <- list(popNames(cats), popNames(cats)) 
  npop <- length(levels(cats@pop))
} else
  
{
  dimnames(eucl.mat) <- list(indNames(cats), indNames(cats)) 
  npop <- length(indNames(cats))
}





#check if fric.raster is a stack or not...
mats <- list()
mats.names<- NA
mats.pathlength<- list()
mats.paths<- list()

pathlength.mat <- NULL
paths <- NULL



n.mats <- dim(fric.raster)[3] #number of rasters in the stack



for (ci in 1:n.mats)
{

plot(fric.raster[[ci]], main=paste(names(fric.raster)[ci],":",pathtype,", NN=",NN,sep=""))
 #image(fric.raster, col=fric.raster@legend@colortable, asp=1)

 mapcolor <-   col2rgb(cats@other$mapdotcolor)/255
 if (is.null(mapcolor)) colX= rgb(0,0,1,0.8) else 
 colX = rgb(mapcolor[1,],mapcolor[2,], mapcolor[3,],cats@other$mapdotalpha)
points(cats@other$xy,cex=1, pch=16, col="blue")
# points(cats@other$xy,cex=cats@other$mapdotsize, pch= cats@other$mapdottype, col=colX)
if (dist.type=="pop")  points(cp,cex=cats@other$mapdotsize*1.5 , pch= 16, col="black")


#create friction matrix
fric.mat <- transition(fric.raster[[ci]],function(x) 1/x[2],NN)

# fric.mat <- transition(fric.raster,function(x) 1/(abs(x[1]-x[2])),8)
#set distances to meters  if no projected already
fric.mat@crs@projargs<- "+proj=merc +units=m"
fric.mat.cor <- geoCorrection(fric.mat)


if (pathtype=="leastcost") cd.mat <-costDistance(fric.mat.cor, cp, cp)

if (pathtype=="rSPDistance") cd.mat <- rSPDistance(fric.mat.cor, cp, cp, theta=1)

if (pathtype=="commute") cd.mat <-as.matrix(commuteDistance(fric.mat.cor, cp))


  dimnames(cd.mat) <- dimnames(eucl.mat) 




if (pathtype=="leastcost" & plotpath==TRUE)   #only show paths if leastcost otherwise not possible
{
comb <- t(combn(1:npop,2))


#pathlength matrix
pathlength.mat <- cd.mat
pathlength.mat[,] <- 0
paths<- list()

cols <- rainbow(dim(comb)[1], alpha=0.5)
for (i in 1:dim(comb)[1])
{

if (dist(rbind(cp[comb[i,1],], cp[comb[i,2],]))==0)
{
 ll <- Line(rbind( cp[comb[i,1],], cp[comb[i,2],]))
 S1 <- Lines(list(ll),ID="Null")
 sPath <- SpatialLines(list(S1))
} else 
{
sPath <- shortestPath(fric.mat.cor, cp[comb[i,1],], cp[comb[i,2],], output="SpatialLines")
}

lines(sPath, lwd=1.5, col=cols[i])
paths[[i]] <- sPath
ll <-  round(SpatialLinesLengths(sPath),3)
pathlength.mat[comb[i,1],comb[i,2]] <- ll
pathlength.mat[comb[i,2],comb[i,1]] <- ll
}

}

mats[[ci]] <- cd.mat
mats.names[[ci]] <- names(fric.raster)[ci]
mats.pathlength[[ci]] <- pathlength.mat
mats.paths[[ci]] <- paths



} #end of ci loop

#mats[[n.mats+1]] <- eucl.mat
#mats.names[n.mats+1]<- "Euclidean"
#
names(mats)  <- names(fric.raster)
#put other calculations here....
# Calculate genetic distances across subpopulations

if (gen.distance=="Gst.Nei")
{
gendist.mat<-as.matrix(pairwise_Gst_Nei2(cats))
}
if (gen.distance=="Gst.Hedrick")
{
gendist.mat<-as.matrix(pairwise_Gst_Hedrick2(cats))
}
if (gen.distance=="D")
{
gendist.mat<-round(as.matrix(pairwise_D2(cats)),4)
}

if (gen.distance=="Smouse")
{
gendist.mat <- as.matrix(gd.smouse(cats,verbose=FALSE))
}
if (gen.distance=="Kosman")
{
gendist.mat <-as.matrix(as.dist(gd.kosman(cats)$geneticdist))
}
if (gen.distance=="propShared")
{
  gendist.mat <-as.matrix(as.dist(propShared(cats)))
}

  dimnames(gendist.mat)<-dimnames(eucl.mat)


return(list( gen.mat=gendist.mat, eucl.mat=eucl.mat,cost.matnames=mats.names, cost.mats=mats,pathlength.mats= mats.pathlength,  paths=mats.paths))
}


