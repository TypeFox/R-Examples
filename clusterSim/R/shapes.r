.numObjects<-function(numObjects,numClusters){
    if(length(numObjects)<numClusters){
        numObjects<-rep(numObjects,numClusters)
     }
    if(length(numObjects)>numClusters){
        numObjects<-numObjects[1:numClusters]
     }
     numObjects
}

.toCsv<-function(file,data,klasy,outputRowNames,outputColNames,csv2){
  if(paste(file,"",sep="")!=""){
    if(csv2){
      write.table(cbind(1:dim(data)[1],klasy,data),file=file,sep=";",dec=",",row.names=outputRowNames,col.names=outputColNames)
    }
    else{
      write.table(cbind(1:dim(data)[1],klasy,data),file=file,row.names=outputRowNames,col.names=outputColNames)
    }
  }
}

shapes.worms<-function(numObjects=180,shape1x1=-2,shape1x2=2,shape2x1=-0.5,shape2x2=2.5,shape2a=1.5,shape2b=5.5,tol=0.1,outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE){
    f1<-function(x){
    x^2
    }


    f2<-function(x,shape2a,shape2b){
    -(x-shape2a)^2+shape2b
    }

    worm<-function(lo,shape1x1,shape1x2,shape2x1,shape2x2,shape2a,shape2b,tol){
      data<-array(0,c(sum(lo),2))
      lim<-seq(shape1x1,shape1x2,length.out=lo[1])
      for(i in 1:lo[1]){
        data[i,1]<-lim[i]+rnorm(1,0,tol)
        data[i,2]<-f1(lim[i])+rnorm(1,0,tol)
      }
      lim<-seq(shape2x1,shape2x2,length.out=lo[2])
      for(i in 1:lo[2]){
        data[i+lo[1],1]<-lim[i]+rnorm(1,0,tol)
        data[i+lo[1],2]<-f2(lim[i],shape2a,shape2b)+rnorm(1,0,tol)
      }
      data
    }
    lo<-.numObjects(numObjects,2)
    klasy<-c(rep(1,lo[1]),rep(2,lo[2]))
    data<-worm(lo,shape1x1,shape1x2,shape2x1,shape2x2,shape2a,shape2b,tol)
    .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
    .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
    list(data=data,clusters=klasy)
}



shapes.circles2<-function(numObjects=180, shape1rFrom=0.75,shape1rTo=0.9,shape2rFrom=0.35,shape2rTo=0.5,outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE){
     lo<-.numObjects(numObjects,2)
     t1 <- 2 * pi * runif(lo[1])
     t2 <- 2 * pi * runif(lo[2])
     y2<-x2<-y1<-x1<-NULL
     for(t in t1) x1 <- c(x1,cos(t)*runif(1,shape1rFrom,shape1rTo))
     for(t in t1) y1 <- c(y1,sin(t)*runif(1,shape1rFrom,shape1rTo))
     X1 <- t(as.matrix(rbind(x1, y1)))
     for(t in t2) x2 <- c(x2,cos(t)*runif(1,shape2rFrom,shape2rTo))
     for(t in t2) y2 <- c(y2,sin(t)*runif(1,shape2rFrom,shape2rTo))
     X2 <- t(as.matrix(rbind(x2, y2)))
     data <- as.matrix(rbind(X1, X2))
     klasy<-c(rep(1,lo[1]),rep(2,lo[2]))
    .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
    .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
    list(data=data,clusters=klasy)
}


shapes.circles3<-function(numObjects=180,shape1rFrom=0.15,shape1rTo=0.3,shape2rFrom=0.55,shape2rTo=0.7,shape3rFrom=1.15,shape3rTo=1.3,outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE){

     lo<-.numObjects(numObjects,3)
     t1 <- 2 * pi * runif(lo[1])
     t2 <- 2 * pi * runif(lo[2])
     t3 <- 2 * pi * runif(lo[3])
       y3<-x3<-y2<-x2<-y1<-x1<-NULL
     for(t in t1) x1 <- c(x1,cos(t)*runif(1,shape1rFrom,shape1rTo))
     for(t in t1) y1 <- c(y1,sin(t)*runif(1,shape1rFrom,shape1rTo))
         X1 <- t(as.matrix(rbind(x1, y1)))
     for(t in t2) x2 <- c(x2,cos(t)*runif(1,shape2rFrom,shape2rTo))
     for(t in t2) y2 <- c(y2,sin(t)*runif(1,shape2rFrom,shape2rTo))
         X2 <- t(as.matrix(rbind(x2, y2)))
     for(t in t3) x3 <- c(x3,cos(t)*runif(1,shape3rFrom,shape3rTo))
     for(t in t3) y3 <- c(y3,sin(t)*runif(1,shape3rFrom,shape3rTo))
         X3 <- t(as.matrix(rbind(x3, y3)))
     data <- as.matrix(rbind(X1,X2,X3))
     klasy<-c(rep(1,lo[1]),rep(2,lo[2]),rep(3,lo[3]))
    .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
    .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
    list(data=data,clusters=klasy)
}


shapes.bulls.eye<-function(numObjects=180, shape1rFrom=0.75,shape1rTo=0.95,shape2rTo=0.45,outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE){
    shapes.circles2(numObjects, shape1rFrom,shape1rTo,shape2rFrom=0,shape2rTo,outputCsv, outputCsv2, outputColNames,outputRowNames)
}

shapes.two.moon<-function(numObjects=180,shape1a=-0.4,shape2b=1,shape1rFrom=0.8, shape1rTo=1.2,shape2rFrom=0.8, shape2rTo=1.2, outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE){

            lo<-.numObjects(numObjects,2)
            x <- matrix(0, nrow=sum(lo), ncol=2)
            
            for(i in 1:sum(lo)){
              alpha<-runif(1,0,2*pi)
              if(i>lo[1]){
                r=runif(1,shape2rFrom,shape2rTo)
              }
              else{
                r=runif(1,shape1rFrom,shape1rTo)
              }
              x[i,1]<-r*cos(alpha)
              x[i,2]<-r*sin(alpha)
              if(i<=lo[1]){
                x[i,1]=shape1a+abs(x[i,1])
              }
              else{
                x[i,1]=-abs(x[i,1])
                x[i,2]=x[i,2]-shape2b
              }
           }
          data<-x
          klasy<-c(rep(1,lo[1]),rep(2,lo[2]))
          .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
          .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
          list(data=data,clusters=klasy)

}

shapes.blocks3d<-function(numObjects=180,shapesUnitSize=0.5, shape2coordinateX=1.2,shape2coordinateY=1.2,shape2coordinateZ=1.2, outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE){
          lo<-.numObjects(numObjects,2)
          x <- matrix(0, nrow=sum(lo), ncol=3)
          for(i in 1:sum(lo)){
          t<-sample(1:4,1)
            if(t==1){
              x[i,1]<-runif(1,0,shapesUnitSize)
              x[i,2]<-runif(1,0,shapesUnitSize)
              x[i,3]<-runif(1,0,shapesUnitSize)
            }
            if(t==2){
              x[i,1]<-runif(1,0,shapesUnitSize)
              x[i,2]<-runif(1,0,shapesUnitSize)
              x[i,3]<-runif(1,shapesUnitSize,shapesUnitSize*2)
            }
            if(t==3){
              x[i,1]<-runif(1,0,shapesUnitSize)
              x[i,2]<-runif(1,shapesUnitSize,shapesUnitSize*2)
              x[i,3]<-runif(1,0,shapesUnitSize)
            }
            if(t==4){
              x[i,1]<-runif(1,shapesUnitSize,shapesUnitSize*2)
              x[i,2]<-runif(1,0,shapesUnitSize)
              x[i,3]<-runif(1,0,shapesUnitSize)
            }
          if(i>lo[1]){
            x[i,]<-c(shape2coordinateX,shape2coordinateY,shape2coordinateZ)-x[i,]
          }
          }
          data<-x
          klasy<-c(rep(1,lo[1]),rep(2,lo[2]))
          .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
          .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
          list(data=data,clusters=klasy)

}


#sw<-shapes.worms(180)
#plot(sw$data,col=rainbow(2)[sw$clusters])

#sw<-shapes.circles2(180)
#plot(sw$data,col=rainbow(2)[sw$clusters])

#sw<-shapes.circles3 (180)
#plot(sw$data,col=rainbow(3)[sw$clusters])

#sw<-shapes.bulls.eye (180)
#plot(sw$data,col=rainbow(2)[sw$clusters])

#sw<-shapes.two.moon (180)
#plot(sw$data,col=rainbow(2)[sw$clusters])

#library(rgl)
#sw<-shapes.blocks3d (300,1,3,3,3)
#plot3d(sw$data,col=rainbow(2)[sw$clusters])