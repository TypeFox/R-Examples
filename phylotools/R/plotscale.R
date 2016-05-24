#### Function plotscale as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


plotscale <-
function( inputdata, len = NULL, wid = NULL, scale = NULL )
{
     if (any(c(is.null(inputdata$tag), is.null(inputdata$gx),is.null(inputdata$gy))))
         stop( "The input data must include one colomn named \"tag\", \"gx\" and \"gy\".")
     
     if (any(c(is.null(len), is.null(wid))))
         stop("You must specify the length and width of the FDP.")
     
     if (is.null(scale))
         stop("Oops! You must specify the scale.")
     
     ######Calc x
     inputdatax<-inputdata[order(inputdata$gx),]
     inputdatax$gx[inputdatax$gx == len] <- len-0.001;
     x1<-seq(0,len, by = scale)
     nox<-NULL
     for (i in 1:(length(x1))){
         ggg <- subset(inputdatax, (inputdatax$gx >= x1[i])&(inputdatax$gx < x1[i+1])); 
         nggg <- nrow(ggg)
         if (nggg==0){
             nox = nox
         }else{
             nox<-append(nox, nggg)
         }
     }
     grid.nox <-rep(1:length(nox), nox)
     resultx<-data.frame(inputdatax, grid.nox)
     
     ########Calc y
     inputdatay<-inputdata[order(inputdata$gy), ]
     inputdatay$gy[inputdatay$gy == wid] <- wid-0.001
     scale = scale
     y1 <- seq(0, wid, by = scale)
     noy <- NULL
     for (i in 1:(length(y1))){
         ggg <- subset(inputdatay, (inputdatay$gy >= y1[i])&(inputdatay$gy < y1[i+1])) 
           nggg<-nrow(ggg)
          if (nggg == 0){
            noy = noy
         }else{
            noy <- append(noy, nggg)
         }
     }
     grid.noy <- rep(1:length(noy),noy)
     resulty <- data.frame(inputdatay, grid.noy)
     ordered.x <- resultx[order(resultx$tag), ]
     ordered.y <- resulty[order(resulty$tag), ]
     result <- data.frame(ordered.x, grid.noy=ordered.y$grid.noy, gridnames=paste("X", ordered.x$grid.nox,"Y", ordered.y$grid.noy, sep=""))
     return(result)
}

