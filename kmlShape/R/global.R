.onLoad <- function(lib,pkg){library.dynam("kmlShape",pkg,lib)}

### A refaire en mieux

#shapeCriterion <- function(Py,Qy,times){
#   return(c(dist(rbind(Py,Qy)/length(times)),
#      distFrechet(Px=times,Py=Py,Qx=times,Qy=Qy,timeScale=1),
#      distFrechet(Px=times,Py=Py,Qx=times,Qy=Qy,timeScale=0.1),
#      dtw(Py,Qy)$distance)
#   )
#}



printLineShort <- function(x){
    if (length(x) <= 10){
        print(x)
    }
    else {
        print(c(as.character(x[1:10]),"..."),quote=FALSE)
    }
}

