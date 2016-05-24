plot.standz.all <-
function(x, which=NULL,...){
if(!is.null(which)) {
        if(which=="w") cual <- x$which.w else cual <- x$which.l
width<-t(x$data.stdz[cual,])} else width<-t(x$data.stdz)
matplot(width,...)
}

