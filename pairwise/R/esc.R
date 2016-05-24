#' @title Expected Score Curves Plots
#' @export esc
#' @description plotting function for plotting expected score curves.
#' @details no details in the moment.
#' @param pers_obj an object of class \code{"pers"} as a result from function \code{\link{pers}}.
#' @param itemnumber an integer, defining the number of the item to plot the respective categoy probability for. This is set to an arbitrary default value of \code{itemnumber = 1} to avoid error messages when you forget to choose an item to plot the expected score curves for.
#' @param integ either an integer defining the number of (ability) groups to integrate the empirical theta vector or the character expression \code{"all"} to plot the empirical theta distribution at the respective item score using symbols (see example).  
#' @param ra an integer, defining the (logit) range for x-axis
#' @param nodes numer of integration nodes 
#' @param lwd see \code{\link{plot}}
#' @param ... arguments passed to plot
#' @examples ########
#' data(bfiN)
#' result <- pers(pair(bfiN))
#' esc(pers_obj=result,1,lwd=2) # plot for first item
#' esc(pers_obj=result,2,lwd=2) # plot for second item
#' for(i in 1:5){esc(pers_obj=result,i,lwd=2)}
#' #########
#' esc(pers_obj=result,2,integ="all",lwd=2) # plot for secod item
#' 
####################################################
####################################################

esc <- function(pers_obj, itemnumber=1, integ=6, ra=4, nodes=100, lwd=2, ... ){

# auslesen pers_obj -----------
responses <- pers_obj$pair$resp
threshold <- na.omit(pers_obj$pair$threshold[itemnumber,]) # theres. des items
itemname <- rownames(pers_obj$pair$threshold)[itemnumber]
d0 <- data.frame(pers_obj$pers$WLE,item.resp=responses[,itemnumber])
d1 <- d0[complete.cases(d0),]
miss <- dim(d0)[1]-dim(d1)[1]
kat <- 0:length(threshold)
theta_wle <- d1[,1]
# range(theta_wle)
npers <- length(theta_wle)
# calculating theoretical theta range from ra and nodes---------------
theta_theo <- seq(-ra,ra,length.out=nodes)
# calculating probabilities---------------
pm <- pvx.matrix(theta_v=theta_theo ,thres=threshold,xm_v=NULL)
# calculating expected score ---------------
expscore <- colSums(kat*pm)

# 1: calculating empirical integrated score and start 'ploting' ----------
if(class(integ)=="numeric"){
group_index <- round(seq(1:npers) / (npers/integ) + 0.5)
d2 <- d1[order(d1[,1]),]
d3 <- aggregate(d2, list(group_index), mean)[,2:3]
# plot-----------
plot(x=range(theta_theo),y=range(kat),type="n",bty="n",ylab="expected Score", xlab= "Theta (ability-logits)",main=paste("Expected Score Curve for Item ",itemname,sep=""),sub=paste("n missing:",miss))
lines(x=theta_theo,y=expscore,type="l", lwd=lwd, ...)
lines(d3, type="o", lwd=lwd+.3, pch=3, ...)
}

# 2: calculating empirical integrated score and start 'ploting' ----------
if(integ=="all"){
d2 <- d1[order(d1[,1]),]
# ff <- apply(as.matrix(d2),1,function(x){paste(x,collapse="_")})
# ff_tab <- table(ff)
scoremean <- aggregate(x=d2[,1], by=list(d2[,2]), FUN=mean)[,c(2,1)]

#fmode <- function(x){  as.numeric(names(table(x))[max(table(x))==table(x)])    }
#scoremode <- aggregate(x=d2[,1], by=list(d2[,2]), FUN=fmode)[,c(2,1)]

ff <- unlist((aggregate(x=d2[,1], by=list(d2[,2]), FUN=table, simplify = F))$x)
# obj1 <- strsplit(x=names(ff_tab), split="_", fixed = TRUE)
# koordinaten <- t(data.frame(lapply(obj1,as.numeric)))
# dim(table(str.pattern(responses)))
#cbind(koordinaten,ff)
koordinaten<- aggregate(d2, d2, mean, na.rm=TRUE)[,c(3,2)] # NA wird dabei rausgeworfen!!!!
# plot-----------
plot(x=range(theta_theo),y=range(kat),type="n",bty="n",ylab="expected Score", xlab= "Theta (ability-logits)",main=paste("Expected Score Curve for Item ",itemname,sep=""),sub=paste("n missing:",miss))
lines(x=theta_theo,y=expscore,type="l", lwd=lwd) # ,...
symbols(koordinaten, circles=ff , inches=.13, add=TRUE,fg="gray60")
lines(scoremean,type="o", lwd=lwd,lty=2) # ,...
#lines(scoremode,type="l", lwd=lwd,lty=3) # ,...
}

}
