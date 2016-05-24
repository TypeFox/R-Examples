#' @title Graphical Item Fit Plots
#' @export gif
#' @description plotting function for plotting empirical and model derived category probability curves.
#' @details no details in the moment.
#' @param pers_obj an object of class \code{"pers"} as a result from function \code{\link{pers}}.
#' @param itemnumber an integer, defining the number of the item to plot the respective categoy probability for. This is set to an arbitrary default value of \code{itemnumber = 1} to avoid error messages when you forget to choose an item to plot the expected score curves for.
#' @param ra an integer, defining the (logit) range for x-axis
#' @param integ either an integer, defining the number of integration points along the (logit) range on the x-axis to integrate the empirical theta values, or the character expression \code{"raw"} (default) which will use the rawscore groups as integration points.
#' @param kat either an integer, defining for which category the empirical category probabilities should be plotted over the model derived category probability curves, or the character expression \code{"all"} (default) which will plot the empirical category probabilities for all categories.
#' @param ... arguments passed to plot
#' @return a plot with category probabilities.
#' @examples ########
#' data(bfiN)
#' pers_obj <- pers(pair(bfiN))
#' #### plot empirical category probabilities
#' gif(pers_obj = pers_obj, itemnumber = 1 )
#' gif(pers_obj = pers_obj, itemnumber = 1 , integ=8) # integration over 8 points
#' gif(pers_obj = pers_obj, itemnumber = 1 , integ=8, kat=1) # only for category number 1
####################################################


####################################################

gif <- function(pers_obj, itemnumber=1, ra=4, integ="raw", kat="all",...  ){
  itemname <- rownames(pers_obj$pair$threshold)[itemnumber]
  WLE <- pers_obj$pers$WLE
  raw <- pers_obj$pers$raw
  responsesi <- pers_obj$pair$resp[,itemnumber]
  d0 <- data.frame(WLE,item.resp=responsesi,raw=raw)
  d1 <- d0[complete.cases(d0),]
  miss <- dim(d0)[1]-dim(d1)[1] # n missing auf WLE und responsesi
  mi <- ((sapply(pers_obj$pair$sb,length))[itemnumber]) # Anz. Katg. -1 für item
  d2 <- d1[order(d1[,1]),] # sortieren nach WLE
  WLEsort <- d2[,1] # sortiert nach WLE
  responsesi_sort <- d2[,2] # sortiert nach WLE
  raw_sort <- d2[,3]
  npers <- length(WLEsort)
  #interv <- round(seq(1:npers) / (npers/integ) + 0.49999)
  if(class(integ)=="numeric"){interv <- as.numeric(cut(WLEsort,integ))}
  if(integ=="raw"){interv <- raw_sort+1}
 
  # start calculating empirical cat probs ------------  
  Lresp <- list()
  for (i in 1:(max(interv))){Lresp[[i]] <- responsesi_sort[interv==i]}
  Lftab <- lapply(Lresp, ftab, catgories = 0:mi) # Kateg. häufig. in bre theta bereichen
  Ln <- lapply(Lftab, function(x){ sum(x[2:(mi+2),1])}) # n Pers. ohne NA in interv theta bereichen 
  ecp <- mapply(function(x,y){x[2:(mi+2),1] / y }, Lftab, Ln,SIMPLIFY = TRUE)
  xkord <- (aggregate(WLEsort,list(interv),FUN=median))[,2] # ggf. noch variabel median oder mean
  

# start 'ploting' function ------------    
p <- round(t(sapply(seq(-ra,ra,length.out=(ra*2+1)*10), function(x){pvx(theta=x,thres=na.omit(pers_obj$pair$threshold[itemnumber,]))} )),9)
rownames(p) <- seq(-ra,ra,length.out=(ra*2+1)*10)
dim(p) 
##### plotting model derived category probabilities-------------------------
plot(y=p[,1], x= seq(-ra,ra,length.out=(ra*2+1)*10), ylim=c(0,1) ,bty="n",type="n",xaxt="n", ylab=("p"), xlab="logits", main=paste("model derived and empirical \n category probability curves for item ", itemname,sep=""),...)
for (i in 1:dim(p)[2]){
  lines(y=p[,i], x=seq(-ra,ra,length.out=(ra*2+1)*10),lty="dotted",col="grey50",...)#, ...
}
#pos <- (1:((ra*2+1)*10))[((1:length(-ra:ra))*10)-ra]
ticks <- unique(round(seq(-ra,ra,length.out=(ra*2+1)*10)))

length(ticks)
length(-ra:ra)
axis(1,ticks,labels=as.character(-ra:ra))
# abline(h = 0, v = na.omit(pair_obj$threshold[itemnumber,]), col = "gray60")
hohe <- apply(p,2,max)
lage <- apply(p,2,function(x){as.numeric(rownames(p)[max(x)==x])  }    )
text(lage,hohe,labels=names(lage),pos=1,cex=.7)

##### plotting empirical category probabilities-------------------------
# matlines(x=xkord,y=t(ecp),col="black")
if(class(kat)=="character"){if(kat=="all"){von=1; bis=(mi+1); ckat=von:bis}}
if(class(kat)=="numeric" | class(kat)=="integer"){
 if(min(kat)<1){stop("there is no category number smaler than 1")} 
 if(max(kat)>(mi+1)){stop("there are only ", mi+1," categories but you requested category number:", max(kat))}
 von <- 1 
 bis <- length(kat)
 ckat <- kat
}
for(i in von:bis){lines(x=xkord,y=ecp[ckat[i],],col="black",...)}

for(i in von:bis){text(xkord, ecp[ckat[i],] ,labels=(names(lage))[ckat[i]],cex=.7,...)}  

}
kat=c(1,3,6,8)


