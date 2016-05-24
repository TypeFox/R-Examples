#' A function to perform Multiple Resolution Goodness of Fit.
#' 
#' Returns the results of a Multiple Resolution Goodness of Fit after the modified method of Kuhnert et al. 2005, originally by Costanza 1989. This function is computationally intensive and consumes a lot of RAM memory.
#'
#' @param year year validation year (by default 1992)
#'
#' @param w1 w1 initial window size (by default 1)
#'
#' @param w2 w2 final window size (by default 27; max. 113)
#'
#' @param k k parameter for weighting Ft with lower/larger windows resolutions
#'
#' @return Fw vector of fits for each window size
#'
#' @return Ft weighted overall fit
#'
#' @keywords Multiple Resolution Goodness Fit
#'
#' @export
#' 
#' @examples
#' ## Not run mrgf(year='1992',w1=1,w2=113,k=0)

mrgf<-function(year='1992',w1=1,w2=27,k=0){ # ,a='',b=''
rpath = paste(find.package('spdynmod'),'/extdata',sep='')
#require(raster)
af<-raster::raster(paste(rpath,'/y',year,'_rs.asc',sep=''))
bf<-raster::raster(paste(rpath,'/y',year,'_mod.asc',sep=''))
af@nrows->nr
af@ncols->nc
print(nr)
print(nc)

ft<<-NULL
fwes<-NULL
ees<-NULL
fw<<-NULL
#fitw<-NULL
g<-0
for(w in seq(w1,w2,2)){
g<-g+1
raster::getValuesFocal(af,ngb=w)->valss1
raster::getValuesFocal(bf,ngb=w)->valss2
print(paste(w,dim(valss1)[2]))
h<-0
fits<-NULL
for(j in 1:dim(valss1)[1]){
ads<-NULL
valsx1<-valss1[j,]
valsx2<-valss2[j,]
vals1<-valsx1[!is.na(valsx1)]
vals2<-valsx2[!is.na(valsx2)]
lvals1<-length(vals1)
cats<-unique(vals1)
lcats<-length(cats)
if(lcats>0){
 h<-h+1
 for(i in 1:lcats){
  length(vals1[vals1==cats[i]])->a1
  length(vals2[vals2==cats[i]])->a2
  ad<-abs(a1-a2)
  ads[i]<-ad
 }
 fit<-1-(sum(ads,na.rm=T)/(2*lvals1))
 fits[h]<-fit
}
}
fitw<-sum(fits,na.rm=T)/length(fits)
fwe<-fitw*exp(-k*(w-1))
ee<-exp(-k*(w-1))
fw[g]<<-fitw
fwes[g]<-fwe
ees[g]<-ee
}
sfwes<-sum(fwes)
sess<-sum(ees)
ft<<-sfwes/sess

ft<-get('ft')
fw<-get('fw')

#png(paste('mrgf_',year,'.png',sep=''))
graphics::plot(seq(w1,w2,2),fw,ylim=c(0,1),xlab='window size',ylab='Fw',ty='l',main='Multiple Resolution Goodness of Fit',sub=paste('Ft =',round(ft,2),'; k = ',k))
mtext(year,side=3)
#dev.off()
return(ft)
}
