library('geostatsp')
data('swissRain')

Ncores = c(1,2)[1+(.Platform$OS.type=='unix')]



sr2 = swissRain
sr2$elev = raster::extract(swissAltitude, sr2)
swissFit = likfitLgm(
    data=sr2, 
    formula=rain~ elev,
    param=c(range=10000,shape=1,nugget=0,boxcox=0.5,anisoRatio=2,anisoAngleDegrees=45),
    paramToEstimate = c("range",'anisoAngleDegrees','anisoRatio'),
    reml=FALSE,
    verbose=FALSE
)


# calculate log-likelihood at the MLE's, but re-estimate variance
sl = loglikLgm(
    swissFit$param[c('range','shape','boxcox', 'anisoRatio', 'anisoAngleRadians')],
    data=sr2, 
    formula=rain~ elev,
    reml=swissFit$model$reml)  


# calculate log-likelihood without re-estimating variance
sigSqHat = attributes(sl)$totalVarHat
sl1 = loglikLgm(
    c(attributes(sl)$param[
            c('boxcox','anisoRatio','anisoAngleRadians','shape', 'range')], 
        variance=sigSqHat),
    data=sr2, 
    formula=rain~ elev,
    reml=swissFit$model$reml)  
 

# re=estimate the anisotropy parameters but not the range
sf2 = likfitLgm(
    data=swissFit$data, 
    formula=swissFit$model$formula,
    param= swissFit$param[c('range','nugget','shape','boxcox', 'anisoRatio', 'anisoAngleRadians')],
    paramToEstimate = c('variance','anisoAngleRadians','anisoRatio'),
    reml=swissFit$model$reml)  

# these should all be the same
as.numeric(sl1)
as.numeric(sl) 
swissFit$optim$logL
sf2$optim$logL

date()
x=profLlgm(swissFit, mc.cores=Ncores,
    range=seq(15000, 55000 , len=12)
)
date()

 
swissInf = informationLgm(swissFit)



if(!interactive()) 
  pdf("profLswissAngle.pdf")

plot(x[[1]],x[[2]], xlab=names(x)[1],
#		yaxt='n',
		ylab='log L',
		ylim=c(min(x[[2]]),x$maxLogL),
		type='n')
lines(x[[1]],x[[2]])
abline(h=x$breaks[-1],
		col=x$col,
		lwd=1.5)
axis(2,at=x$breaks,labels=x$prob,
		line=-1.2,tick=F,
		las=1,padj=1.2,hadj=0,col.axis='red')

abline(v=x$ciLong$par,
		lty=2,
		col=x$col[as.character(x$ciLong$prob)])


axis(1,at=x$ciLong$par,
		labels=x$ciLong$quantile,
		padj= -6,hadj=0.5, 
		tcl=0.5,cex.axis=0.8,
		col=NA,col.ticks='red',col.axis='red')

ciCols = grep("^ci", colnames(swissInf$summary),
		value=TRUE)

ciValues = unlist(swissInf$summary[intersect(names(x), rownames(swissInf$summary)),ciCols])
if(any(!is.na(ciValues)))
  axis(1,at=ciValues,
		labels=gsub("^ci","",ciCols),
		padj= 2,hadj=0.5, 
		tcl=-2,cex.axis=0.7,
		col=NA,col.ticks='blue',col.axis='blue')

lines(x[[1]],x[[2]], type='o')

if(!interactive()) 
  dev.off()


if(interactive()  | Sys.info()['user'] =='patrick') {
  Ncores = c(1,2)[1+(.Platform$OS.type=='unix')]
  date()
  x2d=profLlgm(swissFit, mc.cores=Ncores,
      anisoAngleRadians=seq(22, 60 , len=24)*(2*pi/360),
		anisoRatio =exp(seq(log(1.5),log(20),len=36))
  )
  date()
if(!interactive()) 
  pdf("profLswiss2d.pdf")
image(x2d[[1]],x2d[[2]],x2d[[3]],
		breaks=x2d$breaks,
		col=x2d$col,log='y',
		xlab=names(x2d)[1],
		ylab=names(x2d)[2])

thesevars = c("anisoAngleRadians","log(anisoRatio)")
thisV = swissInf$information[
		thesevars,thesevars]
thisMean= c(x2d$MLE["anisoAngleRadians"],
		log(x2d$MLE['anisoRatio']))


if(requireNamespace("ellipse", quietly=TRUE)) {

for(D in x2d$prob[x2d$prob>0&x2d$prob<1]) {
	thisE = ellipse::ellipse(thisV, centre=thisMean,
			level=D)
  colnames(thisE) = names(thisMean)
	thisE = cbind(thisE,
			anisoRatioExp = exp(thisE[,"anisoRatio"]))
	lines(thisE[,"anisoAngleRadians"],
			thisE[,"anisoRatioExp"],lwd=4)
	lines(thisE[,"anisoAngleRadians"],
			thisE[,"anisoRatioExp"], col=x2d$col[as.character(D)],
			lwd=3)
}
}

points(x2d$MLE[1],x2d$MLE[2],pch=15) 

legend("topright", fill=x2d$col, legend=x2d$prob[-length(x2d$prob)])

if(!interactive()) 
  dev.off()
}