setGeneric("crossVal",
		function(object, x12Parameter=new("x12Parameter"),
				x12BaseInfo=new("x12BaseInfo"),...) { standardGeneric("crossVal")} )
setMethod(
		f='crossVal',
		signature=signature(object = "ts"),
		definition=function(object, x12Parameter,x12BaseInfo,
				showCI=FALSE,main="Cross Validation",
				col_original="black",col_fc="#2020ff",col_bc="#2020ff",
				col_ci="#d1d1ff",col_cishade="#d1d1ff",
				lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,
				lwd_original=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,ytop=1,
				points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
				showLine=TRUE,col_line="grey",lty_line=3,
				ylab="Value",xlab="Date",ylim=NULL,span=NULL) {
		
		fbcastP<-getP(x12Parameter,whichP=list("forecast_years","backcast_years"))
		forecast_years<-fbcastP$forecast_years
		backcast_years<-fbcastP$backcast_years
		
			start_ts <- start(object)
			freq_ts<-frequency(object)
			orig_ts<-object
			end_ts <- end(object)
			
		if(!is.null(forecast_years) && forecast_years!=0){
			forecast_years<-forecast_years*freq_ts
			object <- object[-((length(object)-forecast_years+1):length(object))]			
			object <- ts(object,start=start_ts,frequency=freq_ts)
			end_ts <- end(object)
			forecast=TRUE
			if(points_fc)
				addpoints_fc<-TRUE
		}else{
		forecast=FALSE
		forecast_years=0
		if(points_fc)
			addpoints_fc<-FALSE
		}
		if(!is.null(backcast_years) && backcast_years!=0){	
			backcast_years<-backcast_years*freq_ts
			object <- object[-(1:backcast_years)]
			object <- ts(object,end=end_ts,freq=freq_ts)
			if(points_bc)
				addpoints_bc<-TRUE
		backcast=TRUE}else{
		backcast=FALSE
		backcast_years=0
		if(points_bc)
			addpoints_bc<-FALSE
		}

#		x12Paramter<-setP(x12Parameter,listP=list(forecast_years=forecast_years,backcast_years=backcast_years))
#if(any(file.exists(grep(basename("Rout"),list.files(dirname("Rout")),value=TRUE))))
#			cat(file.exists(grep(basename("Rout"),list.files(dirname("Rout")),value=TRUE)))
#new.env
#			file.copy(from=grep(basename("Rout"),list.files(dirname("Rout")),value=TRUE),to=paste(basename("Rout"),"Temp",gsub("Rout","",grep(basename("Rout"),list.files(dirname("Rout")),value=TRUE)),sep=""))
#print(getwd())
#dir.create("M:/Meraner/Workspace/Saisonbereinigung_Test/x12Test/tmp")
#file.copy(from=grep(basename("Rout"),list.files(dirname("Rout")),value=TRUE)[1],to="M:/Meraner/Workspace/Saisonbereinigung_Test/x12Test")
olddir<-getwd()
setwd(tempdir())
tryCatch(cvout <- x12(object,x12Parameter,x12BaseInfo),
		finally={if(!exists("cvout"))
			cat("=> No cross validation can be performed!\n")
		})
		
setwd(olddir)			
bc<-cvout@backcast@estimate
fc<-cvout@forecast@estimate
object <- object
freq_ts<-freq_ts
ts_plot<- ts(c(bc,object,fc),start=start(bc),end=end(fc),frequency=freq_ts)
ts.lower<-ts(c(cvout@backcast@lowerci,object,cvout@forecast@lowerci),start=start(cvout@backcast@lowerci),end=end(cvout@forecast@lowerci),frequency=freq_ts)
ts.upper<-ts(c(cvout@backcast@upperci,object,cvout@forecast@upperci),start=start(cvout@backcast@upperci),end=end(cvout@forecast@upperci),frequency=freq_ts)

if(is.null(span)){
			xlim <- NULL
			if(!showCI){
				limits.y.lower<-min(orig_ts,bc,fc,na.rm=TRUE)
				limits.y.upper<-max(orig_ts,bc,fc,na.rm=TRUE)
				limits.y<-c(limits.y.lower,limits.y.upper)
			}
			if(showCI){
				limits.y.lower<-min(orig_ts,bc,fc,cvout@forecast@lowerci,cvout@backcast@lowerci,cvout@forecast@upperci,cvout@backcast@upperci,na.rm=TRUE)
				limits.y.upper<-max(orig_ts,bc,fc,cvout@forecast@lowerci,cvout@backcast@lowerci,cvout@forecast@upperci,cvout@backcast@upperci,na.rm=TRUE)
				limits.y<-c(limits.y.lower,limits.y.upper)
			}				
		}else{
			if(length(span)==4){
				if(any(!c(span[2],span[4]))%in%c(1:12))
					stop("Span argument wrong!")
				xlim <- c(span[1]+(span[2]-1)/freq_ts,span[3]+(span[4]-1)/freq_ts)
					limits.y<-c(min(window(ts_plot,span[1:2],span[3:4]),window(orig_ts,span[1:2],span[3:4]),na.rm=TRUE),max(window(ts_plot,span[1:2],span[3:4]),window(orig_ts,span[1:2],span[3:4]),na.rm=TRUE)*ytop)
			if(showCI)				
				limits.y<-c(min(limits.y[1],window(ts.lower,span[1:2],span[3:4]),window(ts.upper,span[1:2],span[3:4]),na.rm=TRUE),
						max(limits.y[2],window(ts.lower,span[1:2],span[3:4]),window(ts.upper,span[1:2],span[3:4]),na.rm=TRUE))
			
				}else if(length(span)==2){
				xlim <- span
				limits.y <- c(min(window(ts_plot,c(span[1],1),c(span[2],1)),window(orig_ts,c(span[1],1),c(span[2],1)),na.rm=TRUE),max(window(ts_plot,span[1],span[2]),window(orig_ts,span[1],span[2]),na.rm=TRUE)*ytop)
				if(showCI)	
				limits.y<-c(min(limits.y[1],window(ts.lower,c(span[1],1),c(span[2],1)),window(ts.upper,c(span[1],1),c(span[2],1)),na.rm=TRUE),
						max(limits.y[2],window(ts.lower,span[1],span[2]),window(ts.upper,span[1],span[2]),na.rm=TRUE))
			
			}else
				stop("Span argument wrong!")
			
		}
		
		if(!is.null(ylim))
		limits.y<-ylim
	
		ts<-plotFbcast(cvout,backcast=backcast,forecast=forecast,
				showCI=showCI,main=main,
				col_original=col_original,col_fc=col_fc,col_bc=col_bc,
				col_ci=col_ci,col_cishade=col_cishade,
				lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,
				lwd_original=lwd_original,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
				ytop=ytop,points_bc=points_bc,points_fc=points_fc,
				showLine=showLine,col_line=col_line,lty_line=lty_line,
				ylab=ylab,xlab=xlab,points_original=points_original,ylim=limits.y,xlim=xlim)
		lines(orig_ts,col=col_original)
		if(!points_fc)
			addpoints_fc<-FALSE
		if(!points_bc)
			addpoints_bc<-FALSE
		if(addpoints_fc)
			points(x=time(cvout@forecast@estimate),y=orig_ts[((length(orig_ts)-forecast_years+1):length(orig_ts))],col=col_original)
		if(addpoints_bc)
			points(time(cvout@backcast@estimate),orig_ts[(1:backcast_years)],col=col_original)
		
		aT <- aL <- axTicks(1)
		if(!is.null(xlim))
			tp <- expand.grid(floor(xlim[1]):ceiling(xlim[2]),(0:(frequency(ts)-1))/frequency(ts))
		else
			tp <- expand.grid(floor(time(ts)[1]):ceiling(time(ts)[length(ts)]),(0:(frequency(ts)-1))/frequency(ts))			
		mm <- round(tp[,2]*frequency(ts))
		yy <- tp[,1]
		tp <- tp[,1]+tp[,2]
		for(i in 1:length(aT)){
			ii <- which.min(abs(tp-aT[i]))
			aT[i] <- tp[ii]
			if(mm[ii]<9)
				aL[i] <- yy[ii]+(mm[ii]+1)/10
			else
				aL[i] <- yy[ii]+(mm[ii]+1)/100
		}
		axis(1,at=aT,labels=aL)	
		
if(backcast){
res.bc<-as.data.frame(rbind(orig_ts[(1:backcast_years)],cvout@backcast@estimate),row.names=c("original","backcast"))
colnames(res.bc)<-(1:backcast_years)
}
if(forecast){
res.fc<-as.data.frame(rbind(orig_ts[((length(orig_ts)-forecast_years+1):length(orig_ts))],cvout@forecast@estimate),row.names=c("original","forecast"))
colnames(res.fc)<-((length(orig_ts)-forecast_years+1):length(orig_ts))
}


#file.remove(grep(basename("Rout"),list.files(dirname("Rout")),value=TRUE,fixed=TRUE))
#if(any(file.exists(grep(basename("RoutTemp"),list.files(dirname("RoutTemp")),value=TRUE))))
#file.rename(from=grep(basename("RoutTemp"),list.files(dirname("RoutTemp")),value=TRUE),to=paste("Rout",gsub("RoutTemp","",grep(basename("RoutTemp"),list.files(dirname("RoutTemp")),value=TRUE)),sep=""))

	crossVal <- new("crossValidation")
	if(backcast)
	crossVal@backcast<-res.bc
	if(forecast)
	crossVal@forecast<-res.fc
	invisible(crossVal)

}
)
setMethod(
		f='crossVal',
		signature=signature(object = "x12Single"),
		definition=function(object,x12BaseInfo=new("x12BaseInfo"),
				showCI=FALSE,main="Cross Validation",
				col_original="black",col_fc="#2020ff",col_bc="#2020ff",
				col_ci="#d1d1ff",col_cishade="#d1d1ff",
				lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,
				lwd_original=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,ytop=1,
				points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
				showLine=TRUE,col_line="grey",lty_line=3,
				ylab="Value",xlab="Date",ylim=NULL,span=NULL) {
			
			crossVal(object@ts,object@x12Parameter,
				showCI=showCI,main=main,
				col_original=col_original,col_fc=col_fc,col_bc=col_bc,
				col_ci=col_ci,col_cishade=col_cishade,
				lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,
				lwd_original=lwd_original,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
				ytop=ytop,points_bc=points_bc,points_fc=points_fc,
				showLine=showLine,col_line=col_line,lty_line=lty_line,
				ylab=ylab,xlab=xlab,points_original=points_original,ylim=ylim,span=span)
		}
)

