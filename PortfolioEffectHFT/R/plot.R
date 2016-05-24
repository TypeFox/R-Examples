setClass("portfolioPlot",
		slots = c(data="data.frame",start.data="data.frame",option="list",bw="logical",breaks="numeric",labels="character")
)

setMethod ("show" , "portfolioPlot",
		function (object){
			p=ggplot() + geom_line(data=object@data, aes_string(x='Time',y='Data',col='legend'),size=object@option$line_size) +
					xlab(NULL) + 
					ylab(NULL) +
					ggtitle(if(!is.null(object@option$subtitle)){
					  bquote(atop(.(object@option$title), atop(italic(.(object@option$subtitle)), "")))
					  }else{
					    object@option$title
					  }) +
					scale_x_continuous(breaks=object@breaks,
							labels=object@labels)+
					util_plotTheme(base_size=object@option$font_size,bw=object@bw,axis.text.size=object@option$axis.text.size,
							title.size=object@option$title.size,has.subtitle=T)+util_colorScheme()
					if(length(unique(object@data$legend))==1){
						p=p+theme( legend.position = "none")
					}
			print(p)
			return(p)
		})

util_ggplot<-function (portfolioPlot){
			p=ggplot() + geom_line(data=portfolioPlot@data, aes_string(x='Time',y='Data',col='legend'),size=portfolioPlot@option$line_size) +
					xlab(NULL) + 
					ylab(NULL) +
					ggtitle(if(!is.null(portfolioPlot@option$subtitle)){
					  bquote(atop(.(portfolioPlot@option$title), atop(italic(.(portfolioPlot@option$subtitle)), "")))
					}else{
					  portfolioPlot@option$title
					}) +
					scale_x_continuous(breaks=portfolioPlot@breaks,
							labels=portfolioPlot@labels)+
					util_plotTheme(base_size=portfolioPlot@option$font_size,bw=portfolioPlot@bw,axis.text.size=portfolioPlot@option$axis.text.size,
							title.size=portfolioPlot@option$title.size,has.subtitle=T)+util_colorScheme()
			if(length(unique(portfolioPlot@data$legend))==1){
				p=p+theme( legend.position = "none")
			}
			return(p)
		}


plotPlot2df<-function (x,y){	
	plotPlot2d(x)
}
#setMethod("plot" ,c(x="portfolio",y="numeric"),plot.ice9portfolio)
setMethod("plot" ,c(x="portfolioPlot",y="missing"),plotPlot2df)

plotPlot2d<-function (portfolioPlot){
			p=ggplot() + geom_line(data=portfolioPlot@data, aes_string(x='Time',y='Data',col='legend'),size=portfolioPlot@option$line_size,col="#004A61") +
					xlab(NULL) + 
					ylab(NULL) +
					ggtitle(if(!is.null(portfolioPlot@option$subtitle)){
					  bquote(atop(.(portfolioPlot@option$title), atop(italic(.(portfolioPlot@option$subtitle)), "")))
					}else{
					  portfolioPlot@option$title
					}) +
					scale_x_continuous(breaks=portfolioPlot@breaks,
							labels=portfolioPlot@labels)+
					util_plotTheme(base_size=portfolioPlot@option$font_size,bw=portfolioPlot@bw,axis.text.size=portfolioPlot@option$axis.text.size,
							title.size=portfolioPlot@option$title.size,has.subtitle=T)+util_colorScheme()
			if(length(unique(portfolioPlot@data$legend))==1){
				p=p+theme( legend.position = "none")
			}
			print(p)
			return(p)
		}

util_plot2d<-function(metric,title=NULL,subtitle = NULL,font_size=10,line_size=1.2,bw=FALSE,legend="",axis.text.size=1.5,title.size=2){
	n=NROW(metric)
	result<-data.frame(value=metric[,2],time=metric[,1],legend=array(legend,dim=n))
	result=result[complete.cases(result),]
	portfolioPlot=util_plot2df(value~time,result,title=title,subtitle = subtitle,font_size=font_size,line_size=line_size,bw=bw,axis.text.size=axis.text.size,title.size=title.size)
	portfolioPlot
}

util_plot2df<-function(formula,data,title=NULL,subtitle=NULL,font_size=10,line_size=1.2,bw=FALSE,axis.text.size=1.5,title.size=2){
  if(is.null(data$legend)){
    stop("Data should have a parameter 'legend'")
  }
  data=data[complete.cases(data),]
  fml = util_parseFormula(formula, data = data)
  if(fml$right.name=="time"){
    Time=sort(unique(fml$right))
    diff<-difftime(as.POSIXlt(max(Time)/1000,origin = "1970-01-01",tz="America/New_York"),
                   as.POSIXlt(min(Time)/1000,origin = "1970-01-01",tz="America/New_York"),units="days")
    diffNum<-which.max(diff> c(3*365,365,90,20,2,0.1,0.03,0.015,0.008,0.0025,0.001,0))
    legends=data$legend
    normTime<-as.POSIXlt(Time/1000,origin = "1970-01-01",tz="America/New_York")
    n<-NROW(normTime)
    TimeStep<-switch(diffNum,
                     "1"= which(c(0,diff(normTime$year))>0),
                     "2"= which((c(0,diff(normTime$mon))!=0)&(c(0,diff(normTime$mon%%3))<0)),
                     "3"= which(c(0,diff(normTime$mon))!=0),
                     "4"= which((c(0,diff(normTime$wday))<0)),
                     "5"= which(c(0,diff(normTime$wday))!=0),
                     "6"= which(c(0,diff(normTime$hour))!=0),
                     "7"= which((c(0,diff(normTime$min))!=0)&(c(0,diff(normTime$min%%10))<0)),
                     "8"= which((c(0,diff(normTime$min))!=0)&(c(0,diff(normTime$min%%5))<0)),
                     "9"= which(c(0,diff(normTime$min))!=0),
                     "10"= which((c(0,diff(normTime$sec))!=0)&(c(0,diff(normTime$sec%%30))<0)),
                     "11"= which((c(0,diff(normTime$sec))!=0)&(c(0,diff(normTime$sec%%10))<0)),
                     "12"= 1:n)
       
		   if(length(TimeStep)==0){
			   tempTimeStep=NULL
			   for(i in 1:12){
				   tempTimeStep=c(tempTimeStep,list(switch(i,
						   "1"= which(c(0,diff(normTime$year))>0),
						   "2"= which((c(0,diff(normTime$mon))!=0)&(c(0,diff(normTime$mon%%3))<0)),
						   "3"= which(c(0,diff(normTime$mon))!=0),
						   "4"= which((c(0,diff(normTime$wday))<0)),
						   "5"= which(c(0,diff(normTime$wday))!=0),
						   "6"= which(c(0,diff(normTime$hour))!=0),
						   "7"= which((c(0,diff(normTime$min))!=0)&(c(0,diff(normTime$min%%10))<0)),
						   "8"= which((c(0,diff(normTime$min))!=0)&(c(0,diff(normTime$min%%5))<0)),
						   "9"= which(c(0,diff(normTime$min))!=0),
						   "10"= which((c(0,diff(normTime$sec))!=0)&(c(0,diff(normTime$sec%%30))<0)),
						   "11"= which((c(0,diff(normTime$sec))!=0)&(c(0,diff(normTime$sec%%10))<0)),
						   "12"= 1:n)))
			   }
			   lengths<-sapply(tempTimeStep,length)
			   s=which(lengths==max(lengths[lengths<10]))[1]
			   TimeStep=tempTimeStep[[s]]
		   }
		   
		   format<-switch(diffNum,
				   "1"= "%Y",
				   "2"= "%m/%Y",
				   "3"= "%m/%Y",
				   "4"="%d %b",
				   "5"= "%d %b",
				   "6"= "%H:%M",
				   "7"= "%H:%M",
				   "8"= "%H:%M",
				   "9"= "%M:%S",
				   "10"= "%M:%S",			
				   "11"= "%M:%S",
				   "12"= "%S")
		   
		   if(((TimeStep[1]-1)/n>0.1)&(format(normTime[TimeStep[1]],format)!=format(normTime[1],format))){
			   TimeStep<-c(1,TimeStep)
		   }
		   if(((n-tail(TimeStep,1))/n>0.06)&(format(normTime[tail(TimeStep,1)],format)!=format(normTime[n],format))){
			   TimeStep<-c(TimeStep,n)
		   }
		   
    
	i=1
	while(i < length(TimeStep)){
		i=i+1
    if((TimeStep[i]-TimeStep[i-1])/n<0.06){
		TimeStep=TimeStep[-i]
		i=1
	}
}
    
    tem=(Time/1000)
	normTimeTemp=normTime
	normTimeTemp$hour=9;
	normTimeTemp$min=30;
	temp=tem-as.numeric(normTimeTemp);
    temp1=which(c(0,diff(temp))<0)
    if(length(temp1)>0){
      for(i in 1:length(temp1)){
        temp[(temp1[i]):length(temp)]=temp[(temp1[i]):length(temp)]+23400
      }
    }
	if(all(diff(temp)==0)){
		temp=1:NROW(Time)
	}
    
    Time=cbind(Time,temp)
    result<-data.frame(Data=fml$left,Time=temp[match(fml$right,Time)],legend=data$legend)
    result2=NULL
    legend=unique(legends)
    for(leg in legend){
      temp1=result[result$legend==leg,]
      temp1=temp1[!is.na(temp1$Data),]
      n=NROW(temp1)
      if(n>=46800){
        s=n%/%23400
        tempSum=c(array(0,dim=s-1),(cumsum(temp1$Data)[-(1:(s-1))]-cumsum(c(0,temp1$Data))[-((n-s+2):(n+1))])/s)
        result2=rbind(result2,data.frame(Data=tempSum[(1:n)%%s==0],Time=temp1$Time[(1:n)%%s==0],legend=leg))	
      }else{
        result2=rbind(result2,temp1)	
      }
    }
    
	portfolioPlot=new("portfolioPlot",data=result2,
			start.data=data.frame(Data=fml$left,Time=fml$right,legend=data$legend),
			option=list(line_size=line_size,axis.text.size=axis.text.size,title.size=title.size,font_size=font_size,title=title,subtitle=subtitle),bw=bw,
			breaks=temp[TimeStep],
			labels=format(normTime[TimeStep],format))
	portfolioPlot
  }else{
    result<-data.frame(x=fml$right,y=fml$left,legend=data$legend)
    p<-ggplot() + geom_line(data=result, aes_string(x='x',y='y',col='legend'),size=line_size) +
      xlab(fml$right.name) + 
      ylab(fml$left.name) +
      ggtitle(if(!is.null(subtitle)){
        bquote(atop(.(title), atop(italic(.(subtitle)), "")))
      }else{
        title
      }) +
      util_plotTheme(base_size=font_size,bw=bw,axis.text.size=axis.text.size,title.size=title.size,has.subtitle=T)+
	  util_colorScheme()
    p
  }
}
		
util_line2d<-function(metric,legend=""){
	
	portfolioPlot=new("portfolioPlot",data=data.frame(),
			start.data=data.frame(Data=metric[,2],Time=metric[,1],legend=array(legend,dim=NROW(metric))),
			option=list(),bw=T,
			breaks=0,
			labels="")
	portfolioPlot
	}

setMethod("+", signature(e1 = "portfolioPlot", e2 = "portfolioPlot"), function(e1, e2) {
				data=rbind(e1@start.data,e2@start.data)
				Time=sort(unique(data$Time))
				diff<-difftime(as.POSIXlt(max(Time)/1000,origin = "1970-01-01",tz="America/New_York"),
						as.POSIXlt(min(Time)/1000,origin = "1970-01-01",tz="America/New_York"),units="days")
				diffNum<-which.max(diff> c(3*365,365,90,20,2,0.1,0.03,0.015,0.008,0.0025,0.001,0))
				
				legends=data$legend
				normTime<-as.POSIXlt(Time/1000,origin = "1970-01-01",tz="America/New_York")
				n<-NROW(normTime)
				TimeStep<-switch(diffNum,
						"1"= which(c(0,diff(normTime$year))>0),
						"2"= which((c(0,diff(normTime$mon))!=0)&(c(0,diff(normTime$mon%%3))<0)),
						"3"= which(c(0,diff(normTime$mon))!=0),
						"4"= which((c(0,diff(normTime$wday))<0)),
						"5"= which(c(0,diff(normTime$wday))!=0),
						"6"= which(c(0,diff(normTime$hour))!=0),
						"7"= which((c(0,diff(normTime$min))!=0)&(c(0,diff(normTime$min%%10))<0)),
						"8"= which((c(0,diff(normTime$min))!=0)&(c(0,diff(normTime$min%%5))<0)),
						"9"= which(c(0,diff(normTime$min))!=0),
						"10"= which((c(0,diff(normTime$sec))!=0)&(c(0,diff(normTime$sec%%30))<0)),
						"11"= which((c(0,diff(normTime$sec))!=0)&(c(0,diff(normTime$sec%%10))<0)),
						"12"= 1:n)
				
				format<-switch(diffNum,
						"1"= "%Y",
						"2"= "%m/%Y",
						"3"= "%m/%Y",
						"4"="%m/%d",
						"5"= "%m/%d",
						"6"= "%H:%M",
						"7"= "%H:%M",
						"8"= "%H:%M",
						"9"= "%M:%S",
						"10"= "%M:%S",			
						"11"= "%M:%S",
						"12"= "%S")
				
				if(((TimeStep[1]-1)/n>0.1)&(format(normTime[TimeStep[1]],format)!=format(normTime[1],format))){
					TimeStep<-c(1,TimeStep)
				}
				if(((n-tail(TimeStep,1))/n>0.06)&(format(normTime[tail(TimeStep,1)],format)!=format(normTime[n],format))){
					TimeStep<-c(TimeStep,n)
				}
				
#				temp=diff(TimeStep)/n>0.06
#				TimeStep=TimeStep[c(temp,TRUE)]
				i=1
				while(i < length(TimeStep)){
					i=i+1
					if((TimeStep[i]-TimeStep[i-1])/n<0.06){
						TimeStep=TimeStep[-i]
						i=1
					}
				}
				
				tem=(Time/1000)
				temp=(tem%%(60*60*24))-48600
				temp1=which(c(0,diff(temp))<0)
				if(any(diff(temp)<0)){
					for(i in 1:length(temp1)){
						temp[(temp1[i]):length(temp)]=temp[(temp1[i]):length(temp)]+23400
					}
				}
				
				Time=cbind(Time,temp)
				result<-data.frame(Data=data$Data,Time=temp[match(data$Time,Time)],legend=legends)
				result2=NULL
				legend=unique(legends)
				for(leg in legend){
					temp1=result[result$legend==leg,]
					temp1=temp1[!is.na(temp1$Data),]
					n=NROW(temp1)
					if(n>46800){
						s=n%/%23400
						tempSum=c(array(0,dim=s-1),(cumsum(temp1$Data)[-(1:(s-1))]-cumsum(c(0,temp1$Data))[-((n-s+2):(n+1))])/s)
						result2=rbind(result2,data.frame(Data=tempSum[(1:n)%%s==0],Time=temp1$Time[(1:n)%%s==0],legend=leg))	
					}else{
						result2=rbind(result2,temp1)	
					}
				}
				
				portfolioPlot=new("portfolioPlot",data=result2,
						start.data=data,
						option=e1@option,bw=e1@bw,
						breaks=temp[TimeStep],
						labels=format(normTime[TimeStep],format))
				portfolioPlot
			})
	



util_summary<-function(portfolio,bw=FALSE){

	portfolioTemp=portfolio_create(portfolio)
	set<-portfolio_getSettings(portfolioTemp)
	Layout <- grid.layout(nrow = 4, ncol = 4,
			widths = unit(c(2, 2, 2,2.27), array("null",dim=4)), 
			heights = unit(c(array(3,dim=2),2.5,3.19), array("null",dim=4)))
	
	grid.newpage()
	pushViewport(viewport(layout = Layout))
	
		tempSet<-set
		symbols<-portfolio_symbols(portfolioTemp)
	.jcall(portfolioTemp@java,returnSig="V", method="createCallGroup",as.integer(2))
	portfolio_startBatch(portfolioTemp)
	portfolio_value(portfolioTemp)
	portfolio_expectedReturn(portfolioTemp)
	portfolio_variance(portfolioTemp)
	portfolio_endBatch(portfolioTemp)
	
	p1<-util_ggplot(util_plot2d(portfolio_value(portfolioTemp),title='Portfolio value ($)',bw=bw))
	p4<-util_ggplot(util_plot2d(portfolio_expectedReturn(portfolioTemp),title="Portfolio Expected Return",bw=bw))+geom_hline(yintercept=0,col='red',size=0.5)
	p5<-util_ggplot(util_plot2d(portfolio_variance(portfolioTemp),title="Portfolio Variance",bw=bw))

	tempSet$resultsSamplingInterval='last'
	portfolio_settings(portfolioTemp,tempSet)
	

	printMatrix1<-matrix(0,ncol=2,nrow=length(symbols))
	
	portfolio_startBatch(portfolioTemp)
	for(symbol in symbols){
		position_weight(portfolioTemp,symbol)
		position_profit(portfolioTemp,symbol)
	}
	if(length(symbols)>5){
	portfolio_profit(portfolioTemp)
}
	portfolio_endBatch(portfolioTemp)
	
	
	
	j<-1
	for(symbol in symbols){
		printMatrix1[j,1]<-round(position_weight(portfolioTemp,symbol)[2]*100,digits =3)
		printMatrix1[j,2]<-round(position_profit(portfolioTemp,symbol)[2],digits =3)
		j<-j+1		}
	rownames(printMatrix1)=symbols
	
	if(length(symbols)>1){
	symbols<-symbols[order(printMatrix1[,1],decreasing=TRUE)]
	printMatrix<-printMatrix1[symbols,]
	}else{
		printMatrix<-printMatrix1
	}
	
	if(length(symbols)>5){
		symbols<-c(symbols[1:4],"Other")
		printMatrix1<-printMatrix1[symbols[1:4],]
		other<-c(0,0)
		other[1]<-100-sum(printMatrix1[1:4,1])
		other[2]<-round(portfolio_profit(portfolioTemp)[2],digits =3)-sum(printMatrix1[,2])
		printMatrix1<-rbind(printMatrix1,other)
	}
	
	result3<-data.frame(data=printMatrix[,1],symbols=symbols,legend="Time2")
	
	
	p2<-ggplot(data=result3, aes_string(x='symbols',y='data')) + geom_bar(stat="identity",position="dodge",fill="#01526D")+ coord_flip() +
			scale_y_continuous(labels = per_for())+
			xlab(NULL) + 
			ylab(NULL) +
			xlim(rev(symbols))+
			ggtitle('Position Weight (%)')+util_plotTheme(axis.text.size=1.5,bw=bw) 
	
	
	
	result4<-data.frame(data=printMatrix[,2],symbols=symbols,legend="Time2")
	
	
	p3<-ggplot(data=result4, aes_string(x='symbols',y='data')) + geom_bar(stat="identity",position="dodge",fill="#00A3DC")+ coord_flip() +
#			scale_y_continuous(labels = per_for(2))+
	        xlab(NULL) + 
			ylab(NULL) +
			xlim(rev(symbols))+
			ggtitle('Position Profit ($)')+util_plotTheme(axis.text.size=1.5,bw=bw) + scale_fill_brewer()
	
	portfolio_settings(portfolioTemp,set)
	
	print(p1, vp = viewport(layout.pos.row = 1,
					layout.pos.col = 1:4))
	print(p2, vp = viewport(layout.pos.row = 2,
					layout.pos.col = 1:2))
	print(p3, vp = viewport(layout.pos.row = 2,
					layout.pos.col = 3:4))
	print(p4, vp = viewport(layout.pos.row = 3,
					layout.pos.col = 1:4))
	print(p5, vp = viewport(layout.pos.row = 4,
					layout.pos.col = 1:4))
}

util_plotTheme<-function (base_size = 10, base_family = "sans", horizontal = TRUE, 
                            dkpanel = FALSE, bw = FALSE,axis.text.size=1.5,title.size=2, has.subtitle = FALSE) 
{
  if(bw){
    bgcolors <- c("white","#B2BFCB","gray")	
  }else{
    bgcolors <- c("#d5e4eb" ,"#c3d6df","white")
  }
  names(bgcolors)<-c("ebg","edkbg","line")
  
  ret <- theme(line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"), 
               rect = element_rect(fill = bgcolors["ebg"],colour = NA, size = 0.5, linetype = 1), 
               text = element_text(family = base_family,face = "plain", colour = "black", size = base_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1), 
               axis.text = element_text(size = rel(axis.text.size)), 
               axis.line = element_line(size = rel(0.8)), axis.line.y = element_blank(), 
               axis.text.x = element_text(vjust = 1), axis.text.y = element_text(hjust = 0), 
               axis.ticks = element_line(), axis.ticks.y = element_blank(), 
               axis.title = element_text(size = rel(1)), axis.title.x = element_text(size=base_size*1.5), 
               axis.title.y = element_text(angle = 90,size=base_size*1.5), axis.ticks.length = unit(base_size * 
                                                                                                      0.5, "points"),
               legend.margin = unit(-0.1, "cm"), legend.key = element_rect(linetype = 0,fill = bgcolors["ebg"]), 
               legend.key.size = unit(1.2, "lines"), legend.key.height = NULL, 
               legend.key.width = NULL, legend.text = element_text(size = rel(1.25)), 
               legend.text.align = NULL,legend.title=element_blank(), legend.title.align = NULL,
               legend.direction = NULL, legend.justification = "center", legend.position = "top",
               panel.background = element_rect(linetype = 0,fill = bgcolors["ebg"]), panel.border = element_blank(), 
               panel.grid.major = element_line(colour = bgcolors["line"], size = rel(1)), 
               panel.grid.minor = element_blank(), panel.margin = unit(0.25, 
                                                                       "lines"), strip.background = element_rect(fill = bgcolors["ebg"], 
                                                                                                                 colour = NA, linetype = 0), strip.text = element_text(size = rel(1.25)), 
               strip.text.x = element_text(), strip.text.y = element_text(angle = -90), 
               plot.background = element_rect(fill = bgcolors["ebg"], 
                                              colour = NA), plot.title = element_text(size = rel(title.size)), plot.margin = unit(c(0, 
                                                                                                                       5, 6, 5) * 2, "points"))
  if (horizontal) {
    ret <- ret + theme(panel.grid.major.x = element_blank())
  }
  else {
    ret <- ret + theme(panel.grid.major.y = element_blank())
  }
  
  if (has.subtitle) {
    ret <- ret + theme(plot.title = element_text(size = rel(title.size), vjust = -1))
  }
  
  if (dkpanel == TRUE) {
    ret <- ret + theme(panel.background = element_rect(fill = bgcolors["edkbg"]), 
                       strip.background = element_rect(fill = bgcolors["edkbg"]))
  }
  ret
}


util_plotDensity<-function(PDF, bw=FALSE){
	df<-data.frame(x=c(PDF$value),y=c(PDF$pdf),col="Return Density")
	p<-ggplot() +util_plotTheme(bw=bw)+scale_colour_manual(values=c("#004A61"))
	if(!is.null(PDF$pdfNormal)){
		dfNormal<-data.frame(x=c(PDF$valueNormal),y=c(PDF$pdfNormal),fill="Normal Density")
		p<-p+geom_area(data=dfNormal, aes_string(x='x',y='y',fill='fill'),size=1,alpha=0.5)+geom_line(data=dfNormal, aes_string(x='x',y='y'),size=0.5,col="red")
	}
	p+geom_line(data=df, aes_string(x='x',y='y',col='col'),size=1.5)
}
per_for<-function(digits=0){
		function(x) {
			x <- round(x,digits=digits)
			paste0(format(x,big.mark = ",", scientific = FALSE, trim = TRUE), "%")
		}
}

util_colorScheme<- function() {
	discrete_scale("colour", "portfolioeffect", portfolio_pal())
}

util_fillScheme<- function() {
	discrete_scale("fill", "portfolioeffect", portfolio_pal())
}

portfolio_pal<-function(){
	function(n) {
		return(c("#014d64","#01a2d9","#00897E","#ee8f71","#7c260b","#adadad","#6794a7","#7ad2f6","#00887d","#76c0c1"))
	}
}

util_multiplot <- function(..., cols=1) {
	
	# Make a list from the ... arguments and plotlist
	plots <- list(...)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))

	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}
util_parseFormula<-function(model,data){
	if(length(model)!=3){
		stop('Failed to parse provided formula.')
	}
	if(!(paste(model[[2]]) %in% colnames(data))){
		stop(paste("Data should have a parameter" ,paste(model[[2]])))
	}
	if(!(paste(model[[3]]) %in% colnames(data))){
		stop(paste("Data should have a parameter" ,paste(model[[3]])))
	}
	result=list(left=data[[paste(model[[2]])]], right=data[[paste(model[[3]])]],  left.name=paste(model[[2]]),  right.name=paste(model[[3]]))
}