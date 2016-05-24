setMethod(f='plot',
    signature=signature(x = "x12Output"),
    definition=function(x,original=TRUE,sa=FALSE,trend=FALSE,
        log_transform=FALSE,
        ylab="Value",xlab="Date",
        main="TS",
        col_original="black",col_sa="blue",col_trend="green",
        lwd_original=1,lwd_sa=1,lwd_trend=1,lty_sa=1,lty_trend=1,
        ytop=1,showAllout=FALSE,showAlloutLines=FALSE,showOut=NULL,annComp=TRUE,annCompTrend=TRUE,
        col_ao="red",col_ls="red",col_tc="red",col_annComp="grey",lwd_out=1,cex_out=1.5,
        pch_ao=4,pch_ls=2,pch_tc=23,plot_legend=TRUE,legend_horiz=TRUE,legend_bty="o",
        ### implement plotFbcast
        forecast=FALSE,backcast=FALSE,showCI=TRUE,
        col_fc="#2020ff",col_bc="#2020ff",col_ci="#d1d1ff",col_cishade="#d1d1ff",
        lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,
        points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
        showLine=FALSE,col_line="grey",lty_line=3,ylim=NULL,span=NULL,...
    ) 
    {
      if(showAllout)
        legend_horiz <- FALSE
      if(is.null(span))
        xlim <- NULL
      else{
        if(length(span)==4){
          if(any(!c(span[2],span[4]))%in%c(1:12))
            stop("Span argument wrong!")
          xlim <- c(span[1]+(span[2]-1)/frequency(x@a1),span[3]+(span[4]-1)/frequency(x@a1))
          if(is.null(ylim)){
            ylim <- c(min(window(x@a1,span[1:2],span[3:4]),na.rm=TRUE),max(window(x@a1,span[1:2],span[3:4]),na.rm=TRUE)*ytop)
			if(log_transform)
			ylim <- c(min(log(window(x@a1,span[1:2],span[3:4])),na.rm=TRUE),max(log(window(x@a1,span[1:2],span[3:4])),na.rm=TRUE)*ytop)	
		}
		}else if(length(span)==2){
          xlim <- span
          ylim <- c(min(window(x@a1,c(span[1],1),c(span[2],1)),na.rm=TRUE),max(window(x@a1,span[1],span[2]),na.rm=TRUE)*ytop)
		  if(log_transform)
			ylim <- c(min(log(window(x@a1,c(span[1],1),c(span[2],1))),na.rm=TRUE),max(log(window(x@a1,span[1],span[2])),na.rm=TRUE)*ytop)		  
	  }else
          stop("Span argument wrong!")
        
      }
      #Achtung: plotFbcast 
      #keine Option log_transform
      #object und nicht x 
      ##
#			c(object,showCI=TRUE,
#			main="Time Series",forecast=TRUE,backcast=TRUE,
#			col_original="black",col_fc="#2020ff",col_bc="#2020ff",
#			col_ci="#d1d1ff",col_cishade="#d1d1ff",
#			lty_original=1,lty_fc=1,lty_bc=1,lty_ci=1,
#			lwd_original=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,ytop=1,
#			points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
#			showLine=FALSE,col_line="grey",lty_line=3,
#			ylab="Value",xlab="Date",ylim=NULL,...)
      showWarnings=TRUE # fuer plotFbcast
      object <- x
      leg.txt <- vector()
      leg.col <- vector()
      leg.lty <- vector()
      user.main.logical<-FALSE
      if(main!="TS"){
        user.main.logical<-TRUE
        user.main<-main}
#			if(!is.null(ylim))
#				span<-ylim
      if(!is.null(showOut))
        showAllout<-FALSE	
      gp<-par()	
      for(i in c("cin","cra","csi","cxy","din","pin")){
        gp <- gp[-which(names(gp)%in%i)]			
      }
      tryCatch({
            if(original){
              if(!log_transform){
                ts <- object@a1	
                main<-main.orig <- "Original Series"
                leg.txt <- c(leg.txt,"Original")
                leg.col <- c(leg.col,col_original)
                leg.lty <- c(leg.lty,lty_original)
                
              }else{
                ts <- log(object@a1)
                main<-main.orig<- "Log transformed Original Series"
                leg.txt <- c(leg.txt,"Original")
                leg.col <- c(leg.col,col_original)
                leg.lty <- c(leg.lty,lty_original)
                
              }}
            if(sa){
              if(!log_transform){
                ts.sa <- object@d11	
                main<-"Seasonally Adjusted Series"	
                leg.txt <- c(leg.txt,"Seasonally Adjusted")
                leg.col <- c(leg.col,col_sa)
                leg.lty <- c(leg.lty,lty_sa)
                
              }else{
                ts.sa <- log(ts.sa <- object@d11)
                main<-"Log transformed Seasonally Adjusted Series"
                leg.txt<- c(leg.txt,"Seasonally Adjusted")
                leg.col <- c(leg.col,col_sa)
                leg.lty <- c(leg.lty,lty_sa)
                
              }}
            if(trend){
              if(!log_transform){
                ts.trend <- object@d12	
                main<-"Trend"	
                leg.txt <- c(leg.txt,"Trend")
                leg.col <- c(leg.col,col_trend)
                leg.lty <- c(leg.lty,lty_trend)
                
              }else{
                ts.trend <- log(ts.trend <- object@d12)
                main<-"Log transformed Trend"
                leg.txt <- c(leg.txt,"Trend")
                leg.col <- c(leg.col,col_trend)
                leg.lty <- c(leg.lty,lty_trend)
                
              }}
            if(sa && trend  &! original){
              if(!log_transform)
                main <- "Seasonally Adjusted Series and Trend"
              else
                main <- "Log transformed Seasonally Adjusted Series and Trend"
            }
            if(original && sa &!trend)	
              main <- paste(main.orig,"and Seasonally Adjusted Series")
            if(original &! sa &&trend)	
              main <- paste(main.orig,"and Trend")
            if(original && sa && trend)	
              main <- paste(main.orig,", Seasonally Adjusted Series and Trend",sep="")
            if(user.main.logical)
              main<-user.main	
            if(forecast && backcast &! is.na(object@forecast@estimate[1]) &! is.na(object@backcast@estimate[1]))
              main<-"Time Series with Back- and Forecasts"
            if(forecast &! backcast  &! is.na(object@forecast@estimate[1]))
              main<-"Time Series with Forecasts"
            if(!forecast && backcast &! is.na(object@backcast@estimate[1]))
              main<-"Time Series with Backcasts"
            if(forecast && backcast && is.na(object@forecast@estimate[1]))
              main<-"Time Series with Backcasts"
            if(forecast && backcast && is.na(object@backcast@estimate[1]))
              main<-"Time Series with Forecasts"
            
#Falls nur SA/nur Trend geplottet werden soll
            if((sa &!original &!trend) | (sa && trend &!original)){
              ts<-ts.sa
              col_original <- col_sa
              lwd_original <- lwd_sa
              lty_original <- lty_sa
            }
            if(trend &!original &!sa){
              ts <- ts.trend
              col_original <- col_trend
              lwd_original <- lwd_trend
              lty_original <- lty_trend
            }
			
			ts.plot<-ts
#if(sa && trend &!original){
#	ts<-ts.sa
#	col_original <- col_sa
#	lwd_original <- lwd_sa
#}
			if(showAllout && object@dg$outlier=="-" && object@dg$autoout=="-")
			showAllout=FALSE

            if(showAllout | !is.null(showOut)){
              if(showAllout){
                if(any(object@dg$outlier!="-")){
                  names.out <- names(object@dg$outlier)
                  names.out <- tolower(gsub("outlier_","",names.out))}
                if(any(object@dg$autoout!="-")){
                if(!exists("names.out"))
				names.out <- tolower(gsub("autooutlier_","",names(object@dg$autoout)))  
				else	
			  names.out <- tolower(c(names.out,gsub("autooutlier_","",names(object@dg$autoout))))
	  			}}
              if(!is.null(showOut)){
                names.out <- tolower(showOut)	
              }		
              months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
              num.months <- sapply(names.out,function(y) which(unlist(sapply(months,function(x) grepl(x,y,fixed=TRUE)))))
              num.months <- as.numeric(num.months)
              years <- suppressWarnings(as.numeric(sapply(strsplit(names.out,""),function(x)paste(as.vector(na.omit(as.numeric(x)))[1:4],collapse=""))))
              rest.months <- suppressWarnings(as.numeric(sapply(strsplit(names.out,""),function(x){
										  l.mon<-5:ifelse(length(unlist(x))==8, 5, 6)
										  paste(as.vector(na.omit(as.numeric(x)))[l.mon],collapse="")})))
              months <- na.omit(c(rbind(num.months,rest.months)))
              out <- cbind(years,months)
              pch.out <- col.out <- names.out2 <- sapply(1:(dim(out)[1]),function(i) unlist(strsplit(names.out[i],as.character(years)[i]))[1])
#				col.out <- names.out2
              col.out <- replace(col.out,which(col.out%in%"ao"),col_ao)
              col.out <- replace(col.out,which(col.out%in%"ls"),col_ls)
              col.out <- replace(col.out,which(col.out%in%"tc"),col_tc)
              pch.out <- replace(pch.out,which(pch.out%in%"ao"),pch_ao)
              pch.out <- replace(pch.out,which(pch.out%in%"ls"),pch_ls)
              pch.out <- as.numeric(replace(pch.out,which(pch.out%in%"tc"),pch_tc))
            }				
            if(plot_legend){
              par(mar = c(4, 4, 4, 2) + 0.1) 
              layout(matrix(c(rep(1,16),2,2),nrow=9,byrow=TRUE))#,heights = c(1, 6), respect = FALSE)
              if(forecast | backcast){
                if(!original){
                  cat("Fore-/Backcasts are only available for 'original' (log transformed) time series!\n")	
                  if(is.null(xlim)){
                    plot(ts,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
                  }else{
                    plot(ts,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
                  }	
                }else{
                  ts.plot <- plotFbcast(object=object,showCI=showCI,
                      main=main,forecast=forecast,backcast=backcast,log_transform=log_transform,
                      col_original=col_original,col_fc=col_fc,col_bc=col_bc,
                      col_ci=col_ci,col_cishade=col_cishade,points_original=points_original,
                      lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,
                      lwd_original=lwd_original,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
                      ytop=ytop,points_bc=points_bc,points_fc=points_fc,
                      showLine=showLine,col_line=col_line,lty_line=lty_line,
                      ylab=ylab,xlab=xlab,ylim=ylim,xlim=xlim,showWarnings=showWarnings,...)	
                }
              }else	
              if(is.null(xlim)){
                plot(ts,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
              }else{
                plot(ts,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
              }
              if(original && sa)
                lines(ts.sa,col=col_sa,type="l",lwd=lwd_sa,lty=lty_sa)
              if(original && trend)
                lines(ts.trend,col=col_trend,type="l",lwd=lwd_trend,lty=lty_trend)
              if(sa && trend &!original)
                lines(ts.trend,col=col_trend,type="l",lwd=lwd_trend,lty=lty_trend)
			if((!forecast || !backcast) && points_original)
				points(ts.plot,col=col_original,lwd=lwd_original)			
			ts<-ts.plot
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
			
              if(is.null(showOut) &! showAllout){
                par(mar = c(0, 0, 0, 0)) 
                if(forecast | backcast){
					ts.plot <- plotFbcast(object=object,showCI=showCI,
                      main=main,forecast=forecast,backcast=backcast,log_transform=log_transform,
                      col_original=NA,col_fc=NA,col_bc=NA,
                      col_ci=NA,col_cishade=NA,
                      ytop=ytop,
                      col_line=NA,
                      ylim=ylim,xlim=xlim,showWarnings=FALSE,type = "n", axes = FALSE, bty="n", ann = FALSE,...)	
                }else{
                  if(is.null(xlim)){
                    plot(ts,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE, xaxt="n", ...)
                  }else{
                    plot(ts,xlim=xlim,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE, xaxt="n", ...)
                  }
                }
                if(length(leg.txt)>1){
                  legend("center",legend=leg.txt,col=leg.col,lty=leg.lty,bg="white",horiz=legend_horiz,bty=legend_bty)
                }
              }
              
              if(!is.null(showOut)){
                out.type <- toupper(unlist(strsplit(names.out,out[1]))[1])
                out.trend <- out[2]+frequency(ts)*c(0:(floor(length(ts)/ frequency(ts))-1))
                out.trend.x <- time(ts)[out.trend]
                out.trend.y <- ts[out.trend]
                col.out.trend <- rep(col_annComp,length(out.trend))
                col.out.trend[grep(out[1],out.trend.x)]<-col.out
                if(annComp && annCompTrend){
                  lines(out.trend.x,out.trend.y,type="l",lty=1,lwd=lwd_out,col=col_annComp)
                  lapply(1:length(out.trend),function(i)abline(v=out.trend.x[i],col=col.out.trend[i],lwd=lwd_out,lty=3))
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=cex_out,lwd=lwd_out)
                }
                if(!annComp &! annCompTrend){
#points(window(ts, out, out, frequency=frequency(ts)),pch=pch.out,cex=cex_out,lwd=lwd_out,col=col.out)
                  abline(v=out.trend.x[grep(out[1],out.trend.x)],col=col.out,lty=3,lwd=lwd_out)
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=cex_out,lwd=lwd_out)
                }
                if(!annComp && annCompTrend){
                  lines(out.trend.x,out.trend.y,type="l",lty=1,lwd=lwd_out,col=col_annComp)
                  abline(v=out.trend.x[grep(out[1],out.trend.x)],col=col.out,lty=3,lwd=lwd_out)
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=cex_out,lwd=lwd_out)
#	points(window(ts, out, out, frequency=frequency(ts)),pch=pch.out,cex=cex_out,lwd=lwd_out,col=col.out)
                }
                if(annComp &! annCompTrend){
#	lines(out.trend.x,out.trend.y,type="p",lty=1,lwd=lwd_out,col=col.out,pch=pch.out.trend,cex=cex_out)
#lapply(1:(dim(out.othermonths)[1]),function(x)points(window(ts, out.othermonths[x,], out.othermonths[x,], frequency=frequency(ts)),pch=4,cex=cex_out,lwd=lwd_out,col=col.out))
                  lapply(1:length(out.trend),function(i)abline(v=out.trend.x[i],col=col.out.trend[i],lwd=lwd_out,lty=3))
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=cex_out,lwd=lwd_out)
                }	
                par(mar = c(0, 0, 0, 0)) 
                if(forecast | backcast){
					ts.plot <- plotFbcast(object=object,showCI=showCI,
                      main=main,forecast=forecast,backcast=backcast,log_transform=log_transform,
                      col_original=NA,col_fc=NA,col_bc=NA,
                      col_ci=NA,col_cishade=NA,
                      ytop=ytop,
                      col_line=NA,showWarnings=FALSE,
                      ylim=ylim,xlim=xlim,type = "n", axes = FALSE, bty="n", ann = FALSE,...)	
                }else{
                  if(is.null(xlim)){
                    plot(ts,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE,xaxt="n",...)
                  }else{
                    plot(ts,xlim=xlim,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE,xaxt="n",...)
                    
                  }
                }
                if(annComp && annCompTrend){
                  if(original && sa && trend)
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"Annual comparison",leg.txt[3]),col=c(leg.col[1],col.out,leg.col[2],col_annComp,leg.col[3]),bg="white",
                        lty=c(leg.lty[1],3,leg.lty[2],3,leg.lty[3]),pch=c(NA,pch.out,NA,NA,NA),ncol=3,horiz=legend_horiz,bty=legend_bty)	
                  else if((original &! sa &! trend) | (!original && sa &! trend) | (!original && trend &! sa))
                    legend("center",legend=c(leg.txt[1],out.type,"","Annual comparison"),col=c(leg.col[1],col.out,NA,col_annComp),
                        bg="white",lty=c(leg.lty[1],3,NA,3),pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                  else if((original && sa &! trend) | (original &! sa && trend) | (sa && trend &!original))
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"Annual comparison"),
                        col=c(leg.col[1],col.out,leg.col[2],col_annComp),bg="white",
                        lty=c(leg.lty[1],3,leg.lty[2],3),pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                  
                }
                
                if(!annComp && annCompTrend){
                  if(original && sa && trend)
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"Annual comparison",leg.txt[3]),col=c(leg.col[1],col.out,leg.col[2],col_annComp,leg.col[3]),
                      bg="white",lty=c(leg.lty[1],3,leg.lty[2],1,leg.lty[3]),
                      pch=c(NA,pch.out,NA,NA,NA),ncol=3,horiz=legend_horiz,bty=legend_bty)	
                  else if((original &! sa &! trend) | (!original && sa &! trend) | (!original && trend &! sa))
                    legend("center",legend=c(leg.txt[1],out.type,"","Annual comparison"),col=c(leg.col[1],col.out,NA,col_annComp),bg="white",
                        lty=c(leg.lty[1],3,NA,1),pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                  else if((original && sa &! trend) | (original &! sa && trend) | (sa && trend &!original))
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"Annual comparison"),
                        col=c(leg.col[1],col.out,leg.col[2],col_annComp),bg="white",
                        lty=c(leg.lty[1],3,leg.lty[2],1),pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                  
                }
                
                if(annComp &! annCompTrend){
                  if(original && sa && trend)
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"Annual comparison",leg.txt[3]),
                        col=c(leg.col[1],col.out,leg.col[2],col_annComp,leg.col[3]),bg="white",
                        lty=c(leg.lty[1],3,leg.lty[2],3,leg.lty[3]),pch=c(NA,pch.out,NA,NA,NA),ncol=3,
                        horiz=legend_horiz,bty=legend_bty)	
                  else if((original &! sa &! trend) | (!original && sa &! trend) | (!original && trend &! sa))
                    legend("center",legend=c(leg.txt[1],out.type,"","Annual comparison"),
                        col=c(leg.col[1],col.out,NA,col_annComp),bg="white",lty=c(leg.lty[1],3,NA,3),
                        pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                  else if((original && sa &! trend) | (original &! sa && trend) | (sa && trend &!original))
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"Annual comparison"),
                        col=c(leg.col[1],col.out,leg.col[2],col_annComp),bg="white",
                        lty=c(leg.lty[1],3,leg.lty[2],3),pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                }
                
                if(!annComp &! annCompTrend){
                  if(original && sa && trend)
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],"",leg.txt[3]),col=c(leg.col[1],col.out,leg.col[2],NA,leg.col[3]),bg="white",lty=c(leg.lty[1],3,leg.lty[2],NA,leg.lty[3]),
                        pch=c(NA,pch.out,NA,NA,NA),ncol=3,horiz=legend_horiz,bty=legend_bty)	
                  else if((original &! sa &! trend) | (!original && sa &! trend) | (!original && trend &! sa))
                    legend("center",legend=c(leg.txt[1],out.type),col=c(leg.col[1],col.out),
                        bg="white",lty=c(leg.lty[1],3),pch=c(NA,pch.out),horiz=legend_horiz,bty=legend_bty)	
                  else if((original && sa &! trend) | (original &! sa && trend) | (sa && trend &!original))
                    legend("center",legend=c(leg.txt[1],out.type,leg.txt[2],""),col=c(leg.col[1],col.out,leg.col[2],NA),
                        bg="white",lty=c(leg.lty[1],3,leg.lty[2],NA),pch=c(NA,pch.out,NA,NA),ncol=2,horiz=legend_horiz,bty=legend_bty)	
                  
                }
              }
              if(showAllout){
                lapply(1:(dim(out)[1]),function(x)points(window(ts, out[x,], out[x,], frequency=frequency(ts)),cex=cex_out,lwd=lwd_out,col=col.out[x],pch=pch.out[x]))
                out.type.ind<-	unlist(lapply(c("ao","ls","tc"),function(x) any(grepl(x,names.out))))
                out.type<-rep("",3)
                out.type.col<-rep(NA,3)
                out.type.pch<-rep(NA,3)
                if(length(which(out.type.ind))>0){
                  out.type[1:length(which(out.type.ind))] <- toupper(c("ao","ls","tc")[out.type.ind])
                  out.type.col<-c(col_ao,col_ls,col_tc)[out.type.ind]
                  out.type.pch<-c(pch_ao,pch_ls,pch_tc)[out.type.ind]
                }
                if(showAlloutLines)
                  lapply(1:(dim(out)[1]),function(x)abline(v=time(window(ts, out[x,], out[x,], frequency=frequency(ts))),col=col.out[x],lwd=lwd_out,lty=3))
#	par(mar = c(0, 0, 0, 0)) 
#	plot(ts,ylim=c(min(ts,na.rm=TRUE),max(ts,na.rm=TRUE)*ytop),type = "n", axes = FALSE, ann = FALSE)
#	if(original && sa && trend)
#		legend("center",legend=c(leg.txt[1],"AO",leg.txt[2],"LS",leg.txt[3],"TC"),col=c(leg.col[1],col_ao,leg.col[2],col_ls,leg.col[3],col_tc),bg="white",lty=c(1,NA,1,NA,1,NA),pch=c(NA,pch_ao,NA,pch_ls,NA,pch_tc),pt.cex=c(NA,2,NA,2,NA,2),ncol=3)	
#	else if((original &! sa &! trend) | (!original && sa &! trend) | (!original && trend &! sa))
#		legend("center",legend=c("AO","LS","TC"),col=c(col_ao,col_ls,col_tc),bg="white",pch=c(pch_ao,pch_ls,pch_tc),pt.cex=2,horiz=legend_horiz)	
#	else if((original && sa &! trend) | (original &! sa && trend) | (sa && trend &!original))
#		legend("center",legend=c(leg.txt[1],"AO",leg.txt[2],"LS","","TC"),col=c(leg.col[1],col_ao,leg.col[2],col_ls,NA,col_tc),bg="white",lty=c(1,NA,1,NA,NA,NA),pch=c(NA,pch_ao,NA,pch_ls,NA,pch_tc),pt.cex=c(NA,2,NA,2,NA,2),ncol=3)	
                
                par(mar = c(0, 0, 0, 0)) 
                if(forecast | backcast){
					ts.plot <- plotFbcast(object=object,showCI=showCI,
                      main=main,forecast=forecast,backcast=backcast,log_transform=log_transform,
                      col_original=NA,col_fc=NA,col_bc=NA,
                      col_ci=NA,col_cishade=NA,
                      ytop=ytop,
                      col_line=NA,showWarnings=FALSE,
                      ylim=ylim,xlim=xlim,type = "n", axes = FALSE, bty="n", ann = FALSE,...)	
                }else{
                  if(is.null(xlim)){	
                    plot(ts,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE,xaxt="n",...)
                  }else{
                    plot(ts,xlim=xlim,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE,xaxt="n",...)
                  }
                }
                if(original && sa && trend)
                  legend("center",legend=c(leg.txt[1],out.type[1],leg.txt[2],out.type[2],leg.txt[3],out.type[3]),
                      col=c(leg.col[1],out.type.col[1],leg.col[2],out.type.col[2],leg.col[3],out.type.col[3]),
                      bg="white",lty=c(leg.lty[1],NA,leg.lty[2],NA,leg.lty[3],NA),pch=c(NA,out.type.pch[1],NA,out.type.pch[2],NA,out.type.pch[3]),pt.cex=c(NA,2,NA,2,NA,2),
                      ncol=3,horiz=legend_horiz,bty=legend_bty)	
                else if((original &! sa &! trend) | (!original && sa &! trend) | (!original && trend &! sa))
                  legend("center",legend=c(out.type[1],out.type[2],out.type[3]),
                      col=c(out.type.col[1],out.type.col[2],out.type.col[3]),bg="white",
                      pch=c(out.type.pch[1],out.type.pch[2],out.type.pch[3]),
                      pt.cex=2,horiz=legend_horiz,bty=legend_bty)	
                else if((original && sa &! trend) | (original &! sa && trend) | (sa && trend &!original))
                  legend("center",legend=c(leg.txt[1],out.type[1],leg.txt[2],out.type[2],"",out.type[3]),
                      col=c(leg.col[1],out.type.col[1],leg.col[2],out.type.col[2],NA,out.type.col[3]),bg="white",lty=c(leg.lty[1],NA,leg.lty[2],NA,NA,NA),pch=c(NA,out.type.pch[1],NA,out.type.pch[2],NA,out.type.pch[3]),pt.cex=c(NA,2,NA,2,NA,2),
                      ncol=3,horiz=legend_horiz,bty=legend_bty)	
              }
#end if plot legend						
            }else{
              if(forecast | backcast){
                if(!original){
                  cat("Fore-/Backcasts are only available for 'original' (log transformed) time series!\n")
                  if(is.null(xlim)){
                    plot(ts,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
                  }else{
                    plot(ts,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
                  }
                }else{
					ts.plot <- plotFbcast(object=object,showCI=showCI,
                      main=main,forecast=forecast,backcast=backcast,log_transform=log_transform,
                      col_original=col_original,col_fc=col_fc,col_bc=col_bc,
                      col_ci=col_ci,col_cishade=col_cishade,points_original=points_original,
                      lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,
                      lwd_original=lwd_original,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
                      ytop=ytop,points_bc=points_bc,points_fc=points_fc,
                      showLine=showLine,col_line=col_line,lty_line=lty_line,
                      ylab=ylab,xlab=xlab,ylim=ylim,xlim=xlim,...)	
                }
              }else	
              if(is.null(xlim)){
                plot(ts,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
              }else{
                plot(ts,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
              }
              if(original && sa)
                lines(ts.sa,col=col_sa,type="l",lwd=lwd_sa,lty=lty_sa)
              if(original && trend)
                lines(ts.trend,col=col_trend,type="l",lwd=lwd_trend,lty=lty_trend)
              if(sa && trend &!original)
                lines(ts.trend,col=col_trend,type="l",lwd=lwd_trend,lty=lty_trend)
              
              if(!is.null(showOut)){	
                out.trend <- out[2]+frequency(ts)*c(0:(floor(length(ts)/ frequency(ts))-1))
                out.trend.x <- time(ts)[out.trend]
                out.trend.y <- ts[out.trend]
                col.out.trend <- rep(col_annComp,length(out.trend))
                col.out.trend[grep(out[1],out.trend.x)]<-col.out
                if(annComp && annCompTrend){
                  lines(out.trend.x,out.trend.y,type="l",lty=1,lwd=lwd_out,col=col_annComp)
                  lapply(1:length(out.trend),function(i)abline(v=out.trend.x[i],col=col.out.trend[i],lwd=lwd_out,lty=3))
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=ceiling(cex_out/2),lwd=lwd_out)
                }
                if(!annComp &! annCompTrend){
                  abline(v=out.trend.x[grep(out[1],out.trend.x)],col=col.out,lty=3,lwd=lwd_out)
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=ceiling(cex_out/2),lwd=lwd_out)
                }
                if(!annComp && annCompTrend){
                  lines(out.trend.x,out.trend.y,type="l",lty=1,lwd=lwd_out,col=col_annComp)
                  abline(v=out.trend.x[grep(out[1],out.trend.x)],col=col.out,lty=3,lwd=lwd_out)
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=ceiling(cex_out/2),lwd=lwd_out)
                }
                if(annComp &! annCompTrend){
                  lapply(1:length(out.trend),function(i)abline(v=out.trend.x[i],col=col.out.trend[i],lwd=lwd_out,lty=3))
                  points(out.trend.x[grep(out[1],out.trend.x)],out.trend.y[grep(out[1],out.trend.x)],col=col.out,pch=pch.out,cex=ceiling(cex_out/2),lwd=lwd_out)
                }
              }else if(showAllout){			
                points.temp<-lapply(1:(dim(out)[1]),function(x)points(window(ts, out[x,], out[x,], frequency=frequency(ts)),cex=ceiling(cex_out/2),lwd=lwd_out,col=col.out[x],pch=pch.out[x]))
                if(showAlloutLines)
                  lines.temp<-lapply(1:(dim(out)[1]),function(x)abline(v=time(window(ts, out[x,], out[x,], frequency=frequency(ts))),col=col.out[x],lwd=lwd_out,lty=3))
                
              }
#else if(is.null(showOut) &! showAllout)	{		
#				
#				#} if showAllout | !is.null(showOut)
#				
#				plot(ts,ylim=c(min(ts,na.rm=TRUE),max(ts,na.rm=TRUE)*ytop),xlab=xlab,ylab=ylab,main=main,lwd=lwd_original,col=col_original)
#				if(original && sa)
#					lines(ts.sa,col=col_sa,type="l",lwd=lwd_sa)
#				if(original && trend)
#					lines(ts.trend,col=col_trend,type="l",lwd=lwd_trend)
#				if(sa && trend &!original)
#					lines(ts.trend,col=col_trend,type="l",lwd=lwd_trend)
#				
#			}
ts<-ts.plot
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


		}
	if((!forecast || !backcast) && points_original && !plot_legend){
			points(ts.plot,col=col_original,lwd=lwd_original)
	ts<-ts.plot
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
	}
	
			gp.new<-par()	
			invisible(gp.new)
      gp <- gp[-which(names(gp)=="page")]
      par(gp)
      
            },finally=par(gp))
			}
)

setMethod(f='plot',
    signature=signature(x = "x12Single"),
    definition=function(x,original=TRUE,sa=FALSE,trend=FALSE,
        log_transform=FALSE,
        ylab="Value",xlab="Date",
        main="TS",
        col_original="black",col_sa="blue",col_trend="green",
        lwd_original=1,lwd_sa=1,lwd_trend=1,lty_sa=1,lty_trend=1,
        ytop=1,showAllout=FALSE,showAlloutLines=FALSE,showOut=NULL,annComp=TRUE,annCompTrend=TRUE,
        col_ao="red",col_ls="red",col_tc="red",col_annComp="grey",lwd_out=1,cex_out=1.5,
        pch_ao=4,pch_ls=2,pch_tc=23,plot_legend=TRUE,legend_horiz=TRUE,legend_bty="o",
        ### implement plotFbcast
        forecast=FALSE,backcast=FALSE,showCI=TRUE,
        col_fc="#2020ff",col_bc="#2020ff",col_ci="#d1d1ff",col_cishade="#d1d1ff",
        lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,
        points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
        showLine=FALSE,col_line="grey",lty_line=3,ylim=NULL,span=NULL,...){
      plot(x@x12Output,original=original,sa=sa,trend=trend,
          log_transform=log_transform,
          ylab=ylab,xlab=xlab,
          main=main,
          col_original=col_original,col_sa=col_sa,col_trend=col_trend,
          lwd_original=lwd_original,lwd_sa=lwd_sa,lwd_trend=lwd_trend,lty_sa=lty_sa,lty_trend=lty_sa,
          ytop=ytop,showAllout=showAllout,showAlloutLines=showAlloutLines,showOut=showOut,annComp=annComp,annCompTrend=annCompTrend,
          col_ao=col_ao,col_ls=col_ls,col_tc=col_tc,col_annComp=col_annComp,lwd_out=lwd_out,cex_out=cex_out,
          pch_ao=pch_ao,pch_ls=pch_ls,pch_tc=pch_tc,plot_legend=plot_legend,
          legend_horiz=legend_horiz,legend_bty=legend_bty,
          ### implement plotFbcast
          forecast=forecast,backcast=backcast,showCI=showCI,
          col_fc=col_fc,col_bc=col_bc,col_ci=col_ci,col_cishade=col_cishade,
          lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
          points_bc=points_bc,points_fc=points_fc,points_original=points_original,
          showLine=showLine,col_line=col_line,lty_line=lty_line,ylim=ylim,span=span,...)      
    }
)

setMethod(f='plot',
    signature=signature(x = "x12Batch"),
    definition=function(x,what="ask",original=TRUE,sa=FALSE,trend=FALSE,
        log_transform=FALSE,
        ylab="Value",xlab="Date",
        main="TS",
        col_original="black",col_sa="blue",col_trend="green",
        lwd_original=1,lwd_sa=1,lwd_trend=1,lty_sa=1,lty_trend=1,
        ytop=1,showAllout=FALSE,showAlloutLines=FALSE,showOut=NULL,annComp=TRUE,annCompTrend=TRUE,
        col_ao="red",col_ls="red",col_tc="red",col_annComp="grey",lwd_out=1,cex_out=1.5,
        pch_ao=4,pch_ls=2,pch_tc=23,plot_legend=TRUE,legend_horiz=TRUE,legend_bty="o",
        ### implement plotFbcast
        forecast=FALSE,backcast=FALSE,showCI=TRUE,
        col_fc="#2020ff",col_bc="#2020ff",col_ci="#d1d1ff",col_cishade="#d1d1ff",
        lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,
        points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
        showLine=FALSE,col_line="grey",lty_line=3,ylim=NULL,span=NULL,...){
      n <- length(x@x12List)
         if(what=="ask"){
           ask=par("ask")
           par(ask=TRUE)
      for(i in 1:n){
        mainB <- paste(xb@x12List[[i]]@tsName,main,sep="-")
        plot(x@x12List[[i]]@x12Output,original=original,sa=sa,trend=trend,
          log_transform=log_transform,
          ylab=ylab,xlab=xlab,
          main=mainB,
          col_original=col_original,col_sa=col_sa,col_trend=col_trend,
          lwd_original=lwd_original,lwd_sa=lwd_sa,lwd_trend=lwd_trend,lty_sa=lty_sa,lty_trend=lty_sa,
          ytop=ytop,showAllout=showAllout,showAlloutLines=showAlloutLines,showOut=showOut,annComp=annComp,annCompTrend=annCompTrend,
          col_ao=col_ao,col_ls=col_ls,col_tc=col_tc,col_annComp=col_annComp,lwd_out=lwd_out,cex_out=cex_out,
          pch_ao=pch_ao,pch_ls=pch_ls,pch_tc=pch_tc,plot_legend=plot_legend,
          legend_horiz=legend_horiz,legend_bty=legend_bty,
          ### implement plotFbcast
          forecast=forecast,backcast=backcast,showCI=showCI,
          col_fc=col_fc,col_bc=col_bc,col_ci=col_ci,col_cishade=col_cishade,
          lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
          points_bc=points_bc,points_fc=points_fc,points_original=points_original,
          showLine=showLine,col_line=col_line,lty_line=lty_line,ylim=ylim,span=span,...)
      }
    }
    par(ask=ask)
    }
)

setGeneric("plotRsdAcf",
    function(x, ...) { standardGeneric("plotRsdAcf")} )


setMethod(
    f='plotRsdAcf',
    signature=signature(x = "x12Output"),definition=function(x,which="acf",
        xlab="Lag",ylab="ACF",
        main="default",col_acf="darkgrey",lwd_acf=4,
        col_ci="blue",lt_ci=2,ylim="default",...)
    {
      x <-x@dg
      if(which=="acf")
        which <- "rsd.acf"
      else if(which=="pacf")
        which <- "rsd.pacf"
      else if(which=="acf2")
        which <- "rsd.acf2"		
      #lwd_bar=4,plot_legend=TRUE){
      if(main=="default"){
        if(which=="rsd.acf"){main <- "Autocorrelations of the Residuals"}
        else if(which=="rsd.pacf"){main <- "Partial Autocorrelations of the Residuals"}        
        else if(which=="rsd.acf2"){main <- "Autocorrelations of the Squared Residuals"}
      }
	  if(!is.null(x[[which]])){
      if(ylim=="default"){
        ylim<-c(-max(abs(x[[which]][[grep("sample",names(x[[which]]),value=TRUE)]]),2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]]),max(abs(x[[which]][[grep("sample",names(x[[which]]),value=TRUE)]]),2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]]))
      }
      if(which=="rsd.pacf")
        ylab="Partial ACF"
      
      plot(x[[which]]$lag,x[[which]][[grep("sample",names(x[[which]]),value=TRUE)]],type="h",xlab=xlab,ylab=ylab,main=main,col=col_acf,ylim=ylim,lwd=lwd_acf,xaxt="n",...)
      if(length(x[[which]]$lag)%in%c(12,24)){
        aT <- c(6,12,18,24)
        axis(side=1,at=aT)  
      }else{
        aT <- c(4,8,12,16)
        axis(side=1,at=aT)
      }
      abline(h=0,col="black")
      lines(x[[which]]$lag,2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]],type="l",col=col_ci,lty=lt_ci)
      lines(x[[which]]$lag,-2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]],type="l",col=col_ci,lty=lt_ci)
     }
	 else{
	   plot(1:10, type = "n", xaxt="n", yaxt="n", xlab="", ylab="", main=main)	 
	   text(5.5,5.5,"Not Available")
 	 }
		 
 }
)

setMethod(f='plotRsdAcf',
    signature=signature(x = "x12Single"),
    definition=function(x,which="acf",
        xlab="Lag",ylab="ACF",
        main="default",col_acf="darkgrey",lwd_acf=4,
        col_ci="blue",lt_ci=2,ylim="default",...){
      plotRsdAcf(x@x12Output,which=which,
          xlab=xlab,ylab=ylab,
          main=main,col_acf=col_acf,lwd_acf=lwd_acf,
          col_ci=col_ci,lt_ci=lt_ci,ylim=ylim,...)      
    }
)

setGeneric("plotSeasFac",
    function(x, ...) { standardGeneric("plotSeasFac")} )


setMethod(
    f='plotSeasFac',
    signature=signature(x = "x12Output"),definition=function(x,SI_Ratios=TRUE,ylab="Value",xlab="",
        lwd_seasonal=1,col_seasonal="black",lwd_mean=1,col_mean="blue",col_siratio="darkgreen",
        col_replaced="red",cex_siratio=.9,cex_replaced=.9,SI_Ratios_replaced=TRUE,plot_legend=TRUE,
        legend_horiz=FALSE,legend_bty="o",
        ...)
    {
      
      if(!SI_Ratios)
        v <- as.vector(x@d10) # Seasonal Factors
      else
        v <- as.vector(x@d10)[1:length(x@d8)] # Seasonal Factors without forecast
      f <- frequency(x@d10)
      dif <- length(v)%%f
      if(dif>0)
        v[length(v)+(1:(f-dif))]<-NA
      out_matrix <- matrix(v,ncol=f,byrow=TRUE)
      if(f==12){
        lab <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      }else if(f==4){
        lab <- c("Qtr1","Qtr2","Qtr3","Qtr4")
      }else if(f==2){
        lab <- c("1st Half","2nd Half")
      }else{
        lab <- 1:f
      }
      if("main"%in%names(list(...))){
        main <- list(...)[["main"]]  
      }else{
        if(SI_Ratios){
          main="Seasonal Factors by period and SI Ratios"
        }else{
          main="Seasonal Factors by period"
        }
      }
      #ylim <- c(min(v,na.rm=TRUE)*.95,max(v,na.rm=TRUE)*1.09)
      if(!"ylim"%in%names(list(...)))
        ylim <- c(min(v,na.rm=TRUE)*.99,max(v,na.rm=TRUE)*1.01)
      xlim <- c(0,f)
      gp<-par()
      for(i in c("cin","cra","csi","cxy","din")){
        gp <- gp[-which(names(gp)%in%i)]			
      }
      tryCatch({
            if(plot_legend){
              par(mar = c(4, 4, 4, 2) + 0.1) 
              layout(matrix(c(rep(1,16),2,2),nrow=9,byrow=TRUE))
              cex_siratio<-cex_siratio*1.5
              cex_replaced <-cex_replaced*1.5
            }
            if(!"ylim"%in%names(list(...))){
              if(!"main"%in%names(list(...)))
                plot(1,type="n",main=main,xlim=xlim,ylim=ylim,xaxt="n",ylab=ylab,xlab=xlab,cex=cex_siratio,...)
              else
                plot(1,type="n",xlim=xlim,ylim=ylim,xaxt="n",ylab=ylab,xlab=xlab,cex=cex_siratio,...)
            }else{
              if(!"main"%in%names(list(...)))
                plot(1,type="n",main=main,xlim=xlim,xaxt="n",ylab=ylab,xlab=xlab,cex=cex_siratio,...)
              else
                plot(1,type="n",xlim=xlim,xaxt="n",ylab=ylab,xlab=xlab,cex=cex_siratio,...)
            }
            axis(1,at=(1:f)-1/2,labels=lab)
            for(i in 0:(f)){    
              abline(v=i,col="grey")
            }
            if(SI_Ratios){
              vv <- as.vector(x@d8) #final unmodified SI Ratios
              dif <- length(vv)%%f
              if(dif>0)
                vv[length(vv)+(1:(f-dif))]<-NA
              out_matrix2 <- matrix(vv,ncol=f,byrow=TRUE)
              vvv <- as.vector(x@d9) # final replacement for SI Ratios
              dif <- length(vvv)%%f
              if(dif>0)
                vvv[length(vvv)+(1:(f-dif))]<-NA
              out_matrix3 <- matrix(vvv,ncol=f,byrow=TRUE)
            }
            for(i in 0:(f-1)){
              s <- seq(.1+i,(i+1)-.1,l=nrow(out_matrix))
              m <- mean(out_matrix[,i+1],na.rm=TRUE)
              points(rep(m,2)~c(s[1],s[length(s)]),type="l",col=col_mean,lwd=lwd_mean)
              points(out_matrix[,i+1]~s,type="l",col=col_seasonal,lwd=lwd_seasonal)
              if(SI_Ratios){
                points(out_matrix2[,i+1]~s,pch=20,cex=cex_siratio,col=col_siratio)
                if(SI_Ratios_replaced)
                  points(out_matrix3[,i+1]~s,pch=20,cex=cex_replaced,col=col_replaced)
              }
            }
            if(plot_legend){
              par(mar = c(0, 0, 0, 0)) 
              if(!"ylim"%in%names(list(...)))
                plot(1,xlim=xlim,ylim=ylim,type = "n", axes = FALSE, bty="n", ann = FALSE, ...)				
              else
                plot(1,xlim=xlim,type = "n", axes = FALSE, bty="n", ann = FALSE, ...)
              if(SI_Ratios){
                if(SI_Ratios_replaced)
                  legend("center",legend=c("Seasonal Factors","Mean","SI Ratio","Replaced SI Ratio"),col=c(col_seasonal,col_mean,col_siratio,col_replaced),pch=c(NA,NA,20,20),
                      lty=c(1,1,NA,NA),bg="white",pt.cex=1.4,horiz=legend_horiz,bty=legend_bty)
                else
                  legend("center",legend=c("Seasonal Factors","Mean","SI Ratio"),
                      col=c(col_seasonal,col_mean,col_siratio),pch=c(NA,NA,20),
                      lty=c(1,1,NA),bg="white",pt.cex=1.4,horiz=legend_horiz,bty=legend_bty)      
              }else
                legend("center",legend=c("Seasonal Factors","Mean"),col=c(col_seasonal,col_mean),
                    lty=c(1,1),bg="white",horiz=legend_horiz,bty=legend_bty)
            }
			gp.new<-par()	
			invisible(gp.new)
      gp <- gp[-which(names(gp)=="page")]
            par(gp)},finally=par(gp))	
	}
)

setMethod(f='plotSeasFac',
    signature=signature(x = "x12Single"),
    definition=function(x,SI_Ratios=TRUE,ylab="Value",xlab="",
        lwd_seasonal=1,col_seasonal="black",lwd_mean=1,col_mean="blue",col_siratio="darkgreen",
        col_replaced="red",cex_siratio=.9,cex_replaced=.9,SI_Ratios_replaced=TRUE,plot_legend=TRUE,legend_horiz=FALSE,legend_bty="o",...){
      plotSeasFac(x@x12Output,SI_Ratios=SI_Ratios,ylab=ylab,xlab=xlab,
          lwd_seasonal=lwd_seasonal,col_seasonal=col_seasonal,lwd_mean=lwd_mean,col_mean=col_mean,
          col_siratio=col_siratio,col_replaced=col_replaced,cex_siratio=cex_siratio,cex_replaced=cex_replaced,
          SI_Ratios_replaced=SI_Ratios_replaced,plot_legend=plot_legend,legend_horiz=legend_horiz,legend_bty=legend_bty,...)      
    }
)

setGeneric("plotSpec",
    function(x, ...) { standardGeneric("plotSpec")} )


setMethod(
    f='plotSpec',
    signature=signature(x = "x12Output"),definition=function(x,which="sa",
        xlab="Frequency",ylab="Decibels",
        main="Spectrum",highlight=TRUE,
        col_bar="darkgrey",col_seasonal="red",col_td="blue",
        lwd_bar=4,lwd_seasonal=4,lwd_td=4,plot_legend=TRUE,
        legend_horiz=TRUE,legend_bty="o",        
        ...)
    {
      if(main=="Spectrum"){
        if(which=="sa")
          main <- "Spectrum of the Seasonally Adjusted Series"
        else if(which=="original")
          main <- "Spectrum of the Original Series"      
        else if(which=="irregular")
          main <- "Spectrum of the Irregular"
        else if(which=="residuals")
          main <- "Spectrum of the RegARIMA Residuals"   
      }
      f<-frequency(x@a1)
      if(which=="sa"){
        which <- "sp1"
        x <- slot(x,which)
      }else if(which=="original"){
        which <- "sp0"
        x <- slot(x,which)
      }else if(which=="irregular"){
        which <- "sp2"
        x <- slot(x,which)
      }else if(which=="residuals"){
        which <- "spr"
        x <- slot(x,which)
      }	
      #out[[which]]$frequency
      plot_spectrum_work(x@frequency,x@spectrum,xlab=xlab,ylab=ylab,f=f,
          main=main,highlight=highlight,
          col_bar=col_bar,col_seasonal=col_seasonal,col_td=col_td,
          lwd_bar=lwd_bar,lwd_seasonal=lwd_seasonal,lwd_td=lwd_td,plot_legend=plot_legend,
          legend_horiz=legend_horiz,legend_bty=legend_bty,
          ...)
      
    })

setMethod(f='plotSpec',
    signature=signature(x = "x12Single"),
    definition=function(x,which="sa",
        xlab="Frequency",ylab="Decibels",
        main="Spectrum",highlight=TRUE,
        col_bar="darkgrey",col_seasonal="red",col_td="blue",
        lwd_bar=4,lwd_seasonal=4,lwd_td=4,plot_legend=TRUE,
        legend_horiz=TRUE,legend_bty="o",  
        ...){
      plotSpec(x@x12Output,which=which,
          xlab=xlab,ylab=ylab,
          main=main,highlight=highlight,
          col_bar=col_bar,col_seasonal=col_seasonal,col_td=col_td,
          lwd_bar=lwd_bar,lwd_seasonal=lwd_seasonal,lwd_td=lwd_td,plot_legend=plot_legend,
          legend_horiz=legend_horiz,legend_bty=legend_bty,
          ...)      
    }
)

setMethod(
    f='plotSpec',
    signature=signature(x = "spectrum"),definition=function(x,frequency,
        xlab="Frequency",ylab="Decibels",
        main="Spectrum",highlight=TRUE,
        col_bar="darkgrey",col_seasonal="red",col_td="blue",
        lwd_bar=4,lwd_seasonal=4,lwd_td=4,plot_legend=TRUE,
        legend_horiz=TRUE,legend_bty="o",
        ...)
    {
#			myfun <- function(x) deparse(substitute(x)) 
#		which=NULL
#			if(main=="Spectrum" && !is.null(which)){
#				if(which=="sa")
#					main <- "Spectrum of the Seasonally Adjusted Series"
#				else if(which=="original")
#					main <- "Spectrum of the Original Series"      
#				else if(which=="irregular")
#					main <- "Spectrum of the Irregular"
#				else if(which=="residuals")
#					main <- "Spectrum of the RegARIMA Residuals"   
#			}
      
      f<-frequency
      #out[[which]]$frequency
      plot_spectrum_work(x@frequency,x@spectrum,xlab=xlab,ylab=ylab,f=f,
          main=main,highlight=highlight,
          col_bar=col_bar,col_seasonal=col_seasonal,col_td=col_td,
          lwd_bar=lwd_bar,lwd_seasonal=lwd_seasonal,lwd_td=lwd_td,plot_legend=plot_legend,
          legend_horiz=legend_horiz,legend_bty=legend_bty,
          ...)
      
    })

setGeneric("plotFbcast",
    function(object, ...) { standardGeneric("plotFbcast")} )


setMethod(
		f='plotFbcast',
		signature=signature(object = "x12Output"),definition=function(object,showCI=TRUE,
				main="Time Series",forecast=TRUE,backcast=TRUE,log_transform=FALSE,
				col_original="black",col_fc="#2020ff",col_bc="#2020ff",
				col_ci="#d1d1ff",col_cishade="#d1d1ff",
				lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,
				lwd_original=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,ytop=1,
				points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
				showLine=FALSE,col_line="grey",lty_line=3,
				ylab="Value",xlab="Date",ylim=NULL,xlim=NULL,showWarnings=TRUE,...)
		{
			
			ts.plot <- ts <- object@a1	
			fc <- object@forecast@estimate
			bc <- object@backcast@estimate
			lci_fc <- object@forecast@lowerci		
			uci_fc <- object@forecast@upperci
			lci_bc <- object@backcast@lowerci
			uci_bc <- object@backcast@upperci
			if(log_transform){
			ts.plot <- ts <- log(ts)
				fc <- log(fc)
				bc <- log(bc)
				lci_fc <- log(lci_fc)
				uci_fc <- log(uci_fc)
				lci_bc <- log(lci_bc)
				uci_bc <- log(uci_bc)
			}
			
			
			
			if((!forecast &! backcast) | (forecast && is.na(fc[1]) &! backcast) | (backcast && is.na(bc[1]) &! forecast) | (forecast && is.na(fc[1]) && backcast && is.na(bc[1]))){
				if(is.null(ylim))
					ylim<-c(min(ts,na.rm=TRUE),max(ts,na.rm=TRUE)*ytop)
				
				plot(ts,xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
				if(showLine)
					abline(v=time(ts)[length(ts)],col=col_line,lty=lty_line)
				if(points_original)
					points(ts,col=col_original,lwd=lwd_original)
				
				if(showWarnings)
					cat("No forecasts or backcasts plotted!\n")
			}				
			if(forecast && is.na(fc[1]) && showWarnings)		
				cat("Warning: No forecasts available for plotting.\n")
			if(backcast && is.na(bc[1]) && showWarnings)
				cat("Warning: No backcasts available for plotting.\n")
			
			if(forecast &! is.na(fc[1]) && (!backcast | is.na(bc[1]))){
				ts.fc<-ts(c(ts[length(ts)],fc),start=end(ts),end=end(fc),frequency=frequency(ts))
				ts.plot<- ts(c(ts,fc),start=start(ts),end=end(fc),frequency=frequency(ts))
				if(main=="Time Series")
					main<-"Time Series with Forecasts"
				if(showCI){
					limits.y<-c(min(ts,ts.fc,lci_fc,uci_fc,na.rm=TRUE),max(ts,ts.fc,lci_fc,uci_fc,na.rm=TRUE)*ytop)
					limits.x<-c(min(time(ts),time(ts.fc),na.rm=TRUE),max(time(ts),time(ts.fc),na.rm=TRUE))
				}else{
					limits.y<-c(min(ts,ts.fc,na.rm=TRUE),max(ts,ts.fc,na.rm=TRUE)*ytop)
					limits.x<-c(min(time(ts),time(ts.fc),na.rm=TRUE),max(time(ts),time(ts.fc),na.rm=TRUE))
				}
#				leg.txt <- c(leg.txt,"Original TS")
#				leg.col <- c(leg.col,col_original)
				if(is.null(ylim)){	
					if(is.null(xlim)){
						plot(ts,ylim=limits.y,xlim=limits.x,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}else{
						plot(ts,xlim=xlim,ylim=limits.y,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}
				}else{
					#ylim<-c(min(limits.y[1],ylim[1]),max(limits.y[2],ylim[2]))	
					if(is.null(xlim)){
						plot(ts,ylim=ylim,xlim=limits.x,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}else{
						plot(ts,xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}
				}
				if(showCI){
					yy <- as.numeric(uci_fc)
					yCI=c(as.numeric(lci_fc),yy[length(time(yy)):1])
					xCI=c(time(lci_fc),time(lci_fc)[length(time(yy)):1])
#		  yCI=c(yCI[length(yCI)],yCI)#,yCI[1])
#		  xCI=c(xCI[length(xCI)],xCI)#,xCI[1])
					polygon(xCI,yCI,col=col_cishade,border=NA)
					lines(lci_fc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					lines(uci_fc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					if(length(lci_fc)==1){
						lines(x=rep(time(ts.fc)[length(ts.fc)],2),y=c(ts.fc[length(ts.fc)],lci_fc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_")	)
						lines(x=rep(time(ts.fc)[length(ts.fc)],2),y=c(ts.fc[length(ts.fc)],uci_fc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_"))	
					}
					
				}
				lines(ts.fc,col=col_fc,lty=lty_fc,lwd=lwd_fc)
				
				if(showLine)
					abline(v=time(ts)[length(ts)],col=col_line,lty=lty_line)
				if(points_fc)
					points(fc,col=col_fc,lwd=lwd_fc)
				if(points_original)
					points(ts,col=col_original,lwd=lwd_original)
				
			}
			if(backcast &! is.na(bc[1]) && (!forecast | is.na(fc[1]))){		
				ts.bc<-ts(c(bc,ts[1]),start=start(bc),end=start(ts),frequency=frequency(ts))
				ts.plot<- ts(c(bc,ts),start=start(bc),end=end(ts),frequency=frequency(ts))
				
				if(main=="Time Series")
					main<-"Time Series with Backcasts"
#				leg.txt <- c(leg.txt,"Original TS")
#				leg.col <- c(leg.col,col_original)
				if(showCI){
					limits.y<-c(min(ts,ts.bc,lci_bc,na.rm=TRUE),max(ts,ts.bc,uci_bc,na.rm=TRUE)*ytop)
					limits.x<-c(min(time(ts),time(ts.bc),na.rm=TRUE),max(time(ts),time(ts.bc),na.rm=TRUE))		
				}else{
					limits.y<-c(min(ts,ts.bc,na.rm=TRUE),max(ts,ts.bc,na.rm=TRUE)*ytop)
					limits.x<-c(min(time(ts),time(ts.bc),na.rm=TRUE),max(time(ts),time(ts.bc),na.rm=TRUE))		
				}
				if(is.null(ylim)){
					if(is.null(xlim)){
						plot(ts,ylim=limits.y,xlim=limits.x,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}else{
						plot(ts,xlim=xlim,ylim=limits.y,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}
				}else{
					#ylim<-c(min(limits.y[1],ylim[1]),max(limits.y[2],ylim[2]))	
					if(is.null(xlim)){
						plot(ts,ylim=ylim,xlim=limits.x,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}else{
						plot(ts,xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}
				}
				if(showLine)
					abline(v=time(ts)[1],col=col_line,lty=lty_line)
				if(showCI){
					yy <- as.numeric(uci_bc)
					yy <- yy[length(yy):1]
					yCI=c(as.numeric(lci_bc),yy)
					xCI=c(time(lci_bc),time(lci_bc)[length(yy):1])
					polygon(xCI,yCI,col=col_cishade,border=NA)
					lines(lci_bc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					lines(uci_bc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					if(length(lci_bc)==1){
						lines(x=rep(time(ts.bc)[length(ts.bc)],2),y=c(ts.bc[length(ts.bc)],lci_bc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_"))
						lines(x=rep(time(ts.bc)[length(ts.bc)],2),y=c(ts.bc[length(ts.bc)],uci_bc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_"))	
					}
					
				}
				lines(ts.bc,col=col_bc,lty=lty_bc,lwd=lwd_bc)
				if(points_bc)
					points(bc,col=col_bc,lwd=lwd_bc)
				if(points_original)
					points(ts,col=col_original,lwd=lwd_original)
				
			}
			
			if(forecast &! is.na(fc[1]) && backcast &! is.na(bc[1])){
				ts.fc<-ts(c(ts[length(ts)],fc),start=end(ts),end=end(fc),frequency=frequency(ts))	
				ts.bc<-ts(c(bc,ts[1]),start=start(bc),end=start(ts),frequency=frequency(ts))
				ts.plot<- ts(c(bc,ts,fc),start=start(bc),end=end(fc),frequency=frequency(ts))
				
				if(main=="Time Series")
					main<-"Time Series with Back- and Forecasts"
#				leg.txt <- c(leg.txt,"Original TS")
#				leg.col <- c(leg.col,col_original)
				if(showCI){
					limits.y<-c(min(ts,ts.fc,ts.bc,lci_bc,lci_bc,lci_fc,uci_fc,na.rm=TRUE),max(ts,ts.fc,ts.bc,lci_bc,lci_bc,lci_fc,uci_fc,na.rm=TRUE)*ytop)
					limits.x<-c(min(time(ts),time(ts.fc),time(ts.bc),na.rm=TRUE),max(time(ts),time(ts.fc),time(ts.bc),na.rm=TRUE))		
				}else{
					limits.y<-c(min(ts,ts.fc,ts.bc,na.rm=TRUE),max(ts,ts.fc,ts.bc,na.rm=TRUE)*ytop)
					limits.x<-c(min(time(ts),time(ts.fc),time(ts.bc),na.rm=TRUE),max(time(ts),time(ts.fc),time(ts.bc),na.rm=TRUE))		
				}
				if(is.null(ylim)){
					if(is.null(xlim)){
						plot(ts,ylim=limits.y,xlim=limits.x,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}else{
						plot(ts,xlim=xlim,ylim=limits.y,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}
				}else{
					#ylim<-c(min(limits.y[1],ylim[1]),max(limits.y[2],ylim[2]))	
					if(is.null(xlim)){
						plot(ts,ylim=ylim,xlim=limits.x,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}else{
						plot(ts,xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,col=col_original,lty=lty_original,xaxt="n",...)
					}
				}
				if(showLine){
					abline(v=time(ts)[length(ts)],col=col_line,lty=lty_line)
					abline(v=time(ts)[1],col=col_line,lty=lty_line)
				}
				if(showCI){
					yy.fc <- as.numeric(uci_fc)
					yy.fc <- yy.fc[length(yy.fc):1]
					yCI.fc=c(as.numeric(lci_fc),yy.fc)
					xCI.fc=c(time(lci_fc),time(lci_fc)[length(yy.fc):1])
					polygon(xCI.fc,yCI.fc,col=col_cishade,border=NA)
					
					yy.bc <- as.numeric(uci_bc)
					yy.bc <- yy.bc[length(yy.bc):1]
					yCI.bc=c(as.numeric(lci_bc),yy.bc)
					xCI.bc=c(time(lci_bc),time(lci_bc)[length(yy.bc):1])
          
					polygon(xCI.bc,yCI.bc,col=col_cishade,border=NA)
					
					lines(lci_fc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					lines(uci_fc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					lines(lci_bc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					lines(uci_bc,col=col_ci,lty=lty_ci,lwd=lwd_ci)
					if(length(lci_bc)==1){
						lines(x=rep(time(ts.bc)[length(ts.bc)],2),y=c(ts.bc[length(ts.bc)],lci_bc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_"))
						lines(x=rep(time(ts.bc)[length(ts.bc)],2),y=c(ts.bc[length(ts.bc)],uci_bc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_"))	
					}
					if(length(lci_fc)==1){
						lines(x=rep(time(ts.fc)[length(ts.fc)],2),y=c(ts.fc[length(ts.fc)],lci_fc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_")	)
						lines(x=rep(time(ts.fc)[length(ts.fc)],2),y=c(ts.fc[length(ts.fc)],uci_fc),col=col_ci,lty=lty_ci,lwd=lwd_ci,type="o",pch=c(NA,"_"))	
					}	
				}
				lines(ts.fc,col=col_fc,lty=lty_fc,lwd=lwd_fc)
				lines(ts.bc,col=col_bc,lty=lty_bc,lwd=lwd_bc)
				
				if(points_fc)
					points(fc,col=col_fc,lwd=lwd_fc)
				if(points_bc)
					points(bc,col=col_bc,lwd=lwd_bc)
				if(points_original)
					points(ts,col=col_original,lwd=lwd_original)
				
			}
			invisible(ts.plot) 
		}			
)

setMethod(f='plotFbcast',
    signature=signature(object = "x12Single"),
    definition=function(object,showCI=TRUE,
        main="Time Series",forecast=TRUE,backcast=TRUE,log_transform=FALSE,
        col_original="black",col_fc="#2020ff",col_bc="#2020ff",
        col_ci="#d1d1ff",col_cishade="#d1d1ff",
        lty_original=1,lty_fc=2,lty_bc=2,lty_ci=1,
        lwd_original=1,lwd_fc=1,lwd_bc=1,lwd_ci=1,ytop=1,
        points_bc=FALSE,points_fc=FALSE,points_original=FALSE,
        showLine=FALSE,col_line="grey",lty_line=3,
        ylab="Value",xlab="Date",ylim=NULL,xlim=NULL,showWarnings=TRUE,...){
      plotFbcast(object@x12Output,showCI=showCI,
          main=main,forecast=forecast,backcast=backcast,log_transform=log_transform,
          col_original=col_original,col_fc=col_fc,col_bc=col_bc,
          col_ci=col_ci,col_cishade=col_cishade,points_original=points_original,
          lty_original=lty_original,lty_fc=lty_fc,lty_bc=lty_bc,lty_ci=lty_ci,
          lwd_original=lwd_original,lwd_fc=lwd_fc,lwd_bc=lwd_bc,lwd_ci=lwd_ci,
          ytop=ytop,points_bc=points_bc,points_fc=points_fc,
          showLine=showLine,col_line=col_line,lty_line=lty_line,
          ylab=ylab,xlab=xlab,ylim=ylim,xlim=xlim,showWarnings=showWarnings,...)	
    }
)

