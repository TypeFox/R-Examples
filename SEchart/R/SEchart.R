SEchart <- function(

	data,								
	idvar,								
	timevar.start,							
	timevar.stop,							
	statevar = NULL,							
	eventvar = NULL,							
	eventtime= NULL,							
	srt1 = NULL,							
	srt2 = NULL,							
	srt3 = NULL,							
	srt4 = NULL,							
	ord.ud = NULL,							
	stratvar = NULL,							
	xlim = NULL,							
	xlab = "Time",							
	ylim = NULL,							
	ylab = NULL,							
	grd = FALSE,							
	grd.col = "grey",							
	grd.lty = "dashed",						
	idnum = 1,								
	idnum.col = "black",						
	idnum.cex = 1,							
	l.base.col = "grey",						
	l.base.lwd = 1,							
	lend = 1,								
	l.state.col = "pal.heat.rev",					
	l.state.lwd = 3,							
	main = "SEchart",						
	bg = "white",							
	legnd = TRUE,							
	legnd.x = "bottomright",					
	legnd.y = NULL,							
	legnd.cex = 1,							
	legnd.txt.round = 2,						
	legnd.bg = "white",						
	strat.cex = 1,							
	strat.txt = NULL,							
	strat.bg.col = "lightgrey",					
	pch = NULL,								
	p.col = NULL,							
	p.cex = 1,								
	p.lwd = 1,								
	return.output = FALSE){					

			if(!is.numeric(data[,idvar])){stop("'idvar' is not numeric")}
			if(!is.numeric(data[,timevar.start])){stop("'timevar.start' is not numeric")}
			if(!is.numeric(data[,timevar.stop])){stop("'timevar.stop' is not numeric")}
			if(!is.null(statevar) & !is.numeric(data[,statevar])){stop("'statevar' is not numeric")}
			if(!is.null(stratvar) & !is.numeric(data[,stratvar])){stop("'stratvar' is not numeric")}
			if(!is.null(ord.ud)){if(any(sort(ord.ud) != sort(unique(data[,idvar])))){stop("incorrect 'ord.ud' specification")}}
			if(!is.null(stratvar)){if(any(is.na(data[,stratvar]))){stop("stratvar contains NAs")}}
			if(!is.null(strat.txt)){if(!is.null(stratvar)){if(length(unique(data[,stratvar])) != length(strat.txt)){stop("length of 'strat.txt' does not correspond to actual number of categories in 'stratvar'")}}}	
			if(!is.null(eventvar)){if(is.null(eventtime)){stop("'eventtime' is not specified")}}
			if(!is.null(statevar) & !is.null(l.state.col)){if(length(unique(statevar)) != length(unique(l.state.col))){stop("length of 'l.state.col' does not correspond to number of unique states")}}

			Z <- data.frame(lapply(data[c(idvar, timevar.start, timevar.stop, statevar, eventvar, stratvar)], as.numeric), stringsAsFactors=FALSE)

			if(!is.null(statevar) & !is.null(eventvar) & !is.null(stratvar)) 	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", "statevar", eventvar, "stratvar")
			if(!is.null(statevar) & !is.null(eventvar) &  is.null(stratvar)) 	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", "statevar", eventvar)
			if(!is.null(statevar) &  is.null(eventvar) & !is.null(stratvar))	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", "statevar", "stratvar")
			if( is.null(statevar) & !is.null(eventvar) & !is.null(stratvar))	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", eventvar, "stratvar")
			if( is.null(statevar) &  is.null(eventvar) & !is.null(stratvar)) 	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", "stratvar")
			if( is.null(statevar) & !is.null(eventvar) &  is.null(stratvar)) 	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", eventvar)
			if(!is.null(statevar) &  is.null(eventvar) &  is.null(stratvar)) 	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop", "statevar")
			if( is.null(statevar) &  is.null(eventvar) &  is.null(stratvar)) 	colnames(Z) <- c("idvar", "timevar.start", "timevar.stop")

			Z <- Z[order(Z$idvar),]	


			srt1a <- srt2a <- srt3a <- srt4a <- NULL

			if(!is.null(srt1)){
				if(srt1 == "start.time") 					srt1a <- mapply(unique(Z$idvar), FUN=function(i) min(Z$timevar.start[Z$idvar==i]))
				if(srt1 == "end.time") 						srt1a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i]))
				if(srt1 == "tot.time") 						srt1a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i])-min(Z$timevar.start[Z$idvar==i]))
				if(srt1 == "midpoint.time")					srt1a <- mapply(unique(Z$idvar), FUN=function(i) (max(Z$timevar.stop[Z$idvar==i])+min(Z$timevar.start[Z$idvar==i]))/2)
				if(!is.null(statevar) & srt1 == "min.state") 		srt1a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) min(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt1 == "max.state") 		srt1a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) max(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt1 == "first.state") 		srt1a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][1]))
				if(!is.null(statevar) & srt1 == "last.state") 		srt1a <- as.numeric(suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][length(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])])))
				if(!is.null(statevar) & srt1 == "average.state") 	srt1a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) sum(Z$statevar[Z$idvar==i & !is.na(Z$statevar)]*(Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T) / sum((Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T)))
				if(any(!is.null(eventvar))) {
					for(w in 1:length(eventvar)) {
						if(srt1 == paste("sum.", eventvar[w], sep="")) 	srt1a <- mapply(unique(Z$idvar), FUN=function(i) sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T))
						if(srt1 == paste("tf.", eventvar[w], sep="")) 	srt1a <- mapply(unique(Z$idvar), FUN=function(i) ifelse(sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T)>0, 1, 0))
						if(srt1 == paste("time.", eventvar[w], sep="")) {
							if(eventtime[w] == "start") {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.start[Z[,eventvar[w]]==1]
								srt1a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "middle")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- (Z$timevar.start[Z[,eventvar[w]]==1]+Z$timevar.stop[Z[,eventvar[w]]==1])/2
								srt1a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "end")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.stop[Z[,eventvar[w]]==1]
								srt1a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}}
								}}
								}


			if(!is.null(srt2)){
				if(srt2 == "start.time") 					srt2a <- mapply(unique(Z$idvar), FUN=function(i) min(Z$timevar.start[Z$idvar==i]))
				if(srt2 == "end.time") 						srt2a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i]))
				if(srt2 == "tot.time") 						srt2a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i])-min(Z$timevar.start[Z$idvar==i]))
				if(srt2 == "midpoint.time")					srt2a <- mapply(unique(Z$idvar), FUN=function(i) (max(Z$timevar.stop[Z$idvar==i])+min(Z$timevar.start[Z$idvar==i]))/2)
				if(!is.null(statevar) & srt2 == "min.state") 		srt2a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) min(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt2 == "max.state") 		srt2a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) max(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt2 == "first.state") 		srt2a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][1]))
				if(!is.null(statevar) & srt2 == "last.state") 		srt2a <- as.numeric(suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][length(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])])))
				if(!is.null(statevar) & srt2 == "average.state") 	srt2a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) sum(Z$statevar[Z$idvar==i & !is.na(Z$statevar)]*(Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T) / sum((Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T)))
				if(any(!is.null(eventvar))) {
					for(w in 1:length(eventvar)) {
						if(srt2 == paste("sum.", eventvar[w], sep="")) 	srt2a <- mapply(unique(Z$idvar), FUN=function(i) sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T))
						if(srt2 == paste("tf.", eventvar[w], sep="")) 	srt2a <- mapply(unique(Z$idvar), FUN=function(i) ifelse(sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T)>0, 1, 0))
						if(srt2 == paste("time.", eventvar[w], sep="")) {
							if(eventtime[w] == "start") {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.start[Z[,eventvar[w]]==1]
								srt2a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar== x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "middle")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- (Z$timevar.start[Z[,eventvar[w]]==1]+Z$timevar.stop[Z[,eventvar[w]]==1])/2
								srt2a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "end")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.stop[Z[,eventvar[w]]==1]
								srt2a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}}
								}}
								}

			if(!is.null(srt3)){
				if(srt3 == "start.time") 					srt3a <- mapply(unique(Z$idvar), FUN=function(i) min(Z$timevar.start[Z$idvar==i]))
				if(srt3 == "end.time") 						srt3a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i]))
				if(srt3 == "tot.time") 						srt3a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i])-min(Z$timevar.start[Z$idvar==i]))
				if(srt3 == "midpoint.time")					srt3a <- mapply(unique(Z$idvar), FUN=function(i) (max(Z$timevar.stop[Z$idvar==i])+min(Z$timevar.start[Z$idvar==i]))/2)
				if(!is.null(statevar) & srt3 == "min.state") 		srt3a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) min(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt3 == "max.state") 		srt3a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) max(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt3 == "first.state") 		srt3a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][1]))
				if(!is.null(statevar) & srt3 == "last.state") 		srt3a <- as.numeric(suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][length(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])])))
				if(!is.null(statevar) & srt3 == "average.state") 	srt3a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) sum(Z$statevar[Z$idvar==i & !is.na(Z$statevar)]*(Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T) / sum((Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T)))
				if(any(!is.null(eventvar))) {
					for(w in 1:length(eventvar)) {
						if(srt3 == paste("sum.", eventvar[w], sep="")) 	srt3a <- mapply(unique(Z$idvar), FUN=function(i) sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T))
						if(srt3 == paste("tf.", eventvar[w], sep="")) 	srt3a <- mapply(unique(Z$idvar), FUN=function(i) ifelse(sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T)>0, 1, 0))
						if(srt3 == paste("time.", eventvar[w], sep="")) {
							if(eventtime[w] == "start") {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.start[Z[,eventvar[w]]==1]
								srt3a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "middle")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- (Z$timevar.start[Z[,eventvar[w]]==1]+Z$timevar.stop[Z[,eventvar[w]]==1])/2
								srt3a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "end")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.stop[Z[,eventvar[w]]==1]
								srt3a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}}
								}}
								}

			if(!is.null(srt4)){
				if(srt4 == "start.time") 					srt4a <- mapply(unique(Z$idvar), FUN=function(i) min(Z$timevar.start[Z$idvar==i]))
				if(srt4 == "end.time") 						srt4a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i]))
				if(srt4 == "tot.time") 						srt4a <- mapply(unique(Z$idvar), FUN=function(i) max(Z$timevar.stop[Z$idvar==i])-min(Z$timevar.start[Z$idvar==i]))
				if(srt4 == "midpoint.time")					srt4a <- mapply(unique(Z$idvar), FUN=function(i) (max(Z$timevar.stop[Z$idvar==i])+min(Z$timevar.start[Z$idvar==i]))/2)
				if(!is.null(statevar) & srt4 == "min.state") 		srt4a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) min(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt4 == "max.state") 		srt4a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) max(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])))
				if(!is.null(statevar) & srt4 == "first.state") 		srt4a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][1]))
				if(!is.null(statevar) & srt4 == "last.state") 		srt4a <- as.numeric(suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) Z$statevar[Z$idvar==i & !is.na(Z$statevar)][length(Z$statevar[Z$idvar==i & !is.na(Z$statevar)])])))
				if(!is.null(statevar) & srt4 == "average.state") 	srt4a <- suppressWarnings(mapply(unique(Z$idvar), FUN=function(i) sum(Z$statevar[Z$idvar==i & !is.na(Z$statevar)]*(Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T) / sum((Z$timevar.stop[Z$idvar==i & !is.na(Z$statevar)] - Z$timevar.start[Z$idvar==i & !is.na(Z$statevar)]), na.rm=T)))
				if(any(!is.null(eventvar))) {
					for(w in 1:length(eventvar)) {
						if(srt4 == paste("sum.", eventvar[w], sep="")) 	srt4a <- mapply(unique(Z$idvar), FUN=function(i) sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T))
						if(srt4 == paste("tf.", eventvar[w], sep="")) 	srt4a <- mapply(unique(Z$idvar), FUN=function(i) ifelse(sum(Z[,eventvar[w]][Z$idvar==i], na.rm=T)>0, 1, 0))
						if(srt4 == paste("time.", eventvar[w], sep="")) {
							if(eventtime[w] == "start") {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.start[Z[,eventvar[w]]==1]
								srt4a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "middle")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- (Z$timevar.start[Z[,eventvar[w]]==1]+Z$timevar.stop[Z[,eventvar[w]]==1])/2
								srt4a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}
							if(eventtime[w] == "end")  {
								Z$ind.timeev <- NA
								Z$ind.timeev[Z[,eventvar[w]]==1] <- Z$timevar.stop[Z[,eventvar[w]]==1]
								srt4a <- sapply(unique(Z$idvar), FUN = function(x) ifelse(all(is.na(Z$ind.timeev[Z$id == x])), NA, Z$ind.timeev[Z$idvar == x & !is.na(Z$ind.timeev)][1]))
								}}
								}}
								}

			stratvar.ord <- NULL
			if(!is.null(stratvar))	stratvar.ord <- sapply(unique(Z$idvar), FUN=function(i) Z$stratvar[Z$idvar==i][1])
		
			ord.id <- NULL
			if(is.null(stratvar) & is.null(srt1)  & is.null(srt2)  & is.null(srt3) & is.null(srt4) & any(is.null(ord.ud)))																	{ord.id <- as.numeric(unique(Z$idvar))}
			if(!is.null(stratvar) & is.null(srt1)  & is.null(srt2)  & is.null(srt3) & is.null(srt4) & any(is.null(ord.ud)) & !is.null(stratvar.ord))													{ord.id <- as.numeric(unique(Z$idvar))[order(stratvar.ord)]}
			if(is.null(stratvar) & !is.null(srt1) & is.null(srt2)  & is.null(srt3) & is.null(srt4)& any(is.null(ord.ud)) & !is.null(srt1a))														{ord.id <- as.numeric(unique(Z$idvar))[order(srt1a)]}
			if(!is.null(stratvar) & !is.null(srt1) & is.null(srt2)  & is.null(srt3) & is.null(srt4)& any(is.null(ord.ud)) & !is.null(srt1a)& !is.null(stratvar.ord))										{ord.id <- as.numeric(unique(Z$idvar))[order(stratvar.ord, srt1a)]}
			if(is.null(stratvar) & !is.null(srt1) & !is.null(srt2) & is.null(srt3) & is.null(srt4)& any(is.null(ord.ud)) & !is.null(srt1a) & !is.null(srt2a))											{ord.id <- as.numeric(unique(Z$idvar))[order(srt1a, srt2a)]}
			if(!is.null(stratvar) & !is.null(srt1) & !is.null(srt2) & is.null(srt3) & is.null(srt4)& any(is.null(ord.ud)) & !is.null(srt1a) & !is.null(srt2a)& !is.null(stratvar.ord))							{ord.id <- as.numeric(unique(Z$idvar))[order(stratvar.ord, srt1a, srt2a)]}
			if(is.null(stratvar) & !is.null(srt1) & !is.null(srt2) & !is.null(srt3) & is.null(srt4) & any(is.null(ord.ud)) & !is.null(srt1a) & !is.null(srt2a) & !is.null(srt3a))								{ord.id <- as.numeric(unique(Z$idvar))[order(srt1a, srt2a, srt3a)]}
			if(!is.null(stratvar) & !is.null(srt1) & !is.null(srt2) & !is.null(srt3) & is.null(srt4) & any(is.null(ord.ud)) & !is.null(srt1a) & !is.null(srt2a) & !is.null(srt3a)& !is.null(stratvar.ord))				{ord.id <- as.numeric(unique(Z$idvar))[order(stratvar.ord, srt1a, srt2a, srt3a)]}
			if(is.null(stratvar) & !is.null(srt1) & !is.null(srt2) & !is.null(srt3) & !is.null(srt4) & any(is.null(ord.ud)) & !is.null(srt1a)  & !is.null(srt2a) & !is.null(srt3a) & !is.null(srt4a))					{ord.id <- as.numeric(unique(Z$idvar))[order(srt1a, srt2a, srt3a, srt4a)]}
			if(!is.null(stratvar) & !is.null(srt1) & !is.null(srt2) & !is.null(srt3) & !is.null(srt4) & any(is.null(ord.ud)) & !is.null(srt1a)  & !is.null(srt2a) & !is.null(srt3a) & !is.null(srt4a)& !is.null(stratvar.ord))	{ord.id <- as.numeric(unique(Z$idvar))[order(stratvar.ord, srt1a, srt2a, srt3a, srt4a)]}

			if(any(!is.null(ord.ud))) ord.id <- ord.ud

			if(is.null(ord.id)) stop("Unable to determine ordering, invalid srt specification.")

			Z$pos <- sapply(1:length(Z$idvar), FUN=function(x) which(ord.id %in% Z$idvar[x]))

			Z$ycoord <- length(unique(Z$idvar))-Z$pos+1

			if(!is.null(stratvar)) {
				add <- sapply(1:nrow(Z), FUN=function(x) which(rev(sort(unique(Z$stratvar))) == unique(Z$stratvar[x]))) - 1
				Z$ycoord <- Z$ycoord + (add*2)
				}

			if(legnd==FALSE & is.null(xlim)) 											xlim <- c(min(data[timevar.start]), max(data[timevar.stop]))
			if(legnd==TRUE & (legnd.x == "topright" | legnd.x=="bottomright") & is.null(xlim)) 			xlim <- c(min(data[timevar.start]), max(data[timevar.stop]) + ((max(data[timevar.stop])-min(data[timevar.start]))/4))
			if(legnd==TRUE & (legnd.x != "topright" & legnd.x!="bottomright") & is.null(xlim)) 			xlim <- c(min(data[timevar.start]), max(data[timevar.stop]))

			if(is.null(ylim) & !is.null(stratvar)) 				ylim <- c(0, max(unique(Z$ycoord))+2)
			if(is.null(ylim) & is.null(stratvar)) 				ylim <- c(0, max(unique(Z$ycoord))+1)
			if(is.null(ylab))								ylab <- idvar

			plot(0, 0, type = "p", col="white", main=main, ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, las = 1, yaxt = "n", yaxs = "i")
			rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg)
			if(grd==TRUE) {grid(nx = NULL, ny = NA, col=grd.col, lty=grd.lty)}
			if(idnum==2) {axis(2, at=unique(Z$ycoord), labels=unique(Z$idvar), las=1, cex.axis=idnum.cex, tck=F)}
			if(!is.null(stratvar)) {
				hbars <- sapply(sort(unique(Z$stratvar)), FUN=function(x) max(Z$ycoord[Z$stratvar==x])+1)
				strat.txt.y <- hbars+0.5
				strat.txt.x <- rep(mean(xlim), 2)
				rect(xleft=rep(xlim[1]-(diff(xlim)/25), length(strat.txt.y)), ybottom=strat.txt.y-0.5, xright=rep(xlim[2]+(diff(xlim)/25), length(strat.txt.y)), ytop=strat.txt.y+0.5, col=strat.bg.col)
				if(!is.null(strat.txt)) strat.txt <- strat.txt
				if(is.null(strat.txt)) strat.txt <- paste(stratvar,"=", sort(unique(Z$stratvar)))
				text(strat.txt, x=strat.txt.x, y=strat.txt.y, cex=strat.cex)
				}

			for(i in unique(Z$idvar)){
				lines(x=c(min(Z$timevar.start[Z$idvar==i]), max(Z$timevar.stop[Z$idvar==i & !is.na(Z$timevar.stop)])), y=c(unique(Z$ycoord[Z$idvar==i]), unique(Z$ycoord[Z$idvar==i])), col = l.base.col, lwd=l.base.lwd, lend=lend)
				if(idnum==1) {text(x = min(Z$timevar.start[Z$idvar==i & !is.na(Z$timevar.start)]), y = unique(Z$ycoord[Z$idvar==i]), labels = as.character(i), pos = 2, col = idnum.col, cex = idnum.cex)}}

			if(!is.null(statevar)){
				if(any(l.state.col=="pal.heat")) 		l.state.col <- heat.colors(length(unique(data[statevar][!is.na(data[statevar])])))
				if(any(l.state.col=="pal.heat.rev")) 	l.state.col <- rev(heat.colors(length(unique(data[statevar][!is.na(data[statevar])]))))
				if(any(l.state.col=="pal.topo")) 		l.state.col <- topo.colors(length(unique(data[statevar][!is.na(data[statevar])])))
				if(any(l.state.col=="pal.topo.rev")) 	l.state.col <- rev(topo.colors(length(unique(data[statevar][!is.na(data[statevar])]))))
				if(any(l.state.col=="pal.cm"))		l.state.col <- cm.colors(length(unique(data[statevar][!is.na(data[statevar])])))
				if(any(l.state.col=="pal.cm.rev"))		l.state.col <- rev(cm.colors(length(unique(data[statevar][!is.na(data[statevar])]))))
				if(any(l.state.col=="pal.terrain")) 	l.state.col <- terrain.colors(length(unique(data[statevar][!is.na(data[statevar])])))
				if(any(l.state.col=="pal.terrain.rev")) 	l.state.col <- rev(terrain.colors(length(unique(data[statevar][!is.na(data[statevar])]))))
				if(any(l.state.col=="pal.gray"))		l.state.col <- gray.colors(length(unique(data[statevar][!is.na(data[statevar])])))
				if(any(l.state.col=="pal.gray.rev"))	l.state.col <- rev(gray.colors(length(unique(data[statevar][!is.na(data[statevar])]))))
				if(any(l.state.col=="pal.rainbow.rev")) 	l.state.col <- rev(rainbow(length(unique(data[statevar][!is.na(data[statevar])]))))
				if(any(l.state.col=="pal.rainbow"))		l.state.col <- rainbow(length(unique(data[statevar][!is.na(data[statevar])])))
				}

			if(!is.null(statevar)){
				for(p in sort(unique(Z$statevar[!is.na(Z$statevar)]))){ # voor alle mogelijke states:
					for(y in 1:length(Z$idvar[!is.na(Z$statevar) & Z$statevar==p])){	#voor elke state=p
						lines(x=c(Z$timevar.start[!is.na(Z$statevar) & Z$statevar==p][y], Z$timevar.stop[!is.na(Z$statevar) & Z$statevar==p][y]), y=c(Z$ycoord[!is.na(Z$statevar) & Z$statevar==p][y],Z$ycoord[!is.na(Z$statevar) & Z$statevar==p][y]), col = l.state.col[which(sort(unique(Z$statevar[!is.na(Z$statevar)]))==p)], lwd=l.state.lwd, lend=lend)
						}}}

			if(!is.null(eventvar)) if(is.null(pch)) 	pch <- c(1:length(eventvar))
			if(!is.null(eventvar)) if(is.null(p.col)) p.col <- rep("black", times=length(eventvar))

			if(!is.null(eventvar)){
				for(q in 1:length(eventtime)){
					if(eventtime[q] == "start" & !is.null(eventvar)){
						points(x=Z$timevar.start[Z[eventvar[q]]==1], y=Z$ycoord[Z[c(paste(eventvar[q]))]==1], pch=pch[q], col=p.col[q], cex=p.cex, lwd=p.lwd)
						}
					if(eventtime[q] == "middle" & !is.null(eventvar)){
						points(x=((Z$timevar.start[Z[eventvar[q]]==1]+Z$timevar.stop[Z[eventvar[q]]==1])/2), y=Z$ycoord[Z[c(paste(eventvar[q]))]==1], pch=pch[q], col=p.col[q], cex=p.cex, lwd=p.lwd)
						}
					if(eventtime[q] == "end" & !is.null(eventvar)){
						points(x=Z$timevar.stop[Z[eventvar[q]]==1], y=Z$ycoord[Z[c(paste(eventvar[q]))]==1], pch=pch[q], col=p.col[q], cex=p.cex, lwd=p.lwd)
						}
						}}

		output <- NULL

		if(!is.null(statevar) & is.null(eventvar)){
			col.statevar <- cbind(l.state.col, sort(unique(Z$statevar[!is.na(Z$statevar)])))
			colnames(col.statevar) <- c("color", statevar)		
			output <- list(ord.id=ord.id, col.statevar=col.statevar)
			}
		if(is.null(statevar)  & !is.null(eventvar)){
			p.inf <- cbind(eventvar, pch, p.col, p.lwd, p.cex)
			output <- list(ord.id=ord.id, p.inf=p.inf)
			}
		if(!is.null(statevar) & !is.null(eventvar)){
			col.statevar <- cbind(l.state.col, sort(unique(Z$statevar[!is.na(Z$statevar)])))
			colnames(col.statevar) <- c("color", statevar)		
			p.inf <- cbind(eventvar, pch, p.col, p.lwd, p.cex)
			output <- list(ord.id=ord.id, col.statevar=col.statevar, p.inf=p.inf)
			}


		if(legnd==TRUE | legnd=="true") {

			if(!is.null(statevar) & !is.null(eventvar)) {
				leg.txt <- c(eventvar, paste(statevar, "=", as.character(round(as.numeric(output$col.statevar[,2]), legnd.txt.round))))
				leg.lwd <- c(rep(p.lwd, length(eventvar)), rep(l.state.lwd,length(output$col.statevar)))
				leg.lty <- c(rep(NA, length(eventvar)), rep(1, length(leg.txt)))
				leg.col <- c(p.col, output$col.statevar[,1])
				leg.pch <- c(pch, rep(NA, length(output$col.statevar)))
				legend(x=legnd.x, y=legnd.y, leg.txt, lty = leg.lty, lwd=leg.lwd, col=leg.col, pch=leg.pch, cex=legnd.cex, bg=legnd.bg)
				}

			if(is.null(statevar) & !is.null(eventvar)) {
				leg.txt <- eventvar
				leg.lwd <- rep(p.lwd, length(eventvar))
				leg.col <- p.col
				leg.pch <- pch
				legend(x=legnd.x, y=legnd.y, leg.txt, col=leg.col, pch=leg.pch, cex=legnd.cex, bg=legnd.bg)
				}

			if(!is.null(statevar) & is.null(eventvar)) {
				leg.txt <- paste(statevar, "=", as.character(round(as.numeric(output$col.statevar[,2]), legnd.txt.round)))
				leg.lwd <- rep(l.state.lwd,length(output$col.statevar))
				leg.lty <- rep(1, length(leg.txt))
				leg.col <- output$col.statevar[,1]
				legend(x=legnd.x, y=legnd.y, leg.txt, lty = leg.lty, lwd=leg.lwd, col=leg.col, cex=legnd.cex, bg=legnd.bg)
				}

				}

			if(return.output == TRUE & (is.null(statevar) | !is.null(eventvar))) {return(output)}


			}













