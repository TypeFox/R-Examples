CIFplot <-
function(	
						x,
						event.code = NULL,
						covar.code = NULL,
						indiv.times = NULL,
						indiv.events = NULL,
						indiv.covar = NULL,
						xlim = c(0, 30),
						ylim = NULL,
						xlab = "Time",
						ylab = "CIF",
						legend = TRUE,
						...
						){
								
    if(class(x) != "irates"){stop("Object needs to be of class irates")}
    object = x
	if(is.null(event.code)) event.code = object$event.code
	if(is.null(covar.code)) covar.code = object$covar.code

  	event.code = as.character(event.code)
  	covar.code = as.character(covar.code)
	
	if(any(!(event.code %in% object$event.code))){
		stop(paste("event.code", paste(event.code[which(!(event.code %in% object$event.code))], collapse = ", "), "is not contained in irates object"))
		}

	ccif = cif(object, t=seq(0, xlim[2], 0.01))
		
	#ltypes = as.data.frame(	matrix(NA, 
	#						nrow=length(covar.code), 
	#						ncol=length(event.code)), 
	#						row.names=covar.code)
	#names(ltypes) = event.code
	
	cols = as.data.frame(	matrix(NA, 
							nrow=length(covar.code), 
							ncol=length(event.code)), 
							row.names=covar.code)
	names(cols) = event.code

	lwds = as.data.frame(	matrix(NA, 
							nrow=length(covar.code), 
							ncol=length(event.code)), 
							row.names=covar.code)
	names(lwds) = event.code

	# estimate non-parametric cif
	if( !is.null(indiv.times) ){
		
		indiv.times = unlist(indiv.times)
		indiv.events = unlist(indiv.events)
  		indiv.covar = unlist(indiv.covar)
		
		indiv.data = as.data.frame(cbind(	id = seq_len(length(indiv.times)),
											from = indiv.covar,
											to = indiv.events,
											times = indiv.times)) 
		names(indiv.data) = c("id", "cov", "event", "time")
		
		origin.events = levels(as.factor(indiv.data$event))
		indiv.data$event = as.factor(indiv.data$event)
		all = levels(indiv.data$event) = paste("ev",levels(indiv.data$event),sep="")
	
		cov = levels(factor(indiv.data$cov))
	
		cens.name = paste("ev", object$no.event.code, sep="")
		events = all[-which(all == cens.name)]
	
		states = c(cov, events)
	
		tra <- matrix(ncol = length(states), nrow = length(states), FALSE)
	
		tra[which(states %in% cov), which(states %in% events)] = TRUE
	    #na <- mvna(tra.data,c("0","1","2","3"),tra2,12)

		names(indiv.data)[which(names(indiv.data) == "event")] = "to"
		names(indiv.data)[which(names(indiv.data) == "cov")] = "from" 
	 	
    	# etm
		tr.prob <- etm(	indiv.data, 
						state.names = states, 
						tra, 
						cens.name = cens.name, 
						s=0
						)

		loop.events = events[which(event.code %in% origin.events)]

		### compute object-cif object
		timep = max(summary(tr.prob)[[paste(covar.code[1], loop.events[1], sep = " ")]]$time)
		
		ccif = cif(object, t=seq(0, timep, 0.1))
		}
		
			
	## for plotting
	vec = rep(NA, length(event.code))
	for(i in 1:length(event.code)){
		vec[i] = max(ccif$cif[[event.code[i]]])
		}
	
	plot(	x=c(0, ccif$t[length(ccif$t)]), 
			y=c(0, ifelse(is.null(ylim), max(vec), ylim[2])), 
			type="n",
			xlab = xlab,
			ylab = ylab,
			...
			)
	
	for(i in 1:length(event.code)){
		for(j in 1:length(covar.code)){
			
			#ltypes[,i] = i
			cols[,i] = gray(0:(length(event.code)-1) / (1*length(event.code)))[i]
			lwds[j,] = j
			lines(	ccif$t, 
					ccif$cif[[event.code[i]]][covar.code[j],], 
					#lty = ltypes[j,i],
					col = cols[j,i],
					lwd = lwds[j,i])
					
			if( !is.null(indiv.times) ){
			lines(	summary(tr.prob)[[paste(covar.code[j], loop.events[i], sep = " ")]]$time, 
			summary(tr.prob)[[paste(covar.code[j], loop.events[i], sep = " ")]]$P,
			#lty = ltypes[j,i],
			col = cols[j,i],
			lwd = lwds[j,i],
			type = "s")
			}

		}
	}
	
	event.lab = object$event.lab[which(object$event.code %in% event.code)]
	covar.lab = object$covar.lab[which(object$covar.code %in% covar.code)]
	
	if(legend){
		legend("topleft", event.lab, lty = rep(1, length(event.code)), lwd=2, col = as.character(cols[1,]))
		legend("bottomright", covar.lab, lwd = as.numeric(lwds[,1]))
		}

	}

