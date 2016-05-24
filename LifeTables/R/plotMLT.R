plotMLT <-
function(..., mx.out, sex, lt.col="All", log=FALSE, age=c(0,1,seq(5,110,5))){
	lt.out <- lt.mx(nmx=mx.out, age=age, sex=sex)$lt
	
	Intervals <- nrow(lt.out)
	Starts <- c(0,1,seq(from=5,to=5*(Intervals-2),by=5))
	n <- c(1,4,rep(5,Intervals-3),500)
	AgeGroups <- rep(as.character(""),Intervals)
	AgeGroups[1] <- as.character(0)
	for(i in 2:(Intervals-1)) {
		AgeGroups[i] <- paste(Starts[i],(Starts[i+1]-1),sep="-")
	}
	AgeGroups[Intervals] <- paste(Starts[Intervals],"+",sep="")
	
	if(log==TRUE){
		y.log <- "y"
		} else {y.log <- ""}
	
	plot.col <- function(to.plot, value){
		col.num <- value+2
		
		if(col.num==3){
			ylab <- expression(""[n]*"m"[x]*"")
			}
		if(col.num==4){
			ylab <- expression(""[n]*"q"[x]*"")
			}
		if(col.num==5){
			ylab <- expression(""[n]*"p"[x]*"")
			}
		if(col.num==6){
			ylab <- expression(""[n]*"d"[x]*"")
			}
		if(col.num==7){
			ylab <- expression(""[]*"l"[x]*"")
			}
		if(col.num==8){
			ylab <- expression(""[n]*"L"[x]*"")
			}
		if(col.num==9){
			ylab <- expression(""[]*"T"[x]*"")
			}
		if(col.num==10){
			ylab <- expression(""[]*"e"[0]*"")
			}
		
		
		plot(to.plot[,col.num], type="l", las=1, ylab=ylab, xlab="", 			xaxt="n", log=y.log)
		axis(1, at=1:length(AgeGroups), labels=AgeGroups, las=2, cex.axis=.85)
		title(xlab="Age (years)", mgp=c(3.5,1,1))
		}
	
	if(lt.col=="All"){
	par(ask=TRUE)
	for(i in 1:8){
	plot.col(to.plot=lt.out, value=i)
		}
	} else {
		par(...)
		plot.col(to.plot=lt.out, value=lt.col)
		}
	}

plot.LifeTable <- function(x, lt.col="nmx", log=TRUE, ...) {
	life.table <- x
	lt.out <- life.table$lt
	Intervals <- nrow(lt.out)
	Starts <- c(0,1,seq(from=5,to=5*(Intervals-2),by=5))
	n <- c(1,4,rep(5,Intervals-3),500)
	AgeGroups <- rep(as.character(""),Intervals)
	AgeGroups[1] <- as.character(0)
	for(i in 2:(Intervals-1)) {
		AgeGroups[i] <- paste(Starts[i],(Starts[i+1]-1),sep="-")
	}
	AgeGroups[Intervals] <- paste(Starts[Intervals],"+",sep="")
	pars <- list(...)
	if (!is.element('oma', names(pars))) pars$oma <- c(0.7,0.7,0,0)
	if (!is.element('cex.axis', names(pars))) pars$cex.axis <- 0.85
	if (!is.element('mgp', names(pars))) pars$mgp <- c(3,0.5,0) # axis labels closer to the axis
	if (!is.element('tcl', names(pars))) pars$tcl <- -0.3 # shorter tick marks
	if(log==TRUE){
		y.log <- "y"
		} else {y.log <- ""}
	l <- length(AgeGroups)
	ylabs <- list(
					nmx=expression(""[n]*"m"[x]*""), nqx=expression(""[n]*"q"[x]*""), 
					npx=expression(""[n]*"p"[x]*""), ndx=expression(""[n]*"d"[x]*""),
					lx=expression(""[]*"l"[x]*""), nLx=expression(""[n]*"L"[x]*""),
					Tx=expression(""[]*"T"[x]*""), ex=expression(""[]*"e"[x]*""))
					
	plot.col <- function(col.name, ...){
		#par(oma=c(0.7,0.7,0,0))
		do.call('plot', c(list(lt.out[,col.name], type="b", las=1, ylab=ylabs[[col.name]], xlab="", 			xaxt="n", log=y.log), pars))
		grid(l, l)
		axis(1, at=1:length(AgeGroups), labels=AgeGroups, las=2, cex.axis=pars$cex.axis)
		title(xlab="Age (years)")
	}
	
	ask <- FALSE
	if(lt.col=="All"){
		ask <- TRUE
		lt.col <- setdiff(colnames(lt.out), c('Age', 'nax'))
	}
	par(ask=ask)
	for(col in lt.col) plot.col(col, ...)
}

