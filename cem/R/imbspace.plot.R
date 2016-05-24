plot.imbalance.space <- function(...) imbspace.plot(...)

imbspace.plot2 <- function(obj,group="1"){
	if(class(obj) != "imbalance.space")
	stop("obj must be of class `imbalance.space'")

	g <- sprintf("G%s",group)
	n <-  obj$space[[g]]
	ML1 <- obj$space$ML1
	Relaxed <- obj$space$Relaxed
	class <- rep("relax", length(n))
	id.raw <- which(Relaxed=="<raw>")
	id.start <- which(Relaxed=="<start>") 
	class[ id.raw ] <- "raw"
	if(length(id.start)>0)
	 class[ id.start ] <- "start"
	ids <- (1:length(n))[-c(id.raw, id.start)]
	name.vars <- names(obj$coars[[1]])
	n.vars <- length(name.vars)
	tab <- data.frame(n=n, ML1=ML1, class=class)
	
	dotcol <- rgb(0.2,0.2,0.8)
	selcol <- rgb(0.8,0.8,0.2)
	main.txt <- sprintf("Total units=%d, ML1=%.3f", n[id.raw], ML1[id.raw])
	if(length(id.start)>0)
	 main.txt <- sprintf("Initial matched units=%d, ML1=%.3f", n[id.start], 
						  ML1[id.start])
	
	plot(1/sqrt(n[ids]), ML1[ids],
         xlab="number matched, scaled as 1/sqrt(matched)", 
		 ylab="median of L1 profile", 
         pch=20, col=dotcol, ylim=range(ML1), xlim=range(1/sqrt(n)),
         main=main.txt,
		 axes=FALSE)
	axis(2)
	x1 <- pretty( 1/sqrt(n), 5)
	axis(1, x1, round(1/x1^2))
	box()
	points(1/sqrt(n[id.raw]), ML1[id.raw],   col="red", pch=20)
	text(1/sqrt(n[id.raw]), ML1[id.raw], "raw",  col="red",adj=-0.5, cex=0.7)
	if(length(id.start)>0){
		points(1/sqrt(n[id.start]), ML1[id.start],   col="green", pch=20)
		text(1/sqrt(n[id.start]), ML1[id.start], "start",  col="green",
			  adj=-0.5, cex=0.7)
	}

}



print.selected.cem <- function(x, ...){
	print(x$breaks)
}






imbspace.plot <- function(obj,group="1", data, explore=TRUE){
	if(!interactive() | !explore){
		imbspace.plot2(obj, group)
		return(NULL)	
	}
	
	
	if(class(obj) != "imbalance.space")
	stop("obj must be of class `imbalance.space'")
	
	haveTCL <- interactive()
	if(!capabilities("tcltk")){
		haveTCL <- FALSE	
		cat("\ntcltk support is absent")
	}
	
	if(haveTCL){
		have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
		if(have_ttk) {
			tkbutton <- ttkbutton
			tkframe <- ttkframe
			tklabel <- ttklabel
			tkradiobutton <- ttkradiobutton
		}
	}
	
	g <- sprintf("G%s",group)
	n <-  obj$space[[g]]
	ML1 <- obj$space$ML1
	Relaxed <- obj$space$Relaxed
	class <- rep("relax", length(n))
	id.raw <- which(Relaxed=="<raw>")
	id.start <- which(Relaxed=="<start>") 
	class[ id.raw ] <- "raw"
	if(length(id.start)>0)
		class[ id.start ] <- "start"
	ids <- (1:length(n))[-c(id.raw, id.start)]
	name.vars <- names(obj$coars[[1]])
	n.vars <- length(name.vars)
	tab <- data.frame(n=n, ML1=ML1, class=class)
	main.txt <- sprintf("Total units=%d, ML1=%.3f", n[id.raw], ML1[id.raw])
	
	if(length(id.start)>0)
		main.txt <- sprintf("Initial matched units=%d, ML1=%.3f", n[id.start], ML1[id.start])
	
	dotcol <- rgb(0.2,0.2,0.8)
	selcol <- rgb(0.8,0.8,0.2)
	plot(1/sqrt(n[ids]), ML1[ids],
         xlab="number matched, scaled as 1/sqrt(matched)", 
		 ylab="median of L1 profile", 
         pch=20, col=dotcol, ylim=range(ML1), xlim=range(1/sqrt(n)),
         main=main.txt, axes=FALSE)
	axis(2)
	x1 <- pretty( 1/sqrt(n), 5)
	axis(1, x1, round(1/x1^2))
	box()
	points(1/sqrt(n[id.raw]), ML1[id.raw],   col="red", pch=20)
	text(1/sqrt(n[id.raw]), ML1[id.raw], "raw",  col="red",adj=-0.5, cex=0.7)
	if(length(id.start)>0){
		points(1/sqrt(n[id.start]), ML1[id.start],   col="green", pch=20)
		text(1/sqrt(n[id.start]), ML1[id.start], "start",  col="green",adj=-0.5, cex=0.7)
	}
	idx.sav <- id.raw
	if(length(id.start)>0)
		idx.sav <- id.start

	old.idx <- NULL
	xy <- xy.coords(1/sqrt(n), ML1)
	
 	tmp.br <- obj$coars[[2]]
	if(length(id.start)>0)
		tmp.br <- obj$coars[[id.start]]
	new.br <- tmp.br
	
	goOn <- TRUE
	
	if(haveTCL){	
		tclServiceMode(FALSE)
		tt <- tktoplevel()
		
		tkwm.title(tt,"Modify CEM solution")
		entries <- list()
		tcvars <- list()
		infoText <- tclVar( sprintf("Total units=%d, ML1=%.3f", n[id.raw], ML1[id.raw]) )
		if(length(id.start)>0)
			infoText <- tclVar( sprintf("Matched units=%d, ML1=%.3f", n[id.start], ML1[id.start]) )
		label1 <- tklabel(tt, textvariable=infoText) 
		tkpack(label1)
		
		n.tmp.br <- length(tmp.br)
		for(i in 1:n.tmp.br){
			
			tcvars[[i]] <- tclVar( deparse( round(tmp.br[[i]], 2), width.cutoff=500) )  
			entries[[i]] <- tkentry(tt, width="100", textvariable=tcvars[[i]])	    
			
			tkpack( tklabel(tt, text=sprintf("Variable: %s", names(tmp.br)[i]) ))
			tkpack(  entries[[i]] )
		}

		tcvars[[n.tmp.br+1]] <- tclVar(  )  
		entries[[n.tmp.br+1]] <- tkentry(tt, width="100", textvariable=tcvars[[n.tmp.br+1]])	    
		
		tkpack( tklabel(tt, text="Additional CEM args:"))
		tkpack(  entries[[n.tmp.br+1]] )
		
		
		OnOK <- function(){
			other.args <- NULL
			cat("\n... running new cem...\n")
			n.tmp.br <- length(tmp.br)
			for(i in 1:n.tmp.br){
				vv <- names(tmp.br)[i]
				tmpc <- tclvalue( tcvars[[i]] )
				new.br[[i]] <<- try(eval(parse(text=tmpc)), silent = TRUE)
				if( class(new.br[[i]]) == "try-error"){
					cat(sprintf("\nError in settings cutpoints of variable << %s >>:\n\n >> %s <<\n\n Using original ones.\n", vv, tmpc) )
					new.br[[i]] <<- tmp.br[[i]] 
				}
			}
			tmpc <- tclvalue( tcvars[[n.tmp.br+1]] )
			other.args <- try(eval(parse(text=tmpc)), silent = TRUE)
			if( class(other.args) == "try-error"){
				cat(sprintf("\nError in additional CEM arguments specification. Ignoring them.\n", tmpc) )
				other.args <- NULL 
			} else 
			 other.args <- tmpc
			tclServiceMode(FALSE)	 
			if(!is.null(other.args))
			 eval(parse(text=sprintf("tmp.mat <- cem(obj$match$treatment, data=data[,c(obj$match$vars,obj$match$treatment) ], cutpoints=new.br, eval.imbalance=FALSE, %s)", other.args)))
		    else
			   tmp.mat <- cem(obj$match$treatment, data=data[,c(obj$match$vars,obj$match$treatment) ], cutpoints=new.br, eval.imbalance=FALSE)
			   
			tclServiceMode(TRUE)
			tmp.ML1 <- L1.meas(obj$match$groups, data=data[,obj$match$vars], breaks=obj$medianCP, weights=tmp.mat$w)$L1
			tmp.n <- tmp.mat$tab["Matched", g] 
			len.c <- length(n)
			n[len.c+1] <<- tmp.n 	 
			ML1[len.c+1] <<- tmp.ML1 	 
			Relaxed[len.c+1] <<- "<new>"	 
			obj$coars[[len.c+1]] <<- new.br	
			xy <<- xy.coords(1/sqrt(n), ML1)
			update.imbplot(len.c+1)	 
		}
		
		
		OK.but <-tkbutton(tt,text="   Run CEM with these coarsening   ",command=OnOK)
		
		
		tkpack(OK.but) 

		tclServiceMode(TRUE)
	}
	
	firstime <- TRUE
	update.imbplot <- function(mT){
		if(length(mT)==0)
		return(FALSE)
		
		
		idx <- which(n>=n[mT])
		idx2 <- which(ML1[idx]<= ML1[mT])
		idx <- idx[idx2]
		
		id.new <- which(Relaxed =="<new>")
		
		points(1/sqrt(n[idx.sav]), ML1[idx.sav],   col=dotcol, pch=20)
		points(1/sqrt(n[idx]), ML1[idx],   col=selcol, pch=20)
		points(1/sqrt(n[id.raw]), ML1[id.raw],   col="red", pch=20)
		points(1/sqrt(n[mT]), ML1[mT],   col="orange", pch=20)
		if(length(id.start)>0){
			points(1/sqrt(n[id.start]), ML1[id.start],   col="green", pch=20)
			text(1/sqrt(n[id.start]), ML1[id.start], "start",  col="green",adj=-0.5, cex=0.7)
		}
		if(length(id.new))
		 points(1/sqrt(n[id.new]), ML1[id.new],   col="cyan", pch=20)
		id.bad <- which(idx == id.raw)
		if(length(id.bad>0))
		 idx <- idx[-id.bad]
		if(length(idx)>0){
 		  y <- lapply(obj$coars[idx], function(a) unlist(lapply(a, length)))
		  x <- matrix(unlist(y), length(y), length(y[[1]]), byrow=TRUE) 
		
		  colnames(x) <- names(y[[1]])
		  tmp <- as.data.frame(x)
				
		  tmp2 <- data.frame(tmp, ML1=ML1[idx])	
		  rownames(tmp2) <- idx
		
		 if(haveTCL & !is.null(obj$coars[[mT]])){
			ttt <- obj$coars[[mT]]
 			for(i in 1:length(tmp.br)){
				tclvalue( tcvars[[i]] )  <-  deparse( round(ttt[[i]], 2), width.cutoff=500)
  	 		}

			tclvalue(infoText) <- sprintf("Matched units=%d, ML1=%.3f", n[mT], ML1[mT]) 

 			tkfocus(tt)	
		 }
		}

		
		old.idx <<- list(breaks = obj$coars[[mT]], n=n[mT], ML1=ML1[mT], medianCP=obj$medianCP)
	
		
		idx.sav <<- idx
		
		return(TRUE)
		
	}
	
	old.idx <- NULL
	update.imbplot( 2 )
	
	while(goOn){
		xy <<- xy.coords(1/sqrt(n), ML1)
		goOn <- update.imbplot(identify(xy, n=1,plot=FALSE)) 
	}
	
    if(haveTCL)
	 tkdestroy(tt)
	
	
	class(old.idx) <- "selected.cem"
  	old.idx
	
}

