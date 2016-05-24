my.strata <- function (..., nameslist,   sep = ", ") 
{									#nameslist = lista imen spremenljivk
    words <- as.character((match.call())[-1])				#ime podatkov
    allf <- list(...)							#podatki
    if (length(allf) == 1 && is.list(ttt <- unclass(allf[[1]]))) 	#so samo eni podatki
        allf <- ttt							#ohranim le podatke (ne listo podatkov), v obliki list
    nterms <- length(allf)						#nterms= st. spremenljivk +1 (row.names)
    if (is.null(names(allf)))						#ce ni imen 
        argname <- words[1:nterms]					#jih dam
    else argname <- ifelse(names(allf) == "", words[1:nterms], 		#ce so prazna jih dam
        names(allf))							#imena so v argname
   
  
   varnames <- names(nameslist)
  	
    #1. iteracija
    what <- allf[[1]]							#prva spremenljivka
   
    for(it in 1:length(varnames)){
       if (length(grep(varnames[it],names(allf)[[1]]))) break		#poiscem ji mesto v svojem poimenovanju
       }
   
   
    if (is.null(levels(what)))     what <- factor(what)			#ce se ni, jo prisilimo v faktorsko
    levs <- unclass(what) - 1						#nastavim prvi level = 0
    wlab <- levels(what)						#imena faktorjev
    labs <- paste(argname[1], wlab, sep = "=")				#prvo ime = 0/1
    
    labsnow <- 1
    allab <- NULL
    
    dd <- length(nameslist[[it]])
    if(dd!=2) {
    	mylabs <- rep(argname[1],length(wlab))
    	mylabs[wlab==0] <- ""
    }
    else mylabs <- labs
   
 
     for (i in (1:nterms)[-1]) {
        if(length(grep(varnames[labsnow],names(allf)[[i]]))==0){					#ce je zdaj to nova spremenljivka, moram najprej ustimat prejsnjo
    		mylabs[mylabs==""] <- nameslist[[labsnow]][1] 
    		if(!any(allab!=""))allab <- paste(allab,mylabs,sep="")					#the first time - do not separate by comma
    		else allab <- paste(allab,mylabs,sep=",")
    		mylabs <- rep("",length(mylabs))
    		labsnow <- labsnow+1
        }
        what <- allf[[i]]
        if (is.null(levels(what))) what <- factor(what)
        wlev <- unclass(what) - 1
        wlab <- levels(what)
        labsnew <- format(paste(argname[i], wlab, sep = "="))
	levs <- wlev + levs * (length(wlab))
        a <- rep(labs, rep(length(wlab), length(labs)))
	b <- rep(wlab, length(labs))
	
	mya <- rep(mylabs, rep(length(wlab), length(labs)))
	allab <- rep(allab,rep(length(wlab), length(labs)))
	myb <- rep(argname[i],length(labs)*length(wlab))
	for(it in 1:length(varnames)){								#it se ustavi pri trenutni spremenljivki
	    if (length(grep(varnames[it],names(allf)[[i]]))) break
    	}
	dd <- length(nameslist[[it]])
	if(dd==2)myb <- paste(myb,rep(wlab,length(labs)),sep="=")
	else myb[rep(wlab,length(labs))==0] <- ""
	mylabs <- paste(mya,myb,sep="")
	
        labs <- paste(a,b, sep = sep)
    }
    mylabs[mylabs==""] <- nameslist[[labsnow]][1] 
    if(!any(allab!=""))allab <- paste(allab,mylabs,sep="")
    else allab <- paste(allab,mylabs,sep=",")
    levs <- levs + 1
    ulevs <- sort(unique(levs[!is.na(levs)]))
    levs <- match(levs, ulevs)
    labs <- labs[ulevs]
    allab <- allab[ulevs]
    factor(levs, labels = allab)
}


