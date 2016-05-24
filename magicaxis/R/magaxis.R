magaxis <-
function(side=1:2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE, unlog='Auto', mgp=c(2,0.5,0), mtline=2, xlab=NULL, ylab=NULL, crunch=TRUE, logpretty=TRUE, prettybase=10, hersh=FALSE, family='sans', frame.plot=FALSE, usepar=FALSE, ...){
majornlist=majorn
minornlist=minorn
unloglist=unlog
labelslist=labels
crunchlist=crunch
logprettylist=logpretty
prettybaselist=prettybase
if(length(majorn)==1 & length(side)>1){majornlist=rep(majorn,length(side))}
if(length(minorn)==1 & length(side)>1){minornlist=rep(minorn,length(side))}
if(length(unlog)==1 & length(side)>1 & (unlog[1]==T | unlog[1]==F | unlog[1]=='Auto')){unloglist=rep(unlog,length(side))}
if(length(labels)==1 & length(side)>1){labelslist=rep(labels,length(side))}
if(length(crunch)==1 & length(side)>1){crunchlist=rep(crunch,length(side))}
if(length(logpretty)==1 & length(side)>1){logprettylist=rep(logpretty,length(side))}
if(length(prettybase)==1 & length(side)>1){prettybaselist=rep(prettybase,length(side))}
if(unlog[1]=='x'){unloglist=rep(FALSE,length(side));unloglist[side %in% c(1,3)]=TRUE}
if(unlog[1]=='y'){unloglist=rep(FALSE,length(side));unloglist[side %in% c(2,4)]=TRUE}
if(unlog[1]=='xy' | unlog[1]=='yx'){unloglist=rep(TRUE,length(side))}

if(length(majornlist) != length(side)){stop('Length of majorn vector mismatches number of axes!')}
if(length(minornlist) != length(side)){stop('Length of minorn vector mismatches number of axes!')}
if(length(unloglist) != length(side)){stop('Length of unlog vector mismatches number of axes!')}
if(length(labelslist) != length(side)){stop('Length of labels vector mismatches number of axes!')}
if(length(crunchlist) != length(side)){stop('Length of crunch vector mismatches number of axes!')}
if(length(logprettylist) != length(side)){stop('Length of logpretty vector mismatches number of axes!')}
if(length(prettybaselist) != length(side)){stop('Length of prettybase vector mismatches number of axes!')}

currentfamily=par('family')
if(hersh & family=='serif'){par(family='HersheySerif')}
if(hersh & family=='sans'){par(family='HersheySans')}
if(hersh==F & family=='serif'){par(family='serif')}
if(hersh==F & family=='sans'){par(family='sans')}

if(usepar){
  tcl=par()$tcl
  mgp=par()$mgp
}

for(i in 1:length(side)){
		currentside=side[i]
    majorn=majornlist[i]
    minorn=minornlist[i]
		unlog=unloglist[i]
		labels=labelslist[i]
		crunch=crunchlist[i]
		logpretty=logprettylist[i]
    prettybase=prettybaselist[i]
  		lims=par("usr")
  		if(currentside %in% c(1,3)){
  		lims=lims[1:2];if(par('xlog')){logged=T}else{logged=F}
  		}else{
  		lims=lims[3:4];if(par('ylog')){logged=T}else{logged=F}
  		}
        lims=sort(lims)

        if(unlog=='Auto'){if(logged){unlog=T}else{unlog=F}}
        if(logged | unlog){usemultloc=(10^lims[2])/(10^lims[1])<50}else{usemultloc=F}
        	  		
  		if(unlog){
        sci.tick=maglab(10^lims,n=majorn,log=T,exptext=T,crunch=crunch,logpretty=logpretty,usemultloc=usemultloc,prettybase=prettybase, hersh=hersh)
        major.ticks = log10(sci.tick$tickat)
  		uselabels = sci.tick$exp
  		labloc = log10(sci.tick$labat)
        if(usemultloc==F){
          minors = log10(pretty(10^major.ticks[1:2],minorn+2))-major.ticks[1]
        }
 		}
 		if(logged & unlog==F){
 		sci.tick=maglab(10^lims,n=majorn,log=T,exptext=F,crunch=crunch,logpretty=logpretty,usemultloc=usemultloc,prettybase=prettybase, hersh=hersh)
 		major.ticks = log10(sci.tick$tickat)
  		uselabels = sci.tick$exp
  		labloc = log10(sci.tick$labat)
 		if(usemultloc==F){minors = log10(pretty(10^major.ticks[1:2],minorn+2))-major.ticks[1]}
 		}
 		if(logged==F & unlog==F){
 		sci.tick=maglab(lims,n=majorn,log=F,exptext=F,prettybase=prettybase, hersh=hersh)
 		major.ticks = sci.tick$tickat
  		uselabels = sci.tick$exp
  		labloc = sci.tick$labat
  		minors = pretty(major.ticks[1:2],minorn+2)-major.ticks[1]
 		}

 		if(logged){axis(side=currentside,at=10^major.ticks,tcl=tcl,labels=FALSE,mgp=mgp,...)}
 		else axis(side=currentside,at=major.ticks,tcl=tcl,labels=FALSE,mgp=mgp,...)
 		
  		if(labels){
  		    if(logged){axis(side=currentside,at=10^labloc,tick=F,labels=uselabels,mgp=mgp,...)}
  		    else axis(side=currentside,at=labloc,tick=F,labels=uselabels,mgp=mgp,...)
  		}
  		
  		if(usemultloc==F){
  		minors = minors[-c(1,length(minors))]
  		minor.ticks = c(outer(minors, major.ticks, `+`))
		if(logged){axis(currentside,at=10^minor.ticks,tcl=tcl*ratio,labels=FALSE,...)}
		else axis(currentside,at=minor.ticks,tcl=tcl*ratio,labels=FALSE,...)
		}
    if(is.null(xlab)==F & currentside==1){mtext(xlab,1,line=mtline)}
    if(is.null(ylab)==F & currentside==2){mtext(ylab,2,line=mtline)}  
}
if(frame.plot){box()}
par(family=currentfamily)
}
