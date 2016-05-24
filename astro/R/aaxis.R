aaxis = function(side = 1, at = NULL, labels = TRUE, majticks = TRUE, minticks = TRUE, nmaj = NULL, nmin = NULL, unlog = FALSE, format = NULL, digits = NULL, las = NULL, mgp = NULL, tcl = NULL, dexcl = 0.2, ...){
    #side = 1; at = NULL; labels = TRUE; majticks = TRUE; minticks = TRUE; nmaj = NULL; nmin = NULL; unlog = FALSE; format = NULL; digits = NULL; las = NULL; mgp = NULL; tcl = NULL; dexcl = 0.2
    
    # check inputs
    if(!is.null(nmaj)){if(nmaj<2){stop("Not enough major tick marks.")}}
    if(length(side)==0){stop("No axes will be created.")}
    if(length(side)>1){warning("Only using 1st 'side' parameter."); side = side[1]}
    
    # beginning tomfoolery
    if(is.null(las)){las=par("las")}
    if(is.null(mgp)){mgp=(par("mgp")-c(0.5,0.5,0))}
    if(is.null(tcl)){tcl=-1*par("tcl")}
    tcl.maj = tcl
    tcl.min = tcl.maj/2
    if(majticks=="FALSE"){tcl.maj=0}
    
    # logged? / majmin/majmax
    if(side==1|side==3){log=par("xlog"); majmin=par("usr")[1]; majmax=par("usr")[2]}
    if(side==2|side==4){log=par("ylog"); majmin=par("usr")[3]; majmax=par("usr")[4]}
    
    # major ticks
    ticks=axTicks(side)
    if(!is.null(nmaj)){ticks = seq(ticks[1],ticks[length(ticks)],len=nmaj)} 
    
    if(log){
#        ticks = ticks[which(log10(ticks)==as.integer(log10(ticks)))]
#        ticks = sort(unique(10^(c((log10(ticks)-4),(log10(ticks)-3),(log10(ticks)-2),(log10(ticks)-1),log10(ticks),(log10(ticks)+1),(log10(ticks)+2),(log10(ticks)+3)))))
        #if(floor(majmin)<=0){floormajmin = 0.01}else{floormajmin = floor(majmin)}
        ticks = 10^seq(floor(majmin),ceiling(majmax),by=1)
    }else if(unlog){
#        ticks = ticks[which(ticks==as.integer(ticks))]
#        ticks = sort(unique((c((ticks-4),(ticks-3),(ticks-2),(ticks-1),ticks,(ticks+1),(ticks+2),(ticks+3)))))
        ticks = seq(floor(majmin),ceiling(majmax),by=1)
    }else{
        ticks = seq(ticks[1]-diff(ticks)[1],ticks[length(ticks)]+diff(ticks)[1],by=diff(ticks)[1])
    }
    
    # user defined at?
    udefat = FALSE
    if(!is.null(at[1])){
        udefat = TRUE
        if(unlog){
            ticks = log10(at)   #log10(c(ticks[1],at,ticks[length(ticks)]))
        }else{
            ticks = at          #c(ticks[1],at,ticks[length(ticks)])
        }
    }
    
    # rounding & backup
    #ticks = signif(ticks,digits=9)
    #ticks = round(ticks,digits=9)
    if(abs((majmax-majmin))>1E-31){if(any(abs(ticks)<1E-31)){ticks[abs(ticks)<1E-31]=0}}
    axticks = ticks
    
    # define labels
    labs = FALSE
    if(is.logical(labels[1])){
        if(labels[1]==TRUE & !unlog){
            #labs = as.character(formatC(ticks,format="g"))
            labs = ticks
        }else if(labels[1]==TRUE & unlog){
            labs = 10^(ticks)
        }
    }else{
        if(!udefat){
            labels = rep(labels,ceiling(length(ticks)/(length(labels)+2)))[1:(length(ticks)-2)]
            #labs = formatC(c("",labels,""),format="g")
            labs = c("",labels,"")
        }else{
            labels = rep(labels,ceiling(length(ticks)/(length(labels))))[1:(length(ticks))]
            #labs = formatC(labels,format="g")
            labs = labels
        }
    }
    
    # format
    if(is.character(format)){
        if(format=="NULL"){
            format = NULL
        }
    }
    fp = FALSE
    powers = {}
    if(!is.null(format) & labels[1]!=FALSE){
        if(is.null(digits)){digits=1}
        if(format!="x" & format!="X" & format!="p"){
            labs = formatC(as.numeric(labs), format=format, digits=digits)
        }else if((format=="x" | format=="X" | format=="p")){
            labs = formatC(as.numeric(labs), format="E", digits=digits)
            if(format=="x"){format="\u00D7"}
            for(ii in 1:length(labs)){
                temp = strsplit(labs[[ii]],"E")[[1]]
                if(as.numeric(temp[1])!=0){
                    powers = c(powers, as.numeric(temp[2]))
                    if(format!="p"){
                        labs[ii] = paste(temp[1],format,"10",sep="")
                    }else{
                        fp = TRUE
                        if(as.numeric(temp[1])<0){
                            labs[ii] = "-10"
                        }else{
                            labs[ii] = "10"
                        }
                    }
                }else{
                    powers = c(powers, "")
                    labs[ii] = "0"
                }
            }
        }
    }
    
    # plot major axis
    if(length(powers)==0){
        axis(side, at=ticks, labels=labs, las=las, mgp=mgp, tcl=tcl.maj, ...)
    }else{
        axis(side, at=ticks, labels=FALSE, las=las, mgp=mgp, tcl=tcl.maj, ...)
        for(jj in 1:length(powers)){
            axis(side, at=ticks[jj], labels=bquote(paste(.(labs[jj])^{.(powers[jj])})), las=las, mgp=mgp, tick=FALSE, ...)
        }
    }
    
    if(minticks=="TRUE"){
        # minor ticks
        #minspace = c(1,2,4,5,10)
        minax = min(c(axticks[1],axticks[2]))
        majax = max(c(axticks[1],axticks[2]))
        minspace = {}
        for(i in 1:10){
            pmin = pretty(c(minax,majax),n=i)
            if(pmin[1]==minax & pmin[length(pmin)]==majax){
                minspace = c(minspace,(length(pmin)-1))
            }
        }
        if(length(minspace)==0){minspace = 1}
        minspace = unique(minspace)
        minsep = 0.1
        
#        # even major tick marks
#        if((unique(diff(axticks))%%2) == 0){
#            minspace = c(1,2,4,10)
#        }
        
        # user defined min tick spacing
        if(!is.null(nmin) & !unlog){minspace = (nmin+1)}
        
        # minor axes side info
        if(side==1 | side==3){
            pinside=1
            parusrlo=par("usr")[1]
            parusrhi=par("usr")[2]
        }else{
            pinside=2
            parusrlo=par("usr")[3]
            parusrhi=par("usr")[4]
        }
        
        # number of desired minor ticks
        mylo = min(parusrlo,parusrhi)
        myhi = max(parusrlo,parusrhi)
        preticks = abs(((par("pin")[pinside]/length(which(axticks>mylo & axticks<myhi)))/minspace)-minsep)
        desticks = minspace[which.min(preticks)]
        
        # minor tick vector
        if(log){
            logvect = seq(0,10,len=desticks+1)[2:desticks]
            logpos = expand.grid(log10(1:9),log10(axticks))
            mticks = 10^(logpos[,1]+logpos[,2])
        }else if(unlog){
            logvect = seq(0,10,len=desticks+1)[2:desticks]
            logpos = expand.grid(log10(1:9),axticks[1:(length(axticks)-1)])
            mticks = logpos[,1]+logpos[,2]
        }else{
            ticksep = diff(seq(axticks[1],axticks[2],len=desticks+1))[1]
            mticks = seq(axticks[1],axticks[length(axticks)],by=ticksep)
        }
        
        # rounding
        if(any(abs(mticks)<9E-31) & abs(mticks[length(mticks)]-mticks[1])>1E-12){
            mticks[which(abs(mticks)<9E-31)]=0
        }
        
        # rounding
        if(!unlog & !log){mticks = signif(mticks,digits=9)}  # rounding

        # remove major ticks from minor tick vector        
        if(majticks!="FALSE"){
            mticks = mticks[-which(as.character(round(mticks,digits=15))%in%as.character(ticks))]
        }
        
        # not enough major tick marks on axis?
        mlabs = FALSE
        if(length(labs)<5 & !fp){
            if(is.logical(labels[1])){
                if(labels[1]==TRUE){
                    if(!unlog){
                        mlabs = as.character(mticks)
                    }else if(unlog){
                        mlabs = as.character(10^(mticks))
                    }
                
                    # remove minor ticks from minor tick label vector which are too close to a major tick
                    dexcl = dexcl  # distance from major ticks within which minor ticks should be removed (inches)
                    ratio = par("pin")[pinside]/(myhi-mylo)
                    dticks = (ticks-min(ticks))*ratio
                    dmticks = (mticks-min(ticks))*ratio
                    
                    for(i in 1:length(dticks)){
                        # centres
                        if(length(which(abs(dmticks-dticks[i])<dexcl))>0){
                            mlabs[which(abs(dmticks-dticks[i])<dexcl)]=""
                        }
                        # centre - text width/2
                        if(length(which(abs(abs(dmticks-dticks[i])-(strwidth(mlabs)/2))<dexcl))>0){
                            mlabs[which(abs(abs(dmticks-dticks[i])-(strwidth(mlabs)/2))<dexcl)]=""
                        }
                        # centre + text width/2
                        if(length(which(abs(abs(dmticks-dticks[i])+(strwidth(mlabs)/2))<dexcl))>0){
                            mlabs[which(abs(abs(dmticks-dticks[i])+(strwidth(mlabs)/2))<dexcl)]=""
                        }
                    }
                }
            }
        }
        
        # format
        powers = {}
        if(!is.null(format) & mlabs[1]!="FALSE" & !(length(unique(mlabs))==1 & unique(mlabs)[1]=="")){
            if(is.null(digits)){digits=1}
            if(format!="x" & format!="X" & format!="\u00D7" & format!="p"){
                mlabs[mlabs!=""] = formatC(as.numeric(mlabs[mlabs!=""]), format=format, digits=digits)
            }else if((format=="x" | format=="X" | format=="\u00D7" | format=="p")){
                omlabs = mlabs
                mlabs = formatC(as.numeric(mlabs[mlabs!=""]), format="E", digits=digits)
                if(format=="x"){format="\u00D7"}
                for(ii in 1:length(mlabs)){
                    temp = strsplit(mlabs[ii],"E")[[1]]
                    if(as.numeric(temp[1])!=0){
                        powers = c(powers, as.numeric(temp[2]))
                        if(format!="p"){
                            mlabs[ii] = paste(temp[1],format,"10",sep="")
                        }else{
                            if(as.numeric(temp[1])<0){
                                mlabs[ii] = "-10"
                            }else{
                                mlabs[ii] = "10"
                            }
                        }
                    }else{
                        powers = c(powers, "")
                        mlabs[ii] = "0"
                    }
                }
                omlabs[omlabs!=""] = mlabs
                mlabs = omlabs
            }
        }
        
        # plot min ticks
        if(length(powers)==0){
            axis(side, at=mticks, labels=mlabs, las=las, mgp=mgp, tcl=tcl.min, ...)
        }else{
            axis(side, at=mticks, labels=FALSE, las=las, mgp=mgp, tcl=tcl.maj, ...)
            for(jj in 1:length(powers)){
                axis(side, at=mticks[jj], labels=bquote(paste(.(mlabs[jj])^{.(powers[jj])})), las=las, mgp=mgp, tick=FALSE, ...)
            }
        }
        
    }
    
}
