aplot = function(x, y = NULL, z = NULL, xlim = NULL, ylim = NULL, zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL, col = NULL, axes = TRUE, side = 1:4, labels = TRUE, majticks = TRUE, minticks = TRUE, nxmaj = NULL, nymaj = NULL, nxmin = NULL, nymin = NULL, xat = NULL, yat = NULL, log = "", unlog = FALSE, xformat = NULL, yformat = NULL, digits = 0, cex = 1, xlabpos = 1, ylabpos = 2, zcol = NULL, cb = FALSE, cbpos = 4, cbsep = 1.5, cbspan = 2, cbinset = 1, cbx1 = NULL, cbx2 = NULL, cby1 = NULL, cby2 = NULL, cblegend = NULL, cbsteps = 250, las = 0, mgp = c(2.5,0.5,0), tcl = 0.5, dexcl = 0.2, cex.axis = 1, cex.cb = 1, zline = mgp[1]+1, col.axes = "black", col.axis = "black", add = FALSE, type = NULL, bgcol = NULL, ...){
    #y = NULL; z = NULL; xlim = NULL; ylim = NULL; zlim = NULL; xlab = NULL; ylab = NULL; zlab = NULL; col = NULL; axes = TRUE; side = 1:4; labels = TRUE; majticks = TRUE; minticks = TRUE; nxmaj = NULL; nymaj = NULL; nxmin = NULL; nymin = NULL; xat = NULL; yat = NULL; log = ""; unlog = FALSE; xformat = NULL; yformat = NULL; digits = 0; cex = 1; xlabpos = 1; ylabpos = 2; zcol = NULL; cb = FALSE; cbpos = 4; cbsep = 1.5; cbspan = 2; cbinset = 1; cbx1 = NULL; cbx2 = NULL; cby1 = NULL; cby2 = NULL; cblegend = NULL; cbsteps = 250; las = 0; mgp = c(2.5,0.5,0); tcl = 0.5; dexcl = 0.2; cex.axis = 1; cex.cb = 1; zline = mgp[1]+1; col.axes = "black"; col.axis = "black"; add = FALSE; type = NULL; bgcol = NULL
    
    # unlog
    if(unlog!="FALSE" & unlog!="F"){
        if(unlog=="TRUE" | unlog=="T"){unlog=1:4
        }else{
            unlog = strsplit(unlog, "")[[1]]
            temp = {}
            if("x"%in%unlog){temp=c(temp, 1, 3)}
            if("y"%in%unlog){temp=c(temp, 2, 4)}
            if("z"%in%unlog){temp=c(temp, 5)}
            unlog = temp
        }
        if(1%in%unlog & 3%in%unlog & !is.null(xlim)){xlim=log10(xlim)}
        if(2%in%unlog & 4%in%unlog & !is.null(ylim)){ylim=log10(ylim)}
        if(5%in%unlog & !is.null(zlim)){zlim=log10(zlim)}
    }
    
    # z axis
    if(is.null(z[1])){
        z = 0
        zcol = "black"      # no single colour, no z data
        zlim = c(0,1)
    }
    if(any(is.na(z))){
        x = x[!is.na(z)]
        y = y[!is.na(z)]
        z = z[!is.na(z)]
    }
    if(is.null(zlim[1])){zlim=range(z, na.rm=TRUE)}
    z[z<zlim[1]] = zlim[1]
    z[z>zlim[2]] = zlim[2]
    base = (z-zlim[1])/(zlim[2]-zlim[1])
    if(is.null(zcol)){
        zcol = hsv(seq(2/3,0,len=200))
    }
    if(is.null(col[1])){col = zcol[(base*(length(zcol)-1))+1]}
    
    # colour bar
    parmarorig = temp = par("mar")
    if(cb){
        #extramar=0
        #if(!is.null(zlab)){extramar=extramar+1}
        #if(cbpos==1){extramar = extramar+0
        #}else if(cbpos==2){extramar = extramar+0
        #}else if(cbpos==3){extramar = extramar+0
        #}else if(cbpos==4){extramar = extramar+2}
        #temp[cbpos] = temp[cbpos]+extramar
        par("mar"=temp)
        if(is.null(cblegend[1])){
            if(length(z)>1){
                if(!5%in%unlog){
                    n = 3
                    cblegend = round(pretty(zlim, n=n), digits=9)
                    while((cblegend[1]!=zlim[1] | cblegend[length(cblegend)]!=zlim[2]) & n < 9){
                        n = n+1
                        cblegend = round(pretty(zlim, n=n), digits=9)
                    }
                    if(cblegend[1]!=zlim[1] | cblegend[length(cblegend)]!=zlim[2]){
                        cblegend = formatC(zlim,digits=2,format="f")
                    }
                }else{
                    bigdist = 10^seq(floor(zlim[1]),ceiling(zlim[2]),by=1)
                    alldist = 10^seq(zlim[1],zlim[2],len=length(zcol))
                    smalldist = unique(signif(alldist, digits=1))
                    if(length(bigdist)<4){samp = smalldist}else{samp=bigdist}
                    cblegend = rep("", length(zcol))
                    for(i in 1:length(samp)){
                        poss = which(signif(alldist, digits=1)==samp[i])
                        poss[which.min(abs(alldist[poss]-samp[i]))]
                        cblegend[poss[which.min(abs(alldist[poss]-samp[i]))]] = samp[i]
                    }
                }
            }else{
                cblegend = "doh!"
            }
        }else{
            if(5%in%unlog){
                bigdist = cblegend
                alldist = 10^seq(zlim[1],zlim[2],len=length(zcol))
                samp=bigdist
                cblegend = rep("", length(zcol))
                for(i in 1:length(samp)){
                    poss = which(signif(alldist, digits=1)==samp[i])
                    poss[which.min(abs(alldist[poss]-samp[i]))]
                    cblegend[poss[which.min(abs(alldist[poss]-samp[i]))]] = samp[i]
                }
            }
        }
    }
    
    # type
    if(is.null(type)){
        if(class(x)=="density" | class(x)=="function"){
            type="l"
        }else{
            type="p"
        }
    }
    
    # background colour
    if(!is.null(bgcol)){
        type2 = type
        type = "n"
    }
    
    # lab positions
    if(xlabpos!=1){
        xlab2 = xlab
        xlab = ""
    }
    if(ylabpos!=2){
        ylab2 = ylab
        ylab = ""
    }
    
    # plot
    if(!add){
        if(class(x)=="density"){
            temp = plot(x, y=NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE, las=las, mgp=mgp, tcl=tcl, col=col, cex=cex, type=type, log=log, ...)
        }else{
            temp = plot(x, y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE, las=las, mgp=mgp, tcl=tcl, col=col, cex=cex, type=type, log=log, ...)
        }
    }else{
        if(class(x)=="density"){
            plot(x, y=NULL, col=col, cex=cex, add=TRUE, type=type, log=log, xlab=xlab, ylab=ylab, ...)
        }else if(class(x)=="function"){
            plot(x, y, col=col, cex=cex, add=TRUE, type=type, log=log, ...)
        }else{
            points(x, y, col=col, cex=cex, type=type, ...)
        }
    }
    
    # background colour
    if(!is.null(bgcol)){
        xrect = c(par("usr")[1],par("usr")[2])
        yrect = c(par("usr")[3],par("usr")[4])
        if(par("xlog")){xrect=10^xrect}
        if(par("ylog")){yrect=10^yrect}
        rect(xrect[1], yrect[1], xrect[2], yrect[2], col=bgcol, border=bgcol)
        if(class(x)=="density"){
            plot(x, y=NULL, col=col, cex=cex, add=TRUE, type=type2, log=log, xlab=xlab, ylab=ylab, ...)
        }else{
            if(!is.null(temp)){x = temp$x; y = temp$y}
            points(x, y, col=col, cex=cex, type=type2, ...)
        }
        add = FALSE
    }
    
    # lab positions
    if(xlabpos!=1){
        mtext(xlab2, side=xlabpos, line=mgp[1], outer=FALSE)
    }
    if(ylabpos!=2){
        mtext(ylab2, side=ylabpos, line=mgp[1], outer=FALSE)
    }
    
    # extras
    if(labels[1]=="TRUE"){labels=1:2}
    if(majticks[1]=="TRUE"){majticks=1:4}
    if(minticks[1]=="TRUE"){minticks=1:4}
    
    # axes
    if(axes & !add){
        if(1%in%side)
            aaxis(side=1, at=xat, labels=1%in%labels, majticks=1%in%majticks, minticks=1%in%minticks, nmaj=nxmaj, nmin=nxmin, unlog=1%in%unlog, format=xformat, digits=digits, las=las, mgp=mgp, tcl=tcl, dexcl=dexcl, cex.axis=cex.axis, col=col.axes, col.axis=col.axis)
        if(2%in%side)
            aaxis(side=2, at=yat, labels=2%in%labels, majticks=2%in%majticks, minticks=2%in%minticks, nmaj=nymaj, nmin=nymin, unlog=2%in%unlog, format=yformat, digits=digits, las=las, mgp=mgp, tcl=tcl, dexcl=dexcl, cex.axis=cex.axis, col=col.axes, col.axis=col.axis)
        if(3%in%side)
            aaxis(side=3, at=xat, labels=3%in%labels, majticks=3%in%majticks, minticks=3%in%minticks, nmaj=nxmaj, nmin=nxmin, unlog=3%in%unlog, format=xformat, digits=digits, las=las, mgp=mgp, tcl=tcl, dexcl=dexcl, cex.axis=cex.axis, col=col.axes, col.axis=col.axis)
        if(4%in%side)
            aaxis(side=4, at=yat, labels=4%in%labels, majticks=4%in%majticks, minticks=4%in%minticks, nmaj=nymaj, nmin=nymin, unlog=4%in%unlog, format=yformat, digits=digits, las=las, mgp=mgp, tcl=tcl, dexcl=dexcl, cex.axis=cex.axis, col=col.axes, col.axis=col.axis)
    }
    
    # colour bar
    if(cb){
        usr = par("usr")
        fin = par("fin")
        fus = c(diff(c(usr[1],usr[2])),diff(c(usr[3],usr[4])))  # usr dimensions
        innerline=strheight("L", units="inches")    # inches per line
        conv = fus/fin
        #outerline = as.numeric(unique(as.character(par("mai")/par("mar"))))
        #cbdist = (cbsep+cbspan)*(outerline/innerline)
        if(cbpos==1){            
            xlower = usr[1]+((conv[1])*(cbinset*innerline))
            xupper = usr[2]-((conv[1])*(cbinset*innerline))
            yupper = usr[3]-((conv[2])*(cbsep*innerline))
            ylower = yupper-((conv[2])*(cbspan*innerline))
            colstep = diff(c(xlower,xupper))/length(cblegend)
            strwidths = strwidth(cblegend)
            align = "rb"
            gradient = "x"
        }else if(cbpos==2){
            xupper = usr[1]-((conv[1])*(cbsep*innerline))
            xlower = xupper-((conv[1])*(cbspan*innerline))
            ylower = usr[3]+((conv[2])*(cbinset*innerline))
            yupper = usr[4]-((conv[2])*(cbinset*innerline))
            colstep = diff(c(ylower,yupper))/length(cblegend)
            strwidths = strheight(cblegend)
            align = "lt"
            gradient = "y"
        }else if(cbpos==3){
            xlower = usr[1]+((conv[1])*(cbinset*innerline))
            xupper = usr[2]-((conv[1])*(cbinset*innerline))
            ylower = usr[4]+((conv[2])*(cbsep*innerline))
            yupper = ylower+((conv[2])*(cbspan*innerline))
            colstep = diff(c(xlower,xupper))/length(cblegend)
            strwidths = strwidth(cblegend)
            align = "lt"
            gradient = "x"
        }else if(cbpos==4){
            xlower = usr[2]+((conv[1])*(cbsep*innerline))
            xupper = xlower+((conv[1])*(cbspan*innerline))
            ylower = usr[3]+((conv[2])*(cbinset*innerline))
            yupper = usr[4]-((conv[2])*(cbinset*innerline))
            colstep = diff(c(ylower,yupper))/length(cblegend)
            strwidths = strheight(cblegend)
            align = "rb"
            gradient = "y"
        }
        if(par("xlog")){xlower=10^xlower; xupper=10^xupper}
        if(par("ylog")){ylower=10^ylower; yupper=10^yupper}
        if(is.null(cbx1)){cbx1 = xlower}
        if(is.null(cbx2)){cbx2 = xupper}
        if(is.null(cby1)){cby1 = ylower}
        if(is.null(cby2)){cby2 = yupper}
        # cblegend
        valpos = try(which(cblegend%in%(unique(cblegend)[unique(cblegend)!=""])),silent=TRUE)
        if(length(grep("Error",valpos))>0){
            valwidths = (1:length(cblegend))*colstep
            srad = 0.5
            upper = valwidths+(srad*strwidths)
            lower = valwidths-(srad*strwidths)
            if(length(valpos)>3){
                for(i in c(1,((length(valpos)):2))){
                    poss = which(upper>lower[valpos[i]] & lower<upper[valpos[i]])
                    poss = poss[-which(poss==valpos[i])]
                    if(cblegend[valpos[i]]!=""){cblegend[poss]=""}
                }
            }
        }
        color.legend(xl=cbx1, yb=cby1, xr=cbx2, yt=cby2, legend=cblegend, rect.col=zcol, cex=cex.cb, align=align, gradient=gradient, pos=cbpos, offset=0.25)
        if(!is.null(zlab)){
            mtext(zlab, side=cbpos, line=zline, outer=FALSE)
        }
        par("mar"=parmarorig)
    }
    
}



