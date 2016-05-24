acb = function(zlim = NULL, zlab = NULL, unlog = FALSE, cex = 1, zcol = NULL, cbpos = 4, cbsep = 0, cbspan = 0, cbinset = 1, cbx1 = NULL, cbx2 = NULL, cby1 = NULL, cby2 = NULL, cblegend = NULL, cex.cb = 1, zline = 2.5, ...){
    #zlim=NULL; zlab=NULL; unlog=FALSE; cex=1; zcol=NULL; cbpos=4; cbsep=2; cbspan=3; cbinset=1; cbx1=NULL; cbx2=NULL; cby1=NULL; cby2=NULL; cblegend=NULL; cex.cb=1; zline=3
    
    # unlog
    if(unlog==TRUE & !is.null(zlim)){zlim=log10(zlim)}
    
    # z axis
    if(is.null(zlim[1])){zlim=c(0,1)}
    if(is.null(zcol)){zcol = hsv(seq(2/3,0,len=200))}
    
    # colour bar
    parmarorig = temp = par("mar")
    extramar=0
    if(!is.null(zlab)){extramar=extramar+1}
    if(cbpos==1){extramar = extramar+0
    }else if(cbpos==2){extramar = extramar+0
    }else if(cbpos==3){extramar = extramar+0
    }else if(cbpos==4){extramar = extramar+2}
    temp[cbpos] = temp[cbpos]+extramar
    par("mar"=temp)
    if(is.null(cblegend[1])){
        if(!unlog){
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
        if(unlog){
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
    
    # dummy plot
    plot(sin, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
    
    # colour bar
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
        ylower = 0.0-((conv[2])*(cbspan*innerline))
        yupper = usr[4]-((conv[2])*(cbsep*innerline))
        colstep = diff(c(xlower,xupper))/length(cblegend)
        strwidths = strwidth(cblegend)
        align = "rb"
        gradient = "x"
    }else if(cbpos==2){
        xlower = 0.0-((conv[1])*(cbspan*innerline))
        xupper = usr[2]-((conv[1])*(cbsep*innerline))
        ylower = usr[3]+((conv[2])*(cbinset*innerline))
        yupper = usr[4]-((conv[2])*(cbinset*innerline))
        colstep = diff(c(ylower,yupper))/length(cblegend)
        strwidths = strheight(cblegend)
        align = "lt"
        gradient = "y"
    }else if(cbpos==3){
        xlower = usr[1]+((conv[1])*(cbinset*innerline))
        xupper = usr[2]-((conv[1])*(cbinset*innerline))
        ylower = usr[3]+((conv[2])*(cbsep*innerline))
        yupper = 1.0+((conv[2])*(cbspan*innerline))
        colstep = diff(c(xlower,xupper))/length(cblegend)
        strwidths = strwidth(cblegend)
        align = "lt"
        gradient = "x"
    }else if(cbpos==4){
        xlower = usr[1]+((conv[1])*(cbsep*innerline))
        xupper = 1.0+((conv[1])*(cbspan*innerline))
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
    valpos = which(cblegend%in%(unique(cblegend)[unique(cblegend)!=""]))
    valwidths = (1:length(cblegend))*colstep
    srad = 1
    upper = valwidths+(srad*strwidths)
    lower = valwidths-(srad*strwidths)
    if(length(valpos)>3){
        for(i in c(1,((length(valpos)):2))){
            poss = which(upper>lower[valpos[i]] & lower<upper[valpos[i]])
            poss = poss[-which(poss==valpos[i])]
            if(cblegend[valpos[i]]!=""){cblegend[poss]=""}
        }
    }
    color.legend(xl=cbx1, yb=cby1, xr=cbx2, yt=cby2, legend=cblegend, rect.col=zcol, cex=cex.cb, align=align, gradient=gradient, pos=cbpos, offset=0.25, ...)
    if(!is.null(zlab)){
        mtext(zlab, side=cbpos, line=zline, outer=FALSE, cex=cex)
    }
    par("mar"=parmarorig)
    
}



