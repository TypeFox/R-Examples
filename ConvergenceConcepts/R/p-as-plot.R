p.as.plot <- function(data,critp,critas,epsilon=0.05,nb.sp=10,mode="p") {
  
# A few functions to create interactivity tools  
  plotleft <- function(...) {
  # To draw nb.sp sample paths on right plot window

    # The moving bar width
    if (mode=="p") {barwidth <- 0.9}
    if (mode=="as") {barwidth <- 0.1}
    

    datafull <- data
    data <- data[,NameVal1:NameVal2]
    
    mini.data <- min(data[1:nb.sp,])
    maxi.data <- max(data[1:nb.sp,])

    if (maxi.data<epsilon) {maxi.data <- epsilon+0.05}
    if (mini.data>-epsilon) {mini.data <- -epsilon-0.05}

    # We create the plot window with well chosen axis
    plot.new()
    plot.window(xlim=c(1,NameVal2-NameVal1+1),ylim=c(mini.data,maxi.data))
    x <- seq(from=1,length=10,to=(NameVal2-NameVal1+1))
    axis(1,at=x,labels=as.vector(round(x+NameVal1-1)))
    axis(2,at=seq(from=mini.data,to=maxi.data,length=10),labels=signif(seq(from=mini.data,to=maxi.data,length=10),2))
    box()
    
    # This is the grey moving bar
    if (mode=="p") {
    rect(nn-NameVal1+1-barwidth, mini.data, nn-NameVal1+1+barwidth, maxi.data, density=NULL,col="lightgrey")
    }
    if (mode=="as") {
    rect(nn-NameVal1+1-barwidth, mini.data, nmax-NameVal1+1+barwidth, maxi.data, density=NULL,col="lightgrey")
    }
    
    # This draws the nb.sp sample paths
    visualize.sp(data,epsilon,nb.sp,plotfunc=points,col="grey")

    # This puts red marks for each sample path going off the band [-epsilon;+epsilon] at (or after if a.s. convergence is considered) n=nn
    if (mode=="p") {
    masque <- abs(data[1:nb.sp,nn-NameVal1+1])>epsilon
    nb.outside <- sum(masque)
    points(rep(nn-NameVal1+1,nb.outside),data[1:nb.sp,nn-NameVal1+1][masque],col="red",pch=95)
     }
    if (mode=="as") {
    masque <- apply(abs(datafull[1:nb.sp,nn:nmax])>epsilon,FUN=function(x){(1:length(x))[x==1][1]},MARGIN=1)
    nb.outside <- sum(!is.na(masque))
    points(nn+masque[!is.na(masque)]-NameVal1,datafull[cbind(as.vector(1:nb.sp)[!is.na(masque)],nn+masque[!is.na(masque)]-1)],col="red",pch=95)
    }

    # The title
   if (mode=="p") {title(bquote(textstyle(bold(atop("Only" ~ .(nb.sp) ~ "sample paths are plotted among which" ~ .(nb.outside)  ~ "go off [" -epsilon ~ ";" +epsilon ~ "] ","in the bar at position n=" ~ .(nn) ~ "(in red).")))))
    title(sub=paste("Bar at position n=",format(nn),sep=""))
                 }
    if (mode=="as") {title(bquote(textstyle(bold(atop("Only" ~ .(nb.sp) ~ "sample paths are plotted among which" ~ .(nb.outside)  ~ "go off [" -epsilon ~ ";" +epsilon ~ "]","in the block beginning at position n=" ~ .(nn) ~ "(1st occurence in red).")))))
    title(sub=paste("Block beginning at position n=",format(nn),sep=""))
                  }
                   
    if (mode=="p") {abline(v=nn-NameVal1+1)}

    abline(h=0,col="black")

    # Some useful information to help reading the plot
    mtext(expression(X[n~','~omega]-X[omega]),side=3,line=-3.5,outer=TRUE,at=0.11,cex=0.8)
    mtext(bquote(bold(n[1] ~ "=")),side=1,line=-4,outer=TRUE,at=0.11)
    mtext(bquote(bold("=" ~ n[2])),side=1,line=-4,outer=TRUE,at=0.97)
    mtext(bquote(epsilon == .(epsilon) ~ scriptstyle(phantom("(step back coefficient)"))),side=1,line=-3,outer=TRUE,at=0.2,col="blue")
    if (mode=="as") {mtext(bquote(K == .(K) ~ scriptstyle("(step back coefficient)")),side=1,line=-2,outer=TRUE,at=0.2)}
   
  }
  
  
  plotright <- function(...) {
  # To draw the convergence criterion
    
    # We create the plot window with well chosen axis
    plot.new()
    plot.window(xlim=c(1,NameVal2-NameVal1+1),ylim=c(0.0,1.0))
    x <- seq(from=1,length=10,to=(NameVal2-NameVal1+1))
    axis(1,at=x,labels=as.vector(round(x+NameVal1-1)))
    axis(2,at=pretty(0:1))
    box()
    
    # This draws the convergence criterion curves (prob and a.s.)
    critp <- critp[NameVal1:NameVal2]
    visualize.crit(critp,plotfunc=points,col="blue")
    critas[(K*nmax+1):nmax] <- NA
    critas <- critas[NameVal1:NameVal2]
    visualize.crit(critas,plotfunc=points,col="red")

    index <- nn-NameVal1+1
    if (index<1) index <- NA
    
    # This adds the little blue or red circle
    if (mode=="p") {text(index,critp[index],labels="O",col="blue")}
    if (mode=="as") {text(index,critas[index],labels="O",col="red")}

    # The title
    if (mode=="p") {title(bquote(textstyle(bold(atop("Criterion value for convergence in probability.","We have " ~ hat(p)[n] ==  .(round(critp[index]*M,3)) ~ "/" ~ .(M) == .(round(critp[index],3)) ~ "(based on all the M=" ~ .(M) ~ "sample paths).")))))}
    if (mode=="as") {title(bquote(textstyle(bold(atop("Criterion value for convergence almost surely.","We have " ~ hat(a)[n] ==  .(round(critas[index]*M,3)) ~ "/" ~ .(M) == .(round(critas[index],3)) ~ "(based on all the M=" ~ .(M) ~ "sample paths).")))))}
    title(sub=paste("n=",format(nn),sep=""))

    # Some useful information to help reading the plot
    mtext(bquote(bold(n[1] ~ "=")),side=1,line=-4,outer=TRUE,at=0.11)
    mtext(bquote(bold("=" ~ n[2])),side=1,line=-4,outer=TRUE,at=0.97)
    mtext(bquote(bold(hat(p)[n])),side=2,line=2.5,xpd=TRUE,col="blue",at=0.35)
    mtext("and",side=2,line=2.5,xpd=TRUE,at=0.5)
    mtext(bquote(bold(hat(a)[n])),side=2,line=2.5,xpd=TRUE,col="red",at=0.65)
  
  }
  
  
  movebar <- function(...) {
  # To move the grey bar
    
    tkfocus(tt)

    n <- as.numeric(tclvalue(sliderbarVal))
    if (n != nn) {
      
      nn <<- n
      tkrreplot(imgleft)
      tkrreplot(imgright)
    }
  }
  
  movewindow <- function(...) {
  # To move the (right and left) windows to the right or left
    
    tkfocus(tt)

    x <- as.numeric(tclvalue(tkget(zoombegin)))
    y <- as.numeric(tclvalue(tkget(zoomend)))
    
    tclvalue(zoombeginVar) <- tclvalue(sliderwindowVal)
    tclvalue(zoomendVar) <- as.numeric(tclvalue(sliderwindowVal))+y-x
    
    zoom()
    
  }



  zoom <- function(...) {
  # To redraw right and left plots from NameVal1 to NameVal2

    if (as.numeric(tclvalue(zoombeginVar))<1) {tclvalue(zoombeginVar) <<- 1;tkmessageBox(message="Bad value")}
    if (as.numeric(tclvalue(zoomendVar))>nmax) {tclvalue(zoomendVar) <<- nmax;tkmessageBox(message="Bad value")}
    if (as.numeric(tclvalue(zoombeginVar))>as.numeric(tclvalue(zoomendVar))) {tclvalue(zoombeginVar) <<- 1;tclvalue(zoomendVar) <<- nmax;tkmessageBox(message="Bad value")}
    
    NameVal1 <<- as.numeric(tclvalue(tkget(zoombegin)))
    NameVal2 <<- as.numeric(tclvalue(tkget(zoomend)))
    
    # This way, the moving grey bar is always visible
    if ( (as.numeric(tclvalue(sliderbarVal)) < NameVal1) )
    {
      if (mode == "p") {tclvalue(sliderbarVal) <- as.character(min(nmax,NameVal1));nn <<- min(nmax,NameVal1)}
      if (mode == "as") {tclvalue(sliderbarVal) <- as.character(min(K*nmax,NameVal1));nn <<- min(K*nmax,NameVal1)}
      
    }
    if ( (as.numeric(tclvalue(sliderbarVal)) > NameVal2))
      { tclvalue(sliderbarVal) <- as.character(NameVal2)
        nn <<- NameVal2
      }

    if (mode == "p") {tkconfigure(sliderbar,showvalue=TRUE,from=min(NameVal1,nmax),to=min(NameVal2,nmax))}
    if (mode == "as") {tkconfigure(sliderbar,showvalue=TRUE,from=min(NameVal1,K*nmax),to=min(NameVal2,K*nmax))}
    tkconfigure(sliderwindow,showvalue=TRUE,from=1,to=(nmax-NameVal2+NameVal1))
    tclvalue(sliderwindowVal) <- as.character(NameVal1)

    tkrreplot(imgleft)
    tkrreplot(imgright)
    
  }

# We initialize the variables used afterwards
  K <- 0.5
  M <- nrow(data)
  nmax <- ncol(data)
  zoombeginVar <- tclVar(1)
  zoomendVar <- tclVar(init="")
  tclvalue(zoomendVar) <- nmax
  sliderbarVal <- tclVar("1")
  sliderwindowVal <- tclVar(1)
  nn <- 1
  NameVal1 <- 1
  NameVal2 <- nmax
  mini.data <- min(data)
  maxi.data <- max(data)
  limhaute <- max(abs(range(data)))
  limbasse <- -max(abs(range(data)))
  if (mode=="p")  rbValue <- tclVar("prob")
  if (mode=="as")  rbValue <- tclVar("as")


  validzoombegin <- function(...) {

    myvarbegin <- as.numeric(tclvalue(zoombeginVar))
    myvarend <- as.numeric(tclvalue(zoomendVar))

    if (myvarbegin<1 | myvarbegin!=as.integer(myvarbegin) | myvarbegin>myvarend) {return(tclvalue(tclVar(FALSE))) }
    else {return(tclvalue(tclVar(TRUE)))}
  }

  invalidzoombegin <- function(...) {

    tclvalue(zoombeginVar) <- as.character(round(as.numeric(tclvalue(zoombeginVar))))

    if (as.numeric(tclvalue(zoombeginVar))<1) {tclvalue(zoombeginVar) <<- 1}
    if (as.numeric(tclvalue(zoombeginVar))>as.numeric(tclvalue(zoomendVar))) {tclvalue(zoombeginVar) <<- 1;tclvalue(zoomendVar) <<- nmax}

    tkconfigure(zoombegin,width="5",textvariable=zoombeginVar,width=9,validate="focus",invalidcommand=invalidzoombegin,validatecommand=validzoombegin)
    tkconfigure(zoomend,width="5",textvariable=zoomendVar,width=9,validate="focus",invalidcommand=invalidzoomend,validatecommand=validzoomend)
  }

  validzoomend <- function(...) {

    myvarbegin <- as.numeric(tclvalue(zoombeginVar))
    myvarend <- as.numeric(tclvalue(zoomendVar))

    if (myvarend>nmax | myvarend!=as.integer(myvarend) | myvarbegin>myvarend) {return(tclvalue(tclVar(FALSE))) }
    else {return(tclvalue(tclVar(TRUE)))}
  }

  invalidzoomend <- function(...) {

    tclvalue(zoomendVar) <- as.character(round(as.numeric(tclvalue(zoomendVar))))

    if (as.numeric(tclvalue(zoomendVar))>nmax) {tclvalue(zoomendVar) <<- nmax}
    if (as.numeric(tclvalue(zoombeginVar))>as.numeric(tclvalue(zoomendVar))) {tclvalue(zoombeginVar) <<- 1;tclvalue(zoomendVar) <<- nmax}
    
    tkconfigure(zoombegin,width="5",textvariable=zoombeginVar,width=9,validate="focus",invalidcommand=invalidzoombegin,validatecommand=validzoombegin)
    tkconfigure(zoomend,width="5",textvariable=zoomendVar,width=9,validate="focus",invalidcommand=invalidzoomend,validatecommand=validzoomend)
  }
  
  rb1func <- function(...) {

    tkconfigure(sliderbar,showvalue=TRUE,from=1,to=nmax)
    mode <<- "p"
    tkrreplot(imgleft)
    tkrreplot(imgright)
    tkconfigure(empty,text=" < =   nmax    ")
    tkconfigure(bartitle,text="Move bar : ")
    zoom()
    tkfocus(tt)
    
  }

  rb2func <- function(...) {

    tkconfigure(sliderbar,showvalue=TRUE,from=1,to=K*nmax)
    mode <<- "as"
    tkrreplot(imgleft)
    tkrreplot(imgright)
    tkconfigure(empty,text="<= K.nmax    ")
    tkconfigure(bartitle,text="Move block : ")
    zoom()
    tkfocus(tt)
  
  }
  
# We create all the Tcl/Tk widgets
  tt <- tktoplevel()
  tkwm.title(tt,"Seeing stochastic convergence in action")
  imgleft <- tkrplot(tt,plotleft,hscale=1.2,vscale=1)
  imgright <- tkrplot(tt,plotright,hscale=1.2,vscale=1)
  if (mode == "p") {bartitle <- tklabel(tt,text="Move bar : ")}
  if (mode == "as") {bartitle <- tklabel(tt,text="Move block : ")}
  if (mode == "p") {empty <- tklabel(tt,text=" < =   nmax    ")}
  if (mode == "as") {empty <- tklabel(tt,text="<= K.nmax    ")}
  if (mode == "p") {maxvalue <- tklabel(tt,text=paste(" <= nmax = ",format(nmax),sep=""))}
  if (mode == "as") {maxvalue <- tklabel(tt,text=paste(" <= nmax = ",format(nmax),sep=""))}
  if (mode == "p") {sliderbar <-tkscale(tt, command=movebar, from=1, to=nmax,showvalue=TRUE, variable=sliderbarVal,resolution=1, orient="horizontal",sliderlength=10,length=68)}
  if (mode == "as") {sliderbar <-tkscale(tt, command=movebar, from=1, to=K*nmax,showvalue=TRUE, variable=sliderbarVal,resolution=1, orient="horizontal",sliderlength=10,length=68)}
  slidewindowtitle <- tklabel(tt,text="Slide window : ")
  sliderwindow <- tkscale(tt, command=movewindow, from=1, to=1,showvalue=TRUE, variable=sliderwindowVal,resolution=1, orient="horizontal",sliderlength=10,length=68)
  titlezoombegin <- tklabel(tt,text="Zoom in/out from 1 <= n1 = ")
  zoombegin <- tkentry(tt,width="5",textvariable=zoombeginVar,width=9,validate="focus",invalidcommand=invalidzoombegin,validatecommand=validzoombegin)
  titlezoomend <- tklabel(tt,text="to n2 = ")
  zoomend <- tkentry(tt,width="5",textvariable=zoomendVar,width=12,validate="focus",invalidcommand=invalidzoomend,validatecommand=validzoomend)
  modetitle <- tklabel(tt,text="Mode of convergence: ")
  rb1 <- tkradiobutton(tt,command=rb1func,variable=rbValue,value="prob")
  rb2 <- tkradiobutton(tt,command=rb2func,variable=rbValue,value="as")
  rb1title <- tklabel(tt,text="Probability")
  rb2title <- tklabel(tt,text="Almost Sure")
  

# We put the Tcl/Tk widgets on the main window tt
  tkgrid(imgleft,imgright,columnspan=6)
  tkgrid(bartitle,sliderbar,empty,slidewindowtitle,sliderwindow,empty,modetitle,rb1title,rb1,rb2title,rb2,columnspan=1)
  tkgrid.configure(bartitle,slidewindowtitle,sticky="es")
  tkgrid(titlezoombegin,zoombegin,titlezoomend,zoomend,maxvalue,columnspan=1)
  tkgrid.configure(titlezoombegin,titlezoomend,rb1title,rb2title,sticky="e")
  tkgrid.configure(zoombegin,zoomend,sliderbar,rb1,rb2,sticky="w")
  tkgrid.configure(empty,sticky="s")
  tkbind(zoombegin, "<Return>",zoom)
  tkbind(zoomend, "<Return>",zoom)
  

  return(tt)
  
}

