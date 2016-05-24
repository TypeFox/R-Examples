######################################################################
#
# visualization.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/28/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains various routines related to graph visualization.
#
# Contents:
#   gplot
#   gplot.arrow
#   gplot.layout.adj
#   gplot.layout.circle
#   gplot.layout.circrand
#   gplot.layout.eigen
#   gplot.layout.fruchtermanreingold
#   gplot.layout.geodist
#   gplot.layout.hall
#   gplot.layout.kamadakawai
#   gplot.layout.mds
#   gplot.layout.princoord
#   gplot.layout.random
#   gplot.layout.rmds
#   gplot.layout.segeo
#   gplot.layout.seham
#   gplot.layout.spring
#   gplot.layout.springrepulse
#   gplot.layout.target
#   gplot.loop
#   gplot.target
#   gplot.vertex
#   gplot3d
#   gplot3d.arrow
#   gplot3d.layout.adj
#   gplot3d.layout.eigen
#   gplot3d.layout.fruchtermanreingold
#   gplot3d.layout.geodist
#   gplot3d.layout.hall
#   gplot3d.layout.kamadakawai
#   gplot3d.layout.mds
#   gplot3d.layout.princoord
#   gplot3d.layout.random
#   gplot3d.layout.rmds
#   gplot3d.layout.segeo
#   gplot3d.layout.seham
#   gplot3d.loop
#   plot.sociomatrix
#   sociomatrixplot
#
######################################################################


#gplot - Two-dimensional graph visualization
gplot<-function(dat,g=1,gmode="digraph",diag=FALSE,label=NULL,coord=NULL,jitter=TRUE,thresh=0,thresh.absval=TRUE,usearrows=TRUE,mode="fruchtermanreingold",displayisolates=TRUE,interactive=FALSE,interact.bycomp=FALSE,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,pad=0.2,label.pad=0.5,displaylabels=!is.null(label),boxed.labels=FALSE,label.pos=0,label.bg="white",vertex.enclose=FALSE,vertex.sides=NULL,vertex.rot=0,arrowhead.cex=1,label.cex=1,loop.cex=1,vertex.cex=1,edge.col=1,label.col=1,vertex.col=NULL,label.border=1,vertex.border=1,edge.lty=NULL,edge.lty.neg=2,label.lty=NULL,vertex.lty=1,edge.lwd=0,label.lwd=par("lwd"),edge.len=0.5,edge.curve=0.1,edge.steps=50,loop.steps=20,object.scale=0.01,uselen=FALSE,usecurve=FALSE,suppress.axes=TRUE,vertices.last=TRUE,new=TRUE,layout.par=NULL,...){
   #Turn the annoying locator bell off, and remove recursion limit
   bellstate<-options()$locatorBell
   expstate<-options()$expression
   on.exit(options(locatorBell=bellstate,expression=expstate))
   options(locatorBell=FALSE,expression=Inf)
   #Create a useful interval inclusion operator
   "%iin%"<-function(x,int) (x>=int[1])&(x<=int[2])
   #Extract the graph to be displayed and obtain its properties
   d<-as.edgelist.sna(dat,force.bipartite=(gmode=="twomode"))
   if(is.list(d))
     d<-d[[g]]
   n<-attr(d,"n")
   if(is.null(label)){
     if(displaylabels!=TRUE)
       displaylabels<-FALSE
     if(!is.null(attr(d,"vnames")))
       label<-attr(d,"vnames")
     else if((gmode=="twomode")&&(!is.null(attr(d,"bipartite"))))
       label<-c(paste("R",1:attr(d,"bipartite"),sep=""), paste("C",(attr(d,"bipartite")+1):n,sep=""))
     else{
       label<-1:n
     }
   }
   #Make adjustments for gmode, if required, and set up other defaults
   if(gmode=="graph"){
      usearrows<-FALSE
   } else if ((gmode=="twomode")&&(!is.null(attr(d,"bipartite")))) {
     #For two-mode graphs, make columns blue and 4-sided (versus 
     #red and 50-sided)
     #If defaults haven't been modified
     Rn <- attr(d,"bipartite")
     if (is.null(vertex.col)) vertex.col <- c(rep(2,Rn),rep(4,n-Rn))
     if (is.null(vertex.sides)) vertex.sides <- c(rep(50,Rn),rep(4,n-Rn))
   } 
   if (is.null(vertex.col)) vertex.col <- 2
   if (is.null(vertex.sides)) vertex.sides <- 50
   #Remove missing edges
   d<-d[!is.na(d[,3]),,drop=FALSE]
   #Set edge line types
   if (is.null(edge.lty)){    #If unset, assume 1 for pos, edge.lty.neg for neg
     edge.lty<-rep(1,NROW(d))
     if(!is.null(edge.lty.neg))   #If NULL, just ignore it
       edge.lty[d[,3]<0]<-edge.lty.neg
   }else{                     #If set, see what we were given...
     if(length(edge.lty)!=NROW(d)){      #Not specified per edge, so modify
       edge.lty<-rep(edge.lty,NROW(d))
       if(!is.null(edge.lty.neg))           #If given neg value, use it
         edge.lty[d[,3]<0]<-edge.lty.neg
     }else{                              #Might modify negative edges
       if(!is.null(edge.lty.neg))
         edge.lty[d[,3]<0]<-edge.lty.neg
     }
   }
   #Save a copy of d, in case values are needed
   d.raw<-d
   #Dichotomize d
   if(thresh.absval)
     d<-d[abs(d[,3])>thresh,,drop=FALSE] #Threshold by absolute value
   else
     d<-d[d[,3]>thresh,,drop=FALSE]      #Threshold by signed value
   attr(d,"n")<-n                    #Restore "n" to d
   #Determine coordinate placement
   if(!is.null(coord)){      #If the user has specified coords, override all other considerations
      x<-coord[,1]
      y<-coord[,2]
   }else{   #Otherwise, use the specified layout function
     layout.fun<-try(match.fun(paste("gplot.layout.",mode,sep="")),silent=TRUE)
     if(class(layout.fun)=="try-error")
       stop("Error in gplot: no layout function for mode ",mode)
     temp<-layout.fun(d,layout.par)
     x<-temp[,1]
     y<-temp[,2]
   }
   #Jitter the coordinates if need be
   if(jitter){
      x<-jitter(x)
      y<-jitter(y)
   }
   #Which nodes should we use?
   use<-displayisolates|(!is.isolate(d,ego=1:n))   
   #Deal with axis labels
   if(is.null(xlab))
     xlab=""
   if(is.null(ylab))
     ylab=""
   #Set limits for plotting region
   if(is.null(xlim))
     xlim<-c(min(x[use])-pad,max(x[use])+pad)  #Save x, y limits
   if(is.null(ylim))
     ylim<-c(min(y[use])-pad,max(y[use])+pad)
   xrng<-diff(xlim)          #Force scale to be symmetric
   yrng<-diff(ylim)
   xctr<-(xlim[2]+xlim[1])/2 #Get center of plotting region
   yctr<-(ylim[2]+ylim[1])/2
   if(xrng<yrng)
     xlim<-c(xctr-yrng/2,xctr+yrng/2)
   else
     ylim<-c(yctr-xrng/2,yctr+xrng/2)     
   baserad<-min(diff(xlim),diff(ylim))*object.scale*
     16/(4+n^(1/2))  #Set the "base radius," letting it shrink for large graphs
   #Create the base plot, if needed
   if(new){  #If new==FALSE, we add to the existing plot; else create a new one
     plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab=xlab,ylab=ylab,asp=1, axes=!suppress.axes,...)
   }
   #Fill out vertex vectors
   vertex.cex <- rep(vertex.cex,length=n)
   vertex.radius<-rep(baserad*vertex.cex,length=n)   #Create vertex radii
   vertex.sides <- rep(vertex.sides,length=n)
   vertex.border <- rep(vertex.border,length=n)
   vertex.col <- rep(vertex.col,length=n)
   vertex.lty <- rep(vertex.lty,length=n)
   vertex.rot <- rep(vertex.rot,length=n)
   loop.cex <- rep(loop.cex,length=n)
#AHM bugfix begin: label attributes weren't being filled out or restricted with [use]
   label.bg <- rep(label.bg,length=n)
   label.border <- rep(label.border,length=n)
   if(!is.null(label.lty)) {label.lty <- rep(label.lty,length=n)}
   label.lwd <- rep(label.lwd,length=n)
   label.col <- rep(label.col,length=n)
   label.cex <- rep(label.cex,length=n)
#AHM bugfix end: label attributes weren't being filled out or restricted with [use]
   #Plot vertices now, if desired
   if(!vertices.last){
#AHM feature start: enclose vertex polygons with circles (makes labels and arrows look better connected) 
     if(vertex.enclose) gplot.vertex(x[use],y[use],radius=vertex.radius[use], sides=50,col="#FFFFFFFF",border=vertex.border[use],lty=vertex.lty[use])
#AHM feature end: enclose vertex polygons with circles (makes labels and arrows look better connected) 
     gplot.vertex(x[use],y[use],radius=vertex.radius[use], sides=vertex.sides[use],col=vertex.col[use],border=vertex.border[use],lty=vertex.lty[use],rot=vertex.rot[use])
     }
   #Generate the edges and their attributes
   px0<-vector()   #Create position vectors (tail, head)
   py0<-vector()
   px1<-vector()
   py1<-vector()
   e.lwd<-vector() #Create edge attribute vectors
   e.curv<-vector()
   e.type<-vector()
   e.col<-vector()
   e.hoff<-vector() #Offset radii for heads
   e.toff<-vector() #Offset radii for tails
   e.diag<-vector() #Indicator for self-ties
   e.rad<-vector()  #Edge radius (only used for loops)
   if(NROW(d)>0){
     if(length(dim(edge.col))==2)   #Coerce edge.col/edge.lty to vector form
       edge.col<-edge.col[d[,1:2]]
     else
       edge.col<-rep(edge.col,length=NROW(d))
     if(length(dim(edge.lty))==2)
       edge.lty<-edge.lty[d[,1:2]]
     else
       edge.lty<-rep(edge.lty,length=NROW(d))
     if(length(dim(edge.lwd))==2){
       edge.lwd<-edge.lwd[d[,1:2]]
       e.lwd.as.mult<-FALSE
     }else{ 
       if(length(edge.lwd)==1)
         e.lwd.as.mult<-TRUE
       else
         e.lwd.as.mult<-FALSE
       edge.lwd<-rep(edge.lwd,length=NROW(d))
     }
     if(!is.null(edge.curve)){
       if(length(dim(edge.curve))==2){
         edge.curve<-edge.curve[d[,1:2]]
         e.curv.as.mult<-FALSE
       }else{ 
         if(length(edge.curve)==1)
           e.curv.as.mult<-TRUE
         else
           e.curv.as.mult<-FALSE
         edge.curve<-rep(edge.curve,length=NROW(d))
       }
     }else
       edge.curve<-rep(0,length=NROW(d))
     dist<-((x[d[,1]]-x[d[,2]])^2+(y[d[,1]]-y[d[,2]])^2)^0.5  #Get the inter-point distances for curves
     tl<-d*dist   #Get rescaled edge lengths
     tl.max<-max(tl)  #Get maximum edge length
     for(i in 1:NROW(d))
       if(use[d[i,1]]&&use[d[i,2]]){  #Plot edges for displayed vertices
         px0<-c(px0,as.double(x[d[i,1]]))  #Store endpoint coordinates
         py0<-c(py0,as.double(y[d[i,1]]))
         px1<-c(px1,as.double(x[d[i,2]]))
         py1<-c(py1,as.double(y[d[i,2]]))
         e.toff<-c(e.toff,vertex.radius[d[i,1]]) #Store endpoint offsets
         e.hoff<-c(e.hoff,vertex.radius[d[i,2]])
         e.col<-c(e.col,edge.col[i])    #Store other edge attributes
         e.type<-c(e.type,edge.lty[i])
         if(edge.lwd[i]>0){
           if(e.lwd.as.mult)
             e.lwd<-c(e.lwd,edge.lwd[i]*d.raw[i,3])
           else
             e.lwd<-c(e.lwd,edge.lwd[i])
         }else
           e.lwd<-c(e.lwd,1)
         e.diag<-c(e.diag,d[i,1]==d[i,2])  #Is this a loop?
         e.rad<-c(e.rad,vertex.radius[d[i,1]]*loop.cex[d[i,1]])
         if(uselen){   #Should we base curvature on interpoint distances?
           if(tl[i]>0){ 
             e.len<-dist[i]*tl.max/tl[i]
             e.curv<-c(e.curv,edge.len*sqrt((e.len/2)^2-(dist[i]/2)^2))
           }else{
             e.curv<-c(e.curv,0)   
           }
         }else{        #Otherwise, use prespecified edge.curve
           if(e.curv.as.mult)    #If it's a scalar, multiply by edge str
             e.curv<-c(e.curv,edge.curve[i]*dist[i])
           else
             e.curv<-c(e.curv,edge.curve[i])
         }
       }
     }
   #Plot loops for the diagonals, if diag==TRUE, rotating wrt center of mass
   if(diag&&(length(px0)>0)&&sum(e.diag>0)){  #Are there any loops present?
     gplot.loop(as.vector(px0)[e.diag],as.vector(py0)[e.diag], length=1.5*baserad*arrowhead.cex,angle=25,width=e.lwd[e.diag]*baserad/10,col=e.col[e.diag],border=e.col[e.diag],lty=e.type[e.diag],offset=e.hoff[e.diag],edge.steps=loop.steps,radius=e.rad[e.diag],arrowhead=usearrows,xctr=mean(x[use]),yctr=mean(y[use]))
   }
   #Plot standard (i.e., non-loop) edges
   if(length(px0)>0){  #If edges are present, remove loops from consideration
     px0<-px0[!e.diag] 
     py0<-py0[!e.diag]
     px1<-px1[!e.diag]
     py1<-py1[!e.diag]
     e.curv<-e.curv[!e.diag]
     e.lwd<-e.lwd[!e.diag]
     e.type<-e.type[!e.diag]
     e.col<-e.col[!e.diag]
     e.hoff<-e.hoff[!e.diag]
     e.toff<-e.toff[!e.diag]
     e.rad<-e.rad[!e.diag]
   }
   if(!usecurve&!uselen){   #Straight-line edge case
     if(length(px0)>0)
       gplot.arrow(as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1), length=2*baserad*arrowhead.cex,angle=20,col=e.col,border=e.col,lty=e.type,width=e.lwd*baserad/10,offset.head=e.hoff,offset.tail=e.toff,arrowhead=usearrows,edge.steps=edge.steps) #AHM edge.steps needed for lty to work
   }else{   #Curved edge case
     if(length(px0)>0){
       gplot.arrow(as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1), length=2*baserad*arrowhead.cex,angle=20,col=e.col,border=e.col,lty=e.type,width=e.lwd*baserad/10,offset.head=e.hoff,offset.tail=e.toff,arrowhead=usearrows,curve=e.curv,edge.steps=edge.steps)
     }
   }
   #Plot vertices now, if we haven't already done so
   if(vertices.last){
#AHM feature start: enclose vertex polygons with circles (makes labels and arrows look better connected) 
     if(vertex.enclose) gplot.vertex(x[use],y[use],radius=vertex.radius[use], sides=50,col="#FFFFFFFF",border=vertex.border[use],lty=vertex.lty[use])
#AHM feature end: enclose vertex polygons with circles (makes labels and arrows look better connected) 
     gplot.vertex(x[use],y[use],radius=vertex.radius[use], sides=vertex.sides[use],col=vertex.col[use],border=vertex.border[use],lty=vertex.lty[use],rot=vertex.rot[use])
     }
   #Plot vertex labels, if needed
   if(displaylabels&(!all(label==""))&(!all(use==FALSE))){
     if (label.pos==0){
       xhat <- yhat <- rhat <- rep(0,n) 
       #Set up xoff yoff and roff when we get odd vertices
       xoff <- x[use]-mean(x[use])
       yoff <- y[use]-mean(y[use])
       roff <- sqrt(xoff^2+yoff^2)
       #Loop through vertices
       for (i in (1:n)[use]){
         #Find all in and out ties that aren't loops
         ij <- unique(c(d[d[,2]==i&d[,1]!=i,1],d[d[,1]==i&d[,2]!=i,2]))
         ij.n <- length(ij)
         if (ij.n>0) {
           #Loop through all ties and add each vector to label direction
           for (j in ij){
             dx <- x[i]-x[j]
             dy <- y[i]-y[j]
             dr <- sqrt(dx^2+dy^2)
             xhat[i] <- xhat[i]+dx/dr
             yhat[i] <- yhat[i]+dy/dr
           }
           #Take the average of all the ties
           xhat[i] <- xhat[i]/ij.n
           yhat[i] <- yhat[i]/ij.n
           rhat[i] <- sqrt(xhat[i]^2+yhat[i]^2)
           if (rhat[i]!=0) { # normalize direction vector
             xhat[i] <- xhat[i]/rhat[i]
             yhat[i] <- yhat[i]/rhat[i]
           } else { #if no direction, make xhat and yhat away from center
             xhat[i] <- xoff[i]/roff[i]
             yhat[i] <- yoff[i]/roff[i]
           }
         } else { #if no ties, make xhat and yhat away from center
           xhat[i] <- xoff[i]/roff[i]
           yhat[i] <- yoff[i]/roff[i]
         }
         if (xhat[i]==0) xhat[i] <- .01 #jitter to avoid labels on points
         if (yhat[i]==0) yhat[i] <- .01
       }
       xhat <- xhat[use]
       yhat <- yhat[use]
     } else if (label.pos<5) {
       xhat <- switch(label.pos,0,-1,0,1)
       yhat <- switch(label.pos,-1,0,1,0)
     } else if (label.pos==6) {
       xoff <- x[use]-mean(x[use])
       yoff <- y[use]-mean(y[use])
       roff <- sqrt(xoff^2+yoff^2)
       xhat <- xoff/roff
       yhat <- yoff/roff
     } else {
       xhat <- 0
       yhat <- 0
     }
#AHM bugfix start: label attributes weren't being filled out or restricted with [use]
#     os<-par()$cxy*label.cex #AHM not used and now chokes on properly filled label.cex
     lw<-strwidth(label[use],cex=label.cex[use])/2
     lh<-strheight(label[use],cex=label.cex[use])/2
     if(boxed.labels){
       rect(x[use]+xhat*vertex.radius[use]-(lh*label.pad+lw)*((xhat<0)*2+ (xhat==0)*1),
         y[use]+yhat*vertex.radius[use]-(lh*label.pad+lh)*((yhat<0)*2+ (yhat==0)*1),
         x[use]+xhat*vertex.radius[use]+(lh*label.pad+lw)*((xhat>0)*2+ (xhat==0)*1),
         y[use]+yhat*vertex.radius[use]+(lh*label.pad+lh)*((yhat>0)*2+ (yhat==0)*1),
         col=label.bg[use],border=label.border[use],lty=label.lty[use],lwd=label.lwd[use])
     }
     text(x[use]+xhat*vertex.radius[use]+(lh*label.pad+lw)*((xhat>0)-(xhat<0)),
          y[use]+yhat*vertex.radius[use]+(lh*label.pad+lh)*((yhat>0)-(yhat<0)),
          label[use],cex=label.cex[use],col=label.col[use],offset=0)         
   }
#AHM bugfix end: label attributes weren't being filled out or restricted with [use]
	 #If interactive, allow the user to mess with things
   if((interactive|interact.bycomp)&&((length(x)>0)&&(!all(use==FALSE)))){ #AHM bugfix: interact.bycomp wouldn't fire without interactive also being set
     #Set up the text offset increment
     os<-c(0.2,0.4)*par()$cxy
     #Get the location for text messages, and write to the screen
     textloc<-c(min(x[use])-pad,max(y[use])+pad)
     tm<-"Select a vertex to move, or click \"Finished\" to end."
     tmh<-strheight(tm)
     tmw<-strwidth(tm)
     text(textloc[1],textloc[2],tm,adj=c(0,0.5)) #Print the initial instruction
     fm<-"Finished"
     finx<-c(textloc[1],textloc[1]+strwidth(fm))
     finy<-c(textloc[2]-3*tmh-strheight(fm)/2,textloc[2]-3*tmh+strheight(fm)/2)
     finbx<-finx+c(-os[1],os[1])
     finby<-finy+c(-os[2],os[2])
     rect(finbx[1],finby[1],finbx[2],finby[2],col="white")
     text(finx[1],mean(finy),fm,adj=c(0,0.5))
     #Get the click location
     clickpos<-unlist(locator(1))
     #If the click is in the "finished" box, end our little game.  Otherwise,
     #relocate a vertex and redraw.
     if((clickpos[1]%iin%finbx)&&(clickpos[2]%iin%finby)){
       cl<-match.call()                #Get the args of the current function
       cl$interactive<-FALSE           #Turn off interactivity
       cl$coord<-cbind(x,y)            #Set the coordinates
       cl$dat<-dat                     #"Fix" the data array
       return(eval(cl))     #Execute the function and return
     }else{
       #Figure out which vertex was selected
       clickdis<-sqrt((clickpos[1]-x[use])^2+(clickpos[2]-y[use])^2)
       selvert<-match(min(clickdis),clickdis)
       #Create usable labels, if the current ones aren't
       if(all(label==""))
         label<-1:n
       #Clear out the old message, and write a new one
       rect(textloc[1],textloc[2]-tmh/2,textloc[1]+tmw,textloc[2]+tmh/2, border="white",col="white")
       if (interact.bycomp) tm <- "Where should I move this component?"
       else tm<-"Where should I move this vertex?"
       tmh<-strheight(tm)
       tmw<-strwidth(tm)
       text(textloc[1],textloc[2],tm,adj=c(0,0.5))
       fm<-paste("Vertex",label[use][selvert],"selected")
       finx<-c(textloc[1],textloc[1]+strwidth(fm))
       finy<-c(textloc[2]-3*tmh-strheight(fm)/2,textloc[2]-3*tmh+ strheight(fm)/2)
       finbx<-finx+c(-os[1],os[1])
       finby<-finy+c(-os[2],os[2])
       rect(finbx[1],finby[1],finbx[2],finby[2],col="white")
       text(finx[1],mean(finy),fm,adj=c(0,0.5))
       #Get the destination for the new vertex
       clickpos<-unlist(locator(1))
       #Set the coordinates accordingly
       if (interact.bycomp) {
         dx <- clickpos[1]-x[use][selvert]
         dy <- clickpos[2]-y[use][selvert]             
         comp.mem <- component.dist(d,connected="weak")$membership
         same.comp <- comp.mem[use]==comp.mem[use][selvert]
         x[use][same.comp] <- x[use][same.comp]+dx
         y[use][same.comp] <- y[use][same.comp]+dy
       } else {
         x[use][selvert]<-clickpos[1]
         y[use][selvert]<-clickpos[2]
       }
       #Iterate (leaving interactivity on)
       cl<-match.call()                #Get the args of the current function
       cl$coord<-cbind(x,y)            #Set the coordinates
       cl$dat<-dat                     #"Fix" the data array
       return(eval(cl))     #Execute the function and return
     }
   }
   #Return the vertex positions, should they be needed
   invisible(cbind(x,y))
}


#gplot.arrow - Custom arrow-drawing method for gplot
gplot.arrow<-function(x0,y0,x1,y1,length=0.1,angle=20,width=0.01,col=1,border=1,lty=1,offset.head=0,offset.tail=0,arrowhead=TRUE,curve=0,edge.steps=50,...){
  if(length(x0)==0)   #Leave if there's nothing to do
    return;
  #Introduce a function to make coordinates for a single polygon
  make.coords<-function(x0,y0,x1,y1,ahangle,ahlen,swid,toff,hoff,ahead,curve,csteps,lty){ 
	if (lty=="blank"|lty==0) return(c(NA,NA)) #AHM leave if lty is "blank"
    slen<-sqrt((x0-x1)^2+(y0-y1)^2)  #Find the total length
#AHM begin code to fix csteps so all dashed lines look the same
	xlenin=(abs(x0-x1)/(par()$usr[2]-par()$usr[1]))*par()$pin[1]
	ylenin=(abs(y0-y1)/(par()$usr[4]-par()$usr[3]))*par()$pin[2]
	csteps=csteps*sqrt(xlenin^2+ylenin^2)
#AHM end code to fix csteps so all dashed lines look the same
#AHM begin code to decode lty (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash)
	if (is.character(lty)){
		lty <- switch (lty,blank=0,solid=1,dashed=2,dotted=3,dotdash=4,longdash=5,twodash=6,lty)
	} else {
		lty  <- as.character(lty)
    }
	if (is.na(as.integer(lty))) lty <- "10"
	if (as.integer(lty)<10) lty <- c("01","10","44", "13", "1343", "73", "2262")[as.integer(lty)+1]
#AHM end code to decode lty
    if(curve==0&lty=="10"){         #Straight, solid edges
      if(ahead){    
        coord<-rbind(                    #Produce a "generic" version w/head
          c(-swid/2,toff),
          c(-swid/2,slen-0.5*ahlen-hoff),
          c(-swid/2-ahlen*sin(ahangle),slen-ahlen*cos(ahangle)-hoff),
          c(0,slen-hoff),
          c(swid/2+ahlen*sin(ahangle),slen-ahlen*cos(ahangle)-hoff),
          c(swid/2,slen-0.5*ahlen-hoff),
          c(swid/2,toff),
          c(NA,NA)
        )
      }else{
        coord<-rbind(                    #Produce a "generic" version w/out head
          c(-swid/2,toff),
          c(-swid/2,slen-hoff),
          c(swid/2,slen-hoff),
          c(swid/2,toff),
          c(NA,NA)
        )
      }
    }else{             #Curved or non-solid edges (requires incremental polygons)
      theta<-atan2(y1-y0,x1-x0)  #Adjust curved arrows to make start/stop points meet at edge of polygon
      x0<-x0+cos(theta)*toff
      x1<-x1-cos(theta)*hoff
      y0<-y0+sin(theta)*toff
      y1<-y1-sin(theta)*hoff
      slen<-sqrt((x0-x1)^2+(y0-y1)^2)
#AHM begin toff/hoff bugfix and simplification of curve code (elimination of toff and hoff)
      if(ahead){    
        inc<-(0:csteps)/csteps
        coord<-rbind(
          cbind(
          -curve*(1-(2*(inc-0.5))^2)-swid/2,inc*(slen-ahlen*0.5)),
          c(-swid/2+ahlen*sin(-ahangle-(curve>0)*pi/16), slen-ahlen*cos(-ahangle-(curve>0)*pi/16)),
          c(0,slen),
          c(swid/2+ahlen*sin(ahangle-(curve>0)*pi/16), slen-ahlen*cos(ahangle-(curve>0)*pi/16)),
          cbind(-curve*(1-(2*(rev(inc)-0.5))^2)+swid/2,rev(inc)*(slen-ahlen*0.5)),
          c(NA,NA)
        )
      }else{
        inc<-(0:csteps)/csteps
        coord<-rbind(
          cbind(-curve*(1-(2*(inc-0.5))^2)-swid/2, inc*slen),
          cbind(-curve*(1-(2*(rev(inc)-0.5))^2)+swid/2, rev(inc)*slen),
          c(NA,NA)
        )
      }
    }
#AHM end bugfix and simplification of curve code
    theta<-atan2(y1-y0,x1-x0)-pi/2     #Rotate about origin
    rmat<-rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
    coord<-coord%*%rmat
    coord[,1]<-coord[,1]+x0            #Translate to (x0,y0)
    coord[,2]<-coord[,2]+y0
#AHM begin code to allow for lty other than 1
    if (lty!="10"){         #Straight, solid edges
      inc <- 1
	  lty.i <- 1
	  lty.n <- nchar(lty)
	  inc.solid=as.integer(substr(lty,lty.i,lty.i))
	  inc.blank=as.integer(substr(lty,lty.i+1,lty.i+1))
      coord.n <- dim(coord)[1]
      coord2 <- NULL
      while (inc<(csteps-inc.solid-inc.blank+1)) {
	    coord2 <- rbind(coord2,coord[inc:(inc+inc.solid),],
	    	coord[(coord.n-inc.solid-inc):(coord.n-inc),],c(NA,NA))
	    inc <- inc+inc.solid+inc.blank
	    lty.i=lty.i+2
        if (lty.i>lty.n) lty.i <- 1
	  }
	  if (inc<(coord.n-inc)) coord2 <- rbind(coord2,coord[inc:(coord.n-inc),],c(NA,NA))
      coord <- coord2
    }
  coord
  }
#AHM end code to allow for lty other than 1
  #"Stretch" the arguments
  n<-length(x0)
  angle<-rep(angle,length=n)/360*2*pi
  length<-rep(length,length=n)
  width<-rep(width,length=n)
  col<-rep(col,length=n)
  border<-rep(border,length=n)
  lty<-rep(lty,length=n)
  arrowhead<-rep(arrowhead,length=n)
  offset.head<-rep(offset.head,length=n)
  offset.tail<-rep(offset.tail,length=n)
  curve<-rep(curve,length=n)
  edge.steps<-rep(edge.steps,length=n)
  #Obtain coordinates
  coord<-vector()
  for(i in 1:n)  
    coord<-rbind(coord,make.coords(x0[i],y0[i],x1[i],y1[i],angle[i],length[i], width[i],offset.tail[i],offset.head[i],arrowhead[i],curve[i],edge.steps[i],lty[i]))
  coord<-coord[-NROW(coord),]
  #Draw polygons
  polygon(coord,col=col,border=border,...) #AHM no longer pass lty, taken care of internally.
}


#gplot.layout.adj - Layout method (MDS of inverted adjacency matrix) for gplot
gplot.layout.adj<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="invadj"
  layout.par$dist="none"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}


#gplot.layout.circle - Place vertices in a circular layout
gplot.layout.circle<-function(d,layout.par){
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-attr(d,"n")
  cbind(sin(2*pi*((0:(n-1))/n)),cos(2*pi*((0:(n-1))/n)))
}


#gplot.layout.circrand - Random circular layout for gplot
gplot.layout.circrand<-function(d,layout.par){ 
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$dist="uniang"
  gplot.layout.random(d,layout.par)
}


#gplot.layout.eigen - Place vertices based on the first two eigenvectors of
#an adjacency matrix
gplot.layout.eigen<-function(d,layout.par){     
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Determine the matrix to be used
  if(is.null(layout.par$var))
    vm<-d
  else
    vm<-switch(layout.par$var,
      symupper=symmetrize(d,rule="uppper"),
      symlower=symmetrize(d,rule="lower"),
      symstrong=symmetrize(d,rule="strong"),
      symweak=symmetrize(d,rule="weak"),
      user=layout.par$mat,
      raw=d
    )
  #Pull the eigenstructure
  e<-eigen(vm)
  if(is.null(layout.par$evsel))
    coord<-Re(e$vectors[,1:2])
  else
    coord<-switch(layout.par$evsel,
      first=Re(e$vectors[,1:2]),
      size=Re(e$vectors[,rev(order(abs(e$values)))[1:2]])
    )
  #Return the result
  coord
}


#gplot.layout.fruchtermanreingold - Fruchterman-Reingold layout routine for #gplot
gplot.layout.fruchtermanreingold<-function(d,layout.par){
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Provide default settings
  n<-attr(d,"n")
  if(is.null(layout.par$niter))
    niter<-500
  else
    niter<-layout.par$niter
  if(is.null(layout.par$max.delta))
    max.delta<-n
  else
    max.delta<-layout.par$max.delta
  if(is.null(layout.par$area))
    area<-n^2
  else
    area<-layout.par$area
  if(is.null(layout.par$cool.exp))
    cool.exp<-3
  else
    cool.exp<-layout.par$cool.exp
  if(is.null(layout.par$repulse.rad))
    repulse.rad<-area*log(n)
  else
    repulse.rad<-layout.par$repulse.rad
  if(is.null(layout.par$ncell))
    ncell<-ceiling(n^0.4)
  else
    ncell<-layout.par$ncell
  if(is.null(layout.par$cell.jitter))
    cell.jitter<-0.5
  else
    cell.jitter<-layout.par$cell.jitter
  if(is.null(layout.par$cell.pointpointrad))
    cell.pointpointrad<-0
  else
    cell.pointpointrad<-layout.par$cell.pointpointrad
  if(is.null(layout.par$cell.pointcellrad))
    cell.pointcellrad<-18
  else
    cell.pointcellrad<-layout.par$cell.pointcellrad
  if(is.null(layout.par$cellcellcellrad))
    cell.cellcellrad<-ncell^2
  else
    cell.cellcellrad<-layout.par$cell.cellcellrad
  if(is.null(layout.par$seed.coord)){
    tempa<-sample((0:(n-1))/n) #Set initial positions randomly on the circle
    x<-n/(2*pi)*sin(2*pi*tempa)
    y<-n/(2*pi)*cos(2*pi*tempa)
  }else{
    x<-layout.par$seed.coord[,1]
    y<-layout.par$seed.coord[,2]
  }
  #Symmetrize the network, just in case
  d<-symmetrize(d,rule="weak",return.as.edgelist=TRUE) 
  #Perform the layout calculation
  layout<-.C("gplot_layout_fruchtermanreingold_R", as.double(d), as.double(n), as.double(NROW(d)), as.integer(niter), as.double(max.delta), as.double(area), as.double(cool.exp), as.double(repulse.rad), as.integer(ncell), as.double(cell.jitter), as.double(cell.pointpointrad), as.double(cell.pointcellrad), as.double(cell.cellcellrad), x=as.double(x), y=as.double(y), PACKAGE="sna")
  #Return the result
  cbind(layout$x,layout$y)
}


#gplot.layout.geodist - Layout method (MDS of geodesic distances) for gplot
gplot.layout.geodist<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="geodist"
  layout.par$dist="none"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}


#gplot.layout.hall - Hall's layout method for gplot
gplot.layout.hall<-function(d,layout.par){
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-NROW(d)
  #Build the Laplacian matrix
  sd<-symmetrize(d)
  laplacian<--sd
  diag(laplacian)<-degree(sd,cmode="indegree")
  #Return the eigenvectors with smallest eigenvalues
  eigen(laplacian)$vec[,(n-1):(n-2)]
}


#gplot.layout.kamadakawai
gplot.layout.kamadakawai<-function(d,layout.par){
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-attr(d,"n")
  if(is.null(layout.par$niter)){
    niter<-1000
  }else
    niter<-layout.par$niter
  if(is.null(layout.par$sigma)){
    sigma<-n/4
  }else
    sigma<-layout.par$sigma
  if(is.null(layout.par$initemp)){
    initemp<-10
  }else
    initemp<-layout.par$initemp
  if(is.null(layout.par$coolexp)){
    coolexp<-0.99
  }else
    coolexp<-layout.par$coolexp
  if(is.null(layout.par$kkconst)){
    kkconst<-n^2
  }else
    kkconst<-layout.par$kkconst
  if(is.null(layout.par$edge.val.as.str))
    edge.val.as.str<-TRUE
  else
    edge.val.as.str<-layout.par$edge.val.as.str
  if(is.null(layout.par$elen)){
    d<-symmetrize(d,return.as.edgelist=TRUE)
    if(edge.val.as.str)
      d[,3]<-1/d[,3]
    elen<-geodist(d,ignore.eval=FALSE)$gdist
    elen[elen==Inf]<-max(elen[is.finite(elen)])*1.25
  }else
    elen<-layout.par$elen
  if(is.null(layout.par$seed.coord)){
    x<-rnorm(n,0,n/4)
    y<-rnorm(n,0,n/4)
  }else{
    x<-layout.par$seed.coord[,1]
    y<-layout.par$seed.coord[,2]
  }
  #Obtain locations
  pos<-.C("gplot_layout_kamadakawai_R",as.integer(n),as.integer(niter), as.double(elen),as.double(initemp),as.double(coolexp),as.double(kkconst),as.double(sigma), x=as.double(x),y=as.double(y), PACKAGE="sna")
  #Return to x,y coords
  cbind(pos$x,pos$y)
}


#gplot.layout.mds - Place vertices based on metric multidimensional scaling
#of a distance matrix
gplot.layout.mds<-function(d,layout.par){     
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Determine the raw inputs for the scaling
  if(is.null(layout.par$var))
    vm<-cbind(d,t(d))
  else
    vm<-switch(layout.par$var,
      rowcol=cbind(d,t(d)),
      col=t(d),
      row=d,
      rcsum=d+t(d),
      rcdiff=t(d)-d,
      invadj=max(d)-d,
      geodist=geodist(d,inf.replace=NCOL(d))$gdist,
      user=layout.par$vm
    )
  #If needed, construct the distance matrix
  if(is.null(layout.par$dist))
    dm<-as.matrix(dist(vm))
  else
    dm<-switch(layout.par$dist,
      euclidean=as.matrix(dist(vm)),
      maximum=as.matrix(dist(vm,method="maximum")),
      manhattan=as.matrix(dist(vm,method="manhattan")),
      canberra=as.matrix(dist(vm,method="canberra")),
      none=vm
    )
  #Transform the distance matrix, if desired
  if(is.null(layout.par$exp))
    dm<-dm^2
  else
    dm<-dm^layout.par$exp
  #Perform the scaling and return
  cmdscale(dm,2)
}


#gplot.layout.princoord - Place using the eigenstructure of the correlation 
#matrix among concatenated rows/columns (principal coordinates by position
#similarity)
gplot.layout.princoord<-function(d,layout.par){     
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Determine the vectors to be related
  if(is.null(layout.par$var))
    vm<-rbind(d,t(d))
  else
    vm<-switch(layout.par$var,
      rowcol=rbind(d,t(d)),
      col=d,
      row=t(d),
      rcsum=d+t(d),
      rcdiff=d-t(d),
      user=layout.par$vm
    )
  #Find the correlation/covariance matrix
  if(is.null(layout.par$cor)||layout.par$cor)
    cd<-cor(vm,use="pairwise.complete.obs")
  else    
    cd<-cov(vm,use="pairwise.complete.obs")
  cd<-replace(cd,is.na(cd),0)
  #Obtain the eigensolution
  e<-eigen(cd,symmetric=TRUE)
  x<-Re(e$vectors[,1])
  y<-Re(e$vectors[,2])
  cbind(x,y)
}


#gplot.layout.random - Random layout for gplot
gplot.layout.random<-function(d,layout.par){     
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-attr(d,"n")
  #Determine the distribution
  if(is.null(layout.par$dist))
    temp<-matrix(runif(2*n,-1,1),n,2)
  else if (layout.par$dist=="unif")
    temp<-matrix(runif(2*n,-1,1),n,2)
  else if (layout.par$dist=="uniang"){
    tempd<-rnorm(n,1,0.25)
    tempa<-runif(n,0,2*pi)
    temp<-cbind(tempd*sin(tempa),tempd*cos(tempa))
  }else if (layout.par$dist=="normal")
    temp<-matrix(rnorm(2*n),n,2)
  #Return the result
  temp
}


#gplot.layout.rmds - Layout method (MDS of euclidean row distances) for gplot
gplot.layout.rmds<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="row"
  layout.par$dist="euclidean"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}


#gplot.layout.segeo - Layout method (structural equivalence in geodesic 
#distances) for gplot
gplot.layout.segeo<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="geodist"
  layout.par$dist="euclidean"
  gplot.layout.mds(d,layout.par)
}


#gplot.layout.seham - Layout method (structural equivalence under Hamming
#metric) for gplot
gplot.layout.seham<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="rowcol"
  layout.par$dist="manhattan"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}


#gplot.layout.spring - Place vertices using a spring embedder
gplot.layout.spring<-function(d,layout.par){
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Set up the embedder params
  ep<-vector()
  if(is.null(layout.par$mass))  #Mass is in "quasi-kilograms"
    ep[1]<-0.1
  else
    ep[1]<-layout.par$mass
  if(is.null(layout.par$equil)) #Equilibrium extension is in "quasi-meters"
    ep[2]<-1
  else
    ep[2]<-layout.par$equil
  if(is.null(layout.par$k)) #Spring coefficient is in "quasi-Newtons/qm"
    ep[3]<-0.001
  else
    ep[3]<-layout.par$k
  if(is.null(layout.par$repeqdis)) #Repulsion equilibrium is in qm
    ep[4]<-0.1
  else
    ep[4]<-layout.par$repeqdis
  if(is.null(layout.par$kfr)) #Base coef of kinetic friction is in qn-qkg
    ep[5]<-0.01
  else
    ep[5]<-layout.par$kfr
  if(is.null(layout.par$repulse))
    repulse<-FALSE
  else
    repulse<-layout.par$repulse
  #Create initial condidions
  n<-dim(d)[1]
  f.x<-rep(0,n)       #Set initial x/y forces to zero
  f.y<-rep(0,n)
  v.x<-rep(0,n)       #Set initial x/y velocities to zero
  v.y<-rep(0,n)
  tempa<-sample((0:(n-1))/n) #Set initial positions randomly on the circle
  x<-n/(2*pi)*sin(2*pi*tempa)
  y<-n/(2*pi)*cos(2*pi*tempa)
  ds<-symmetrize(d,"weak")            #Symmetrize/dichotomize the graph
  kfr<-ep[5]                          #Set initial friction level
  niter<-1                            #Set the iteration counter
  #Simulate, with increasing friction, until motion stops    
  repeat{
    niter<-niter+1                    #Update the iteration counter
    dis<-as.matrix(dist(cbind(x,y)))  #Get inter-point distances
    #Get angles relative to the positive x direction
    theta<-acos(t(outer(x,x,"-"))/dis)*sign(t(outer(y,y,"-"))) 
    #Compute spring forces; note that we assume a base spring coefficient
    #of ep[3] units ("pseudo-Newtons/quasi-meter"?), with an equilibrium
    #extension of ep[2] units for all springs
    f.x<-apply(ds*cos(theta)*ep[3]*(dis-ep[2]),1,sum,na.rm=TRUE)
    f.y<-apply(ds*sin(theta)*ep[3]*(dis-ep[2]),1,sum,na.rm=TRUE)
    #If node repulsion is active, add a force component for this
    #as well.  We employ an inverse cube law which is equal in power
    #to the attractive spring force at distance ep[4]
    if(repulse){
      f.x<-f.x-apply(cos(theta)*ep[3]/(dis/ep[4])^3,1,sum,na.rm=TRUE)
      f.y<-f.y-apply(sin(theta)*ep[3]/(dis/ep[4])^3,1,sum,na.rm=TRUE)
    }
    #Adjust the velocities (assume a mass of ep[1] units); note that the
    #motion is roughly modeled on the sliding of flat objects across
    #a uniform surface (e.g., spring-connected cylinders across a table).
    #We assume that the coefficients of static and kinetic friction are
    #the same, which should only trouble you if you are under the 
    #delusion that this is a simulation rather than a graph drawing
    #exercise (in which case you should be upset that I'm not using
    #Runge-Kutta or the like!).
    v.x<-v.x+f.x/ep[1]         #Add accumulated spring/repulsion forces
    v.y<-v.y+f.y/ep[1]
    spd<-sqrt(v.x^2+v.y^2)     #Determine frictional forces
    fmag<-pmin(spd,kfr)  #We can't let friction _create_ motion!
    theta<-acos(v.x/spd)*sign(v.y)  #Calculate direction of motion
    f.x<-fmag*cos(theta)        #Decompose frictional forces
    f.y<-fmag*sin(theta)
    f.x[is.nan(f.x)]<-0         #Correct for any 0/0 problems
    f.y[is.nan(f.y)]<-0
    v.x<-v.x-f.x                #Apply frictional forces (opposing motion -
    v.y<-v.y-f.y                #note that mass falls out of equation)
    #Adjust the positions (yep, it's primitive linear updating time!)
    x<-x+v.x
    y<-y+v.y
    #Check for cessation of motion, and increase friction
    mdist<-mean(dis)
    if(all(v.x<mdist*1e-5)&&all(v.y<mdist*1e-5))
      break
    else
      kfr<-ep[5]*exp(0.1*niter)
  }
  #Return the result
  cbind(x,y)
}


#gplot.layout.springrepulse
gplot.layout.springrepulse<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$repulse<-TRUE
  gplot.layout.spring(d,layout.par)
}


#gplot.layout.target
gplot.layout.target<-function(d,layout.par){
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-NROW(d)
  if(is.null(layout.par$niter)){
    niter<-1000
  }else
    niter<-layout.par$niter
  if(is.null(layout.par$radii)){
    temp<-degree(d)
    offset<-min(sum(temp==max(temp))/(n-1),0.5)
    radii<-1-(temp-min(temp))/(diff(range(temp))+offset)
  }else
    radii<-layout.par$radii
  if(is.null(layout.par$minlen)){
    minlen<-0.05
  }else
    minlen<-layout.par$minlen
  if(is.null(layout.par$initemp)){
    initemp<-10
  }else
    initemp<-layout.par$initemp
  if(is.null(layout.par$coolexp)){
    coolexp<-0.99
  }else
    coolexp<-layout.par$coolexp
  if(is.null(layout.par$maxdelta)){
    maxdelta<-pi
  }else
    maxdelta<-layout.par$maxdelta
  if(is.null(layout.par$periph.outside)){
    periph.outside<-FALSE
  }else
    periph.outside<-layout.par$periph.outside
  if(is.null(layout.par$periph.outside.offset)){
    periph.outside.offset<-1.2
  }else
    periph.outside.offset<-layout.par$periph.outside.offset
  if(is.null(layout.par$disconst)){
    disconst<-1
  }else
    disconst<-layout.par$disconst
  if(is.null(layout.par$crossconst)){
    crossconst<-1
  }else
    crossconst<-layout.par$crossconst
  if(is.null(layout.par$repconst)){
    repconst<-1
  }else
    repconst<-layout.par$repconst
  if(is.null(layout.par$minpdis)){
    minpdis<-0.05
  }else
    minpdis<-layout.par$minpdis
  theta<-runif(n,0,2*pi)
  #Find core/peripheral vertices (in the sense of Brandes et al.)
  core<-apply(d&t(d),1,any)
  #Adjust radii if needed
  if(periph.outside)
    radii[!core]<-periph.outside.offset
  #Define optimal edge lengths
  elen<-abs(outer(radii,radii,"-"))
  elen[elen<minlen]<-(outer(radii,radii,"+")/sqrt(2))[elen<minlen]
  elen<-geodist(elen*d,inf.replace=n)$gdist
  #Obtain thetas
  pos<-.C("gplot_layout_target_R",as.integer(d),as.double(n), as.integer(niter),as.double(elen),as.double(radii),as.integer(core), as.double(disconst),as.double(crossconst),as.double(repconst), as.double(minpdis),as.double(initemp),as.double(coolexp),as.double(maxdelta), theta=as.double(theta),PACKAGE="sna")
  #Transform to x,y coords
  cbind(radii*cos(pos$theta),radii*sin(pos$theta))
}


#gplot.loop - Custom loop-drawing method for gplot
gplot.loop<-function(x0,y0,length=0.1,angle=10,width=0.01,col=1,border=1,lty=1,offset=0,edge.steps=10,radius=1,arrowhead=TRUE,xctr=0,yctr=0,...){
  if(length(x0)==0)   #Leave if there's nothing to do
    return;
  #Introduce a function to make coordinates for a single polygon
  make.coords<-function(x0,y0,xctr,yctr,ahangle,ahlen,swid,off,rad,ahead){
    #Determine the center of the plot
    xoff <- x0-xctr
    yoff <- y0-yctr
    roff <- sqrt(xoff^2+yoff^2)
    x0hat <- xoff/roff
    y0hat <- yoff/roff
    r0.vertex <- off
    r0.loop <- rad
    x0.loop <- x0hat*r0.loop
    y0.loop <- y0hat*r0.loop
    ang <- (((0:edge.steps)/edge.steps)*(1-(2*r0.vertex+0.5*ahlen*ahead)/ (2*pi*r0.loop))+r0.vertex/(2*pi*r0.loop))*2*pi+atan2(-yoff,-xoff)
    ang2 <- ((1-(2*r0.vertex)/(2*pi*r0.loop))+r0.vertex/(2*pi*r0.loop))*2*pi+ atan2(-yoff,-xoff)
    if(ahead){
      x0.arrow <- x0.loop+(r0.loop+swid/2)*cos(ang2)
      y0.arrow <- y0.loop+(r0.loop+swid/2)*sin(ang2)
      coord<-rbind(
        cbind(x0.loop+(r0.loop+swid/2)*cos(ang), 
          y0.loop+(r0.loop+swid/2)*sin(ang)),
        cbind(x0.arrow+ahlen*cos(ang2-pi/2),
          y0.arrow+ahlen*sin(ang2-pi/2)),
        cbind(x0.arrow,y0.arrow),
        cbind(x0.arrow+ahlen*cos(-2*ahangle+ang2-pi/2),
          y0.arrow+ahlen*sin(-2*ahangle+ang2-pi/2)),
        cbind(x0.loop+(r0.loop-swid/2)*cos(rev(ang)),
          y0.loop+(r0.loop-swid/2)*sin(rev(ang))),
        c(NA,NA)
      )
    }else{
      coord<-rbind(
        cbind(x0.loop+(r0.loop+swid/2)*cos(ang),
          y0.loop+(r0.loop+swid/2)*sin(ang)),
        cbind(x0.loop+(r0.loop-swid/2)*cos(rev(ang)),
          y0.loop+(r0.loop-swid/2)*sin(rev(ang))),
        c(NA,NA)
      )
    }
    coord[,1]<-coord[,1]+x0            #Translate to (x0,y0)
    coord[,2]<-coord[,2]+y0
    coord
  }
  #"Stretch" the arguments
  n<-length(x0)
  angle<-rep(angle,length=n)/360*2*pi
  length<-rep(length,length=n)
  width<-rep(width,length=n)
  col<-rep(col,length=n)
  border<-rep(border,length=n)
  lty<-rep(lty,length=n)
  rad<-rep(radius,length=n)
  arrowhead<-rep(arrowhead,length=n)
  offset<-rep(offset,length=n)
  #Obtain coordinates
  coord<-vector()
  for(i in 1:n)  
    coord<-rbind(coord,make.coords(x0[i],y0[i],xctr,yctr,angle[i],length[i], width[i],offset[i],rad[i],arrowhead[i]))
  coord<-coord[-NROW(coord),]
  #Draw polygons
  polygon(coord,col=col,border=border,lty=lty,...)
}


#gplot.target - Draw target diagrams using gplot
gplot.target<-function(dat,x,circ.rad=(1:10)/10,circ.col="blue",circ.lwd=1,circ.lty=3,circ.lab=TRUE,circ.lab.cex=0.75,circ.lab.theta=pi,circ.lab.col=1,circ.lab.digits=1,circ.lab.offset=0.025,periph.outside=FALSE,periph.outside.offset=1.2,...){
  #Transform x
  offset<-min(0.5,sum(x==max(x))/(length(x)-1))
  xrange<-diff(range(x))
  xmin<-min(x)
  x<-1-(x-xmin)/(xrange+offset)
  circ.val<-(1-circ.rad)*(xrange+offset)+xmin
  #Check for a layout.par, and set radii
  cl<-match.call()
  if(is.null(cl$layout.par))
    cl$layout.par<-list(radii=x)
  else
    cl$layout.par$radii<-x
  cl$layout.par$periph.outside<-periph.outside
  cl$layout.par$periph.outside.offset<-periph.outside.offset
  cl$x<-NULL
  cl$circ.rad<-NULL
  cl$circ.col<-NULL
  cl$circ.lwd<-NULL
  cl$circ.lty<-NULL
  cl$circ.lab.theta<-NULL
  cl$circ.lab.col<-NULL
  cl$circ.lab.cex<-NULL
  cl$circ.lab.digits<-NULL
  cl$circ.lab.offset<-NULL
  cl$periph.outside<-NULL
  cl$periph.outside.offset<-NULL
  cl$mode<-"target"
  cl$xlim=c(-periph.outside.offset,periph.outside.offset)
  cl$ylim=c(-periph.outside.offset,periph.outside.offset)
  cl[[1]]<-match.fun("gplot")
  #Perform the plotting operation
  coord<-eval(cl)
  #Draw circles
  if(length(circ.col)<length(x))
    circ.col<-rep(circ.col,length=length(x))
  if(length(circ.lwd)<length(x))
    circ.lwd<-rep(circ.lwd,length=length(x))
  if(length(circ.lty)<length(x))
    circ.lty<-rep(circ.lty,length=length(x))
  for(i in 1:length(circ.rad))
    segments(circ.rad[i]*sin(2*pi/100*(0:99)), circ.rad[i]*cos(2*pi/100*(0:99)),circ.rad[i]*sin(2*pi/100*(1:100)), circ.rad[i]*cos(2*pi/100*(1:100)),col=circ.col[i], lwd=circ.lwd[i],lty=circ.lty[i])
  if(circ.lab)
    text((circ.rad+circ.lab.offset)*cos(circ.lab.theta), (circ.rad+circ.lab.offset)*sin(circ.lab.theta), round(circ.val,digits=circ.lab.digits),cex=circ.lab.cex,col=circ.lab.col)
  #Silently return the resulting coordinates
  invisible(coord)
}


#gplot.vertex - Routine to plot vertices, using polygons
gplot.vertex<-function(x,y,radius=1,sides=4,border=1,col=2,lty=NULL,rot=0,...){
  #Introduce a function to make coordinates for a single polygon
  make.coords<-function(x,y,r,s,rot){
    ang<-(1:s)/s*2*pi+rot*2*pi/360
    rbind(cbind(x+r*cos(ang),y+r*sin(ang)),c(NA,NA))  
  }
  #Prep the vars
  n<-length(x)
  radius<-rep(radius,length=n)
  sides<-rep(sides,length=n)
  border<-rep(border,length=n)
  col<-rep(col,length=n)
  lty<-rep(lty,length=n)
  rot<-rep(rot,length=n)
  #Obtain the coordinates
  coord<-vector()
  for(i in 1:length(x))
    coord<-rbind(coord,make.coords(x[i],y[i],radius[i],sides[i],rot[i]))
  #Plot the polygons
  polygon(coord,border=border,col=col,lty=lty,...)
}


#gplot3d - Three-dimensional graph visualization
gplot3d<-function(dat,g=1,gmode="digraph",diag=FALSE,label=NULL,coord=NULL,jitter=TRUE,thresh=0,mode="fruchtermanreingold",displayisolates=TRUE,displaylabels=!missing(label),xlab=NULL,ylab=NULL,zlab=NULL,vertex.radius=NULL,absolute.radius=FALSE,label.col="gray50",edge.col="black",vertex.col=NULL,edge.alpha=1,vertex.alpha=1,edge.lwd=NULL,suppress.axes=TRUE,new=TRUE,bg.col="white",layout.par=NULL){
   #Require that rgl be loaded
   require(rgl)
   #Extract the graph to be displayed
   d<-as.edgelist.sna(dat,force.bipartite=(gmode=="twomode"))
   if(is.list(d))
     d<-d[[g]]
   n<-attr(d,"n")
   if(is.null(label)){
     if(displaylabels!=TRUE)
       displaylabels<-FALSE
     if(!is.null(attr(d,"vnames")))
       label<-attr(d,"vnames")
     else if((gmode=="twomode")&&(!is.null(attr(d,"bipartite"))))
       label<-c(paste("R",1:attr(d,"bipartite"),sep=""), paste("C",(attr(d,"bipartite")+1):n,sep=""))
     else{
       label<-1:n
     }
   }
   #Make adjustments for gmode, if required
   if(gmode=="graph"){
      usearrows<-FALSE
   }else if(gmode=="twomode"){
     if(is.null(vertex.col))
       vertex.col<-rep(c("red","blue"),times=c(attr(d,"bipartite"), n-attr(d,"bipartite")))
   }
   if(is.null(vertex.col))
     vertex.col<-"red"
   #Remove missing edges
   d<-d[!is.na(d[,3]),,drop=FALSE]
   #Save a copy of d, in case values are needed
   d.raw<-d
   #Dichotomize d
   d<-d[d[,3]>thresh,,drop=FALSE]
   attr(d,"n")<-n                    #Restore "n" to d
   #Determine coordinate placement
   if(!is.null(coord)){      #If the user has specified coords, override all other considerations
      x<-coord[,1]
      y<-coord[,2]
      z<-coord[,3]
   }else{   #Otherwise, use the specified layout function
     layout.fun<-try(match.fun(paste("gplot3d.layout.",mode,sep="")), silent=TRUE)
     if(class(layout.fun)=="try-error")
       stop("Error in gplot3d: no layout function for mode ",mode)
     temp<-layout.fun(d,layout.par)
     x<-temp[,1]
     y<-temp[,2]
     z<-temp[,3]
   }
   #Jitter the coordinates if need be
   if(jitter){
      x<-jitter(x)
      y<-jitter(y)
      z<-jitter(z)
   }
   #Which nodes should we use?
   use<-displayisolates|(!is.isolate(d,ego=1:n))   
   #Deal with axis labels
   if(is.null(xlab))
     xlab=""
   if(is.null(ylab))
     ylab=""
   if(is.null(zlab))
     zlab=""
   #Create the base plot, if needed
   if(new){  #If new==FALSE, we add to the existing plot; else create a new one
     rgl.clear()
     if(!suppress.axes)      #Plot axes, if desired
       rgl.bbox(xlab=xlab,ylab=ylab,zlab=zlab);
   }
   rgl.bg(color=bg.col)  
   #Plot vertices
   temp<-as.matrix(dist(cbind(x[use],y[use],z[use])))
   diag(temp)<-Inf
   baserad<-min(temp)/5
   if(is.null(vertex.radius)){
     vertex.radius<-rep(baserad,n)
   }else if(absolute.radius)
     vertex.radius<-rep(vertex.radius,length=n)
   else
     vertex.radius<-rep(vertex.radius*baserad,length=n)
   vertex.col<-rep(vertex.col,length=n)
   vertex.alpha<-rep(vertex.alpha,length=n)
   if(!all(use==FALSE))
     rgl.spheres(x[use],y[use],z[use],radius=vertex.radius[use], color=vertex.col[use], alpha=vertex.alpha[use])
   #Generate the edges and their attributes
   pt<-vector()   #Create position vectors (tail, head)
   ph<-vector()
   e.lwd<-vector() #Create edge attribute vectors
   e.col<-vector()
   e.alpha<-vector()
   e.diag<-vector() #Indicator for self-ties
   if(length(dim(edge.col))==2)   #Coerce edge.col/edge.lty to vector form
     edge.col<-edge.col[d[,1:2]]
   else
     edge.col<-rep(edge.col,length=NROW(d))
   if(is.null(edge.lwd)){
     edge.lwd<-0.5*apply(cbind(vertex.radius[d[,1]],vertex.radius[d[,2]]),1, min) + vertex.radius[d[,1]]*(d[,1]==d[,2])
   }else if(length(dim(edge.lwd))==2){
     edge.lwd<-edge.lwd[d[,1:2]]
   }else{
     if(edge.lwd==0)
       edge.lwd<-0.5*apply(cbind(vertex.radius[d[,1]],vertex.radius[d[,2]]),1, min) + vertex.radius[d[,1]]*(d[,1]==d[,2])
     else
       edge.lwd<-rep(edge.lwd,length=NROW(d))
   }
   if(length(dim(edge.alpha))==2){
     edge.alpha<-edge.alpha[d[,1:2]]
   }else{ 
     edge.alpha<-rep(edge.alpha,length=NROW(d))
   }
   for(i in 1:NROW(d))
     if(use[d[i,1]]&&use[d[i,2]]){    #Plot edges for displayed vertices
       pt<-rbind(pt,as.double(c(x[d[i,1]],y[d[i,1]],z[d[i,1]]))) #Store endpoint coordinates
       ph<-rbind(ph,as.double(c(x[d[i,2]],y[d[i,2]],z[d[i,2]])))
         e.col<-c(e.col,edge.col[i])    #Store other edge attributes
         e.alpha<-c(e.alpha,edge.alpha[i])
         e.lwd<-c(e.lwd,edge.lwd[i])
         e.diag<-c(e.diag,d[i,1]==d[i,2])  #Is this a loop?
       }
   m<-NROW(pt)  #Record number of edges
   #Plot loops for the diagonals, if diag==TRUE
   if(diag&&(m>0)&&sum(e.diag>0)){  #Are there any loops present?
     gplot3d.loop(pt[e.diag,],radius=e.lwd[e.diag],color=e.col[e.diag], alpha=e.alpha[e.diag])
   }
   #Plot standard (i.e., non-loop) edges
   if(m>0){  #If edges are present, remove loops from consideration
     pt<-pt[!e.diag,] 
     ph<-ph[!e.diag,]
     e.alpha<-e.alpha[!e.diag]
     e.lwd<-e.lwd[!e.diag]
     e.col<-e.col[!e.diag]
   }
   if(length(e.alpha)>0){
     gplot3d.arrow(pt,ph,radius=e.lwd,color=e.col,alpha=e.alpha)
   }
   #Plot vertex labels, if needed
   if(displaylabels&(!all(label==""))&(!all(use==FALSE))){
     rgl.texts(x[use]-vertex.radius[use],y[use],z[use],label[use], color=label.col)
   }
   #Return the vertex positions, should they be needed
   invisible(cbind(x,y,z))
}


#gplot3d.arrow- Draw a three-dimensional "arrow" from the positions in a to
#the positions in b, with specified characteristics.
gplot3d.arrow<-function(a,b,radius,color="white",alpha=1){
  #First, define an internal routine to make triangle coords
  make.coords<-function(a,b,radius){
    alen<-sqrt(sum((a-b)^2))
    xos<-radius*sin(pi/8)
    yos<-radius*cos(pi/8)
    basetri<-rbind(         #Create single offset triangle, pointing +z
      c(-xos,-yos,0), 
      c(0,0,alen), 
      c(xos,-yos,0)
    )
    coord<-vector()
    for(i in (1:8)/8*2*pi){  #Rotate about z axis to make arrow
      rmat<-rbind(c(cos(i),sin(i),0),c(-sin(i),cos(i),0), c(0,0,1))
      coord<-rbind(coord,basetri%*%rmat)
    }
    #Rotate into final angle (spherical coord w/+z polar axis...I know...)
    phi<--atan2(b[2]-a[2],a[1]-b[1])-pi/2
    psi<-acos((b[3]-a[3])/alen)
    coord<-coord%*%rbind(c(1,0,0),c(0,cos(psi),sin(psi)), c(0,-sin(psi),cos(psi)))
    coord<-coord%*%rbind(c(cos(phi),sin(phi),0),c(-sin(phi),cos(phi),0), c(0,0,1))
    #Translate into position
    coord[,1]<-coord[,1]+a[1]
    coord[,2]<-coord[,2]+a[2]
    coord[,3]<-coord[,3]+a[3]
    #Return the matrix
    coord
  }
  #Expand argument vectors if needed
  if(is.null(dim(a))){
    a<-matrix(a,ncol=3)
    b<-matrix(b,ncol=3)
  }  
  n<-NROW(a)
  radius<-rep(radius,length=n)
  color<-rep(color,length=n)
  alpha<-rep(alpha,length=n)
  #Obtain the joint coordinate matrix
  coord<-vector()
  for(i in 1:n)
    coord<-rbind(coord,make.coords(a[i,],b[i,],radius[i]))
  #Draw the triangles
  rgl.triangles(coord[,1],coord[,2],coord[,3],color=rep(color,each=24), alpha=rep(alpha,each=24))
}


#gplot3d.layout.adj - Layout method (MDS of inverse adjacencies) for gplot3d
gplot3d.layout.adj<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="invadj"
  layout.par$dist="none"
  layout.par$exp=1
  gplot3d.layout.mds(d,layout.par)
}


#gplot3d.layout.eigen - Place vertices based on the first three eigenvectors of
#an adjacency matrix
gplot3d.layout.eigen<-function(d,layout.par){     
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Determine the matrix to be used
  if(is.null(layout.par$var))
    vm<-d
  else
    vm<-switch(layout.par$var,
      symupper=symmetrize(d,rule="uppper"),
      symlower=symmetrize(d,rule="lower"),
      symstrong=symmetrize(d,rule="strong"),
      symweak=symmetrize(d,rule="weak"),
      user=layout.par$mat,
      raw=d
    )
  #Pull the eigenstructure
  e<-eigen(vm)
  if(is.null(layout.par$evsel))
    coord<-Re(e$vectors[,1:3])
  else
    coord<-switch(layout.par$evsel,
      first=Re(e$vectors[,1:3]),
      size=Re(e$vectors[,rev(order(abs(e$values)))[1:3]])
    )
  #Return the result
  coord
}


#gplot3d.layout.fruchtermanreingold - Fruchterman-Reingold layout method for
#gplot3d
gplot3d.layout.fruchtermanreingold<-function(d,layout.par){
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-attr(d,"n")
  #Provide default settings
  if(is.null(layout.par$niter))
    niter<-300
  else
    niter<-layout.par$niter
  if(is.null(layout.par$max.delta))
    max.delta<-n
  else
    max.delta<-layout.par$max.delta
  if(is.null(layout.par$volume))
    volume<-n^3
  else
    volume<-layout.par$volume
  if(is.null(layout.par$cool.exp))
    cool.exp<-3
  else
    cool.exp<-layout.par$cool.exp
  if(is.null(layout.par$repulse.rad))
    repulse.rad<-volume*n
  else
    repulse.rad<-layout.par$repulse.rad
  if(is.null(layout.par$seed.coord)){
    tempa<-runif(n,0,2*pi) #Set initial positions randomly on the sphere
    tempb<-runif(n,0,pi)
    x<-n*sin(tempb)*cos(tempa)
    y<-n*sin(tempb)*sin(tempa)
    z<-n*cos(tempb)
  }else{
    x<-layout.par$seed.coord[,1]
    y<-layout.par$seed.coord[,2]
    z<-layout.par$seed.coord[,3]
  }
  #Symmetrize the graph, just in case
  d<-symmetrize(d,return.as.edgelist=TRUE)
  #Set up positions
  #Perform the layout calculation
  layout<-.C("gplot3d_layout_fruchtermanreingold_R", as.double(d), as.integer(n), as.integer(NROW(d)), as.integer(niter), as.double(max.delta), as.double(volume), as.double(cool.exp), as.double(repulse.rad), x=as.double(x), y=as.double(y), z=as.double(z),PACKAGE="sna")
  #Return the result
  cbind(layout$x,layout$y,layout$z)
}


#gplot3d.layout.geodist - Layout method (MDS of geodesic distances) for gplot3d
gplot3d.layout.geodist<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="geodist"
  layout.par$dist="none"
  layout.par$exp=1
  gplot3d.layout.mds(d,layout.par)
}


#gplot3d.layout.hall - Hall's layout method for gplot3d
gplot3d.layout.hall<-function(d,layout.par){
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-NCOL(d)
  #Build the Laplacian matrix
  sd<-symmetrize(d)
  laplacian<--sd
  diag(laplacian)<-degree(sd,cmode="indegree")
  #Return the eigenvectors with smallest eigenvalues
  eigen(laplacian)$vec[,(n-1):(n-3)]
}


#gplot3d.layout.kamadakawai
gplot3d.layout.kamadakawai<-function(d,layout.par){
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-attr(d,"n")
  if(is.null(layout.par$niter)){
    niter<-1000
  }else
    niter<-layout.par$niter
  if(is.null(layout.par$sigma)){
    sigma<-n/4
  }else
    sigma<-layout.par$sigma
  if(is.null(layout.par$initemp)){
    initemp<-10
  }else
    initemp<-layout.par$initemp
  if(is.null(layout.par$coolexp)){
    coolexp<-0.99
  }else
    coolexp<-layout.par$coolexp
  if(is.null(layout.par$kkconst)){
    kkconst<-n^3
  }else
    kkconst<-layout.par$kkconst
  if(is.null(layout.par$edge.val.as.str))
    edge.val.as.str<-TRUE
  else
    edge.val.as.str<-layout.par$edge.val.as.str
  if(is.null(layout.par$elen)){
    d<-symmetrize(d,return.as.edgelist=TRUE)
    if(edge.val.as.str)
      d[,3]<-1/d[,3]
    elen<-geodist(d,ignore.eval=FALSE)$gdist
    elen[elen==Inf]<-max(elen[is.finite(elen)])*1.5
  }else
    elen<-layout.par$elen
  if(is.null(layout.par$seed.coord)){
    x<-rnorm(n,0,n/4)
    y<-rnorm(n,0,n/4)
    z<-rnorm(n,0,n/4)
  }else{
    x<-layout.par$seed.coord[,1]
    y<-layout.par$seed.coord[,2]
    z<-layout.par$seed.coord[,3]
  }
  #Obtain locations
  pos<-.C("gplot3d_layout_kamadakawai_R",as.double(n), as.integer(niter),as.double(elen),as.double(initemp),as.double(coolexp), as.double(kkconst),as.double(sigma),x=as.double(x),y=as.double(y), z=as.double(z),PACKAGE="sna")
  #Return to x,y coords
  cbind(pos$x,pos$y,pos$z)
}


#gplot3d.layout.mds - Place vertices based on metric multidimensional scaling
#of a distance matrix
gplot3d.layout.mds<-function(d,layout.par){     
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Determine the raw inputs for the scaling
  if(is.null(layout.par$var))
    vm<-cbind(d,t(d))
  else
    vm<-switch(layout.par$var,
      rowcol=cbind(d,t(d)),
      col=t(d),
      row=d,
      rcsum=d+t(d),
      rcdiff=t(d)-d,
      invadj=max(d)-d,
      geodist=geodist(d,inf.replace=NROW(d))$gdist,
      user=layout.par$vm
    )
  #If needed, construct the distance matrix
  if(is.null(layout.par$dist))
    dm<-as.matrix(dist(vm))
  else
    dm<-switch(layout.par$dist,
      euclidean=as.matrix(dist(vm)),
      maximum=as.matrix(dist(vm,method="maximum")),
      manhattan=as.matrix(dist(vm,method="manhattan")),
      canberra=as.matrix(dist(vm,method="canberra")),
      none=vm
    )
  #Transform the distance matrix, if desired
  if(is.null(layout.par$exp))
    dm<-dm^2
  else
    dm<-dm^layout.par$exp
  #Perform the scaling and return
  cmdscale(dm,3)
}


#gplot3d.layout.princoord - Place using the eigenstructure of the correlation 
#matrix among concatenated rows/columns (principal coordinates by position
#similarity)
gplot3d.layout.princoord<-function(d,layout.par){     
  d<-as.sociomatrix.sna(d)
  if(is.list(d))
    d<-d[[1]]
  #Determine the vectors to be related
  if(is.null(layout.par$var))
    vm<-rbind(d,t(d))
  else
    vm<-switch(layout.par$var,
      rowcol=rbind(d,t(d)),
      col=d,
      row=t(d),
      rcsum=d+t(d),
      rcdiff=d-t(d),
      user=layout.par$vm
    )
  #Find the correlation/covariance matrix
  if(is.null(layout.par$cor)||layout.par$cor)
    cd<-cor(vm,use="pairwise.complete.obs")
  else    
    cd<-cov(vm,use="pairwise.complete.obs")
  cd<-replace(cd,is.na(cd),0)
  #Obtain the eigensolution
  e<-eigen(cd,symmetric=TRUE)
  x<-Re(e$vectors[,1])
  y<-Re(e$vectors[,2])
  z<-Re(e$vectors[,3])
  cbind(x,y,z)
}


#gplot3d.layout.random - Layout method (random placement) for gplot3d
gplot3d.layout.random<-function(d,layout.par){     
  d<-as.edgelist.sna(d)
  if(is.list(d))
    d<-d[[1]]
  n<-attr(d,"n")
  #Determine the distribution
  if(is.null(layout.par$dist))
    temp<-matrix(runif(3*n,-1,1),n,3)
  else if (layout.par$dist=="unif")
    temp<-matrix(runif(3*n,-1,1),n,3)
  else if (layout.par$dist=="uniang"){
    tempd<-rnorm(n,1,0.25)
    tempa<-runif(n,0,2*pi)
    tempb<-runif(n,0,pi)
    temp<-cbind(tempd*sin(tempb)*cos(tempa),tempd*sin(tempb)*sin(tempa), tempd*cos(tempb))
  }else if (layout.par$dist=="normal")
    temp<-matrix(rnorm(3*n),n,3)
  #Return the result
  temp
}


#gplot3d.layout.rmds - Layout method (MDS of euclidean row distances) for 
#gplot3d
gplot3d.layout.rmds<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="row"
  layout.par$dist="euclidean"
  layout.par$exp=1
  gplot3d.layout.mds(d,layout.par)
}


#gplot3d.layout.segeo - Layout method (structural equivalence on geodesic
#distances) for gplot3d
gplot3d.layout.segeo<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="geodist"
  layout.par$dist="euclidean"
  gplot3d.layout.mds(d,layout.par)
}


#gplot3d.layout.seham - Layout method (structural equivalence under Hamming
#metric) for gplot3d
gplot3d.layout.seham<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="rowcol"
  layout.par$dist="manhattan"
  layout.par$exp=1
  gplot3d.layout.mds(d,layout.par)
}


#gplot3d.loop - Draw a three-dimensional "loop" at position a, with specified 
#characteristics.
gplot3d.loop<-function(a,radius,color="white",alpha=1){
  #First, define an internal routine to make triangle coords
  make.coords<-function(a,radius){
    coord<-rbind(
      cbind(
        a[1]+c(0,-radius/2,0), 
        a[2]+c(0,radius/2,radius/2), 
        a[3]+c(0,0,radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,-radius/2,0), 
        a[2]+c(0,radius/2,radius/2), 
        a[3]+c(0,0,-radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,radius/2,0), 
        a[2]+c(0,radius/2,radius/2), 
        a[3]+c(0,0,radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,radius/2,0), 
        a[2]+c(0,radius/2,radius/2), 
        a[3]+c(0,0,-radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,-radius/2,0), 
        a[2]+c(radius,radius/2,radius/2), 
        a[3]+c(0,0,radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,-radius/2,0), 
        a[2]+c(radius,radius/2,radius/2), 
        a[3]+c(0,0,-radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,radius/2,0), 
        a[2]+c(radius,radius/2,radius/2), 
        a[3]+c(0,0,radius/4),
        c(NA,NA,NA)
      ),
      cbind(
        a[1]+c(0,radius/2,0), 
        a[2]+c(radius,radius/2,radius/2), 
        a[3]+c(0,0,-radius/4),
        c(NA,NA,NA)
      )
    )
  }
  #Expand argument vectors if needed
  if(is.null(dim(a))){
    a<-matrix(a,ncol=3)
  }  
  n<-NROW(a)
  radius<-rep(radius,length=n)
  color<-rep(color,length=n)
  alpha<-rep(alpha,length=n)
  #Obtain the joint coordinate matrix
  coord<-vector()
  for(i in 1:n)
    coord<-rbind(coord,make.coords(a[i,],radius[i]))
  #Plot the triangles
  rgl.triangles(coord[,1],coord[,2],coord[,3],color=rep(color,each=24), alpha=rep(alpha,each=24))
}


#plot.sociomatrix - An odd sort of plotting routine; plots a matrix (e.g., a 
#Bernoulli graph density, or a set of adjacencies) as an image.  Very handy for 
#visualizing large valued matrices...
plot.sociomatrix<-function(x,labels=NULL,drawlab=TRUE,diaglab=TRUE,drawlines=TRUE,xlab=NULL,ylab=NULL,cex.lab=1,...){       
   #Begin preprocessing
   if((!(class(x)%in%c("matrix","array","data.frame")))||(length(dim(x))>2))
     x<-as.sociomatrix.sna(x)
   if(is.list(x))
     x<-x[[1]]
   #End preprocessing
   n<-dim(x)[1]
   o<-dim(x)[2]
   if(is.null(labels))
     labels<-list(NULL,NULL)
   if(is.null(labels[[1]])){  #Set labels, if needed
     if(is.null(rownames(x)))
       labels[[1]]<-1:dim(x)[1]
     else
       labels[[1]]<-rownames(x)
   }
   if(is.null(labels[[2]])){ 
     if(is.null(colnames(x)))
       labels[[2]]<-1:dim(x)[2]
     else
       labels[[2]]<-colnames(x)
   }
   d<-1-(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
   if(is.null(xlab))
     xlab<-""
   if(is.null(ylab))
     ylab<-""
   plot(1,1,xlim=c(0,o+1),ylim=c(n+1,0),type="n",axes=FALSE,xlab=xlab,ylab=ylab, ...)
   for(i in 1:n)
      for(j in 1:o)
         rect(j-0.5,i+0.5,j+0.5,i-0.5,col=gray(d[i,j]),xpd=TRUE, border=drawlines)
   rect(0.5,0.5,o+0.5,n+0.5,col=NA,xpd=TRUE)
   if(drawlab){
      text(rep(0,n),1:n,labels[[1]],cex=cex.lab)
      text(1:o,rep(0,o),labels[[2]],cex=cex.lab)
   }
   if((n==o)&(drawlab)&(diaglab))
      if(all(labels[[1]]==labels[[2]]))
         text(1:o,1:n,labels[[1]],cex=cex.lab)
}


#sociomatrixplot - an alias for plot.sociomatrix
sociomatrixplot<-plot.sociomatrix
