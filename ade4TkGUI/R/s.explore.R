#### ==============================
#### Explorateur de graphique ade4 
#### en tcl\tk                     
#### ==============================
#### par Stephane Dray             
#### 1/12/2004 Montreal            
#### ==============================

explore <- function(call.graph, scale.graph=1.3){
  
#  require(tcltk)
#  require(tkrplot)
#  
#  usrCoords <- NULL
#  usrCoords.small <- NULL
#  usrCoords.ori <- NULL
#  parPlotSize <- NULL
#  parPlotSize.small <- NULL
#  tl.select.tt2 <- tclVar()
#  rbValue.tt2 <- tclVar()
#  wh.include.row.tt2 <- tclVar()
#  rbValue.ttlabval <- tclVar()
#  xinf <- tclVar()
#  xsup <- tclVar()
#  yinf <- tclVar()
#  ysup <- tclVar()
#  xy1.rect <- NULL
#  zz <- NULL
#  
#  ## -----------------------------------------------------------------------
#  ##    Fonctionnement general
#  ## -----------------------------------------------------------------------
#
#  ## Pour enlever des warn inexpliques sur Linux
#  opar<-options(warn=-1)
#  on.exit(options(opar))
#  ## Recuperation de l'appel    
#  appel<-as.list(match.call())
#  if(!is.call(appel$call.graph))
#    appel$call.graph <- eval(appel$call.graph)
#  appel.plot <- as.call(appel$call.graph)
#  appel.list <- as.list(match.call(eval(appel.plot[[1]]),call=as.call(appel.plot)))
#  
#  type <- deparse(appel.list[[1]])
#  formals.s <- formals(eval(appel.plot[[1]]))
#  appel.list.small <- appel.list
#  
#  type.available <- c("s.arrow", "s.class","s.label","s.value")
#  
#  if(!(type%in%type.available)) {
#    stop(paste("Not implemented for",type,"\n"))
#  }
#  
#  ## Recuperation des parametres de l'appel  
#  ## differents type de s. , differents parametres
#
#  appel.list <- c(appel.list[1],lapply(appel.list[-1], eval, envir = sys.frame(0)))
#  appel.list <- c(appel.list,formals.s[which(is.na(pmatch(names(formals.s),names(appel.list))))])
#  appel.list <- c(appel.list[1],lapply(appel.list[-1],eval,envir=appel.list))
#  
#  if((appel.list$add.plot)) {
#    stop("set add.plot to FALSE")
#  }
#  
#  dftoplot <- appel.list$dfxy[,c(appel.list$xax,appel.list$yax)]
#  add.rectzoom.plot1 <- NULL
#  add.labels.plot1 <- NULL
#  add.find.plot1 <- NULL  
#  zf <- 100
#
#  if(type=="s.class"){
#    dftoplot <- cbind(dftoplot,fac=appel.list$fac)
#    clabel.ori <- appel.list$clabel
#    labelled.points <- NULL
#    label.ori <- row.names(dftoplot)
#    appel.list.small$clabel <- appel.list$clabel/3
#  }
#  
#  if(type=="s.arrow"){
#    label.ori <- appel.list$label 
#    clabel.ori <- appel.list$clabel
#    labelled.points <- NULL
#    origin.ori <- appel.list$origin
#    appel.list.small$clabel <- appel.list$clabel/3
#  }
#  
#  if(type=="s.label"){
#    label.ori <- appel.list$label
#    clabel.ori <- appel.list$clabel
#    appel.list.small$clabel <- appel.list$clabel/3
#    labelled.points <- NULL
#  }
#  
#  
#  if(type=="s.value"){
#    dftoplot <- cbind(dftoplot,Z=appel.list$z)
#    appel.list.small$csize <- appel.list$csize
#    appel.list.small$clegend <- appel.list$clegend/2
#  }
#  
#  
#  ## -----------------------------------------------------------------------
#  ##    Fonctions Graphiques generales
#  ## -----------------------------------------------------------------------
#  
#  
#  
#  plot.s.dens <- function(dens,x1,x2,type="horizontal"){
#    wh.x1 <- which.min(abs(dens$x-x1))
#    wh.x2 <- which.min(abs(dens$x-x2))
#    params <- par(mar = c(0, 0, 0, 0),bg="white")
#    if(type=="horizontal"){
#      plot(dens$x,dens$y,axes=FALSE,xlab="",ylab="",main="",ty='l',xlim=usrCoords.ori[1:2])
#      polygon(c(dens$x[wh.x1],dens$x[wh.x1:wh.x2],dens$x[wh.x2]),c(0,dens$y[wh.x1:wh.x2],0),col="red")
#    } else {
#      plot(dens$y,dens$x,axes=FALSE,xlab="",ylab="",main="",lwd=2,ty='l',ylim=usrCoords.ori[3:4])
#      polygon(c(0,dens$y[wh.x1:wh.x2],0),c(dens$x[wh.x1],dens$x[wh.x1:wh.x2],dens$x[wh.x2]),col="red")
#    }
#    box()
#    par(params)
#  }
#  
#  plot.s.densx <- function(){
#    plot.s.dens(densx,usrCoords[1],usrCoords[2])
#  }
#  
#  plot.s.densy <- function(){
#    plot.s.dens(densy,usrCoords[3],usrCoords[4],type="vertical")
#  }
#  
#  plot.s <- function() {
#    params <- par(bg="white",mar=rep(0.1,4))
#    eval(as.call(appel.list), sys.frame(0))
#    if(is.null(usrCoords.ori)) {usrCoords.ori <<- par("usr")}
#    if(!is.null(add.labels.plot1)){
#      eval(parse(text=add.labels.plot1), sys.frame(0))
#      if(as.numeric(tclvalue(Show.Value))==1){
#        cat(add.labels.plot1,"\n")
#      }
#    }
#    if(!is.null(add.rectzoom.plot1)){
#      eval(parse(text=add.rectzoom.plot1), sys.frame(0))
#      add.rectzoom.plot1 <<- NULL
#    }
#
#    parPlotSize <<- par("plt")
#    usrCoords   <<- par("usr")
#    par(params)
#    if(as.numeric(tclvalue(Show.Value))==1){
#      cat(deparse(as.call(appel.list)),"\n")
#    }
#    if(!is.null(add.find.plot1)){
#      eval(parse(text=add.find.plot1), sys.frame(0))
#      add.find.plot1 <<- NULL
#    }
#  }
#  
#  plot.s.small <- function () {
#    params=par(bg="white",mar = c(0.1, 0.1, 0.1, 0.1))
#    eval(as.call(appel.list.small), sys.frame(0))
#    rect(usrCoords[1],usrCoords[3],usrCoords[2],usrCoords[4],border="red")
#    parPlotSize.small <<- par("plt")
#    usrCoords.small   <<- par("usr")
#    par(params)
#  }
#  
#  plot.s.change.origin <- function(xnew,ynew){
#    xrange <- diff(usrCoords[1:2])
#    yrange <- diff(usrCoords[3:4])
#    appel.list$xlim[1] <<- xnew-xrange/2
#    appel.list$xlim[2] <<- xnew+xrange/2
#    appel.list$ylim[1] <<- ynew-yrange/2
#    appel.list$ylim[2] <<- ynew+yrange/2
#    tkrreplot(plot1)
#    tkrreplot(plot2)
#    tkrreplot(plot3)
#    tkrreplot(plot4)
#    refresh.textinfo()
#  }
#  
#  plot.s.change.origin.find <- function(xnew,ynew,label){
#    xrange <- diff(usrCoords[1:2])
#    yrange <- diff(usrCoords[3:4])
#    appel.list$xlim[1] <<- xnew-xrange/2
#    appel.list$xlim[2] <<- xnew+xrange/2
#    appel.list$ylim[1] <<- ynew-yrange/2
#    appel.list$ylim[2] <<- ynew+yrange/2
#    if (type =="s.label"){
#      add.find.plot1 <<- paste("scatterutil.eti(",xnew,",",ynew,",'",label,"',",clabel.ori,", boxes = TRUE, coul ='red')",sep="")
#    }
#    
#    
#    if (type =="s.class"){
#      add.find.plot1 <<- paste("scatterutil.eti(",xnew,",",ynew,",'",label,"',",clabel.ori,", boxes = TRUE, coul ='red')",sep="")
#    }
#    if (type =="s.arrow"){
#      add.find.plot1 <<- paste("scatterutil.eti.circ(",xnew,",",ynew,",'",label,"',",clabel.ori,", origin =c(",paste(origin.ori,collapse=","),"))",sep="")
#    }
#    tkrreplot(plot1)
#    tkrreplot(plot2)
#    tkrreplot(plot3)
#    tkrreplot(plot4)
#    refresh.textinfo()
#  }
#
#  
#  plot.s.refresh <- function() {
#    appel.list$xlim <<- usrCoords.ori[1:2]
#    appel.list$ylim <<- usrCoords.ori[3:4]
#    add.labels.plot1 <<- NULL
#    if(type=="s.label" | type=="s.arrow" | type=="s.class"){
#      appel.list$clabel <<- clabel.ori
#    }
#    tkrreplot(plot1)
#    tkrreplot(plot2)
#    tkrreplot(plot3)
#    tkrreplot(plot4)
#    tclvalue(zoomfactor)<<-100
#    refresh.textinfo()
#  }
#  
#  redo.plot.s.zoomfac <- function(...) {
#    if(as.numeric(tclObj(zoomfactor))!=zf){
#      zf <<- as.numeric(tclvalue(zoomfactor))
#      centerx <- mean(usrCoords[1:2])
#      centery <- mean(usrCoords[3:4])
#      xrange <- diff(usrCoords.ori[1:2])/2/zf*100
#      yrange <- diff(usrCoords.ori[3:4])/2/zf*100
#      appel.list$xlim <<- rep(centerx,2) + c(-xrange,xrange) 
#      appel.list$ylim <<- rep(centery,2) + c(-yrange,yrange)
#      tkrreplot(plot1)
#      tkrreplot(plot2)
#      tkrreplot(plot3)
#      tkrreplot(plot4)
#      refresh.textinfo()
#      
#    }
#  }
#  
#  refresh.textinfo<-function(){
#    wh.include.row <- which(dftoplot[,1]>usrCoords[1] & dftoplot[,1]<usrCoords[2] & dftoplot[,2]>usrCoords[3] & dftoplot[,2]<usrCoords[4],arr.ind=TRUE)
#    tkconfigure(txtdata, state="normal")
#    tkconfigure(txtlabrow, state="normal")
#    tkdelete(txtlabrow,"0.0","end")
#    tkdelete(txtdata,"0.0","end")
#    if(nrow(dftoplot)>0){
#      namesrowchar <- row.names(dftoplot)[wh.include.row]
#      tkinsert(txtlabrow,"end", paste(namesrowchar, collapse = "\n"))
#      tkinsert(txtdata,"end", paste(datachar[wh.include.row], collapse = "\n"))
#    }
#    tkconfigure(txtdata, state="disabled")
#    tkconfigure(txtlabrow, state="disabled")
#    
#  }
#  find.key <- function()
#    {
#      point.name <- tclvalue(point.found)
#      tkdelete(tl.find,0,"end")
#      indx.point <- grep(point.name,substr(name.sites,1,nchar(point.name)))
#      name.sites.new <<- name.sites[indx.point]
#        for (i in 1:length(name.sites.new))
#          {
#            tkinsert(tl.find,"end",name.sites.new[i])
#          }
#    }
#  
#  find.return <- function()
#    {
#      point.name <- tclvalue(point.found)
#      indx.point <- which(point.name==name.sites)
#      if(length(indx.point)>0) {
#      pointChoice <- as.vector(as.matrix(dftoplot[indx.point,c(1,2)]))
#      plot.s.change.origin.find(pointChoice[1],pointChoice[2],name.sites[indx.point])
#      }
#    }
#  
#  recenterplot <- function(){
#    indx.point <- which(name.sites.new[as.numeric(tkcurselection(tl.find))+1]==name.sites)
#    pointChoice <- as.vector(as.matrix(dftoplot[indx.point,c(1,2)]))
#    if (length(pointChoice)>0){
#      plot.s.change.origin.find(pointChoice[1],pointChoice[2],name.sites[indx.point])
#    }
#  }
#
#  nothing <- function(x,y){}
#  
#  getxy <- function(x,y){
#    xy1.rect <<- change.systeme.coordonnes(x,y)        
#  }
#  
#  change.systeme.coordonnes <- function(xClick,yClick,plot=plot1,parplosize=parPlotSize,usrcor=usrCoords){
#    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",plot)))
#    height <- as.numeric(tclvalue(tkwinfo("reqheight",plot)))
#    xMin <- parplosize[1] * width
#    xMax <- parplosize[2] * width
#    yMin <- parplosize[3] * height
#    yMax <- parplosize[4] * height
#    rangeX <- usrcor[2] - usrcor[1]
#    rangeY <- usrcor[4] - usrcor[3]
#    xClick <- as.numeric(xClick)
#    yClick <- as.numeric(yClick)
#    yClick <- height - yClick
#    xPlotCoord <- usrcor[1]+(xClick-xMin)*rangeX/(xMax-xMin)
#    yPlotCoord <- usrcor[3]+(yClick-yMin)*rangeY/(yMax-yMin)
#    return(c(xPlotCoord,yPlotCoord))
#  }
#  
#  findclosestpoint <- function(xClick,yClick,df){
#      squared.Distance <- (xClick-df[,1])^2 + (yClick-df[,2])^2
#      indexClosest <- which.min(squared.Distance)
#      return(indexClosest)
#    }
#  
#  identifyclosestpoint <- function(x,y){
#      newxy <- change.systeme.coordonnes(x,y)
#      index <- findclosestpoint(newxy[1],newxy[2],dftoplot)
#      msg <- paste("Id: ",row.names(dftoplot)[index],"\n",sep="")
#      msg <- paste(msg,paste(names(dftoplot),dftoplot[index,],sep=": ",collapse="\n"),collapse="",sep="")
#      tkmessageBox(title="Identify",message=msg,type="ok",icon="info")
#    }
#  
#
#  ## -----------------------------------------------------------------------
#  ##    Fonctions pour le mode Zoom 
#  ## -----------------------------------------------------------------------
#
#
#  
#  draw.rect <- function(x,y){
#    xy2.rect <- change.systeme.coordonnes(x,y)
#    add.rectzoom.plot1 <<- paste("rect(",xy1.rect[1],",",xy1.rect[2],",",xy2.rect[1],",",xy2.rect[2],",border=\"red\",lty=3)",sep="")
#    tkrreplot(plot1)        
#  }
#  
#  zoom.rect <- function(x,y){
#    xy2.rect <- change.systeme.coordonnes(x,y)
#    xx <- c(xy1.rect[1],xy2.rect[1])
#    yy <- c(xy1.rect[2],xy2.rect[2])
#    if((abs(diff(xx)/max(xx)))<1e-3 & (abs(diff(yy)/max(yy)))<1e-3){
#      ## si zone trop petite, zoomer et centrer sur le point
#      if(as.numeric(tclvalue(zoomfactor))<=(maxzoom-100)){
#        tclvalue(zoomfactor) <<- as.numeric(tclvalue(zoomfactor))+100
#        redo.plot.s.zoomfac()
#      } 
#      plot.s.change.origin(xy1.rect[1],xy1.rect[2])  
#      
#    } else {
#      zoom.usr <- diff(usrCoords.ori[1:2])*100/max(c(diff(xx),diff(yy)))
#      if(zoom.usr < maxzoom){
#        realzoom <- seq(50,maxzoom, by = 50)[which.min(abs(seq(50,maxzoom, by = 50)-zoom.usr))]
#        tclvalue(zoomfactor) <<- realzoom
#      } else {
#        tclvalue(zoomfactor) <<- maxzoom
#      }
#      redo.plot.s.zoomfac()
#      plot.s.change.origin(mean(c(xy1.rect[1],xy2.rect[1])), mean(c(xy1.rect[2],xy2.rect[2])))  
#    }  
#  }
#
#  ## -----------------------------------------------------------------------
#  ##    Fonctions pour le mode Label
#  ## -----------------------------------------------------------------------
#  
#  addonelabel <- function(x,y) {
#      newxy <- change.systeme.coordonnes(x,y)
#      index <- findclosestpoint(newxy[1],newxy[2],dftoplot)
#      msg <- paste("Id: ",row.names(dftoplot)[index],"\n",sep="")
#      msg <- paste("Add label for this point ?\n\n",msg,sep="",collapse="")
#      msg <- paste(msg,paste(names(dftoplot),dftoplot[index,],sep=": ",collapse="\n"),collapse="",sep="")
#      mbval<- tkmessageBox(title="Label",
#                           message=msg,type="yesno",icon="question",default="yes")
#      if (tclvalue(mbval)=="yes") {
#        labelled.points <<- unique(c(labelled.points,index))
#        if (type =="s.label"){
#          appel.list$clabel <<- 0
#          add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                     paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#        }
#        
#        
#        if (type =="s.class"){
#          add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                     paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#        }
#        if (type =="s.arrow"){
#          appel.list$clabel <<- 0
#          add.labels.plot1 <<- paste("scatterutil.eti.circ(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", origin =c(",paste(origin.ori,collapse=","),"))",sep="")
#        }
#        tkrreplot(plot1) 
#      }
#    }
#  
#  removelabel <- function()
#    {   labelled.points <<- NULL
#        add.labels.plot1 <<- NULL
#        if(type!="s.class" ) appel.list$clabel <<- 0
#        if (!(is.null(appel.list$cpoint))){
#          if (appel.list$cpoint == 0) appel.list$cpoint <<- 1
#        }
#        tkrreplot(plot1) 
#      }
#  
#  addalllabel <- function()
#    {   labelled.points <<- NULL
#        add.labels.plot1 <<- NULL
#        appel.list$clabel <<- clabel.ori
#        if (type =="s.class"){
#          labelled.points <<- 1:nrow(dftoplot)
#          add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                     paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#        }
#        tkrreplot(plot1) 
#      }
#  
#  addselectlabel <- function()
#    {   
#      okbutton.f1.tt2 <- function(){
#        cursel<-(as.numeric(tkcurselection(tl.select.tt2))+1)
#        if(length(cursel)==0) {tkdestroy(tt2)}
#        else {
#          if(tclvalue(rbValue.tt2)=="current"){labelled.points <<- unique(c(labelled.points,wh.include.row.tt2[cursel]))}
#          if(tclvalue(rbValue.tt2)=="whole"){labelled.points <<- unique(c(labelled.points,cursel))}
#          
#          if (type =="s.label"){
#            appel.list$clabel <<- 0
#            add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                       paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                       clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#          }
#          
#          if (type =="s.class"){
#            add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                       paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                       clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#          }
#          
#          
#          if (type =="s.arrow"){
#            appel.list$clabel <<- 0
#            add.labels.plot1 <<- paste("scatterutil.eti.circ(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                       paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                       clabel.ori,", origin =c(",paste(origin.ori,collapse=","),"))",sep="")
#          }
#          tkrreplot(plot1)
#          tkdestroy(tt2) 
#        }
#      }        
#      change.df.affichage <- function(){
#        tkdelete(tl.select.tt2,0,"end")
#        if(tclvalue(rbValue.tt2)=="whole"){
#          name.sites.df <- row.names(dftoplot)
#        }
#        else if(tclvalue(rbValue.tt2)=="current"){
#          wh.include.row.tt2 <<- which(dftoplot[,1]>usrCoords[1] & dftoplot[,1]<usrCoords[2] & dftoplot[,2]>usrCoords[3] & dftoplot[,2]<usrCoords[4],arr.ind=TRUE)
#          name.sites.df <- row.names(dftoplot)[wh.include.row.tt2]
#        }
#        if(length(name.sites)>0){    
#          for (i in 1:length(name.sites))
#            {
#              tkinsert(tl.select.tt2,"end",name.sites[i])
#            }
#        }
#      }           
#      
#      tt2 <- tktoplevel()
#      tkwm.title(tt2,"Add selected labels")
#      rb1.tt2 <- tkradiobutton(tt2,text="Current plot")
#      rb2.tt2 <- tkradiobutton(tt2,text="Whole plot")
#      rbValue.tt2 <<- tclVar("current")
#      tkconfigure(rb1.tt2,variable=rbValue.tt2,value="current",command=change.df.affichage)
#      tkconfigure(rb2.tt2,variable=rbValue.tt2,value="whole",command=change.df.affichage)
#      tkpack(tklabel(tt2,text="Select from:"))
#      tkpack(rb1.tt2,anchor="w")
#      tkpack(rb2.tt2,anchor="w")
#      scr2.tt2 <- tkscrollbar(tt2, repeatinterval=5,
#                              command=function(...)tkyview(tl.select.tt2,...))
#      tl.select.tt2<<-tklistbox(tt2,height=15,selectmode="extended",yscrollcommand=function(...)tkset(scr2.tt2,...),background="white")
#      
#      bt.cancel.tt2 <- tkbutton(tt2, text="  Cancel  ",command=function() tkdestroy(tt2))
#      bt.ok.tt2 <- tkbutton(tt2, text="  OK  ",command=okbutton.f1.tt2)
#      tkpack(bt.cancel.tt2,bt.ok.tt2,side="bottom")
#      tkpack(tl.select.tt2, side="left", fill="both", expand=TRUE)
#      tkpack(scr2.tt2, side="right", fill="y")
#      change.df.affichage()
#      
#    }
#  
#  labelbyvalues <- function(){
#
#    ok.button.labval <- function(){
#      
#      cursel <- which(dftoplot[,1] >= as.numeric(tclvalue(xinf)) & dftoplot[,1] <= as.numeric(tclvalue(xsup)) & dftoplot[,2] >= as.numeric(tclvalue(yinf)) & dftoplot[,2] <= as.numeric(tclvalue(ysup)))
#      if(length(cursel)==0) {tkdestroy(tt.labval)}
#      else {
#        
#        if(tclvalue(rbValue.ttlabval)=="add"){labelled.points <<- unique(c(labelled.points,cursel))}
#        if(tclvalue(rbValue.ttlabval)=="new"){labelled.points <<- cursel}
#        
#        if (type =="s.label"){
#          appel.list$clabel <<- 0
#          add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                     paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#          
#        }
#        
#        if (type =="s.class"){
#          add.labels.plot1 <<- paste("scatterutil.eti(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                     paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", boxes = TRUE, coul =rep(1,",paste(length(labelled.points)),"))",sep="")
#        }
#        
#        
#        if (type =="s.arrow"){
#          appel.list$clabel <<- 0
#          add.labels.plot1 <<- paste("scatterutil.eti.circ(c(",paste(dftoplot[labelled.points,1],collapse=","),"),c(",
#                                     paste(dftoplot[labelled.points,2],collapse=","),"),"," c('", paste(label.ori[labelled.points],collapse="','"),"'),",
#                                     clabel.ori,", origin =c(",paste(origin.ori,collapse=","),"))",sep="")
#        }
#        tkrreplot(plot1)
#        
#      }
#    }
#    
#    xinf <<- tclVar(usrCoords.ori[1])
#    xsup <<- tclVar(usrCoords.ori[2])
#    yinf <<- tclVar(usrCoords.ori[3])
#    ysup <<- tclVar(usrCoords.ori[4])
#    
#    tt.labval <- tktoplevel()
#    tt.x <- tkframe(tt.labval)
#    tt.y <- tkframe(tt.labval)
#    tkwm.title(tt.labval,"Add labels by values")
#    entry.xinf <-tkentry(tt.x,width="20",textvariable=xinf,background="white")
#    entry.xsup <-tkentry(tt.x,width="20",textvariable=xsup,background="white")
#    entry.yinf <-tkentry(tt.y,width="20",textvariable=yinf,background="white")
#    entry.ysup <-tkentry(tt.y,width="20",textvariable=ysup,background="white")
#    rb1.ttlabval <- tkradiobutton(tt.labval,text="add to current")
#    rb2.ttlabval <- tkradiobutton(tt.labval,text="new")
#    rbValue.ttlabval <<- tclVar("new")
#    tkconfigure(rb1.ttlabval,variable=rbValue.ttlabval,value="add")
#    tkconfigure(rb2.ttlabval,variable=rbValue.ttlabval,value="new")
#    tkpack(tklabel(tt.labval,text="Selection :"))
#    tkpack(rb1.ttlabval,anchor="w")
#    tkpack(rb2.ttlabval,anchor="w")
#    bt.cancel.tt.labval <- tkbutton(tt.labval, text="  Dismiss  ",command=function() {tkdestroy(tt.labval)})
#    bt.ok.tt.labval <- tkbutton(tt.labval, text="  OK  ",command=ok.button.labval)
#    tkpack(entry.xinf,tklabel(tt.x,text=" < x < "),entry.xsup,side="left")
#    tkpack(entry.yinf,tklabel(tt.y,text=" < y < "),entry.ysup,side="left",anchor="w")
#    tkpack(tt.x)
#    tkpack(tt.y)
#    tkpack(bt.cancel.tt.labval,side="left")
#    tkpack(bt.ok.tt.labval,side="right")
#    
#    
#  }
#  
#
#  ## -----------------------------------------------------------------------
#  ##    Fonctions generales associees a plot1
#  ## -----------------------------------------------------------------------
#
#  copytoclip <- function() {tkrreplot(plot1)}
#  
#  outgraph.explore <- function()
#    {
#      ##
#      ## Main dialog window with title and frames
#      ##
#      tf <- tktoplevel()
#      tkwm.title(tf,"Save graphic")
#      ##
#      ## Frames
#      ##
#      frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
#      frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
#      frame3 <- tkframe(tf, relief="groove", borderwidth=2)	
#      devframe <- tkframe(frame2, relief="groove", borderwidth=2)
#      ##
#      ## Tcl/Tk variables
#      ##
#      done <- tclVar(0)
#      formatvar <- tclVar(1)
#      widthvar <- tclVar(6)
#      heightvar <- tclVar(6)
#      ##
#      ## Save function
#      ##
#      savefic <- function(formatvar, widthvar, heightvar)
#	{
#          outform <- tclvalue(formatvar)
#          width <- as.numeric(tclvalue(widthvar))
#          height <- as.numeric(tclvalue(heightvar))
#          
#          if (outform == 1) { # postcript
#            filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.ps", defaultextension=".ps",
#                                               title="Save graph...", filetypes="{PostScript {.ps .eps}} {{All Files} {*.*}}"))
#            if (filename != "") {
#              postscript(file=filename, width=width, height=height)
#            }
#          } else if (outform == 2) { # pdf
#            filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.pdf", defaultextension=".pdf",
#                                               title="Save graph...", filetypes="{PDF {.pdf}} {{All Files} {*.*}}"))
#            if (filename != "") {
#              pdf(file=filename, width=width, height=height)
#            }
#          } else if (outform == 3) { # pictex
#            filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.tex", defaultextension=".tex",
#                                               title="Save graph...", filetypes="{PicTeX {.tex}} {{All Files} {*.*}}"))
#            if (filename != "") {
#              pictex(file=filename, width=width, height=height)
#            }
#          } else if (outform == 4) { # xfig
#            filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.fig", defaultextension=".fig",
#                                               title="Save graph...", filetypes="{XFig {.fig}} {{All Files} {*.*}}"))
#            if (filename != "") {
#              xfig(file=filename, width=width, height=height)
#            }
#          } else if (outform == 5) { # png
#            filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.png", defaultextension=".png",
#                                               title="Save graph...", filetypes="{PNG {.png}} {{All Files} {*.*}}"))
#            if (filename != "") {
#              png(filename=filename, width=width, height=height)
#            }
#          } else if (outform == 6) { # jpeg
#            filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.jpeg", defaultextension=".jpeg",
#                                               title="Save graph...", filetypes="{JPEG {.jpeg .jpg}} {{All Files} {*.*}}"))
#            if (filename != "") {
#              jpeg(filename=filename, width=width, height=height)
#            }
#          }
#          plot.s()
#          dev.off()
#          tkdestroy(tf)
#	}
#      ##
#      ## Frames setup
#      ##
#      tkgrid(tklabel(tf,text="Save current graphic", font="Times 18"), columnspan=2)
#
#      tkgrid(tklabel(frame2,text="Output format : "), sticky="n")
#      tkgrid(tkradiobutton(frame2, text="postscript", value=1, variable=formatvar), sticky="w")
#      tkgrid(tkradiobutton(frame2, text="pdf", value=2, variable=formatvar), sticky="w")
#      tkgrid(tkradiobutton(frame2, text="pictex", value=3, variable=formatvar), sticky="w")
#      tkgrid(tkradiobutton(frame2, text="xfig", value=4, variable=formatvar), sticky="w")
#      tkgrid(tkradiobutton(frame2, text="png", value=5, variable=formatvar), sticky="w")
#      tkgrid(tkradiobutton(frame2, text="jpeg", value=6, variable=formatvar), sticky="w")
#      tkgrid(frame2, rowspan=2, sticky="n")
#      
#      tkgrid(tklabel(frame3,text="Output size : "))
#      width.entry <- tkentry(frame3, textvariable=widthvar, width=10)
#      height.entry <- tkentry(frame3, textvariable=heightvar, width=10)
#      tkgrid(tklabel(frame3,text="Width : "), width.entry)
#      tkgrid(tklabel(frame3,text="Height : "), height.entry)
#      tkgrid(frame3, column=1, row=1, sticky="n")
#
#      save.but <- tkbutton(frame1, text="Save", command=function() savefic(formatvar, widthvar, heightvar))
#      cancel.but <- tkbutton(frame1, text="Dismiss", command=function() tkdestroy(tf))
#      tkgrid(save.but, cancel.but)
#      tkgrid(frame1, column=1, row=2, sticky="n")
#      
#      tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
#      tkbind(tf, "<KeyPress-Return>", function() savefic(formatvar, widthvar, heightvar))
#      tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
#      tkwait.variable(done)
#      if(tclvalue(done) == "2") return(0)
#      tkdestroy(tf)
#    }
#
#
#  
#  menurightbutton <- function(x,y){
#    menuright <-tkmenu(plot1, tearoff=FALSE)
#    ## Copy seulement pour windows (a verifier)
#    if(.Platform$OS.type=="windows") {
#      tkadd(menuright,"command",label="Copy",command=copytoclip)
#      tkadd(menuright,"separator")
#    }
#    
#    tkadd(menuright,"command",label="Save",command=outgraph.explore)
#    tcl("tk_popup",menuright,x,y)
#    
#  }
#
#
#  ## -----------------------------------------------------------------------
#  ##    Configuration des fonctions associees a plot1 en fonction du mode
#  ## -----------------------------------------------------------------------
#
#  modezoom <- function(){
#    tkbind(plot1, "<Button-1>",getxy)
#    tkbind(plot1, "<B1-Motion>",draw.rect)
#    tkbind(plot1, "<ButtonRelease-1>",zoom.rect)
#    tkconfigure(plot1,cursor="crosshair")
#  }
#
#  modelabel <- function(){
#    tkbind(plot1, "<Button-1>",addonelabel)
#    tkbind(plot1, "<B1-Motion>",nothing)
#    tkbind(plot1, "<ButtonRelease-1>",nothing)
#    tkconfigure(plot1,cursor="top_tee")
#  }
#  
#  modeidentify <- function(){
#    tkbind(plot1, "<Button-1>",identifyclosestpoint)
#    tkbind(plot1, "<B1-Motion>",nothing)
#    tkbind(plot1, "<ButtonRelease-1>",nothing)
#    tkconfigure(plot1,cursor="question_arrow")
#  }
#  
#  modespan <- function(){
#    tkbind(plot1, "<Button-1>",plot.s.span)
#    tkbind(plot1, "<B1-Motion>",nothing)
#    tkbind(plot1, "<ButtonRelease-1>",nothing)
#    tkconfigure(plot1,cursor="hand1")
#    
#  }
#  
#  ## --------------------------------------------------
#  ##    Configuration des fonctions associees a plot2 
#  ## --------------------------------------------------
#
#  
#  plot.s.span <- function (x,y) {
#    newxy <-change.systeme.coordonnes(x,y)
#    plot.s.change.origin(newxy[1],newxy[2])
#  }
#
#  plot.s.span.small <- function (x,y) {
#    newxy <-change.systeme.coordonnes(x,y,plot=plot2,parplosize=parPlotSize.small,usrcor=usrCoords.small)
#    plot.s.change.origin(newxy[1],newxy[2])
#  }
#
#  ## -----------------------------------------------------------------------
#  ##    Fonction de  changement de mode 
#  ## -----------------------------------------------------------------------
#
#  changermode <- function(){
#    rbVal <- as.character(tclvalue(rbValue))
#    if(rbVal!=rbValold){
#      rbValold <<- rbVal
#      if(rbVal=="zoom") {
#        tkconfigure(b1.lab,state="disabled")
#        tkconfigure(b2.lab,state="disabled")
#        tkconfigure(b3.lab,state="disabled")
#        tkconfigure(b4.lab,state="disabled")
#        tkconfigure(titrelab,state="disabled")
#        modezoom()
#      }
#      else if(rbVal=="span") {
#        tkconfigure(b1.lab,state="disabled")
#        tkconfigure(b2.lab,state="disabled")
#        tkconfigure(b3.lab,state="disabled")
#        tkconfigure(b4.lab,state="disabled")
#        tkconfigure(titrelab,state="disabled")
#        modespan()
#      }
#      else if(rbVal=="identify") {
#        tkconfigure(b1.lab,state="disabled")
#        tkconfigure(b2.lab,state="disabled")
#        tkconfigure(b3.lab,state="disabled")
#        tkconfigure(b4.lab,state="disabled")
#        tkconfigure(titrelab,state="disabled")
#        modeidentify()
#      }
#      else if(rbVal=="label") {
#        tkconfigure(b1.lab,state="normal")
#        tkconfigure(b2.lab,state="normal")
#        tkconfigure(b3.lab,state="normal")
#        tkconfigure(b4.lab,state="normal")
#        tkconfigure(titrelab,state="normal")
#        modelabel()
#      }
#    }
#  }
#  
#  ## -----------------------------------------------------------------------
#  ##    tctk : interface globale
#  ## -----------------------------------------------------------------------
#
#  tt <- tktoplevel()
#  tkwm.title(tt,"Explore.s")
#  zoomfactor  <-  tclVar(100)
#  Show.Value <- tclVar("0")
#  
#  ## 3 frame 2 en haut, 1 en bas
#  global.frm <- tkframe(tt,borderwidth=2)
#  top.frm <- tkframe(global.frm)
# 
#  top.left.frm <- tkframe(top.frm)
#  top.right.frm <- tkframe(top.frm)
#  ## Top Left frame : plot
#  frame1 <- tkframe(top.left.frm, relief="groove", borderwidth=2)
#  usrCoords.ori <- NULL
#  plot1 <- tkrplot(frame1, fun=plot.s,hscale=scale.graph,vscale=scale.graph )
#  tkpack(plot1,side="top")
#  tkbind(plot1, "<Button-3>",menurightbutton)
#  densx <- density(dftoplot[,1],from=usrCoords[1],to=usrCoords[2])
#  densy <- density(dftoplot[,2],from=usrCoords[3],to=usrCoords[4])
#
#  ## Top Right frame : small plot + text
#  frame2b <- tkframe(top.right.frm, relief="groove", borderwidth=2)
#  plot2 <- tkrplot(frame2b, fun=plot.s.small,hscale=.4,vscale=.4 )
#  plot3 <- tkrplot(frame2b, fun=plot.s.densx,hscale=.4,vscale=.1 )
#  plot4 <- tkrplot(frame2b, fun=plot.s.densy,hscale=.1,vscale=.4 )
#  
#  tkpack(plot3,side="top",anchor="w")
#  tkpack(plot4,side="right",anchor="center")
#  tkpack(plot2,side="left",anchor="center")
#  tkconfigure(plot2,cursor="hand1")
#  tkbind(plot2, "<Button-1>",plot.s.span.small)
#
#  frame2 <- tkframe(top.right.frm, relief="groove", borderwidth=2)
#
#  ##  label des points
#  ## fenetre texte
#  oldwidth <- options()$width
#  options(width = 10000)
#  conn <- textConnection("zz", "w", local=TRUE)
#  sink(conn)
#  print(dftoplot)
#  sink()
#  close(conn)
#  options(width = oldwidth)
#
#  namesrowchar <- row.names(dftoplot)
#  wnamesrowchar <- max(nchar(namesrowchar))
#  namescolchar <- substring(zz[1], 2)
#  datachar <- substring(zz, 2 + wnamesrowchar)
##  rm(zz, envir = .GlobalEnv)
#  datachar <- datachar[-1]
#  wdatachar <- max(nchar(datachar))
#  wnamesrowchar <- max(nchar(namesrowchar))
#  hdatachar <- min(nrow(dftoplot),15)
#  scrv <- tkscrollbar(frame2, repeatinterval=5,command=function(...) {tkyview(txtdata,...);tkyview(txtlabrow,...);tkyview(txtlabcol,...);})
#  txtlabrow <- tktext(frame2,bg="grey",font="courier 12",yscrollcommand=function(...)tkset(scrv,...),width = wnamesrowchar,height=hdatachar)
#  txtlabcol <- tktext(frame2,bg="blue",font="courier 12",yscrollcommand=function(...)tkset(scrv,...),width = wnamesrowchar+wdatachar,height=1)
#  txtdata <- tktext(frame2,bg="white",font="courier 12",yscrollcommand=function(...)tkset(scrv,...),width = max(wdatachar,21),height=hdatachar)
#  tkinsert(txtlabcol,"end", namescolchar)
#  tkinsert(txtlabrow,"end", paste(namesrowchar, collapse = "\n"))
#  tkinsert(txtdata,"end", paste(datachar, collapse = "\n"))
#  tkconfigure(txtdata, state="disabled")
#  tkconfigure(txtlabrow, state="disabled")
#  tkconfigure(txtlabcol, state="disabled")
#  tkpack(tklabel(frame2,text="Information:"))
#  tkpack(scrv, side="right",fill="y")
#  tkpack(txtlabcol, side="top",fill="x",anchor="e")
#  tkpack(txtlabrow, side="left")
#  tkpack(txtdata,side="left")
#  tkpack(frame1, side="top",anchor="w")
#  tkpack(frame2b,side="top",anchor="center")
#
#  ## -----------------------------------------------------------------------
#  ##    tcltk : interface pour le choix du Mode 
#  ## -----------------------------------------------------------------------
#  
#  frame6 <- tkframe(top.right.frm, relief="groove", borderwidth=2,background = "white")
#  rb1 <- tkradiobutton(frame6,text="Zoom",foreground = "darkgreen", background = "white")
#  rb2 <- tkradiobutton(frame6,text="Pan",foreground = "darkgreen", background = "white")
#  rb3 <- tkradiobutton(frame6,text="Labels",foreground = "darkgreen", background = "white")
#  rb4 <- tkradiobutton(frame6,text="Identify",foreground = "darkgreen", background = "white")
#  rbValue <- tclVar("zoom")
#  rbValold <- as.character(tclvalue(rbValue))
#  tkconfigure(rb1,variable=rbValue,value="zoom",command=changermode)
#  tkconfigure(rb2,variable=rbValue,value="span",command=changermode)
#  if(type=="s.value"){
#    tkconfigure(rb3,variable=rbValue,value="label",command=changermode,state="disabled")
#  }
#  else{
#    tkconfigure(rb3,variable=rbValue,value="label",command=changermode)
#  }
#  tkconfigure(rb4,variable=rbValue,value="identify",command=changermode)
#  tkpack(tklabel(frame6,text="Mode options:",foreground = "darkgreen", background = "white"))
#  tkpack(rb1,rb2,rb3,rb4,side="left",fill="x")
#  cb.show <- tkcheckbutton(top.right.frm)
#  tkconfigure(cb.show,variable=Show.Value,text="Show commands in R console")
#  
#  tkpack(frame6,side="top",anchor="center")
#  tkpack(cb.show,anchor="center")
#  modezoom()
#  ## Find menu
#  frame74 <- tkframe(top.right.frm, relief="groove", borderwidth=2)
#  frame4 <- tkframe(frame74, relief="groove", borderwidth=2)
#  scr2 <- tkscrollbar(frame4, repeatinterval=5,
#                      command=function(...)tkyview(tl.find,...))
#  tl.find<-tklistbox(frame4,height=7,width=14,selectmode="single",yscrollcommand=function(...)tkset(scr2,...),background="white")
#  tkpack(tklabel(frame4,text="Find a point"))
#  point.found <- tclVar("")
#  entry.found <- tkentry(frame4,width="14",textvariable=point.found,background="white")
#
#  tkbind(entry.found, "<KeyRelease>",find.key)
#  tkbind(entry.found, "<Return>",find.return)
#  
#  tkpack(entry.found, side="top", fill="x", expand=TRUE)
#  tkpack(tl.find, side="left", fill="both", expand=TRUE)
#  tkpack(scr2, side="right", fill="y")
#  name.sites <- row.names(dftoplot)
#  if(exists("label.ori"))
#    name.sites <- label.ori
#    name.sites.new <- name.sites  
#  for (i in 1:length(name.sites))
#    {
#      tkinsert(tl.find,"end",name.sites[i])
#    }
# 
#  tkbind(tl.find, "<ButtonRelease-1>",recenterplot)
#  ## -----------------------------------------------------------------------
#  ##    tcltk : Menu Label 
#  ## -----------------------------------------------------------------------
#  
#  frame7 <- tkframe(frame74, relief="groove", borderwidth=2)
#  b1.lab <- tkbutton(frame7,text=" Add all ")
#  b2.lab <- tkbutton(frame7,text=" Remove all ")
#  b3.lab <- tkbutton(frame7,text=" Select... ")
#  b4.lab <- tkbutton(frame7,text=" By values... ")
#  tkconfigure(b1.lab,command=addalllabel,state="disabled")
#  tkconfigure(b2.lab,command=removelabel,state="disabled")
#  tkconfigure(b3.lab,command=addselectlabel,state="disabled")
#  tkconfigure(b4.lab,command=labelbyvalues,state="disabled")
#  titrelab<-tklabel(frame7,text="Label \n options: \n",state="disabled")
#  tkpack(titrelab)
#  tkpack(b1.lab,b2.lab,anchor="center",fill="x")
#  tkpack(b3.lab,b4.lab,anchor="center",fill="x")
#  tkpack(frame4,anchor="e",side="left",fill="x")
#  tkpack(frame7,anchor="w",side="right")
#  tkpack(frame74,fill="x")
#
#  tkpack(frame2,fill="both")
# 
#
#  ## -----------------------------------------------------------------------
#  ##    tcltk : Boutons Quit, Refresh... 
#  ## -----------------------------------------------------------------------
#  low.frm <- tkframe(top.left.frm,relief="groove", borderwidth=2)
#  quit.fn <- function(){
#    list.obj <- ls(pos=sys.frame())
#    obj.glob <- c("add.rectzoom.plot1","usrCoords.ori","parPlotSize","usrCoords","parPlotSize.small","usrCoords.small","appel.list","add.labels.plot1","zoomfactor","xy1.rect","labelled.points","wh.include.row.tt2","xinf","yinf","xsup","ysup","rbValue.ttlabval","rbValold","name.sites","name.sites.new","point.found","rbValue.tt2","tl.select.tt2")
#    to.remove <- obj.glob[pmatch(list.obj,obj.glob)[!is.na(pmatch(list.obj,obj.glob))]]
##    rm(list=eval(to.remove),pos=sys.frame())
#    tkdestroy(tt)
#  }
#
#  frame5 <- tkframe(low.frm)
#  button.quit <- tkbutton(frame5, text="   Dismiss   ",command=function() quit.fn(),foreground = "darkgreen", background = "white")
#  button.refresh <- tkbutton(frame5, text="   Redraw   ",command=plot.s.refresh,foreground = "darkgreen", background = "white")
#  button.save <- tkbutton(frame5, text="   Save   ",command=outgraph.explore,foreground = "darkgreen", background = "white")
#  tkgrid(button.quit,button.save,button.refresh)
# 
#
#  ## -----------------------------------------------------------------------
#  ## tcltk : zoom scroll
#  ## -----------------------------------------------------------------------
#  ## Zoom
#  maxzoom <- 1000
#  frame3 <- tkframe(low.frm)
#  tkpack(tklabel(frame3,text="Zoom factor "),side="left",anchor="n")
#  scaleview <- tkscale(frame3, command=redo.plot.s.zoomfac, from=100, to=maxzoom, showvalue=FALSE, variable=zoomfactor, resolution=10, orient="horizontal", length=300, tick=100)
#  tkpack(scaleview, anchor="center")
#  tkpack(frame5,side="left",fill="both")
#  tkpack(frame3,side="right")
#
#  ## -----------------------------------------------------------------------
#  ## tcltk : mise en page globale
#  ## -----------------------------------------------------------------------
#
#
## tkpack(low.frm, fill="both", side="bottom",expand=TRUE)
#  tkpack(low.frm, anchor="s",fill="x")
#  tkpack(top.left.frm, side="left", anchor="n", fill="both", expand=TRUE)
#  tkpack(top.right.frm, side="right", anchor="n", fill="both", expand=TRUE)
#
#  tkpack(top.frm, fill="both", expand=TRUE)
#  tkpack(global.frm, fill="both", expand=TRUE)
#  done <- tclVar(0)
#  tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
##  tkwait.variable(done)
#
#  if(tclvalue(done) == "2") quit.fn()

}
