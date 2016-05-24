viewR<-function(data=NULL,xyz=NULL,ret=FALSE){
  ###############################
  ###### DATA CHECK #############
  ###############################
  if(is.null(data)){
    file<-selectR()
    if(length(file)==0){stop("You must select a file")}
    data<-readNii(file)
  }
  if(is.character(data)){
    file<-data
    if(length(file)==0){stop("You must select a file")}
    data<-readNii(file)
  }
  if(!is.array(data)){
    stop("This is not an array with 3 or more dimensions")
  }
  orig<-data
  ##############################
  ##### DESIRED HEIGHT #########
  ##############################
  screenHeight<-as.numeric(strsplit(tclvalue(tkwm.maxsize("."))," ")[[1]])[2]
  desHeight<-round(screenHeight/5*3)
  ##############################
  ##### ASPECT RATIO ###########
  ##############################
  d<-dim(data)
  d1<-d[1]
  d2<-d[2]
  d3<-d[3]
  xi<-xyz[1]
  yi<-xyz[2]
  zi<-xyz[3]
  xx <- ifelse(is.null(xi)||xi>d1||xi<1||!is.numeric(xi), round(d1/2), xi)
  yy <- ifelse(is.null(yi)||yi>d2||yi<1||!is.numeric(xi), round(d2/2), yi)
  zz <- ifelse(is.null(zi)||zi>d3||zi<1||!is.numeric(xi), round(d3/2), zi)
  
  if(length(d)>3){
    d4<-d[4]
  }else{
    dim(data)<-c(dim(data),1)
    d4<-1
  }
  
  xCoords<-1:d1
  yCoords<-1:d2
  zCoords<-1:d3
  tCoords<-1:d4
  at<-attr(data,"pixdim")
  if(is.null(at)){pixdim<-c(1,1,1)}else{pixdim<-at}
  if(length(pixdim)==8){pixdim<-pixdim[2:4]}
  data<-zeroNa(input = data)
  w1<-d1*pixdim[1]
  w2<-d2*pixdim[2]
  w3<-d3*pixdim[3]
  #########################
  ### IMAGE VARIABLES #####
  #########################
  r<-round(range(data),digits = 2)
  time<-tclVar("1")
  xl<-tclVar(as.character(xx))
  yl<-tclVar(as.character(yy))
  zl<-tclVar(as.character(zz))
  X<-tclVar(tclvalue(xl))
  Y<-tclVar(tclvalue(yl))
  Z<-tclVar(tclvalue(zl))
  parPlotSize1<-c()
  usrCoords1<-c()
  low<-tclVar(as.character(r[1]))
  high<-tclVar(as.character(r[2]))
  intens<-tclVar(as.character(round(data[xx,yy,zz,1],digits = 2)))
  crossHairsOn<-TRUE
  tmovie<-tclVar(FALSE)
  after_ID <-""
  someFlag<-TRUE
  #########################
  ####### FRAMES ##########
  #########################
  tclServiceMode(FALSE)
  base<-tktoplevel()
  tktitle(base)<-"Display"
  master<-tkframe(parent=base)
  f1<-tkframe(parent = master)
  #########################
  ##### FUNCTIONS #########
  #########################
  plotf<-function(){
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    mat<-matrix(1:4,nrow = 2,ncol = 2,byrow = TRUE)
    nf<-layout(mat, widths = c(w1,w2),heights = c(w3,w2), respect = FALSE)
    xl<-as.numeric(tclvalue(xl))
    yl<-as.numeric(tclvalue(yl))
    zl<-as.numeric(tclvalue(zl))
    x<-round(xl)
    y<-round(yl)
    z<-round(zl)
    tclvalue(X)<<-x
    tclvalue(Y)<<-y
    tclvalue(Z)<<-z
    t<-as.numeric(tclvalue(time))
    if(is.na(t)){t<-1}
    tclvalue(time)<<-t
    lim<-c(as.numeric(tclvalue(low)),as.numeric(tclvalue(high)))
    if(any(is.na(lim))){lim<-r;tclvalue(low)<<-r[1];tclvalue(high)<<-r[2]}
    tclvalue(intens)<<-round(data[x,y,z,t],digits = 2)
    im1<-data[,y,,t]
    im1[im1<lim[1]]<-lim[1]
    im1[im1>lim[2]]<-lim[2]
    image(1:d1,1:d3,z=im1,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = zl,v = xl,col="green")
    im2<-data[x,,,t]
    im2[im2<lim[1]]<-lim[1]
    im2[im2>lim[2]]<-lim[2]
    image(1:d2,1:d3,z=im2,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = zl,v = yl,col="green")
    im3<-data[,,z,t]
    im3[im3<lim[1]]<-lim[1]
    im3[im3>lim[2]]<-lim[2]
    image(1:d1,1:d2,z=im3,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = yl,v = xl,col="green")
    if(d4>1){
    tseries<-data[x,y,z,]
    plot(1:d4,tseries,type="l",axes=FALSE,col="grey",xlim=c(1,d4))
    points(x = t,y = tseries[t],col="green")
    }
    parPlotSize1 <<- par("plt")
    usrCoords1   <<- par("usr")
  }
  OnLeftClick1 <- function(x,y){
    xClick<-x
    yClick<-y
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img1)))-2
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img1)))-2
    xBorder<-w1/(w1+w2)*(width)
    yBorder<-w3/(w3+w2)*(height)
    xClick <- as.numeric(xClick)
    yClick <- as.numeric(yClick)
    if(xClick<=xBorder && xClick>0 && yClick>0 && yClick<=height ){
      xl<-xClick/xBorder*d1+.5
      if(xl>d1){xl<-d1}
      if(xl<1){xl<-1}
      tclvalue(xl)<<-xl
      tclvalue(X)<<-round(xl)
    }
    if(yClick<=yBorder && yClick>0&& xClick>0 &&xClick<width){
      zl<-d3-yClick/yBorder*d3+.5
      if(zl>d3){zl<-d3}
      if(zl<1){zl<-1}
      tclvalue(zl)<<-zl
      tclvalue(Z)<<-round(zl)
    }
    if(xClick>xBorder && xClick<width && yClick<=yBorder && yClick>0){
      yl<-(xClick-xBorder)/(width-xBorder)*d2+.5
      if(yl>d2){yl<-d2}
      if(yl<1){yl<-1}
      tclvalue(yl)<<-yl
      tclvalue(Y)<<-round(yl)
    }
    if(yClick>yBorder && xClick<xBorder && yClick<=height && xClick>0){
      yl<-d2-(yClick-yBorder)/(height-yBorder)*d2+.5
      if(yl>d2){yl<-d2}
      if(yl<1){yl<-1}
      tclvalue(yl)<<-yl
      tclvalue(Y)<<-round(yl)
    }
    if(yClick>yBorder&&xClick>xBorder &&d4>1){
      first<-1.04-.04*d4
      last<- -.04+1.04*d4
      s<-round(first):round(last)
      rangeUsr<-last-first
      rangePix<-width-xBorder
      click<-xClick-xBorder
      tind<-round(click/rangePix*length(s))
      t<-s[tind]
      if(t<1){t<-1}
      if(t>d4){t<-d4}
      tclvalue(time)<<-t
      }
    tkrreplot(img1)
  }
  changeTemp<-function(){
    if(is.na(strtoi(tclvalue(time)))){
      tclvalue(time)<<-as.character(1)
    }
    t<-strtoi(tclvalue(time))
    if(t>d4){tclvalue(time)<<-as.character(d4)}
    if(t<1){tclvalue(time)<<-as.character(1)}
    tkrreplot(img1)
  }
  tkspinbox <- function(parent, ...) {tkwidget(parent, "tk::spinbox", ...)}
  onSpin<-function(){
    if(is.na(strtoi(tclvalue(X)))||strtoi(tclvalue(X))==0){
      tclvalue(xl)<<-as.character(xx)
      tclvalue(X)<<-as.character(xx)}else{
        tclvalue(xl)<<-tclvalue(X)
      }
    
    if(is.na(strtoi(tclvalue(Y)))||strtoi(tclvalue(Y))==0){
      tclvalue(yl)<<-as.character(yy)
      tclvalue(Y)<<-as.character(yy)}else{
        tclvalue(yl)<<-tclvalue(Y)
      }
    
    if(is.na(strtoi(tclvalue(Z)))||strtoi(tclvalue(Z))==0){
      tclvalue(zl)<<-as.character(zz)
      tclvalue(Z)<<-as.character(zz)}else{
        tclvalue(zl)<<-tclvalue(Z)
      }
    if(is.na(strtoi(tclvalue(time)))||strtoi(tclvalue(time))==0){
      tclvalue(time)<<-as.character(1)}else{
        tclvalue(time)<<-tclvalue(time)
      }
    x<-round(as.numeric(tclvalue(xl)))
    y<-round(as.numeric(tclvalue(yl)))
    z<-round(as.numeric(tclvalue(zl)))
    if(x<1){tclvalue(xl)<<-1}
    if(y<1){tclvalue(yl)<<-1}
    if(z<1){tclvalue(zl)<<-1}
    if(x>d1){tclvalue(xl)<<-d1}
    if(y>d2){tclvalue(yl)<<-d2}
    if(z>d3){tclvalue(zl)<<-d3}
    tkrreplot(img1)
  }
  maxminChange<-function(){
    if(suppressWarnings(is.na(as.numeric(tclvalue(high))))||tclvalue(high)==""){
      tclvalue(high)<<-r[2]
    }
    if(suppressWarnings(is.na(as.numeric(tclvalue(low))))||tclvalue(low)==""){
      tclvalue(low)<<-r[1]
    }
    max<-as.numeric(tclvalue(high))
    min<-as.numeric(tclvalue(low))
    if(max<min){tclvalue(high)<<-r[2];tclvalue(low)<<-r[1]}
    tkrreplot(img1)  
  }
  crossHairs<-function(){
    crossHairsOn<<-!crossHairsOn
    tkrreplot(img1)
  }
  repeat_call<-function(ms = 200 , f) {
    after_ID <<- tcl( "after" , ms,function(){ 
      if(someFlag){
        f()
        after_ID<<-repeat_call(ms,f)
      }else{
        tcl("after" , "cancel" , after_ID)
      }
    })
  }
  movieT<-function(){
    if(tclvalue(tmovie)==1){
      someFlag<<-TRUE
      repeat_call(1,function() {tkinvoke(spin,"buttonup")})
    }else{someFlag<<-FALSE}
  }
  ##########################
  ###### THE PLOTS #########
  ##########################
  testImg<-tkrplot(parent = f1,fun =plotf)
  testHeight<-as.numeric(tkwinfo("reqheight",testImg))
  testWidth<-as.numeric(tkwinfo("reqwidth",testImg))
  scaleFactor<-(desHeight)/(testHeight)
  asp<-(w1+w2)/(w2+w3)
  if(asp>1){
  img1<-tkrplot(parent = f1,fun = plotf,hscale=scaleFactor,vscale=scaleFactor/asp)
  }else{
    img1<-tkrplot(parent = f1,fun = plotf,hscale=scaleFactor*asp,vscale=scaleFactor)
  }
  f5<-tkframe(parent=master,
              width=as.numeric(tkwinfo("reqwidth",img1)),
              height=100,
              borderwidth=0,
              relief="groove")
  ##########################
  ##### Widgets ############
  ##########################
  spin<-tkspinbox(f5,textvariable=time,from=1,to=d4,command=changeTemp,increment=1,repeatdelay=200,width=5,wrap=TRUE)
  coX<-tkspinbox(f5,textvariable=X,from=1,to=d1,command=function()onSpin(),increment=1,repeatdelay=10,width=5)
  coY<-tkspinbox(f5,textvariable=Y,from=1,to=d2,command=function()onSpin(),increment=1,repeatdelay=10,width=5)
  coZ<-tkspinbox(f5,textvariable=Z,from=1,to=d3,command=function()onSpin(),increment=1,repeatdelay=10,width=5)
  MAX<-tkentry(f5,textvariable=high,width=7)
  MIN<-tkentry(f5,textvariable=low,width=7)
  intensity<-tkentry(f5,textvariable=intens,state="readonly",readonlybackground="white",width=7)
  coXLab<-tklabel(parent = f5,text="X")
  coYLab<-tklabel(parent = f5,text="Y")
  coZLab<-tklabel(parent = f5,text="Z")
  coTLab<-tklabel(parent = f5,text="Time")
  maxLab<-tklabel(parent = f5,text="Max")
  minLab<-tklabel(parent = f5,text="Min")
  intensLab<-tklabel(parent=f5,text="Intensity")
  tmovieBut<-tkcheckbutton(f5,variable=tmovie, command=movieT,text="Movie")
  ##########################
  ###### GEOMETRY ##########
  ##########################
  tkgrid(master,columnspan=2)
  tkgrid(f1)
  tkgrid(f5)
  tkgrid(img1)
  tkgrid(f5,columnspan=2,rowspan=1)
  tkplace(coXLab,relx=0,rely=1/5,y=-10)
  tkplace(coYLab,relx=0,rely=2/5,y=-10)
  tkplace(coZLab,relx=0,rely=3/5,y=-10)
  tkplace(coTLab,relx=0,rely=4/5,y=-10)
  tkplace(coX,'in'=coXLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coY,'in'=coYLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coZ,'in'=coZLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(spin,'in'=coTLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  if(d4>1){tkplace(tmovieBut,"in"=spin,rely=.5,anchor="w",x=as.numeric(tkwinfo("reqwidth",spin)),height=100/4)}
  wiBox<-as.numeric(tkwinfo("reqwidth",intensity))+10
  tkplace(intensity,relx=1,x=-wiBox,rely=0/4,y=10)
  tkplace(MAX,relx=1,x=-wiBox,rely=1/4,y=10)
  tkplace(MIN,relx=1,x=-wiBox,rely=2/4,y=10)
  tkplace(intensLab,'in'=intensity,x=-1,rely=.5,anchor="e",height=100/4)
  tkplace(maxLab,'in'=MAX,x=-1,rely=.5,anchor="e",height=100/4)
  tkplace(minLab,'in'=MIN,x=-1,rely=.5,anchor="e",height=100/4)
  ###########################
  ##### BINDINGS ############
  ###########################
  tkbind(MAX, "<Return>",function()maxminChange())
  tkbind(MIN, "<Return>",function()maxminChange())
  tkbind(coX, "<Return>",function()onSpin())
  tkbind(coY, "<Return>",function()onSpin())
  tkbind(coZ, "<Return>",function()onSpin())
  tkbind(spin, "<Return>",function()changeTemp())
  tkbind(img1, "<Button-3>",crossHairs)
  tkbind(img1, "<Button-1>",OnLeftClick1)
  tkbind(img1, "<B1-Motion>",OnLeftClick1)
  tkconfigure(img1,cursor="crosshair")
  tclServiceMode(TRUE)
  tkwm.resizable(base,FALSE,FALSE)
  tkfocus(base)
  if(ret){return(orig)}
}

