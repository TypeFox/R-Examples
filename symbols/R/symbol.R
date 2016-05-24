symbol <-
function(df,type="star",colin=NULL,colout=NULL,colin2=NA,colout2=1,ssize=NULL,labels=0,labelsize=0.6,scheme=1,sortby=0,descending=FALSE,coorx=NULL,coory=NULL,lty=1,
                 main=NULL,sub=NULL,xlab=NULL,ylab=NULL,add=FALSE,xlim=NULL,ylim=NULL,
                 facew=0.5,faceh=0.5,eyes=0.5,eyed=0.5,mouthw=0.5,mouthc=0.5,brows=0.5,browp=0.5,nosel=0.5,nosew=0.5,ears=0.5,pupils=0.5,
                 body=0.5,limb1=0.5,limb2=0.5,limb3=0.5,limb4=0.5,defcol=NA){
 n<-dim(df)[1]
 m<-dim(df)[2]
 if (is.null(coorx) & !is.null(coory)) coorx<-1:n
 if (!is.null(coorx) & is.null(coory)) coory<-1:n
 # Plot range
 if (is.null(coorx) & is.null(coory)){
  minx<-0
  maxx<-10
  miny<-0
  maxy<-10
 } else
 {
  sortby<-0
  minx<-min(coorx)-(max(coorx)-min(coorx))*0.05
  maxx<-max(coorx)+(max(coorx)-min(coorx))*0.05
  miny<-min(coory)-(max(coory)-min(coory))*0.05
  maxy<-max(coory)+(max(coory)-min(coory))*0.05
 }
 if (!add) plot.new() else{
  if (!is.null(xlim) & !is.null(ylim)){
   minx<-xlim[1]
   maxx<-xlim[2]
   miny<-ylim[1]
   maxy<-ylim[2]
  }
 }
 lengthx<-maxx-minx
 lengthy<-maxy-miny
 plot.window(c(minx,maxx),c(miny,maxy))
 if (!is.null(coorx) & !is.null(coory) & !add){
  axis(1)
  axis(2)
  box()
 }
 title(main=main,sub=sub,xlab=xlab,ylab=ylab)
 # Undefined faces
 if (type=="face"){
  if (is.na(colin2)) colin2<-"white"
  if (facew==0.5 & faceh==0.5 & eyes==0.5 & eyed==0.5 & mouthw==0.5 & mouthc==0.5 & brows==0.5 & browp==0.5 & nosel==0.5 & nosew==0.5 & ears==0.5 & pupils==0.5){
   for (i in 1:m){
    if (is.numeric(df[,i])){
     if (facew==0.5) facew<-i else if (faceh==0.5) faceh<-i else if (eyes==0.5) eyes<-i else if (eyed==0.5) eyed<-i else if (mouthw==0.5) mouthw<-i else if (mouthc==0.5) mouthc<-i else if (brows==0.5) brows<-i else if (browp==0.5) browp<-i else if (nosel==0.5) nosel<-i else if (nosew==0.5) nosew<-i else if (ears==0.5) ears<-i else if (pupils==0.5) pupils<-i
    }
   }
  }
 }
 # Undefined sticks
 if (type=="stick"){
  if (body==0.5 & limb1==0.5 & limb2==0.5 & limb3==0.5 & limb4==0.5 ){
   for (i in 1:m){
    if (is.numeric(df[,i])){
     if (body==0.5) body<-i else if (limb1==0.5) limb1<-i else if (limb2==0.5) limb2<-i else if (limb3==0.5) limb3<-i else if (limb4==0.5) limb4<-i
    }
   }
  }
 }
 numattr<-nnumattr(df)
 # Color segments of color icons
 if (type=="icon"){
  colorseg<-rep(defcol,times=8)
  if (numattr==1){
   for (i in 1:m) if (is.numeric(df[,i])){
    colorseg<-rep(i,times=8)
    break
   }
  } else
  if (numattr==2){
   a<-0
   for (i in 1:m) if (is.numeric(df[,i])){
    colorseg[(1+a*4):(4+a*4)]<-i
    a<-a+1
    if (a==2) break
   }
  } else
  if (numattr<=4){
   a<-0
   for (i in 1:m) if (is.numeric(df[,i])){
    colorseg[(1+a*2):(2+a*2)]<-i
    a<-a+1
    if (a==4) break
   }
  } else
  {
   a<-0
   for (i in 1:m){
    if (is.numeric(df[,i])){
     colorseg[a+1]<-i
     a<-a+1
    }
    if (a==8) break
   }
  }
 }
 # Symbol area, sorting, labels
 squares<-ceiling(sqrt(n))
 if (sortby>0 & sortby<=m) if (!descending) {df<-df[order(df[,sortby]),]} else {df<-df[order(-df[,sortby]),]}
 if (labels>0 & labels<=m) labels<-df[,labels] else labels=NULL
 border=colout2
 col=colin2
 # Cathegorical colors
 colrs<-1:n
 colrs2<-1:n
 if (!is.null(colout)){
  if (!is.numeric(df[,colout])) 
   for (i in 1:n){
    for (a in 1:nlevels(df[,colout])) if (levels(df[,colout])[a]==df[i,colout]) colrs[i]<-rainbow(nlevels(df[,colout]))[a]
   }
 }
 if (!is.null(colin)){
  if (!is.numeric(df[,colin])) 
   for (i in 1:n){
    for (a in 1:nlevels(df[,colin])) if (levels(df[,colin])[a]==df[i,colin]) colrs2[i]<-rainbow(nlevels(df[,colin]))[a]
   }
 }
 # Coordinates and size of symbols
 if (!is.null(coorx) & !is.null(coory)){
  if (is.null(ssize)){
   if (type=="stick") ssize<-0.02 else if (type=="icon") ssize=0.01 else ssize<-0.04
  }
 }else
 {
  coorx<-1:n
  coory<-1:n
  if (scheme==3){
   coorx=NULL
   coory=NULL
   for (i in 0:(squares%/%2)){
    coorx<-c(coorx,-i:i)
    coorx<-c(coorx,seq(i+1,i+1,length=(i+1)*2))
    coorx<-c(coorx,i:-i)
    coorx<-c(coorx,seq(-i-1,-i-1,length=(i+1)*2))
    coory<-c(coory,seq(i,i,length=i*2+1))
    coory<-c(coory,i:(-i-1))
    coory<-c(coory,seq(-i-1,-i-1,length=i*2+1))
    coory<-c(coory,(-i-1):i)
   }
   if (squares%/%2==squares/2) shiftxy<-0.5 else shiftxy<-0
   coorx<-coorx*lengthx/squares+lengthx/2-lengthx/squares*shiftxy+minx
   coory<-coory*lengthy/squares+lengthy/2+lengthy/squares*shiftxy+miny
  } else
  {
   for (i in 1:n){
     if (scheme==2 & ((i-1)%/%squares)/2!=((i-1)%/%squares)%/%2) coorx[i]<-maxx-(((i-1)%%squares)*(lengthx/squares)+lengthx/squares/2)
      else coorx[i]<-((i-1)%%squares)*(lengthx/squares)+lengthx/squares/2+minx
     coory[i]<-maxy-(((i-1)%/%squares)*(lengthy/squares)+lengthy/squares/2)
   }
  }
  if (is.null(ssize)) ssize<-1/(squares*3)
 }
 # Normalization
 for (i in 1:m) df[,i]<-normalize(df[,i])
 polx<-1:numattr
 poly<-1:numattr
 angle<-2*pi/numattr
 for (i in 1:n){
  numberofnum<-numattr
  polx<-1:numattr
  poly<-1:numattr
  # Colors
  if (!is.null(colout)){
   if (is.numeric(df[,colout]))
    border<-hsv(0,1,df[i,colout]) else
    border<-colrs[i]
  }
  if (!is.null(colin)){
   if (is.numeric(df[,colin]))
    col<-hsv(0,0,df[i,colin]) else
    col<-colrs2[i]
  }
  if (type=="icon"){
   # Color icons
   rect(coorx[i]-ssize*lengthx,coory[i]-ssize*lengthy,coorx[i]+ssize*lengthx,coory[i]+ssize*lengthy,col=defcol,border=NA)
   polygon(c(coorx[i]-ssize*lengthx,coorx[i],coorx[i]),c(coory[i]+ssize*lengthy,coory[i]+ssize*lengthy,coory[i]),col=rainbow(133)[1+100*df[i,colorseg[1]]],border=NA)
   polygon(c(coorx[i]+ssize*lengthx,coorx[i],coorx[i]),c(coory[i]+ssize*lengthy,coory[i]+ssize*lengthy,coory[i]),col=rainbow(133)[1+100*df[i,colorseg[2]]],border=NA)
   polygon(c(coorx[i]+ssize*lengthx,coorx[i]+ssize*lengthx,coorx[i]),c(coory[i]+ssize*lengthy,coory[i],coory[i]),col=rainbow(133)[1+100*df[i,colorseg[3]]],border=NA)
   polygon(c(coorx[i]+ssize*lengthx,coorx[i]+ssize*lengthx,coorx[i]),c(coory[i]-ssize*lengthy,coory[i],coory[i]),col=rainbow(133)[1+100*df[i,colorseg[4]]],border=NA)
   polygon(c(coorx[i]+ssize*lengthx,coorx[i],coorx[i]),c(coory[i]-ssize*lengthy,coory[i]-ssize*lengthy,coory[i]),col=rainbow(133)[1+100*df[i,colorseg[5]]],border=NA)
   polygon(c(coorx[i]-ssize*lengthx,coorx[i],coorx[i]),c(coory[i]-ssize*lengthy,coory[i]-ssize*lengthy,coory[i]),col=rainbow(133)[1+100*df[i,colorseg[6]]],border=NA)
   polygon(c(coorx[i]-ssize*lengthx,coorx[i]-ssize*lengthx,coorx[i]),c(coory[i]-ssize*lengthy,coory[i],coory[i]),col=rainbow(133)[1+100*df[i,colorseg[7]]],border=NA)
   polygon(c(coorx[i]-ssize*lengthx,coorx[i]-ssize*lengthx,coorx[i]),c(coory[i]+ssize*lengthy,coory[i],coory[i]),col=rainbow(133)[1+100*df[i,colorseg[8]]],border=NA)
  } else
  if (type=="stick"){
   # Stick figures
   if (body==0.5) {bd<-body} else {bd<-df[i,body]}
   if (limb1==0.5) {l1<-limb1} else {l1<-df[i,limb1]}
   if (limb2==0.5) {l2<-limb2} else {l2<-df[i,limb2]}
   if (limb3==0.5) {l3<-limb3} else {l3<-df[i,limb3]}
   if (limb4==0.5) {l4<-limb4} else {l4<-df[i,limb4]}
   bx<-ssize*lengthx*sin(pi/2*(bd-0.5))*0.5
   by<-ssize*lengthy*sin(pi/2-pi/2*(bd-0.5))*0.5
   lines(c(coorx[i]+bx,coorx[i]-bx),c(coory[i]+by,coory[i]-by),col=border)
   lines(c(coorx[i]+bx-ssize*lengthx*sin(pi/2*(1-l1))*0.5,coorx[i]+bx),c(coory[i]+by+ssize*lengthy*sin(pi/2-pi/2*(1-l1))*0.5,coory[i]+by),col=border)
   lines(c(coorx[i]+bx+ssize*lengthx*sin(pi/2*l2)*0.5,coorx[i]+bx),c(coory[i]+by+ssize*lengthy*sin(pi/2-pi/2*l2)*0.5,coory[i]+by),col=border)
   lines(c(coorx[i]-bx+ssize*lengthx*sin(pi/2*(1-l3))*0.5,coorx[i]-bx),c(coory[i]-by-ssize*lengthy*sin(pi/2-pi/2*(1-l3))*0.5,coory[i]-by),col=border)
   lines(c(coorx[i]-bx-ssize*lengthx*sin(pi/2*l4)*0.5,coorx[i]-bx),c(coory[i]-by-ssize*lengthy*sin(pi/2-pi/2*l4)*0.5,coory[i]-by),col=border)
  } else
  if (type=="face"){
   # Chernoff faces
   if (facew==0.5) {fw<-ssize*lengthx/2*(1+facew)} else {fw<-ssize*lengthx/2*(1+df[i,facew])}
   if (faceh==0.5) {fh<-ssize*lengthy/2*(1+faceh)} else {fh<-ssize*lengthy/2*(1+df[i,faceh])}
   if (eyes==0.5) {es<-eyes} else {es<-df[i,eyes]+0.05}
   if (eyed==0.5) {ed<-eyed} else {ed<-df[i,eyed]}
   if (mouthw==0.5) {mw<-mouthw} else {mw<-df[i,mouthw]}
   if (mouthc==0.5) {mc<-mouthc} else {mc<-df[i,mouthc]}
   if (brows==0.5) {bs<-brows} else {bs<-df[i,brows]}
   if (browp==0.5) {bp<-browp} else {bp<-df[i,browp]}
   if (nosel==0.5) {nl<-nosel} else {nl<-df[i,nosel]}
   if (nosew==0.5) {nw<-nosew} else {nw<-df[i,nosew]}
   if (ears==0.5) {ea<-ears} else {ea<-df[i,ears]}
   if (pupils==0.5) {ps<-pupils} else {ps<-df[i,pupils]}
   # Ears
   plotellipse(mid=c(coorx[i]-fw,coory[i]),lwd=1,rx=ssize*lengthx*0.5*ea,ry=ssize*lengthy*0.5*ea,dr=0.1,lcol=border,col=col)
   plotellipse(mid=c(coorx[i]+fw,coory[i]),lwd=1,rx=ssize*lengthx*0.5*ea,ry=ssize*lengthy*0.5*ea,dr=0.1,lcol=border,col=col)
   # Face
   plotellipse(mid=c(coorx[i],coory[i]),lwd=1,rx=fw,ry=fh,dr=0.1,lcol=border,col=col)
   # Eyes
   ewx<-ssize*lengthx*es*0.25
   ewy<-ssize*lengthy*es*0.25
   ex<-ewx+(fw-2*ewx)*ed
   ey<-ssize*lengthy/7+(fh-(ssize*lengthy/2))/3
   plotellipse(mid=c(coorx[i]+ex,coory[i]+ey),lwd=1,rx=ewx,ry=ewy*0.8,dr=0.1,lcol=border,col=col)
   plotellipse(mid=c(coorx[i]-ex,coory[i]+ey),lwd=1,rx=ewx,ry=ewy*0.8,dr=0.1,lcol=border,col=col)
   # Pupils
   plotellipse(mid=c(coorx[i]+ex-ewx+ewx*0.25+(2*ewx-2*ewx*0.25)*ps,coory[i]+ey),lwd=1,rx=ewx*0.25,ry=ewy*0.25,dr=0.1,lcol=border,col=border)
   plotellipse(mid=c(coorx[i]-ex-ewx+ewx*0.25+(2*ewx-2*ewx*0.25)*ps,coory[i]+ey),lwd=1,rx=ewx*0.25,ry=ewy*0.25,dr=0.1,lcol=border,col=border)
   # Brows
   lines(c(coorx[i]+ex-ewx,coorx[i]+ex+ewx*1.1),c(coory[i]+ey+ewy*0.8+fh*0.4*bp-(ewy/2-ewy*bs),coory[i]+ey+ewy*0.8+fh*0.4*bp+(ewy/2-ewy*bs)),col=border)
   lines(c(coorx[i]-ex-ewx,coorx[i]-ex+ewx*1.1),c(coory[i]+ey+ewy*0.8+fh*0.4*bp+(ewy/2-ewy*bs),coory[i]+ey+ewy*0.8+fh*0.4*bp-(ewy/2-ewy*bs)),col=border)   
   # Mouth
   if (mc>=0.5){
    mc<-abs(mc-0.5)
    plotellipse(mid=c(coorx[i],coory[i]-fh/2+mc*ssize*lengthy*mw*0.7),rx=fw*mw*0.6,ry=mc*ssize*lengthy*mw,from=pi,to=pi*2,lwd=1,dr=0.1,lcol=border)
   } else
   {
    mc<-abs(mc-0.5)
    plotellipse(mid=c(coorx[i],coory[i]-fh/2-mc*ssize*lengthy*mw*0.7),rx=fw*mw*0.6,ry=mc*ssize*lengthy*mw,from=2*pi,to=3*pi,lwd=1,dr=0.1,lcol=border)
   }
   # Nose
   polygon(c(coorx[i]-ssize*lengthx*nw*0.18,coorx[i],coorx[i]+ssize*lengthx*nw*0.18),c(coory[i]-ssize*lengthy*nl*0.25,coory[i]+ssize*lengthy*nl*0.25,coory[i]-ssize*lengthy*nl*0.25),col=col,border=border)
  } else
  for (a in 1:m){
   if (is.numeric(df[,a])){
    if (type=="profile"){
     # Profiles
     polx[numattr+1-numberofnum]<-coorx[i]-lengthx*ssize+(numattr-numberofnum)*lengthx*ssize*2/(numattr-1)
     poly[numattr+1-numberofnum]<-coory[i]-lengthy*ssize/1.5+df[i,a]*lengthy*ssize*1.6
    } else
    if (type=="bar"){
     # Bars
     rect(coorx[i]-lengthx*ssize+(numattr-numberofnum)*lengthx*ssize*2/numattr,coory[i]-lengthy*ssize/1.5,coorx[i]-lengthx*ssize+(numattr-numberofnum)*lengthx*ssize*2/numattr+lengthx*ssize*2/numattr,coory[i]-lengthy*ssize/1.5+df[i,a]*lengthy*ssize*1.6,col=col,border=border)
    } else
    {
     # Stars, suns, polygons
     polx[numattr+1-numberofnum]<-coorx[i]+sin(angle*(numattr-numberofnum))*df[i,a]*lengthx*ssize
     poly[numattr+1-numberofnum]<-coory[i]+sin(pi/2-angle*(numattr-numberofnum))*df[i,a]*lengthy*ssize
    }
    numberofnum<-numberofnum-1
   }
  }
  if (type=="profile"){
   polx<-c(coorx[i]-lengthx*ssize,polx,coorx[i]+lengthx*ssize)
   poly<-c(coory[i]-lengthy*ssize/1.5,poly,coory[i]-lengthy*ssize/1.5)
  }
  if (type=="star" | type=="sun" | type=="polygon" | type=="profile") polygon(polx,poly,col=col,border=border)
  if (type=="sun"){
   for (a in 1:numattr){
    polx[a]<-coorx[i]+sin(angle*(a-1))*lengthx*ssize
    poly[a]<-coory[i]+sin(pi/2-angle*(a-1))*lengthy*ssize
   }
  }
  if (type=="star" | type=="sun") {for (a in 1:numattr) lines(c(coorx[i],polx[a]),c(coory[i],poly[a]),col=border,lty=lty)}
  if (!is.null(labels)){
   text(coorx[i],coory[i]-lengthy*ssize*1.2,labels[i],cex=labelsize,col=border)
  }
 }
}

