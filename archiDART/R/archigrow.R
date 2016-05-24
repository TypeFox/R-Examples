archigrow<-function(inputlie, inputtps, res=NULL, unitlength="px", unittime, unitangle="d", rotation=0, numdate=NULL, finalscale=NULL, coldyn, GRscale=NULL, main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, ...){
  
  # Errors interception
   
  if (mode(inputlie)!="character"){stop("mode(inputlie) must be character")}
  
  if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}
  
  if (is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
  if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
  if (res<=0){stop("res must be a positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (mode(unittime)!="character"){stop("mode(unittime) must be character")}
  
  if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
  if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
  
  if (mode(rotation)!="numeric"){stop("mode(rotation) must be numeric")}
  if (rotation<0){stop("rotation must be a positive value")}
  
  if (is.null(numdate)==FALSE)
  {if (mode(numdate)!="numeric"){stop("mode(numdate) must be numeric")}
   for (i in 1:length(numdate)){if (numdate[i]<=0){stop("numdate must be either a positive value or a vector of positive values")}}
   numdate.sort<-sort(numdate)
   for (i in 1:length(numdate)) {if (numdate[i]!=numdate.sort[i]){stop("Numeric elements in numdate must be sorted by increasing values")}}}
  
  if (is.null(finalscale)==FALSE) {if (mode(finalscale)!="logical"){stop("mode(finalscale) must be logical")}}
  
  if (is.null(GRscale)==FALSE){
    if (mode(GRscale)!="numeric") {stop("mode(GRscale) must be numeric")}
    if (length(GRscale)!=2|diff(GRscale)==0){stop("length(GRscale) must be equal to 2: c(min, max)")}
    if (GRscale[1]<0|GRscale[2]<0){stop("GRscale must be a vector of positive values")}}
  
  # Reading of DART output files
  
  filenames.tps<-list.files(path=inputtps, pattern="\\.tps$")
  path.tps<-rep(inputtps, length.out=length(filenames.tps))
  filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
  TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
  message(paste("Number of DART tps files in inputtps:", length(TIME), sep=" "))
  
  filenames.lie<-list.files(path=inputlie, pattern="\\.lie$")
  path.lie<-rep(inputlie, length.out=length(filenames.lie))
  filenameslie<-sub(x=filenames.lie, pattern="\\.lie$", replacement="")
  LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
  message(paste("Number of DART lie files in inputlie:", length(LIE), sep=" "))
  
  if (length(TIME)==1) {
    age<-list()
    for (i in 1:length(LIE)) {age[[i]]<-TIME[[1]]$Date}
    for (i in 1:length(LIE)) {if(length(age[[i]])!=max(LIE[[i]]$Date)){message(paste("Note: The number of observation dates in", filenames.tps[[1]], "is not equal to max(Date) in", filenames.lie[[i]], sep=" "))}}}
  else {
    if (length(TIME)!=length(LIE)) {stop("If there is more than one tps file in inputtps, the number of lie files in inputlie and tps files in inputtps must be equal")}
      for (i in 1:length(LIE)) {if(filenameslie[i]!=filenamestps[i]) {stop("Input lie files and their corresponding tps files must have the same name")}}
      age<-list()
      for (i in 1:length(LIE)){age[[i]]<-TIME[[i]]$Date}
      for (i in 1:length(LIE)) {if (length(age[[i]])!=max(LIE[[i]]$Date)) {message(paste("Note: The number of observation dates in", filenames.tps[[i]], "is not equal to max(Date) in", filenames.lie[[i]], sep=" "))}}}
  
  # Unit conversion and rotation
  
  if (unitangle=="r") {
    cunitangle<-1
    rotation<-rotation}
  if (unitangle=="d") {
    cunitangle<-180/pi
    rotation<-rotation*(1/cunitangle)}
  
  rot.matrix<-matrix(c(cos(rotation), sin(rotation), -sin(rotation), cos(rotation)), nrow=2, ncol=2)
  
  if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
  if (unitlength=="cm") {cunit<-(cm(1)/res)}
  if (unitlength=="px") {cunit<-1}
  for (i in 1:length(LIE)){
    LIE[[i]]$X<-LIE[[i]]$X*cunit
    LIE[[i]]$Y<-LIE[[i]]$Y*cunit
    newcoord<-rot.matrix%*%t(as.matrix(data.frame(LIE[[i]]$X, LIE[[i]]$Y)))
    LIE[[i]]$X<-newcoord[1,]
    LIE[[i]]$Y<-newcoord[2,]}
  
  GRunit<-paste("(", unitlength, "/", unittime, ")", sep="")
  
  # Calculation of a standardized growth rate matrix and a corresponding color matrix
  
  GRlist<-list()
  
  for (i in 1:length(LIE)){
    
    growthrate<-matrix(nrow=sum(LIE[[i]]$Bran=="true"), ncol=length(age[[i]]))
    r<-0
    
    for (l in 1:nrow(LIE[[i]])){
      
      if (LIE[[i]]$Bran[l]=="true"){
        
        if (LIE[[i]]$Date[l]==0) {
          r<-r+1}
        if (LIE[[i]]$Date[l]!=0) {
          r<-r+1
          c<-LIE[[i]]$Date[l]
          prec<-LIE[[i]]$Prec[l]
          if (c==1) {growthrate[r,c]<-sqrt((LIE[[i]]$X[l]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[l]-LIE[[i]]$Y[prec])^2)/age[[i]][c]}
          if (c>1) {growthrate[r,c]<-sqrt((LIE[[i]]$X[l]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[l]-LIE[[i]]$Y[prec])^2)/(age[[i]][c]-age[[i]][c-1])}}
      
        m<-LIE[[i]]$Suiv[l]
        
        if (m!=0){
        
        while (LIE[[i]]$Apic[m]=="false") {
          
          if (LIE[[i]]$Date[m]==LIE[[i]]$Date[LIE[[i]]$Prec[m]]){
          
              c<-LIE[[i]]$Date[m]
              prec<-LIE[[i]]$Prec[m]
              if (c==1) {growthrate[r,c]<-growthrate[r,c] + sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/age[[i]][c]}
              if (c>1) {growthrate[r,c]<-growthrate[r,c] + sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/(age[[i]][c]-age[[i]][c-1])}}
          
          if (LIE[[i]]$Date[m]!=LIE[[i]]$Date[LIE[[i]]$Prec[m]]){
            
            c<-LIE[[i]]$Date[m]
            prec<-LIE[[i]]$Prec[m]
            if (c==1) {growthrate[r,c]<-sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/age[[i]][c]}
            if (c>1) {growthrate[r,c]<-sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/(age[[i]][c]-age[[i]][c-1])}}
        
          m<-LIE[[i]]$Suiv[m]}
      
          if (LIE[[i]]$Apic[m]=="true"){
            
            if (LIE[[i]]$Date[m]==LIE[[i]]$Date[LIE[[i]]$Prec[m]]){
              
              c<-LIE[[i]]$Date[m]
              prec<-LIE[[i]]$Prec[m]
              if (c==1) {growthrate[r,c]<-growthrate[r,c] + sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/age[[i]][c]}
              if (c>1) {growthrate[r,c]<-growthrate[r,c] + sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/(age[[i]][c]-age[[i]][c-1])}}
            
            if (LIE[[i]]$Date[m]!=LIE[[i]]$Date[LIE[[i]]$Prec[m]]){
              
              c<-LIE[[i]]$Date[m]
              prec<-LIE[[i]]$Prec[m]
              if (c==1) {growthrate[r,c]<-sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/age[[i]][c]}
              if (c>1) {growthrate[r,c]<-sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)/(age[[i]][c]-age[[i]][c-1])}}}}}}
    
    growthrate[is.na(growthrate)]<-0
    if (is.null(GRscale)==TRUE) {
      maxi<-matrix(max(growthrate), nrow=nrow(growthrate), ncol=ncol(growthrate))
      mini<-matrix(min(growthrate), nrow=nrow(growthrate), ncol=ncol(growthrate))}
    else {
      maxi<-max(GRscale)
      mini<-min(GRscale)
      if (maxi<max(growthrate)){stop(paste("max(GRscale) must be positive and superior or equal to", max(growthrate), sep=" "))}
      if (mini>min(growthrate)){stop(paste("min(GRscale) must belong to [0,", min(growthrate), "]", sep=""))}}
    growthratesd<-(growthrate-mini)/(maxi-mini)
    pal<-colorRamp(coldyn)
    colorlegend<-colorRampPalette(coldyn)(30)
    colors<-matrix(nrow=nrow(growthratesd), ncol=ncol(growthratesd))
    for (c in 1:ncol(growthratesd)){
        colors[,c]<-rgb(pal(growthratesd[,c]), maxColorValue=255)}    
    
    # Creating an output list for the growth rate matrices
    
      Num<-c(1:nrow(growthrate)-1)
      outputresults<-data.frame(Num, growthrate)
      if (length(TIME)==1){num<-TIME[[1]]$Num} else {num<-TIME[[i]]$Num}
      colnames(outputresults)<-c("Root", paste("GR.Date", t(num), sep=""))
      GRlist[[i]]<-outputresults
      names(GRlist)[i]<-filenameslie[i]
    
    # X-Y plotting of vectorized root systems
    
    if (is.null(numdate)==TRUE){
      LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
      LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
      minx<-min(LIE[[i]]$X)
      maxx<-max(LIE[[i]]$X)
      miny<-min(LIE[[i]]$Y)
      maxy<-max(LIE[[i]]$Y)
      layout(matrix(2:1,ncol=2), widths = c(3,1),heights = c(1,1))
      legendimage<-as.raster(matrix(rev(colorlegend), ncol=1))
      par(mar=c(1,1,2,1))
      plot(c(0,2),c(0,1),type="n", axes=FALSE, xlab="", ylab="", main=paste("Growth rate",GRunit, sep="\n"), cex.main=0.7)
      if (is.null(GRscale)==TRUE) {text(x=1.5, y=seq(0.5,1,l=5), labels=round(seq(round(min(growthrate),1),round(max(growthrate),1),l=5),1),...)}
      else {text(x=1.5, y=seq(0.5,1,l=5), labels=round(seq(round(min(GRscale),1),round(max(GRscale),1),l=5),1),...)}
      rasterImage(legendimage,xleft=0.5,ybottom=0.5, xright=1,ytop=1)
      par(mar=c(5,4,4,2)+0.1)
      if (is.null(main)==TRUE){main1<-filenameslie[i]} else {main1<-main}
      if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
      if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
      if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
      if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
      plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
      r<-0
      for (k in 1:nrow(LIE[[i]])){
        if(LIE[[i]]$Bran[k]=="true"){
          r<-r+1
          a<-LIE[[i]]$Prec[k]
          b<-LIE[[i]]$Date[k]
          if (a!=0) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=colors[r,b],...)}
          m<-LIE[[i]]$Suiv[k]
          if (m!=0){
          while (LIE[[i]]$Apic[m]=="false"){
            a<-LIE[[i]]$Prec[m]
            b<-LIE[[i]]$Date[m]
            segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r,b],...)
            m<-LIE[[i]]$Suiv[m]}
          if (LIE[[i]]$Apic[m]=="true"){
            a<-LIE[[i]]$Prec[m]
            b<-LIE[[i]]$Date[m]
            segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r,b],...)}}}}
      layout(1)}
    
    if (is.null(numdate)==FALSE) {
      for (l in 1:length(numdate)){if (numdate[l]>max(LIE[[i]]$Date)){message(paste("Warning: numdate contains a numerical value greater than max(Date) in", filenames.lie[[i]], sep=" "))}}
      if (finalscale==TRUE){
        LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
        LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
        minx<-min(LIE[[i]]$X)
        maxx<-max(LIE[[i]]$X)
        miny<-min(LIE[[i]]$Y)
        maxy<-max(LIE[[i]]$Y)
        for (j in 1:length(numdate)){
          layout(matrix(2:1,ncol=2), widths = c(3,1),heights = c(1,1))
          legendimage<-as.raster(matrix(rev(colorlegend), ncol=1))
          par(mar=c(1,1,2,1))
          plot(c(0,2),c(0,1),type="n", axes=FALSE, xlab="", ylab="", main=paste("Growth rate",GRunit, sep="\n"), cex.main=0.7)
          if (is.null(GRscale)==TRUE) {text(x=1.5, y=seq(0.5,1,l=5), labels=round(seq(round(min(growthrate),1),round(max(growthrate),1),l=5),1),...)}
          else {text(x=1.5, y=seq(0.5,1,l=5), labels=round(seq(round(min(GRscale),1),round(max(GRscale),1),l=5),1),...)}
          rasterImage(legendimage,xleft=0.5,ybottom=0.5, xright=1,ytop=1)
          par(mar=c(5,4,4,2)+0.1)
          if (is.null(main)==TRUE){main1<-paste(filenameslie[i], "-numdate=", numdate[j], sep="")} else {main1<-main}
          if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
          if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
          if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
          if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
          plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, xlab=xlab1, ylab=ylab1,...)
          r<-0
          for (k in 1:nrow(LIE[[i]])){
            if(LIE[[i]]$Bran[k]=="true"){
              r<-r+1
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0 & b<=numdate[j]) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=colors[r,b],...)}
              m<-LIE[[i]]$Suiv[k]
              if (m!=0){
              while (LIE[[i]]$Apic[m]=="false"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                if (b<=numdate[j]) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r,b],...)}
                m<-LIE[[i]]$Suiv[m]}
              if (LIE[[i]]$Apic[m]=="true"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                if (b<=numdate[j]) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r,b],...)}}}}}
          layout(1)}}
      
      if (finalscale==FALSE){
        for (j in 1:length(numdate)){
          Date=NULL
          Num=NULL
          Suiv=NULL
          Prec=NULL
          X=NULL
          Y=NULL
          subset.LIE<-subset(LIE[[i]], Date<=numdate[j], select=c(Num, Date, Suiv, Prec, X, Y))
          LIE[[i]]$X<-LIE[[i]]$X-min(subset.LIE$X)
          LIE[[i]]$Y<-LIE[[i]]$Y-min(subset.LIE$Y)
          subset.LIE$X<-subset.LIE$X-min(subset.LIE$X)
          subset.LIE$Y<-subset.LIE$Y-min(subset.LIE$Y)
          minx<-min(subset.LIE$X)
          maxx<-max(subset.LIE$X)
          miny<-min(subset.LIE$Y)
          maxy<-max(subset.LIE$Y)
          layout(matrix(2:1,ncol=2), widths = c(3,1),heights = c(1,1))
          legendimage<-as.raster(matrix(rev(colorlegend), ncol=1))
          par(mar=c(1,1,2,1))
          plot(c(0,2),c(0,1),type="n", axes=FALSE, xlab="", ylab="", main=paste("Growth rate",GRunit, sep="\n"), cex.main=0.7)
          if (is.null(GRscale)==TRUE) {text(x=1.5, y=seq(0.5,1,l=5), labels=round(seq(round(min(growthrate),1),round(max(growthrate),1),l=5),1),...)}
          else {text(x=1.5, y=seq(0.5,1,l=5), labels=round(seq(round(min(GRscale),1),round(max(GRscale),1),l=5),1),...)}
          rasterImage(legendimage,xleft=0.5,ybottom=0.5, xright=1,ytop=1)
          par(mar=c(5,4,4,2)+0.1)
          if (is.null(main)==TRUE){main1<-paste(filenameslie[i], "-numdate=", numdate[j], sep="")} else {main1<-main}
          if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
          if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
          if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
          if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
          plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, xlab=xlab1, ylab=ylab1,...)
          r<-0
          for (k in 1:nrow(LIE[[i]])){
            if(LIE[[i]]$Bran[k]=="true"){
              r<-r+1
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0 & b<=numdate[j]) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=colors[r,b],...)}
              m<-LIE[[i]]$Suiv[k]
              if (m!=0){
              while (LIE[[i]]$Apic[m]=="false"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                if (b<=numdate[j]) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r,b],...)}
                m<-LIE[[i]]$Suiv[m]}
              if (LIE[[i]]$Apic[m]=="true"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                if (b<=numdate[j]) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r,b],...)}}}}}
          layout(1)}}}}
        
          return(GRlist)}