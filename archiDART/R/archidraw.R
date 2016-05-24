archidraw<-function(inputlie=NULL, inputrsml=NULL, res=NULL, unitlength="px", unitangle="d", rotation=0, rsml.connect=FALSE, numdate=NULL, finalscale=NULL, coldate=par("col"), ltydate=par("lty"), lwddate=par("lwd"), main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,...){
    
    # Errors interception
    
    if (is.null(inputlie)==TRUE & is.null(inputrsml)==TRUE){stop("inputlie and/or inputrsml must be provided")}
  
    if (is.null(inputlie)==FALSE) {if (mode(inputlie)!="character"){stop("mode(inputlie) must be character")}}
    
    if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
    
    if (is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
    if (is.null(res)==FALSE){
      if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
      if (res<=0){stop("res must be a positive value")}}
    
    if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
    if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
    
    if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
    if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
    
    if (mode(rotation)!="numeric"){stop("mode(rotation) must be numeric")}
    if (rotation<0){stop("rotation must be a positive value")}
    
    if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
    
    if (is.null(numdate)==FALSE)
    {if (mode(numdate)!="numeric"){stop("mode(numdate) must be numeric")}
     for (i in 1:length(numdate)){if (numdate[i]<=0){stop("numdate must be either a positive value or a vector of positive values")}}
     numdate.sort<-sort(numdate)
     for (i in 1:length(numdate)) {if (numdate[i]!=numdate.sort[i]){stop("Numeric elements in numdate must be sorted by increasing values")}}}
    
    if (is.null(finalscale)==FALSE) {if (mode(finalscale)!="logical"){stop("mode(finalscale) must be logical")}}
    
    # Reading of DART and rsml files
    
    if (is.null(inputlie)==FALSE){
      filenames.lie<-list.files(path=inputlie, pattern="\\.lie$")
      path.lie<-rep(inputlie, length.out=length(filenames.lie))
      filenameslie<-sub(x=filenames.lie, pattern="\\.lie$", replacement="")
      message(paste("Number of DART lie files in inputlie:", length(filenames.lie), sep=" "))}
    
    if (is.null(inputrsml)==FALSE) {
      filenames.rsml<-list.files(path=inputrsml, pattern="\\.rsml$")
      path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
      filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
      message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
    
    if (is.null(inputrsml)==TRUE){
      if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}}
    else {
      if (is.null(inputlie)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
      else{
        if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}
        if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}
    
    if (is.null(inputrsml)==TRUE){ # Only DART files

    LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)}
    
    else {
      
      if (is.null(inputlie)==TRUE){ # Only rsml files
        
        LIE<-list()
        filenameslie<-c()
        RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=1, connect=rsml.connect)
        for (i in 1:length(RSML)){
          LIE<-append(LIE, RSML[[i]]$lie)
          length1<-length(RSML[[i]]$lie)
          if (length1>1){
            num<-c(1:length1)
            filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
          if (length1==1){
            filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}}
      
      else { # DART and rsml files
        
        LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
        RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=1, connect=rsml.connect)
        for (i in 1:length(RSML)){
          LIE<-append(LIE, RSML[[i]]$lie)
          length1<-length(RSML[[i]]$lie)
          if (length1>1){
            num<-c(1:length1)
            filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
          if (length1==1){
            filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}}}
    
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
    
    # Drawing the root system architecture for each DART output file
  
    for (i in 1:length(LIE)){
  
    if (is.null(numdate)==TRUE){
      if (length(coldate)>max(LIE[[i]]$Date)){message(paste("Note: The number of colours in coldate is greater than max(Date) in ", filenameslie[i], ".lie", sep=""))}
      if (length(coldate)<max(LIE[[i]]$Date)){message(paste("Note: The number of colours in coldate is lower than max(Date) in ", filenameslie[i], ".lie", sep=""))}
      coldate1<-rep(coldate, length.out=max(LIE[[i]]$Date))
      if (length(ltydate)>max(LIE[[i]]$Date)){message(paste("Note: The number of elements in ltydate is greater than max(Date) in ", filenameslie[i], ".lie", sep=""))}
      if (length(ltydate)<max(LIE[[i]]$Date)){message(paste("Note: The number of elements in ltydate is lower than max(Date) in ", filenameslie[i], ".lie", sep=""))}
      ltydate1<-rep(ltydate, length.out=max(LIE[[i]]$Date))
      if (length(lwddate)>max(LIE[[i]]$Date)){message(paste("Note: The number of elements in lwddate is greater than max(Date) in ", filenameslie[i], ".lie", sep=""))}
      if (length(lwddate)<max(LIE[[i]]$Date)){message(paste("Note: The number of elements in lwddate is lower than max(Date) in ", filenameslie[i], ".lie", sep=""))}
      lwddate1<-rep(lwddate, length.out=max(LIE[[i]]$Date))
      LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
      LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
      minx<-min(LIE[[i]]$X)
      maxx<-max(LIE[[i]]$X)
      miny<-min(LIE[[i]]$Y)
      maxy<-max(LIE[[i]]$Y)
      if (is.null(main)==TRUE){main1<-filenameslie[i]} else {main1<-main}
      if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
      if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
      if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
      if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
      plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
      for (k in 1:nrow(LIE[[i]])){
        a<-LIE[[i]]$Prec[k]
        b<-LIE[[i]]$Date[k]
        if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b], lty=ltydate1[b], lwd=lwddate1[b],...)}}}
    
    if (is.null(numdate)==FALSE) {
      for (l in 1:length(numdate)){if (numdate[l]>max(LIE[[i]]$Date)){message(paste("Note: numdate contains a numerical value greater than max(Date) in ", filenameslie[i], ".lie", sep=""))}}
      if (length(coldate)>max(numdate)){message("Note: The number of colours in coldate is greater than max(numdate)")}
      if (length(coldate)<max(numdate)){message("Note: The number of colours in coldate is lower than max(numdate)")}
      coldate1<-rep(coldate, length.out=max(numdate))
      if (length(ltydate)>max(numdate)){message("Note: The number of elements in ltydate is greater than max(numdate)")}
      if (length(ltydate)<max(numdate)){message("Note: The number of elements in ltydate is lower than max(numdate)")}
      ltydate1<-rep(ltydate, length.out=max(numdate))
      if (length(lwddate)>max(numdate)){message("Note: The number of elements in lwddate is greater than max(numdate)")}
      if (length(lwddate)<max(numdate)){message("Note: The number of elements in lwddate is lower than max(numdate)")}
      lwddate1<-rep(lwddate, length.out=max(numdate))
      
      if (finalscale==TRUE){
        LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
        LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
        minx<-min(LIE[[i]]$X)
        maxx<-max(LIE[[i]]$X)
        miny<-min(LIE[[i]]$Y)
        maxy<-max(LIE[[i]]$Y)
        for (j in 1:length(numdate)){
          if (is.null(main)==TRUE){main1<-paste(filenameslie[i], "-numdate=", numdate[j], sep="")} else {main1<-main}
          if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
          if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
          if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
          if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
          plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, xlab=xlab1, ylab=ylab1,...)
          for (k in 1:nrow(LIE[[i]])){
            if (LIE[[i]]$Date[k]<=numdate[j]){
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k],col=coldate1[b],lty=ltydate1[b],lwd=lwddate1[b],...)}}}}}
      
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
          if (is.null(main)==TRUE){main1<-paste(filenameslie[i], "-numdate=", numdate[j], sep="")} else {main1<-main}
          if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
          if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
          if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
          if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
          plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, xlab=xlab1, ylab=ylab1,...)
          for (k in 1:nrow(LIE[[i]])){
            if (LIE[[i]]$Date[k]<=numdate[j]){
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k],col=coldate1[b],lty=ltydate1[b],lwd=lwddate1[b],...)}}}}}}}}