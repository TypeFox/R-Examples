# Creation of a function to calculate the distance between two points

distance<-function(x1,y1,x2,y2){
  a<-sqrt((x1-x2)^2+(y1-y2)^2)
  return(a)}

# Creation of a function to calculate the X and Y coordinates of an unknown point (Xn) on a line knowing two points (X1 and X2) 
# of that line and the distance between X1 and Xn

XYcoord<-function(x1,x2,y1,y2,d){
  if (x1!=x2){
    distx1x2<-distance(x1=x1, y1=y1, x2=x2, y2=y2)
    
    if (d>=(distx1x2/2)){
      dd<-d
      invslope<-(x2-x1)/(y2-y1)
      b<-(-2*x1)-((2*x1)/(invslope^2))
      a<-1+(1/(invslope^2))
      c<-(x1^2)+((x1^2)/(invslope^2))-(dd^2)
      delta<-(b^2)-4*a*c
      Xn1<-((-b)+sqrt(delta))/(2*a)
      Xn2<-((-b)-sqrt(delta))/(2*a)
      Yn1<-(1/invslope)*Xn1+(y1-(1/invslope)*x1)
      Yn2<-(1/invslope)*Xn2+(y1-(1/invslope)*x1)
      if (Xn1>=min(x1,x2) & Xn1<=max(x1,x2)){result<-c(Xn1,Yn1)} else {result<-c(Xn2,Yn2)}}
    
    if (d<(distx1x2/2)){
      dd<-distx1x2-d
      invslope<-(x2-x1)/(y2-y1)
      b<-(-2*x2)-((2*x2)/(invslope^2))
      a<-1+(1/(invslope^2))
      c<-(x2^2)+((x2^2)/(invslope^2))-(dd^2)
      delta<-(b^2)-4*a*c
      Xn1<-((-b)+sqrt(delta))/(2*a)
      Xn2<-((-b)-sqrt(delta))/(2*a)
      Yn1<-(1/invslope)*Xn1+(y2-(1/invslope)*x2)
      Yn2<-(1/invslope)*Xn2+(y2-(1/invslope)*x2)
      if (Xn1>=min(x1,x2) & Xn1<=max(x1,x2)){result<-c(Xn1,Yn1)} else {result<-c(Xn2,Yn2)}}
    
    return(result)}
  
  else{
    Xn1<-x1
    Yn1<-y1+d
    Yn2<-y1-d
    if (Yn1>=min(y1,y2) & Yn1<=max(y1,y2)){result<-c(Xn1,Yn1)} else {result<-c(Xn1,Yn2)}
    return(result)}}

# Calculation of a normal vector to a line (ax+by+c=0)
# u is a direction vector of ax+by+c=0

normal<-function(u){
  n<-c(u[2], -u[1])
  return(n)}

################################################################################################

trajectory<-function(inputrac=NULL, inputlie=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px", unitangle="d", rotation=0, l.brangle, l.curv, l.tipangle, rsml.date=NULL){
  
  # Errors interception
  
  if (is.null(inputrac)==TRUE & is.null(inputlie)==TRUE & is.null(inputtps)==TRUE & is.null(inputrsml)==TRUE){stop("inputrac/inpulie/inputtps and/or inputrsml must be provided")}
  
  if (is.null(inputrac)==FALSE) {if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}}
  
  if (is.null(inputlie)==FALSE) {if (mode(inputlie)!="character"){stop("mode(inputlie) must be character")}}
  
  if (is.null(inputtps)==FALSE) {if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}}
  
  if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
  
  if (is.null(inputrac)==FALSE|is.null(inputtps)==FALSE|is.null(inputlie)==FALSE){
    if (is.null(inputrac)==TRUE|is.null(inputtps)==TRUE|is.null(inputlie)==TRUE){stop("If inputrac/inputlie/inputtps is not NULL, inputrac/inputlie/inputtps must be provided")}}
  
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
  
  if (mode(l.brangle)!="numeric"){stop("mode(l.brangle) must be numeric")}
  if (l.brangle<=0) {stop("l.brangle must be a positive value")}
  
  if (mode(l.curv)!="numeric"){stop("mode(l.curv) must be numeric")}
  if (l.curv<=0) {stop("l.curv must be a positive value")}
  
  if (mode(l.tipangle)!="numeric"){stop("mode(l.tipangle) must be numeric")}
  if (l.tipangle<=0) {stop("l.tipangle must be a positive value")}
  
  if (is.null(inputrsml)==FALSE & is.null(rsml.date)==TRUE) {stop("If inputrsml is not NULL, rsml.date must be a positive numeric value")}
  
  if (is.null(rsml.date)==FALSE){
    if (rsml.date<=0|length(rsml.date)>1){stop("rsml.date must be a single positive value")}}
  
  # Reading of DART and rsml files
  
  if (is.null(inputrac)==FALSE){
    filenames.rac<-list.files(path=inputrac, pattern="\\.rac$")
    path.rac<-rep(inputrac, length.out=length(filenames.rac))
    filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
    message(paste("Number of DART rac files in inputrac:", length(filenames.rac), sep=" "))}
  
  if (is.null(inputtps)==FALSE){
    filenames.tps<-list.files(path=inputtps, pattern="\\.tps$")
    path.tps<-rep(inputtps, length.out=length(filenames.tps))
    filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
    message(paste("Number of DART tps files in inputtps:", length(filenames.tps), sep=" "))}
  
  if (is.null(inputlie)==FALSE){
    filenames.lie<-list.files(path=inputlie, pattern="\\.lie$")
    path.lie<-rep(inputlie, length.out=length(filenames.lie))
    filenameslie<-sub(x=filenames.lie, pattern="\\.lie$", replacement="")
    message(paste("Number of DART lie files in inputlie:", length(filenames.lie), sep=" "))}
  
  if (is.null(inputrsml)==FALSE){
    filenames.rsml<-list.files(path=inputrsml, pattern="\\.rsml$")
    path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
    filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
    message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
  
  if (is.null(inputrsml)==TRUE){
    if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
    if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
    if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}}
  else {
    if (is.null(inputrac)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
    else{
      if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
      if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
      if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}
      if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}
  
  if (is.null(inputrsml)==TRUE){ # Only DART files
    
    LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
    
    DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
    for (i in 1:length(DATA)) {
      colnames(DATA[[i]])<-c()
      colnames(DATA[[i]])[1]<-"Root"
      colnames(DATA[[i]])[2]<-"Mother"
      colnames(DATA[[i]])[3]<-"Ord"
      colnames(DATA[[i]])[4]<-"DBase"
      colnames(DATA[[i]])[5]<-"DApp"
      for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
    
    TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
  
    if (length(LIE)!=length(DATA)) {stop("The number of rac files in inputrac and lie files in inputlie must be equal")}
    else {
      for (i in 1:length(DATA)) {if(filenamesrac[i]!=filenameslie[i]) {stop("Input rac files and their corresponding lie files must have the same name")}}}  	
    
    if (length(TIME)==1) {
      for (i in 1:length(DATA)) {if(length(TIME[[1]]$Date)!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
    else {
      if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac/lie files in inputrac/inputlie and tps files in inputtps must be equal")}
      else {
        for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac/lie files and their corresponding tps files must have the same name")}}
        for (i in 1:length(DATA)) {if (length(TIME[[i]]$Date)!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}} 
  
  else {
    
    if (is.null(inputrac)==TRUE){ # Only rsml files
      
      LIE<-list()
      DATA<-list()
      TIME<-list()
      time<-list(data.frame(Num=1, Date=rsml.date, CoulR=0, CoulG=0, CoulB=0))
      filenameslie<-c()
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=TRUE)
      for (i in 1:length(RSML)){
        for (j in 1:length(RSML[[i]]$rac)){colnames(RSML[[i]]$rac[[j]])[6]<-"Lengths1"}
        DATA<-append(DATA, RSML[[i]]$rac)
        LIE<-append(LIE, RSML[[i]]$lie)
        length1<-length(RSML[[i]]$rac)
        TIME[(length(TIME)+1):(length(TIME)+length1)]<-time[1]
        if (length1>1){
          num<-c(1:length1)
          filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}}
    
    else { # DART and rsml files
      
      LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
      
      if (length(LIE)!=length(DATA)) {stop("The number of rac files in inputrac and lie files in inputlie must be equal")}
      else {
        for (i in 1:length(DATA)) {if(filenamesrac[i]!=filenameslie[i]) {stop("Input rac files and their corresponding lie files must have the same name")}}}  	
      
      if (length(TIME)==1) {
        for (i in 1:length(DATA)) {if(length(TIME[[1]]$Date)!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}
        TIME[1:length(DATA)]<-TIME[1]}
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac/lie files in inputrac/inputlie and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac/lie files and their corresponding tps files must have the same name")}}
          for (i in 1:length(DATA)) {if (length(TIME[[i]]$Date)!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}
      
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=TRUE)
      
      time<-list(data.frame(Num=1, Date=rsml.date, CoulR=0, CoulG=0, CoulB=0))
      
      for (i in 1:length(RSML)){
        for (j in 1:length(RSML[[i]]$rac)){colnames(RSML[[i]]$rac[[j]])[6]<-"Lengths1"}
        DATA<-append(DATA, RSML[[i]]$rac)
        LIE<-append(LIE, RSML[[i]]$lie) 
        length1<-length(RSML[[i]]$rac)
        TIME[(length(TIME)+1):(length(TIME)+length1)]<-time[1]
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
    
  # Creating vectors and matrices for root architecture parameters calculation
  
    filenames<-c()
    finallength.lie<-list()
    Root<-c()
    Mother<-c()
    Ord<-c()
    DBase<-c()
    DApp<-c()
    Finalrootlength<-c()
    Orientation<-c()
    Br.Angle<-c()
    MeanAngleVar<-c()
    SDAngleVar<-c()
    rac<-list()
    tip<-list()
  
  # Storing the final root length of each root constituting the vectorized root systems
  
  for (i in 1:length(LIE)){
    
    finallength<-c()
    k<-0
    
    for (j in 1:nrow(LIE[[i]])){
      
      if (LIE[[i]]$Bran[j]=="true"){
        k<-k+1
        if (j==1){
          finallength[k]<-0
          m<-LIE[[i]]$Suiv[j]}
        
        if (j>1){
          prec<-LIE[[i]]$Prec[j]
          finallength[k]<-sqrt((LIE[[i]]$X[j]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[j]-LIE[[i]]$Y[prec])^2)
          m<-LIE[[i]]$Suiv[j]}
        
          while (m!=0){
            prec<-LIE[[i]]$Prec[m]
            finallength[k]<-finallength[k]+sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)
            m<-LIE[[i]]$Suiv[m]}}}
    
    finallength.lie[[i]]<-finallength}
    
  # Calculating the X and Y coordinates of interpolated points
  
  t<-1
  for (i in 1:length(LIE)){
    
    k<-0
    XYcurv<-list()
    XYangle<-list()
    orientation<-c()
    tortuosity<-c()
    if (length(TIME)==1){
      num<-TIME[[1]]$Num
      date<-TIME[[1]]$Date}
    else {
      num<-TIME[[i]]$Num
      date<-TIME[[i]]$Date}
    tipangle<-matrix(nrow=nrow(DATA[[i]]), ncol=length(date))
    
    for (j in 1:nrow(LIE[[i]])){
      
      if (LIE[[i]]$Bran[j]=="true"){
        
        k<-k+1
        l<-1
        distangle<-0
        if (finallength.lie[[i]][k]<l.tipangle) {tipangle[k,]<-NA}
        
        # Points used for root curvature calculation
        
        if (finallength.lie[[i]][k]<l.curv){XYcurv[[k]]<-NA}
        else {
        XYcurv[[k]]<-matrix(ncol=2, nrow=floor(finallength.lie[[i]][k]/l.curv)+1)
        if (j==1){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[j]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[j]
              suiv<-LIE[[i]]$Suiv[j]
              l<-l+1}
        
        if (j>1){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[LIE[[i]]$Prec[j]]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[LIE[[i]]$Prec[j]]
              suiv<-j
              l<-l+1}
                
        D<-distance(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
        distcurv<-D
        
        while (l<=floor(finallength.lie[[i]][k]/l.curv)+1){
          
          # First situation
          if (distcurv<l.curv){
            while (distcurv<l.curv){
              suiv<-LIE[[i]]$Suiv[suiv]
              D<-distance(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
              distcurv<-distcurv+D}
            
            if (distcurv==l.curv){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              distcurv<-0}
            
            if (distcurv>l.curv){
              xycoord<-XYcoord(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv-distcurv+D)
              XYcurv[[k]][l,1]<-xycoord[1]
              XYcurv[[k]][l,2]<-xycoord[2]
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              distcurv<-distance(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
              
              if (distcurv==l.curv){
                XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                l<-l+1
                if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                distcurv<-0}
              
              while (distcurv>l.curv){
                xycoord<-XYcoord(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv)
                XYcurv[[k]][l,1]<-xycoord[1]
                XYcurv[[k]][l,2]<-xycoord[2]
                l<-l+1
                if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                distcurv<-distance(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
                if (distcurv==l.curv){
                  XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                  XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                  l<-l+1
                  if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                  distcurv<-0}}
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}}
            suiv<-LIE[[i]]$Suiv[suiv]}
          
          # Second situation
          if (distcurv>l.curv){
            xycoord<-XYcoord(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv-distcurv+D)
            XYcurv[[k]][l,1]<-xycoord[1]
            XYcurv[[k]][l,2]<-xycoord[2]
            l<-l+1
            if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
            distcurv<-distance(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
            
            if (distcurv==l.curv){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              distcurv<-0}
            
            while (distcurv>l.curv){
              xycoord<-XYcoord(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv)
              XYcurv[[k]][l,1]<-xycoord[1]
              XYcurv[[k]][l,2]<-xycoord[2]
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              distcurv<-distance(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
              if (distcurv==l.curv){
                XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                l<-l+1
                if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                distcurv<-0}}
            
            if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
            
          suiv<-LIE[[i]]$Suiv[suiv]}
          
          # Third situation
          if (distcurv==l.curv){
            XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
            XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
            l<-l+1
            if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
            distcurv<-0
            suiv<-LIE[[i]]$Suiv[suiv]}
        
          D<-distance(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])
          distcurv<-distcurv+D}}
      
      # Points used to calculate the branching angles of daughter roots on mother roots
        
        if (j==1){
          XYangle[[k]]<-NA
          orientation[k]<-NA}
        else{
          if (finallength.lie[[i]][k]<l.brangle|(finallength.lie[[i]][DATA[[i]]$Mother[k]+1]-(cunit*DATA[[i]]$DBase[k]))<l.brangle){XYangle[[k]]<-NA}
          else{
            XYangle[[k]]<-matrix(nrow=3, ncol=2)
            XYangle[[k]][1,1]<-LIE[[i]]$X[LIE[[i]]$Prec[j]]
            XYangle[[k]][1,2]<-LIE[[i]]$Y[LIE[[i]]$Prec[j]]
            
            # For daughter roots
            m<-j
            
            while(distangle<l.brangle){
              prec<-LIE[[i]]$Prec[m]
              D<-distance(x1=LIE[[i]]$X[prec], x2=LIE[[i]]$X[m], y1=LIE[[i]]$Y[prec], y2=LIE[[i]]$Y[m])
              distangle<-distangle+D
              if (distangle<l.brangle) {m<-LIE[[i]]$Suiv[m]}}
            
            if (distangle==l.brangle){
              XYangle[[k]][3,1]<-LIE[[i]]$X[m]
              XYangle[[k]][3,2]<-LIE[[i]]$Y[m]}
            
            if (distangle>l.brangle){
              xycoord<-XYcoord(x1=LIE[[i]]$X[prec], x2=LIE[[i]]$X[m], y1=LIE[[i]]$Y[prec], y2=LIE[[i]]$Y[m], d=l.brangle-distangle+D)
              XYangle[[k]][3,1]<-xycoord[1]
              XYangle[[k]][3,2]<-xycoord[2]}
            
            # For mother roots
            distangle<-0
            m<-LIE[[i]]$Prec[j]
            
            while(distangle<l.brangle){
              suiv<-LIE[[i]]$Suiv[m]
              D<-distance(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[suiv], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[suiv])
              distangle<-distangle+D
              if (distangle<l.brangle) {m<-suiv}}
            
            if (distangle==l.brangle){
              XYangle[[k]][2,1]<-LIE[[i]]$X[suiv]
              XYangle[[k]][2,2]<-LIE[[i]]$Y[suiv]}
            
            if (distangle>l.brangle){
              xycoord<-XYcoord(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[suiv], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[suiv], d=l.brangle-distangle+D)
              XYangle[[k]][2,1]<-xycoord[1]
              XYangle[[k]][2,2]<-xycoord[2]}}
          
          # Calculating the Orientation of lateral roots
          
          prec<-LIE[[i]]$Prec[j]
          suivMR<-LIE[[i]]$Suiv[prec]
          suivLR<-LIE[[i]]$Suiv[j]
          u<-c(LIE[[i]]$X[suivMR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivMR]-LIE[[i]]$Y[prec])
          while (u[1]==0 & u[2]==0){
            suivMR<-LIE[[i]]$Suiv[suivMR]
            u<-c(LIE[[i]]$X[suivMR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivMR]-LIE[[i]]$Y[prec])}
          n<-normal(u)
          lateral<-c(LIE[[i]]$X[j]-LIE[[i]]$X[prec], LIE[[i]]$Y[j]-LIE[[i]]$Y[prec])
          while (lateral[1]==0 & lateral[2]==0){
            lateral<-c(LIE[[i]]$X[suivLR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivLR]-LIE[[i]]$Y[prec])
            suivLR<-LIE[[i]]$Suiv[suivLR]}
          
          if (n%*%lateral>0){orientation[k]<-"Left"}
          if (n%*%lateral<0){orientation[k]<-"Right"}
          
          if (n%*%lateral==0 & suivLR==0){
            if (runif(1, min=-1, max=1)>0) {orientation[k]<-"Left"} else {orientation[k]<-"Right"}}
          
          if (n%*%lateral==0 & suivLR!=0){
            
            while(n%*%lateral==0){
              lateral<-c(LIE[[i]]$X[suivLR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivLR]-LIE[[i]]$Y[prec])
              suivLR<-LIE[[i]]$Suiv[suivLR]
              if (suivLR==0 & n%*%lateral==0){
                if (runif(1, min=-1, max=1)>0) {orientation[k]<-"Left"} else {orientation[k]<-"Right"}
                break}}
            
            if (n%*%lateral>0){orientation[k]<-"Left"}
            if (n%*%lateral<0){orientation[k]<-"Right"}}}
      
      # Calculating the points used for the calculation of tip angles
      
      if (finallength.lie[[i]][k]>=l.tipangle){
       
        m<-j
        
        while (LIE[[i]]$Apic[m]!="true"){
          
          a<-LIE[[i]]$Date[m]
          b<-LIE[[i]]$Date[LIE[[i]]$Suiv[m]]
          
          if (a!=0 & a<b){
            
            if (DATA[[i]][k,5+a]*cunit>=l.tipangle) {
            
            disttip<-distance(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]])
            if (disttip>l.tipangle){xycoord<-XYcoord(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], d=l.tipangle)}
            if (disttip==l.tipangle){xycoord<-c(LIE[[i]]$X[LIE[[i]]$Prec[m]], LIE[[i]]$Y[LIE[[i]]$Prec[m]])}
            if (disttip<l.tipangle){
              prec1<-LIE[[i]]$Prec[m]
              while(disttip<l.tipangle){
                prec2<-LIE[[i]]$Prec[prec1]
                D<-distance(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2])
                disttip<-disttip+D
                if (disttip<l.tipangle) {prec1<-prec2}}
              if (disttip==l.tipangle){xycoord<-c(LIE[[i]]$X[prec2], LIE[[i]]$Y[prec2])}
              if (disttip>l.tipangle) {xycoord<-XYcoord(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], d=l.tipangle-disttip+D)}}
            
            tipangle[k,a:(b-1)]<-acos((c(0,1)%*%c(LIE[[i]]$X[m]-xycoord[1], LIE[[i]]$Y[m]-xycoord[2]))/(sqrt((LIE[[i]]$X[m]-xycoord[1])^2+(LIE[[i]]$Y[m]-xycoord[2])^2)))*cunitangle}}
          
          m<-LIE[[i]]$Suiv[m]}
      
      if (LIE[[i]]$Apic[m]=="true"){
                
          a<-LIE[[i]]$Date[m]
        
          disttip<-distance(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]])
          if (disttip>l.tipangle){xycoord<-XYcoord(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], d=l.tipangle)}
          if (disttip==l.tipangle){xycoord<-c(LIE[[i]]$X[LIE[[i]]$Prec[m]], LIE[[i]]$Y[LIE[[i]]$Prec[m]])}
          if (disttip<l.tipangle){
            prec1<-LIE[[i]]$Prec[m]
            while(disttip<l.tipangle){
              prec2<-LIE[[i]]$Prec[prec1]
              D<-distance(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2])
              disttip<-disttip+D
              if (disttip<l.tipangle) {prec1<-prec2}}
            if (disttip==l.tipangle){xycoord<-c(LIE[[i]]$X[prec2], LIE[[i]]$Y[prec2])}
            if (disttip>l.tipangle) {xycoord<-XYcoord(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], d=l.tipangle-disttip+D)}}
          
          tipangle[k,(a:ncol(tipangle))]<-acos((c(0,1)%*%c(LIE[[i]]$X[m]-xycoord[1], LIE[[i]]$Y[m]-xycoord[2]))/(sqrt((LIE[[i]]$X[m]-xycoord[1])^2+(LIE[[i]]$Y[m]-xycoord[2])^2)))*cunitangle}}
      
      # Tortuosity
      
      m<-j
      
      while (LIE[[i]]$Apic[m]!="true"){m<-LIE[[i]]$Suiv[m]}
      
      if (k==1) {tortuosity[k]<-finallength.lie[[i]][k]/sqrt((LIE[[i]]$X[j]-LIE[[i]]$X[m])^2+(LIE[[i]]$Y[j]-LIE[[i]]$Y[m])^2)} else {tortuosity[k]<-finallength.lie[[i]][k]/sqrt((LIE[[i]]$X[LIE[[i]]$Prec[j]]-LIE[[i]]$X[m])^2+(LIE[[i]]$Y[LIE[[i]]$Prec[j]]-LIE[[i]]$Y[m])^2)}}}
    
    tip[[i]]<-data.frame(DATA[[i]]$Root, tipangle)
    colnames(tip[[i]])<-c("Root", paste("Ang.Date", t(num), sep=""))
  
  # Calculating branching angles and curvatures
  
    br.angle<-c()
    meananglevar<-c()
    sdanglevar<-c()
    
    for (j in 1:length(XYangle)){
      
      if (class(XYangle[[j]])=="matrix"){
        
        VECTangle<-matrix(nrow=2, ncol=2)
        VECTangle[1,]<-XYangle[[j]][2,]-XYangle[[j]][1,] # For mother roots
        VECTangle[2,]<-XYangle[[j]][3,]-XYangle[[j]][1,] # For daugther roots
        normVECTangle<-sqrt(VECTangle[,1]^2+VECTangle[,2]^2)
        br.angle[j]<-acos((VECTangle[1,]%*%VECTangle[2,])/(normVECTangle[1]*normVECTangle[2]))*cunitangle} 
      
      else {br.angle[j]<-NA}
      
      if (class(XYcurv[[j]])=="matrix") {
        
        if (nrow(XYcurv[[j]])>2){
          
          VECTcurv<-diff(XYcurv[[j]])
          angle<-c()
          for (k in 1:(nrow(XYcurv[[j]])-2)){
            ratio<-(VECTcurv[k,]%*%VECTcurv[k+1,])/(sqrt(VECTcurv[k,1]^2+VECTcurv[k,2]^2)*sqrt(VECTcurv[k+1,1]^2+VECTcurv[k+1,2]^2))
            if (ratio>1) {ratio<-1} # The value of a cosinus must be < or equal to 1
            if (ratio<(1*(-1))) {ratio<-1*(-1)} # The value of a cosinus must be > or equal to -1
            angle[k]<-acos(ratio)*cunitangle}
          meananglevar[j]<-mean(angle)
          sdanglevar[j]<-sd(angle)}
      
        else {
          
          meananglevar[j]<-NA
          sdanglevar[j]<-NA}} 
      
      else {
        
        meananglevar[j]<-NA
        sdanglevar[j]<-NA}}
  
  rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Mother=DATA[[i]]$Mother, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase*cunit, DApp=DATA[[i]]$DApp, FinalRootLength=finallength.lie[[i]], Tortuosity=tortuosity, Orientation=orientation, Branching.Angle=br.angle, Mean.Curv=meananglevar, SD.Curv=sdanglevar)}
  
  names(rac)<-filenameslie
  names(tip)<-filenameslie
  outputresults<-list(root=rac, tip=tip)
  return(outputresults)}