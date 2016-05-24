architect<-function(inputrac=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px", rsml.date=NULL, rsml.connect=FALSE, rootdiv=1){
  
  # Errors interception
  
  if (is.null(inputrac)==TRUE & is.null(inputtps)==TRUE & is.null(inputrsml)==TRUE){stop("inputrac/inputps and/or inputrsml must be provided")}
  
  if (is.null(inputrac)==FALSE) {if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}}
  
  if (is.null(inputtps)==FALSE) {if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}}
  
  if (is.null(inputrac)==FALSE|is.null(inputtps)==FALSE){
    if (is.null(inputrac)==TRUE|is.null(inputtps)==TRUE){stop("If inputrac/inputtps is not NULL, inputtps/inputrac must be provided")}}
  
  if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
  
  if (is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0|length(res)>1){stop("res must be a single positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(inputrsml)==FALSE & is.null(rsml.date)==TRUE) {stop("If inputrsml is not NULL, rsml.date must be a positive numeric value")}
  
  if (is.null(rsml.date)==FALSE){
    if (rsml.date<=0|length(rsml.date)>1){stop("rsml.date must be a single positive value")}}
  
  if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
  
  if (mode(rootdiv)!="numeric"){stop("mode(rootdiv) must be numeric")}
  for (i in 1:length(rootdiv)){if (rootdiv[i]<0){stop("rootdiv must be either a positive value or a vector of positive values")}}
  rootdiv.sort<-sort(rootdiv)
  for (i in 1:length(rootdiv)) {if (rootdiv[i]!=rootdiv.sort[i]){stop("Numeric elements in rootdiv must be sorted by increasing values")}}
  
  # Reading of DART and rsml files
  
  if (is.null(inputtps)==FALSE){
    filenames.tps<-list.files(path=inputtps, pattern="\\.tps$")
    path.tps<-rep(inputtps, length.out=length(filenames.tps))
    filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
    message(paste("Number of DART tps files in inputtps:", length(filenames.tps), sep=" "))}
  
  if (is.null(inputrac)==FALSE){
    filenames.rac<-list.files(path=inputrac, pattern="\\.rac$")
    path.rac<-rep(inputrac, length.out=length(filenames.rac))
    filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
    message(paste("Number of DART rac files in inputrac:", length(filenames.rac), sep=" "))}
  
  if (is.null(inputrsml)==FALSE) {
    filenames.rsml<-list.files(path=inputrsml, pattern="\\.rsml$")
    path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
    filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
    message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
  
  if (is.null(inputrsml)==TRUE){
    if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
    if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}}
  else {
    if (is.null(inputrac)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
    else{
      if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
      if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
      if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}

  if (is.null(inputrsml)==TRUE){ # Only DART files
  
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      if (length(TIME)==1) {
        age<-list()
        obstot<-0
        for (i in 1:length(DATA)) {
          age[[i]]<-TIME[[1]]$Date
          obstot<-obstot+length(age[[i]])}
        for (i in 1:length(DATA)) {if(length(age[[i]])!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac files in inputrac and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac files and their corresponding tps files must have the same name")}}
          age<-list()
          obstot<-0
          for (i in 1:length(TIME)) {
            age[[i]]<-TIME[[i]]$Date
            obstot<-obstot+length(age[[i]])}
          for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}}
  
  else {
    if (is.null(inputrac)==TRUE){ # Only RSML files
      
      DATA<-list()
      age<-list()
      filenamesrac<-c()
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect) # Read RSML files
      obstot<-0
      for (i in 1:length(RSML)){
        for (j in 1:length(RSML[[i]]$rac)){colnames(RSML[[i]]$rac[[j]])[6]<-"Lengths1"}
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        obstot<-obstot+length1
        age[(length(age)+1):(length(age)+length1)]<-rsml.date
        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}}
    
    else { # DART and RSML files
      
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE) # Read TPS files
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1) # Read RAC files
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      if (length(TIME)==1) {
        age<-list()
        obstot<-0
        for (i in 1:length(DATA)) {
          age[[i]]<-TIME[[1]]$Date
          obstot<-obstot+length(age[[i]])}
        for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac files in inputrac and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac files and their corresponding tps files must have the same name")}}
          age<-list()
          obstot<-0
          for (i in 1:length(DATA)) {
            age[[i]]<-TIME[[i]]$Date
            obstot<-obstot+length(age[[i]])}
          for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}
      
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect) # Read RSML files
      
      for (i in 1:length(RSML)){
        for (j in 1:length(RSML[[i]]$rac)){colnames(RSML[[i]]$rac[[j]])[6]<-"Lengths1"}
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        obstot<-obstot+length1
        age[(length(age)+1):(length(age)+length1)]<-rsml.date
        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}}}
  
  # Creating vectors and matrices for root architecture parameters calculation
  
  FileNames<-c()
  Time<-c()
  FirstOrderRootLength<-c()
  TotalRootLength<-c()
  TotalLateralRootNumber<-c()
  TotalLateralRootLength<-c()
  LateralRootDensity<-c()
  GrowthRateFirstOrderRoot<-c()
  GrowthRateTotalRoot<-c()
  
  if ((length(rootdiv)==1 & rootdiv[1]!=1)|(length(rootdiv)>1)) {
    ColLengthsNumber<-0
    for (i in 1:length(DATA)) {ColLengthsNumber<-ColLengthsNumber+(ncol(DATA[[i]])-5)}
    
    if (length(rootdiv)==1){
      LRNDroot<-matrix(nrow=ColLengthsNumber,ncol=rootdiv)
      LRLDroot<-matrix(nrow=ColLengthsNumber,ncol=rootdiv)
      LRDDroot<-matrix(nrow=ColLengthsNumber,ncol=rootdiv)}
    else {
      LRNDroot<-matrix(nrow=ColLengthsNumber,ncol=length(rootdiv)-1)
      LRLDroot<-matrix(nrow=ColLengthsNumber,ncol=length(rootdiv)-1)
      LRDDroot<-matrix(nrow=ColLengthsNumber,ncol=length(rootdiv)-1)}}		
  
  # Calculation of root architecture parameters
  
  k<-0
  maxord<-max(sapply(DATA, function(x){max(x[[3]])}))
  if (maxord>1) {latroot<-matrix(ncol=4*(maxord-1), nrow=obstot)}
  
  for(i in 1:length(DATA)){
    
    LastFirstOrderRootLength<-DATA[[i]][[ncol(DATA[[i]])]][DATA[[i]]$Ord==1]
    
    for (t in 1:(ncol(DATA[[i]])-5)){
      
      k<-k+1
      
      data.split<-split(DATA[[i]][[paste("Lengths",t,sep="")]],as.factor(DATA[[i]]$Ord))
      
      # File names
      
      FileNames[k]<-filenamesrac[i]
      
      # Time
      Time[k]<-age[[i]][t]
      
      # Unit conversion
      
      if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
      if (unitlength=="cm") {cunit<-(cm(1)/res)}
      if (unitlength=="px") {cunit<-1}	
      
      # Total root length
      
      TotalRootLength[k]<-sum(DATA[[i]][[paste("Lengths",t,sep="")]])*cunit
      
      # Growth rate of the root system
      
      if (t==1) {GrowthRateTotalRoot[k]<-TotalRootLength[k]/age[[i]][t]}
      if (t>1) {GrowthRateTotalRoot[k]<-(TotalRootLength[k]-TotalRootLength[k-1])/(age[[i]][t]-age[[i]][t-1])}
      
      # First-order root length
      
      FirstOrderRootLength[k]<-sum(data.split$'1')*cunit
      
      # Growth rate of the first-order root
      
      if (t==1) {GrowthRateFirstOrderRoot[k]<-FirstOrderRootLength[k]/age[[i]][t]}
      if (t>1) {GrowthRateFirstOrderRoot[k]<-(FirstOrderRootLength[k]-FirstOrderRootLength[k-1])/(age[[i]][t]-age[[i]][t-1])}
      
      # Total number of lateral roots
      
      if (FirstOrderRootLength[k]==0) {TotalLateralRootNumber[k]<-0}  else {TotalLateralRootNumber[k]<-sum(DATA[[i]][[paste("Lengths",t,sep="")]]!=0)-length(data.split$'1')}
      
      # Total length of lateral roots
      
      TotalLateralRootLength[k]<-TotalRootLength[k]-FirstOrderRootLength[k]
      
      # Number, length and growth rate of lateral roots (by branching order)		
      
      if (maxord>1){
      
      for (l in 1:(maxord-1)){
        if (l<=max(DATA[[i]]$Ord)-1) {latroot[k,l]<-sum(data.split[[l+1]]!=0)} else {latroot[k,l]<-0}
        if (l<=max(DATA[[i]]$Ord)-1) {latroot[k,(l+(maxord-1))]<-sum(data.split[[l+1]])*cunit} else {latroot[k,(l+(maxord-1))]<-0}
        if (latroot[k,l]==0){latroot[k,(l+2*(maxord-1))]<-0} else {latroot[k,(l+2*(maxord-1))]<-latroot[k,(l+(maxord-1))]/latroot[k,l]}
        if (t==1) {latroot[k,(l+3*(maxord-1))]<-latroot[k,(l+(maxord-1))]/age[[i]][t]} else {latroot[k,(l+3*(maxord-1))]<-(latroot[k,(l+(maxord-1))]-latroot[k-1,(l+(maxord-1))])/(age[[i]][t]-age[[i]][t-1])}}}
      
      # Density of secondary roots on the first-order root
      
      if (FirstOrderRootLength[k]==0|maxord==1) {LateralRootDensity[k]<-0} else {LateralRootDensity[k]<-latroot[k,1]/FirstOrderRootLength[k]} 
      
      # Rootdiv
      
      if ((length(rootdiv)==1 & rootdiv[1]!=1)|(length(rootdiv)>1)) {
        
        Ord=NULL
        DBase=NULL
        
        # Secondary root number, length and density distribution on the first-order root
        
        if (length(rootdiv)==1){
          rootdiv1<-(0:rootdiv)*LastFirstOrderRootLength*cunit/rootdiv}
        else {rootdiv1<-rootdiv}
        
        for (j in 1:(length(rootdiv1)-1)){
          
         sub<-subset(DATA[[i]], (DBase!=0) & (DBase >= rootdiv1[j]/cunit) & (DBase < rootdiv1[j+1]/cunit) & (Ord==2))
         LRNDroot[k,j]<-sum(sub[[paste("Lengths",t,sep="")]]!=0)
         LRLDroot[k,j]<-sum(sub[[paste("Lengths",t,sep="")]])*cunit
         LRDDroot[k,j]<-LRNDroot[k,j]/(rootdiv1[j+1]-rootdiv1[j])}}}}	
  
  # Summarying results in a data frame
  
  if ((length(rootdiv)==1 & rootdiv[1]!=1)|(length(rootdiv)>1)) {
    
    if (maxord>1) {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, TotalLateralRootNumber, TotalLateralRootLength, latroot, LateralRootDensity, LRNDroot, LRLDroot, LRDDroot)}
    else {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, TotalLateralRootNumber, TotalLateralRootLength, LateralRootDensity, LRNDroot, LRLDroot, LRDDroot)}}
  
  else {
    if (maxord>1) {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, TotalLateralRootNumber, TotalLateralRootLength, latroot, LateralRootDensity)}
    else {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, TotalLateralRootNumber, TotalLateralRootLength, LateralRootDensity)}}
  
  if (maxord>1){
  LRnumberheading<-c()
  LRlengthheading<-c()
  LRmeanlengthheading<-c()
  LRgrowthrateheading<-c()
  
  for (h in 2:maxord){
    LRnumberheading[h-1]<-paste("N", h, "LR", sep="")
    LRlengthheading[h-1]<-paste("L", h, "LR", sep="")
    LRmeanlengthheading[h-1]<-paste("ML", h, "LR", sep="")
    LRgrowthrateheading[h-1]<-paste("GR", h, "L", sep="")}}
  
  if ((length(rootdiv)==1 & rootdiv[1]!=1)|(length(rootdiv)>1)) {
    
    LRNDheading<-c()
    if (length(rootdiv)==1){
      for (h in 1:rootdiv)
      {LRNDheading[h]<-paste("N2LR.Layer", h, sep="")}}
    else {
      for (h in 2:length(rootdiv)-1)
      {LRNDheading[h]<-paste("N2LR.", rootdiv[h], "to", rootdiv[h+1], sep="")}}
    
    LRLDheading<-c()
    if (length(rootdiv)==1){
      for (h in 1:rootdiv)
      {LRLDheading[h]<-paste("L2LR.Layer", h, sep="")}}
    else {
      for (h in 2:length(rootdiv)-1)
      {LRLDheading[h]<-paste("L2LR.", rootdiv[h], "to", rootdiv[h+1], sep="")}}
    
    LRDDheading<-c()
    if (length(rootdiv)==1){
      for (h in 1:rootdiv)
      {LRDDheading[h]<-paste("D2LR.Layer", h, sep="")}}
    else {
      for (h in 2:length(rootdiv)-1)
      {LRDDheading[h]<-paste("D2LR.", rootdiv[h], "to", rootdiv[h+1], sep="")}}}
  
  if ((length(rootdiv)==1 & rootdiv[1]!=1)|(length(rootdiv)>1)) {
    
    if (maxord>1) {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TNLR", "TLRL", t(LRnumberheading), t(LRlengthheading), t(LRmeanlengthheading), t(LRgrowthrateheading), "D2LR", t(LRNDheading), t(LRLDheading), t(LRDDheading))}
    else {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TNLR", "TLRL", "D2LR", t(LRNDheading), t(LRLDheading), t(LRDDheading))}}
  
  else {
    if (maxord>1) {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TNLR", "TLRL", t(LRnumberheading), t(LRlengthheading), t(LRmeanlengthheading), t(LRgrowthrateheading), "D2LR")}
    else {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TNLR", "TLRL", "D2LR")}}
  
  return(outputresults)}