latdist<-function(inputrac=NULL, inputrsml=NULL, output=c("lrd","dtp"), res=NULL, unitlength="px", int.length=NULL, interpol=NULL, rsml.date=NULL, rsml.connect=FALSE){
  
  # Errors interception
  
  if (is.null(inputrac)==TRUE & is.null(inputrsml)==TRUE){stop("inputrac and/or inputrsml must be provided")}
  
  if (is.null(inputrac)==FALSE) {if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}}
  
  if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
  
  if (output=="lrd"|output=="dtp"){} else {stop("The character string in output must be lrd or dtp")}
  
  if (is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0){stop("res must be a positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(interpol)==FALSE) {
    if (mode(interpol)!="numeric"){stop("mode(interpol) must be numeric")}
    if (interpol<=0){stop("interpol must be a positive value")}}
  
  if (is.null(int.length)==FALSE) {
    if (mode(int.length)!="numeric"){stop("mode(int.length) must be numeric")}
    if (int.length<=0){stop("int.length must be a positive value")}}
  
  if (is.null(inputrsml)==FALSE & is.null(rsml.date)==TRUE) {stop("If inputrsml is not NULL, rsml.date must be a positive numeric value")}
  
  if (is.null(rsml.date)==FALSE){
    if (rsml.date<=0|length(rsml.date)>1){stop("rsml.date must be a single positive value")}}
  
  if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
  
  # Reading of DART and rsml files
  
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
    if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}}
  else {
    if (is.null(inputrac)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
    else{
      if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
      if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}
  
  if (is.null(inputrsml)==TRUE){ # Only DART files
    
    DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
    
    for (i in 1:length(DATA)) {
      colnames(DATA[[i]])<-c()
      colnames(DATA[[i]])[1]<-"Root"
      colnames(DATA[[i]])[2]<-"Mother"
      colnames(DATA[[i]])[3]<-"Ord"
      colnames(DATA[[i]])[4]<-"DBase"
      colnames(DATA[[i]])[5]<-"DApp"
      for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}}
  
  else {
    
    if (is.null(inputrac)==TRUE){ # Only rsml files
      
      DATA<-list()
      filenamesrac<-c()
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect)
      for (i in 1:length(RSML)){
        for (j in 1:length(RSML[[i]]$rac)){colnames(RSML[[i]]$rac[[j]])[6]<-"Lengths1"}
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}}
    
    else { # DART and rsml files
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
      
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect)
      
      for (i in 1:length(RSML)){
        for (j in 1:length(RSML[[i]]$rac)){colnames(RSML[[i]]$rac[[j]])[6]<-"Lengths1"}
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}}}
 
  # Unit conversion
  
  if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
  if (unitlength=="cm") {cunit<-(cm(1)/res)}
  if (unitlength=="px") {cunit<-1}
  for (i in 1:length(DATA)){
    DATA[[i]]$DBase<-DATA[[i]]$DBase*cunit
    for (j in 6:ncol(DATA[[i]])) {DATA[[i]][,j]<-DATA[[i]][,j]*cunit}}
    
  # Calculating lateral root distributions

  rac<-list()
  res<-list()
  
  if (output=="lrd"){
  
    if (is.null(interpol)==TRUE) {
      
      # Calculating lateral root length and density distributions alons each mother root without interpolation
    
        for (i in 1:length(DATA)){
          
          if (nrow(DATA[[i]])>1){
          
          gen<-list()
          LRdensity<-list()
          latrootnum<-c()
          
          for (j in 1:nrow(DATA[[i]])){
            Root<-NULL
            DBase<-NULL
            sub<-subset(DATA[[i]], Mother==j-1, select=c(Root, DBase))
            latrootnum[j]<-nrow(sub)
            if (nrow(sub)==0) {gen[[j]]<-NULL} else {gen[[j]]<-sub}}
          
          finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
          global.latrootdensity<-latrootnum/finalrootlength
          
          rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
          
          for (j in 1:length(gen)){
            
            if (is.null(gen[[j]])==FALSE){
              
              lrdensity<-c()
              dbase<-c()
              lrlength<-c()
              t<-0
            
              for (k in 1:nrow(gen[[j]])){
                
                if (gen[[j]]$DBase[k]>=(int.length/2) & gen[[j]]$DBase[k]<=finalrootlength[j]-(int.length/2)){
                  
                  t<-t+1
                
                  min<-gen[[j]]$DBase[k]-(int.length/2)
                  max<-gen[[j]]$DBase[k]+(int.length/2)
                  
                  dbase[t]<-gen[[j]]$DBase[k]
                  lrdensity[t]<-sum(gen[[j]]$DBase>=min & gen[[j]]$DBase<=max)/int.length
                  lrlength[t]<-sum(finalrootlength[1+gen[[j]]$Root[gen[[j]]$DBase>=min & gen[[j]]$DBase<=max]])/int.length}}
            
              if (length(dbase)!=0) {LRdensity[[j]]<-data.frame(DBase=dbase, LRD=lrdensity, LRL=lrlength)} else {LRdensity[[j]]<-NULL}}
          
            else {
              
              LRdensity[[j]]<-NULL}}
          
          res[[i]]<-LRdensity}
        
        else {
          
          latrootnum<-0
          finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
          global.latrootdensity<-latrootnum/finalrootlength
          
          rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
          res[[i]]<-NULL}}
        
          names(rac)<-filenamesrac
          names(res)<-filenamesrac
          outputresults<-list(root=rac, res=res)}
  
    else {
      
      # Calculating lateral root length and density distributions alons each mother root with interpolation
    
      for (i in 1:length(DATA)){
        
        if (nrow(DATA[[i]])>1){
        
        gen<-list()
        LRdensity<-list()
        latrootnum<-c()
        
        for (j in 1:nrow(DATA[[i]])){
          Mother<-NULL
          DBase<-NULL
          Root<-NULL
          sub<-subset(DATA[[i]], Mother==j-1, select=c(Root, DBase))
          latrootnum[j]<-nrow(sub)
          if (nrow(sub)==0) {gen[[j]]<-NULL} else {gen[[j]]<-sub}}
        
        finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
        global.latrootdensity<-latrootnum/finalrootlength
        
        rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
        
        for (j in 1:length(gen)){
          
          if (is.null(gen[[j]])==FALSE){
            
            lrdensity<-c()
            dbase<-c()
            lrlength<-c()
            seqDBase<-seq(from=0, to=finalrootlength[j], by=finalrootlength[j]/(interpol-1))
            t<-0
            
            for (k in 1:length(seqDBase)){
              
              if (seqDBase[k]>=(int.length/2) & seqDBase[k]<=finalrootlength[j]-(int.length/2)){
                
                t<-t+1
                
                min<-seqDBase[k]-(int.length/2)
                max<-seqDBase[k]+(int.length/2)
                
                dbase[t]<-seqDBase[k]
                lrdensity[t]<-sum(gen[[j]]$DBase>=min & gen[[j]]$DBase<=max)/int.length
                lrlength[t]<-sum(finalrootlength[1+gen[[j]]$Root[gen[[j]]$DBase>=min & gen[[j]]$DBase<=max]])/int.length}}
            
            if (length(dbase)!=0) {LRdensity[[j]]<-data.frame(DBase=dbase, LRD=lrdensity, LRL=lrlength)} else {LRdensity[[j]]<-NULL}}
          
          else {
            
            LRdensity[[j]]<-NULL}}
        
        res[[i]]<-LRdensity}
        
        else{
          
          latrootnum<-0
          finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
          global.latrootdensity<-latrootnum/finalrootlength
          
          rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
          res[[i]]<-NULL}}
      
      names(rac)<-filenamesrac
      names(res)<-filenamesrac
      outputresults<-list(root=rac, res=res)}}

  
  if (output=="dtp"){
    
    # Calculating the distance between neighbouring roots along each mother root
    
    for (i in 1:length(DATA)){
      
      if (nrow(DATA[[i]])>1){
      
      gen<-list()
      dtp<-list()
      latrootnum<-c()
      
      for (j in 1:nrow(DATA[[i]])){
        DBase<-NULL
        sub<-subset(DATA[[i]], Mother==j-1, select=DBase)
        latrootnum[j]<-nrow(sub)
        if (nrow(sub)<2) {gen[[j]]<-NULL} else {
          gen[[j]]<-sub[order(sub$DBase),]}}
      
      finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
      global.latrootdensity<-latrootnum/finalrootlength
      
      rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
      
      if (length(gen)>0){
      
        for (j in 1:length(gen)){
        
          if (is.null(gen[[j]])==FALSE) {
        
            difftoprec<-diff(gen[[j]])
            dbase<-gen[[j]][2:length(gen[[j]])]
          
            dtp[[j]]<-data.frame(DBase=dbase, DTP=difftoprec)}
        
          else {
          
            dtp[[j]]<-NULL}}
      
        res[[i]]<-dtp}
      
      if (length(gen)==0) {res[[i]]<-NULL}}
      
      else {
        
        latrootnum<-0
        finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
        global.latrootdensity<-latrootnum/finalrootlength
        
        rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
        res[[i]]<-NULL}}
    
    names(rac)<-filenamesrac
    names(res)<-filenamesrac
    outputresults<-list(root=rac, res=res)}

return(outputresults)}