#==========================================================#
# Name: AMGET.R
# Purpose: Post processing tool for ADAPT 5
# Author: Benjamin. G
# Last Edited: Aug 1st 2013
#==========================================================#

AMGET <- function() {
  
#-----------------------Create environment
ENV  <- new.env()


#-----------------------SetWD
SetWD <- function() {
  setwd(paste(.libPaths(),'AMGET',sep='/'))
  if(file.exists('MyDir.Rdata')==F) {
    if(.Platform$OS.type == "unix") {    # On Mac
      MyDir   <- readline(cat('\n\nIt appears to be the first time you are running AMGET. The settings spreadsheet will be created, please indicate the destination folder below:'))
      }else{                             # On PC
      cat('\nIt appears to be the first time you are running AMGET.\nThe settings spreadsheet will be created, please indicate the destination folder in the window that just promped\n')
      MyDir   <- choose.dir(caption='Indicate the settings destination folder')
      }
    
    if(file.exists(MyDir)==F) {SetWD()
                               }else{
     save(MyDir,file='MyDir.Rdata')
     FILES   <- c('AMGET_Settings.csv','AMGET_User_Manual.pdf')
     file.copy(FILES,paste(MyDir,FILES,sep='/'))
     cat(paste('\nThe settings spreadsheet and the user manual were successfully created in:\n',MyDir,sep=''))}}
  
  load('MyDir.Rdata') 
  
  if(file.exists(paste(MyDir,'AMGET_Settings.csv',sep='/'))==F) {
    file.remove('MyDir.Rdata')
    SetWD()
  }else{
    setwd(MyDir)
    settings   <- read.csv('AMGET_Settings.csv', as.is=T)  
    MacUSER    <- toupper(settings[grep('#DIR',settings[,1])+2,2])
    OUT.DIR    <- settings[grep('#DIR',settings[,1])+1,2]
    if(OUT.DIR=='') {
      
      if(.Platform$OS.type == "unix") {    # On Mac
        OUT.DIR  <- readline(cat('\n\nDefine the folder where all the graphics will be saved (this can be redefined in the settings file under output directory later on for enhanced organization of additional projects/runs):'))
        }else{                             # On PC
          OUT.DIR   <- choose.dir(caption='Indicate the plots destination folder')    
          cat('\n\nEnter the path to the folder where all the graphics will be saved (this can be redefined in the settings file under output directory later on for enhanced organization of additional projects/runs):')
        }
            
      settings[grep('#DIR',settings[,1])+1,2] <- OUT.DIR
      write.table(settings,file='AMGET_Settings.csv',row.names=F,quote=F,sep=',')
      cat('\n\n** The set up was successful, you can now start using AMGET **\n\n')}
    
    cat('\n\nWelcome to AMGET')    
    if(grepl('N',MacUSER)){
       DIR   <- readline(cat('\n\nPlease indicate the path to the folder containing the model(s) files you want to evaluate:')) 
    }else{
       cat('\n(a window just promped to allow you to define your WD)')
       DIR   <- choose.dir(default=getwd(),caption='Select the folder containing your model')}
  
    if(toupper(DIR)=='EXIT') {stop('Au revoir!',call. = FALSE)
        }else{
          if(file.exists(DIR)==F) {
             cat('\n\nSorry but this directory does not appear to exist...\n\n')
             SetWD()
             }else{
      
            assign('WORK.DIR', DIR, envir = ENV)
            assign('MyDir', MyDir,  envir = ENV)
            Main.menu()}}}
  }



# -----------------------Main Menu
Main.menu <- function() {
  setwd(ENV$WORK.DIR)
  switch(menu(c('Indiv. DV vs. IV profiles',
                'Goodness of fit plots',
                'Posthoc fit plots',
                'Parameter distribution profiles',
                'Visual Predictive Check',
                'Change WD',
                'Quit'),
         title= cat(paste('\n\nAMGET Main Menu\n\nWD:',
                          ENV$WORK.DIR,'\n\nSelect one of the following options:',sep=''))) +1,
         cat('Au revoir!'), 
         IDV.menu(),          # Option 1
         GOF(),               # Option 2
         PHF(),               # Option 3
         PRM(),               # Option 4
         VPC(),               # Option 5
         SetWD(),             # Option 6
         cat('Au revoir!'))}  # Option 7



# -----------------------IDV Menu
IDV.menu <- function() {
  switch(menu(c('NONMEM-shaped dataset','ADAPT population dataset','ADAPT individualized datasets',
                'Main Menu','Quit'), 
              title= cat(paste('\n\nIndividual DV vs. IV\n',
                               '\nSelect one of the following options:', sep='')))+1,
         cat('Au revoir!'), 
         IDV.NONMEM(),        # Option 1
         IDV('POP'),          # Option 2
         IDV.IND(),           # Option 3
         Main.menu(),         # Option 4
         cat('Au revoir!'))}  # Option 5



# -----------------------IDV GRAPHS  
IDV.Graphs <- function(DIR,table,dat.obs,dat.bolus,dat.inf,CMT.name,IV.coln,ID.coln,DV.coln,AMT,
                       CMTs,C,IDs,Totsub,LOG,LAXIS,LABEL,SOFT){
  # Get the setting options
  settings   <- read.csv(paste(DIR,'AMGET_Settings.csv',sep='/'), as.is=T)
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  DV.name    <- settings[grep('#IDV',settings[,1])+1,2]
  DV.units   <- paste('(',settings[grep('#IDV',settings[,1])+2,2],')',sep='')
  IV.name    <- settings[grep('#IDV',settings[,1])+3,2]
  IV.units   <- paste('(',settings[grep('#IDV',settings[,1])+4,2],')',sep='')
  Dose.units <- settings[grep('#IDV',settings[,1])+5,2]
  Mcol       <- c(col2rgb(settings[grep('#IDV',settings[,1])+6,2]),settings[grep('#IDV',settings[,1])+6,2])
  Mch        <- as.numeric(settings[grep('#IDV',settings[,1])+7,2])
  Alpha      <- as.numeric(settings[grep('#IDV',settings[,1])+8,2])
  Lty        <- as.numeric(settings[grep('#IDV',settings[,1])+9,2])
  Lwd        <- as.numeric(settings[grep('#IDV',settings[,1])+10,2])
  Dcol       <- settings[grep('#IDV',settings[,1])+11,2]
  BQL        <- as.numeric(settings[grep('#IDV',settings[,1])+12,2])
  BQLcol     <- settings[grep('#IDV',settings[,1])+13,2]
  BQLty      <- as.numeric(settings[grep('#IDV',settings[,1])+14,2])
  if(length(IDs)>1){Layout <- as.numeric(settings[grep('#IDV',settings[,1])+15,2])
  }else{Layout  <- 1}
  Scatter    <- toupper(settings[grep('#IDV',settings[,1])+16,2])
  PosLegend  <- settings[grep('#IDV',settings[,1])+17,2]
  Title      <- toupper(settings[grep('#IDV',settings[,1])+18,2])
     
  Population  <- function(Name,Type,PCH,LT,LW){   # Spaghetti or Scatter plots
    par(mfrow=c(1,1),cex=1,las=1)
    plot.new()
    mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
    mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
    if (LOG=='LOG') {mtext(paste(Name,'plots (Log Scale)'),line=-16, side=3, cex=2)
    }else{mtext(paste(Name,'plots'),line=-16, side=3, cex=2)}  
    mtext(paste('Individuals', DV.name,'vs.',IV.name,sep=' '),line=-18.2, side=3, cex=1.4)
    mtext(paste('Y(',CMTs[C],'): ',CMT.name[C],sep=''),line=-20.1, side=3, cex=1.2)
    mtext(paste(Totsub,'subjects |',nrow(dat.obs),'observations of which',
                length(dat.obs[dat.obs[,DV.coln]<BQL ,DV.coln]),'BQL',sep=' '),line=-21.7, 
          side=3, cex=1) 
    
    if(nrow(dat.bolus)!=0){
      dat.bolus <- as.data.frame(dat.bolus)
      dat.dose  <- dat.bolus[!duplicated(dat.bolus[,ID.coln]),]
    }else{dat.inf   <- as.data.frame(dat.inf)
          dat.dose  <- dat.inf[!duplicated(dat.inf[,ID.coln]),]} 
    
    AMTs  <- sort(unique(dat.dose[,AMT]))
    
    for (d in 1:length(AMTs)) { 
      dat.AMT <- dat.obs[dat.obs[,ID.coln]%in%unique(dat.dose[dat.dose[,AMT]==AMTs[d],ID.coln]),]
      nobs    <- length(dat.AMT[,DV.coln])
      nID     <- unique(dat.AMT[,ID.coln])
      nsub    <- length(nID)
      nBQL    <- length(dat.AMT[dat.AMT[,DV.coln]<BQL ,DV.coln])
      ivlim   <- as.numeric(c(0,max(dat.AMT[,IV.coln])))    
      if (LOG=='LOG') {limind  <- c(min(dat.AMT[,DV.coln]),max(dat.AMT[,DV.coln]))
      }else{limind <- c(0,max(dat.AMT[,DV.coln]))}
      if (Title=='Y') {Main <- paste('Dose:',AMTs[d],Dose.units,'|',nsub,'subj. |', nobs, 
                                     'obs. of which',nBQL,'BQL |',CMT.name[C],LABEL, sep=' ')
      }else{Main <- NULL}
      if (LOG=='LOG') {
        plot(x=1,y=1,type='n',xlab=paste(IV.name,IV.units,sep=' '),ylab=paste(DV.name,DV.units,sep=' '), 
             xlim=ivlim,ylim=limind,log='y',yaxt='n',main=Main)
        axis(2,as.vector(1:10 %o% 10^(-5:5)),labels=F,tcl = -0.25)
        axis(2,as.vector(10*10^(-5:5)),labels=as.character(as.vector(10*10^(-5:5))))
      }else{
        plot(x=1,y=1,type='n',xlab=paste(IV.name,IV.units,sep=' '),ylab=paste(DV.name,DV.units,sep=' '), 
             xlim=ivlim,ylim=limind,main=Main)}
      grid(col='grey90',lty=1,lwd=0.01,equilogs=(LOG!='LOG'))
      abline(h=as.numeric(BQL),col=BQLcol,lty=BQLty)
      for(s in 1:nsub){
        lines(x=dat.AMT[dat.AMT[,ID.coln]==nID[s],IV.coln],y=dat.AMT[dat.AMT[,ID.coln]==nID[s],DV.coln],
              type=Type,pch=PCH,lwd=LW,col=rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),
              lty=LT,cex=1.2)}
      if(Title=='Y'){mtext(paste(table,'| Individuals', DV.name,'vs.',IV.name,'profiles',sep=' '),
                           cex=0.9,line=0.3, side=3)
      }else{mtext(paste('Dose:',AMTs[d],Dose.units),cex=0.8,side=1,line=-1.2)}
      text(x=max(dat.AMT[,IV.coln]),y=as.numeric(BQL),pos=1,cex=0.8,offset=0.4,
           labels=paste(BQL,substring(DV.units,2,nchar(DV.units)-1),sep=' '))
      legend (PosLegend,legend=c('Observations','Quantification Lim.'),lty=c(LT,BQLty),
              col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),BQLcol),pch=c(PCH,NA),
              lwd=c(LW,1),bty='n',cex=0.9)
    }} # End dose levels loop
  
  #Population graphs
  if(Scatter!='NONE' & length(IDs)>1){
  switch(Scatter,
         SPAGHETTI = Population('Spaghetti','l',NA,Lty,Lwd),
         SCATTER   = Population('Scatter','p',Mch,0,0),
         BOTH      = {Population('Spaghetti','l',NA,Lty,Lwd);Population('Scatter','p',Mch,0,0)})}

  # Individual profiles
  par(mfrow=c(1,1),cex=1,las=1)
  plot.new()
  mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
  mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
  if (LOG=='LOG') {mtext('Individual profiles (Log Scale)',line=-16, side=3, cex=2)
  }else{mtext('Individual profiles',line=-16, side=3, cex=2)}
  mtext(paste('Individuals', DV.name,'vs.',IV.name,sep=' '),line=-18.2, side=3, cex=1.4)
  mtext(paste('Y(',CMTs[C],'): ',CMT.name[C],sep=''),line=-20.1, side=3, cex=1.2)
  mtext(paste(Totsub,'sub. |',nrow(dat.obs),'obs.of which',
              length(dat.obs[dat.obs[,DV.coln]<BQL ,DV.coln]),'BQL',sep=' '),line=-21.7, 
              side=3, cex=1)
  
  switch(Layout,'1'=par(las=1,mfrow=c(1,1),mar=c(5,4,4,2),cex=0.9),
         '2'=par(las=1,mfrow=c(2,2),mar=c(3.6,4,4,0.6),cex=0.8),
         '3'=par(las=1,mfrow=c(3,3),mar=c(3.6,4,2,0.6),cex=0.55))
  for(a in 1:Totsub) {
      nobs   <- length(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])
      if(nobs!=0){
      nBQL   <- length(dat.obs[dat.obs[,ID.coln]==IDs[a] & dat.obs[,DV.coln]<BQL ,DV.coln])
      ivlim  <- c(0,max(dat.obs[dat.obs[,ID.coln]==IDs[a],IV.coln]))
      if (LOG=='LOG') {limind <- c(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln]),
                                   max(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln]))
      }else{limind <- c(0,max(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln]))}
      if(SOFT=='NONMEM' & Title=='Y') {Main <- paste('Sub.',IDs[a],'|', nobs, 'obs. of which',nBQL,
                                                     'BQL |',CMT.name[C],LABEL, sep=' ')}
      if(SOFT!='NONMEM' & Title=='Y') {Main <- paste(ifelse(substr(IDs[a],nchar(IDs[a]),
                                                     nchar(IDs[a]))==',',unlist(strsplit(IDs[a],',',
                                                     fixed=T)),IDs[a]) ,'|', nobs, 'obs. of which',
                                                     nBQL,'BQL |',CMT.name[C],LABEL, sep=' ')}
      if (LOG=='LOG') {
      plot(x=1,y=1,type='n',xlab=paste(IV.name,IV.units,'\n',sep=' '), 
           ylab=paste(DV.name,DV.units,sep=' '),yaxt='n',log='y',xlim=ivlim,ylim=limind,main=Main)
      axis(2,as.vector(1:10 %o% 10^(-5:5)),labels=F,tcl = -0.25)
      axis(2,as.vector(10*10^(-5:5)),labels=as.character(as.vector(10*10^(-5:5))))
      }else{
      plot(x=1,y=1,type='n',xlab=paste(IV.name,IV.units,'\n',sep=' '),
           ylab=paste(DV.name,DV.units,sep=' '),xlim=ivlim,ylim=limind,main=Main)}     
      grid(col='grey90',lty=1,lwd=0.01,equilogs=(LOG!='LOG'))  
      if(Layout!=3){ if(Title=='Y') {mtext(paste(table,'| Individuals', DV.name,'vs.',IV.name,
                                           'profiles',sep=' '),line=0.3, side=3,cex=0.8)
      }else{ if(SOFT=='NONMEM'){mtext(paste('ID',IDs[a]),cex=0.6,side=1, line=-1.2)
                          }else{mtext(IDs[a],cex=0.8,side=1, line=-1.2)}}}
      lines(x=dat.obs[dat.obs[,ID.coln]==IDs[a],IV.coln],y=dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln],
            col=rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),type='o',lty=Lty,lwd=Lwd,pch=Mch,cex=1.2)
      
      # Add Bolus info
      if(nrow(dat.bolus)>0){
        if(IDs[a] %in% dat.bolus[,ID.coln]){
          abline(v=dat.bolus[dat.bolus[,ID.coln]==IDs[a] & dat.bolus[,AMT]>0,IV.coln], col=Dcol,lty=2)
          text(x=dat.bolus[dat.bolus[,ID.coln]==IDs[a] & dat.bolus[,AMT]>0,IV.coln],
               y=limind[2]-limind[2]*0.05,pos=2,col=Dcol,srt=90,
               labels=paste('B',dat.bolus[dat.bolus[,ID.coln]==IDs[a] & dat.bolus[,AMT]>0,'CMT'],':',
                            dat.bolus[dat.bolus[,ID.coln]==IDs[a] & dat.bolus[,AMT]>0,AMT],
                            Dose.units, sep=''))}}
      
      # Add Inf info for ADAPT
      if(nrow(dat.inf)>0 & SOFT=='ADAPT'){
        if(IDs[a] %in% dat.inf[,ID.coln]){
          for(I in 1:(length(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln])/2)) {
            arrow.line <- seq(0,(length(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln])/2),by=2)
            arrows(x0 =as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln][arrow.line[I]+1]),
                   y0 =(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])+(limind[2]-limind[1])*0.04*I),
                   x1 =as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln][arrow.line[I]+2]),
                   y1 =(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])+(limind[2]-limind[1])*0.04*I),
                   col=Dcol,length=0.1,lwd=1)
            text(x=((as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln][arrow.line[I]+1])+
              as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln][arrow.line[I]+2]))/2),
              y=(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])+(limind[2]-limind[1])*0.04*I),
              pos=1,col=Dcol, offset = 0.2,
              labels=paste('R',dat.inf[dat.inf[,ID.coln]==IDs[a],'CMT'][arrow.line[I]+1],':',
                           dat.inf[dat.inf[,ID.coln]==IDs[a],AMT][arrow.line[I]+1],Dose.units,'/',
                           (as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln][arrow.line[I]+2])-
                           as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln][arrow.line[I]+1])),
                           substring(IV.units,2,nchar(IV.units)-1),sep=''))}}}
      
      # Add Inf info for NONMEM
      if(nrow(dat.inf)>0 & SOFT=='NONMEM'){
        if(IDs[a] %in% dat.inf[,ID.coln]){
          for(I in 1:length(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln])) { # Nb of inf per subject
            start_inf <- as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],IV.coln])
            end_inf   <- start_inf+(as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],AMT])/
                         as.numeric(dat.inf[dat.inf[,ID.coln]==IDs[a],'RATE']))
            arrows(x0 =start_inf,
                   y0 =(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])+(limind[2]-limind[1])*0.04*I),
                   x1 =end_inf,
                   y1 =(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])+(limind[2]-limind[1])*0.04*I),
                   col=Dcol,length=0.1,lwd=1)
            text(x=(end_inf+start_inf)/2,
                 y=(min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])+(limind[2]-limind[1])*0.04*I),
                 pos=1,col=Dcol, offset = 0.2,
                 labels=paste('R',dat.inf[dat.inf[,ID.coln]==IDs[a],'CMT'],':',
                              dat.inf[dat.inf[,ID.coln]==IDs[a],AMT],Dose.units,'/',
                              (end_inf-start_inf),substring(IV.units,2,nchar(IV.units)-1),sep=''))}}}
      
      # Add Cmax and Cmin info
      Cmax <- max(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])
      Tmax <- dat.obs[dat.obs[,ID.coln]==IDs[a] & dat.obs[,DV.coln]==Cmax,IV.coln]
      Cmin <- min(dat.obs[dat.obs[,ID.coln]==IDs[a],DV.coln])
      Tmin <- dat.obs[dat.obs[,ID.coln]==IDs[a] & dat.obs[,DV.coln]==Cmin,IV.coln]
      text(x=Tmax,y=Cmax,pos=4,col=Mcol[4],offset=0.3,
           labels=paste(signif(Cmax,digits=4),substring(DV.units,2,nchar(DV.units)-1),' at ',
                        signif(Tmax,digits=3),substring(IV.units,2,nchar(IV.units)-1),sep=''))
      text(x=Tmin,y=Cmin,pos=4,col=Mcol[4],offset=0.3,
           labels=paste(signif(Cmin,digits=4),substring(DV.units,2,nchar(DV.units)-1),' at ',
                        signif(Tmin,digits=3),substring(IV.units,2,nchar(IV.units)-1),sep=''))      
      abline(h=as.numeric(BQL), col=BQLcol,lty=BQLty)
      text(x=max(dat.obs[dat.obs[,ID.coln]==IDs[a],IV.coln])*0.97,y=as.numeric(BQL),pos=1,
           offset = 0.4, labels=paste(BQL,substring(DV.units,2,nchar(DV.units)-1),sep=' '))
      legend(PosLegend,legend=c('Observations','Dosing Times','Quantification Lim.'),
             lty=c(Lty,2,BQLty),pch=c(Mch,NA,NA),col=c(rgb(Mcol[1],Mcol[2],Mcol[3],
             alpha=255*Alpha,maxColorValue=255),Dcol,BQLcol),bty='n',cex=0.9)
      }else{cat(paste('\nSkipped subject',IDs[a],': no observations in compartment',CMTs[C]))}
      } # End of loop for each subject
  dev.off()}



# -----------------------IDV Nonmem
IDV.NONMEM <- function() {
  # Import the data file  
  setwd(ENV$WORK.DIR)
  List.csv    <- list.files(pattern='.csv$')
  
  if(length(List.csv)>1) {
    table     <-  select.list(List.csv,title=cat('\nSelect the NONMEM-shaped dataset you want to evaluate:'))
  }else{table <- List.csv[as.numeric(1)]}
  
  dat         <- read.csv(table,header=T,as.is=T)
  
  # Get the setting options
  settings   <- read.csv(paste(ENV$MyDir,'AMGET_Settings.csv',sep='/'), as.is=T)
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  DV.name    <- settings[grep('#IDV',settings[,1])+1,2]
  IV.name    <- settings[grep('#IDV',settings[,1])+3,2]
  IV.coln    <- toupper(IV.name)
  ID.coln    <- c('ID')
  DV.coln    <- c('DV')
  AMT        <- c('AMT')
  LOG_PLOT   <- toupper(settings[grep('#IDV',settings[,1])+19,2])
  
  # Check the column names
  if(!'CMT'%in%colnames(dat)){
    dat$CMT<-1
    CMT.name<-c('Default')
  }else{CMT.name<-settings[grep('#IDV',settings[,1])+20,2:(length(unique(dat[dat$MDV==0,'CMT']))+1)]} 
  if(!'ID'%in%colnames(dat)){
    ID.coln <- select.list(colnames(dat),title=cat('\nThe name ID could not be found in the column names,\nselect a substitution column:'))
    if(!ID.coln%in%colnames(dat)){stop('Au revoir!',call.=FALSE)}}
  if(!'DV'%in%colnames(dat)){
    DV.coln <- select.list(colnames(dat),title=cat('\nThe name DV could not be found in the column names,\nselect a substitution column:'))
    if(!DV.coln%in%colnames(dat)){stop('Au revoir!',call.=FALSE)}}
  if(!IV.coln%in%colnames(dat)){
    IV.coln <- select.list(colnames(dat),title=cat(paste('\nThe name ',IV.coln,' could not be found in the column names,\nselect a substitution column:',sep='')))
    if(!IV.coln%in%colnames(dat)){stop('Au revoir!',call.=FALSE)}}
  if(!AMT%in%colnames(dat)){
    AMT <- select.list(colnames(dat),title=cat('\nThe name AMT could not be found in the column names,\nselect a substitution column:'))
    if(!AMT%in%colnames(dat)){stop('Au revoir!',call. = FALSE)}}
  
  # Subset the data
  if('RATE' %in% colnames(dat)) {
    dat.inf     <- dat[dat$MDV==1 & dat[,AMT]>0 & dat$RATE>0,]
    dat.bolus   <- dat[dat$MDV==1 & dat[,AMT]>0 & (dat$RATE==0 | dat$RATE=='.'),]
  }else{
    dat.inf     <- dat[dat$MDV==999,]  # Create an empty file to avoid crash
    dat.bolus   <- dat[dat$MDV==1 & dat[,AMT]>0,]}
  
  dat           <- dat[dat$MDV==0,]
  dat[,DV.coln] <- as.numeric(dat[,DV.coln])
  dat[,IV.coln] <- as.numeric(dat[,IV.coln])
  CMTs          <- unique(dat$CMT)
  
  dir.create(path=paste(output.dir,'/',substring(table,1,nchar(table)-4),'_PLOTS',sep=''),
             showWarnings=F)
  setwd(paste(output.dir,'/',substring(table,1,nchar(table)-4),'_PLOTS', sep=''))
  cat('\nWorking...')
  
  # Plots
  for (C in 1:length(CMTs)) {  
    # Linear Scale
    dat.obs <- dat[dat$CMT==CMTs[C] & dat[,DV.coln]>=0,]
    dat.obs <- dat.obs[order(dat.obs[,ID.coln],dat.obs[,IV.coln]),]
    IDs     <- unique(dat.obs[,ID.coln])
    Totsub  <- length(IDs)
    pdf(file=paste(substring(table,1,nchar(table)-4),'_',DV.name,'_vs_',IV.name,'_Y',CMTs[C],
                   '.pdf',sep=''),width=11.7,height=8.3)
    IDV.Graphs(ENV$MyDir,table,dat.obs,dat.bolus,dat.inf,CMT.name,IV.coln,ID.coln,DV.coln,AMT,
               CMTs,C,IDs,Totsub,'LIN','','','NONMEM')
    
    # Log Scale
    dat.obs <- dat[dat$CMT==CMTs[C] & dat[,DV.coln]>0,]
    dat.obs <- dat.obs[order(dat.obs[,ID.coln],dat.obs[,IV.coln]),]
    IDs    <- unique(dat.obs[,ID.coln])
    Totsub <- length(IDs)
    pdf(file=paste(substring(table,1,nchar(table)-4),'_',DV.name,'_vs_',IV.name,'_Y',CMTs[C],
                   '_(Log_Scale).pdf',sep=''),width=11.7,height=8.3)
    IDV.Graphs(ENV$MyDir,table,dat.obs,dat.bolus,dat.inf,CMT.name,IV.coln,ID.coln,DV.coln,AMT,
               CMTs,C,IDs,Totsub,'LOG','y','| Log Scale','NONMEM')
  }
  
  switch(menu(c('Generate more individual profiles',
                'Go back to the individual DV vs. IV Menu',
                'Go back to the Main Menu','Quit'), 
              title= cat(paste('\n\nThe individuals ',DV.name,' vs. ',IV.name,
                               ' plots have successfully\nbeen created in the folder:\n',
                               getwd(),'\n\nTo continue select one of the following options:',
                               sep='')))+1,
         cat('Au revoir!'), 
         IDV.NONMEM(),        # Option 1
         IDV.menu(),          # Option 2
         Main.menu(),         # Option 3
         cat('Au revoir!'))}  # Option 4



# -----------------------IDV Ind and Pop
IDV.IND  <- function(){
  setwd(ENV$WORK.DIR)
  List.dat  <- list.files(pattern='.dat$')
  if(length(List.dat)==1) {IDV('IND')}else{
  switch(menu(c('Yes',
                'Go back to the individual DV vs. IV Menu',
                'Go back to the Main Menu','Quit'), 
              title=cat('\nThe',length(List.dat),'following datasets have been found and will be imported:',
                        paste('\n  ',List.dat[1:length(List.dat)],sep=''),
                        '\n\nDo you wish to continue?'))+1,
         stop('Au revoir!',call. = FALSE),
         IDV('IND'),                         # Option 1
         IDV.menu(),                         # Option 2
         Main.menu(),                        # Option 3
         stop('Au revoir!',call. = FALSE))}} # Option 4

IDV <- function(METHOD) {
  setwd(ENV$WORK.DIR)
  List.dat  <- list.files(pattern='.dat$')
  if(METHOD=='POP'){
    if(length(List.dat)>1) {
      table   <-  select.list(List.dat,title= cat('\nSelect the ADAPT population dataset you want to evaluate:'))
    }else{table <- List.dat[as.numeric(1)]}
    # Import dataset
    dat    <- read.table(table,stringsAsFactors=F,header=F,fill=T,col.names=1:100)
    dat    <- dat[,colSums(is.na(dat[,1:100]))!=nrow(dat)]  
    CMTs   <- 1:(ncol(dat)-1)
    if(is.numeric(dat[1,1])){
        METHOD <- 'IND'
        cat('\nWarning:\nThe selected file does not look like a population dataset, AMGET is trying using individual instead.')}}  
  if(METHOD=='IND'){
    # Import Ind dataset and compile pop dataset
    dat <- c()
    for(A in 1:length(List.dat)) {
      dat.temp  <- read.table(List.dat[A],stringsAsFactors=F,header=F,fill=T,col.names=1:100)
      dat.temp  <- rbind(c(substring(List.dat[A],1,nchar(List.dat[A])-4),rep(NA,99)),dat.temp)
      dat       <- rbind(dat,dat.temp)}
    dat         <- dat[,colSums(is.na(dat[,1:100]))!=nrow(dat)]
    dat.raw     <- dat # Save for export later
    CMTs        <- 1:(ncol(dat)-1)
    if(length(List.dat)==1){
    table <- List.dat[as.numeric(1)]
    }else{
    table       <- strsplit(getwd(),'/')
    table       <- table[[1]][length(table[[1]])]
    table       <- paste('Pooled datasets from .../',table,sep='')}}

  # Compile the data
  cat('\n\nWorking...')
  dat.bolus <- dat.obs <- dat.inf <- c()
  rows <- NumID <- Ninf <- Nbolus <- 0
  
  while((rows)<nrow(dat)) {
    IDs       <- dat[(1+rows),1]
    NumID     <- NumID+1
    Ndose     <- as.numeric(dat[(4+rows),1])
    Ninf      <- as.numeric(dat[(2+rows),1])
    Nbolus    <- as.numeric(dat[(3+rows),1])
    dos.Times <- dat[(1:Ndose)+(4+rows),1]
    
    # Infusion
    if(Ninf>0){AMT.inf <- dat[(1:Ndose)+(4+rows),(1:Ninf)+1]
               if(Ninf==1) {dim(AMT.inf) <- c(Ndose,Ninf)}
               for(R in 1:Ninf){
                 temp    <- cbind(rep(IDs,times=length(dos.Times)),
                                  rep(NumID,times=length(dos.Times)),dos.Times,AMT.inf[,R],R)
                 dat.inf <- rbind(dat.inf,temp)}
               colnames(dat.inf) <- c('ID','COUNT','TIME','AMT','CMT')}
    
    # Bolus
    if(Nbolus>0){AMT.bolus <- dat[(1:Ndose)+(4+rows),(1:Nbolus)+1+Ninf]
                 if(Nbolus==1) {dim(AMT.bolus) <- c(Ndose,Nbolus)}
                 for(B in 1:Nbolus){
                   temp      <- cbind(rep(IDs,times=length(dos.Times)),
                                      rep(NumID,times=length(dos.Times)),dos.Times,AMT.bolus[,B],B)
                   dat.bolus <- rbind(dat.bolus,temp)}
                 colnames(dat.bolus) <- c('ID','COUNT','TIME','AMT','CMT')}
    
    # Observation
    obs.Times <- dat[(1:dat[(6+Ndose+rows),1])+(6+Ndose+rows),1]
    obs       <- dat[(1:dat[(6+Ndose+rows),1])+(6+Ndose+rows),CMTs+1]
    if(dat[(6+Ndose+rows),1]==1){obs <- t(obs)} # If only one obs avoid it to be put in column
    temp.obs  <- data.frame(stringsAsFactors=F,
                            ID   = rep(IDs,times=length(obs.Times)),
                            COUNT= rep(NumID,times=length(obs.Times)),
                            TIME = obs.Times,
                            DV   = obs)
    colnames(temp.obs) <- c('ID','COUNT','TIME',paste('DV',CMTs,sep=''))
    dat.obs   <- rbind(dat.obs,temp.obs) 
    rows = (rows + 6 + length(obs.Times) + length(dos.Times))}
  
  dat       <- dat.obs  # To match the other scripts
  dat.inf   <- as.data.frame(dat.inf,stringsAsFactors=F)
  dat.bolus <- as.data.frame(dat.bolus,stringsAsFactors=F)
  
  # Get the setting options
  settings   <- read.csv(paste(ENV$MyDir,'AMGET_Settings.csv',sep='/'), as.is=T)
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  DV.name    <- settings[grep('#IDV',settings[,1])+1,2]
  IV.name    <- settings[grep('#IDV',settings[,1])+3,2]
  LOG_PLOT   <- toupper(settings[grep('#IDV',settings[,1])+19,2])
  CMT.name   <- settings[grep('#IDV',settings[,1])+20,CMTs+1]
  IV.coln    <- c('TIME')
  ID.coln    <- c('ID')
  AMT        <- c('AMT')
  
  if(METHOD=='POP' | length(List.dat)==1){
    dir.create(path=paste(output.dir,'/',substring(table,1,nchar(table)-4),'_PLOTS',sep=''),
               showWarnings=F)
    setwd(paste(output.dir,'/',substring(table,1,nchar(table)-4),'_PLOTS', sep=''))
  }else{
    dir.create(path=paste(output.dir,'/',length(List.dat),'_pooled_datasets_PLOTS',sep=''),
               showWarnings=F)
    setwd(paste(output.dir,'/',length(List.dat),'_pooled_datasets_PLOTS',sep=''))
    # Write table
    dat.raw[is.na(dat.raw)] <- c('')
    write.table(dat.raw, paste(length(List.dat),'_pooled_datasets.dat', sep=''),
                row.names=F, col.names=F, quote = F, sep='\t')}
  
  if(nrow(dat.bolus)!=0){
    dat[,IV.coln]       <- as.numeric(dat[,IV.coln])
    dat.bolus[,IV.coln] <- as.numeric(dat.bolus[,IV.coln])
    dat.bolus[,AMT]     <- as.numeric(dat.bolus[,AMT])
  } else {
    dat[,IV.coln]       <- as.numeric(dat[,IV.coln])
    dat.inf[,IV.coln]   <- as.numeric(dat.inf[,IV.coln])
    dat.inf[,AMT]       <- as.numeric(dat.inf[,AMT])}
  
  # Plots  
  for (C in 1:length(CMTs)) {  
    DV.coln       <- paste('DV',CMTs[C],sep='') # Change DV column for each CMT
    dat[,DV.coln] <- as.numeric(dat[,DV.coln])
    
    # Linear Scale
    dat.obs <- dat[dat[,DV.coln]>=0 & !is.na(dat[,DV.coln]),]  # Eliminate missing values
    IDs     <- unique(dat.obs[,ID.coln])
    Totsub  <- length(IDs)
    if(METHOD=='POP'| length(List.dat)==1){
      pdf(file=paste(substring(table,1,nchar(table)-4),'_',DV.name,'_vs_',IV.name,'_Y',CMTs[C],
                     '.pdf',sep=''),width=11.7,height=8.3)
    }else{
      pdf(file=paste(length(List.dat),'_pooled_datasets_',DV.name,'_vs_',IV.name,'_Y',CMTs[C],
                     '.pdf',sep=''),width=11.7,height=8.3)}   
    IDV.Graphs(ENV$MyDir,table,dat.obs,dat.bolus,dat.inf,CMT.name,IV.coln,ID.coln,DV.coln,AMT,
               CMTs,C,IDs,Totsub,'LIN','','','ADAPT')
    
    # Log Scale
    dat.obs <- dat[dat[,DV.coln]>0 & !is.na(dat[,DV.coln]),]  # Eliminate missing values
    IDs     <- unique(dat.obs[,ID.coln])
    Totsub  <- length(IDs)
    if(METHOD=='POP'| length(List.dat)==1){
      pdf(file=paste(substring(table,1,nchar(table)-4),'_',DV.name,'_vs_',IV.name,'_Y',CMTs[C],
                     '_(Log_Scale).pdf',sep=''),width=11.7,height=8.3)
    }else{
      pdf(file=paste(length(List.dat),'_pooled_datasets_',DV.name,'_vs_',IV.name,'_Y',CMTs[C],
                     '_(Log_Scale).pdf',sep=''),width=11.7,height=8.3)}
    IDV.Graphs(ENV$MyDir,table,dat.obs,dat.bolus,dat.inf,CMT.name,IV.coln,ID.coln,DV.coln,AMT,
               CMTs,C,IDs,Totsub,'LOG','y','| Log Scale','ADAPT')
  }
  
  switch(menu(c('Generate more individual profiles',
                'Go back to the individual DV vs. IV Menu',
                'Go back to the Main Menu','Quit'), 
              title= cat(paste('\n\nThe individuals ',DV.name,' vs. ',IV.name,
                               ' plots have successfully been created in the folder:\n',
                               getwd(),'\n\nTo continue select one of the following options:',
                               sep='')))+1,
         cat('Au revoir!'), 
         IDV(METHOD),        # Option 1
         IDV.menu(),         # Option 2
         Main.menu(),        # Option 3
         cat('Au revoir!'))} # Option 4



# -----------------------GOF
GOF <- function(){
  Precision <- function(Xs,Ys,Limx,Limy,Labx,Laby,LOG,LOGlab,Algo) {
    if(Title=='Y'){Main <- paste(Laby,' vs. ',Labx,LOGlab,' | Y(',CMTs[C],'): ',
                                 Sorts.names[1,CMTs[C]+1],sep='')}else{Main <- NULL}
    if(LOG=='xy'){
    plot(x=1,y=1,type='n',xlab=paste(Labx,DV.units,'\n'),ylab=paste(Laby,DV.units), 
         xlim=Limx,ylim=Limy,log='xy',yaxt='n',xaxt='n',main=Main)
    axis(1,as.vector(1:10 %o% 10^(-5:5)),labels=F,tcl=-0.25)
    axis(1,as.vector(10*10^(-5:5)),labels=as.character(as.vector(10*10^(-5:5))))
    axis(2,as.vector(1:10 %o% 10^(-5:5)),labels=F,tcl=-0.25)
    axis(2,as.vector(10*10^(-5:5)),labels=as.character(as.vector(10*10^(-5:5))))    
    RR <- round(summary(lm(log(dat.CMT[,Ys])~log(dat.CMT[,Xs]),data=dat.CMT))$r.squared,
                digits=4)
    }else{
    plot(x=1,y=1,type='n',xlab=paste(Labx,DV.units,'\n'),ylab=paste(Laby,DV.units), 
         xlim=Limx,ylim=Limy,main=Main)
    RR <- round(summary(lm(dat.CMT[,Ys]~dat.CMT[,Xs],data=dat.CMT))$r.squared,digits=4)}
    grid(col='grey90',lty=1,lwd=0.01,equilogs=(LOG!='xy'))
    abline (a=0,b=1)
    legend('bottomright',legend=c(Laby,'Unity Line','Loess Line'),lty=c(0,1,Loess.lty),
           lwd=c(NA,1,Loess.lwd),pch=c(Mch,NA,NA),,bty='n',cex=0.9,
           col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),'black',
                 rgb(Loess.col[1],Loess.col[2],Loess.col[3],alpha=255*Alpha,maxColorValue=255)))
    points(x=dat.CMT[,Xs],y=dat.CMT[,Ys],type='p',pch=Mch,cex=1.2,
           col=rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255))
    
    if(Title=='Y'){mtext(paste(model.name,'|',niter,Algo,'iter. on', nsub, 'sub. |',subtitle),
                         line=0.3,side=3,cex=0.8)}
    mtext(paste('OFV:',OFV),cex=0.72,line=-1.1,adj=0.01, side=3)
    mtext(paste('AIC:',AIC),cex=0.72,line=-2,adj=0.01, side=3)    
    mtext(paste('BIC:',BIC),cex=0.72,line=-2.9,adj=0.01, side=3)
    mtext(substitute(paste(R^2,':',r2),list(r2=RR)),cex=0.72,line=-3.9,adj=0.01, side=3)
 
    if(Algo!='ID'){                            # Loess Line
     order.obs <- dat.CMT[order(dat.CMT[,Xs]),]
     y.loess   <- loess(order.obs[,Ys] ~ order.obs[,Xs], span=Loess.span)
     y.loess   <- predict(y.loess, order.obs[,Xs])
     lines(x=order.obs[,Xs],y=y.loess,lty=Loess.lty,lwd=Loess.lwd,
           col=rgb(Loess.col[1],Loess.col[2],Loess.col[3],alpha=255*Alpha,maxColorValue=255))}}
  
  Residuals <- function(Xs,Ys,Limx,Labx,Unitx,Laby,Algo) {
    if(Title=='Y'){Main <- paste(Laby,' vs. ',Labx,' | Y(',CMTs[C],'): ',Sorts.names[1,CMTs[C]+1],
                                 sep='')}else{Main<-NULL}
    plot(x=1,y=1,type='n',xlim=Limx,ylim=c(min(dat.CMT[,Ys],-2),max(dat.CMT[,Ys],2)),
         xlab=paste(Labx,Unitx,'\n'),ylab=Laby,main=Main)
    grid(col='grey90',lty=1,lwd=0.01)
    polygon(x=c(min(-10,Limx[1]*1.5),Limx[2]*1.2,Limx[2]*1.2,min(-10,Limx[1]*1.5)),border=NA,
            y=c(2,2,-2,-2),density=Sden,lty=1,col=rgb(Scol[1],Scol[2],Scol[3],alpha=255*0.15,maxColorValue=255))
    points(x=dat.CMT[,Xs],y=dat.CMT[,Ys],type='p',pch=Mch,cex=1.2,
           col=rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255))
    if(Sden==0){abline(h=c(-2,2),lwd=1,lty=1,
                       col=rgb(Scol[1],Scol[2],Scol[3],alpha=255*1,maxColorValue=255))
      legend(PosLegend,legend=c(Laby,'Unity Line','Loess Line'),lty=c(0,1,Loess.lty),
             lwd=c(NA,1,Loess.lwd),pch=c(Mch,NA,NA),bty='n',cex=0.9,
             col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),'black',
                   rgb(Loess.col[1],Loess.col[2],Loess.col[3],alpha=255*Alpha,maxColorValue=255)))
    }else{legend(PosLegend,legend=c(Laby,'Unity Line','Loess Line','95%CI'),
            lty=c(0,1,Loess.lty,0),lwd=c(0,1,Loess.lwd,0),pch=c(Mch,NA,NA,15),
            col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),'black',
                  rgb(Loess.col[1],Loess.col[2],Loess.col[3],alpha=255*Alpha,maxColorValue=255),
                  rgb(Scol[1],Scol[2],Scol[3],alpha=255*0.15,maxColorValue=255)),
                  pt.cex=c(1,1,1,2.4),bty='n',cex=0.9)}
    abline(h=0)
    if(Title=='Y'){
    mtext(paste(model.name,'|',niter,Algo,'iter. on',nsub,'sub. |',subtitle,sep=' '),
          line=0.3,side=3,cex=0.8)}
    mtext(paste('OFV: ',OFV,'\nAIC: ',AIC,'\nBIC: ',BIC,sep=''),cex=0.72,
          line=-3.2,adj=0.01, side=3)
    
    if(Algo!='ID'){                            # Loess Line
     order.obs <- dat.CMT[order(dat.CMT[,Xs]),]
     y.loess   <- loess(order.obs[,Ys]~order.obs[,Xs],span=Loess.span)
     y.loess   <- predict(y.loess,order.obs[,Xs])
     lines(x=order.obs[,Xs],y=y.loess,lty=Loess.lty,lwd=Loess.lwd,
           col=rgb(Loess.col[1],Loess.col[2],Loess.col[3],alpha=255*Alpha,maxColorValue=255))}}
  
  # Set up script
  setwd(ENV$WORK.DIR)  
  List.csv   <- list.files(pattern='RSD.csv$')
  List.csv   <- sapply(strsplit(List.csv,'RSD.csv', fixed=T),`[`,1)
  if(length(List.csv)>1) {
    table    <- select.list(List.csv,title=cat('\nSelect the model you want to evaluate:'))
  } else {table <- List.csv[as.numeric(1)]}
  model.name <- table
  table      <- paste(table,'RSD.csv',sep='')
  run.name   <- paste(model.name,'.run',sep='')
  
  # Get info on the model
  if(file.exists(run.name)==F) {cat('\n\nSorry but it appears that some required files are missing to perform this analyis');Main.menu()
  }else{run.file   <- readLines(run.name,warn=F)
  Algo       <- strsplit(run.file[2],'ADAPT', fixed=T)
  Algo       <- strsplit(Algo[[1]][2],' ', fixed=T)
  Algo       <- Algo[[1]][min(grep('[A-Z]',Algo[[1]]))]
  subtitle   <- run.file[max(grep('Model: ', run.file))]
  subtitle   <- strsplit(subtitle,'Model:', fixed=T)
  subtitle   <- unlist(subtitle[[1]][2])
  subtitle   <- substr(subtitle,min(c(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                                      unlist(gregexpr('[0-9]',subtitle,ignore.case=T)))
                                   [c(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                                      unlist(gregexpr('[0-9]',subtitle,ignore.case=T)))>0]),
                          max(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                              unlist(gregexpr('[0-9]',subtitle,ignore.case=T))))      
  if(nchar(model.name)+nchar(subtitle)>35){
    subtitle <- paste(substr(subtitle,1,(32-nchar(model.name))),'...',sep='')}
  niter      <- run.file[max(grep('iterations: ', run.file))]
  niter      <- strsplit(niter,':', fixed=T)
  niter      <- as.numeric(niter[[1]][grep('[0-9]',niter[[1]])]) 
  OFV        <- run.file[max(grep('Likelihood: ', run.file))]
  OFV        <- strsplit(OFV,'Likelihood:', fixed=T)
  ifelse(length(grep('[0-9]',OFV[[1]][1]))>0,
         OFV        <- as.numeric(OFV[[1]][length(OFV[[1]])]),
         OFV        <- as.numeric(OFV[[1]][length(OFV[[1]])])*2)  # Negative Log Likelihood x2
  AIC        <- run.file[max(grep('AIC: ', run.file))]
  AIC        <- strsplit(AIC,':', fixed=T)
  AIC        <- as.numeric(AIC[[1]][length(AIC[[1]])]) 
  BIC        <- run.file[max(grep('BIC: ', run.file))]
  BIC        <- strsplit(BIC,':', fixed=T)
  BIC        <- as.numeric(BIC[[1]][length(BIC[[1]])])
  
  # Import the dataset
  dat.sim    <- read.csv (table, header=T, as.is=T)
  dat.sim    <- dat.sim[dat.sim$Data!=-1,]
  nsub       <- length(unique(dat.sim[,'Individ.']))
  
  # Get the setting options
  settings   <- read.csv(paste(ENV$MyDir,'AMGET_Settings.csv',sep='/'), as.is=T)
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  DV.name    <- settings[grep('#GOF',settings[,1])+1,2]
  DV.units   <- paste('(',settings[grep('#GOF',settings[,1])+2,2],')',sep='')
  IV.name    <- settings[grep('#GOF',settings[,1])+3,2]
  IV.units   <- paste('(',settings[grep('#GOF',settings[,1])+4,2],')',sep='')
  Mcol       <- c(col2rgb(settings[grep('#GOF',settings[,1])+5,2]),settings[grep('#GOF',settings[,1])+5,2])
  Mch        <- as.numeric(settings[grep('#GOF',settings[,1])+6,2])
  Alpha      <- as.numeric(settings[grep('#GOF',settings[,1])+7,2])
  Loess.span <- as.numeric(settings[grep('#GOF',settings[,1])+8,2])
  Loess.col  <- col2rgb(settings[grep('#GOF',settings[,1])+9,2])
  Loess.lty  <- as.numeric(settings[grep('#GOF',settings[,1])+10,2])
  Loess.lwd  <- as.numeric(settings[grep('#GOF',settings[,1])+11,2])
  Scol       <- col2rgb(settings[grep('#GOF',settings[,1])+12,2])
  Sden       <- as.numeric(settings[grep('#GOF',settings[,1])+13,2])
  Layout     <- as.numeric(settings[grep('#GOF',settings[,1])+14,2])
  PosLegend  <- settings[grep('#GOF',settings[,1])+15,2]
  Title      <- toupper(settings[grep('#GOF',settings[,1])+16,2])
  Sorts.names<- settings[grep('#GOF',settings[,1])+17,]
  
  CMTs       <- unique(dat.sim[,'Output.'])
  
  # Destination Folder
  dir.create(path=paste(output.dir,'/',model.name,'_PLOTS',sep=''),showWarnings=F)
  setwd(paste(output.dir,'/',model.name,'_PLOTS',sep=''))
  
  for (C in 1:length(CMTs)) {
    pdf(file=paste('GOF_[',model.name,']_Y',CMTs[C],'.pdf', sep=''),width=11.7,height=8.3)  
    plot.new()
    mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
    mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
    mtext('Goodness of fit plots',line=-12.5, side=3, cex=2)
    mtext(paste('Run:', model.name,sep=' '),line=-14.3, side=3, cex=1.4)
    mtext(paste('Y(',CMTs[C],'): ',Sorts.names[CMTs[C]+1],sep=''),line=-15.8, side=3, cex=1.2)
    mtext('Caution: In LOG Scale Null and Negative values are excluded',line=-17, side=3, cex=0.8)
    
    switch(Layout,'1'=par(las=1,mfrow=c(1,1),mar=c(5,4,4,2),cex=1),
           '2'=par(las=1,mfrow=c(2,2),mar=c(3.6,4,4,0.6),cex=0.7))
    dat.CMT <- dat.sim[dat.sim[,'Output.']==CMTs[C],]
    
    # Define X and Y Limits
    ivlim   <- range(dat.CMT$Obser.Time)                                  # Ind Variable limits
    logind  <- c(min(c(dat.CMT$Data,dat.CMT$ModelPred.)[c(dat.CMT$Data,dat.CMT$ModelPred.)>0]), 
                 max(dat.CMT$Data,dat.CMT$ModelPred.))                    # Log scale IPRED limits
    limind  <- c(0, max(dat.CMT$Data,dat.CMT$ModelPred.))                 # Linear scale DV and IPRED limits
 
    if(Algo!='ID'){logpop  <-c(min(c(dat.CMT$Data,dat.CMT$PopModelPred.)[c(dat.CMT$Data,
                               dat.CMT$PopModelPred.)>0]), max(dat.CMT$Data,dat.CMT$PopModelPred. )) # Log scale PopModelPred. limits
                   limpop  <-c(0, max(dat.CMT$Data,dat.CMT$PopModelPred.))# Linear scale DV and PRED limits
                   limpop2 <-c(0, max(dat.CMT$PopModelPred.))}            # Linear scale PRED only limits
    
    if(ivlim[1]!=ivlim[2]){
    #  Graph DV vs. IPRED
    Precision('ModelPred.','Data',limind,limind,'Ind. Model Pred.',DV.name, '','',Algo)

    #  Graph DV vs. IPRED (Log Scale)
    if(length(dat.CMT[dat.CMT[,'Data']<=0,'Data'])!=0 | length(dat.CMT[dat.CMT[,'ModelPred.']<=0,
                                                                       'ModelPred.'])!=0) {
       cat(paste('\nNull data were removed in the Log scale DV vs. IPRED for Y(',CMTs[C],')',
                 sep=''))}
    dat.CMT <- dat.CMT[dat.CMT[,'Data']>0 & dat.CMT[,'ModelPred.']>0,] # Remove Null or neg data
    Precision('ModelPred.','Data',logind,logind,'Ind. Model Pred.',DV.name, 'xy',' (Log Scale)',Algo)
    dat.CMT <- dat.sim[dat.sim[,'Output.']==CMTs[C],] # Restaure initial dataset
    
    #  Graph DV vs. PRED
    if('PopModelPred.' %in% colnames(dat.CMT)){
      Precision('PopModelPred.','Data',limpop,limpop,'Model Pred.',DV.name, '','',Algo)
    
    #  Graph DV vs. PRED (Log Scale)
    if(length(dat.CMT[dat.CMT[,'Data']<=0,'Data'])!=0 | length(dat.CMT[dat.CMT[,'PopModelPred.']<=0,
                                                                       'PopModelPred.'])!=0) {
       cat(paste('\nNull data were removed in the Log scale DV vs. PRED for Y(',CMTs[C],')',sep=''))}
    dat.CMT <- dat.CMT[dat.CMT[,'Data']>0 & dat.CMT[,'PopModelPred.']>0,] # Remove Null or neg data
    Precision('PopModelPred.','Data',logpop,logpop,'Model Pred.',DV.name, 'xy',' (Log Scale)',Algo)
    dat.CMT <- dat.sim[dat.sim[,'Output.']==CMTs[C],]} # Restaure initial dataset
    
    #  Std. Residual vs. IV
    Residuals('Obser.Time','Std.Resid.',ivlim,IV.name,IV.units,'Std. Residuals',Algo)
    
    #  Std. Residual vs. IPRED
    Residuals('ModelPred.','Std.Resid.',limind,'Ind. Model Pred',DV.units,'Std. Residuals',Algo)

    #  Std. Residual vs. PRED
    if('PopModelPred.' %in% colnames(dat.CMT)) {
    Residuals('PopModelPred.','Std.Resid.',limpop2,'Model Pred',DV.units,'Std. Residuals',Algo)}
    
    } else {cat(paste('\nSkipped Y(',CMTs[C],'): No valid data',sep=''))}
    
    dev.off()}  # End of outputs loop 
  
  switch(menu(c('Generate more Goodness of fit plots',
                'Posthoc fit plots',
                'Parameter distribution profiles',
                'Go back to the Main Menu','Quit'), 
              title= cat(paste('\n\nThe Goodness of fit plots have successfully been created in the folder:\n',
                               getwd(),'\n\nTo continue select one of the following options:',
                               sep='')))+1,
         cat('Au revoir!'), 
         GOF(),              # Option 1
         PHF(),              # Option 2
         PRM(),              # Option 3
         Main.menu(),        # Option 4
         cat('Au revoir!'))}}# Option 5



# -----------------------PHF
PHF <- function (){
  setwd(ENV$WORK.DIR)
  List.csv   <- list.files(pattern='PLT.csv$')
  List.csv   <- sapply(strsplit(List.csv,'PLT.csv', fixed=T),`[`,1)
  
  if(length(List.csv)>1) {
    table       <-  select.list(List.csv,title=cat('\nSelect the model you want to evaluate:'))
  } else {table <- List.csv[as.numeric(1)]}
  
  model.name <- table
  table      <- paste(table,'PLT.csv',sep='')
  run.name   <- paste(model.name,'.run',sep='')
  
  # Get info on the model
  if(file.exists(run.name)==F | file.exists(table)==F) {
    cat('\n\nSorry but it appears that some required files are missing to perform this analyis')
    Main.menu()
    }else{
  cat('\n\nPloting...')
  run.file   <- readLines(run.name,warn=F)
  Algo       <- run.file[2]
  Algo       <- strsplit(Algo,'ADAPT', fixed=T)
  Algo       <- strsplit(Algo[[1]][2],' ', fixed=T)
  Algo       <- Algo[[1]][min(grep('[A-Z]',Algo[[1]]))]
  subtitle   <- run.file[max(grep('Model: ', run.file))]
  subtitle   <- strsplit(subtitle,'Model:', fixed=T)
  subtitle   <- unlist(subtitle[[1]][2])
  subtitle   <- substr(subtitle,min(c(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                                      unlist(gregexpr('[0-9]',subtitle,ignore.case=T)))
                                    [c(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                                       unlist(gregexpr('[0-9]',subtitle,ignore.case=T)))>0]),
                       max(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                           unlist(gregexpr('[0-9]',subtitle,ignore.case=T))))
  if(nchar(model.name)+nchar(subtitle)>45){subtitle <- paste(substr(subtitle,
                                                       1,(42-nchar(model.name))),'...',sep='')}
  niter      <- run.file[max(grep('iterations: ', run.file))]
  niter      <- strsplit(niter,':', fixed=T)
  niter      <- as.numeric(niter[[1]][grep('[0-9]',niter[[1]])])
  if(Algo=='ITS') {OFV <- run.file[grep('MAP Objective Function: ', run.file)]
    OFV      <- strsplit(OFV,'MAP Objective Function:', fixed=T)
  }else{OFV  <- run.file[grep('Likelihood: ', run.file)]
    OFV      <- strsplit(OFV,'Likelihood:', fixed=T)}
  
  # Import the dataset
  dat.sim    <- read.csv(table, header=T, as.is=T)
  if(Algo=='ID') {dat.sim[,'Individ.'] <- '1'
                  dat.sim[,'IndividID'] <- '_Subject_1'}
  nsub       <- length(unique(dat.sim[,'Individ.']))
  
  # Get the number of outputs
  CMTs       <- unique(grep('Y.[0-99].$',colnames(dat.sim),value=T))
  CMTs       <- strsplit(CMTs,'.', fixed=T)
  CMTs       <- unlist(CMTs)
  CMTs       <- CMTs[seq(2,length(CMTs), by=2)]
  
  # Get the setting options
  settings   <- read.csv(paste(ENV$MyDir,'AMGET_Settings.csv',sep='/'), as.is=T)
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  DV.name    <- settings[grep('#PHF',settings[,1])+1,2]
  DV.units   <- paste('(',settings[grep('#PHF',settings[,1])+2,2],')',sep='')  
  IV.name    <- settings[grep('#PHF',settings[,1])+3,2]
  IV.units   <- paste('(',settings[grep('#PHF',settings[,1])+4,2],')',sep='')
  Mcol       <- col2rgb(settings[grep('#PHF',settings[,1])+5,2])
  Mch        <- as.numeric(settings[grep('#PHF',settings[,1])+6,2])
  Alpha      <- as.numeric(settings[grep('#PHF',settings[,1])+7,2])
  OBS.col    <- col2rgb(settings[grep('#PHF',settings[,1])+8,2])
  OBS.lty    <- as.numeric(settings[grep('#PHF',settings[,1])+9,2])
  OBS.lwd    <- as.numeric(settings[grep('#PHF',settings[,1])+10,2])
  PRED.col   <- col2rgb(settings[grep('#PHF',settings[,1])+11,2])
  PRED.lty   <- as.numeric(settings[grep('#PHF',settings[,1])+12,2])
  PRED.lwd   <- as.numeric(settings[grep('#PHF',settings[,1])+13,2])
  POP.line   <- toupper(settings[grep('#PHF',settings[,1])+14,2])
  POP.col    <- col2rgb(settings[grep('#PHF',settings[,1])+15,2])
  POP.lty    <- as.numeric(settings[grep('#PHF',settings[,1])+16,2])
  POP.lwd    <- as.numeric(settings[grep('#PHF',settings[,1])+17,2])  
  Layout     <- settings[grep('#PHF',settings[,1])+18,2]
  PosLegend  <- settings[grep('#PHF',settings[,1])+19,2]
  Log        <- toupper(settings[grep('#PHF',settings[,1])+20,2])
  Param      <- toupper(settings[grep('#PHF',settings[,1])+21,2])
  Title      <- toupper(settings[grep('#PHF',settings[,1])+22,2])  
  CMT.names  <- settings[grep('#PHF',settings[,1])+23,2:(length(CMTs)+1)]
  
  # Import IND.csv and get PRM names
  if(Algo!='ID' & Param=='Y'){
  dat.prm    <- read.csv(paste(model.name,'IND.csv',sep=''), header=F,skip=3, as.is=T)  
  nprm       <- as.numeric(dat.prm[1,2])
  prm.names  <- dat.prm[2,(1:nprm)+2]
  for(p in 1:nprm) {
    temp  <- unlist(strsplit(prm.names[1,p],' '))
    prm.names[1,p] <- temp[grep('[A-Z]',temp,ignore.case=T)]
    prm.names[1,p] <- substring(prm.names[1,p],2,nchar(prm.names[1,p]))}}
  
  # Import RSD.csv
  if(Algo!='ID'){
  dat.pop    <- read.csv(paste(model.name,'RSD.csv',sep=''), header=T, as.is=T)
  dat.pop    <- dat.pop[,c(1,3,5,11)]
  dat.pop    <- dat.pop[order(dat.pop[,2]),]}

  # Destination Folder
  dir.create(path=paste(output.dir,'/',model.name,'_PLOTS',sep=''),showWarnings=F)
  setwd(paste(output.dir,'/',model.name,'_PLOTS',sep=''))
  
  # Extract the subject ID
  IDs <- strsplit(unique(dat.sim$IndividID),' ', fixed=T)
  IDs <- unlist(IDs)
  IDs <- grep('[A-Z]',IDs,ignore.case=T,value=T)
  IDs <- substring(IDs,2,nchar(IDs))
  
  # Loop for each output
  for (C in 1:length(CMTs)){
    OBS.column  <- paste('Z.',CMTs[C],'.',sep='')
    PRED.column <- paste('Y.',CMTs[C],'.',sep='')
    dat.obs     <- dat.sim[,c('Individ.','Obser..Time',OBS.column)]  
    dat.obs     <- dat.obs[!is.na(dat.obs[,OBS.column]),]
    dat.obs     <- dat.obs[dat.obs[,OBS.column]!=-1,]
    R2          <- run.file[grep('R-squared', run.file)+C]
    
    if(Log=='Y'){
      pdf(file=paste('Posthoc_fits_[',model.name,']_Y',CMTs[C],'_(Log_Scale).pdf',
                      sep=''),width=11.7,height=8.3)
      plot.new()
      mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
      mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
      mtext('Posthoc fit plots (Log Scale)',line=-12.5, side=3, cex=2)
      mtext(paste('Run:', model.name,sep=' '),line=-14.3, side=3, cex=1.4)
      mtext(paste('Y(',CMTs[C],'): ',CMT.names[as.numeric(CMTs[C])],sep=''),line=-15.8, 
            side=3, cex=1.2)
      mtext('Caution: In LOG Scale Null and Negative values are excluded',line=-17, side=3, cex=0.8)  
      }else{
      pdf(file=paste('Posthoc_fits_[',model.name,']_Y',CMTs[C],'.pdf', sep=''),width=11.7,height=8.3) 
      plot.new()
      mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
      mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
      mtext('Posthoc fit plots',line=-12.5, side=3, cex=2)
      mtext(paste('Run:', model.name,sep=' '),line=-14.3, side=3, cex=1.4)
      mtext(paste('Y(',CMTs[C],'): ',CMT.names[as.numeric(CMTs[C])],sep=''),line=-15.8, 
            side=3, cex=1.2)}
    
    if(Algo!='ID' & Param=='Y'){
     switch(Layout,'1'={par(las=1,mfrow=c(1,1),mar=c(5,4,4,7),cex=1);Cex=1},
            '2'={par(las=1,mfrow=c(2,2),mar=c(3.6,4,4,7),cex=0.72);Cex=0.72},
            '3'={par(las=1,mfrow=c(3,3),mar=c(3.6,4,2,7),cex=0.55);Cex=0.55})
    }else{
     switch(Layout,'1'={par(las=1,mfrow=c(1,1),mar=c(5,4,4,1),cex=1);Cex=1},
            '2'={par(las=1,mfrow=c(2,2),mar=c(3.6,4,4,1),cex=0.72);Cex=0.72},
            '3'={par(las=1,mfrow=c(3,3),mar=c(3.6,4,2,1),cex=0.55);Cex=0.55})}
      
  # Prepare the graphics
  PHF.graph  <- function(Log,Log.Main){  
  if(Algo=='ID'){ivlim   <- range(ID.sim$Plot.Time,ID.obs$Obser..Time)
                 limind  <- range(ID.sim[,PRED.column],ID.obs[,OBS.column])
  }else{
    ivlim   <- range(ID.sim$Plot.Time,ID.obs$Obser..Time,
                     dat.pop[dat.pop[,1]==a & dat.pop[,2]==CMTs[C],3])
    limind  <- range(ID.sim[,PRED.column],ID.obs[,OBS.column],
                     dat.pop[dat.pop[,1]==a & dat.pop[,2]==CMTs[C],4])}
  if(Title=='Y'){Main<-paste(IDs[a], ' | ', nobs, ' observations ',Log.Main,'| Y(',CMTs[C],'): ',
                             CMT.names[C],sep='')}
  if(Log=='Y'){ 
    plot(x=1,y=1,type='n',xlim=ivlim,ylim=limind,log='y',xlab=paste(IV.name, IV.units,'\n'), 
         ylab=paste(DV.name, DV.units),main=Main,yaxt='n')
    axis(2,as.vector(1:10 %o% 10^(-5:5)),labels=F,tcl=-0.25)
    axis(2,as.vector(10*10^(-5:5)),labels=as.character(as.vector(10*10^(-5:5))))
  }else{
    plot(x=1,y=1,type='n',xlim=ivlim,ylim=limind,xlab=paste(IV.name, IV.units,'\n'), 
         ylab=paste(DV.name, DV.units),main=Main)}
  grid(col='grey90',lty=1,lwd=0.01,equilogs=(Log!='Y'))
  lines(x=ID.sim$Plot.Time,y=ID.sim[,PRED.column],lty=PRED.lty,lwd=PRED.lwd,
        col=rgb(PRED.col[1],PRED.col[2],PRED.col[3],alpha=255*Alpha,maxColorValue=255))
  if(Title=='Y'){if(Layout!=3){mtext(paste(model.name,'|',niter,Algo,'iter. |',
                                           subtitle),cex=Cex,line=0.3,side=3)}
  }else{mtext(IDs[a],cex=Cex,side=1,line=-1)}
  if(Algo=='ITS'){mtext(paste('MAP OFV:',iOFV),cex=Cex*0.9,line=-1,side=3)
                  mtext(substitute(paste(R^2,': ',r2),list(r2=iR2)),cex=Cex*0.9,line=-1.8,side=3)          
  }else{mtext(paste('OFV:',iOFV),cex=Cex*0.9,line=-1,side=3)
        mtext(substitute(paste(R^2,': ',r2),list(r2=iR2)),cex=Cex*0.9,line=-1.8,side=3)}
  if(Algo!='ID' & Param=='Y'){mtext(parameters,line=0.3,side=4,cex=Cex)}
  if(Algo!='ID' & POP.line=='Y'){
    lines(x=dat.pop[dat.pop[,1]==a & dat.pop[,2]==CMTs[C] & dat.pop[,4]!=0,3],
          y=dat.pop[dat.pop[,1]==a & dat.pop[,2]==CMTs[C] & dat.pop[,4]!=0,4],
          col=rgb(POP.col[1],POP.col[2],POP.col[3],alpha=255*Alpha,maxColorValue=255),
          lty=POP.lty,lwd=POP.lwd)
    legend(PosLegend,legend=c('Observations','Indiv. Predictions','Pop Predictions'),
           lty=c(NA,PRED.lty,POP.lty),lwd=c(1,PRED.lwd,POP.lwd),pch=c(Mch,NA,NA),
           col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),
                 rgb(PRED.col[1],PRED.col[2],PRED.col[3],alpha=255*Alpha,maxColorValue=255),
                 rgb(POP.col[1],POP.col[2],POP.col[3],alpha=255*Alpha,maxColorValue=255)),bty='n',cex=0.9)
  }else{
    legend(PosLegend,legend=c('Observations','Indiv. Predictions'),lty=c(NA,PRED.lty),
           lwd=c(1,PRED.lwd),pch=c(Mch,NA),bty='n',cex=0.9,
           col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),
                 rgb(PRED.col[1],PRED.col[2],PRED.col[3],alpha=255*Alpha,maxColorValue=255)))}
     lines(x=ID.obs$Obser..Time,y=ID.obs[,OBS.column],lty=OBS.lty,lwd=OBS.lwd,
           col=rgb(OBS.col[1],OBS.col[2],OBS.col[3],alpha=255*Alpha,maxColorValue=255))
     points(x=ID.obs$Obser..Time,y=ID.obs[,OBS.column],pch=Mch,cex=1.2,
            col=rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255))}    
    
    # Loop for each subject
    for (a in 1:length(IDs)){
      ID.sim <- dat.sim[dat.sim$Individ.==a,]
      ID.obs <- dat.obs[dat.obs$Individ.==a,]
      nobs   <- length(ID.obs[,OBS.column])
      iR2    <- strsplit(R2[[a]],' ',fixed=T)
      iR2    <- unlist(iR2)
      iR2    <- iR2[grep('[0-9]',iR2)]
      if((Algo=='ITS' & length(iR2)<3)|(Algo!='ITS' & length(iR2)<2)){iR2='n/a'
      }else{iR2 <- as.numeric(iR2[2])}
      if(Algo=='ITS'){iOFV <- as.numeric(OFV[[a]][length(OFV[[1]])])
      }else{   ifelse(length(grep('[0-9]',OFV[[a]][1]))>0,
               iOFV        <- as.numeric(OFV[[a]][length(OFV[[1]])]),
               iOFV        <- as.numeric(OFV[[a]][length(OFV[[1]])])*2)} # Negative Log Likelihood x2
      
      if(Algo!='ID' & Param=='Y'){
      parameters <- paste(prm.names,':',signif(as.numeric(dat.prm[dat.prm[,1]==a,(1:nprm)+2]),
                                               digits=3),'\n',collapse='')
      parameters <- paste('Parameters:\n',parameters,sep='')}
      Main       <- NULL
      
      if(Log=='Y'){
        ID.obs  <- ID.obs[ID.obs[,OBS.column]>0,]
        ID.sim  <- ID.sim[ID.sim[,PRED.column]>0,]
        PHF.graph('Y','(Log scale) ')
      }else{PHF.graph('N','')}}
  
    dev.off()}
  
  switch(menu(c('Generate more Posthoc fit plots',
                'Goodness of fit plots',
                'Parameter distribution profiles',
                'Go back to the main menu','Quit'), 
              title= cat(paste('\n\nThe Posthoc fit plots have successfully been created in the folder:\n',
                               getwd(),'\n\nTo continue select one of the following options:', sep='')))+1,
         cat('Au revoir!'), 
         PHF(),              # Option 1
         GOF(),              # Option 2
         PRM(),              # Option 3
         Main.menu(),        # Option 4        
         cat('Au revoir!'))}}# Option 5



# -----------------------PRM
PRM <- function (){
  setwd(ENV$WORK.DIR)
  List.csv   <- list.files(pattern='IND.csv$')
  List.csv   <- sapply(strsplit(List.csv,'IND.csv', fixed=T),`[`,1)
  
  if(length(List.csv)>1) {
    table        <-  select.list(List.csv,title=cat('\nSelect the model you want to evaluate:'))
  } else {table <- List.csv[as.numeric(1)]}
  
  model.name <- table
  table      <- paste(table,'IND.csv',sep='')
  run.name   <- paste(model.name,'.run',sep='')
  
  # Get info on the model
  if(file.exists(run.name)==F) {
    cat('\n\nSorry but it appears that some required files are missing to perform this analyis')
    Main.menu()
  }else{
  
  run.file   <- readLines(run.name,warn=F)
  Algo       <- run.file[2]
  Algo       <- strsplit(Algo,'ADAPT', fixed=T)
  Algo       <- strsplit(Algo[[1]][2],' ', fixed=T)
  Algo       <- Algo[[1]][min(grep('[A-Z]',Algo[[1]]))]
  subtitle   <- run.file[max(grep('Model: ', run.file))]
  subtitle   <- strsplit(subtitle,'Model:', fixed=T)
  subtitle   <- unlist(subtitle[[1]][2])
  subtitle   <- substr(subtitle,min(c(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                                      unlist(gregexpr('[0-9]',subtitle,ignore.case=T)))
                                    [c(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                                       unlist(gregexpr('[0-9]',subtitle,ignore.case=T)))>0]),
                       max(unlist(gregexpr('[A-Z]',subtitle,ignore.case=T)),
                           unlist(gregexpr('[0-9]',subtitle,ignore.case=T))))      
  if(nchar(subtitle)>30) {subtitle <- paste(substr(subtitle,1,27),'...',sep='')}
  niter      <- run.file[max(grep('iterations: ', run.file))]
  niter      <- strsplit(niter,':', fixed=T)
  niter      <- as.numeric(niter[[1]][grep('[0-9]',niter[[1]])]) 
  
  # Import the dataset
  dat.prm    <- read.csv(table,header=F,skip=3, as.is=T)
  nsub       <- length(dat.prm[,1])-3
  nprm       <- as.numeric(dat.prm[1,2])
  prm.names  <- dat.prm[2,(1:nprm)+2]
  for(p in 1:nprm) {
    temp  <- strsplit(prm.names[1,p],' ')
    temp  <- unlist(temp)
    prm.names[1,p] <- temp[grep('[A-Z]',temp,ignore.case=T)]
    prm.names[1,p] <- substring(prm.names[1,p],2,nchar(prm.names[1,p]))}
  dat.prm  <- dat.prm[c(-1:-3),c(-1,-2)]

  # Get the data for PRM vs. Iteration
  IT.file   <- read.csv(paste(model.name,'IT.csv',sep=''),header=T,skip=3, as.is=T)
  IT.file   <- IT.file[,colSums(is.na(IT.file[,1:ncol(IT.file)]))!=nrow(IT.file)]
  IT.file   <- IT.file[,c(1,ncol(IT.file),2:(ncol(IT.file)-1))] 
  IT.file.name <- colnames(IT.file)
  IT.file.name[2:ncol(IT.file)] <- substr(IT.file.name[2:ncol(IT.file)],2,
                                          nchar(IT.file.name[2:ncol(IT.file)]))
  IT.file.name <- strsplit(IT.file.name,'.',fixed=T)
  colnames(IT.file) <- sapply(IT.file.name, paste,collapse='')
  
  # Get the setting options
  settings   <- read.csv(paste(ENV$MyDir,'AMGET_Settings.csv',sep='/'), as.is=T)
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  Layout     <- as.numeric(settings[grep('#PRM',settings[,1])+1,2])
  nBreak     <- as.numeric(settings[grep('#PRM',settings[,1])+2,2])
  HISTO      <- toupper(settings[grep('#PRM',settings[,1])+3,2])
  Lcol       <- col2rgb(settings[grep('#PRM',settings[,1])+4,2])
  Lty        <- as.numeric(settings[grep('#PRM',settings[,1])+5,2])
  Lwd        <- as.numeric(settings[grep('#PRM',settings[,1])+6,2])
  Alpha      <- as.numeric(settings[grep('#PRM',settings[,1])+7,2])
  Title      <- toupper(settings[grep('#PRM',settings[,1])+8,2])
  
  # Destination Folder
  dir.create(path=paste(output.dir,'/',model.name,'_PLOTS',sep=''),showWarnings=F)
  setwd(paste(output.dir,'/',model.name,'_PLOTS',sep=''))
  
  # Plots
  pdf(file=paste('PRM_Distribution_[',model.name,'].pdf', sep=''),width=11.7,height=8.3)
  plot.new()
  mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
  mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
  mtext('Model parameters distribution',line=-12.5, side=3, cex=2)
  mtext(paste(model.name,'|',subtitle,sep=' '),line=-14.3, side=3, cex=1.4)
  mtext(paste(niter,Algo,'iterations on', nsub, 'subjects |',nprm,'model parameters',sep=' '),
        line=-16, side=3, cex=1.2)
  
  if(HISTO=='Y') {
    switch(Layout,'1'=par(las=1,mfcol=c(1,1),mar=c(5,4,4,2),cex=0.9),
           '2'=par(las=1,mfcol=c(2,2),mar=c(4,4,4,0.6),cex=0.8))
  }else{switch(Layout,'1'=par(las=1,mfrow=c(1,1),mar=c(5,4,4,2),cex=0.9),
               '2'=par(las=1,mfrow=c(2,2),mar=c(4,4,4,0.6),cex=0.8))}
  Main <- ''

  # Density
  for (P in 1:nprm) {
    dat.prm[,P] <- as.numeric(dat.prm[,P])
    if(Title=='Y'){Main<-paste('Distribution of ',prm.names[1,P],sep='')}
    plot(density(dat.prm[,P]),type='n',main=Main)
    grid(col='grey90',lty=1,lwd=0.01)
    lines(density(dat.prm[,P]),lwd=Lwd,lty=Lty,
          col=rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255))
    if(Title=='Y'){mtext(paste('Range:',signif(min(dat.prm[,P]),digits=4),'-',
                        signif(max(dat.prm[,P]),digits=4),'| Mean:',
                        signif(mean(dat.prm[,P]),digits=4),'| Median:', 
                        signif(median(dat.prm[,P]),digits=4), sep=' '),
                         cex=0.8,line=0.4, side=3)
    }else{mtext(prm.names[1,P],cex=0.8,side=1, line=-1.2)}
    
    # Histogram
    if(HISTO=='Y'){BREAKS <- seq(min(dat.prm[,P]),max(dat.prm[,P]),length.out=nBreak)
    if(Title=='Y'){Main<-paste('Distribution of ',prm.names[1,P],sep='')}
    hist(dat.prm[,P], breaks=BREAKS, col=rainbow(nBreak, start=0.5, end=0.63),freq=T,main=Main,
         xlab=paste(prm.names[1,P],'\n'), axes=T, labels=T)
    if(Title=='Y'){mtext(paste('Range:',signif(min(dat.prm[,P]),digits=4),'-',
                         signif(max(dat.prm[,P]),digits=4),'| Mean:',
                         signif(mean(dat.prm[,P]),digits=4),'| Median:',
                         signif(median(dat.prm[,P]),digits=4)),cex=0.8,line=0.4, side=3)}
    }}
  dev.off()
  
  # PRM vs. Iteration Plots.
  pdf(file=paste('PRM_Iterations_[',model.name,'].pdf', sep=''),width=11.7,height=8.3)
  plot.new()
  mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
  mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
  mtext('Model parameters value vs. Iteration',line=-12.5, side=3, cex=2)
  mtext(paste(model.name,'|',subtitle,sep=' '),line=-14.3, side=3, cex=1.4)
  mtext(paste(niter,Algo,'iterations on', nsub, 'subjects |',nprm,'model parameters',sep=' '),
        line=-16, side=3, cex=1.2)
  
  switch(Layout,'1'=par(las=1,mfrow=c(1,1),mar=c(5,4,4,2),cex=0.9),
         '2'=par(las=1,mfrow=c(2,2),mar=c(3.6,4,4,0.6),cex=0.8))
  for(I in 2:ncol(IT.file)){
    if(Title=='Y'){Main<-paste(colnames(IT.file)[I],' values vs. Iteration',sep='')
    }else{Main <- NULL}
    plot(IT.file[,1],IT.file[,I],type='n',xlab='Iterations\n',
         ylab=paste(colnames(IT.file)[I],'values'),main=Main)
    grid(col='grey90',lty=1,lwd=0.01)
    lines(IT.file[,1],IT.file[,I],lwd=Lwd,lty=Lty,
          col=rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255))
    if(Title=='Y'){mtext(paste('Init: ',signif(IT.file[length(IT.file[is.na(IT.file[,I]),I])+1,I],
                         digits=4),' | Range: ',
                               signif(min(IT.file[,I],na.rm=T),digits=4),' to ',
                               signif(max(IT.file[,I],na.rm=T),digits=4),' | Final: ',
                               signif(IT.file[nrow(IT.file),I],digits=4),sep=''),
                               cex=0.8,line=0.4, side=3)}
    }
  dev.off()
  
  switch(menu(c('Generate more Parameter distribution profiles',
                'Goodness of fit plots','Posthoc fits',
                'Go back to the main menu','Quit'), 
              title= cat(paste('\n\nThe Parameter distribution profiles have successfully been created in the folder:\n',
                               getwd(),'\n\nTo continue select one of the following options:',
                               sep='')))+1,
         cat('Au revoir!'), 
         PRM(),              # Option 1
         GOF(),              # Option 2
         PHF(),              # Option 3
         Main.menu(),        # Option 4        
         cat('Au revoir!'))}}# Option 5



# -----------------------VPC
VPC <- function() {
  setwd(ENV$WORK.DIR)
  List.csv   <- list.files(pattern='POP.csv$')
  List.csv   <- sapply(strsplit(List.csv,'POP.csv', fixed=T),`[`,1) 
  
  if(length(List.csv)==0){
    cat('\n\nSorry but it appears that some required files are missing to perform this analyis')
    Main.menu()  
    }else{
  if(length(List.csv)>1) {
    table        <-  select.list(List.csv,title=cat('\nSelect the model you want to evaluate:'))
  }else{table <- List.csv[1]}
  
  table      <- paste(table,'POP.csv',sep='')
  dat.sim    <- read.csv(table, header=T,as.is=T)
  dat.sim    <- dat.sim[,grep('Y.[0-99].',colnames(dat.sim))]
  
  # Get info on the model
  model.name <- substring(table,1,nchar(table)-7)
  nsim       <- nrow(dat.sim)
  CMTs       <- unique(grep('Y.[0-99].$',colnames(dat.sim),value=T))
  CMTs       <- unlist(strsplit(CMTs,'.', fixed=T))
  CMTs       <- CMTs[seq(2,length(CMTs), by=2)]
  
  # Percentiles calculation
  settings   <- read.csv(paste(ENV$MyDir,'AMGET_Settings.csv',sep='/'), as.is=T)
  LowP       <- as.numeric(settings[grep('#VPC',settings[,1])+15,2])
  HighP      <- as.numeric(settings[grep('#VPC',settings[,1])+16,2])
  Percentiles<- LPercentiles <- c()
  
  for(a in 1:length(CMTs)) {
    col.dat <- grep(paste('Y.',a,'.',sep=''),colnames(dat.sim))
    dat.tmp <- dat.sim[,col.dat]
    out <- out.L <- temp <- temp.L <- c()
    for(b in 1:length(col.dat)) {
      #Lin
      temp      <- quantile(dat.tmp[,b], c(LowP,0.5,HighP), names=F)
      out       <- rbind(out,cbind(temp[1],temp[2],temp[3]))
    
      #Log
      temp.L    <- quantile(dat.tmp[dat.tmp[,b]>0,b], c(LowP,0.5,HighP), names=F)
      out.L     <- rbind(out.L,cbind(temp.L[1],temp.L[2],temp.L[3]))}
    
    Percentiles <- cbind(Percentiles,out)
    LPercentiles<- cbind(LPercentiles,out.L)}
  
  # Get the data time and observations
  List.dat <- list.files(pattern='.dat$')
  List.dat <- List.dat[List.dat!=paste(substring(model.name,1,nchar(model.name)-4),'.dat',sep='')]
  if(length(List.dat)>1) {
    table2 <- select.list(List.dat,title=cat(paste(
                          '\nSelect the datafile containing the OBSERVATIONS used by ',
                           model.name,':',sep='')))
  }else{table2 <- List.dat[1]} 
  List.dat <- table2
  table2   <- read.table(table2,stringsAsFactors=F,header=F,fill=T,col.names=1:100)
  table2   <- table2[,colSums(is.na(table2[,1:100]))!=nrow(table2)]  
  
  dat.obs <- c()
  rows <- NumID <- 0
  while((rows)<nrow(table2)) {
    IDs       <- table2[(1+rows),1]
    NumID     <- NumID+1
    Ndose     <- as.numeric(table2[(4+rows),1])
    dos.Times <- table2[(1:Ndose)+(4+rows),1]
    obs.Times <- table2[(1:table2[(6+Ndose+rows),1])+(6+Ndose+rows),1]
    obs       <- table2[(1:table2[(6+Ndose+rows),1])+(6+Ndose+rows),(1:CMTs)+1]
    temp.obs  <- data.frame(stringsAsFactors=F,
                            ID   = rep(IDs,times=length(obs.Times)),
                            COUNT= rep(NumID,times=length(obs.Times)),
                            TIME = obs.Times, DV = obs)
    colnames(temp.obs) <- c('ID','COUNT','TIME',paste('DV',1:CMTs,sep=''))
    dat.obs   <- rbind(dat.obs,temp.obs) 
    rows = (rows + 6 + length(obs.Times) + length(dos.Times))}

  dat.obs[,3] <- as.numeric(dat.obs[,3])
  IVar        <- unique(dat.obs[,'TIME'])
  dat.sim     <- cbind(IVar,Percentiles)
  dat.sim.L   <- cbind(IVar,LPercentiles)
  
  # Get the setting options
  output.dir <- settings[grep('#DIR',settings[,1])+1,2]
  DV.name    <- settings[grep('#VPC',settings[,1])+1,2]
  DV.units   <- paste('(',settings[grep('#VPC',settings[,1])+2,2],')',sep='')
  IV.name    <- settings[grep('#VPC',settings[,1])+3,2]
  IV.units   <- paste('(',settings[grep('#VPC',settings[,1])+4,2],')',sep='')
  Mcol       <- col2rgb(settings[grep('#VPC',settings[,1])+5,2])
  Mch        <- as.numeric(settings[grep('#VPC',settings[,1])+6,2])
  Alpha      <- as.numeric(settings[grep('#VPC',settings[,1])+7,2])
  Lcol       <- col2rgb(settings[grep('#VPC',settings[,1])+8,2])
  Lty        <- as.numeric(settings[grep('#VPC',settings[,1])+9,2])
  Lwd        <- as.numeric(settings[grep('#VPC',settings[,1])+10,2])
  Scol       <- col2rgb(settings[grep('#VPC',settings[,1])+11,2])
  Sden       <- as.numeric(settings[grep('#VPC',settings[,1])+12,2])
  Sty        <- as.numeric(settings[grep('#VPC',settings[,1])+13,2])
  Layout     <- as.numeric(settings[grep('#VPC',settings[,1])+14,2])
  Interval   <- as.numeric(settings[grep('#VPC',settings[,1])+17,2])
  PosLegend  <- settings[grep('#VPC',settings[,1])+18,2]
  Title      <- toupper(settings[grep('#VPC',settings[,1])+19,2])
  Sorts.names<- unlist(settings[grep('#VPC',settings[,1])+20,2:(length(CMTs)+1)])
  
  # Create output dir
  dir.create(path=paste(output.dir,'/',model.name,'_PLOTS',sep=''),showWarnings=F)
  setwd(paste(output.dir,'/',model.name,'_PLOTS',sep=''))
  
  VPC.graph  <- function(DAT.SIM,Log,Name1,Name2,Name3){
    pdf(file=paste('VPC_[',model.name,']',Name1,'.pdf', sep=''),width=11.7,height=8.3)
    plot.new()
    mtext('AMGET',line=-36.5, side=3, cex=0.8,adj=-0.065)
    mtext(format(Sys.time(), '%d %b %Y %H:%M'),line=-36.5, side=3, cex=0.8,adj=1.02)
    mtext(paste('Visual Predictive Check',Name2),line=-12.5, side=3, cex=2)
    mtext(paste(DV.name,'vs.',IV.name,sep=' '),line=-14.3, side=3, cex=1.4)
    mtext(paste(nsim,'simulations | Sim:', table,'| Obs:',List.dat,sep=' '),line=-15.8,side=3,cex=1)
    
    switch(Layout,'1'=par(las=1,mfrow=c(1,1),mar=c(5,4,4,2),cex=0.9),
           '2'=par(las=1,mfrow=c(2,2),mar=c(3.6,4,4,0.6),cex=0.7),
           '3'=par(las=1,mfrow=c(3,3),mar=c(3.6,4,2,0.6),cex=0.55))
    for (C in 1:length(CMTs)) {
      OBS.column  <- paste('DV',CMTs[C],sep='') 
      if(Log=='Y'){dat.cmt <- dat.obs[!is.na(dat.obs[,OBS.column]) & dat.obs[,OBS.column]>0,]  
      }else{dat.cmt <- dat.obs[!is.na(dat.obs[,OBS.column]) & dat.obs[,OBS.column]!=-1,]}
      timeclass   <- cut(dat.cmt[,'TIME'],breaks=seq(0,max(dat.cmt[,'TIME']),by=Interval))
      dat.med     <- data.frame(Times=tapply(dat.cmt[,'TIME'],factor(timeclass),mean),
                                Medians=tapply(dat.cmt[,OBS.column],factor(timeclass),median),
                                LowP=tapply(dat.cmt[,OBS.column],factor(timeclass),quantile,LowP),
                                HighP=tapply(dat.cmt[,OBS.column],factor(timeclass),quantile,HighP))
      if(Title=='Y') {Main<-paste('VPC | ',DV.name,' vs. ',IV.name,' | Y(',CMTs[C],'): ',
                                  Sorts.names[C],Name3,sep='')}else{Main <- NULL}
      if(Log=='Y'){
        limind    <- range(dat.cmt[,OBS.column])
        ivlim     <- range(dat.sim.L[,1])
        plot(x=1,y=1,type='n',xlim=ivlim,ylim=limind,xlab=paste(IV.name,IV.units,'\n'), 
             ylab=paste(DV.name,DV.units),log='y',yaxt='n',main=Main)
        axis(2,as.vector(1:10 %o% 10^(-5:5)),labels=F,tcl=-0.25)
        axis(2,as.vector(10*10^(-5:5)),labels=as.character(as.vector(10*10^(-5:5))))
      }else{
        limind    <- c(0,max(dat.cmt[,OBS.column]))
        ivlim     <- range(dat.sim[,1])
        plot(x=1,y=1,type='n',xlim=ivlim,ylim=limind,xlab=paste(IV.name,IV.units,'\n'), 
             ylab=paste(DV.name,DV.units),main=Main)}
      grid(col='grey90',lty=1,lwd=0.01,equilogs=(Log!='Y'))
      PolyS <- cbind(DAT.SIM[,1],DAT.SIM[,(C*3)-1])
      PolyS <- rbind(cbind(DAT.SIM[,1],DAT.SIM[,(C*3)+1]),PolyS[order(-PolyS[,1]),])
      polygon(x=PolyS[,1],y=PolyS[,2],density=Sden,lty=1,border=NA,
              col=rgb(Scol[1],Scol[2],Scol[3],alpha=255*(Alpha/3),maxColorValue=255))   # Surface Simulated data
      PolyO <- cbind(dat.med[,1],dat.med[,3])
      PolyO <- rbind(cbind(dat.med[,1],dat.med[,4]),PolyO[order(-PolyO[,1]),])
      polygon(x=PolyO[,1],y=PolyO[,2],density=Sden,lty=1,border=NA,
              col=rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*(Alpha/3),maxColorValue=255))   # Surface Observed data
      points(x=dat.cmt[,'TIME'],y=dat.cmt[,OBS.column], pch=Mch,
             col=rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255))  # Observations
      if(Title=='Y'){if(Layout!=3){mtext(paste('Sim:', table,'| Obs:',List.dat,'|',nsim,'simulations'), 
                           cex=0.8, line=0.3, side=3)}
      }else{mtext(Sorts.names[C],cex=0.8,side=1, line=-1.2)}
      lines(x=dat.med[,1],y=dat.med[,2],lty=Lty,lwd=Lwd,
            col=rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255))   #  Median obs.
      lines(x=DAT.SIM[,1],y=DAT.SIM[,(C*3)],lty=Sty,lwd=Lwd,
            col=rgb(Scol[1],Scol[2],Scol[3],alpha=255*Alpha,maxColorValue=255))   #  Median sim.
      if(Sden==0){
        lines(x=DAT.SIM[,1],y=DAT.SIM[,(C*3)-1],lty=Sty,lwd=Lwd,
              col=rgb(Scol[1],Scol[2],Scol[3],alpha=255*Alpha,maxColorValue=255)) # Low percentile sim.
        lines(x=DAT.SIM[,1],y=DAT.SIM[,(C*3)+1],lty=Sty,lwd=Lwd,
              col=rgb(Scol[1],Scol[2],Scol[3],alpha=255*Alpha,maxColorValue=255)) # High percentile sim.
        lines(x=dat.med[,1],y=dat.med[,3],lty=Lty,lwd=Lwd,
              col=rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255)) # Low percentile Obs.
        lines(x=dat.med[,1],y=dat.med[,4],lty=Lty,lwd=Lwd,
              col=rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255)) # High percentile Obs.
        legend(PosLegend,legend=c('Observations',
               paste(LowP*100,'th, 50th and ',HighP*100,'th %iles of Obs.',sep=''),
               paste(LowP*100,'th, 50th and ',HighP*100,'th %iles of Sim.',sep='')),
               lty=c(NA,Lty,Sty),lwd=c(0,Lwd,Lwd),pch=c(Mch,NA,NA),bty='n',cex=0.9,
               col=c(rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255),
                     rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255),
                     rgb(Scol[1],Scol[2],Scol[3],alpha=255*Alpha,maxColorValue=255)))
      }else{legend(PosLegend,legend=c('Observations','Median of Obs.','Median of Sim.',
                                      paste((HighP-LowP)*100,'% PI of Obs.',sep=''),
                                      paste((HighP-LowP)*100,'% PI of Sim.',sep='')),
                   pch=c(Mch,NA,NA,15,15),pt.cex=c(1,1,1,2.4,2.4),lwd=c(0,Lwd,Lwd,0,0),
                   lty=c(NA,Lty,Sty,NA,NA),bty='n',cex=0.9,
                   col=c(rgb(Mcol[1],Mcol[2],Mcol[3],alpha=255*Alpha,maxColorValue=255),
                         rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*Alpha,maxColorValue=255),
                         rgb(Scol[1],Scol[2],Scol[3],alpha=255*Alpha,maxColorValue=255),
                         rgb(Lcol[1],Lcol[2],Lcol[3],alpha=255*(Alpha/3),maxColorValue=255),
                         rgb(Scol[1],Scol[2],Scol[3],alpha=255*(Alpha/3),maxColorValue=255)))}}
    dev.off()}

    VPC.graph(dat.sim,'N','','','')                          # Linear
    if(Sden<=0){VPC.graph(dat.sim.L,'Y','_(Log)','(Log Scale)',' | Log')} # Log
    
  switch(menu(c('Generate more Visual Predictive Check',
                'Go back to the Main Menu','Quit'), 
              title= cat(paste('\n\nThe Visual Predictive Check have successfully been created in the folder:\n',
                               getwd(),'\n\nTo continue select one of the following options:', 
                               sep='')))+1,
         cat('Au revoir!'), 
         VPC(),               # Option 1
         Main.menu(),         # Option 2
         cat('Au revoir!'))}} # Option 3

 SetWD()}
