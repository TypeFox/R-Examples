#Load the required packages. They need to be installed first.
require(sem)
require(igraph)
#require(Cairo)
#require(RJDBC) #only for database!
require(tcltk)
tclRequire("BWidget")

############################################## DATA functions ########################################################
## To change the working directory ##
dpa.data.setWorkingDirectory <- function()
{
 baseDir <-tclvalue(tkchooseDirectory())
 setwd(baseDir)
}

## To load data either from the disk or database ##
dpa.data.loadData<- function()
{
 ak1 <- tktoplevel()
 tktitle(ak1) <- "  Load data  "
 assign("ak1",ak1,env=.GlobalEnv)
 DISK.but <-tkbutton(ak1,text="   Load data from disk  ",command=dpa.data.loadDataFromDisk )
 tkconfigure(DISK.but,width=20,foreground="black",background="white")
 tkgrid(DISK.but,row=0,column=0,padx=10,pady=10)
 DATABASE.but <-tkbutton(ak1,text="   Load data from database  ",command=dpa.data.loadDataFromDatabase)
 tkconfigure(DATABASE.but,width=20,foreground="black",background="white") 
 tkgrid(DATABASE.but,row=1,column=0,padx=10,pady=10)
 CANCEL.but <-tkbutton(ak1,text="   Cancel  ",command=dpa.data.loadCancel)
 tkgrid(CANCEL.but,row=2,column=1,padx=10,pady=10)
}

## To cancel the loading of data ##
dpa.data.loadCancel<-function()
{ 
 try(tkdestroy(ak1),TRUE)
}

## To load data file available on the disk ##
dpa.data.loadDataFromDisk <- function(dataFileName=NULL)
{
 if (is.null(dataFileName)){
   try(tkconfigure(ak1,cursor="watch"),TRUE)
   dataFileName <- tclvalue(tkgetOpenFile())
 }
 assign("dataFileName",dataFileName,env=.GlobalEnv)
 
 if (!nchar(dataFileName)){
   tkconfigure(ak1,cursor="arrow")
   tkmessageBox(message = "No file was selected.")
 } else {
   if(is.data.frame(e)==FALSE){

     #sort out the extension of the file
     filesplit<-unlist(strsplit(dataFileName, ".", fixed = TRUE))
     extension<-filesplit[length(filesplit)]
     #print(extension)
     
     if(extension=="rda"){
       #print(paste("Loading rda file"))
       #load .RDA file
	 load(dataFileName, .GlobalEnv)
     }
     if(extension=="csv"){
       #print(paste("Loading csv file"))
       #load .csv file
       e <- read.csv(file=dataFileName)  
       assign("e",e,env=.GlobalEnv)
     }
     try(tkconfigure(ak1,cursor="arrow"),TRUE)
     try(tkdestroy(ak1),TRUE)
     #}
    }
   else{
     message <- tkmessageBox(message="There is already a data loaded. Do you want to over-write it?",icon="question",type="yesnocancel",default="yes")
     if(as.character(message)=="yes"){
       if(extension=="rda"){
         #load .RDA file
         print(paste("Loading rda file"))
  	   load(dataFileName, .GlobalEnv)
       }
       if(extension=="csv"){
         #load .csv file
         print(paste("Loading csv file"))
         e <- read.csv(file=dataFileName)  
         assign("e",e,env=.GlobalEnv)
       }
       assign("e",e,env=.GlobalEnv)
       try(tkconfigure(ak1,cursor="arrow"),TRUE)
       try(tkdestroy(ak1),TRUE)
       }  
    }
  }
dataFileName<-NULL;	
}

## To load the data from a database after getting required information from user ##
dpa.data.loadDataFromDatabase<- function()
{
 ak2<-tktoplevel()
 tktitle(ak2) <- "  Authentication Required  "
 assign("ak2",ak2,env=.GlobalEnv)
 dbServerType<-tclVar()
 assign("dbServerType",dbServerType,env=.GlobalEnv)
 dbServer   <- tclVar()
 assign("dbServer",dbServer,env=.GlobalEnv)
 dbName<-tclVar()
 assign("dbName",dbName,env=.GlobalEnv)
 dbTable    <- tclVar()
 assign("dbTable",dbTable,env=.GlobalEnv)
 UserName <- tclVar()
 assign("UserName",UserName,env=.GlobalEnv)
 Password <- tclVar()
 assign("Password",Password,env=.GlobalEnv)
  
 entry.dbServerType <- tkentry(ak2,width="30",textvariable=dbServerType)
 entry.dbServer     <- tkentry(ak2,width="30",textvariable=dbServer)
 entry.dbName       <- tkentry(ak2,width="30",textvariable=dbName)
 entry.dbTable      <- tkentry(ak2,width="30",textvariable=dbTable)
 entry.UserName     <- tkentry(ak2,width="30",textvariable=UserName)
 entry.Password     <- tkentry(ak2,width="30",textvariable=Password)
  
 tkgrid(tklabel(ak2,text="Server Type"),row=0,column=0)
 tkgrid(entry.dbServerType,padx=10,pady=10,row=0,column=1)
 tkgrid(tklabel(ak2,text="Server"),row=1,column=0)
 tkgrid(entry.dbServer,padx=10,pady=10,row=1,column=1)
 tkgrid(tklabel(ak2,text="Data Base Name"),row=2,column=0)
 tkgrid(entry.dbName,padx=10,pady=10,row=2,column=1)
 tkgrid(tklabel(ak2,text="Data Base Table"),row=3,column=0)
 tkgrid(entry.dbTable,padx=10,pady=10,row=3,column=1)
 tkgrid(tklabel(ak2,text="User Name"),row=4,column=0)
 tkgrid(entry.UserName,padx=10,pady=10,row=4,column=1)
 tkgrid(tklabel(ak2,text="Password"),row=5,column=0)
 tkgrid(entry.Password,padx=10,pady=10,row=5,column=1)
  
 NoteFont<-tkfont.create(family="times",size=12)
 tkgrid(tklabel(ak2,text="Note : This is a time consuming step. Consider saving data in the next step in order to allow quick loading later.",font=NoteFont),row=6,columnspan=4) 
  
 SUBMIT.but <-tkbutton(ak2,text="   Submit   ",command=dpa.data.authenticationSubmit)
 tkgrid(SUBMIT.but,row=7,column=1,padx=10,pady=10)
 
 AuthenticationCancel.but <-tkbutton(ak2,text="   Cancel  ",command=dpa.data.authentificationCancel)
 tkgrid(AuthenticationCancel.but,row=7,column=2,padx=10,pady=10)
}

## To submit the required information in order to access the database ##
dpa.data.authenticationSubmit <- function()
{
 dbServerTypeVal <- tclvalue(dbServerType)
 dbServerVal <- tclvalue(dbServer)
 dbName <- tclvalue(dbName)
 dbTableVal <- tclvalue(dbTable)
 print(dbTableVal)
 UserNameVal <- tclvalue(UserName)
 PasswordVal <- tclvalue(Password)  
 
 if(dbServerTypeVal=="PostGres") {
   drv <- JDBC("org.postgresql.Driver","c:/postgresql-8.3-604.jdbc4.jar",identifier.quote="`")
   #Make a connection!
   #conn <- dbConnect(drv, "jdbc:postgresql://test.eeni.tbm.tudelft.nl:5432/emiledb", usr, pwd)
   servertable <- paste("jdbc:postgresql://",dbServerVal,"/",dbName);
   conn <- dbConnect(drv, servertable, UserNameVal ,PasswordVal)
   #query <- paste("select * from ",dbTableVal LIMIT 10 ,";")
   dbGetQuery(conn,"select 10 from dbTableVal")
   e<-dbReadTable(conn,"dbTableVal")
  }
}

## To cancel the authentication screen ##
dpa.data.authenticationCancel<- function()
 { tkdestroy(ak2)}

## To edit the data to either to have a look or to modify it ##
dpa.data.viewOrEditData<- function()
{
 e <- fix(e)
}

## To first sort the data and then find missing rows and terms in the original data ##
dpa.data.checkData <- function()
{
 ak3 <- tktoplevel()
 tktitle(ak3) <- " Sort and check the Data  "
 dpa.data.checkData.sort<-function()
 { 
 #the columns with whose respect data is ordered to be entered by user.
 dpa.sort.data(e,"tick","job") 
 }
 Sort.but <-tkbutton(ak3,text=" sort the data  ",command=dpa.data.checkData.sort)
 tkconfigure(Sort.but,width=30,foreground="black",background="white")
 tkgrid(Sort.but,row=0,column=0)

 dpa.data.missingRow <- function()
{
 dpa.find.missingRow(e,"tick","job")
 }
 MissingRow.but <-tkbutton(ak3,text=" Find missing rows from the data ",command=dpa.data.missingRow)
 tkconfigure(MissingRow.but,width=30,foreground="black",background="white")
 tkgrid(MissingRow.but,row=1,column=0)
  
 dpa.data.missingData <- function()
 {
 dpa.locate.missing(e,"tick","job")
 }
 MissingData.but <-tkbutton(ak3,text=" Find missing data from the existing rows ",command=dpa.data.missingData)
 tkconfigure(MissingData.but,width=30,foreground="black",background="white")
 tkgrid(MissingData.but,row=2,column=0) 
}
 
## To save the data to the disk ##
dpa.data.saveDataToDisk <- function(dataFileName=NULL)
{
 if (is.null(dataFileName)){
 dataFileName<-tclvalue(tkgetSaveFile())
  }
 if (!nchar(dataFileName)){
   tkmessageBox(message="No file was selected!")
  }
 else{
   save(e, file=dataFileName)
  }
}

#########################################  Relations functions  ########################################
#To load the relations from the disk
dpa.relations.loadRelations<- function(loadRelFileName=NULL)
{
 if(is.null(loadRelFileName)){
   loadRelFileName <- tclvalue(tkgetOpenFile())
 }
 if (!nchar(loadRelFileName)) {
   tkmessageBox(message = "No file was selected.")}
 else {
  if(is.data.frame(relations)==FALSE){
   #load relations from disk.
   relations <-read.table(file=loadRelFileName,header=TRUE)
   assign("relations",relations,env=.GlobalEnv)
  }
  else{ 
    messageRel <- tkmessageBox(message="There is already a relation loaded. Do you want to over-write it?",icon="question",type="yesnocancel",default="yes")
    if(as.character(messageRel)=="yes"){
    #load relations from disk.
    relations <-read.table(file=loadRelFileName,header=TRUE)
    assign("relations",relations,env=.GlobalEnv)
    }
   }
 }	
}

dpa.relations.saveRelations<- function(saveRelFileName=NULL)
{
 if (is.null(saveRelFileName)){
   saveRelFileName<-tclvalue(tkgetSaveFile())
 }
 if (!nchar(saveRelFileName)){
   tkmessageBox(message="No file was selected!")}
 else{
   write.table(relations,file=saveRelFileName,row.names=F) 
  }
}

#######  Lag generating function  #######
dpa.generate.lag <-function(dataframe=NULL,tickColumn=NULL,sourceColumn=NULL,minLag=1,maxLag=1)
{
 attach(dataframe)
 on.exit(detach(dataframe))
# tick <- eval(parse(text=tickColumn))
#print(tick)
 source <- eval(parse(text=sourceColumn))
 newColumn <- source
 n=(nrow(dataframe))/(max(dataframe[tickColumn]))   
 for(i in as.numeric(minLag):as.numeric(maxLag)){
  newName <- paste(sourceColumn,"_L",i,sep="")
  if(is.na(match(newName,names(dataframe)))){
  print("Generating lag...");
  for(k in 1:n)
  {
  limitMax=k*max(dataframe[tickColumn])
  limitMin=(k-1)*max(dataframe[tickColumn])+i+1       
   for(j in limitMax:limitMin)
    { newColumn[j]<- newColumn[j-1] }
   }
 #To add new column to the original data
 dataframe<-data.frame(dataframe,newColumn)      
 #To rename the column just added to the dataframe
 names( dataframe )[ which( names( dataframe ) == "newColumn" ) ] <- newName
 }
 }
 assign("dataframe",dataframe,env=.GlobalEnv)
}

#######   A function to Sort the original data    ########
dpa.sort.data<-function(dataframe=NULL,tickColumn=NULL,runColumn=NULL)
{
  attach(dataframe)
  on.exit(detach(dataframe))
  tick <- eval(parse(text=tickColumn))
  run <- eval(parse(text=runColumn)) 
  dataframe<-dataframe[order(run,tick),]   
}

#A function to find the missing data in existing rows
dpa.locate.missing<-function(dataframe=NULL,tickColumn=NULL,jobColumn=NULL)
{
   attach(dataframe)
   on.exit(detach(dataframe))
   tick <- eval(parse(text=tickColumn))
   job <- eval(parse(text=jobColumn))

   dataframe[1:10,19]<-NA

   missingData<-dataframe[!complete.cases(dataframe),]
   print(missingData)
   n=nrow(missingData)
   m=ncol(missingData)
   names(missingData)<-names(dataframe)
   print(n)
   print(m)
   for(i in 1:n){
      for(j in 1:m){
         if(is.na(missingData[i,j])==T)
          { 
            cat("jobNumber",job[i],"\t")
            cat("tickNumber",tick[i],"\t")
            cat("value for",names(missingData[j]),"is missing","\n")
           }
        }
      } 
}
#dpa.locate.missing(e,"tick","job")



#A function to find a missing row.
dpa.find.missingRow<-function(dataframe=NULL,tickColumn=NULL,jobColumn=NULL)
{
   attach(dataframe)
   on.exit(detach(dataframe))
   tick <- eval(parse(text=tickColumn))
   job <- eval(parse(text=jobColumn))

#To sort the data first. 
   dataframe<-dataframe[order(job,tick),]

#some of the rows are deleted
   dataframe<-dataframe[-3,]
   dataframe<-dataframe[-3,]
   dataframe<-dataframe[-1983,]
   dataframe<-dataframe[-49967,]
   dataframe<-dataframe[-which(dataframe$job==10),]
    
   n=max(job)
   m=max(tick)
   #print(n)
   #print(m)

   l=1
   #To initialize the number of missing jobs and missing ticks.
   mj=0
   mt=0
   
   for(i in 1:n){
      if(is.na(dataframe[l,jobColumn])==T)
        {
         cat("jobNumber",i,"missing","\n")
         }

      if(dataframe[l,jobColumn]!=i)
        { 
         cat("jobNumber",i,"missing","\n")
         mj=mj+1
         }
      else{
         j=1
         while(j<=m){
         k=m*(i-1-mj)+j-mt
         if(dataframe[k,tickColumn]!=j)
         {
         d<-dataframe[k,tickColumn]-dataframe[k-1,tickColumn]
         if(d==2)
         {
         cat("jobNumber",i,"\t")
         cat("tickNumber",j,"missing","\n")
         j=j+1
         mt=mt+1
         }
         if(d>2)
         {
         cat("jobNumber",i,"\t")
         cat("tickNumber",j,"to",j+(d-2),"missing","\n")
         j=j+(d-1)
         mt=mt+(d-1)
         }
        }
        j=j+1
       }
       l=1+m*(i-mj)-mt
      }          
     }
     print(mj)
     print(mt)
}
#dpa.find.missingRow(e,"tick","job")


#######  TO create a dataframe of relations  #######
i=0
relations<-NULL
dpa.relations.addRelations<-function(From_column=NULL,To_column=NULL,Lag_in=NULL,minLag=NULL,maxLag=NULL,Direction=NULL)
{
 relations <- rbind(relations,data.frame(From_column=From_column,To_column=To_column,Lag_in=Lag_in,minLag=minLag,maxLag=maxLag,Direction=Direction,stringsAsFactors = FALSE))
 assign("relations",relations,env=.GlobalEnv)
 #print(relations)
  ## Call the lag function and include lags in the original data ##
  if(as.numeric(maxLag)!=0)
  { 
  if(Lag_in=="From")
  {source <- From_column
   #print(source)
  }
  else{source <- To_column
   }
  if(is.na(match("time_column",ls()))){
    print("Warning! Time column was not set, now set to tick")
    time_column="tick"
    assign("time_column",time_column,env=.GlobalEnv)
  }
  dpa.generate.lag(e,time_column,source,minLag,maxLag)
  e<-dataframe
  assign("e",e,env=.GlobalEnv)
  #print(names(e)) 
  #print("Done creating lags...");
 }
}
#A function to make i global and increase its value 
dpa.incrementValue<-function(i)
{
 i=i+1
 assign("i",i,env=.GlobalEnv)
 print(i)
}  
#To edit and add relations
dpa.relations.editRelations<- function()
{ 
tt <- tktoplevel()
Font1<-tkfont.create(family="times",size=10)
Font2<-tkfont.create(family="times",size=12)
tktitle(tt) <- "  DATAFRAME OF RELATIONS  "

#To submit the entries of a relation
dpa.relations.submitRelations<- function()
{ 
  dpa.incrementValue(i)
  From_column <- source[as.numeric(tclvalue(tcl(FromBox,"getvalue")))+1]
  print(From_column)
  tkgrid(tklabel(tt,text=paste(From_column),font=Font1),row=i+1,column=0,padx=5,pady=5)
  To_column <- destination[as.numeric(tclvalue(tcl(ToBox,"getvalue")))+1]
  print(To_column)
  tkgrid(tklabel(tt,text=paste(To_column),font=Font1),row=i+1,column=1,padx=5,pady=5)
  Lag_in <- lagsIn[as.numeric(tclvalue(tcl(LagIn,"getvalue")))+1]
  print(Lag_in)
  tkgrid(tklabel(tt,text=paste(Lag_in),font=Font1),row=i+1,column=2,padx=5,pady=5)
  minLag <- lags[as.numeric(tclvalue(tcl(LagFrom,"getvalue")))+1]
  print(minLag)
  tkgrid(tklabel(tt,text=paste(minLag),font=Font1),row=i+1,column=3,padx=5,pady=5)
  maxLag <- lags[as.numeric(tclvalue(tcl(LagTo,"getvalue")))+1]
  print(maxLag)
  tkgrid(tklabel(tt,text=paste(maxLag),font=Font1),row=i+1,column=4,padx=5,pady=5)
  Direction <- Directions[as.numeric(tclvalue(tcl(Dir,"getvalue")))+1]
  print(Direction)
  tkgrid(tklabel(tt,text=paste(Direction),font=Font1),row=i+1,column=5,padx=5,pady=5)
  dpa.relations.addRelations(From_column,To_column,Lag_in,minLag,maxLag,Direction)
  #edit(relations)
  dpa.relations.deleteRelations<-function()
 {
  #rows.Remove(i+1)
  tkdelete(tt,"rows",tclvalue(tcl(tt,"active","row")),"1")
 }
  DELETE.but <-tkbutton(tt,text="   delete   ",command=dpa.relations.deleteRelations)
  tkgrid(DELETE.but,row=i+1,column=6,padx=5,pady=0)

}
  SUBMIT.but <-tkbutton(tt,text="   submit   ",command=dpa.relations.submitRelations)
  tkgrid(SUBMIT.but,row=1,column=6,padx=5,pady=0) 

  tkgrid(tklabel(tt,text="From",font=Font2),row=0,column=0,padx=5,pady=5)
  source <- names(e)
  FromBox <- tkwidget(tt,"ComboBox",editable=FALSE,values=source)
  tkconfigure(FromBox)
  tkgrid(FromBox,row=1,column=0,padx=5,pady=5)
        
  tkgrid(tklabel(tt,text="To",font=Font2),row=0,column=1,padx=5,pady=5)
  destination <- names(e)
  ToBox <- tkwidget(tt,"ComboBox",editable=FALSE,values=destination)
  tkconfigure(ToBox) 
  tkgrid(ToBox,row=1,column=1,padx=5,pady=5)

  tkgrid(tklabel(tt,text="Create lags for",font=Font2),row=0,column=2,padx=5,pady=5)
  lagsIn <- c("From","To")
  LagIn <- tkwidget(tt,"ComboBox",editable=FALSE,values=lagsIn)
  tkconfigure(LagIn)  
  tkgrid(LagIn,row=1,column=2,padx=5,pady=5)

  tkgrid(tklabel(tt,text="Minimum lag",font=Font2),row=0,column=3,padx=5,pady=5)
  lags <- c("0","1","2","3","4","5","6","7","8","9","10","11")
  LagFrom <- tkwidget(tt,"ComboBox",editable=TRUE,values=lags)
  tkconfigure(LagFrom)  
  tkgrid(LagFrom,row=1,column=3,padx=5,pady=5)

  tkgrid(tklabel(tt,text="Maximum lag",font=Font2),row=0,column=4,padx=5,pady=5)
  lags <- c("0","1","2","3","4","5","6","7","8","9","10","11")
  LagTo <- tkwidget(tt,"ComboBox",editable=TRUE,values=lags)
  tkconfigure(LagTo)
  tkgrid(LagTo,row=1,column=4,padx=5,pady=5)
   
  tkgrid(tklabel(tt,text="Direction",font=Font2),row=0,column=5,padx=5,pady=5)
  Directions <- c("UniDirectional","BiDirectional")
  Dir <- tkwidget(tt,"ComboBox",editable=FALSE,values=Directions)
  tkconfigure(Dir)
  tkgrid(Dir,row=1,column=5,padx=5,pady=5)
}

################################   Analysis functions ################################

dpa.analysis.options <- function()
{
opt <- tktoplevel()
tktitle(opt) <- " Options "
tkgrid(tklabel(opt,text="Select Time column  "),row=0,column=0,padx=5,pady=5)
columns <- names(e)
Time <- tkwidget(opt,"ComboBox",editable=FALSE,values=columns) 
tkgrid(Time,row=0,column=1,padx=5,pady=5)
tkgrid(tklabel(opt,text="Select an option for model run  "),row=1,column=0,padx=5,pady=5)

#To add radiobuttons for selecting one of the analysis option.
rb1 <- tkradiobutton(opt)
rb2 <- tkradiobutton(opt)
rb3 <- tkradiobutton(opt)
rbValue <- tclVar("time_irrespective")
tkconfigure(rb1,variable=rbValue,value="time_irrespective")
tkconfigure(rb2,variable=rbValue,value="every_timeStep")
tkconfigure(rb3,variable=rbValue,value="time_interval")
tkgrid(tklabel(opt,text="Select time interval per model run  "),row=1,column=0,padx=5,pady=5)

tkgrid(tklabel(opt,text="All data irrespective of time"),rb1)
tkgrid(tklabel(opt,text="For each time step in the data"),rb2)
tkgrid(tklabel(opt,text="For a time interval you want"),rb3)

OnOK <- function()
{
rbVal <- as.character(tclvalue(rbValue))
assign("rbVal",rbVal,env=.GlobalEnv)
if (rbVal=="time_interval")
 {  
 interval <- tclVar()
 entry.interval <-tkentry(opt,width="20",textvariable=interval)
 tkgrid(tklabel(opt,text=" enter time interval "))
 tkgrid(entry.interval)
 }


OnSave <- function()
 {
  time_column <- columns[as.numeric(tclvalue(tcl(Time,"getvalue")))+1]
  print(time_column)
  assign("time_column",time_column,env=.GlobalEnv)
  print(rbVal)

  if(rbVal=="time_interval"){
  intervalVal <- tclvalue(interval)
  print(intervalVal)
  assign("intervalVal",intervalVal,env=.GlobalEnv) 

######################### codes to average data over an interval entered by user  #########################
 #avgData <- NULL
 #numberOfRuns = as.numeric(nrow(e))/as.numeric(intervalVal)
 #for(i in 1:numberOfRuns)
 #{
 # avgRow <- NULL
 # for(j in 1:intervalVal)
 #{
 # rowNumber <- j+(i-1)*as.numeric(intervalVal)
 # if(j==1)
 #{avgRow <- e[rowNumber, ]}
 #else
 #{avgRow<- avgRow + e[rowNumber, ]}
 #}
 # avgData<- rbind(avgData,avgRow)
 #}
 # avgData<-avgData/as.numeric(intervalVal)
 #e<-avgData
 #assign("e",e,env=.GlobalEnv) 
 }
   
 if (rbVal=="time_interval")
 {
  time <- tktoplevel()
  Font1<-tkfont.create(family="times",size=10)
  Font2<-tkfont.create(family="times",size=12)
  tktitle(time) <- "  tick for analysis  "
  tkgrid(tklabel(time,text="enter tick number for analysis",font=Font2),row=0,column=0,padx=5,pady=5)
   
  tkgrid(tklabel(time,text="tick number",font=Font2),row=0,column=1,padx=5,pady=5)
  numTick <- c(1:max(e["tick"]))
  tickNum <- tkwidget(time,"ComboBox",editable=FALSE,values=numTick)
  tkconfigure(tickNum)
  tkgrid(tickNum,row=1,column=1,padx=5,pady=5)
  }
  dpa.analysis.options.ticks.enter<-function()
  { 
  dpa.incrementValue(i)
  NumTick <- numTick[as.numeric(tclvalue(tcl(tickNum,"getvalue")))+1]
  assign("NumTick",NumTick,env=.GlobalEnv)
  print(NumTick)
  tkgrid(tklabel(time,text=paste(NumTick),font=Font1),row=i+1,column=1,padx=5,pady=5)
  }
  enterTick.but <-tkbutton(time,text="   enter   ",command=dpa.analysis.options.ticks.enter)
  tkgrid(enterTick.but,row=1,column=2,padx=5,pady=0) 
  #tkdestroy(opt)
  }
  Save.but <-tkbutton(opt,text="   save   ",command=OnSave)
  tkgrid(Save.but,row=9,column=1)

 OnCancel <- function()
 { 
  tkdestroy(opt)
  }

 Cancel.but <-tkbutton(opt,text="   cancel   ",command=OnCancel)
 tkgrid(Cancel.but,row=9,column=2)
 }

OK.but <- tkbutton(opt,text="OK",command=OnOK)
tkgrid(OK.but)
}
##### Save dpa result to a text file on disk  #####
dpa.analysis.saveDPA <- function(dpaFileName=NULL)
{
 if(is.null(dpaFileName)){
	dpaFileName<-tclvalue(tkgetSaveFile())
}
 assign("dpaFileName",dpaFileName,env=.GlobalEnv)
 if (!nchar(dataFileName)) {
    tkmessageBox(message = "No file was selected.")}
 else {
 sink(file =dpaFileName, append = FALSE, type = c("output", "message"),split=FALSE)
 print("#######################################################")
 print("                         RELATIONS                     ")
 print(relations)
 print("                       ANALYSIS SUMMARY                ")
 print("                       COEFFICIENTS                    ")
 print(sem.results.coefficients)
 sink() 
 }
}


################################# START FUNCTION FOR THE STARTING SCREEN ###########################################

dpa.start <- function(){
#To design the main window of the DPA tool.
dpa <- tktoplevel()
mainFont<-tkfont.create(family="times",size=13)
tktitle(dpa) <- "Dynamic Path Approach"

#reset all the variables, so we can restart the tool.
e<-NULL
listOfTicks<-NULL
relations<-NULL
relevantData<-NULL
row<-NULL
sem.DPA<-NULL
sem.results.coefficients<-NULL
sem.results.parameters<-NULL
sem.results.statistics<-NULL
sem.standardized<-NULL
variables<-NULL


assign("e",e,env=.GlobalEnv)
assign("listOfTicks",listOfTicks,env=.GlobalEnv)
assign("relations",relations,env=.GlobalEnv)
assign("row",row,env=.GlobalEnv)
assign("sem.DPA",sem.DPA,env=.GlobalEnv)
assign("sem.results.coefficients",sem.results.coefficients,env=.GlobalEnv)
assign("sem.results.parameters",sem.results.parameters,env=.GlobalEnv)
assign("sem.results.statistics",sem.results.statistics,env=.GlobalEnv)
assign("sem.standardized",sem.standardized,env=.GlobalEnv)
assign("variables",variables,env=.GlobalEnv)
assign("dpa",dpa,env=.GlobalEnv)


tkgrid(tklabel(dpa,text="  Data  ",font=mainFont),row=0,column=0,padx=5,pady=5)
tkgrid(tklabel(dpa,text="  Relations",font=mainFont),row=0,column=1,padx=5,pady=5)
tkgrid(tklabel(dpa,text="  Analysis",font=mainFont),row=0,column=2,padx=5,pady=5)
tkgrid(tklabel(dpa,text="  Results ",font=mainFont),row=0,column=3,padx=5,pady=5)

#####################    DATA  buttons   ###################################
WorkingDirectory.but <-tkbutton(dpa,text=" Change working directory  ",command=dpa.data.setWorkingDirectory)
tkconfigure(WorkingDirectory.but,width=20,foreground="black",background="white")
tkgrid(WorkingDirectory.but,row=1,column=0,padx=10,pady=10)

LoadData.but <-tkbutton(dpa,text="   Load data  ",command=dpa.data.loadData)
tkconfigure(LoadData.but,width=20,foreground="black",background="white")
tkgrid(LoadData.but,row=2,column=0,padx=10,pady=10)
 
ViewData.but <-tkbutton(dpa,text="   View / Edit data  ",command=dpa.data.viewOrEditData)
tkconfigure(ViewData.but,width=20,foreground="black",background="white")
tkgrid(ViewData.but,row=3,column=0,padx=10,pady=10)

CheckData.but <-tkbutton(dpa,text="  Check data  ",command=dpa.data.checkData)
tkconfigure(CheckData.but,width=20,foreground="black",background="white")
tkgrid(CheckData.but,row=4,column=0,padx=10,pady=10)

SaveData.but <-tkbutton(dpa,text="   Save data  ",command=dpa.data.saveDataToDisk)
tkconfigure(SaveData.but,width=20,foreground="black",background="white")
tkgrid(SaveData.but,row=5,column=0,padx=10,pady=10)

#####################    RELATIONS  buttons   ###############################
LoadRelations.but <-tkbutton(dpa,text="   Load relations  ",command=dpa.relations.loadRelations)
tkconfigure(LoadRelations.but,width=20,foreground="black",background="white")
tkgrid(LoadRelations.but,row=1,column=1,padx=10,pady=10)


EditRelations.but <-tkbutton(dpa,text="   Edit relations   ",command=dpa.relations.editRelations)
tkconfigure(EditRelations.but,width=20,foreground="black",background="white")
tkgrid(EditRelations.but,row=2,column=1,padx=10,pady=10)

SaveRelations.but <-tkbutton(dpa,text="   Save relations  ",command=dpa.relations.saveRelations)
tkconfigure(SaveRelations.but ,width=20,foreground="black",background="white")
tkgrid(SaveRelations.but,row=3,column=1,padx=10,pady=10)

DPA.but <-tkbutton(dpa,text="   Perform DPA  ",command=dpa.analysis.performDPA)
tkconfigure(DPA.but ,width=20,foreground="black",background="white")
tkgrid(DPA.but,row=2,column=2,padx=10,pady=10)

Options.but <-tkbutton(dpa,text="   Options  ",command=dpa.analysis.options)
tkconfigure(Options.but ,width=20,foreground="black",background="white")
tkgrid(Options.but,row=1,column=2,padx=10,pady=10)

saveDPA.but <-tkbutton(dpa,text="   Save DPA O/P  ",command=dpa.analysis.saveDPA)
tkconfigure(saveDPA.but ,width=20,foreground="black",background="white")
tkgrid(saveDPA.but,row=3,column=2,padx=10,pady=10)

############################### RESULTS buttons ################################

graphDir.but <-tkbutton(dpa,text="   Change results directory  ",command=dpa.results.setGraphDir)
tkconfigure(graphDir.but,width=20,foreground="black",background="white")
tkgrid(graphDir.but,row=1,column=3,padx=10,pady=10)

NodesPlot.but <-tkbutton(dpa,text="   Generate nodes plot  ",command=dpa.results.viewNodePlots)
tkconfigure(NodesPlot.but,width=20,foreground="black",background="white")
tkgrid(NodesPlot.but,row=2,column=3,padx=10,pady=10)

RelationsPlot.but <-tkbutton(dpa,text="   Generate relations plot  ",command=dpa.results.viewRelationsPlots)
tkconfigure(RelationsPlot.but,width=20,foreground="black",background="white")
tkgrid(RelationsPlot.but,row=3,column=3,padx=10,pady=10)
}

dpa.exit <- function(){
#devs<-dev.off()
tkdestroy(dpa)
}
