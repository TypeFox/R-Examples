mkSEER<-function(df,seerHome="~/data/SEER",outDir="mrgd",outFile="cancDef",
                 indices = list(c("sex","race"), c("histo3","seqnum"),  "ICD9"),
                 writePops=TRUE,writeRData=TRUE,writeDB=TRUE){
  #   require(dplyr); require(LaF); seerHome="~/data/SEER"
  #   outDir="mrgd";outFile="cancDef";writePops=T;writeRData=TRUE;writeDB=TRUE  # for debugging
  #   indices = list(c("sex","race"), "histo2", "histo3", "ICD9")
  
  # gimic to get rid of unwanted notes in R CMD check
  db=reg=race=sex=age=agerec=year=py=agedx=age19=age86=NULL
#   
#   mkPopsae=function(popsa) {
# #     data(stdUS) #load us 2000 standard population up to 99+ and use it for the extrapolation
#     props=SEERaBomb::stdUS$prop[86:100] #get rid of check Note by directly getting stdUS
#     props=props/sum(props)
#     elders=popsa%>%filter(age86==90)%>%group_by(db,reg,race,sex,year) 
#     #      head(elders,2)
#     #      head(popsa,2)
#     grow=function(x) data.frame(x,age=85.5:99.5,PY=x$py*props) 
#     elders=elders%>%do_(~grow(.))%>%select(-age86,-py) # no note on check about .
# #     elders=elders%>%do(grow(.))%>%select(-age86,-py) # => Note on check about .
#     nms=names(elders)
#     names(elders)[which(nms=="PY")]<-c("py")
#     nelders=popsa%>%filter(age86<90)
#     nms=names(nelders)
#     names(nelders)[which(nms=="age86")]<-c("age")
#     elders=elders[,names(nelders)] # reorder columns to match
#     rbind(nelders,elders) # this now becomes popsa extended (i.e. popsae) 
#   }

mkPopsae=function(popsa) { # this replaces the version above it
  rates=SEERaBomb::nvsr01[86:100,] #get rid of check Note by directly getting nvsr
  props=cumprod(1-rates[,c("pm","pf")])
  tots=cumsum(props)["100",]
  for (i in 1:dim(props)[2]) props[,i]=props[,i]/tots[,i]
  elders=popsa%>%filter(age86==90)%>%group_by(db,reg,race,sex,year) 
  grow=function(x) {if (x$sex=="male") y=data.frame(x,age=85.5:99.5,PY=x$py*props$pm) else
    y=data.frame(x,age=85.5:99.5,PY=x$py*props$pf)
    y
  }
  elders=elders%>%do_(~grow(.))%>%select(-age86,-py) # no note on check about .
  nms=names(elders)
  names(elders)[which(nms=="PY")]<-c("py")
  nelders=popsa%>%filter(age86<90)
  nms=names(nelders)
  names(nelders)[which(nms=="age86")]<-c("age")
  elders=elders[,names(nelders)] # reorder columns to match
  rbind(nelders,elders) # this now becomes popsa extended (i.e. popsae) 
}

#  mkPopsae(popsa) # for testing

  seerHome=path.expand(seerHome)
  outD=file.path(seerHome,outDir) 
  outF=file.path(seerHome,outDir,paste0(outFile,".RData")) 
  outDB=file.path(seerHome,outDir,paste0(outFile,".db")) 
  dir.create(outD,showWarnings = FALSE)   # OK if already exists, so suppress warning
  
  if(writePops|writeDB) {
    cat("Making population file data.tables\n")
    dirs=list.dirs(file.path(seerHome,"populations"))
    dirs=dirs[grep("yr",dirs)]
    dirs=dirs[-grep("yr2005",dirs)]
    
    colTypes=c("integer","string",rep("integer",5),"double")    # double, integer, categorical and string
    colWidths=c(4,7,2,1,1,1,2, 10) 
    # colNames1 = c('year','X0','reg','race','origin','sex', 'age19','py')
    colNames2 = c('year','X0','reg','race','origin','sex', 'age86','py')
    
    ptm <- proc.time()
    # popgaL=vector(mode="list",3)
    popsaL=vector(mode="list",3)
    ii=1
    for (i in dirs) {
#       f=file.path(i,"19agegroups.txt") 
#       laf<-laf_open_fwf(f,column_types=colTypes,column_widths=colWidths,column_names = colNames1)
#       popgaL[[ii]]=tbl_df(laf[,colNames1[-c(2,5)]] )
      
      f=file.path(i,"singleages.txt") 
      laf<-laf_open_fwf(f,column_types=colTypes,column_widths=colWidths,column_names = colNames2)
      popsaL[[ii]]=tbl_df(laf[,colNames2[-c(2,5)]] )
      
      if (grepl("plus",i)) {#5/21/14 fixed systematic lower 73 incidence due to 73 PY being also in 92
        # popgaL[[ii]]=popgaL[[ii]]%>%filter(reg%in%c(29,31,35,37)) 
        popsaL[[ii]]=popsaL[[ii]]%>%filter(reg%in%c(29,31,35,37)) 
        cat("Removing SEER 9 person years from:\n",i,"\nbefore pooling into one file.\n")
      }
      ii=ii+1
    }
    # popga=rbind_all(popgaL)
    # popsa=rbind_all(popsaL)
    popsa=bind_rows(popsaL)
    # popga$reg=as.integer(popga$reg+1500)
    popsa$reg=as.integer(popsa$reg+1500)
    
#     popga=popga%>%
#       mutate(race=cut(race,labels=c("white","black","other"),breaks=c(1,2,3,100), right=F))  %>%
#       mutate(db=cut(reg,labels=c("73","92","00"),breaks=c(1500,1528,1540,1550), right=F))  %>%
#       mutate(reg=mapRegs(reg)) %>%
#       mutate(age19=c(0.5,3,seq(7.5,82.5,5),90)[age19+1]) %>%
#       mutate(sex=factor(sex,labels=c("male","female"))) %>%
#       group_by(db,reg,race,sex,age19,year) %>%
#       #       summarise(py=sum(py))%>%   #summing here  over counties and hispanic origin or not, which are in popga as rows but not as columns
#       #       group_by(add=F) # clear grouping
#       summarise(py=sum(py))
#     popga=as.data.frame(popga) # clears everything, down to the data.frame
    
    #     class(popga)
    
    popsa=popsa%>%    # note: here age groups indices 0-18 are replaced by actual ages 1-85, so no need to remap
      mutate(race=cut(race,labels=c("white","black","other"),breaks=c(1,2,3,100), right=F))  %>%
      mutate(db=cut(reg,labels=c("73","92","00"),breaks=c(1500,1528,1540,1550), right=F))  %>%
      mutate(reg=mapRegs(reg)) %>%    
      mutate(age86=age86+0.5) %>%    
      mutate(sex=factor(sex,labels=c("male","female"))) %>%
      group_by(db,reg,race,sex,age86,year) %>%
      #       summarise(py=sum(py))%>%   #summing here  over counties and hispanic origin or not, which are in popga as rows but not as columns
      #       group_by(add=F) # clear grouping
      summarise(py=sum(py))
    popsa=as.data.frame(popsa) # clears everything down to the data.frame
    popsa[popsa$age86==85.5,"age86"]=90 # touch up to better guess of average age of >85 group
    popsae=mkPopsae(popsa)  # extend popsa to fill in ages 85-99 based on the US STD 2000 population  
    
    delT=proc.time() - ptm  
    cat("The population files of SEER were processed in ",delT[3]," seconds.\n",sep="")
    
  } #if writePops or writeDB
  
  
  if(writeRData|writeDB) {
    dirs=list.dirs(file.path(seerHome,"incidence"))
    dirs=dirs[grep("yr",dirs)]
    dirs=dirs[-grep("yr2005",dirs)]
    
    y=df[which(df$names!=" "),"names"] 
    cat("Cancer Data:\nThe following fields will be written:\n");  print(y)
    cancers=c('breast','digothr','malegen','femgen','other','respir','colrect','lymyleuk','urinary') 
    files=paste0(cancers,".txt")
    DFL=vector(mode="list",length(dirs)*length(files))
    ptm <- proc.time()
    ii=1
    for (i in dirs) 
      for (j in files) {
        f=file.path(i,toupper(j))
        print(f)
        laf<-laf_open_fwf(f,column_types=df$type, column_widths=df$width)
        DFL[[ii]]=tbl_df(laf[,which(df$names!=" ")])
        ii=ii+1
      }
    print("using bind_rows() to make DF canc")
    # canc=rbind_all(DFL)
    canc=bind_rows(DFL)
    colnames(canc)<-y
    
    canc=canc%>%
      filter(agedx<200)%>%
      mutate(race=cut(race,labels=c("white","black","other"),breaks=c(1,2,3,100), right=F))  %>%
      mutate(db=cut(reg,labels=c("73","92","00"),breaks=c(1500,1528,1540,1550), right=F))  %>%
      mutate(reg=mapRegs(reg)) %>%
      # mutate(age19=c(0.5,3,seq(7.5,82.5,5),90)[agerec+1]) %>%
      mutate(age86=as.numeric(as.character(cut(agedx,c(0:85,150),right=F,labels=c(0.5:85,90)))))%>%
      #       select(-agerec)%>%
      mutate(sex=factor(sex,labels=c("male","female")))
    canc=mapCancs(canc)
    canc=mapTrts(canc)
    delT=proc.time() - ptm  
    cat("Cancer files were processed in ",delT[3]," seconds.\n")
  } #if writeRData or writeDB
  
  if (writeDB) {
    print("Deleting old version of this SQLite database file.")
    unlink(outDB)
    print("Creating new SQLite database file.")
    ptm <- proc.time()
    #     m=dbDriver("SQLite")
    #     con <- dbConnect(m, dbname = outDB)
    my_db <- src_sqlite(outDB, create = T)
    #     copy_to(my_db,DF, name="canc",temporary = FALSE, indexes = indices,overwrite=TRUE) 
    copy_to(my_db,canc, temporary = FALSE, indexes = indices,overwrite=TRUE) 
    # copy_to(my_db,popga, temporary = FALSE,overwrite=TRUE)
    copy_to(my_db,popsa, temporary = FALSE,overwrite=TRUE)
    copy_to(my_db,popsae, temporary = FALSE,overwrite=TRUE)
    #     dbWriteTable(con, "canc", DF,overwrite=TRUE)
    #     dbWriteTable(con, "popga", popga,overwrite=TRUE)
    #     dbWriteTable(con, "popsa", popsa,overwrite=TRUE)
    delT=proc.time() - ptm  
    cat("\ndata.tables were written to ",outDB," in ",delT[3]," seconds.\n")
    cat("tables in this SQLite file are:\n")
    #     print(dbListTables(con))   
    #     dbDisconnect(con)
  }
  
  if (writeRData) {
    print("save()-ing DF canc to disk")
    canc=tbl_df(canc)
    save(canc,file=outF)
    cat("Cancer data has been written to: ",outF,"\n")
  }
  
  if (writePops) {
    # popga=tbl_df(popga)
    popsa=tbl_df(popsa)
    popsae=tbl_df(popsae)
    # save(popga,file=file.path(outD,"popga.RData"))  
    save(popsa,file=file.path(outD,"popsa.RData"))  
    save(popsae,file=file.path(outD,"popsae.RData"))  
  } #only if writePops
  #if(.Platform$OS.type=="unix") {
  s=dir(outD,full.names=T)
  d=file.info(s)[,c("size","mtime")] 
  d$size=paste0(round(d$size/1e6)," M")
  print(d) 
  #}
}
