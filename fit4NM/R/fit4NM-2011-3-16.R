fit4NM<-function()
{  options(guiToolkit="RGtk2")
   require(gWidgets)
   require(tcltk)
   require(tkrplot)
   require(RGtk2)
   require(gWidgetsRGtk2)
   require(cairoDevice)   
   
   dir.create("c:/fit4NM",showWarnings=FALSE)
   if(length(dir("c:/fit4NM",all.files=TRUE))>2)
   {   load("c:/fit4NM/.Rdata",.GlobalEnv)
   }
   
   PRED.list.7<<-c("PRED","IPRED", "NPRED" ,"PREDI","CPRED","CPREDI","EPRED")
   RES.list.7<<-c("RES","IRES","WRES","IWRES","NWRES","NRES","RESI","WRESI","CRES",
               "CWRES","CRESI","CWRESI","ERES","EWRES","NPDE")
   PRED.list.6<<-c("PRED","IPRE")
   RES.list.6<<-c("RES","WRES","IWRE","IRES")
 
   ############################################################################# 
   # Configuration
   #############################################################################

   #### Notes before configuration #############################################
   ConfigNotes<-function(h,...)
   {  gmessage("*** Notes before configuration ***
                \nNONMEM path
                \n      Copy \'c:\\NONMEM folder\\run or util\\nmfe6.bat or nmfe7.bat\' to other folder.
                \n      To run NONMEM in R, add REM to the lines containing \'del\' in \'nmfe6.bat\' and \'nmfe7.bat\'.
                \n      Create \'c:\\NONMEM folder\\fit4NM\\' and copy  \'these modified nmfe6.bat or nmfe7.bat' to 'c:\\NONMEM folder\\fit4NM'.
                \n      Use \'these modified nmfe7.bat and nmfe6.bat\' as default and alternative NONMEM pathways.
                \nNONMEM help
                \n      Find \'c:\\NONMEM folder\\html\\index.htm\\' for both of the default and alternative NONMEM.",cont=TRUE,width=600)   
   }
   
   #### Set NONMEM path-default ################################################
   NMpath1<-function(h,...)
   {  Default.NMpath<<-gfile(text="Choose default NONMEM",type="open")
   }

   #### Set NONMEM help-default ################################################  
   NMpathHelp1<-function(h,...)
   {  Default.Helppath<<-gfile(text="Choose default NONMEM help - index.htm",type="open")
   }
    
   #### Set NONMEM path-alternative ############################################
   NMpath2<-function(h,...)
   {  Alternative.NMpath<<-gfile(text="Choose alternative NONMEM",type="open")
   }
    
   #### Set NONMEM help-alternative ############################################
   NMpathHelp2<-function(h,...)
   {  Alternative.Helppath<<-gfile(text="Choose alternative NONMEM help - index.htm",type="open")
   }
     
   #### Set external editor ####################################################
   Editorpath<-function(h,...)
   {  Editor.path<<-gfile(text="Choose external editor",type="open")
   }
    
   #### Save configuration #####################################################
   saveConfig<-function(h,...)
   {  dir.create("c:/fit4NM",showWarnings=FALSE)
      dir.list<-ls(pos=1)
      flag.D.H<-sum(dir.list=="Default.Helppath")
      flag.D.N<-sum(dir.list=="Default.NMpath")
      flag.A.H<-sum(dir.list=="Alternative.Helppath")
      flag.A.N<-sum(dir.list=="Alternative.NMpath")
      flag.E<-sum(dir.list=="Editor.path")
      save.list<-c(ifelse(flag.D.H==1,"Default.Helppath",NA),ifelse(flag.A.H==1,"Alternative.Helppath",NA),
                   ifelse(flag.D.N==1,"Default.NMpath",NA),ifelse(flag.A.N==1,"Alternative.NMpath",NA),
                   ifelse(flag.E==1,"Editor.path",NA))
      save.list<-save.list[!is.na(save.list)]
      save(list = save.list, file = "c:/fit4NM/.Rdata")
   }

   ############################################################################# 
   # Data
   #############################################################################

   #### Data manipulation ###################################################### 
   ###### Calculate elapsed time ###############################################  
   CalcTime<-function(h,...)
   {  openStandard<-function(h,...)
      {  standard.file<<-gfile(text="Choose standard time file",type="open")
         temp<-strsplit(standard.file,split="\\\\")[[1]]
         file.id<-temp[length(temp)]
         dir.name<-strsplit(standard.file,split=file.id)[[1]]
         setwd(dir.name)
         svalue(standard.t)<-standard.file
         ref.data<-read.csv(standard.file)
         ind.list<-colnames(ref.data)
         tmp<-gframe("",cont=BBgroup,height=50)
         button4<-gbutton("Select indicator",width=20,height=10)
         add(tmp,button4)
         indc.t<<-gdroplist(c("NONE",ind.list),width=200,height=200)
         add(tmp,indc.t)        
         tmp<-gframe("",cont=BBgroup)
         button2<-gbutton("Select data folder #ID, DATE (yyyy-mm-dd), TIME (h:mm:ss),... after data split by ID",handler=openIndiv1,width=20,height=10)
         add(tmp,button2)
         indiv.t<<-gdroplist(c(" "),width=200,height=200)
         add(tmp,indiv.t)
         tmp<-gframe("",cont=BBgroup,height=50)
         button4<-gbutton("Select time unit",width=20,height=10)
         add(tmp,button4)
         time.t<<-gdroplist(c("secs","mins","hours","days","weeks"),width=200,height=200)
         add(tmp,time.t)
         tmp<-gframe("",cont=BBgroup)
         button3<-gbutton("Calculate and save as csv",handler=ElapseTime,width=20,height=10)
         add(tmp,button3)         
      }
      
      openIndiv1<-function(h,...)
      {  indiv.dir<<-gfile(text="Choose folder for individual data",type="selectdir")
         indiv.file<-dir(indiv.dir)
         for(i in 1:length(indiv.file))
         {  indiv.t[i]<-indiv.file[i]
         }
      }
  
      ElapseTime<-function(h,...)
      {       
         setTIME<-svalue(standard.t)
         setDIR<-indiv.dir 
         UNIT<-svalue(time.t)
         ind.VAR<-svalue(indc.t)
         New.Time<-CalcTIME0(setTIME,setDIR,UNIT,ind.VAR)
         temp<-gfile(text="Save calculated elapsed time as csv",
           type="save") #,filter=list("csv files"=list(patterns=c("*.csv"))))
         file.name<-strsplit(temp,split="\\.")[[1]][1]
         write.csv(New.Time,paste(file.name,".csv",sep=""),row.names=F,quote=F)
         dispose(timewin1)  
      }

      CalcTIME0<-function(setTIME,setDIR,UNIT,ind.VAR)
      {  file.list<-dir(setDIR)
         setTIME<-read.csv(standard.file)
         tot.data<-NULL
         for(i in 1:length(file.list))
         {  data.tt<-read.csv(paste(setDIR,"\\",file.list[i],sep=""))
            if(ind.VAR!="NONE")
            {  ind.list<-names(table(data.tt[,ind.VAR]))
               for(j in 1:length(ind.list))
               {  sel.ind<-which(as.character(data.tt[,ind.VAR])==ind.list[j])
                  data.t<-data.tt[sel.ind,]
                  if( sum(setTIME$X.ID==data.t$X.ID[1] )!=0)
                  {  orig.colname<-colnames(data.t)
                     colnames(data.t)<-toupper(orig.colname)
                     DATETIME<-paste(data.t$DATE,data.t$TIME,sep=" ")
                     data.temp<-strptime(DATETIME,"%Y-%m-%d %H:%M:%S")
                     standard<-paste(as.character(setTIME$DATE0[setTIME$X.ID==data.t$X.ID[1] & setTIME[,ind.VAR]==ind.list[j]]),
                             as.character(setTIME$TIME0[setTIME$X.ID==data.t$X.ID[1]& setTIME[,ind.VAR]==ind.list[j]]),sep=" ")
                     TIME<-as.numeric(difftime(strptime(data.temp,"%Y-%m-%d %H:%M:%S"),
                             strptime(standard,"%Y-%m-%d %H:%M:%S"), units=UNIT))
                     data.t$TIME<-TIME
                     data.t<-data.t[,-which(toupper(orig.colname)=="DATE")]
                     colnames(data.t)<-orig.colname[-which(toupper(orig.colname)=="DATE")]
                     tot.data<-rbind(tot.data,data.t)
                  }
               }        
            } else
            {  data.t<-data.tt
               if( sum(setTIME$X.ID==data.t$X.ID[1])!=0)
               {  orig.colname<- colnames(data.t)
                  colnames(data.t)<-toupper(orig.colname)
                  DATETIME<-paste(data.t$DATE,data.t$TIME,sep=" ")
                  data.temp<-strptime(DATETIME,"%Y-%m-%d %H:%M:%S")
                  standard<-paste(as.character(setTIME$DATE0[setTIME$X.ID==data.t$X.ID[1]]),
                          as.character(setTIME$TIME0[setTIME$X.ID==data.t$X.ID[1]]),sep=" ")
                  TIME<-as.numeric(difftime(strptime(data.temp,"%Y-%m-%d %H:%M:%S"),
                          strptime(standard,"%Y-%m-%d %H:%M:%S"), units=UNIT))
                  data.t$TIME<-TIME
                  data.t<-data.t[,-which(toupper(orig.colname)=="DATE")]
                  colnames(data.t)<-orig.colname[-which(toupper(orig.colname)=="DATE")]
                  tot.data<-rbind(tot.data,data.t)
               }   
            }   
         }
         ID.sort<-sort.list(tot.data$X.ID)
         tot.data<-tot.data[ID.sort,]
         temp<-colnames(tot.data)         
         temp[temp=="X.ID"]<-"#ID"
         colnames(tot.data)<-temp
         return(tot.data)
      }

      timewin1<<-gwindow("Caclulate elapsed time")
      BBgroup<<-ggroup(cont=timewin1,horizontal=FALSE)
      Bgroup1<<-ggroup(cont=BBgroup,horizontal=TRUE)
      tmp<-gframe("",cont=BBgroup)
      standard.t<-gedit(" ",width=50)
      button1<-gbutton("Open reference time file #ID, DATE0 (yyyy-mm-dd), TIME0 (h:mm:ss)",handler=openStandard)
      add(tmp,button1)
      add(tmp,standard.t)
   }

   ###### Join data ############################################################  
   DataJoinhandler<-function(h,...)
   {  DataJoin<-function(h,...)
      {  tmp<-DataJoin.totDATA
         dir.name<<-gfile(text="Data file folder",type="selectdir")
         tclvalue(DataJoin.Name)<<-dir.name
         setwd(dir.name)
         indiv.file<-dir(dir.name)
         for(i in 1:length(indiv.file))
         {  tmp[[i]]<-read.csv(indiv.file[i],na.strings=".")
         }
         DataJoin.totDATA<<-tmp
         DJk<<-length(indiv.file)
      }

      JoinData<-function(h,...)
      {  merge.data<-DataJoin.totDATA[[1]]
            if(DJk>=2)
            {  for(i in 2:DJk)
               {  #merge.data<-merge(merge.data,DataJoin.totDATA[[i]],by=by.var,suffixes=c("",""))
                  merge.data<-merge(merge.data,DataJoin.totDATA[[i]],all=T)
               }
               temp<-colnames(merge.data)
               temp[temp=="X.ID"]<-"#ID"
               colnames(merge.data)<-temp
               file.new<-gfile(text="Save as csv",type="save")
               write.csv(merge.data,paste(file.new,".csv",sep=""),quote=FALSE,row.names=FALSE)
            } else
            {  warnings()
            }
            tkdestroy(Toptt)    
      }

      DataJoin.totDATA<<-list()
      Toptt<<-tktoplevel()
      tkwm.title(Toptt,"Data join")
      tt<-tkframe(Toptt)
      ttg<-tkframe(tt)
      tkgrid(ttg)      	
      OK.but2 <-tkbutton(tt,text="Select csv file folder",command=DataJoin)
      Join.but<-tkbutton(tt,text="Join",command=JoinData)
      tkgrid(tt)
      DataJoin.Name<-tclVar("")
      DataJoinName<-tkentry(tt,width="60",textvariable=DataJoin.Name) 
      tkgrid(tklabel(tt,text=""),
                    OK.but2,DataJoinName,tklabel(tt,text=""),tklabel(tt,text=""),Join.but)
      tkgrid(tt)	
   }
      
   ###### Split data ###########################################################  

   DataSplit<-function(h,...)
   {  splitD<-function(h,...)
      {  splitDIR<-gfile("Choose folder",type="selectdir")
         setwd(splitDIR)
         cond1<-svalue(VarList.g1)
         cond2<-svalue(VarList.g2)
         temp<- strsplit(whole.file,"\\\\")[[1]]
         temp<-temp[length(temp)]
         header<-strsplit(temp,"\\.")[[1]][1]
         if(cond2=="NONE")
         {  cond<-paste(cond1,as.character(whole.data[,cond1]),sep="_")
         } else
         {  cond<-paste(cond1,as.character(whole.data[,cond1]),cond2,as.character(whole.data[,cond2]),sep="_")
         }  
         sel.id<-names(table(cond))
     
         for(i in sel.id)
         {  select.id<-which(cond==i)
            select.data<-whole.data[select.id,]
            sel.filename<-paste(header,"_",i,".csv",sep="")
            temp<-colnames(select.data)
            temp[temp=="X.ID"]<-"#ID"
            colnames(select.data)<-temp
            write.csv(select.data,sel.filename,row.names=F)
         }  
         dispose(gDS)
      }

      whole.file<-gfile(text="Open control file",type="open")
      whole.data<-read.csv(whole.file)
      Var.Name<-colnames(whole.data)
      VarList.g1<-gdroplist(Var.Name)
      VarList.g2<-gdroplist(c("NONE",Var.Name))
      gDS<-gwindow("Data split",width=100)
      ggDS<-ggroup(horizontal=FALSE,cont=gDS)
      tmp<-gframe("Level 1",cont=ggDS)
      add(tmp,VarList.g1)
      tmp<-gframe("Level 2",cont=ggDS)
      add(tmp,VarList.g2)
      Button1<-gbutton("OK",handler=splitD)
      tmp<-gframe("Split",cont=ggDS)
      add(tmp,Button1)
   }
      
   ###### Convert data #########################################################  
   ######## colume to line #####################################################  
   ColtoLine<-function(h,...)
   {  convert.CL<-function(h,...)
      {  CL.name.list<-NULL
         if(sum(toupper(colnames(conv.data))=="MDV")!=0)
            conv.data<-conv.data[conv.data$MDV!=1,]        
         for(i in 1:3)
            CL.name.list<-c(CL.name.list, svalue(CL.list[[i]]))
         conv.temp.data<-conv.data[,CL.name.list]
         colnames(conv.temp.data)<-c("ID","X","Y")
         time<-sort(unique(as.character(conv.temp.data$X), replace=T))
         ntime<- length(time)
         subj<-unique(as.character(conv.temp.data$ID), replace=T)
         nsubj<-length(subj)
         prep <- matrix(nrow=ntime, ncol=1)
         colnames(prep) <- "X"
         prep[,1] <-time
         TIME.t<-prep
         prep.i<-prep
         for(i in subj)
         {  d<-conv.temp.data[conv.temp.data$ID == i & !is.na(conv.temp.data$Y),]
            conc<-data.frame(X=d$X,Y=d$Y)
            temp<-merge(TIME.t, conc, by="X",all=T)
            prep.i<-cbind(prep.i,temp$Y)
         }
         colnames(prep.i)<- c("X",subj)
         CL.final<<-prep.i
         tmp<-gframe("Save as csv",cont=Bgroup1)
         gbutton("Save",handler=save.CL,cont=tmp)
      }   

      save.CL<-function(h,...)
      {  save.name<-gfile(text="save as csv",type="save")
         write.table(CL.final, paste(save.name,".csv",sep=""), sep=",", row.names=F,na=".")
         dispose(CtoL.win)
      } 

      select.CL<-function(h,...)
      {  sum.name<-gfile(text="csv file with column data (eg: ID, TIME, DV)",type="open")
         svalue(file.CL)<-sum.name
         conv.data<<-read.csv(sum.name,na.strings=".")
         CL.g<-list()
         CL.list<<-list()
         CLparam.input<-c("ID","X","Y")
         for(i in 1:3)
         {  CL.g[[i]]<-ggroup(cont=Bgroup1)
            glabel(CLparam.input[i],cont=CL.g[[i]])
            temp<- gdroplist(colnames(conv.data),cont=CL.g[[i]])
            CL.list[[i]]<<-temp
         }
         tmp<-gframe("Convert",cont=Bgroup1)
         gbutton("Convert",handler=convert.CL,cont=tmp)
      }  
 
      CtoL.win<<-gwindow("Convert data : column to line")
      BBgroup<-ggroup(cont=CtoL.win,horizontal=TRUE)
      Bgroup1<-ggroup(cont=BBgroup, horizontal=FALSE)
      tmp<-gframe("Open csv file with column data (eg: ID, TIME, DV)",cont=Bgroup1)
      file.CL<<-gedit(" ",cont=tmp)
      glabel("(missing=\".\")",cont=tmp)        
      gbutton("Open",handler=select.CL,cont=tmp)
   }
      
   ######## line to colume #####################################################  
   LinetoCol<-function(h,...)
   {  convert.CL<-function(h,...)
      {  ID.n<-ncol(conv.data)-1
         tot.data<-NULL
         ID.list<-as.character(conv.data[1,-1])
         name.by.ID<-as.character(conv.data[,1])[-1]
         for(i in 1:ID.n)
         {  temp<-cbind(rep(ID.list[i],length(name.by.ID)),name.by.ID,conv.data[-1,i+1])
            tot.data<-rbind(tot.data,temp)
         }
         colnames(tot.data)<-c("ID",as.character(conv.data[1,1])," ")
         LC.final<<-tot.data
         tmp<-gframe("Save as csv",cont=Bgroup1)
         gbutton("Save",handler=save.LC,cont=tmp)
      }   

      save.LC<-function(h,...)
      {  save.name<-gfile(text="save as csv",type="save")
         write.table(LC.final, paste(save.name,".csv",sep=""), sep=",", row.names=F,na=".")
         dispose(CtoL.win)
      } 

      select.CL<-function(h,...)
      {  sum.name<-gfile(text="colume to line file",type="open")
         svalue(file.CL)<-sum.name
         conv.data<<-read.csv(sum.name,na.strings=".",header=F)
         tmp<-gframe("Convert",cont=Bgroup1)
         gbutton("Convert",handler=convert.CL,cont=tmp)
      }  

      CtoL.win<<-gwindow("Convert data : line to column")
      BBgroup<-ggroup(cont=CtoL.win,horizontal=TRUE)
      Bgroup1<-ggroup(cont=BBgroup, horizontal=FALSE)
      tmp<-gframe("Open csv file with line data",cont=Bgroup1)
      file.CL<<-gedit(" ",cont=tmp)
      glabel("(missing=\".\")",cont=tmp)        
      gbutton("Open",handler=select.CL,cont=tmp)
   }
   ###### Biosignal data selection##############################################  
   #### Notes before selection #############################################
   BiosignalNotes<-function(h,...)
   {  gmessage("*** Notes before biosignal data preparation ***
                \nPlease use elapsed time instead of actual time.",cont=TRUE,width=600)   
   }   
   ######## Batch process ######################################################  
   BDS.batch<-function(h,...)
   {
      openBDS<-function(h,...)
      {   BDS.dir<<-gfile(text="Open biosignal data file folder",type="selectdir")
          svalue(control.t)<-BDS.dir
      }
      
      BDScalc<-function(h,...)
      {  N.BDS<-length(indiv.file)
         saveDIR<-gfile(text="Select folder for saving data",type="selectdir")
         win<-gwindow(paste("Biosignal data : Selection : Barch progressing : N=",N.BDS,sep=""),width=300,height=50)
         BDS.progress<-gslider(from=0,to=N.BDS,by=1,value=0,cont=win)               
         for(i in 1:length(indiv.file))
         {   svalue(BDS.progress)<-i
             temp.tot<-read.csv(indiv.file[i],na.string=".")
             var.name<-colnames(temp.tot)
             X.id<-which(var.name==svalue(BDS.list[[1]]))
             Y.id<-which(var.name==svalue(BDS.list[[2]]))
             X<-temp.tot[,X.id]
             Y<-temp.tot[,Y.id]
             N<-length(Y)
             diff.Y<-(Y[-1]-Y[-length(Y)])/(X[-1]-X[-length(X)])
             cut.rate<-as.numeric(svalue(sel.rate))/100
             cut.n<-round(N*cut.rate)

             sel.id<-sort(sort.list(abs(diff.Y),decreasing=T)[-c(1:cut.n)])
             select.flag<-rep(FALSE,length(X))
             select.flag[c(1,sel.id+1)]<-TRUE
             select.data<-cbind(temp.tot,select.flag)
             colnames(select.data)<-c(colnames(temp.tot),paste("selected",svalue(BDS.list[[2]]),sep="."))
             file.name<-paste(strsplit(indiv.file[i],split="\\.")[[1]][1],"-",
                              svalue(BDS.list[[2]]),".csv",sep="")
             full.file<-paste(saveDIR, file.name,sep="\\")                
             write.csv(select.data,file.name)           
          }  
          dispose(win)
#          dispose(BDSwin)  
      }      
      
      TDselect<-function(h,...)
      {  setwd(BDS.dir)
         indiv.file<<-dir(BDS.dir)
         temp<-read.csv(indiv.file[1],na.string=".")
         var.name.temp<-colnames(temp)
         BDSparam.input<-c("TIME  ","DV ")
            
         BDS.g<-list()
         BDS.list<-list()
         for(i in 1:2)
         {  BDS.g[[i]]<-ggroup(cont=BBgroup)
            temp<- gdroplist(var.name.temp,cont=BDS.g[[i]])
            glabel(BDSparam.input[i],cont=BDS.g[[i]])
             BDS.list[[i]]<-temp
         }
         BDS.list<<-BDS.list           
         
         Button2<-gbutton("Select data",handler= BDScalc)
         tmp=gframe("",container=BBgroup)
         add(tmp,Button2)
         tmp<-gframe("",container=BBgroup)
         glabel("Difference[t] = (DV[t] - DV[t-1])/(TIME[t]-TIME[t-1]),discard DV[t] if abs(difference[t]) ranked in descending order > cutoff percentile",
                   container=tmp)       
        
#         add(BigGroup,ggraphics())
         
      }  
    
      BDSwin<<-gwindow("Biosignal data selection : Batch processing")
      BigGroup<<-ggroup(cont=BDSwin,horizontal=TRUE)
      BBgroup<<-ggroup(cont=BigGroup,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open data folder",handler=openBDS)
      add(tmp,button1)
      add(tmp,control.t)
      sel.rate<<-gedit(" ",width=50)   
      tmp=gframe("Cutoff percentile of difference",container=BBgroup)
      add(tmp,sel.rate)                   
      Button<-gbutton("Select TIME and DV",handler= TDselect)
      tmp=gframe("",container=BBgroup)
      add(tmp,Button)
   }  

   ######## Individual process##################################################  
   BDS.indiv<-function(h,...)
   {
      BDSsave<-function(h,...)
      {  temp.tot<-read.csv(file.name,na.string=".")    
         select.data<-cbind(temp.tot,sel.FINAL)
         colnames(select.data)<-c(colnames(temp.tot),paste("selected",svalue(BDS.list[[2]]),sep="."))
         write.csv(select.data,paste(gfile(text="Save case deletion raw data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)         
      }
              
      BDS.indiv.GUI<-function(h,...)
      {  file.name<<-svalue(FileList)
         temp.tot<-read.csv(file.name,na.string=".")
         var.name.temp<-c("NONE",colnames(temp.tot))
         BDSparam.input<-c("TIME  ","DV ","Previous selection")
      
         BDSINDwin<<-gwindow(paste("Biosignal data : Selection : Individual processing",
         file.name,sep="---"))
         BigGroup<<-ggroup(cont=BDSINDwin,horizontal=TRUE)
         BBgroup<<-ggroup(cont=BigGroup,horizontal=FALSE)
         sel.rate<<-gedit(" ",width=20)   
         tmp=gframe("Cutoff percentile of difference",container=BBgroup)
         add(tmp,sel.rate)  
         BDS.g<-list()
         BDS.list<-list()
         for(i in 1:3)
         {  BDS.g[[i]]<-ggroup(cont=BBgroup)
            temp<- gdroplist(var.name.temp,cont=BDS.g[[i]])
            glabel(BDSparam.input[i],cont=BDS.g[[i]])
             BDS.list[[i]]<-temp
         }
         BDS.list<<-BDS.list           
#         checkB<-gcheckbox("Use  previous selection")
#         tmp=gframe("",container=BBgroup)
#         add(tmp,checkB)
         Button2<-gbutton("Select data",handler= BDScalc)
         tmp<-gframe("",container=BBgroup)
         add(tmp,Button2)
         Button3<-gbutton("Save data",handler= BDSsave)
         tmp<-gframe("",container=BBgroup)
         add(tmp,Button3)
         tmp<-gframe("",container=BBgroup)
         glabel("Difference[t] = (DV[t] - DV[t-1])/(TIME[t]-TIME[t-1]),discard DV[t] if abs(difference[t]) ranked in descending order > cutoff percentile",
                   container=tmp)       
         
         flag.start<<-TRUE
         iter.N<<-0
    }
 
      openBDS<-function(h,...)
      {   BDS.dir<<-gfile(text="Open biosignal data file folder",type="selectdir")
          svalue(control.t)<-BDS.dir
          setwd(BDS.dir)
          indiv.file<<-dir(BDS.dir)
          tmp<-gframe("Files",container=BBgroup)
          Files<-matrix(indiv.file,ncol=1)
          colnames(Files)<-"Files"
          FileList<<-gtable(Files,multiple=T,handler=BDS.indiv.GUI)
          size(FileList)<-c(500,200)
          add(tmp,FileList)
      }
      
      BDScalc<-function(h,...)
      {   iter.N<<-iter.N+1
      
          if(flag.start)
          {  flag.start<<-FALSE 
             temp.tot<-read.csv(file.name,na.string=".")
             N1<-N2<-nrow(temp.tot)
             if(svalue(BDS.list[[3]])!="NONE")
             {  N2<-sum(temp.tot[,svalue(BDS.list[[3]])])         
             }
#             tmp<-gframe("",cont=BBgroup,horizontal=FALSE)
#             glabel(paste("total number of observation = ", N1,"\n",sep=""),cont=tmp)
#             glabel(paste("total number of previous selection = ", N2,"\n",sep=""),cont=tmp)

             niter<-matrix(c("total N","0",N1,N2),ncol=2)
             colnames(niter)<-c("iteration","N")
             NList<<-gtable(niter,multiple=T)
             size(NList)<-c(150,200)             
             add(BigGroup,NList)
             add(BigGroup,ggraphics())                  
             var.name<-colnames(temp.tot)
             X.id<-which(var.name==svalue(BDS.list[[1]]))
             Y.id<-which(var.name==svalue(BDS.list[[2]]))
             presel.id<-which(var.name==svalue(BDS.list[[3]]))
             Y<-temp.tot[,Y.id]
             X<-temp.tot[,X.id]
             if(length(presel.id)!=0)
             {  presel<-temp.tot[,presel.id]
                sel.NULL<-ifelse(presel==1,TRUE,FALSE)
             } else
             {  presel<-NULL            
                sel.NULL<-rep(TRUE,length(Y))
             }
             sel.FINAL<-sel.NULL
             Y<<-Y
             X<<-X    
             flag.start<<-FALSE        
          }
          sel.NULL<-sel.FINAL
          flag.start<<-FALSE  
          Y.new<-Y[sel.NULL]
          X.new<-X[sel.NULL]     
          N<-length(Y.new)
          diff.Y<-(Y.new[-1]-Y.new[-length(Y.new)])/(X.new[-1]-X.new[-length(X.new)])
          cut.rate<-as.numeric(svalue(sel.rate))/100
          cut.n<-round(N*cut.rate)

          sel.id<-sort(sort.list(abs(diff.Y),decreasing=T)[-c(1:cut.n)])
          select.flag.temp<-rep(FALSE,length(X.new))
          select.flag.temp[c(1,sel.id+1)]<-TRUE
          sel.FINAL[sel.NULL]<-select.flag.temp
          sel.FINAL<<-sel.FINAL          
          plot(X.new[sel.id],Y.new[sel.id],type='l',col=3,xlab="TIME",ylab="DV",main=paste("N=",length(sel.id),sep=""))
          NList[]<-rbind(NList[],c(iter.N,length(sel.id)))
      }      
    
      BDSwin<<-gwindow("Biosignal data : Selection : Individual processing")
      BigGroup<<-ggroup(cont=BDSwin,horizontal=TRUE)
      BBgroup<<-ggroup(cont=BigGroup,horizontal=FALSE)
      NITER<<-1
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open data folder",handler=openBDS)
      add(tmp,button1)
      add(tmp,control.t)
   }    
   ######## Batch smoothing ##################################################  
   BDS.smooth.batch<-function(h,...)
   {  openBDS<-function(h,...)
      {   BDS.dir<<-gfile(text="Open biosignal data file folder",type="selectdir")
          svalue(control.t)<-BDS.dir
      }
      
      BDScalc<-function(h,...)
      {  N.BDS<-length(indiv.file)
         saveDIR<<-gfile(text="Select folder for data saving",type="selectdir")
         win<-gwindow(paste("Biosignal data : Central moving average : Barch progressing : N=",N.BDS,sep=""),width=300,height=50)
         BDS.progress<-gslider(from=0,to=N.BDS,by=1,value=0,cont=win)               
         for(i in 1:length(indiv.file))
         {   svalue(BDS.progress)<-i
             file.name<<-paste(strsplit(indiv.file[i],split="\\.")[[1]][1],"-",
                              svalue(BDS.list[[2]]),".csv",sep="")                
             full.file<-paste(saveDIR, file.name,sep="\\")                          
             temp.tot<-read.csv(indiv.file[i],na.string=".")
             var.name<-colnames(temp.tot)
             X.id<-which(var.name==svalue(BDS.list[[1]]))
             Y.id<-which(var.name==svalue(BDS.list[[2]]))             
             presel.id<-which(var.name==svalue(BDS.list[[3]]))
             if(length(presel.id)!=0)
             {  Y<-temp.tot[temp.tot[,presel.id],Y.id]
                X<-temp.tot[temp.tot[,presel.id],X.id]
             } else
             {  Y<-temp.tot[,Y.id]
                X<-temp.tot[,X.id]
             }               

             N<-length(Y)
             MovingW<-as.numeric(svalue(MW))
             addMW<-floor(MovingW/2)
             Y.new<-NULL
             for(i in 1:N)
             {  if(i < addMW)
                {  Y.new<-c(Y.new,mean(Y[1:(i+addMW)]))
                } else if(i >N-addMW)
                {  Y.new<-c(Y.new,mean(Y[(i-addMW):N]))
                } else
                {  Y.new<-c(Y.new,mean(Y[(i-addMW):(i+addMW)]))
                }  
             }
             Y<-rep(NA,nrow(temp.tot))
             if(length(presel.id)!=0)
             {  Y[temp.tot[,presel.id]]<-Y.new
             } else
             {  Y<-Y.new
             } 
             tot.data<-cbind(temp.tot,Y)

             colnames(tot.data)<-c(colnames(temp.tot),paste(svalue(BDS.list[[2]]),".smooth",sep=""))
             tot.data<<-tot.data            
             write.csv(tot.data,full.file)           
          }  
          dispose(win)
#          dispose(BDSwin)  
      }      
      
      TDselect<-function(h,...)
      {  setwd(BDS.dir)
         indiv.file<<-dir(BDS.dir)
         temp<-read.csv(indiv.file[1],na.string=".")
         var.name.temp<-colnames(temp)
         BDSparam.input<-c("TIME  ","DV ","Previous selection")
            
         BDS.g<-list()
         BDS.list<-list()
         for(i in 1:3)
         {  BDS.g[[i]]<-ggroup(cont=BBgroup)
            temp<- gdroplist(var.name.temp,cont=BDS.g[[i]])
            glabel(BDSparam.input[i],cont=BDS.g[[i]])
             BDS.list[[i]]<-temp
         }
         BDS.list<<-BDS.list           
         
         Button2<-gbutton("Calculate central moving average",handler= BDScalc)
         tmp=gframe("",container=BBgroup)
         add(tmp,Button2)
        
      }  
    
      BDSwin<<-gwindow("Biosignal data : Selection : Batch processing")
      BigGroup<<-ggroup(cont=BDSwin,horizontal=TRUE)
      BBgroup<<-ggroup(cont=BigGroup,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open data folder",handler=openBDS)
      add(tmp,button1)
      add(tmp,control.t)
      MW<<-gedit(" ",width=20)   
      tmp=gframe("Moving window (# of obs, odd number)",container=BBgroup)
      add(tmp,MW)                  
      Button<-gbutton("Select TIME and DV",handler= TDselect)
      tmp=gframe("",container=BBgroup)
      add(tmp,Button)
   }  
      
   ######## Individual smoothing ##################################################  
   BDS.smooth<-function(h,...)
   {
      BDSsave<-function(h,...)
      {    write.csv(tot.data,paste(gfile(text="Save central moving average data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)         
      }
              
      BDS.indiv.GUI<-function(h,...)
      {  file.name<<-svalue(FileList)
         temp.tot<-read.csv(file.name,na.string=".")
         var.name.temp<-c("NONE",colnames(temp.tot))
         BDSparam.input<-c("TIME  ","DV ","Previous selection")
      
         BDSINDwin<<-gwindow(paste("Biosignal data : Central moving average : Individual process",
         file.name,sep="---"))
         BigGroup<<-ggroup(cont=BDSINDwin,horizontal=TRUE)
         BBgroup<<-ggroup(cont=BigGroup,horizontal=FALSE)
         MW<<-gedit(" ",width=20)   
         tmp=gframe("Moving window (# of obs, odd number)",container=BBgroup)
         add(tmp,MW)  
         BDS.g<-list()
         BDS.list<-list()
         for(i in 1:3)
         {  BDS.g[[i]]<-ggroup(cont=BBgroup)
            temp<- gdroplist(var.name.temp,cont=BDS.g[[i]])
            glabel(BDSparam.input[i],cont=BDS.g[[i]])
             BDS.list[[i]]<-temp
         }
         BDS.list<<-BDS.list           
         Button2<-gbutton("Calculate central moving average",handler= BDSsmoothcalc)
         tmp<-gframe("",container=BBgroup)
         add(tmp,Button2)
         Button21<-gbutton("Draw plot",handler= BDSplot)
         tmp<-gframe("",container=BBgroup)
         add(tmp,Button21)         
         Button3<-gbutton("Save data",handler= BDSsave)
         tmp<-gframe("",container=BBgroup)
         add(tmp,Button3)
         flag.start<<-TRUE
      }
      
      BDSplot<-function(h,...)
      {   plot(X.plot,Y.plot,type='l',lty=2,col="grey",xlab="TIME",ylab="DV")
          lines(X.plot,Y.new.plot,col=2,lwd=2)
       }
      
      openBDS<-function(h,...)
      {   BDS.dir<<-gfile(text="Open biosignal data file folder",type="selectdir")
          svalue(control.t)<-BDS.dir
          setwd(BDS.dir)
          indiv.file<<-dir(BDS.dir)
          tmp<-gframe("Files",container=BBgroup)
          FileList<<-gtable(indiv.file,handler=BDS.indiv.GUI)
          size(FileList)<-c(500,200)
          add(tmp,FileList)
      }
 
      BDSsmoothcalc<-function(h,...)
      {   if(start.flag)
          {  add(BigGroup,ggraphics())  
          } 
          start.flag<<-FALSE
          temp.tot<-read.csv(file.name,na.string=".")
          var.name<-colnames(temp.tot)
          X.id<-which(var.name==svalue(BDS.list[[1]]))
          Y.id<-which(var.name==svalue(BDS.list[[2]]))
          presel.id<-which(var.name==svalue(BDS.list[[3]]))
          if(length(presel.id)!=0)
          {  Y<-temp.tot[temp.tot[,presel.id],Y.id]
             X<-temp.tot[temp.tot[,presel.id],X.id]
          } else
          {  Y<-temp.tot[,Y.id]
             X<-temp.tot[,X.id]
          }             
          N<-length(Y)
          MovingW<-as.numeric(svalue(MW))
          addMW<-floor(MovingW/2)
          Y.new<-NULL
          for(i in 1:N)
          {  if(i < addMW)
             {  Y.new<-c(Y.new,mean(Y[1:(i+addMW)]))
             } else if(i >N-addMW)
             {  Y.new<-c(Y.new,mean(Y[(i-addMW):N]))
             } else
             {  Y.new<-c(Y.new,mean(Y[(i-addMW):(i+addMW)]))
             }  
          }
          Y.t<-rep(NA,nrow(temp.tot))
          if(length(presel.id)!=0)
          {  Y.t[temp.tot[,presel.id]]<-Y.new
          } else
          {  Y.t<-Y.new
          } 
          tot.data<-cbind(temp.tot,Y.t)
          colnames(tot.data)<-c(colnames(temp.tot),paste(svalue(BDS.list[[2]]),".smooth",sep=""))
          tot.data<<-tot.data     
          X.plot<<-X
          Y.plot<<-Y
          Y.new.plot<<-Y.new        
          plot(X,Y,type='l',lty=2,col="grey",xlab="TIME",ylab="DV")
          lines(X,Y.new,col=2,lwd=2)
      }      
      

    
      BDSwin<<-gwindow("Biosignal data selection : Central moving average")
      BigGroup<<-ggroup(cont=BDSwin,horizontal=TRUE)
      BBgroup<<-ggroup(cont=BigGroup,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open data folder",handler=openBDS)
      add(tmp,button1)
      add(tmp,control.t)
      start.flag<<-TRUE
   } 
     
   #### NONMEM data ############################################################  
   ###### Notes before NONMEM data creation ####################################  
   dataNotes<-function(h,...)
   {  gmessage("*** Notes before NONMEM data creation ***
                \nIDs should be successive in ascending order.
                \n   1, 2, 3, 4, 5, ... (O)
                \n   1, 2, 5, 7, 8, ... (X)",cont=TRUE,width=600)   
   }
   
   ###### Create NONMEM data ###################################################  
   DataPrep<-function(h,...)
   {  DemogOK <- function()
      {  fileName1<<-gfile(text="Choose demographic data file",type="open")
         tclvalue(DemogName)<-fileName1
         dir.name<-strsplit(fileName1,split="\\.")[[1]][1]
         temp<-strsplit(dir.name,"\\\\")[[1]]    
         dataname<-temp[length(temp)]
         dir.name<-strsplit(dir.name,dataname)[[1]]  
         setwd(dir.name)       
      }
      
      AdmOK <- function()
      {  fileName2<<-gfile(text="Choose dosing data file",type="open")
         tclvalue(AdmName)<-fileName2
      }
      
      DVOK <- function()
      {  fileName3<<-gfile(text="Choose DV data file",type="open")
         tclvalue(DVName)<-fileName3
      }
      
      IPKOK<-function()
      {  fileName5<<-gfile(text="Choose IPK data file",type="open")
         tclvalue(IPKName)<-fileName5
      }

      displayInTable <- function(tclarray,title="",height=-1,width=-1,nrow=-1,ncol=-1)
      {  tt <- tktoplevel()
  	 tclRequire("Tktable")
  	 tkwm.title(tt,title)
  	 table1 <- tkwidget(tt,"table",rows=nrow,cols=ncol,titlerows=1,titlecols=0,
          	height=height+1,width=width+1,xscrollcommand=function(...) tkset(xscr,...),
	        yscrollcommand=function(...) tkset(yscr,...))
  	 xscr <-tkscrollbar(tt,orient="horizontal", command=function(...)tkxview(table1,...))
  	 yscr <- tkscrollbar(tt,command=function(...)tkyview(table1,...))
  	 tkgrid(table1,yscr)
  	 tkgrid.configure(yscr,sticky="nsw")
  	 tkgrid(xscr,sticky="new")
  	 tkconfigure(table1,variable=tclarray,background="white",selectmode="extended")
  	 return (table1)
      }

      openSpread<-function()
      {  NM.data.temp<-matrix(as.character(NM.data),ncol=ncol(NM.data))
   	 NM.data.temp<-rbind(colnames(NM.data),NM.data.temp)
    	 tclArray1 <- tclArray()
    	 for(i in (1:nrow(NM.data.temp)))
    	    for(j in (1:ncol(NM.data.temp)))
       	       tclArray1[[i-1,j-1]] <- NM.data.temp[i,j]
	 table1 <- displayInTable(tclArray1,nrow=nrow(NM.data.temp),ncol=ncol(NM.data.temp))
      }

      Combine<-function()
      {  Demog<<-read.csv(fileName1)
    	 Adm<<-read.csv(fileName2)
    	 DV<<-read.csv(fileName3)
    	 fileName5<<-tclvalue(IPKName)
    	 if(fileName5!="") IPK<<-read.csv(fileName5)
    	 ID.list<-unique(Demog$X.ID)
    	 if(fileName5!="")
    	 {  Demog.temp<-matrix(0,nrow=length(ID.list),ncol=(ncol(Demog)+ncol(IPK)-1))
       	    colnames(Demog.temp)<-c("X.ID",colnames(Demog)[-1],colnames(IPK)[-1])
       	    Demog.temp<-data.frame(Demog.temp)
       	    Demog.temp[,1]<-ID.list
       	    for(i in ID.list)
               Demog.temp[Demog.temp$X.ID==i,]<-c(Demog[Demog$X.ID==i,],IPK[IPK$X.ID==i,-1])
       	    Demog<<-Demog.temp
    	 }
    	 tot.data<-NULL
    	 colname.final1<-c(colnames(Adm),"DV","MDV")
    	 colname.final2<-c(colnames(Adm),"DV","MDV",colnames(Demog)[-1])
    	 for(i in ID.list) 
    	 {  DV.temp<-DV[which(DV$X.ID==i),]
    	    MDV<-rep(10,nrow(DV.temp))
    	    DV.temp<-cbind(DV.temp,MDV)
    	    temp.list<-which(Adm$X.ID==i)
    	    temp<-Adm[temp.list,]
    	    MDV<-rep(1,nrow(temp))
    	    temp<-cbind(temp,MDV)
     	    temp<-merge(temp,DV.temp,all=TRUE)   	      
     	    temp$MDV[temp$MDV==10]<-0
     	    temp<-temp[,colname.final1]
     	    Demog.temp<-matrix(rep(Demog[which(Demog$X.ID==i),-1],nrow(temp)),nrow=nrow(temp),byrow=T)
    	    temp<-cbind(temp,Demog.temp)
            colnames(temp)<-colname.final2
            tot.data<-rbind(tot.data,temp)
    	 }
    	 tot.data<-as.matrix(tot.data)
    	 temp<-as.character(tot.data)
    	 temp[which(temp=="NA")]<-"."
    	 tot.data<-matrix(temp,ncol=ncol(tot.data))
    	 colname.final2[1]<-"#ID"
    	 tot.data<-rbind(colname.final2,tot.data)
    	 NM.data<<-tot.data
    	 openSpread()
      }
      
      Save<-function()
      {  fileName4<-tclvalue(tkgetSaveFile(filetypes="{{CSV Files} {.csv}}")) 
    	 fileName4<-paste(fileName4,".csv",sep="")
    	 write.table(NM.data,fileName4,sep=",",quote=FALSE,row.names=FALSE, col.names=FALSE)
    	 tkdestroy(Toptt)
      }

      Demog<<-NULL
      Adm<<-NULL
      DV<<-NULL
      IPK<<-NULL
      Toptt<<-tktoplevel()
      tkwm.title(Toptt,"NM data preparation for PREDPP")
      tt<-tkframe(Toptt)
      DemogName <- tclVar("")
      Demog.Name <-tkentry(tt,width="60",textvariable=DemogName)
      IPKName <- tclVar("")
      IPK.Name <-tkentry(tt,width="60",textvariable=IPKName)
      AdmName <- tclVar("")
      Adm.Name <-tkentry(tt,width="60",textvariable=AdmName)
      DVName <- tclVar("")
      DV.Name <-tkentry(tt,width="60",textvariable=DVName)
      OK.but1 <-tkbutton(tt,text="   Select   ",width=7,height=1,command=DemogOK)
      OK.but2 <-tkbutton(tt,text="   Select   ",width=7,height=1,command=AdmOK)
      OK.but3 <-tkbutton(tt,text="   Select   ",width=7,height=1,command=DVOK)
      OK.but4 <-tkbutton(tt,text="   Select   ",width=7,height=1,command=IPKOK)
      CombineOK<-tkbutton(tt,text="Combine",width=7,height=1,command=Combine)
      SaveOK<-tkbutton(tt,text="    Save   ",width=7,height=1,command=Save)
      tkgrid(tklabel(tt,text="Demographics (#ID, covariates)"),Demog.Name,OK.but1)
      tkgrid(tklabel(tt,text="IPK (#ID, IPK)"),IPK.Name,OK.but4)
      tkgrid(tklabel(tt,text="Dosing (#ID, TIME (elapse), AMT,RATE)"),Adm.Name,OK.but2)
      tkgrid(tklabel(tt,text="DV (#ID, TIME (elapse), DV)"),DV.Name,OK.but3)
      tkgrid(tklabel(tt,text=" "),tklabel(tt,text=" "),CombineOK)
      tkgrid(tklabel(tt,text=" "),tklabel(tt,text=" "),SaveOK)
      tkgrid(tt)
   }   
   
   ###### Create successive ID #################################################  
   data.ID<-function(h,...)
   {  data.file<-gfile(text="Choose data file",type="open")
      DD.data<- read.csv(data.file,na.string=".")
      rep.n<-table(DD.data$X.ID)
      seq.id<-1:length(rep.n)
      new.id<-rep(seq.id,rep.n)
      new.data<-cbind(new.id,DD.data)
      colnames(new.data)[1:2]<-c("#ID","OID")
      dir.name<-strsplit(data.file,split="\\.")[[1]][1]
      temp<-strsplit(dir.name,"\\\\")[[1]]    
      dataname<-temp[length(temp)]
      dir.name<-strsplit(dir.name,dataname)[[1]]
      setwd(dir.name)   
      write.csv(new.data,paste(gfile(text="Save as csv",
                 type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
   }
   #### Explore NONMEM data ####################################################  
   ###### Select data file #####################################################  
   OpenEDAData<-function(h,...)
   {  EDAfileName<<-gfile(text="Choose EDA data file",type="open")
      EDA.data<<-read.csv(EDAfileName,na.string=".")
      Var.Name<<-colnames(EDA.data)
      dir.name<-strsplit(EDAfileName,split="\\.")[[1]][1]
      temp<-strsplit(dir.name,"\\\\")[[1]]    
      dataname<-temp[length(temp)]
      dir.name<-strsplit(dir.name,dataname)[[1]]
      setwd(dir.name)
   }
   
   ###### Summary ##############################################################  
   ######## Summary statistics-continuous ######################################  
   Summary.stat<-function(h,...)
   {  calc.summary<-function()
      {  DA.data<-as.matrix(EDA.data[,Con.list])
         colnames(DA.data)<-Con.list
         summary.stat1<-rbind(apply(DA.data,2,function(x) mean(x,na.rm=T)),
         apply(DA.data,2,function(x) sd(x,na.rm=T)), 
         apply(DA.data,2,function(x) min(x,na.rm=T)),
         apply(DA.data,2,function(x) quantile(x,na.rm=T,probs=0.25)),
         apply(DA.data,2,function(x) quantile(x,na.rm=T,probs=0.5)),
         apply(DA.data,2,function(x) quantile(x,na.rm=T,probs=0.75)),
         apply(DA.data,2,function(x) max(x,na.rm=T)))
         summary.stat1<-t(summary.stat1)
         colnames(summary.stat1)<-c("Mean","SD","Mininum","Q1","Median","Q3","Maximum")
         summary.stat<-round(summary.stat1,3)
         summary.stat<-cbind(c(rownames(summary.stat1)),summary.stat)
         savesummarydata<-function(h,...)
         {  dir.name<-strsplit(EDAfileName,split="\\.")[[1]][1]
            temp<-strsplit(dir.name,"\\\\")[[1]]    
            dataname<-temp[length(temp)]
            dir.name<-strsplit(dir.name,dataname)[[1]]
            setwd(dir.name)
            write.csv(summary.stat,paste(gfile(text="Save as csv",
                 type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         }
         summary.w<-gwindow("Summary statistics")
         gsummary<-ggroup(cont=summary.w,horizontal=F)
         summary.table<-gtable(summary.stat,do.subset=TRUE,cont=gsummary)
         size(summary.table)<-c(20,200)
         tmp<-gframe(cont=gsummary,spacing=10000)
         Button1<-gbutton("Save",handler=savesummarydata,spacing=2000)
         size(Button1)<-c(50,30)
         add(tmp,Button1)
      }

      saveCat<-function(h,...)
      {  Con.list<<-Var.Name[svalue(catcheck)]
         dispose(checkg)
         calc.summary()
      }

      Var.Name<-colnames(EDA.data)
      checkg<-gwindow("Select continuous variable for summary statistics")
      catcheck<-gcheckboxgroup(Var.Name,use.table=TRUE,cont=checkg)
      Button1<-gbutton("OK",type="OK",handler=saveCat,cont=checkg)
   }   
   ######## Summary statistics-continuous by ID ################################
   Summary.stat.ind<-function(h,...)
   {  calc.summary<-function()
      {  DA.data<-c(EDA.data[,Con.list])
         if(length(Cat.list)==0)
         {  Con.data<-as.matrix(EDA.data[,c("X.ID")]) 
         } else
         {  Con.data<-as.matrix(EDA.data[,c("X.ID",Cat.list)])
         }
         fp<-function(x)
         { n<-length(x)
           name<-x[1]
           if(n>=2)
             for(i in 2:n)
               name<-paste(name,x[i],sep="-")
           return(name)
         }      
         ID<-apply(Con.data,1,fp)
         summary.stat1<-cbind(tapply(DA.data,ID,function(x) length(x[!is.na(x)])),
         tapply(DA.data,ID,function(x) mean(x,na.rm=T)),
         tapply(DA.data,ID,function(x) sd(x,na.rm=T)), 
         tapply(DA.data,ID,function(x) min(x,na.rm=T)),
         tapply(DA.data,ID,function(x) quantile(x,na.rm=T,probs=0.25)),
         tapply(DA.data,ID,function(x) quantile(x,na.rm=T,probs=0.5)),
         tapply(DA.data,ID,function(x) quantile(x,na.rm=T,probs=0.75)),
         tapply(DA.data,ID,function(x) max(x,na.rm=T)))
         colnames(summary.stat1)<-c("N","Mean","SD","Mininum","Q1","Median","Q3","Maximum")
         summary.stat<-round(summary.stat1,3)
         
         cat.name<-matrix(unlist(strsplit(rownames(summary.stat1),split="-")),byrow=T,nrow=nrow(summary.stat))
         if(length(Cat.list)==0)
         {  colnames(cat.name)<-"ID"
         } else
         { colnames(cat.name)<-c("ID",Cat.list)
         }
         summary.stat<-cbind(cat.name,summary.stat)
         savesummarydata<-function(h,...)
         {  dir.name<-strsplit(EDAfileName,split="\\.")[[1]][1]
            temp<-strsplit(dir.name,"\\\\")[[1]]    
            dataname<-temp[length(temp)]
            dir.name<-strsplit(dir.name,dataname)[[1]]
            setwd(dir.name)
            write.csv(summary.stat,paste(gfile(text="Save as csv",
                 type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         }
         summary.w<-gwindow(paste("Summary statistics by ID : ",Con.list ))
         gsummary<-ggroup(cont=summary.w,horizontal=F)
         summary.table<-gtable(summary.stat,do.subset=TRUE,cont=gsummary)
         size(summary.table)<-c(20,200)
         tmp<-gframe(cont=gsummary,spacing=10000)
         Button1<-gbutton("Save",handler=savesummarydata,spacing=2000)
         size(Button1)<-c(50,30)
         add(tmp,Button1)
      }

      saveCat<-function(h,...)
      {  Con.list<<-svalue(catcheck)
         Cat.list<<-Var.Name[svalue(catcheck1)]
         dispose(checkg)
         calc.summary()
      }

      Var.Name<-colnames(EDA.data)[-1]
      checkg<-gwindow("Summary statistics by ID")
      tmp<-gframe("Select continuous variable", cont=checkg)      
      catcheck<-gradio(Var.Name,use.table=TRUE,cont=tmp) 
      gg<-ggroup(cont=checkg)
      tmp<-gframe("Select categorical variables for levels", cont=gg,horizontal=FALSE)      
      catcheck1<-gcheckboxgroup(Var.Name,use.table=TRUE)
      size(catcheck1)<-c(200,300)
      add(tmp,catcheck1)
      Button1<-gbutton("OK",type="OK",handler=saveCat,cont=checkg)
   }
   
   ######## Summary statistics-categorical by ID ################################
   Summary.stat.cat.ind<-function(h,...)
   {  calc.summary<-function()
      {  DA.data<-data.frame(EDA.data[,c(Con.list,"X.ID")])
         sum.cat<-as.data.frame(table(DA.data))
         p<-ncol(sum.cat)
         summary.stat<-sum.cat[,c(p-1,1:(p-2),p)]
         savesummarydata<-function(h,...)
         {  dir.name<-strsplit(EDAfileName,split="\\.")[[1]][1]
            temp<-strsplit(dir.name,"\\\\")[[1]]    
            dataname<-temp[length(temp)]
            dir.name<-strsplit(dir.name,dataname)[[1]]
            setwd(dir.name)
            write.csv(summary.stat,paste(gfile(text="Save as csv",
                 type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         }
         title.N<-"Summary table by ID :"
         for(i in 1:length(Con.list))
           title.N<-paste(title.N,Con.list[i])         
         summary.w<-gwindow(title.N)
         gsummary<-ggroup(cont=summary.w,horizontal=F)
         summary.table<-gtable(summary.stat,do.subset=TRUE,cont=gsummary)
         size(summary.table)<-c(20,200)
         tmp<-gframe(cont=gsummary,spacing=10000)
         Button1<-gbutton("Save",handler=savesummarydata,spacing=2000)
         size(Button1)<-c(50,30)
         add(tmp,Button1)
      }

      saveCat<-function(h,...)
      {  Con.list<<-Var.Name[svalue(catcheck)]
         dispose(checkg)
         calc.summary()
      }

      Var.Name<-colnames(EDA.data)[-1]
      checkg<-gwindow("Select categorical variables for summary statistics by ID")
      catcheck<-gcheckboxgroup(Var.Name,use.table=TRUE,cont=checkg)
      Button1<-gbutton("OK",type="OK",handler=saveCat,cont=checkg)
   }   
        
   ######## Summary statistics-categorical-single level per person #############  
   Summary.cat<-function(h,...)
   {  saveCat<-function(h,...)
      {  savesummarycdata<-function(h,...)
         {  dir.name<-strsplit(EDAfileName,split="\\.")[[1]][1]
            temp<-strsplit(dir.name,"\\\\")[[1]]    
            dataname<-temp[length(temp)]
            dir.name<-strsplit(dir.name,dataname)[[1]]
            setwd(dir.name)
            write.csv(Cat.summary,paste(gfile(text="Save as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         }
         
         Cat.list<<-Var.Name[svalue(catcheck)]
         dispose(checkg)
         n.table<-tapply(rep(1,nrow(EDA.data)),EDA.data$X.ID,sum)
         temp<-EDA.data[1,]   
         for(i in 2:length(n.table)) 
            temp<-rbind(temp,c(EDA.data[sum(n.table[1:(i-1)])+1,]))
         cont.list<-list()
         for(i in 1:length(Cat.list))
         {  cont.list[[i]]<-table(temp[,Cat.list[i]])
         }
         names(cont.list)<-Cat.list
         p<-max(unlist(lapply(cont.list,length)))
         n<-length(Cat.list)*2
         Cat.summary<-matrix(" ",nrow=n,ncol=(p+1))
         colnames(Cat.summary)<-c("Variable",paste("level",1:p,sep=""))
         for(i in 1:length(Cat.list))
         {  Cat.summary[(i-1)*2+1,1]<-Cat.list[i]
            Cat.summary[(i-1)*2+1,2:(length(cont.list[[i]])+1)]<-names(cont.list[[i]])
            Cat.summary[(i-1)*2+2,2:(length(cont.list[[i]])+1)]<-cont.list[[i]]
         }    
         summary.cw<-gwindow("Summary Categorical variable Statistics")
         gsummary<-ggroup(cont=summary.cw,horizontal=F)
         summary.table<-gtable(Cat.summary,do.subset=TRUE,cont=gsummary)
         size(summary.table)<-c(20,200)
         tmp<-gframe(cont=gsummary,spacing=10000)
         Button1<-gbutton("Save",handler=savesummarycdata,spacing=2000)
         size(Button1)<-c(50,30)
         add(tmp,Button1)
      }
      
      Var.Name<-colnames(EDA.data)
      checkg<-gwindow("Select categorical variable for summary statistics")
      catcheck<-gcheckboxgroup(Var.Name,use.table=TRUE,cont=checkg)
      Button1<-gbutton("OK",type="OK",handler=saveCat,cont=checkg)
   }
      
   ######## Summary statistics-categorical-multiple levels per person ##########  
   Summary.cat1<-function(h,...)
   {  saveCat<-function(h,...)
      {  savesummarycdata<-function(h,...)
         {  dir.name<-strsplit(EDAfileName,split="\\.")[[1]][1]
            temp<-strsplit(dir.name,"\\\\")[[1]]    
            dataname<-temp[length(temp)]
            dir.name<-strsplit(dir.name,dataname)[[1]]
            setwd(dir.name)
            write.csv(Cat.summary,paste(gfile(text="Save as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         }
         
         Cat.list<<-Var.Name[svalue(catcheck)]
         dispose(checkg)
         n.table<-tapply(rep(1,nrow(EDA.data)),EDA.data$X.ID,sum)
         temp<-EDA.data
         cont.list<-list()
         for(i in 1:length(Cat.list))
         {  IDs<-EDA.data$X.ID
            X<-temp[,Cat.list[i]]
            temp.A<-paste(IDs,X,sep="-")
            cont.list[[i]]<- table(X[!duplicated(temp.A)])
         }
         names(cont.list)<-Cat.list
         p<-max(unlist(lapply(cont.list,length)))
         n<-length(Cat.list)*2
         Cat.summary<-matrix(" ",nrow=n,ncol=(p+1))
         colnames(Cat.summary)<-c("Variable",paste("level",1:p,sep=""))
         for(i in 1:length(Cat.list))
         {  Cat.summary[(i-1)*2+1,1]<-Cat.list[i]
            Cat.summary[(i-1)*2+1,2:(length(cont.list[[i]])+1)]<-names(cont.list[[i]])
            Cat.summary[(i-1)*2+2,2:(length(cont.list[[i]])+1)]<-cont.list[[i]]
         }
         summary.cw<-gwindow("Summary Categorical variable Statistics")
         gsummary<-ggroup(cont=summary.cw,horizontal=F)
         summary.table<-gtable(Cat.summary,do.subset=TRUE,cont=gsummary)
         size(summary.table)<-c(20,200)
         tmp<-gframe(cont=gsummary,spacing=10000)
         Button1<-gbutton("Save",handler=savesummarycdata,spacing=2000)
         size(Button1)<-c(50,30)
         add(tmp,Button1)
      }
      
      Var.Name<-colnames(EDA.data)
      checkg<-gwindow("Select categorical variable for summary statistics")
      catcheck<-gcheckboxgroup(Var.Name,use.table=TRUE,cont=checkg)
      Button1<-gbutton("OK",type="OK",handler=saveCat,cont=checkg)
   }
   
   ###### Plot #################################################################  
   ######## XY plot ############################################################  

   XY.plot<-function(h,...)
   {  D.data<-EDA.data
      ID.id<-NULL
      if(sum(colnames(D.data)=="X.ID")!=0)
         ID.id<-which(colnames(D.data)=="X.ID")
      updatePlot<-function(h,...)
      {  condX.V <-svalue(VarList.X,index=T)
         condY.V<-svalue(VarList.Y,index=T)
         select.data<-D.data
         if(!is.na(condX.V)& !is.na(condY.V))
         {  X<-select.data[,condX.V]
            Y<-select.data[,condY.V]
            dev.set(which=Dev.XYplot)
            plot(X,Y,xlab=Var.Name[condX.V],ylab=Var.Name[condY.V])
         }
      }
      saveData<-function(h,...)
      {  condX.V <-svalue(VarList.X,index=T)
         condY.V<-svalue(VarList.Y,index=T)
         select.data<-D.data[,c(ID.id,condX.V,condY.V)]
         if(is.null(ID.id))
         {  colnames(select.data)<-Var.Name[c(condX.V,condY.V)]
         } else
         {  colnames(select.data)<-Var.Name[c(ID.id,condX.V,condY.V)]
         }
         write.csv(select.data,paste(gfile(text="Save as csv",
             type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)          
      }
      VarList.X<-gdroplist(Var.Name)
      VarList.Y<-gdroplist(Var.Name)
      Button1<-gbutton("OK",handler=updatePlot)
      Button2<-gbutton("Save",handler=saveData) 
      win<-gwindow("XY plot")
      BigGroup<-ggroup(cont=win)
      group<-ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe(" X variable",container=group)
      add(tmp,VarList.X)
      tmp<-gframe(" Y variable",container=group)
      add(tmp,VarList.Y)
      tmp<-gframe("Plot",container=group)
      add(tmp,Button1,expand=TRUE)
      tmp<-gframe("Save",container=group)
      add(tmp,Button2,expand=TRUE)
      add(BigGroup,ggraphics())
      Dev.XYplot<-dev.cur()
   }   
   
   ######## DV vs TIME by ID ################################################### 
   ID.plot<-function(h,...)
   {  D.data<-EDA.data
      ID.id<-NULL
      if(sum(colnames(D.data)=="X.ID")!=0)
         ID.id<-which(colnames(D.data)=="X.ID")
      TIME.id<-which(tolower(colnames(D.data))=="time")
      DV.id<-which(tolower(colnames(D.data))=="dv")      
      if(is.null(ID.id))  
      {  gconfirm("Need ID information",icon="warning")
      } else
      {  ID<-sort(unique(D.data[,ID.id]))
         ID<-matrix(ID,nrow=length(ID))
         colnames(ID)<-c("ID")
         updatePlot<-function(h,...) 
         {  dev.set(which=Dev.IDplot)
            plot(D.data[,TIME.id],D.data[,DV.id],type='n',xlab="TIME",ylab="DV")
            id<-names(table(D.data[,ID.id]))
            for(i in 1:length(id))
            {  data<-D.data[which(D.data[,ID.id]==as.numeric(id[i])),]
               data<-data[!is.na(data[,DV.id]),]
               lines(data[,TIME.id],data[,DV.id],lty=3)
            }
            select.id<-svalue(IDlist)
            for(kk in 1:length(select.id))
            {  data<-D.data[which(D.data[,ID.id]==select.id[kk]),]
               data<-data[!is.na(data[,DV.id]),]
               lines(data[,TIME.id],data[,DV.id],col=2,lwd=2)
            }   
         }
   
         savedata<-function(h,...)
         {  header<-strsplit(EDAfileName,"\\.")[[1]]
            select.id<-svalue(IDlist)
            sel.data<-NULL
            for(i in select.id)
            {  sel.id<-which(EDA.data$X.ID==i)
               sel.data<-rbind(sel.data,EDA.data[sel.id,]) 
            }   
            write.csv(sel.data,paste(gfile(text="Save as csv",
               type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         }
         IDlist <- gtable(ID,multiple=T,handler=updatePlot) 
         window <- gwindow("DV vs TIME by ID")
         Biggroup<-ggroup(cont=window,horizontal=TRUE)
         group<-ggroup(cont=Biggroup,horizontal=FALSE)
         tmp<-gframe("ID",container=group)
         size(IDlist)<-c(100,200)
         add(tmp,IDlist)
         button1<-gbutton("OK",handler=updatePlot)
         tmp<-gframe("Display selected IDs",container=group)
         add(tmp,button1)
         button2<-gbutton("Save",handler=savedata)     
         tmp<-gframe("Save selected IDs",container=group)
         add(tmp,button2)
         add(Biggroup,tmp)
         add(Biggroup, ggraphics())
         Dev.IDplot<-dev.cur()
         par(mfrow=c(1,1))
         plot(D.data[,TIME.id],D.data[,DV.id],type='n',xlab="TIME",ylab="DV")
         id<-names(table(D.data[,ID.id]))
         for(i in 1:length(id))
         {  data<-D.data[which(D.data[,ID.id]==as.numeric(id[i])),]
            data<-data[!is.na(data[,DV.id]),]
            lines(data[,TIME.id],data[,DV.id])
         }
      }
   }
       
   ######## DV vs TIME by covariates ###########################################  
   IDCOV.plot<-function(h,...)
   {  D.data<-EDA.data
      ID.id<-NULL
      if(sum(colnames(D.data)=="X.ID")!=0)
         ID.id<-which(colnames(D.data)=="X.ID")
      TIME.id<-which(tolower(colnames(D.data))=="time")
      DV.id<-which(tolower(colnames(D.data))=="dv")
      Var.Name<-colnames(D.data)
      cov.list<-1:ncol(D.data)
      COV.data<-as.matrix(D.data[,c(ID.id,cov.list)])
      colnames(COV.data)<-c("ID",Var.Name[cov.list])
      sample.data<-D.data
      sample.data[is.na(sample.data)]<-"NA"
      colnames(sample.data)<-Var.Name

      updatePlot<-function(h,...)
      {  cond1.V<-svalue(VarList.g1)
         cond1.C<-svalue(VarType.g1)
         cond1.F<-as.numeric(svalue(From.id.g1))
         cond1.T<-as.numeric(svalue(To.id.g1))
         cond2.V<-svalue(VarList.g2)
         cond2.C<-svalue(VarType.g2)
         cond2.F<-as.numeric(svalue(From.id.g2))
         cond2.T<-as.numeric(svalue(To.id.g2))
         select.data<-D.data
         if(cond1.V!="NONE     ")
         {  if(!is.na(cond1.F))
            {  selected.id<-which(select.data[,which(Var.Name==cond1.V)]>=cond1.F)
               select.data<-select.data[selected.id,]
            }
            if(!is.na(cond1.T))
            {  selected.id<-which(select.data[,which(Var.Name==cond1.V)]<=cond1.T)
               select.data<-select.data[selected.id,]
            }
         }
         if(cond2.V!="NONE     ")
         {  if(!is.na(cond2.F))
            {  selected.id<-which(select.data[,which(Var.Name==cond2.V)]>=cond2.F)
               select.data<-select.data[selected.id,]
            }
            if(!is.na(cond2.T))
            {  selected.id<-which(select.data[,which(Var.Name==cond2.V)]<=cond2.T)
               select.data<-select.data[selected.id,]
            }
         }
         if(nrow(select.data)==0)
         {  gmessage("No data!",title="Error",icon="error") 
         } else if(svalue(VarType.g1)!="Categorical-multiple levels per person")
         {  id<-names(table(select.data[,1]))
            data<-D.data[which(D.data[,ID.id]==id[1]),c(ID.id,DV.id,TIME.id)]
            data<-data[!is.na(data[,2]),]
            dev.set(which=Dev.IDCOVplot)
            plot(data[,3],data[,2],type='l',xlab="TIME",ylab="DV",
               xlim=range(D.data[,TIME.id],na.rm=T),ylim=range(D.data[,DV.id],na.rm=T))
            if(length(id)!=1)
               for(i in 2:length(id))
               {  data<-D.data[which(D.data[,ID.id]==id[i]),c(ID.id,DV.id,TIME.id)]
                  data<-data[!is.na(data[,2]),]
                  lines(data[,3],data[,2])
               }
         } else
         {  id<-names(table(select.data[,1]))
            data.D<-EDA.data
            data.D<-data.D[!is.na(data.D[,"DV"]),]
            dev.set(which=Dev.IDCOVplot)
            plot(data.D[,"TIME"],data.D[,"DV"],type='n',xlab="TIME",ylab="DV",
                  xlim=range(data.D[,"TIME"],na.rm=T),ylim=range(data.D[,"DV"],na.rm=T))
            for(i in 1:length(id))
            {  data<-select.data[which(select.data[,"X.ID"]==id[i]),c(ID.id,DV.id,TIME.id)]
               data<-data[!is.na(data[,"DV"]),]
               lines(data[,"TIME"],data[,"DV"])
            }
         }    
      }
    
      saveData<-function(h,...)
      {  cond1.V<-svalue(VarList.g1)
         cond1.C<-svalue(VarType.g1)
         cond1.F<-as.numeric(svalue(From.id.g1))
         cond1.T<-as.numeric(svalue(To.id.g1))
         cond2.V<-svalue(VarList.g2)
         cond2.C<-svalue(VarType.g2)
         cond2.F<-as.numeric(svalue(From.id.g2))
         cond2.T<-as.numeric(svalue(To.id.g2))
         select.data<-D.data
         if(cond1.V!="NONE")
         {  if(!is.na(cond1.F))
            {  selected.id<-which(select.data[,which(Var.Name==cond1.V)]>=cond1.F)
               select.data<-select.data[selected.id,]
            }
            if(!is.na(cond1.T))
            {  selected.id<-which(select.data[,which(Var.Name==cond1.V)]<=cond1.T)
               select.data<-select.data[selected.id,]
            }
         }
         if(cond2.V!="NONE")
         {  if(!is.na(cond2.F))
            {  selected.id<-which(select.data[,which(Var.Name==cond2.V)]>=cond2.F)
               select.data<-select.data[selected.id,]
            }
            if(!is.na(cond2.T))
            {  selected.id<-which(select.data[,which(Var.Name==cond2.V)]<=cond2.T)
               select.data<-select.data[selected.id,]
            }
         }
         if(nrow(select.data)==0)
         {  gmessage("No data!",title="Error",icon="error") 
         } else  if(svalue(VarType.g1)!="Categorical-multiple levels per person")
         {  id<-names(table(select.data[,1]))
            data<-NULL
            for(i in 1:length(id))
               data<-rbind(data,D.data[which(D.data[,ID.id]==id[i]),])
            write.csv(data,paste(gfile(text="Save as csv",
               type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
         } else
         {  write.csv(select.data,paste(gfile(text="Save as csv",
               type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)

         }      
      }

      From.label.1<-glabel("From")
      To.label.1<-glabel("To")
      VarList.g1<-gdroplist(c("NONE     ",Var.Name))
      VarType.g1<-gradio(c("Categorical-single level per person","Categorical-multiple levels per person","Continuous"))
      From.id.g1<-gedit(" ",width=10)
      To.id.g1<-gedit(" ",width=10)
      From.label.2<-glabel("From")
      To.label.2<-glabel("To")
      VarList.g2<-gdroplist(c("NONE     ",Var.Name))
      VarType.g2<-gradio(c("Categorical-single level per person","Categorical-multiple levels per person","Continuous"))
      From.id.g2<-gedit(" ",width=10)
      To.id.g2<-gedit(" ",width=10)
      Button1<-gbutton("OK",handler=updatePlot)
      Button2<-gbutton("Save",handler=saveData)
      win<-gwindow("DV vs TIME by covariates")
      BigGroup<-ggroup(cont=win)
      group<-ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe("Covariate 1 : ",container=group)
      add(tmp,VarList.g1)
      add(tmp,VarType.g1)
      tmp<-gframe("Select  : ",container=group)
      add(tmp,From.label.1);add(tmp,From.id.g1)
      add(tmp,To.label.1);add(tmp,To.id.g1)
      tmp<-gframe("Covariate 2 : ",container=group)
      add(tmp,VarList.g2)
      add(tmp,VarType.g2)
      tmp<-gframe("Select :",container=group)
      add(tmp,From.label.2);add(tmp,From.id.g2)
      add(tmp,To.label.2);add(tmp,To.id.g2)
      tmp<-gframe("Plot",container=group)
      add(tmp,Button1,expand=TRUE)
      tmp<-gframe("Save",container=group)
      add(tmp,Button2,expand=TRUE)   
      add(BigGroup,ggraphics())
      Dev.IDCOVplot<-dev.cur()
   }  
    
   ######## Covariates vs covariates ###########################################  
   COVvsCOV.plot<-function(h,...)
   {  CovPlot<-function()
      {  updatePlot<-function(h,...)
         {  ID.id<-NULL
            Var.Name<-Cov.list
            if(sum(colnames(EDA.data)=="X.ID")!=0)
               ID.id<-which(colnames(EDA.data)=="X.ID")
            D.data<-EDA.data
            DD1<-svalue(VarList.cond1)
            DD1.F<-svalue(From.id.cond1)
            DD1.T<-svalue(To.id.cond1)
            DD2<-svalue(VarList.cond2)
            DD2.F<-svalue(From.id.cond2)
            DD2.T<-svalue(To.id.cond2)

            if(DD1!="NONE     ")
            {  if(DD2!="NONE     ")
               {  if(DD1.F==" " | DD2.F==" " | DD1.T ==" " |DD2.T==" ")
                  {  gmessage("Enter select condition",icon="error")
                  } else
                  {  sel.id<-which(D.data[,DD1]>=as.numeric(svalue(DD1.F)) & D.data[,DD1]<=as.numeric(svalue(DD1.T)) &
                              D.data[,DD2]>=as.numeric(svalue(DD2.F)) & D.data[,DD2]<=as.numeric(svalue(DD2.T)) )
                  }
               } else
               {  sel.id<-which(D.data[,DD1]>=as.numeric(svalue(DD1.F)) & D.data[,DD1]<=as.numeric(svalue(DD1.T)))
               }    
            } else
            {  sel.id<-1:nrow(D.data)
            }

            temp.data<-EDA.data[sel.id,]
            n.table<-tapply(rep(1,nrow(temp.data)),temp.data$X.ID,sum)
            temp<-apply(temp.data[1:n.table[1],],2,function(x) median(x,na.rm=T))
      
            for(i in 2:length(n.table)) 
               temp<-rbind(temp,apply(temp.data[(sum(n.table[1:(i-1)])+1):sum(n.table[1:i]),],2,function(x) median(x,na.rm=T)))

            D.data<-temp[,c("X.ID",Cov.list)]
            DD<<-D.data
      
            ID.id<-NULL
            Var.Name11<-colnames(D.data)
            if(sum(colnames(D.data)=="X.ID")!=0)
               ID.id<-which(colnames(D.data)=="X.ID")
            condX.V <-svalue(VarList.X)
            condY.V<-svalue(VarList.Y)
            select.data<-D.data
            if(!is.na(condX.V)& !is.na(condY.V))
            {  X<-select.data[,condX.V]; if(svalue(VarType.g1)=="Categorical") X<-factor(X)
               Y<-select.data[,condY.V]; if(svalue(VarType.g2)=="Categorical") Y<-factor(Y)
               dev.set(which=Dev.COVvsCOVplot)  
               plot(X,Y,xlab=condX.V,ylab=condY.V)
               if(svalue(VarType.g1)!="Categorical" & svalue(VarType.g2)!="Categorical")
                  lines(lowess(X,Y),col=2,lwd=2)
            } 
         }

         saveData<-function(h,...)
         {  condX.V <-svalue(VarList.X,index=T)+1
            condY.V<-svalue(VarList.Y,index=T)+1
            select.data<-D.data[,c(ID.id,condX.V,condY.V)]
            write.csv(select.data,paste(gfile(text="Save as csv",
                  type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)          
         }

         Var.Name<-Cov.list
         From.label.1<-glabel("From")
         To.label.1<-glabel("To")
         VarList.cond1<-gdroplist(c("NONE     ",Var.Name))
         From.id.cond1<-gedit(" ",width=10)
         To.id.cond1<-gedit(" ",width=10)

         From.label.2<-glabel("From")
         To.label.2<-glabel("To")
         VarList.cond2<-gdroplist(c("NONE     ",Var.Name))
         From.id.cond2<-gedit(" ",width=10)
         To.id.cond2<-gedit(" ",width=10)

         VarList.X<-gdroplist(Var.Name)
         VarType.g1<-gradio(c("Categorical","Continuous"))
         VarList.Y<-gdroplist(Var.Name)
         VarType.g2<-gradio(c("Categorical","Continuous"))
         Button1<-gbutton("OK",handler=updatePlot)
         Button2<-gbutton("SAVE",handler=saveData)
    
         win<-gwindow("Covariate vs covariate plot")
         BigGroup<-ggroup(cont=win)
         group<-ggroup(horizontal=FALSE,cont=BigGroup)
         tmp<-gframe(" Select ",container=group)
         add(tmp,VarList.cond1)
         add(tmp,From.label.1);add(tmp,From.id.cond1)
         add(tmp,To.label.1);add(tmp,To.id.cond1)
         tmp<-gframe(" Select ",container=group)
         add(tmp,VarList.cond2)
         add(tmp,From.label.2);add(tmp,From.id.cond2)
         add(tmp,To.label.2);add(tmp,To.id.cond2)

         tmp<-gframe(" Covariate 1 : X-axis",container=group)
         add(tmp,VarList.X)
         add(tmp,VarType.g1)
         tmp<-gframe(" Covariate 2 : Y-axis",container=group)
         add(tmp,VarList.Y)
         add(tmp,VarType.g2)

         tmp<-gframe("Plot",container=group)
         add(tmp,Button1,expand=TRUE)
         tmp<-gframe("Save",container=group)
         add(tmp,Button2,expand=TRUE)
         add(BigGroup,ggraphics())
         Dev.COVvsCOVplot<-dev.cur()
      }

      saveCov<-function(h,...)
      {  Cov.list<<-Var.Name1[svalue(groupcheck)]
         dispose(checkg)
         CovPlot()
      }

      Var.Name1<-colnames(EDA.data)
      checkg<-gwindow("Covariate selection",width=200,height=400)
      groupcheck<-gcheckboxgroup(Var.Name1,cont=checkg,use.table=TRUE)
      Button1<-gbutton("OK",type="OK",handler=saveCov,cont=checkg)
   }

   ############################################################################# 
   # Control stream
   #############################################################################

   #### Edit with default editor ###############################################
   EditEditor<-function(h,...)
   {  file.ctl<-gfile(text="Choose NONMEM control file",type="open")
      D.temp<-readLines(file.ctl)
      Current.CTL<<-D.temp
      NONMEM.CTL<<-D.temp
      if(is.null(Current.CTL))
      {  edit.txt<-" "
      } else
      {  edit.txt<-Current.CTL[1]
         for(i in 2:length(Current.CTL))
            edit.txt<-paste(edit.txt,Current.CTL[i],sep="\n")
      }
      dir.name<-strsplit(file.ctl,split="\\.")[[1]][1]
      temp<-strsplit(file.ctl,"\\\\")[[1]]
      file.name<-temp[length(temp)]    
      setwd(strsplit(file.ctl,split=file.name)[[1]][1])
      edit.win<-gwindow(temp[length(temp)])
      g<-ggroup(horizontal=FALSE,cont=edit.win)
     	tmp<-gframe("Editor",container=g)
      a<-gtext(edit.txt,width=600,height=300,
                font.attr=c(sizes="large",family="monospace"))
      add(tmp,a)
      save.b<-gbutton("Save",handler=function(h,...){ 
               write.table(svalue(a),file.ctl,quote=FALSE,row.names=FALSE,col.names=FALSE)})
      saveas.b<-gbutton("Save as",handler=function(h,...){ 
                 file.new<-gfile(text="Save as",type="save")
                 write.table(svalue(a),file.new,quote=FALSE,row.names=FALSE,col.names=FALSE)
                 dispose(edit.win)
                 Editor1(file.new)})
      gb<-ggroup(horizontal=TRUE,cont=g)
      tmp<-gframe("",container=gb)  
      add(tmp,save.b)
      tmp<-gframe("",container=gb)  
      add(tmp,saveas.b)
   }
      
   #### Edit with external editor ##############################################
   ExternalEditor<-function(h,...)
   {  control.path<-gfile(text="Choose control file",type="open")
      temp<-strsplit( Editor.path,split="\\\\")[[1]]
      Editor.name<-temp[length(temp)]
      Editor.path.t<-strsplit(Editor.path,split=Editor.name)[[1]]
      setwd(Editor.path.t)
      system(paste(Editor.name,control.path),invisible=F,wait=F)
   }
   
   ############################################################################# 
   # NONMEM run
   #############################################################################

   #### Run ####################################################################
   ###### Notes before run #####################################################
   BeforeRun<-function(h,...)
   {  gmessage("*** Notes before run ***
           \nThe file name of a control stream should be runnumber.ctl.
           \nID should be #ID.
           \nIPRED and IWRES should be defined in the control file.           
           \nIndividual predictions = IPRED  
           \nIndividual weighted residuals = IWRES 
           \nSemicolons used for comments in a NONMEM control stream should be preceded by one space as follows.
           \n  $THETA ; #8     (O)
           \n  $THETA; #8      (X)
           \n$TABLE statements to generate runnumber.eta (ID, ETA(1), ETA(2), ...) and runnumber.par (ID, V, CL, ...) should be included in a control stream.
           \nPrior to the beginning of a run,
           \n        $TABLE statements to generate runnumber.noh is inserted automatically in a control stream.
           \n        Control and data input files are copied automatically to a runnumber subfolder. 
           \n        The data input file copied to a runnumber subfolder is duplicated automatically to runnumber.csv.
           \nThe number of data records in several ouput files should be equal to that of data input file (ACCEPT argument in $DATA is not allowed).
           \nModel description should be written in English and should not include comma(,) separators",cont=T)
   }   
   
   ###### From default editor ##################################################
   Editor<-function(h,...)
   {  file.ctl<<-gfile(text="Choose NONMEM control file",type="open")
      Editor1<-function(file.ctl)
      {	 D.temp<-readLines(file.ctl)
      	 Current.CTL<-D.temp
      	 NONMEM.CTL<-D.temp
  	 if(is.null(Current.CTL))
      	 {  edit.txt<-" "
      	 } else
      	 {  edit.txt<-Current.CTL[1]
            for(i in 2:length(Current.CTL))
               edit.txt<-paste(edit.txt,Current.CTL[i],sep="\n")
      	 }
         dir.name<-strsplit(file.ctl,split="\\.")[[1]][1]
         temp<-strsplit(file.ctl,"\\\\")[[1]]
         file.name<-temp[length(temp)]    
         setwd(strsplit(file.ctl,split=file.name)[[1]][1])

      	 edit.win<-gwindow(temp[length(temp)])
      	 g<-ggroup(horizontal=FALSE,cont=edit.win)
     	 tmp<-gframe("Editor",container=g)
      	 a<-gtext(edit.txt,width=600,height=300,
                font.attr=c(sizes="large",family="monospace"))
      	 add(tmp,a)
      	 save.b<-gbutton("Save",handler=function(h,...){ 
           	write.table(svalue(a),file.ctl,quote=FALSE,row.names=FALSE,col.names=FALSE)})
      	 saveas.b<-gbutton("Save as",handler=function(h,...){ 
             	file.new<-gfile(text="Save As",type="save")
             	write.table(svalue(a),file.new,quote=FALSE,row.names=FALSE,col.names=FALSE)
             	dispose(edit.win)
             	Editor1(file.new)
             	file.ctl<<-file.new})
      	 gb<-ggroup(horizontal=TRUE,cont=g)
      	 tmp<-gframe("",container=gb)  
      	 add(tmp,save.b)
      	 tmp<-gframe("",container=gb)  
      	 add(tmp,saveas.b)
     	 RunNM.1<-gbutton("Run default NONMEM",handler=function(h,...){ 
                dir.create(dir.name,showWarnings=F)
                file.copy(file.ctl,dir.name,overwrite=T)
      	        temp<-strsplit(data.file,"\\\\")[[1]]
                data.name<-temp[length(temp)] 
                temp<-strsplit(dir.name,split="\\\\")[[1]]
             	file.id<-temp[length(temp)]
                file.copy(data.file,dir.name,overwrite=T)
                file.copy(data.file,paste(dir.name,"\\",file.id,".csv",sep=""),overwrite=T)      
                NONMEM.command<-paste(Default.NMpath, file.name, paste(strsplit(file.name,split="\\.")[[1]][1],"res",sep="."),
                                     paste(">",strsplit(file.name,split="\\.")[[1]][1],".console",sep=""))
                setwd(dir.name)
                param.num<-as.numeric(svalue(num.param))
                Description<-svalue(model.description)
                tt<-strsplit(Description,"\n")[[1]]
                Description<-""
                for(i in 1:length(tt))
                   Description<-paste(Description,tt[i])

                temp<-strsplit(Default.NMpath,split="\\\\")[[1]]
                NMversion<-ifelse(temp[length(temp)]=="nmfe7.bat","NM7","NM6")
                write.table(NMversion,"NM.version")

                RemakeCTL(paste(dir.name,"\\",file.id,".ctl",sep=""))
                write.table(" ",paste(dir.name,"\\",file.id,".console",sep=""))                
                system(NONMEM.command,wait=F) 
          	OpenNMConsole(file.id,paste(dir.name,"\\",file.id,".console",sep=""))
                alarm()
                TOT.temp<-TOT.RUN
                TOT.temp$num<-TOT.temp$num+1
                TOT.temp$data<-rbind(TOT.temp$data,c(file.id,dir.name,svalue(parent)))
                colnames(TOT.temp$data)<-c("ID","path","parents")
                TOT.RUN<<-TOT.temp
                OpenResult(file.id,paste(dir.name,"\\",file.id,".res",sep=""))
                D.LST<-readLines(paste(dir.name,"\\",file.id,".res",sep=""))
                if(sum(dir()==paste(file.id,".ETA",sep=""))!=0)
             	{   ETA<-read.table(paste(dir.name,"\\",file.id,".ETA",sep=""),skip=1,header=T)
                   temp.ETA<-colnames(ETA)
                   ET.id<-NULL
                   for(i in 1:length(temp.ETA))
                   {  if(length(strsplit(temp.ETA[i],"ET")[[1]])>1)
                         ET.id<-c(ET.id,i)
                   }
                   ETA<-ETA[,ET.id]       
                } else
                { ETA<-NULL
                }                        
                data.tempt<-read.csv(data.file,na.string=".")   
                if(sum(toupper(colnames(data.tempt))=="MDV")!=0)
                {  data.i<-which(toupper(colnames(data.tempt))=="MDV")
                   temp<-data.tempt[,data.i]
                   temp<-temp[temp==0]
                   data.n<-length(temp)         
                } else
                {  data.n<-nrow(data.tempt)
                }
                ShowResult1(D.LST,param.num,data.n,Description,ETA,file.id,dir.name)
         })
                
      	 RunNM.2<-gbutton("Run alternative NONMEM",handler=function(h,...){ 
                dir.create(dir.name,showWarnings=F)
                file.copy(file.ctl,dir.name,overwrite=T)
      	        temp<-strsplit(data.file,"\\\\")[[1]]
                data.name<-temp[length(temp)] 
                temp<-strsplit(dir.name,split="\\\\")[[1]]
            	file.id<-temp[length(temp)]
                file.copy(data.file,dir.name,overwrite=T)
                file.copy(data.file,paste(dir.name,"\\",file.id,".csv",sep=""),overwrite=T)
                
                NONMEM.command<-paste(Alternative.NMpath, file.name, paste(strsplit(file.name,split="\\.")[[1]][1],"res",sep="."),
                                     paste(">",strsplit(file.name,split="\\.")[[1]][1],".console",sep=""))
                setwd(dir.name)
                param.num<-as.numeric(svalue(num.param))
                Description<-svalue(model.description)
                tt<-strsplit(Description,"\n")[[1]]
                Description<-""
                for(i in 1:length(tt))
                   Description<-paste(Description,tt[i])
                
                temp<-strsplit(Alternative.NMpath,split="\\\\")[[1]]
                NMversion<-ifelse(temp[length(temp)]=="nmfe7.bat","NM7","NM6")
                write.table(NMversion,"NM.version")

                RemakeCTL(paste(dir.name,"\\",file.id,".ctl",sep=""))
                write.table(" ",paste(dir.name,"\\",file.id,".console",sep=""))                
                system(NONMEM.command,wait=F) 
             	OpenNMConsole(file.id,paste(dir.name,"\\",file.id,".console",sep=""))
                alarm()
                TOT.temp<-TOT.RUN
                TOT.temp$num<-TOT.temp$num+1
                TOT.temp$data<-rbind(TOT.temp$data,c(file.id,dir.name,svalue(parent)))
                colnames(TOT.temp$data)<-c("ID","path","parents")
                TOT.RUN<<-TOT.temp
                OpenResult(file.id,paste(dir.name,"\\",file.id,".res",sep=""))
                D.LST<-readLines(paste(dir.name,"\\",file.id,".res",sep=""))
                if(sum(dir()==paste(file.id,".ETA",sep=""))!=0)
             	{  ETA<-read.table(paste(dir.name,"\\",file.id,".ETA",sep=""),skip=1,header=T)
                   temp.ETA<-colnames(ETA)
                   ET.id<-NULL
                   for(i in 1:length(temp.ETA))
                   {  if(length(strsplit(temp.ETA[i],"ET")[[1]])>1)
                         ET.id<-c(ET.id,i)
                   }
                   ETA<-ETA[,ET.id]                              
                } else
                {  ETA<-NULL
                }   
                data.tempt<-read.csv(data.file,na.string=".")   
                if(sum(toupper(colnames(data.tempt))=="MDV")!=0)
                {  data.i<-which(toupper(colnames(data.tempt))=="MDV")
                   temp<-data.tempt[,data.i]
                   temp<-temp[temp==0]
                   data.n<-length(temp)         
                } else
                {  data.n<-nrow(data.tempt)
                }  
                ShowResult1(D.LST,param.num,data.n,Description,ETA,file.id,dir.name)               
         })  
                  
         openrundata<-function(h,...)
         {  data.file<<-gfile(text="Choose data file",type="open")
            svalue(data.set)<-data.file
         }

         grunNM<-ggroup(horizontal=F,cont=g) 
         tmp<-gframe("Model description",cont=grunNM)
         model.description<-gtext(" ",width=250,height=20)
         add(tmp,model.description)
         tmp<-gframe("Choose parents",cont=grunNM)
         item<-unique(c("ROOT",TOT.RUN$data[,"ID"]))
         parent<-gdroplist(item)
         add(tmp,parent)         
         tmp<-gframe("Number of parameters\n # of theta + # of omega and sigma(except fixed)",cont=grunNM)
         num.param<-gtext("",width=100,height=15)
         add(tmp,num.param)         
         tmp<-gframe("NONMEM data file",cont=grunNM)
         data.set<-gtext(" ",width=400,height=30)
         add(tmp,data.set)
         open.but<-gbutton("Open",handler=openrundata)
         add(tmp,open.but)                         
     	 gc<-ggroup(horizontal=TRUE,cont=g)
      	 tmp<-gframe("",container=gc)  
      	 add(tmp,RunNM.1)
      	 tmp<-gframe("",container=gc)  
      	 add(tmp,RunNM.2)
      }
      Editor1(file.ctl)
   }

   ###### From default editor -additional function##############################
   RemakeCTL<-function(file.name)
   {  
      D.CTL<-readLines(file.name)
      temp<-strsplit(file.name,split="\\\\")[[1]]
      file.id<-strsplit(temp[length(temp)],split="\\.")[[1]][1]
      D.temp<-matrix(D.CTL)
      NMversion<-read.table("NM.version")
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      if(sum(indicator==";runnumber.ctl",na.rm=T)==0)
      {  table.id<-which(tolower(indicator)=="$table")[1]
         add.CTL1<-";runnumber.ctl"
         if(NMversion=="NM7")
         {  add.CTL2<-paste("$TABLE ID TIME DV IPRED IRES IWRES NPRED NRES NWRES PREDI RESI WRESI", 
                   " CPRED CRES CWRES CPREDI CRESI CWRESI EPRED ERES EWRES NPDE FILE=",file.id,".NOH NOPRINT ONEHEADER",sep="")
         } else
         {  add.CTL2<-paste("$TABLE ID TIME DV IPRED IRES IWRES PRED RES WRES FILE=",file.id,".NOH NOPRINT ONEHEADER",sep="")
         }                    
         D.new<-c(D.CTL[1:(table.id-1)],add.CTL1,add.CTL2,D.CTL[table.id:length(D.CTL)])
      } else
      {  D.new<-D.CTL
      }  
      write.table(D.new,file.name,quote=FALSE,row.names=FALSE,col.names=FALSE)
   }

   OpenResult<-function(NonmemRunID,filename)
   {  NonmemRes<-gwindow(paste(as.character(NonmemRunID),".res",sep=""))
      a<-gtext("",cont=NonmemRes,font.attr=c(family="korea1deb"),width=700)
      old<-0
      console.text<-readLines(filename)
      svalue(a)<-console.text
   }  
 
   OpenNMConsole<-function(NonmemRunID,filename)
   {  for(k in 1:10000)
      { n<-1+1
      }
      NonmemConsole<-gwindow(as.character(NonmemRunID))
      a<-gtext("",cont=NonmemConsole,font.attr=c(family="korea1deb"))
      diff.t<-TRUE
      old<-0
      console.text<-readLines(filename)
      svalue(a)<-console.text
      while(diff.t)
      {  new<-length(console.text)
         if(new!=old & new!=0) add(a,console.text[(old+1):new])
         old<-new
         diff.t<-ifelse(new>2, console.text[new-1]!="OUTPUT",TRUE)
         console.text<-readLines(filename)
      } 
   } 

   CatchNMversion<-function(D.LST)
   {  D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      NM.id<- which(indicator=="1NONLINEAR")
      temp<- strsplit(D.LST[NM.id],split=" ")[[1]]
      NMversion<-ifelse(temp[which(temp=="VERSION")+1]=="VI",6,7)
      return(NMversion)
   }

   ShowResult1<-function(D.LST,param.num,data.n,Description,ETA,file.id,dir.name)
   {  RunID<-TOT.RUN$data[TOT.RUN$num,1]
      n.lst<-length(D.LST)
      Date<-D.LST[n.lst-1]
      Time<-D.LST[n.lst]

      D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      min.id<- which(indicator=="0MINIMIZATION")
      min.id<-min.id[length(min.id)]
      Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]

      indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])

      cond.id<-grep(" EIGENVALUES ",D.LST)
      if(length(cond.id)!=0)
      {  
         flag<-T
         id.current<-cond.id+5
         cond.line<-0
         while(flag)
         {  ttemp<-D.LST[id.current]
            flag<-ttemp!=" "
            if(flag)
            {  cond.line<-cond.line+1
               id.current<-id.current+1
            }
         }
         temp<-NULL
         for(i in 1:cond.line)
            temp<-c(temp,strsplit(D.LST[cond.id+5+cond.line+i],split=" ")[[1]])
         temp<-as.numeric(temp[temp!=""&temp!="+"])
         cond.num<-round(max(temp)/min(temp),3)        
      } else
      {  cond.num<-NA
      }   
      
      obj.id<-which(indicator==" #OBJV")
      if(length(obj.id)!=0)
      { obj.id<-obj.id[length(obj.id)]
      } else
      { obj.id<-9+which(indicator==" ********************                           MINIMUM VALUE OF OBJECTIVE FUNCTION                  ********************" )
      } 
      temp<-strsplit(D.LST[obj.id],split=" ")[[1]]
      temp<-temp[3:(length(temp)-3)]
      temp<-as.numeric(temp[temp!=""])
      Obj<-temp[!is.na(temp)]
      AIC<-Obj+2*param.num
      AICc<-round(Obj+2*param.num+2*param.num*(param.num+1)/(data.n-param.num-1),3)
      SBC<-round(Obj+param.num*log(data.n),3)
      parent<-TOT.RUN$data[TOT.RUN$num,"parents"]

# choose THETA,seTHETA,OMEGA,seOMEGA,SIGMA,seSIGMA

      final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
      final.start.id<-final.start.id[length(final.start.id)]
      Result.LST<-D.LST[final.start.id:length(D.LST)]

      theta.id<-grep("THETA",Result.LST)
      theta.line<-0
      theta.flag<-TRUE
      while(theta.flag)
      {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
         {  theta.line<-theta.line+1
         } else
         {  theta.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:theta.line)
         temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
      temp<-as.numeric(temp[temp!=""])
      THETA<-temp
      
      seTHETA<-rep(NA,length(THETA))  
      if(length(theta.id)!=1)
      {  temp<-NULL
         for(i in 1:theta.line)
            temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
         temp<-as.numeric(temp[temp!=""])
         seTHETA<-temp
      } 

      omega.id<-grep("OMEGA",Result.LST)
      omega.line<-0
      omega.flag<-TRUE
      while(omega.flag)
      {  if(Result.LST[omega.id[1]+3+omega.line]!=" ")
         {  omega.line<-omega.line+1
         } else
         {  omega.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:omega.line)
         temp<-c(temp,unlist(strsplit(Result.LST[omega.id[1]+2+i],split=" ")))
      temp<-temp[temp!=""]

      N.eta<-length(temp)
      OMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      seOMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      omega.name<-NULL
      id.current<-omega.id[1]+4+omega.line  
      id.current1<-omega.id[2]+4+omega.line  
      
      for(i in 1:N.eta)
      {  temp<-NULL;temp1<-NULL
         flag<-TRUE
         while(flag)
         {  id.current<-id.current+1;id.current1<-id.current1+1
            flag<-Result.LST[id.current]!=" "
            if(flag)
            {  temp<-c(temp,unlist(strsplit(Result.LST[id.current],split="+ ")))
               temp1<-c(temp1,unlist(strsplit(Result.LST[id.current1],split="+ ")))            
            }
         }
         temp<-temp[temp!=""];temp1<-temp1[temp1!=""]
         temp<-temp[temp!="+"] ;temp1<-temp1[temp1!="+"] 
         temp1[ temp1=="........."]<-NA 
         
         for(j in 1:N.eta)
         {  OMEGA[i,j]<-as.numeric(temp[j])
            seOMEGA[i,j]<-as.numeric(temp1[j])
            omega.name<-c(omega.name,paste("OMEGA(",i,"/",j,")",sep=""))
         }  
         id.current<-id.current+1 
         id.current1<-id.current1+1              
      }     

      sigma.id<-grep("SIGMA",Result.LST)
      temp<-unlist(strsplit(Result.LST[sigma.id[1]+6],split="  "))
      temp<-temp[temp!=""]; temp<-temp[temp!="+"]
      SIGMA<-as.numeric(temp)
      seSIGMA<-rep(NA,length(SIGMA))
      if(length(theta.id)!=1)
      {  temp<-unlist(strsplit(Result.LST[sigma.id[2]+6],split="  "))
         temp<-temp[temp!=""]; temp<-temp[temp!="+"]
         seSIGMA<-as.numeric(temp)      
      } 

      names.est<-c(paste("TH",1:length(THETA),sep=""),omega.name,paste("SIGMA",1:length(SIGMA)))      
      EST<-c(THETA,OMEGA,SIGMA)
      SE<-c(seTHETA,seOMEGA,seSIGMA)

      RSE<-round(SE/EST*100,4)
      Lower<-round(EST-1.96*SE,4)
      Upper<-round(EST+1.96*SE,4)

      if(sum(!is.na(seTHETA))==0)
      {  COV<-"NONE"
         cond.num<-NA
      } else
      {  COV<-"OK"
      }

      temp<-c(RunID,Date,Time,Min,COV,Obj,AIC,AICc,SBC,cond.num,parent,Description,param.num)
      if(TOT.RUN$num>2)
      {  run.table[]<-rbind(run.table[],temp)
      } else
      {  run.table[][TOT.RUN$num,]<-temp
      }
      shrinkage.ETA<-matrix(NA,ncol=N.eta,nrow=N.eta)    
      if(!is.null(ETA))
      {  se.ETA<-apply(ETA,2,sd) 
         diag(shrinkage.ETA)<-(1-se.ETA/sqrt(diag(OMEGA)))*100
      }
      shrinkage<-c(rep("NA",length(THETA)),round(shrinkage.ETA,3),rep("NA",length(SIGMA)))   
      CV.ETA<-matrix(NA,ncol=N.eta,nrow=N.eta)
      diag(CV.ETA)<-sqrt(diag(OMEGA))*100
      CV<-c(rep("NA",length(THETA)),round(CV.ETA,3),rep("NA",length(SIGMA)))
      
      tot.res<-cbind(names.est,EST,SE,RSE,Lower,Upper,shrinkage,CV)
      colnames(tot.res)<-c("Parameters","Estimates","SE","%RSE","Lower","Upper","%Shrinkage","%CV")
      tot.res[is.na(tot.res)| is.nan(tot.res) | tot.res=="NA"| tot.res=="NaN"]<-" "
      tot.res<-tot.res[-which(apply(tot.res,1,function(x) sum(x==" "))==7),]
      gtable(tot.res, cont=gwindow(paste(file.id,".sum",sep="")),do.subset=TRUE,width=150)
      write.csv(tot.res,paste(dir.name,"\\",file.id,".sum",sep=""),quote=F)
      write.csv(tot.res,paste(dir.name,"\\",file.id,".sum.csv",sep=""),quote=F)      
   }
      
   ###### From external editor #################################################

   ExternalRun<-function(h,...)
   {  opendata<-function(h,...)
      {  data.file<<-gfile(text="Open data file",type="open")
         svalue(data.t)<-data.file      
      }

      openControl<-function(h,...)
      {  control.file<<-gfile(text="Open control file (runnumber.ctl)",type="open")
         svalue(control.t)<-control.file
         dir.name<-strsplit(control.file,split="\\.")[[1]][1]
         svalue(dir.t)<-dir.name

         temp<-strsplit(dir.name,split="\\\\")[[1]]
         file.id<-temp[length(temp)]
         setwd(strsplit(dir.name,split=file.id)[[1]])
       }

      openEdtRun<-function(h,...)
      {  dir.name<-strsplit(control.file,split="\\.")[[1]][1]
         dir.create(dir.name,showWarnings=F)
         temp<-strsplit(control.file,split="\\\\")[[1]]
         file.id<-temp[length(temp)]
         file.id1<-strsplit(file.id,split="\\.")[[1]][1]
         setwd(dir.name)
         file.copy(control.file,dir.name,overwrite=T)
         file.copy(data.file,dir.name,overwrite=T)
         file.copy(data.file,paste(dir.name,"\\",file.id1,".csv",sep=""),overwrite=T)
         
         temp<-strsplit(Default.NMpath,split="\\\\")[[1]]
         NMversion<-ifelse(temp[length(temp)]=="nmfe7.bat","NM7","NM6")
         write.table(NMversion,"NM.version")
         RemakeCTL(file.id)

         temp<-strsplit( Editor.path,split="\\\\")[[1]]
         Editor.name<-temp[length(temp)]
         Editor.path.t<-strsplit(Editor.path,split=Editor.name)[[1]]
         setwd(Editor.path.t)
         system(paste(Editor.name,paste(dir.name,"\\",file.id,sep="")),wait=F)
         setwd(dir.name)
      }
 
      Outerwin<-gwindow("Run from external editor")
      Bgroup<-ggroup(cont=Outerwin,horizontal=TRUE)
      BBgroup<-ggroup(cont=Bgroup,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit("",width=50)
      button1<-gbutton("Open control file",handler=openControl)
      add(tmp,button1)
      add(tmp,control.t)

      tmp<-gframe("",cont=BBgroup)
      button2<-gbutton("Open data files",handler=opendata,width=20,height=10)
      data.t<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,data.t) 
      tmp<-gframe("NONMEM run number directory",cont=BBgroup)
      dir.t<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,dir.t) 
 
      tmp<-gframe("",cont=BBgroup)
      button3<-gbutton("Open external editor for run",handler=openEdtRun,width=20,height=10)
      add(tmp,button3)   
   }
      
   ###### Direct run ###########################################################
   DirectRun<-function(h,...)
   {  RunNONMEM<-function(id)
      {  i<-id
         data.file<-tclvalue(DataFile.Name[[i]])
         file.ctl<-tclvalue(ControlFile.Name[[i]])
         param.num<-as.numeric(tclvalue(Param.Num[[i]]))
         Description<-tclvalue(Description.N[[i]])
         parents<-toupper(tclvalue(Parent.Num[[i]]) )
         dir.name<-strsplit(file.ctl,split="\\.")[[1]][1]
         temp<-strsplit(file.ctl,"\\\\")[[1]]
         file.name<-temp[length(temp)]              
         dir.create(dir.name,showWarnings=F)
         file.copy(file.ctl,paste(dir.name,file.name,sep="\\"),overwrite=T)
         temp<-strsplit(data.file,"\\\\")[[1]]
         data.name<-temp[length(temp)] 
         file.id<-strsplit(data.name,"\\.")[[1]][1]
         file.copy(data.file,paste(dir.name,data.name,sep="\\"),overwrite=T)
         file.copy(data.file,paste(dir.name,"\\",file.id,".csv",sep=""),overwrite=T)

         NONMEM.command<-paste(Default.NMpath, file.name, paste(strsplit(file.name,split="\\.")[[1]][1],"res",sep="."),
                               paste(">",strsplit(file.name,split="\\.")[[1]][1],".console",sep=""))
         setwd(dir.name)
         temp<-strsplit(Default.NMpath,split="\\\\")[[1]]
         NMversion<-ifelse(temp[length(temp)]=="nmfe7.bat","NM7","NM6")
         write.table(NMversion,"NM.version")

         RemakeCTL(paste(dir.name,file.name,sep="\\"))
         write.table(" ",paste(dir.name,"\\",file.id,".console",sep=""))
         system(NONMEM.command,invisible=F,show.output.on.console=F,wait=F) 
         OpenNMConsole(file.id,paste(dir.name,"\\",file.id,".console",sep=""))
         alarm()
         temp<-strsplit(dir.name,split="\\\\")[[1]]
         file.id<-temp[length(temp)]
         TOT.temp<-TOT.RUN
         TOT.temp$num<-TOT.temp$num+1
         TOT.temp$data<-rbind(TOT.temp$data,c(file.id,dir.name,parents))
         colnames(TOT.temp$data)<-c("ID","path","parents")
         TOT.RUN<<-TOT.temp

         OpenResult(file.id,paste(dir.name,"\\",file.id,".res",sep=""))
         D.LST<-readLines(paste(dir.name,"\\",file.id,".res",sep=""))
         if(sum(dir()==paste(file.id,".ETA",sep=""))!=0)
         {  ETA<-read.table(paste(dir.name,"\\",file.id,".ETA",sep=""),skip=1,header=T)
            temp.ETA<-colnames(ETA)
            ET.id<-NULL
            for(i in 1:length(temp.ETA))
            {  if(length(strsplit(temp.ETA[i],"ET")[[1]])>1)
                  ET.id<-c(ET.id,i)
            }
            ETA<-ETA[,ET.id]           
         } else
         {  ETA<-NULL
         }   
         data.tempt<-read.csv(data.file,na.string=".")   
         if(sum(toupper(colnames(data.tempt))=="MDV")!=0)
         {  data.i<-which(toupper(colnames(data.tempt))=="MDV")
            temp<-data.tempt[,data.i]
            temp<-temp[temp==0]
            data.n<-length(temp)         
         } else
         {  data.n<-nrow(data.tempt)
         }
         ShowResult1(D.LST,param.num=param.num,data.n,Description=Description,ETA,file.id,dir.name)               
          
      }

      SeqRun<-function()
      {  for(i in 1:k)
         {  data.file<-tclvalue(DataFile.Name[[i]])
            file.ctl<-tclvalue(ControlFile.Name[[i]])
            param.num<-as.numeric(tclvalue(Param.Num[[i]]))
            Description<-tclvalue(Description.N[[i]])
            parents<-toupper(tclvalue(Parent.Num[[i]]) )

            dir.name<-strsplit(file.ctl,split="\\.")[[1]][1]
            temp<-strsplit(file.ctl,"\\\\")[[1]]
            file.name<-temp[length(temp)]              
            dir.create(dir.name,showWarnings=F)
            file.copy(file.ctl,paste(dir.name,file.name,sep="\\"),overwrite=T)
            temp<-strsplit(data.file,"\\\\")[[1]]
            data.name<-temp[length(temp)] 
            file.copy(data.file,paste(dir.name,data.name,sep="\\"),overwrite=T)
            file.id<-strsplit(file.name,split="\\.")[[1]][1]
            file.copy(data.file,paste(dir.name,"\\",file.id,".csv",sep=""),overwrite=T)
            NONMEM.command<-paste(Default.NMpath, paste(file.id,".ctl",sep=""), paste(file.id,"res",sep="."),
                                     paste(">",file.id,".console",sep=""))
            setwd(dir.name)
            temp<-strsplit(Default.NMpath,split="\\\\")[[1]]
            NMversion<-ifelse(temp[length(temp)]=="nmfe7.bat","NM7","NM6")
            write.table(NMversion,"NM.version")

            RemakeCTL(paste(dir.name,file.name,sep="\\"))
            write.table(" ",paste(dir.name,"\\",file.id,".console",sep=""))            
            system(NONMEM.command,invisible=F,show.output.on.console=F,wait=F)
            OpenNMConsole(file.id,paste(dir.name,"\\",file.id,".console",sep=""))           
            alarm()
            temp<-strsplit(dir.name,split="\\\\")[[1]]
            file.id<-temp[length(temp)]
            TOT.temp<-TOT.RUN
            TOT.temp$num<-TOT.temp$num+1
            TOT.temp$data<-rbind(TOT.temp$data,c(file.id,dir.name,parents))
            colnames(TOT.temp$data)<-c("ID","path","parents")
            TOT.RUN<<-TOT.temp

            OpenResult(file.id,paste(dir.name,"\\",file.id,".res",sep=""))      
            D.LST<-readLines(paste(dir.name,"\\",file.id,".res",sep=""))
            if(sum(dir()==paste(file.id,".ETA",sep=""))!=0)
            {  ETA<-read.table(paste(dir.name,"\\",file.id,".ETA",sep=""),skip=1,header=T)
               temp.ETA<-colnames(ETA)
               ET.id<-NULL
               for(i in 1:length(temp.ETA))
               {  if(length(strsplit(temp.ETA[i],"ET")[[1]])>1)
                    ET.id<-c(ET.id,i)
               }
               ETA<-ETA[,ET.id]              
            } else
            {  ETA<-NULL
            }   
            data.tempt<-read.csv(data.file,na.string=".")   
            if(sum(toupper(colnames(data.tempt))=="MDV")!=0)
            {  data.i<-which(toupper(colnames(data.tempt))=="MDV")
               temp<-data.tempt[,data.i]
               temp<-temp[temp==0]
               data.n<-length(temp)         
            } else
            {  data.n<-nrow(data.tempt)
            }
            ShowResult1(D.LST,param.num=param.num,data.n,Description=Description,ETA,file.id,dir.name)                                       
         }
      }

      Add<-function()
      {  AddLine(k+1)
      }

      AddLine<-function(id)
      {  k<<-id
         ControlFile.Name[[id]]<<-tclVar("")
         ControlFileName[[id]]<<-tkentry(tt,width="20",textvariable=ControlFile.Name[[id]])
         DataFile.Name[[id]]<<-tclVar("")
         DataFileName[[id]]<<-tkentry(tt,width="20",textvariable=DataFile.Name[[id]])
         Run.Num[[id]]<<-tclVar("")
         RunNum[[id]]<<-tkentry(tt,width="15",textvariable=Run.Num[[id]])
         Description.N[[id]]<<-tclVar("")
         DescriptionN[[id]]<<-tkentry(tt,width="60",textvariable=Description.N[[id]])
         Param.Num[[id]]<<-tclVar("")
         ParamNum[[id]]<<-tkentry(tt,width="15",textvariable=Param.Num[[id]])
         Parent.Num[[id]]<<-tclVar("")
         ParentNum[[id]]<<-tkentry(tt,width="15",textvariable=Parent.Num[[id]])
         tkgrid(ControlFileName[[id]],RunNum[[id]],tklabel(tt,text=" "),ParentNum[[id]],
                   tklabel(tt,text=" "),ParamNum[[id]], tklabel(tt,text=" "),
                      DescriptionN[[id]],DataFileName[[id]],tklabel(tt,text=" "))
         tkgrid(tt)
      }

      ConFile<-function()
      {  kk<-k
         tclvalue(ControlFile.Name[[kk]])<<-gfile(text="Open control file (runnumber.ctl)",type="open")
         file.ctl<-tclvalue(ControlFile.Name[[kk]])
         temp<-strsplit(file.ctl,"\\\\")[[1]]
         RunNumber<-strsplit(tolower(temp[length(temp)]),split="\\.ctl")[[1]][1]  
         tclvalue(Run.Num[[kk]])<<-RunNumber
         DirectRunNum<<-c(DirectRunNum,RunNumber)
      }

      DataFile<-function()
      {  kk<-k
         file.ctl<-tclvalue(ControlFile.Name[[kk]])
         temp<-strsplit(file.ctl,"\\\\")[[1]]
         file.name<-temp[length(temp)]
         dir.name<-strsplit(file.ctl,split=file.name)[[1]] 
         setwd(dir.name)
         tclvalue(DataFile.Name[[kk]])<<-gfile(text="Data file",type="open")
      }
      
      SaveRun<-function(h,...)
      {  TOT.table<-NULL
         for(i in 1:k)
         {  data.file<-tclvalue(DataFile.Name[[i]])
            file.ctl<-tclvalue(ControlFile.Name[[i]])
            param.num<-tclvalue(Param.Num[[i]])
            Description<-tclvalue(Description.N[[i]])
            parents<-toupper(tclvalue(Parent.Num[[i]]))
            Runnum<-tclvalue(Run.Num[[i]])
            TOT.table<-rbind(TOT.table,c(file.ctl,Runnum,parents,param.num,Description,data.file))
         }  
         colnames(TOT.table)<-c("ControlFile","Runnumber","Parents","paramnum","Description","DataFile")
         run.file<-paste(gfile(text="Save as csv",type="save"),".csv",sep="")
         write.csv(TOT.table,run.file)
      }
      
      DirectRunNum<<-NULL
      Toptt<<-tktoplevel()
      tkwm.title(Toptt,"Direct run")
      tt<-tkframe(Toptt)
      ttg<-tkframe(tt)
      OK.but3 <-tkbutton(ttg,text="Sequential runs",command=SeqRun)
      OK.but4 <-tkbutton(ttg,text="Save as csv",command=SaveRun)
      tkgrid(tklabel(ttg,text=""),tklabel(ttg,text=""),tklabel(ttg,text=""),
                     OK.but3,tklabel(ttg,text=""),OK.but4)
      tkgrid(ttg)
      OK.but1 <-tkbutton(tt,text="Control file (runnumber.ctl)",command=ConFile)
      OK.but2 <-tkbutton(tt,text="Data file",command=DataFile)
      Add.but<-tkbutton(tt,text="Add",command=Add)
      
      tkgrid(OK.but1,tklabel(tt,text="Run number"),tklabel(tt,text=""),
                    tklabel(tt,text="Parents"),tklabel(tt,text=""),
                    tklabel(tt,text="# of parameters"),tklabel(tt,text=""),tklabel(tt,text="Description (English only)"),OK.but2,Add.but)

      tkgrid(tt)
      ControlFile.Name<<-list()
      ControlFileName<<-list()
      DataFile.Name<<-list()
      DataFileName<<-list()
      RunNum<<-list()
      Run.Num<<-list()
      DescriptionN<<-list()
      Description.N<<-list()
      ParamNum<<-list()
      Param.Num<<-list()      
      ParentNum<<-list()
      Parent.Num<<-list()
            
      k<<-1
      ControlFile.Name[[k]]<<-tclVar("")
      ControlFileName[[k]]<<-tkentry(tt,width="20",textvariable=ControlFile.Name[[k]])
      DataFile.Name[[k]]<<-tclVar("")
      DataFileName[[k]]<<-tkentry(tt,width="20",textvariable=DataFile.Name[[k]]) 
      Run.Num[[k]]<<-tclVar("")
      RunNum[[k]]<<-tkentry(tt,width="15",textvariable=Run.Num[[k]])
      Description.N[[k]]<<-tclVar("")
      DescriptionN[[k]]<<-tkentry(tt,width="60",textvariable=Description.N[[k]])
      Param.Num[[k]]<<-tclVar("")
      ParamNum[[k]]<<-tkentry(tt,width="15",textvariable=Param.Num[[k]])
      Parent.Num[[k]]<<-tclVar("")
      ParentNum[[k]]<<-tkentry(tt,width="15",textvariable=Parent.Num[[k]])

      tkgrid(ControlFileName[[k]],RunNum[[k]],tklabel(tt,text=" "),ParentNum[[k]],
                   tklabel(tt,text=" "),ParamNum[[k]], tklabel(tt,text=" "),
                      DescriptionN[[k]],DataFileName[[k]],tklabel(tt,text=" "))
      tkgrid(tt)	
   }
      
   #### Run table ##############################################################
   ###### Make run table from runnumber subfolders ############################
   AddRunTable<-function(h,...)
   {  AddTable<-function(h,...)
      {  for(ii in 1:kk)
         {  dir.name<-tclvalue(ControlFile.Dir[[ii]])
            file.id<-tclvalue(Run.Num[[ii]])
            param.num<-as.numeric(tclvalue(Param.Num[[ii]]))
            Description<-tclvalue(Description.N[[ii]])
            TOT.temp<-TOT.RUN
            TOT.temp$num<-TOT.temp$num+1
            TOT.temp$data<-rbind(TOT.temp$data,c(file.id,dir.name,toupper(tclvalue(Parent.Num[[ii]]))))
            colnames(TOT.temp$data)<-c("ID","path","parents")
            TOT.RUN<<-TOT.temp

            D.LST<-readLines(paste(dir.name,"\\",file.id,".res",sep=""))
            ETA<-read.table(paste(dir.name,"\\",file.id,".ETA",sep=""),skip=1,header=T)
            temp.ETA<-colnames(ETA)
            ET.id<-NULL
            for(i in 1:length(temp.ETA))
            {  if(length(strsplit(temp.ETA[i],"ET")[[1]])>1)
                  ET.id<-c(ET.id,i)
            }
            ETA<-ETA[,ET.id]              
            data.file<-paste(dir.name,"\\",file.id,".csv",sep="")
            data.tempt<-read.csv(data.file,na.string=".")   
            if(sum(toupper(colnames(data.tempt))=="MDV")!=0)
            {  data.i<-which(toupper(colnames(data.tempt))=="MDV")
               temp<-data.tempt[,data.i]
               temp<-temp[temp==0]
               data.n<-length(temp)         
            } else
            {  data.n<-nrow(data.tempt)
            }

            RunID<-TOT.RUN$data[TOT.RUN$num,1]
            n.lst<-length(D.LST)
            Date<-D.LST[n.lst-1]
            Time<-D.LST[n.lst]

            D.temp<-matrix(D.LST)
            indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
            min.id<- which(indicator=="0MINIMIZATION")
            min.id<-min.id[length(min.id)]
            Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]

            indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])
            cond.id<-grep(" EIGENVALUES ",D.LST)
            if(length(cond.id)!=0)
            {  flag<-T
               id.current<-cond.id+5
               cond.line<-0
               while(flag)
               {  ttemp<-D.LST[id.current]
                  flag<-ttemp!=" "
                  if(flag)
                  {  cond.line<-cond.line+1
                     id.current<-id.current+1
                  }
               }
               temp<-NULL
               for(i in 1:cond.line)
                  temp<-c(temp,strsplit(D.LST[cond.id+5+cond.line+i],split=" ")[[1]])
               temp<-as.numeric(temp[temp!=""&temp!="+"])
               cond.num<-round(max(temp)/min(temp),3)        
            } else
            {  cond.num<-NA
            }
            
            obj.id<-which(indicator==" #OBJV")
            if(length(obj.id)!=0)
            {  obj.id<-obj.id[length(obj.id)]
            } else
            {  obj.id<-9+which(indicator==" ********************                           MINIMUM VALUE OF OBJECTIVE FUNCTION                  ********************" )
            } 
            temp<-strsplit(D.LST[obj.id],split=" ")[[1]]
            temp<-temp[3:(length(temp)-3)]
            temp<-as.numeric(temp[temp!=""])
            Obj<-temp[!is.na(temp)]
            AIC<-Obj+2*param.num
            AICc<-round(Obj+2*param.num+2*param.num*(param.num+1)/(data.n-param.num-1),3)
            SBC<-round(Obj+param.num*log(data.n),3)
            parent<-TOT.RUN$data[TOT.RUN$num,"parents"]

            final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
            final.start.id<-final.start.id[length(final.start.id)]
            Result.LST<-D.LST[final.start.id:length(D.LST)]

            theta.id<-grep("THETA",Result.LST)
            theta.line<-0
            theta.flag<-T
            while(theta.flag)
            {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
               {  theta.line<-theta.line+1
               } else
               {  theta.flag<-FALSE
               }  
            }
            temp<-NULL
            for(i in 1:theta.line)
              temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
            temp<-as.numeric(temp[temp!=""])
            THETA<-temp
      
            seTHETA<-rep(NA,length(THETA))  
            if(length(theta.id)!=1)
            {  temp<-NULL
               for(i in 1:theta.line)
                  temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
               temp<-as.numeric(temp[temp!=""])
               seTHETA<-temp
            }     

            if(sum(!is.na(seTHETA))==0)
            {  COV<-"NONE"
               cond.num<-NA
            } else
            {  COV<-"OK"
            }

            temp<-c(RunID,Date,Time,Min,COV,Obj,AIC,AICc,SBC,cond.num,parent,Description,param.num)

            if(TOT.RUN$num>2)
            {  run.table[]<-rbind(run.table[],temp)
            } else
            {  run.table[][TOT.RUN$num,]<-temp
            }
         } 
         tkdestroy(Toptt)
      }

      Add<-function()
      {  AddLine(kk+1)
      }

      AddLine<-function(id)
      {  kk<<-id
         ControlFile.Dir[[id]]<<-tclVar("")
         ControlFileDir[[id]]<<-tkentry(tt,width="20",textvariable=ControlFile.Dir[[id]])
         Run.Num[[id]]<<-tclVar("")
         RunNum[[id]]<<-tkentry(tt,width="15",textvariable=Run.Num[[id]])
         Description.N[[id]]<<-tclVar("")
         DescriptionN[[id]]<<-tkentry(tt,width="60",textvariable=Description.N[[id]])
         Param.Num[[id]]<<-tclVar("")
         ParamNum[[id]]<<-tkentry(tt,width="15",textvariable=Param.Num[[id]])
         Parent.Num[[id]]<<-tclVar("")
         ParentNum[[id]]<<-tkentry(tt,width="15",textvariable=Parent.Num[[id]])

         tkgrid(ControlFileDir[[id]],RunNum[[id]],tklabel(tt,text=" "),ParentNum[[id]],
                   tklabel(tt,text=" "),ParamNum[[id]], tklabel(tt,text=" "),
                      DescriptionN[[id]],tklabel(tt,text=" "))
         tkgrid(tt)
      }

      ConFile<-function()
      {  kkk<-kk
         tclvalue(ControlFile.Dir[[kkk]])<<-gfile(text="Choose run number subfolder",type="selectdir")
         file.ctl<-tclvalue(ControlFile.Dir[[kkk]])
         temp<-strsplit(file.ctl,"\\\\")[[1]]
         RunNumber<-temp[length(temp)]
         tclvalue(Run.Num[[kkk]])<<-RunNumber
         setwd(strsplit(file.ctl,split=RunNumber)[[1]])
      }

      Toptt<<-tktoplevel()
      tkwm.title(Toptt,"Make run table from runnumber subfolder")
      tt<-tkframe(Toptt)
      ttg<-tkframe(tt)
      OK.but3 <-tkbutton(ttg,text="OK",command=AddTable)

      tkgrid(tklabel(ttg,text=""),tklabel(ttg,text=""),tklabel(ttg,text=""),
                     OK.but3,tklabel(ttg,text=""))
      tkgrid(ttg)
      OK.but1 <-tkbutton(tt,text="Run number subfolder",command=ConFile)
      Add.but<-tkbutton(tt,text="Add",command=Add)
      
      tkgrid(OK.but1,tklabel(tt,text="Run number"),tklabel(tt,text=""),
                    tklabel(tt,text="Parents"),tklabel(tt,text=""),
                    tklabel(tt,text="# of parameters"),tklabel(tt,text=""),tklabel(tt,text="Description (English only)"),Add.but)
      tkgrid(tt)
      ControlFile.Dir<<-list()
      ControlFileDir<<-list()
      RunNum<<-list()
      Run.Num<<-list()
      DescriptionN<<-list()
      Description.N<<-list()
      ParamNum<<-list()
      Param.Num<<-list()      
      ParentNum<<-list()
      Parent.Num<<-list()
            
      kk<<-1
      ControlFile.Dir[[kk]]<<-tclVar("")
      ControlFileDir[[kk]]<<-tkentry(tt,width="20",textvariable=ControlFile.Dir[[kk]])
      Run.Num[[kk]]<<-tclVar("")
      RunNum[[kk]]<<-tkentry(tt,width="15",textvariable=Run.Num[[kk]])
      Description.N[[kk]]<<-tclVar("")
      DescriptionN[[kk]]<<-tkentry(tt,width="60",textvariable=Description.N[[kk]])
      Param.Num[[kk]]<<-tclVar("")
      ParamNum[[kk]]<<-tkentry(tt,width="15",textvariable=Param.Num[[kk]])
      Parent.Num[[kk]]<<-tclVar("")
      ParentNum[[kk]]<<-tkentry(tt,width="15",textvariable=Parent.Num[[kk]])

      tkgrid(ControlFileDir[[kk]],RunNum[[kk]],tklabel(tt,text=" "),ParentNum[[kk]],
                   tklabel(tt,text=" "),ParamNum[[kk]], tklabel(tt,text=" "),
                      DescriptionN[[kk]],tklabel(tt,text=" "))
      tkgrid(tt)
   }
      
   ###### Save run table as csv ################################################
   saveRUNTABLE.handler<-function(h,...)
   {  runtable<-run.table[]
      runtable<-runtable[runtable$Run!="",]
      temp.table<-cbind(runtable,TOT.RUN$data[,2])
      colnames(temp.table)<-c(colnames(runtable),"Path")
      write.csv(temp.table,paste(gfile(text="Save run table as csv",
           type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F,quote=F)
   }
    
   ###### Load run table #######################################################
   loadRUNTABLE.handler<-function(h,...)
   {  run.filename<-gfile("Open run table file",type="open")
      run.f<-read.csv(run.filename)

      for(ii in 1:nrow(run.f))
      {  temp<-as.character(run.f[ii,1:13])
         for(j in c(1,2,3,4,5,11,12))
            temp[j]<-as.character(run.f[ii,j])
         if(TOT.RUN$num>2)
         {  run.table[]<-rbind(run.table[],temp)
         } else
         {  run.table[][TOT.RUN$num+ii,]<-temp
         }
      }    
      TOT.temp<-TOT.RUN
      TOT.temp$data<-cbind(as.character(run.f$Run.number),as.character(run.f$Path),as.character(run.f$Parents))          
      colnames(TOT.temp$data)<-c("ID","path","parents")
      TOT.temp$data<-rbind(TOT.RUN$data,TOT.temp$data)
      TOT.temp$num<-TOT.RUN$num+nrow(run.f)   
      TOT.RUN<<-TOT.temp      
    }
   
   #### Model tree #############################################################
   Tree.handler<-function(h,...)
   {  Tree.make<-function(Tree.data,Root.name)
      {  id.tree<-Tree.data[,1]
         parent.tree<-Tree.data[,2]
         Tree.struct<-list()
         Tree.struct$ROOT<-NULL
         Tree.struct$Child<-list()
         for(i in 1:nrow(Tree.data))
         {  if(parent.tree[i]==Root.name)
            {  Tree.struct$ROOT<-c(Tree.struct$ROOT,id.tree[i])
            }
         }
         if(!is.null(Tree.struct$ROOT))
         {  for(j in 1:length(Tree.struct$ROOT))
            {  Tree.struct$Child[[j]]<-id.tree[which(parent.tree==Tree.struct$ROOT[j])]
               names(Tree.struct$Child)[j]<-paste(Root.name,Tree.struct$ROOT[j],sep="-")
            }
            return(Tree.struct)
         } else
         {  return(NULL)
         }
      }
 
      expandTree<-function(Tree.data,Root.list,TREE.display)
      {  empty.col<-ncol(TREE.display)
         Tree.Flag<-FALSE
         for(i in 1:length(Root.list))
         {  temp<-Tree.make(Tree.data,Root.list[i])

            if(length(temp$ROOT)!=0)
            {  empty.matrix<-matrix(0,nrow=nrow(TREE.display)+length(temp$ROOT)-1,ncol=empty.col+1)
               ind.i<-which(TREE.display[,ncol(empty.matrix)-1]==Root.list[i])
               empty.matrix[1:ind.i,1:ncol(TREE.display)]<-TREE.display[1:ind.i,]
               empty.matrix[ind.i:(ind.i+length(temp$ROOT)-1),empty.col+1]<-temp$ROOT
               if(i != length(Root.list))
                  empty.matrix[(ind.i+length(temp$ROOT)):nrow(empty.matrix),1:ncol(TREE.display)]<-TREE.display[(ind.i+1):nrow(TREE.display),]
               TREE.display<-empty.matrix
               Tree.Flag<-TRUE
            }
          
         }
         TREE<-list(Display=TREE.display,Flag=Tree.Flag)
         return(TREE)  
      }   

      Tree.data <- as.matrix(cbind(as.character(TOT.RUN$data[,1]),as.character(TOT.RUN$data[,3])))
      colnames(Tree.data)<-c("ID","parents")      
      Tree.data<-matrix(c(Tree.data[!duplicated(Tree.data[,"ID"],fromLast=T),]),ncol=2)
       
      temp<-Tree.make(Tree.data,"ROOT")
      Root.list<-temp$ROOT    
      TREE.display<-matrix(Root.list,ncol=1)      
      flag<-T
      while(flag)
      {  TREE<-expandTree(Tree.data,Root.list,TREE.display)
         TREE.display<-TREE$Display
         Root.list<-TREE.display[,ncol(TREE.display)]
         flag<-TREE$Flag
      }   
      A<-TREE.display
      colnames(A)<-c("ROOT",rep(" ",ncol(A)-1))
      MT<-gwindow("Model tree")
      gtable(A,cont=MT)
   }
   
   #### Explore output #########################################################
   ###### Select output data ###################################################
   outputselect<-function(h,...)
   {  selRUNnum<-function(h,...)
      {  selectIDs<-function(h,...)
         {  D.data<-Orig.Data
            DD1<-svalue(VarList.cond1)
            DD1.F<-svalue(From.id.cond1)
            DD1.T<-svalue(To.id.cond1)
            DD2<-svalue(VarList.cond2)
            DD2.F<-svalue(From.id.cond2)
            DD2.T<-svalue(To.id.cond2)
            temp.i<-which(toupper(colnames(Orig.Data))=="MDV")
            if(DD1!="NONE     ")
            {  if(DD2!="NONE     ")
               {  if(DD1.F==" " | DD2.F==" " | DD1.T ==" " |DD2.T==" ")
                  {  gmessage("Enter select condition",icon="error")
                  } else
                  {  
                     if(length(temp.i)!=0)
                     {  sel.id<-which(D.data[,DD1]>=as.numeric(svalue(DD1.F)) & D.data[,DD1]<=as.numeric(svalue(DD1.T)) &
                              D.data[,DD2]>=as.numeric(svalue(DD2.F)) & D.data[,DD2]<=as.numeric(svalue(DD2.T))&D.data[,temp.i]==0 )
                     } else
                     {  sel.id<-which(D.data[,DD1]>=as.numeric(svalue(DD1.F)) & D.data[,DD1]<=as.numeric(svalue(DD1.T)) &
                              D.data[,DD2]>=as.numeric(svalue(DD2.F)) & D.data[,DD2]<=as.numeric(svalue(DD2.T)) )
                     }                              
                  }
               } else
               {  if(length(temp.i)!=0)
                  { sel.id<-which(D.data[,DD1]>=as.numeric(svalue(DD1.F)) & D.data[,DD1]<=as.numeric(svalue(DD1.T))&D.data[,temp.i]==0)
                  } else
                  { sel.id<-which(D.data[,DD1]>=as.numeric(svalue(DD1.F)) & D.data[,DD1]<=as.numeric(svalue(DD1.T)))
                  }                     
               } 
            } else
            { if(length(temp.i)!=0)
              {  sel.id<-(1:nrow(D.data))[D.data[,temp.i]==0]
              } else
              {  sel.id<-(1:nrow(D.data))
              }               
            }
            SEL.ID<<-sel.id
            TEMP<-read.table(paste(file.id,".noh",sep=""),skip=1,header=T)   
            if(length(sel.id)!=nrow(D.data))
            {  OUTPUT.file<<- paste(gfile(text="Save selected data as csv",
                  type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep="")        
               write.csv(TEMP[SEL.ID,],OUTPUT.file,row.names=F,quote=F)     
            } else
            {  OUTPUT.file<<-paste(svalue(id.sel),".noh",sep="")
            }  
            dispose(win)
         }
         file.id<-svalue(id.sel)
         FILE.ID<<-file.id
         runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==file.id),"path"]
         setwd(runnum.path)
         Orig.Data<-read.csv(paste(file.id,".csv",sep=""))
         Var.Name<-colnames(Orig.Data)
         temp.id<-which(toupper(Var.Name)=="MDV")
         if(length(temp.id)!=0)
           Var.Name<-Var.Name[-temp.id]
         From.label.1<-glabel("From")
         To.label.1<-glabel("To")
         VarList.cond1<-gdroplist(c("NONE     ",Var.Name))
         From.id.cond1<-gedit(" ",width=10)
         To.id.cond1<-gedit(" ",width=10)
         From.label.2<-glabel("From")
         To.label.2<-glabel("To")
         VarList.cond2<-gdroplist(c("NONE     ",Var.Name))
         From.id.cond2<-gedit(" ",width=10)
         To.id.cond2<-gedit(" ",width=10)
         VarList.X<-gdroplist(Var.Name)
         VarType.g1<-gradio(c("Categorical","Continuous"))
         VarList.Y<-gdroplist(Var.Name)
         VarType.g2<-gradio(c("Categorical","Continuous"))
         Button1<-gbutton("OK",handler=selectIDs)

         tmp<-gframe(" Select level 1 ",container=group)
         add(tmp,VarList.cond1)
         add(tmp,From.label.1);add(tmp,From.id.cond1)
         add(tmp,To.label.1);add(tmp,To.id.cond1)
         tmp<-gframe(" Select level 2",container=group)
         add(tmp,VarList.cond2)
         add(tmp,From.label.2);add(tmp,From.id.cond2)
         add(tmp,To.label.2);add(tmp,To.id.cond2)
         tmp<-gframe("  ",container=group)
         add(tmp,Button1)
      }
      Button<-gbutton("OK",handler=selRUNnum)
      runid.list<-unique(TOT.RUN$data[,"ID"])
      id.sel<-gdroplist(runid.list)

      win<-gwindow("Select output data")
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe("Run number",container=group)
      add(tmp,id.sel)
      add(tmp,Button)
   }
      
   ###### View run summary #####################################################
   showRES<-function(h,...)
   {   runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
       tot.res<-read.csv(paste(runnum.path,"\\",FILE.ID,".sum",sep=""))[,-1]
       temp<-matrix(as.character(as.matrix(tot.res)),ncol=ncol(tot.res))
       temp[is.na(temp)]<-""
       tot.res<-temp
       colnames(tot.res)<-c("Parameters","Estimates","SE","%RSE","Lower","Upper","%Shrinkage","%CV")
       gtable(tot.res, cont=gwindow(paste(FILE.ID,".sum",sep="")),do.subset=TRUE,width=150)
    }
   
   ###### Measures of performance 1 ##############################################
   Measure.Performance1<-function(h,...)
   {  MeasureP<-function(Ex.data)
      {  ID.list.t<-sort(as.numeric(unique(Ex.data$ID)))
         MP.Indiv<-NULL
        for(i in ID.list.t)
         {  sel.id<-which(Ex.data$ID==i)
            Ni<-length(sel.id)
            temp.X<-Ex.data$X[sel.id]
            temp.Y<-Ex.data$Y[sel.id]
            temp.Y.hat<-Ex.data$Y.hat[sel.id]
            R2i<-sum((temp.Y-temp.Y.hat)^2)/sum((temp.Y-mean(temp.Y))^2)
            PEi<-(temp.Y-temp.Y.hat)/temp.Y.hat*100
            MDAWR.i<-median(abs(PEi),na.rm=T)
            MDWR.i<-median(PEi,na.rm=T)
            MAWR.i<-mean(abs(PEi),na.rm=T)
            MWR.i<-mean(PEi,na.rm=T)
            RMSWR.i<-sqrt(mean((PEi)^2,na.rm=T))
            
            MP.Indiv<-rbind(MP.Indiv,c(i,MDWR.i,MDAWR.i,MAWR.i,MWR.i))            
         }
         PE<-(Ex.data$Y-Ex.data$Y.hat)/Ex.data$Y.hat*100
         MDAWR<-median(abs(PE),na.rm=T)
         MDWR<-median(PE,na.rm=T)
         MAWR<-mean(abs(PE),na.rm=T)
         MWR<-mean(PE,na.rm=T)      
         RMSWR<-sqrt(mean((PE)^2,na.rm=T))          
         MP.pop<-matrix(c(MDWR,MDAWR,MAWR,MWR),nrow=1)
         colnames(MP.Indiv)<-c("ID","MDWR","MDAWR","MAWR","MWR")
         colnames(MP.pop)<-c("MDWR","MDAWR","MAWR","MWR")
         MP.pop.temp<-round(MP.pop,5)

         g1<-gwindow("Performance errors - Population")
         gtable(MP.pop.temp,cont=g1)

         g2<-gwindow("Performance errors - Individual")
         gtable(MP.Indiv,cont=g2)

         MP.Indiv<<-MP.Indiv
         MP.pop<<-MP.pop       
      }

      save.MP.pop<-function(h,...)
      {  MP.pop.name<-gfile(text="Save performance errors - population as csv",type="save")
         write.csv(MP.pop,paste(MP.pop.name,".csv",sep=""))
      }

      save.MP.indiv<-function(h,...)
      {  MP.indiv.name<-gfile(text="Save performance errors - individual as csv",type="save")
         write.csv(MP.Indiv,paste(MP.indiv.name,".csv",sep=""))
      }

      calc.MP<-function(h,...)
      {  ID.list<-NULL
         for(i in 1:3)
            ID.list<-c(ID.list,svalue(MP.list[[i]]))
         EX.data.temp<-EX.data.T[,ID.list]
         colnames(EX.data.temp)<-c("ID","Y","Y.hat")
         EX.data<-data.frame(EX.data.temp)
         MeasureP(EX.data)

         tmp<-gframe("Save performance errors - Population",cont=Bgroup1)
         gbutton("OK",handler=save.MP.pop,cont=tmp)
         tmp<-gframe("Save performance errors - Individual",cont=Bgroup1)
         gbutton("OK",handler=save.MP.indiv,cont=tmp)
      }
      
      select.MP<-function(h,...)
      {  noh.name<-gfile(text="data file",type="open")
         svalue(file.MP)<-noh.name
         EX.data.T<<-read.csv(noh.name,na.string=".")
         temp.list<-colnames(EX.data.T)  
         MPparam.input<-c("ID      ","Observation(Y)   ","Prediction(Y_hat)    ")
         MP.g<-list()
         MP.list<-list()
         for(i in 1:3)
         {  MP.g[[i]]<-ggroup(cont=Bgroup1)
            temp<- gdroplist(temp.list,cont=MP.g[[i]])
            glabel(MPparam.input[i],cont=MP.g[[i]])
            MP.list[[i]]<-temp
         }
         MP.list<<-MP.list
         tmp<-gframe("Calculate performance error",cont=Bgroup1)
         gbutton("OK",handler=calc.MP,cont=tmp)
         tmp<-gframe("Formulae",cont=Bgroup1,horizontal=FALSE)        
         glabel("Weighted residual(WR) = (observation-prediction)/prediction x 100(%)\n",cont=tmp)                
         glabel("Median weighted residual(MDWR) = median(WR)\n",cont=tmp)        
         glabel("Median absolute weighted residual(MDAWR) = median(|WR|)\n",cont=tmp)        
         glabel("Mean weighted residual(MWR) = mean(WR)\n",cont=tmp)        
         glabel("Mean absolute weighted residual(MAWR) = mean(|WR|)\n",cont=tmp)              
      }
      MP.win<<-gwindow("Performance errors")
      BBgroup<-ggroup(cont=MP.win,horizontal=TRUE)
      Bgroup1<-ggroup(cont=BBgroup, horizontal=FALSE)
      tmp<-gframe("Open csv file with ID, observation and prediction",cont=Bgroup1)
      file.MP<<-gedit(" ",cont=tmp)
      gbutton("Open",handler=select.MP,cont=tmp)
   }
   
   ###### Measures of performance 2 ##############################################
   Measure.Performance2<-function(h,...)
   {  Assoc.m<-function(T)
      {  C<-0 ; D<-0; keep<-NULL
         n1<-nrow(T); n2<-ncol(T)
         for(i in 1:(n1-1))
         for(j in 1:(n2-1))
         {  C<-C+T[i,j]*sum(T[(i+1):n1,(j+1):n2])
            D<-D+sum(T[i,(j+1):n2])*sum(T[(i+1):n1,j])
         }

         Z<-sum(T*(T-1)/2)
         P<-sum(T)*(sum(T)-1)/2
         X0<-0
         for(i in 1:n1)
         for(j in 1:(n2-1))
         {  X0<-X0+T[i,j]*sum(T[i,(j+1):n2])
         }

         Y0<-0
         for(j in 1:n2)
         for(i in 1:(n1-1))
         {  Y0<-Y0+T[i,j]*sum(T[(i+1):n1,j])
         }

         assoc.m<-data.frame(C,D,X0,Y0,Z,P)
         return(assoc.m)
      }

      Pred.Prob<-function(X,Y,XY.positive=TRUE)
      {  X.c<-round(X)
         Y.c<-round(Y)
         X.Y<-table(X.c,Y.c)
         if(!XY.positive)
            X.Y<-X.Y[,ncol(X.Y):1]
         temp.A<-Assoc.m(X.Y)
         Sommer.d<-(temp.A$C-temp.A$D)/(temp.A$C+temp.A$D+temp.A$Y0)
         Pk<- (1+Sommer.d)/2
         return(Pk)
      }

      MeasureP<-function(Ex.data)
      {  ID.list.t<-sort(as.numeric(unique(Ex.data$ID)))
         MP.Indiv<-NULL
         XY.positive<-ifelse(svalue(XY)=="Same",TRUE,FALSE)
         for(i in ID.list.t)
         {  sel.id<-which(Ex.data$ID==i)
            Ni<-length(sel.id)
            temp.X<-Ex.data$X[sel.id]
            temp.Y<-Ex.data$Y[sel.id]
            temp.Y.hat<-Ex.data$Y.hat[sel.id]
            Pki<-Pred.Prob(temp.X,temp.Y.hat,XY.positive)
            
            MP.Indiv<-rbind(MP.Indiv,c(i,Pki))            
         }
         Pk<-Pred.Prob(Ex.data$X,Ex.data$Y.hat,XY.positive)
        MP.pop<-matrix(Pk,nrow=1)
         colnames(MP.Indiv)<-c("ID","PredictionProbability")
         colnames(MP.pop)<-c("PredictionProbability")
         MP.pop.temp<-round(MP.pop,5)

         g1<-gwindow("Prediciton probability - Population")
         gtable(MP.pop.temp,cont=g1)

         g2<-gwindow("Prediction probability - Individual")
         gtable(MP.Indiv,cont=g2)

         MP.Indiv<<-MP.Indiv
         MP.pop<<-MP.pop       
      }

      save.MP.pop<-function(h,...)
      {  MP.pop.name<-gfile(text="Save prediction probability - population as csv",type="save")
         write.csv(MP.pop,paste(MP.pop.name,".csv",sep=""))
      }

      save.MP.indiv<-function(h,...)
      {  MP.indiv.name<-gfile(text="Save prediction probability - individual as csv",type="save")
         write.csv(MP.Indiv,paste(MP.indiv.name,".csv",sep=""))
      }

      calc.MP<-function(h,...)
      {  ID.list<-NULL
         for(i in 1:4)
            ID.list<-c(ID.list,svalue(MP.list[[i]]))
         EX.data.temp<-EX.data.T[,ID.list]
         colnames(EX.data.temp)<-c("ID","X","Y","Y.hat")
         EX.data<-data.frame(EX.data.temp)
         MeasureP(EX.data)

         tmp<-gframe("Save prediction probability - Population",cont=Bgroup1)
         gbutton("OK",handler=save.MP.pop,cont=tmp)
         tmp<-gframe("Save prediction probability - Individual",cont=Bgroup1)
         gbutton("OK",handler=save.MP.indiv,cont=tmp)
      }
      select.MP<-function(h,...)
      {  noh.name<-gfile(text="data file",type="open")
         svalue(file.MP)<-noh.name
         EX.data.T<<-read.csv(noh.name)
         temp.list<-colnames(EX.data.T)  
         MPparam.input<-c("ID      ","Predictor(X)    ","Observation(Y)   ","Prediction(Y_hat)    ")
         MP.g<-list()
         MP.list<-list()
         for(i in 1:4)
         {  MP.g[[i]]<-ggroup(cont=Bgroup1)
            temp<- gdroplist(temp.list,cont=MP.g[[i]])
            glabel(MPparam.input[i],cont=MP.g[[i]])
            MP.list[[i]]<-temp
         }
         MP.list<<-MP.list
         tmp<-gframe("Effect of predictor on prediction",cont=Bgroup1)
         XY<<-gradio(c("Stimulation","Inhibition"),selected=1,cont=tmp,horizontal=TRUE)
         gbutton("OK",handler=calc.MP,cont=tmp)
         tmp<-gframe("Formulae",cont=Bgroup1,horizontal=FALSE)        
         glabel("Prediction probability(Pk) = (Somers' d+1)/2\n",cont=tmp) 
      }
      MP.win<<-gwindow("Prediction probability")
      BBgroup<-ggroup(cont=MP.win,horizontal=TRUE)
      Bgroup1<-ggroup(cont=BBgroup, horizontal=FALSE)
      tmp<-gframe("Open csv file with ID, predictor, observation and prediction",cont=Bgroup1)
      file.MP<<-gedit(" ",cont=tmp)
      gbutton("Open",handler=select.MP,cont=tmp)
   }
        
   ###### Plot #################################################################
   ######## XY plot ############################################################
   postXY.plot<-function(h,...)
   {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
      datafile.name1<-paste(FILE.ID,".csv",sep="")
      D.data1<-read.csv(datafile.name1,na.string=".")
      colnames(D.data1)<-paste(colnames(D.data1),"-EDA  ",sep="")
      datafile.name2<-paste(FILE.ID,".noh",sep="")
      D.data2<-read.table(datafile.name2,skip=1,header=T)
      D.data<-cbind(D.data2,D.data1)[SEL.ID,]
      
      Var.Name.post<-colnames(D.data)
      updatePlot<-function(h,...)
      {  condX.V <-svalue(VarList.X,index=T)
         condY.V<-svalue(VarList.Y,index=T)
         print(condX.V)
         print(condY.V)
         select.data<-D.data
         if(!is.na(condX.V)& !is.na(condY.V))
         {  X<-select.data[,condX.V]
            Y<-select.data[,condY.V]
            dev.set(which=Dev.XYplot)
            plot(X,Y,xlab=Var.Name.post[condX.V],ylab=Var.Name.post[condY.V])
            lines(lowess(X,Y),col=2,lwd=2)                
         }
      }
      
      saveData<-function(h,...)
      {  condX.V <-svalue(VarList.X,index=T)
         condY.V<-svalue(VarList.Y,index=T)
         select.data<-D.data[,c(ID.id,condX.V,condY.V)]
         if(is.null(ID.id))
         {  colnames(select.data)<-Var.Name.post[c(condX.V,condY.V)]
         } else
         {  colnames(select.data)<-Var.Name.post[ c(ID,id,condX.V,condY.V)]                                                                                                                                            
         }
         write.csv(select.data,paste(gfile(text="Save as csv",
             type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)          
      }
      
      VarList.X<-gdroplist(Var.Name.post)
      VarList.Y<-gdroplist(Var.Name.post)
      Button1<-gbutton("OK",handler=updatePlot)
      Button2<-gbutton("Save",handler=saveData)   
      win<-gwindow("XY plot")
      BigGroup<-ggroup(cont=win)
      group<-ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe(" X variable",container=group)
      add(tmp,VarList.X)
      tmp<-gframe(" Y variable",container=group)
      add(tmp,VarList.Y)

      tmp<-gframe("Plot",container=group)
      add(tmp,Button1,expand=TRUE)
      tmp<-gframe("Save",container=group)
      add(tmp,Button2,expand=TRUE)
      add(BigGroup,ggraphics())
      Dev.XYplot<-dev.cur()
   }  
    
   ######## PRED vs DV #########################################################
   DVvsPRED.plot<-function(h,...)
   {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
      outfile.name<-paste(FILE.ID,".noh",sep="")
      Output.Table<-read.table(outfile.name,skip=1,header=T)
      select.data<-Output.Table[SEL.ID,]
      ID<-as.character(sort(unique(select.data$ID)))
      ID<-matrix(ID,nrow=length(ID))
      colnames(ID)<-c("ID")

      updatePlot<-function(h,...)
      {  condX.V <-"DV"
         condY.V<-svalue(PRED.sel)
         X<-select.data[,3]
         del.id<-which(X==0)
         if(length(del.id)!=0)
            X<-X[-del.id]
         if(is.null(X))
            X<-rep(NA,nrow(select.data)-length(del.id))
         if(length(del.id)!=0)
         {  Y<-select.data[-del.id,condY.V]
         } else
         {  Y<-select.data[,condY.V]
         }         
           
         if(is.null(Y))
            Y<-rep(NA,nrow(select.data)-length(del.id))
         plot(X,Y,xlab=condX.V,ylab=condY.V)
         abline(a=0,b=1,col=2)        
      }
      
      updatePlot1<-function(h,...)
      {  condX.V <-"DV"
         condY.V<-svalue(PRED.sel)
         X<-select.data[,3]
         del.id<-which(X==0)
         if(length(del.id)!=0)
            X<-X[-del.id]
         if(is.null(X))
            X<-rep(NA,nrow(select.data)-length(del.id))
         if(length(del.id)!=0)
         {  Y<-select.data[-del.id,condY.V]
         } else
         {  Y<-select.data[,condY.V]
         }         
         if(is.null(Y))
            Y<-rep(NA,nrow(select.data)-length(del.id))
         plot(X,Y,xlab=condX.V,ylab=condY.V)
         abline(a=0,b=1,col=2)    
         sel.id<-which(select.data$ID==svalue(IDlist))
         X<-select.data[sel.id,3]
         del.id<-which(X==0)
         if(length(del.id)!=0)
            X<-X[-del.id]
         if(is.null(X))
            X<-rep(NA,nrow(select.data)-length(del.id))
         Y<-select.data[sel.id,condY.V]
         if(length(del.id)!=0)
            Y<-Y[-del.id]
         points(X,Y,col=2,pch=16)        
      }
      
      Button<-gbutton("OK",handler=updatePlot)
      win<-gwindow(paste("DV vs Predictions plot (",OUTPUT.file,")",sep=""))
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
      tmp=gframe(" Y variable",container=group)
      if(read.table("NM.version")=="NM6")
      {  PRED.list<-PRED.list.6
      } else
      {  PRED.list<-PRED.list.7
      }
      PRED.sel<-gdroplist(PRED.list)
      add(tmp,PRED.sel)

      tmp=gframe("Plot",container=group)
      add(tmp,Button,expand=TRUE)
      tmp<-gframe("ID",container=group)
      IDlist <- gtable(ID,multiple=T,handler=updatePlot1) 
      size(IDlist)<-c(100,200)
      add(tmp,IDlist)
      add(BigGroup,ggraphics())
   }
   
   ######## RES vs DV ##########################################################
   DVvsRES.plot<-function(h,...)
   {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
      outfile.name<-paste(runnum.path,"\\",FILE.ID,".noh",sep="")
      Output.Table<-read.table(outfile.name,skip=1,header=T)
      select.data<-Output.Table[SEL.ID,]
      ID<-as.character(sort(unique(select.data$ID)))
      ID<-matrix(ID,nrow=length(ID))
      colnames(ID)<-c("ID")
    
      updatePlot<-function(h,...)
      {  condX.V <-"DV"
         condY.V<-svalue(RES.sel)
         X<-select.data[,3]
         del.id<-which(X==0)
         if(length(del.id)!=0)
            X<-X[-del.id]
         if(is.null(X))
            X<-rep(NA,nrow(select.data)-length(del.id))
         if(length(del.id)!=0)
         {  Y<-select.data[-del.id,condY.V]
         } else
         {  Y<-select.data[,condY.V]
         }         
         if(is.null(Y))
            Y<-rep(NA,nrow(select.data)-length(del.id))
         
         par(mfrow=c(1,1))
         plot(X,Y,xlab=condX.V,ylab=condY.V)
         abline(h=0,col=2)
         WRES.list<-c("WRES","IWRES","NWRES","WRESI","CWRES","CRESI","CWRESI","EWRES","NPDE")
         if(sum(WRES.list==condY.V)!=0)
         {  abline(h=2,col=2,lty=2)
            abline(h=-2,col=2,lty=2)
         }        
      }
      
      updatePlot1<-function(h,...)
      {  condX.V <-"DV"
         condY.V<-svalue(RES.sel)
         X<-select.data[,3]
         del.id<-which(X==0)
         if(length(del.id)!=0)
            X<-X[-del.id]
         if(is.null(X))
            X<-rep(NA,nrow(select.data)-length(del.id))
         if(length(del.id)!=0)
         {  Y<-select.data[-del.id,condY.V]
         } else
         {  Y<-select.data[,condY.V]
         }
         if(is.null(Y))
            Y<-rep(NA,nrow(select.data)-length(del.id))
         
         par(mfrow=c(1,1))
         plot(X,Y,xlab=condX.V,ylab=condY.V)
         abline(h=0,col=2)
         WRES.list<-c("WRES","IWRES","NWRES","WRESI","CWRES","CRESI","CWRESI","EWRES","NPDE")
         if(sum(WRES.list==condY.V)!=0)
         {  abline(h=2,col=2,lty=2)
            abline(h=-2,col=2,lty=2)
         }   
          sel.id<-which(select.data$ID==svalue(IDlist))
         X<-select.data[sel.id,3]
         del.id<-which(X==0)
         if(length(del.id)!=0)
             X<-X[-del.id]
         if(is.null(X))
           X<-rep(NA,nrow(select.data)-length(del.id))
         Y<-select.data[sel.id,condY.V]
         if(length(del.id)!=0)
            Y<-Y[-del.id]
         points(X,Y,col=2,pch=16)        
              
      }

      Button<-gbutton("OK",handler=updatePlot)
      win<-gwindow(paste("DV vs Residuals plot (",OUTPUT.file,")",sep=""))
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)

      tmp=gframe(" Y variable",container=group)
      if(read.table("NM.version")=="NM6")
      {  RES.list<-RES.list.6
      } else
      {  RES.list<-RES.list.7
      }

      RES.sel<-gdroplist(RES.list)
      add(tmp,RES.sel)
      tmp=gframe("Plot",container=group)
      add(tmp,Button,expand=TRUE)
      tmp<-gframe("ID",container=group)
      IDlist <- gtable(ID,multiple=T,handler=updatePlot1) 
      size(IDlist)<-c(100,200)
      add(tmp,IDlist)
      add(BigGroup,ggraphics())
   }   
   
   ######## RES vs TIME ########################################################
   TIMEvsRES.plot<-function(h,...)
   {  updatePlot<-function(h,...)
      {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
         outfile.name<-paste(runnum.path,"\\",FILE.ID,".NOH",sep="")
         Output.Table<-read.table(outfile.name,skip=1,header=T)
         condX.V <-"TIME"
         condY.V<-svalue(RES.sel)
   
         select.data<-Output.Table[SEL.ID,]
         del.id<-which(select.data[,3]==0)
         if(length(del.id)!=0)
         {   X<-select.data[-del.id,condX.V]
         } else
         {   X<-select.data[,condX.V]
         }          
         if(is.null(X))
           X<-rep(NA,nrow(select.data)-length(del.id))
         if(length(del.id)!=0)
         {  Y<-select.data[-del.id,condY.V]
         } else
         {  Y<-select.data[,condY.V]
         }          
         if(is.null(Y))
           Y<-rep(NA,nrow(select.data)-length(del.id))

         plot(X,Y,xlab=condX.V,ylab=condY.V,type='n')
         id<-as.numeric(names(table(Output.Table$ID)))
         for(i in 1:length(id))
         {  if(length(del.id)!=0)
            {  X.data<-X[which(select.data[-del.id,"ID"]==id[i])]
               Y.data<-Y[which(select.data[-del.id,"ID"]==id[i])]
            } else
            {  X.data<-X[which(select.data[,"ID"]==id[i])]
               Y.data<-Y[which(select.data[,"ID"]==id[i])]
            }               
            lines(X.data,Y.data,lty=2)
         }
         abline(h=0,col=2)
         
         WRES.list<-c("WRES","IWRES","NWRES","WRESI","CWRES","CRESI","CWRESI","EWRES","NPDE")

         if(sum(WRES.list==condY.V)!=0)
         {  abline(h=2,col=2,lty=2)
            abline(h=-2,col=2,lty=2)
         }
       }

      Button<-gbutton("OK",handler=updatePlot)

      win<-gwindow(paste("TIME vs Residuals plot (",OUTPUT.file,")",sep=""))
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)

      tmp=gframe(" Y variable",container=group)
      if(read.table("NM.version")=="NM6")
      {  RES.list<-RES.list.6
      } else
      {  RES.list<-RES.list.7
      }
      RES.sel<-gdroplist(RES.list)
      add(tmp,RES.sel)

      tmp=gframe("Plot",container=group)
      add(tmp,Button,expand=TRUE)
      add(BigGroup,ggraphics())
   }
   
   ######## Prediction and DV vs TIME ##########################################
   TIMEvsDVandPRED.plot<-function(h,...)
   {  updatePlot<-function(h,...) 
      {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
         outfile.name<-paste(runnum.path,"\\",FILE.ID,".noh",sep="")
         Output.Table<-read.table(outfile.name,skip=1,header=T)
         D.data<-read.csv(paste(FILE.ID,".csv",sep=""))
         condX.V <-"TIME"
         select.data<-Output.Table[SEL.ID,]
         condY.V<-svalue(PRED.sel)
         del.id<-which(select.data[,3]==0)
         if(length(del.id)!=0)
         {  X<-select.data[-del.id,condX.V]
            Y<-select.data[-del.id,3]
            ID<-select.data[-del.id,1]
            Y1<-select.data[-del.id,condY.V]  
         } else
         {  X<-select.data[,condX.V]
            Y<-select.data[,3]
            ID<-select.data[,1]
            Y1<-select.data[,condY.V]  
         }   
         plot(X,Y,type='n',xlab=condX.V,ylab=condY.V,
            ylim=range(c(Y,Y1)))#,xlim=range(Output.Table[,condX.V]),ylim=range(Output.Table[,condY.V]))
         id.list<-unique(Output.Table$ID)
         
         for(i in id.list)
         {  sel.data.id<-which(ID==i)
            lines(X[sel.data.id],Y[sel.data.id],lwd=0.1,lty=2,col="grey")
         }
         
         for(i in id.list)
         {  sel.data.id<-which(ID==i)
            lines(X[sel.data.id],Y1[sel.data.id],lwd=0.5,lty=1,col="grey50")
         }
      }
      
      Button<-gbutton("OK",handler=updatePlot)
      win<-gwindow(paste("Predictions and DV vs TIME plot (",OUTPUT.file,")",sep=""))
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
 
      tmp=gframe(" Y variable",container=group)
      if(read.table("NM.version")=="NM6")
      {  PRED.list<-PRED.list.6
      } else
      {  PRED.list<-PRED.list.7
      }

      PRED.sel<-gdroplist(PRED.list)
      add(tmp,PRED.sel)

      tmp=gframe("Plot",container=group)
      add(tmp,Button,expand=TRUE)
      add(BigGroup,ggraphics())
   }
   
   ######## Prediction and DV vs TIME by ID ####################################
   TIMEvsDVandPREDID.plot<-function(h,...)
   {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
      outfile.name<-paste(runnum.path,"\\",FILE.ID,".noh",sep="")
      Output.Table<-read.table(outfile.name,skip=1,header=T,na.string=".")
      D.data<-Output.Table[SEL.ID,]
      n<-nrow(D.data)
      p<-names(table(D.data$ID))

      updatePlot1<-function(h,...) 
      {  select.pred<-svalue(PRED.sel)
         select.id<-svalue(IDlist)
         id.list<-which(D.data$ID==select.id)
         if(length(id.list)!=0)
         {  y.lim<-range(D.data$DV[id.list],D.data$IPRED[id.list],D.data[id.list,select.pred],na.rm=T)
            plot(D.data$TIME[id.list],D.data[id.list,"DV"],ylim=y.lim,col="grey10",cex=1.2,xlab="TIME",ylab="",main=paste("ID",select.id,sep=" "))
            lines(D.data$TIME[id.list],D.data[id.list,"DV"],lwd=1)
            lines(D.data$TIME[id.list],D.data[id.list,select.pred],col=4,lwd=1)
            lines(D.data$TIME[id.list],D.data[id.list,"IPRED"],col=2,lwd=1)
         }
      }
      
      ID<-as.character(sort(unique(D.data$ID)))
      ID<-matrix(ID,nrow=length(ID))
      colnames(ID)<-c("ID")

      window<-gwindow("Predictions and DV vs TIME by ID plot")
      Biggroup<-ggroup(cont=window)
      group=ggroup(horizontal=FALSE,cont=Biggroup)
      tmp<-gframe(" Y variable",container=group)
      if(read.table("NM.version")=="NM6")
      {  PRED.list<-PRED.list.6
      } else
      {  PRED.list<-PRED.list.7
      }
      PRED.list<-PRED.list[PRED.list!="IPRED"]
      PRED.sel<-gdroplist(PRED.list)
      add(tmp,PRED.sel)
      tmp<-gframe("ID",container=group)
      IDlist <- gtable(ID,multiple=T,handler=updatePlot1) 
      size(IDlist)<-c(100,200)
      add(tmp,IDlist)
      tmp<-gframe("DV   : black",container=group)
      tmp<-gframe("IPRED : red",container=group)
      tmp<-gframe("Selected prediction : blue",container=group)   
      add(Biggroup,tmp)
      add(Biggroup, ggraphics())       
   }   
   
   ######## Covariate vs parameter #############################################
   EBEvsCOV.plot<-function(h,...)
   {  updatePlot<-function(h,...)
      {  runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
         outfile.name<-paste(runnum.path,"\\",FILE.ID,".noh",sep="")
         Output.Table<-read.table(outfile.name,skip=1,header=T)
         condX.V <-svalue(EBE.sel)
         condY.V<-svalue(COV.sel)
         X<-TOT.data[,condX.V]
         Y<-TOT.data[,condY.V]
         plot(Y,X,xlab=condY.V,ylab=condX.V)
         lines(lowess(Y,X),col=2,lwd=2)        
      }
      
      savePlot<-function(h,...)
      {  condX.V <-svalue(EBE.sel)
         condY.V<-svalue(COV.sel)
         write.csv(T.data[,c("X.ID",condX.V,condY.V)],paste(gfile(text="Save selected data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }     
      
      runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==FILE.ID),"path"]
      runnum.path<-runnum.path[length(runnum.path)]
      outfile.name<-paste(runnum.path,"\\",FILE.ID,".noh",sep="")
      datafile.name<-paste(runnum.path,"\\",FILE.ID,".csv",sep="")
      EBEfile.name<-paste(runnum.path,"\\",FILE.ID,".eta",sep="")
      PARfile.name<-paste(runnum.path,"\\",FILE.ID,".par",sep="")
    
      Orig.data<-read.csv(datafile.name,na.string=".")
      EBE.data<-read.table(EBEfile.name,skip=1,header=T)
      PAR.data<-read.table(PARfile.name,skip=1,header=T)   
      if(nrow(EBE.data)!=nrow(PAR.data))
      {  tot.data<-merge(EBE.data,PAR.data,all=T)
      } else
      {  tot.data<-cbind(EBE.data,PAR.data[,-1])
      }   
      EBE.data<-merge(EBE.data,PAR.data,all=T)
      colnames(EBE.data)[1]<-"X.ID"
      Var.Name<-colnames(Orig.data)
      EBE.name<-colnames(EBE.data);#EBE.name<-EBE.name[EBE.name!="X.ID"]
      if(nrow(EBE.data)!=nrow(Orig.data))
      {  tot.data<-merge(Orig.data,EBE.data,all=T)
      } else
      {  tot.data<-cbind(Orig.data,EBE.data[,-1])
      }      
      
      TOT.data<-tot.data[SEL.ID,]
      T.data<<-TOT.data
      Button2<-gbutton("OK",handler=updatePlot)
      Button3<-gbutton("OK",handler=savePlot)

      win<-gwindow("Parameter vs covariate")
      BigGroup<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
      tmp=gframe(" Parameter",container=group)
      EBE.sel<-gdroplist(EBE.name)
      add(tmp,EBE.sel)
 
      tmp=gframe(" Covariate",container=group)
      COV.sel<-gdroplist(Var.Name)
      add(tmp,COV.sel)

      tmp=gframe("Plot",container=group)
      add(tmp,Button2,expand=TRUE)
 
      tmp=gframe("Save",container=group)
      add(tmp,Button3,expand=TRUE)     
      add(BigGroup,ggraphics())
   }
       
   #### Convert PK parameters ##################################################
   PKparam.converter<-function(h,...)
   {  save.convert.PK<-function(h,...)
      {  PK.convert.name<-gfile(text="Save converted PK parameters as csv",type="save")
         write.csv(convert.T,paste(PK.convert.name,".csv",sep=""),row.names=FALSE)
      }

      convert.V<-function(h,...)
      {  Method<-svalue(Method.t); Comp<-svalue(Comp.t)
         PK.name.list<-NULL
         for(i in 1:(length(PKparam.input)+1))
            PK.name.list<-c(PK.name.list, svalue(PK.list[[i]]))
         input.PK.T<-PK.data[,PK.name.list]
         colnames(input.PK.T)<-c("ID",PKparam.input)
         convert.T<-NULL
         for(i in 1:nrow(input.PK.T))
         {  input.PK<-data.frame(input.PK.T[i,-1])
            names(input.PK)<-PKparam.input
            if(Method=="Volumes and clearances" & Comp=="1-comp")
            {  convert.P<-convert.comp1.m1(input.PK)
            } else  if(Method=="Volumes and clearances" & Comp=="2-comp")
            {  convert.P<-convert.comp2.m1(input.PK)
            } else  if(Method=="Volumes and clearances" & Comp=="3-comp")
            {  convert.P<-convert.comp3.m1(input.PK)
            } else  if(Method=="V1,Rate Constants" & Comp=="1-comp")
            {  convert.P<-convert.comp1.m2(input.PK)
            } else  if(Method=="V1,Rate Constants" & Comp=="2-comp")
            {  convert.P<-convert.comp2.m2(input.PK)
            } else  if(Method=="V1,Rate Constants" & Comp=="3-comp")
            {  convert.P<-convert.comp3.m2(input.PK)
            } else  if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="1-comp")
            {  convert.P<-convert.comp1.m3(input.PK)
            } else  if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="2-comp")
            {  convert.P<-convert.comp2.m3(input.PK)
            } else  if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="3-comp")
            {  convert.P<-convert.comp3.m3(input.PK)
            } else  if(Method=="Coefficients and Exponents" & Comp=="1-comp")
            {  convert.P<-convert.comp1.m4(input.PK)
            } else  if(Method=="Coefficients and Exponents" & Comp=="2-comp")
            {  convert.P<-convert.comp2.m4(input.PK)
            } else  if(Method=="Coefficients and Exponents" & Comp=="3-comp")
            {  convert.P<-convert.comp3.m4(input.PK)
            } else  if(Method=="V1,exponents, K21,K31" & Comp=="1-comp")
            {  convert.P<-convert.comp1.m5(input.PK)
            } else  if(Method=="V1,exponents, K21,K31" & Comp=="2-comp")
            {  convert.P<-convert.comp2.m5(input.PK)
            } else  if(Method=="V1,exponents, K21,K31" & Comp=="3-comp")
            {  convert.P<-convert.comp3.m5(input.PK)
            }  
            ID<-input.PK.T$ID[i]
            convert.T<-rbind(convert.T,data.frame(ID,t(convert.P)))
         }
         convert.P.print<-matrix(as.character(round(unlist(convert.T),6)),ncol=ncol(convert.T))
         convert.P.print[which(is.na(convert.P.print))]<-"  "
         colnames(convert.P.print)<-colnames(convert.T)
         gtable(convert.P.print,cont=gwindow())
         convert.T<<-convert.T
      }

      select.V1<-function(h,...)
      {  sum.name1<-gfile(text="Individual PK parameters(runnumber.par)",type="open")
         svalue(file.n1)<-sum.name1
         print("In V1")       
         PK.data<-read.table(sum.name1,header=T,skip=1)
         Method<-svalue(Method.t); Comp<-svalue(Comp.t)
         list.1.1<-c("V1" , "CL1");list.1.2<-c( "V1" , "V2" , "CL1" , "CL2" )
         list.1.3<-c( "V1" ,"V2" , "V3" , "CL1" , "CL2" , "CL3" )

         list.2.1<-c( "V1" ,"K10");list.2.2<-c( "V1" ,"K10" , "K12" , "K21")
         list.2.3<-c( "V1" ,"K10" , "K12" , "K13" , "K21" , "K31")

         list.3.1<-c(  "CL1" ,"t.0.5.alpha" );list.3.2<-c( "V1" , "CL1" , "t.0.5.alpha" , "t.0.5.beta" )
         list.3.3<-c( "V1" , "CL1" , "Vdss" , "t.0.5.alpha" , "t.0.5.beta" , "t.0.5.gamma")

         list.4.1<-c( "A" , "alpha" );list.4.2<-c( "A" , "B" , "alpha" , "beta" )
         list.4.3<-c( "A" , "B" , "C" , "alpha" , "beta" , "gamma")

         list.5.1<-c( "V1" , "alpha");list.5.2<-c( "V1" , "K21" , "alpha" , "beta" )
         list.5.3<-c( "V1" , "K21" , "K31" , "alpha" , "beta" , "gamma" )

         if(Method=="Volumes and clearances" & Comp=="1-comp")
         {  PKparam.input<-list.1.1
         } else if(Method=="Volumes and clearances" & Comp=="2-comp")
         {  PKparam.input<-list.1.2
         } else if(Method=="Volumes and clearances" & Comp=="3-comp")
         {  PKparam.input<-list.1.3
         } else if(Method=="V1,Rate Constants" & Comp=="1-comp")
         {  PKparam.input<-list.2.1
         } else if(Method=="V1,Rate Constants" & Comp=="2-comp")
         {  PKparam.input<-list.2.2
         } else if(Method=="V1,Rate Constants" & Comp=="3-comp")
         {  PKparam.input<-list.2.3
         } else if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="1-comp")
         {  PKparam.input<-list.3.1
         } else if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="2-comp")
         {  PKparam.input<-list.3.2
         } else if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="3-comp")
         {  PKparam.input<-list.3.3
         } else if(Method=="Coefficients and Exponents" & Comp=="1-comp")
         {  PKparam.input<-list.4.1
         } else if(Method=="Coefficients and Exponents" & Comp=="2-comp")
         {  PKparam.input<-list.4.2
         } else if(Method=="Coefficients and Exponents" & Comp=="3-comp")
         {  PKparam.input<-list.4.3
         } else if(Method=="V1,exponents, K21,K31" & Comp=="1-comp")
         {  PKparam.input<-list.5.1
         } else if(Method=="V1,exponents, K21,K31" & Comp=="2-comp")
         {  PKparam.input<-list.5.2
         } else if(Method=="V1,exponents, K21,K31" & Comp=="3-comp")
         {  PKparam.input<-list.5.3
         } 
         PK.g<-list()
         PK.list<<-list()
         PK.g[[1]]<-ggroup(cont=Bgroup1)
         glabel("ID   ",cont=PK.g[[1]])
         temp<- gdroplist(colnames(PK.data),cont=PK.g[[1]])
         PK.list[[1]]<<-temp
       
         for(i in 2:(length(PKparam.input)+1))
         {  PK.g[[i]]<-ggroup(cont=Bgroup1)
            glabel(PKparam.input[i-1],cont=PK.g[[i]])
            temp<- gdroplist(colnames(PK.data),cont=PK.g[[i]])
            PK.list[[i]]<<-temp
         }
         PKparam.input<<-PKparam.input
         tmp<-gframe("Convert and save",cont=Bgroup1)
         PK.data<<-PK.data
         gbutton("Convert",handler=convert.V,cont=tmp)
         gbutton("Save",handler=save.convert.PK,cont=tmp)
      }
    
      select.V2<-function(h,...)
      {  sum.name2<-gfile(text="Individual PK parameters(csv file)",type="open")
         svalue(file.n2)<-sum.name2
         PK.data<-read.csv(sum.name2)
         Method<-svalue(Method.t); Comp<-svalue(Comp.t)
         list.1.1<-c("V1" , "CL1");list.1.2<-c( "V1" , "V2" , "CL1" , "CL2" )
         list.1.3<-c( "V1" ,"V2" , "V3" , "CL1" , "CL2" , "CL3" )

         list.2.1<-c( "V1" ,"K10");list.2.2<-c( "V1" ,"K10" , "K12" , "K21")
         list.2.3<-c( "V1" ,"K10" , "K12" , "K13" , "K21" , "K31")

         list.3.1<-c(  "CL1" ,"t.0.5.alpha" );list.3.2<-c( "V1" , "CL1" , "t.0.5.alpha" , "t.0.5.beta" )
         list.3.3<-c( "V1" , "CL1" , "Vdss" , "t.0.5.alpha" , "t.0.5.beta" , "t.0.5.gamma")

         list.4.1<-c( "A" , "alpha" );list.4.2<-c( "A" , "B" , "alpha" , "beta" )
         list.4.3<-c( "A" , "B" , "C" , "alpha" , "beta" , "gamma")

         list.5.1<-c( "V1" , "alpha");list.5.2<-c( "V1" , "K21" , "alpha" , "beta" )
         list.5.3<-c( "V1" , "K21" , "K31" , "alpha" , "beta" , "gamma" )

         if(Method=="Volumes and clearances" & Comp=="1-comp")
         {  PKparam.input<-list.1.1
         } else if(Method=="Volumes and clearances" & Comp=="2-comp")
         {  PKparam.input<-list.1.2
         } else if(Method=="Volumes and clearances" & Comp=="3-comp")
         {  PKparam.input<-list.1.3
         } else if(Method=="V1,Rate Constants" & Comp=="1-comp")
         {  PKparam.input<-list.2.1
         } else if(Method=="V1,Rate Constants" & Comp=="2-comp")
         {  PKparam.input<-list.2.2
         } else if(Method=="V1,Rate Constants" & Comp=="3-comp")
         {  PKparam.input<-list.2.3
         } else if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="1-comp")
         {  PKparam.input<-list.3.1
         } else if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="2-comp")
         {  PKparam.input<-list.3.2
         } else if(Method=="V1,Vdss,Cl,Half-lives" & Comp=="3-comp")
         {  PKparam.input<-list.3.3
         } else if(Method=="Coefficients and Exponents" & Comp=="1-comp")
         {  PKparam.input<-list.4.1
         } else if(Method=="Coefficients and Exponents" & Comp=="2-comp")
         {  PKparam.input<-list.4.2
         } else if(Method=="Coefficients and Exponents" & Comp=="3-comp")
         {  PKparam.input<-list.4.3
         } else if(Method=="V1,exponents, K21,K31" & Comp=="1-comp")
         {  PKparam.input<-list.5.1
         } else if(Method=="V1,exponents, K21,K31" & Comp=="2-comp")
         {  PKparam.input<-list.5.2
         } else if(Method=="V1,exponents, K21,K31" & Comp=="3-comp")
         {  PKparam.input<-list.5.3
         } 
         PK.g<-list()
         PK.list<<-list()
         PK.g[[1]]<-ggroup(cont=Bgroup1)
         glabel("ID   ",cont=PK.g[[1]])
         temp<- gdroplist(colnames(PK.data),cont=PK.g[[1]])
         PK.list[[1]]<<-temp
       
         for(i in 2:(length(PKparam.input)+1))
         {  PK.g[[i]]<-ggroup(cont=Bgroup1)
            glabel(PKparam.input[i-1],cont=PK.g[[i]])
            temp<- gdroplist(colnames(PK.data),cont=PK.g[[i]])
            PK.list[[i]]<<-temp
         }
         PKparam.input<<-PKparam.input
         tmp<-gframe("Convert and save",cont=Bgroup1)
         PK.data<<-PK.data
         gbutton("Convert",handler=convert.V,cont=tmp)
         gbutton("Save",handler=save.convert.PK,cont=tmp)
      }
       
      ## method 1 
      convert.comp1.m1<-function(input)
      {  V1<-input$V1 ;  CL1<-input$CL1
         Vdss<- V1;   K10 <- CL1 /V1  ;      l1 <- CL1 /V1  ;      alpha <- l1  
         t.0.5.alpha <- log(2)/l1  ;      C1 <- 1/V1;      T.A <- C1
         F.A <- C1 *V1
         Result.P<-rbind(V1,CL1,NA,Vdss,T.A,F.A,alpha,t.0.5.alpha,K10)
         return(Result.P)
      }

      convert.comp2.m1<-function(input)
      {  V1<-input$V1 ;      V2<-input$V2
         CL1<-input$CL1;      CL2<-input$CL2
         Vdss <- V1 +V2 ;      K10 <- CL1 /V1  
         K12 <- CL2 /V1  ;      K21 <- CL2 /V2  
         a0<- K10*K21;      a1<- -(K10    + K12  + K21)
         l1  <- (-a1 + sqrt(a1*a1 - 4* a0)) / 2
         l2  <- (-a1 - sqrt(a1*a1 - 4* a0)) / 2
         alpha <- l1  
         beta <- l2  
         t.0.5.alpha <- log(2)/l1  
         t.0.5.beta  <- log(2)/l2  
         C1  <- (K21 - l1) / (l2 - l1)/V1
         C2  <- (K21 - l2) / (l1 - l2)/V1
         T.A <- C1
         T.B <- C2
         F.A <- C1 *V1
         F.B <- C2 *V1
         Result.P<-rbind(V1,V2,CL1,CL2,NA,Vdss,T.A,T.B,F.A,F.B,alpha,beta,t.0.5.alpha,t.0.5.beta,K10,K12,K21)
         return(Result.P)
      }

      convert.comp3.m1<-function(input)
      {  V1<-input$V1 ;         V2<-input$V2;         V3<-input$V3
         CL1<-input$CL1;         CL2<-input$CL2;         CL3<-input$CL3
         Vdss <- V1 +V2 +V3 
         K10 <- CL1 /V1  ;         K12 <- CL2 /V1  ;         K13 <- CL3 /V1  
         K21 <- CL2 /V2  ;         K31 <- CL3 /V3  ;         Vdss<-V1+V2+V3
         a0<-K10*K21*K31
         a1<- K10  * K31  + K21  * K31  + K21  * K13  + K10  * K21  + K31  * K12  
         a2<- K10  + K12  + K13  + K21  + K31  
         p <-a1 - (a2  * a2  / 3) 
         q <- (2 * a2  * a2  * a2  / 27) - (a1  * a2  / 3)+ a0 
         r1 <- sqrt(-(p * p * p) / 27) 
         phi <-acos((-q/ 2) / r1)/3 
         r2  <- 2 *exp(log(r1) / 3) 
         root1 <- -(cos(phi) * r2- a2  / 3) 
         root2 <- -(cos(phi+ 2*pi/3) * r2 - a2  / 3) 
         root3 <- -(cos(phi + 4*pi/3)* r2 - a2  / 3) 
         l1 <- max(c(root1,root2,root3) )
         l2 <- median(c(root1,root2,root3) )
         l3 <- min(c(root1,root2,root3) )
         alpha <- l1;       beta <- l2  ;         gamma <- l3  
         t.0.5.alpha <- log(2)/l1  
         t.0.5.beta  <- log(2)/l2  
         t.0.5.gamma  <- log(2)/l3  	
         C1  <- (K21  - l1 ) * (K31  - l1 ) / (l1  - l2 ) / (l1  - l3 )/V1  
         C2  <-(K21  - l2 ) * (K31  - l2 ) / (l2  - l1 ) / (l2  - l3 )/V1  
         C3  <-(K21  - l3 ) * (K31  - l3 ) / (l3  - l2 ) / (l3  - l1 )/V1  
         T.A<-C1;         T.B<-C2;         T.C<-C3
         F.A <-C1 *V1;         F.B <-C2 *V1;         F.C <-C3  *V1
         Result.P<-rbind(V1,V2,V3,CL1,CL2,CL3,NA,Vdss,T.A,T.B,T.C,F.A,F.B,F.C,
              alpha,beta,gamma,t.0.5.alpha,t.0.5.beta,t.0.5.gamma,K10,K12,K13,K21,K31)
         return(Result.P)
      }

      ## method 2 
      convert.comp1.m2<-function(input)
      {  V1<-input$V1 
         K10<-input$K10
         Vdss <-V1  
         CL1 <- V1  *K10   
 
         l1 <-CL1/V1
         C1 <- 1 /V1     
         T.A <- C1 	
         F.A <- C1*V1   	
         alpha <- l1   
         t.0.5.alpha <- log(2)/l1   
         Result.P<-rbind(V1,K10,NA,Vdss,CL1,T.A,F.A,alpha,t.0.5.alpha)
         return(Result.P)
      }

      convert.comp2.m2<-function(input)
      {  V1<-input$V1 
         K10<-input$K10
         K12<-input$K12
         K21<-input$K21
         V2 <- V1 *K12  /K21    
         Vdss <-V1  +V2 	
         CL1 <- V1  *K10   
         CL2 <- V1  *K12   
         a0 <- K10  * K21   
         a1 <- -(K10 + K12 + K21)
         l1 <- (-a1 + sqrt(a1*a1 - 4* a0)) / 2
         l2 <- (-a1 - sqrt(a1*a1 - 4* a0)) / 2
         C1 <- (K21   - l1  )  / (l2   - l1  ) /V1   
         C2 <-(K21   - l2  )  / (l1   - l2 ) /V1   
         T.A <- C1 
         T.B <- C2 
         F.A <- C1*V1   
         F.B <- C2*V1   
         alpha <- l1   
         beta <- l2    
         t.0.5.alpha <- log(2)/l1   
         t.0.5.beta <- log(2)/l2     
         Result.P<-rbind(V1,K10,K12,K21,NA,V2,Vdss,CL1,CL2,T.A,T.B,F.A,F.B,
                  alpha,beta,t.0.5.alpha,t.0.5.beta)   
         return(Result.P)
      }

      convert.comp3.m2<-function(input)
      {  V1<-input$V1 
         K10<-input$K10;         K12<-input$K12;         K13<-input$K13
         K21<-input$K21;         K31<-input$K31
         V2 <- V1 *K12  /K21  
         V3 <- V1  *K13  /K31   
         Vdss <-V1  +V2  +V3  	
         CL1 <- V1  *K10   
         CL2 <- V1  *K12   
         CL3 <- V1  *K13   
         a0 <- K10  * K21  *K31   
         a1 <- K10   * K31   + K21   * K31   + K21   * K13   + K10   * K21   + K31  *K12   
         a2 <- K10   + K12   + K13   + K21  +K31   	
         p <- a1   - (a2  * a2  / 3) 
         q <- (2 * a2   * a2   * a2   / 27) - (a1   * a2   / 3)+ a0 
         r1 <- sqrt(-(p* p * p) / 27) 
         phi <-acos((-q / 2)/ r1)/3 
         r2 <- 2 * exp(log(r1) / 3) 
         root1 <--(cos(phi) * r2- a2   / 3) 
         root2 <- -(cos(phi + 2*pi/3) * r2 -a2   / 3) 
         root3 <--(cos(phi + 4*pi/3) * r2 -a2   / 3) 
         l1 <-max(root1,root2,root3) 
         l2 <-median(c(root1,root2,root3)) 
         l3 <-min(root1,root2,root3) 
         C1 <- (K21   - l1  ) * (K31   - l1  ) / (l1   - l2  ) / (l1   - l3  )/V1   
         C2 <-(K21   - l2  ) * (K31   - l2  ) / (l2   - l1  ) / (l2   - l3  )/V1   
         C3 <-(K21   - l3  ) * (K31   - l3  ) / (l3   - l2  ) / (l3   - l1  )/V1   
         T.A <- C1 ;         T.B <- C2 ;         T.C <- C3 	
         F.A <- C1*V1;           F.B <- C2*V1   ;         F.C <- C3*V1   	
         alpha <- l1 ;         beta <- l2   ;         gamma <- l3   
         t.0.5.alpha <- log(2)/l1;           t.0.5.beta <- log(2)/l2   ;         t.0.5.gamma <- log(2)/l3   

         Result.P<-rbind(V1,K10,K12,K13,K21,K31,NA,V2,V3,Vdss,CL1,CL2,CL3,T.A,T.B,T.C,F.A,F.B,F.C,
               alpha,beta,gamma,t.0.5.alpha,t.0.5.beta,t.0.5.gamma)
         return(Result.P)
      }

      ## method 3 
      convert.comp1.m3<-function(input)
      {  CL1<-input$CL1
         t.0.5.alpha<-input$t.0.5.alpha
         alpha <-log(2)/t.0.5.alpha 
         l1<-alpha
         V1 <-CL1/l1	
         K10 <-CL1 /V1  
         Vdss <-V1 
         C1 <- 1/V1
         T.A <-C1 
         F.A <-C1*V1  
         Result.P<-rbind(CL1,t.0.5.alpha,NA,
                        alpha, V1,Vdss,K10,T.A,F.A)
         return(Result.P)
     }

     convert.comp2.m3<-function(input)
     {   V1<-input$V1 
         CL1<-input$CL1
         t.0.5.alpha<-input$t.0.5.alpha
         t.0.5.beta<-input$t.0.5.beta
         alpha <-log(2)/t.0.5.alpha 
         beta <-log(2)/t.0.5.beta	
         K10 <-CL1 /V1  
         l1<-alpha
         l2<-beta 	
         K21 <-l1*l2/K10 
         K12 <-l1+l2-K21-K10	
         CL2 <-V1 *K12  	
         C1 <- (K21 - l1) / (l2- l1)/V1
         C2 <-(K21 - l2) / (l1 - l2)/V1
         T.A <-C1 
         T.B <-C2 	
         F.A <-C1*V1  
         F.B <-C2*V1  
         V1 <-V1  
         V2 <-V1 *K12 /K21     
         Vdss <-V1 +V2 	
         Result.P<-rbind(V1,CL1,t.0.5.alpha,t.0.5.beta,NA,
                        alpha,beta, V2,Vdss,CL2,K10,K12,K21,T.A,T.B,F.A,F.B)
         return(Result.P)
     }
     
     convert.comp3.m3<-function(input)
     {   V1<-input$V1 
         Vdss<-input$Vdss
         CL1<-input$CL1
         t.0.5.alpha<-input$t.0.5.alpha
         t.0.5.beta<-input$t.0.5.beta
         t.0.5.gamma<-input$t.0.5.gamma
         alpha <-log(2)/t.0.5.alpha 
         beta <-log(2)/t.0.5.beta
         gamma <-log(2)/t.0.5.gamma	
         K10 <-CL1 /V1  
         l1<-alpha
         l2<-beta
         l3<-gamma	
         g0 <- l1  * l2  * l3  / K10  
         g1 <- (Vdss - V1 ) / V1  * g0 
         f <- (l1  * l2  + l1  * l3  + l2  * l3  - g1  - g0) / K10  
         root <- sqrt(f * f - g0 * 4) 
         K21w <- (f + root) / 2 
         K31w <- (f - root) / 2 
         h <- (Vdss / V1  - 1) * g0 
         i <- l1  + l2  + l3  - K10  - K21w  - K31w  
         K13w <- (h - i * K31w ) / (K21w  - K31w )  
         K12w <- i - K13w  
         K21 <-K21w 
         K31 <-K31w 
         K12 <-K12w 
         K13 <-K13w	
         CL2 <-V1 *K12  
         CL3 <-V1 *K13  	
         C1 <- (K21  - l1 ) * (K31  - l1 ) / (l1  - l2 ) / (l1  - l3 )/V1  
         C2 <-(K21  - l2 ) * (K31  - l2 ) / (l2  - l1 ) / (l2  - l3 )/V1  
         C3 <-(K21  - l3 ) * (K31  - l3 ) / (l3  - l2 ) / (l3  - l1 )/V1  
         T.A <-C1 
         T.B <-C2 
         T.C <-C3 	
         F.A <-C1*V1  
         F.B <-C2*V1  
         F.C <-C3*V1  
         V1 <-V1  
         V2 <-V1 *K12 /K21  
         V3 <-V1 *K13 /K31  
         Vdss <-V1 +V2 +V3  	
         Result.P<-rbind(V1,Vdss,CL1,t.0.5.alpha,t.0.5.beta,t.0.5.gamma,NA,
                        alpha,beta,gamma, V1,V2,V3,Vdss,CL2,CL3,K10,K12,K13,K21,K31,T.A,T.B,T.C,F.A,F.B,F.C)
         return(Result.P)
      }

      ## method 4
      convert.comp1.m4<-function(input)
      {  T.A<-input$A
         alpha<-input$alpha
         l1<-alpha	
         asum <-T.A
         K10w <- l1 
         K10 <-K10w 
         V1 <-1/asum 
         Vdss <-V1 	
         CL1 <-V1 *K10  
         F.A <-T.A*V1  
         t.0.5.alpha <-log(2)/l1  		
         Result.P<-rbind(T.A,alpha,NA,
                                V1,Vdss,CL1,F.A,t.0.5.alpha,K10)
         return(Result.P)
      }

      convert.comp2.m4<-function(input)
      {  T.A<-input$A
         T.B<-input$B
         alpha<-input$alpha
         beta<-input$beta
         l1<-alpha
         l2<-beta	
         asum <-T.A+T.B 
         K21w <- (T.A * l2 + T.B * l1) / (T.A + T.B)
         K10w <- l1 * l2/ K21w
         K12w <- l1 + l2 - K21w - K10w
         K10 <-K10w 
         K12 <-K12w 
         K21 <-K21w 
         V1 <-1/asum 
         V2 <-V1 *K12 /K21  
         Vdss <-V1 +V2 	
         CL1 <-V1 *K10  
         CL2 <-V1 *K12   	
         F.A <-T.A*V1  
         F.B <-T.B*V1  
         t.0.5.alpha <-log(2)/l1  
         t.0.5.beta	 <-log(2)/l2  		
         Result.P<-rbind(T.A,T.B,alpha,beta,NA,
                         V1,V2,Vdss,CL1,CL2,F.A,F.B,t.0.5.alpha,t.0.5.beta,K10,K12,K21)
         return(Result.P)
     }

     convert.comp3.m4<-function(input)
     {   T.A<-input$A
         T.B<-input$B
         T.C<-input$C
         alpha<-input$alpha
         beta<-input$beta
         gamma<-input$gamma
         l1<-alpha
         l2<-beta
         l3<-gamma	
         asum <-T.A+T.B+T.C  
         btemp <- -(l1  * T.C  + l1  *T.B + l3  * T.A + l3  * T.B + l2  * T.A + l2  * T.C )/asum 
         ctemp <- (l1  * l2  * T.C  + l1  * l3  * T.B + l2  * l3  * T.A) / asum 
         K21w <- (-btemp + sqrt(btemp * btemp - 4 * ctemp)) / 2 
         K31w <-(-btemp - sqrt(btemp * btemp - 4 * ctemp)) / 2 
         K10w <- l1  * l2  * l3  / K21w / K31w 
         K12w <- ((l2  * l3  + l1  * l2  + l1  * l3 ) - K21w * (l1  + l2  + l3 ) - K10w * K31w + K21w * K21w) / (K31w - K21w) 
         K13w <- l1  + l2  + l3  - (K10w + K12w + K21w + K31w) 
         K10 <-K10w 
         K12 <-K12w 
         K13 <-K13w 
         K21 <-K21w 
         K31 <-K31w 
         V1 <-1/asum 
         V2 <-V1 *K12 /K21  
         V3 <-V1 *K13 /K31  
         Vdss <-V1 +V2 +V3  	
         CL1 <-V1 *K10  
         CL2 <-V1 *K12  
         CL3 <-V1 *K13  	
         F.A <-T.A*V1  
         F.B <-T.B*V1  
         F.C <-T.C*V1  
         t.0.5.alpha <-log(2)/l1  
         t.0.5.beta	 <-log(2)/l2  
         t.0.5.gamma <-log(2)/l3  		
         Result.P<-rbind(T.A,T.B,T.C,alpha,beta,gamma,NA,
                         V1,V2,V3,Vdss,CL1,CL2,CL3,F.A,F.B,F.C,t.0.5.alpha,t.0.5.beta,t.0.5.gamma,K10,K12,K13,K21,K31)  
         return(Result.P)
      }

      ## method 5 
      convert.comp1.m5<-function(input)
      {  V1<-input$V1
         alpha<-input$alpha
         l1<-alpha
         K10 <- l1  
         Vdss <-V1
         CL1 <-V1 *K10  
         T.A <- 1/V1  
         F.A <-T.A*V1  
         t.0.5.alpha <-log(2)/l1  
         Result.P<-rbind(V1,alpha,NA, Vdss,CL1,T.A,F.A, t.0.5.alpha,K10)
         return(Result.P)
      }
    
      convert.comp2.m5<-function(input)
      {  V1<-input$V1
         alpha<-input$alpha
         beta<-input$beta
         K21<-input$K21
         l1<-alpha
         l2<-beta
         K10	 <- l1  * l2   / K21  
         K12	 <- l1 + l2 - K21 - K10
         V2	 <-V1 *K12 /K21  
         Vdss	 <-V1+V2	
         CL1	 <-V1 *K10  
         CL2	 <-V1 *K12  	
         T.A	 <- (K21  - l1 ) / (l2  - l1 )/V1  
         T.B	 <-(K21  - l2 ) / (l1  - l2 )/V1  	
         F.A	 <-T.A*V1  
         F.B	 <-T.B*V1  
         t.0.5.alpha <-log(2)/l1  
         t.0.5.beta<-log(2)/l2  	
         Result.P<<-rbind(V1,alpha,beta,K21,NA,
                      V2,Vdss,CL1,CL2,T.A,T.B,F.A,F.B,
                      t.0.5.alpha,t.0.5.beta,K10,K12)
         return(Result.P)
      }

      convert.comp3.m5<-function(input)
      {  V1<-input$V1
         alpha<-input$alpha
         beta<-input$beta
         gamma<-input$gamma
         K21<-input$K21
         K31<-input$K31
         l1<-alpha
         l2<-beta
         l3<-gamma
         K10	 <- l1  * l2  * l3  / K21  / K31 
         K12	 <- ((l2  * l3  + l1  * l2  + l1  * l3 ) - K21  * (l1  + l2  + l3 ) - K10  * K31  + K21  * K21 ) / (K31  - K21 ) 
         K13	 <- l1  + l2  + l3  - (K10  + K12  + K21  + K31 ) 
         V2	 <-V1 *K12 /K21  
         V3	 <-V1 *K13 /K31  
         Vdss	 <-V1+V2+V3	
         CL1	 <-V1 *K10  
         CL2	 <-V1 *K12  
         CL3	 <-V1 *K13  	
         T.A	 <- (K21  - l1 ) * (K31  - l1 ) / (l1  - l2 ) / (l1  - l3 )/V1  
         T.B	 <-(K21  - l2 ) * (K31  - l2 ) / (l2  - l1 ) / (l2  - l3 )/V1  
         T.C	 <-(K21  - l3 ) * (K31  - l3 ) / (l3  - l2 ) / (l3  - l1 )/V1  	
         F.A	 <-T.A*V1  
         F.B	 <-T.B*V1  
         F.C	 <-T.C*V1  
         t.0.5.alpha <-log(2)/l1  
         t.0.5.beta<-log(2)/l2  
         t.0.5.gamma <-log(2)/l3  
		
         Result.P<<-rbind(V1,alpha,beta,gamma,K21,K31,NA,
                      V2,V3,Vdss,CL1,CL2,CL3,T.A,T.B,T.C,F.A,F.B,F.C,
                      t.0.5.alpha,t.0.5.beta,t.0.5.gamma,K10,K12,K13)
         return(Result.P)
      }

      convert.PK<<-gwindow("Convert PK parameters")
      BBgroup<-ggroup(cont=convert.PK,horizontal=TRUE)
      Bgroup1<-ggroup(cont=BBgroup, horizontal=FALSE)
      tmp<-gframe("Compartment",cont=Bgroup1)
      Comp.t<<-gradio(c("1-comp","2-comp","3-comp"),selected=1,cont=tmp,horizontal=TRUE)    
      tmp<-gframe("Parameterization in NONMEM",cont=Bgroup1)
      Method.t<<-gradio(c("Volumes and clearances","V1,Rate Constants","V1,Vdss,Cl,Half-lives",
                        "Coefficients and Exponents","V1,exponents, K21,K31"),selected=1,cont=tmp,horizontal=FALSE)
      tmp<-gframe("Individual PK parameters(runnumber.par)",cont=Bgroup1)
      file.n1<<-gedit(" ",cont=tmp)
      gbutton("Open",handler=select.V1,cont=tmp)
      tmp<-gframe("Individual PK parameters(csv file)",cont=Bgroup1)
      file.n2<<-gedit(" ",cont=tmp)
      gbutton("Open",handler=select.V2,cont=tmp)
   }
   
   #### Explore simulated data #################################################
   simulationD<-function(h,...)
   {  sim2GUI<-function()
      {  win<-gwindow("Summaries of Simulation")
         BigGroup<<-ggroup(cont=win)
         group<<-ggroup(horizontal=FALSE,cont=BigGroup)    
         dir.g1<<-gedit("")
         button.g1<-gbutton("Open simulated data",handler=openF)
         tmp<-gframe("",container=group)
         add(tmp,button.g1)
         add(tmp,dir.g1)
         N.g1<<-gedit("")
         tmp<-gframe("Number of simulations",container=group)
         add(tmp,N.g1)
      }
      
      openF<-function(h,...)
      { SIM.file<<-gfile("Select simulated data",type="open")
         svalue(dir.g1)<-SIM.file
         SIM.var.temp<-colnames(read.table(SIM.file,skip=1,nrows=1,header=T))
         SIM.var<-SIM.var.temp[which(tolower(SIM.var.temp)!="id" &tolower(SIM.var.temp)!="time")]
         tmp<-gframe("Stratification",cont=group)
         str.list<<-gdroplist(c("NONE",SIM.var))
         str.v<<-gedit("")  
         add(tmp,str.list)
         add(tmp,str.v)          
         tmp<-gframe("Simulation summaries",cont=group)
         var.list<<-gdroplist(SIM.var)  
         add(tmp,var.list)       
         Button1<<-gbutton("Calculate summaries",handler=CalcSIM)
         add(tmp,Button1)
         Button2<<-gbutton("Save summaries",handler=SaveSIM)
         add(tmp,Button2)

         tmp<-gframe("Load summaries of simulation",cont=group)
         Button14<<-gbutton("Open",handler=OpenSIM)
         edit14<<-gedit("")
         add(tmp,Button14)
         add(tmp,edit14)
         CI.list<<-gdroplist(c("95%","90%","80%"))
         tmp<-gframe("PI",cont=group)
         add(tmp,CI.list)
         Button<<-gbutton("OK",handler=updatePlot)
         tmp<-gframe("Plot",cont=group)
         add(tmp,Button)  
         tmp<-gframe("median : red",container=group)
         tmp<-gframe("PI : blue",container=group)   
         add(BigGroup,tmp)
         add(BigGroup, ggraphics())       
      }

      CalcSIM<-function(h,...)
      {  SIM.win<-gwindow("Summaries of Simulation progress",width=300,height=50)
         Bgroup<-ggroup(cont=SIM.win,horizontal=FALSE)
         SIM.progress<-gslider(from=0,to=100,by=1,value=0,cont=Bgroup)
         tmp<-gframe("Read line",cont=Bgroup)
         SIM.line<-gedit(" ",width=50)
         add(tmp,SIM.line)
         svalue(SIM.progress)<-0
         n.sim<-as.numeric(svalue(N.g1))
         flag<-TRUE
         i<-0
         while(flag)
         {  i<-i+1
            temp<-readLines(SIM.file,i*1000)
            t.grep<-grep("TABLE",temp)
            if(length(t.grep)>=2)
             { flag<-FALSE
               N.obs<-t.grep[2]-3
             }
         }
         D.data<-read.table(SIM.file,nrows=N.obs,skip=1,header=T)
         
#         if(sum(colnames(D.data)=="MDV")!=0)
#           D.data<-D.data[D.data$MDV!=1,]
#         if(sum(colnames(D.data)=="DV")!=0)
#           D.data<-D.data[D.data$DV!=0,]
   
         if(svalue(str.list)!="NONE")
            D.data<-D.data[D.data[,svalue(str.list)]==as.numeric(svalue(str.v)),]
         TIME.id<-as.numeric(as.character(names(table(round(D.data[,"TIME"],3)))))
         Quantile.tot1<-NULL
         n.data<-N.obs
         NN<-0
         for(i in 1:length(TIME.id))
         {  b<-round(i/length(TIME.id)*100)
            svalue(SIM.progress)<-b

            ID.list<-which(round(D.data[,"TIME"],3)==TIME.id[i]); NN<-NN+length(ID.list)
            con<-file(SIM.file,"r")
            Line<-0
            TT<-list()
            Line<-Line+1+ID.list[1]
            svalue(SIM.line)<-Line
            temp<-scan(con,skip=1+ID.list[1],nlines=1,quiet=TRUE)
            if(length(ID.list)>1)
            {  for(ki in 2:length(ID.list))
               {  Line<-Line+ID.list[ki]-ID.list[ki-1]
                  svalue(SIM.line)<-Line                
                  temp<-rbind(temp,scan(con,skip=ID.list[ki]-ID.list[ki-1]-1,nlines=1,quiet=TRUE))
               }   
            }
            for(k in 2:n.sim)
            {  Line<-Line+(n.data+1)-ID.list[length(ID.list)]+ID.list[1]+1
               svalue(SIM.line)<-Line
               temp<-rbind(temp,scan(con,skip=(n.data+1)-ID.list[length(ID.list)]+ID.list[1],nlines=1,quiet=TRUE))
               if(length(ID.list)>1)
               {  for(ki in 2:length(ID.list))
                  {  Line<-Line+ID.list[ki]-ID.list[ki-1]
                     svalue(SIM.line)<-Line                     
                     temp<-rbind(temp,scan(con,skip=ID.list[ki]-ID.list[ki-1]-1,nlines=1,quiet=TRUE))
                  }   
               }
            }
            close(con)            
            colnames(temp)<-colnames(D.data)
            temp<-data.frame(temp)
            if(sum(colnames(D.data)=="MDV")!=0)
              temp<-temp[temp$MDV!=1,]
            if(sum(colnames(D.data)=="DV")!=0)
              temp<-temp[temp$DV!=0,]            
            sort.DV<-sort(unlist(temp[, svalue(var.list)]))
            n.DV<-length(sort.DV)
            Q02.5<-sort.DV[max(round(n.DV*0.025),1)]
            Q05<-sort.DV[max(round(n.DV*0.05),1)]
            Q10<-sort.DV[max(round(n.DV*0.10),1)]
            Q25<-sort.DV[max(round(n.DV*0.25),1)]
            Q50<-sort.DV[round(n.DV*0.50)]
            Q75<-sort.DV[round(n.DV*0.75)]
            Q90<-sort.DV[round(n.DV*0.90)]   
            Q95<-sort.DV[round(n.DV*0.95)]
            Q97.5<-sort.DV[round(n.DV*0.975)]
            if(sum(colnames(D.data)=="DV")!=0)
               temp.data<-D.data[D.data$DV!=0,]
            N<-length(which(round(temp.data[,"TIME"],3)==TIME.id[i]))
            q.temp<-c(TIME.id[i],N,Q02.5,Q05,Q10,Q25,
                     Q50,Q75,Q90,Q95,Q97.5)
            Quantile.tot1<-rbind(Quantile.tot1,unlist(q.temp))            
         }

         dispose(SIM.win)
         colnames(Quantile.tot1)<-c("TIME","N","Q02.5","Q05","Q10","Q25","Q50","Q75","Q90","Q95","Q97.5")
         Quantile.keep<<-data.frame(Quantile.tot1)
      }

      OpenSIM<-function(h,...)
      {  file.SIM<<-gfile(text="Open summaries of simulation",
              type="open",filter=list("csv files"=list(patterns=c("*.csv"))))
         svalue(edit14)<-file.SIM       
         Quantile.keep<<-read.csv(file.SIM)
      }

      SaveSIM<-function(h,...)
      {  file.SIM<<-gfile(text="Save Predictive check calculation as csv",
                type="save",filter=list("csv files"=list(patterns=c("*.csv"))))
         write.csv(Quantile.keep,paste(file.SIM,".csv",sep=""),row.names=F)
         svalue(edit14)<-paste(file.SIM,".csv",sep="")  
      }

      updatePlot<-function(h,...)
      {  Quantile.tot<-Quantile.keep
         colnames(Quantile.tot)<-NA
         plot.data<-NULL
         for(j in c(3,4,5,7,9,10,11))
            plot.data<-rbind(plot.data,Quantile.tot[,c(1,j)])
         Quantile.tot<-Quantile.keep
         CI.range<-svalue(CI.list)
         plot(plot.data[,1],plot.data[,2],type='n',
                   xlab="TIME",ylab=svalue(var.list),ylim=range(plot.data[,2]),pch=16,cex=0.7,col="grey")
         lines(Quantile.tot$TIME,Quantile.tot$Q50,lty=1,col=2,lwd=2)                 
         if(CI.range=="95%")
         {  lines(Quantile.tot$TIME,Quantile.tot$Q02.5,lty=1,col=4,lwd=2)
            lines(Quantile.tot$TIME,Quantile.tot$Q97.5,lty=1,col=4,lwd=2) 
         } else if(CI.range=="90%")   
         {  lines(Quantile.tot$TIME,Quantile.tot$Q05,lty=1,col=4,lwd=2)
            lines(Quantile.tot$TIME,Quantile.tot$Q95,lty=1,col=4,lwd=2) 
         } else if(CI.range=="80%")   
         {  lines(Quantile.tot$TIME,Quantile.tot$Q10,lty=1,col=4,lwd=2)
            lines(Quantile.tot$TIME,Quantile.tot$Q90,lty=1,col=4,lwd=2)
         }         
      }
      
      sim2GUI()
   } 
      
   #### Open R terminal with Xpose #############################################
   OpenXpose<-function(h,...)
   {  setwd("c:/fit4NM")
      system("Rgui")
   } 
 
   ############################################################################# 
   # Model evaluation
   #############################################################################

   #### Notes before model evaluations #########################################
   VPCNote<-function(h,...)
   {  gmessage("*** Notes before model evaluations  ***
                \nRadomization test, predictive checks and bootstrap should be conducted within a runnumber subfolder.
                \nRandomization test permutes the median value of a time invariant covariate across every time point in an individual.
                \nSimilary, the median value of a time varying covariate in an individual is permuted. 
                \nThe result of randomization test for time varying covariates should not be adopted. 
                \nFor bootstrap!
                \nSemicolons used for comments in a NONMEM control stream should be preceded by one space as follows.           
                \n   $THETA ; #8     (O)
                \n   $THETA; #8      (X)
                \nFor, predictive checks!
                \n     Please copy the original control file to other folder.
                \n     Rename it to anyname.ctl and move it to runnumber subfolder.             
                \n     All model parameters (THETA, ETA, EPSILON) in anyname.ctl should be replaced by final estimates.          
                \n     Otherwise, leave it as was in the original control file.                
                \n     After starting predictive checks, $EST, $COV, $TABLE blocks are removed and $TABLE block for PC.sim is inserted automatically.                   
                \n     Press \"calculate predictive checks\" button for every stratum when you use stratification.
                \n     All the variables used for stratification should be present at every time point. Please remember this, especially for the stratification based on the dosing data which are not present at every time point in data input file.
                \nFor the menus of open randomization test result, predictive checks with PC.sim and summary data from joined bootstrap raw data file.
                \n     All the conditions but seed number, under which the results of randomization test, predictive checks and bootstrap were obtained, should not be changed.
                \n     Each result file of the randomization test performed in separate PC contains information of the covariate model at the record of prob 0. When you join these result files, leave only one prob 0 data record.",cont=TRUE,width=500)   
   }  
   
   #### Randomization test #####################################################
   RandomTest<-function(h,...)
   {  openControl<-function(h,...)
      {  control.file<-gfile(text="Open control file",type="open")
         current.ctl<<-readLines(control.file)
         svalue(control.t)<-control.file
         temp<-strsplit(control.file,split="\\\\")[[1]]
         Random.RUN<-temp[length(temp)]
         setwd(strsplit(control.file,split=Random.RUN)[[1]])
         Random.RUN<<-strsplit(Random.RUN,split="\\.")[[1]][1]
         file.id<<-strsplit(Random.RUN,split="\\.")[[1]][1]
      }
 
      randomsave<-function(h,...)
      {  file.name<-paste(gfile(text="Save as csv",
                 type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep="")
         write.csv(RAN.result,file.name,row.names=F)
         svalue(edit5)<-file.name                 
      }

      randomopen<-function(h,...)
      {  file.name<-gfile(text="Open randomization result file", type="open",filter=list("csv files"=list(patterns=c("*.csv"))))
         temp<-read.csv(file.name)
         sort.id<-sort.list(as.numeric(as.character(temp$Delta[-1])),decreasing=T)
         temp<-temp[c(1,sort.id+1),]
         n.row<-nrow(temp)-1
         temp$Prob[-1]<-1:n.row
         temp$Quantile[-1]<-c(1:n.row)/n.row
         RAN.result<<-temp
         svalue(edit5)<-file.name
         Rplot()
      }
      
      opendata<-function(h,...)
      {  data.file<<-gfile(text="Open data file",type="open")
         D.data<<-read.csv(data.file,na.string=".")
         svalue(data.t)<-data.file
         Var.list<-colnames(D.data)
         tmp<-gframe("Covariates",cont=BBgroup)
         cov.t<<-gdroplist(Var.list)
         add(tmp,cov.t)
         tmp<-gframe("OBJ for reference model",cont=BBgroup)
         base.obj.t<<-gedit(" ",width=50)
         add(tmp,base.obj.t)
         tmp<-gframe("# of iterations / Seed number",cont=BBgroup)
         RT.label<-glabel("# of replicates")
         iteration.n<<-gedit("1000",width=10)
         add(tmp,RT.label)
         add(tmp,iteration.n)
         Seed.label<-glabel("Seed number")
         Seed.input.RT<<-gedit("0",width=10)
         add(tmp,Seed.label)
         add(tmp,Seed.input.RT) 
         tmp<-gframe("Randomization test",cont=BBgroup)
         button3<-gbutton("Start",handler=randomstart,width=20,height=10)
         add(tmp,button3)
         tmp<-gframe("Save randomization test result",cont=BBgroup)
         button4<-gbutton("Save",handler=randomsave,width=20,height=10)
         add(tmp,button4)
          tmp<-gframe("Open randomization test result",cont=BBgroup)
         button5<-gbutton("Open",handler=randomopen,width=20,height=10)
         edit5<<-gedit("",width=50)
         add(tmp,button5)
         add(tmp,edit5)                  
      }

      randomstart<-function(h,...)
      {  cov.T<-svalue(cov.t)
         print(cov.T)
         temp.CTL<-current.ctl
         temp<-strsplit(temp.CTL,split=" ")
         indicator<-NULL
         for(i in 1:length(temp))
            indicator<-rbind(indicator,temp[[i]][1])
         id<-which(indicator=="$DATA")
         temp.CTL[id]<-"$DATA RT.csv"
         write.table(temp.CTL,"RT.ctl",quote=FALSE,row.names=FALSE,col.names=FALSE)
         id.table<-table(D.data$X.ID)
         COV<-NULL
         for(i in 1:length(id.table))
         {  start.id<-ifelse(i==1,1, (id.table[1:(i-1)]+1))
            COV<-c(COV,median(D.data[start.id:sum(id.table[1:i]),cov.T],na.rm=T))
         }
         iteration.N<-as.numeric(svalue(iteration.n))            
         RT.win<-gwindow("Randomization test progress",width=300,height=50)
         RT.progress<-gslider(from=0,to=iteration.N,by=1,value=0,cont=RT.win)
         svalue(RT.progress)<-0

         random.table<-NULL     
         seed<-as.numeric(svalue(Seed.input.RT))
         set.seed(seed)
         for( k in 0:iteration.N)
         {  svalue(RT.progress)<-k
            D.data.t<-D.data
            if(k !=0)
            {  COV.random<-sample(COV,replace=F)
               cov.random<-COV.random[D.data.t$X.ID]
               D.data.t[,cov.T]<-cov.random
               colnames(D.data.t)[1]<-"#ID"
               write.csv(D.data.t,"RT.csv",quote=FALSE,row.names=FALSE,na=".")
            } else
            {  colnames(D.data.t)[1]<-"#ID"
               write.csv(D.data.t,"RT.csv",quote=FALSE,row.names=FALSE,na=".")
            }
            random.command<-paste(Default.NMpath," RT.ctl RT.res")
            system(random.command,invisible=F,show.output.on.console=F)    
            D.LST<-readLines("RT.res")
            D.temp<-matrix(D.LST)
            indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
            min.id<- which(indicator=="0MINIMIZATION")
            min.id<-min.id[length(min.id)]
            Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]
            indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])
            temp<-strsplit(D.LST[length(D.LST)-4],split=" ")[[1]]
            temp<-as.numeric(temp[temp!=""&temp!="+"])
            cond.num<-round(max(temp)/min(temp),3)
            obj.id<-which(indicator==" #OBJV")
            if(length(obj.id)!=0)
            {  obj.id<-obj.id[length(obj.id)]
            } else
            {  obj.id<-9+which(indicator==" ********************                           MINIMUM VALUE OF OBJECTIVE FUNCTION                  ********************" )
            } 
            temp<-strsplit(D.LST[obj.id],split=" ")[[1]]
            temp<-temp[3:(length(temp)-3)]
            temp<-as.numeric(temp[temp!=""])
            Obj<-temp[!is.na(temp)]            
            final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
            final.start.id<-final.start.id[length(final.start.id)]
            Result.LST<-D.LST[final.start.id:length(D.LST)]
            theta.id<-grep("THETA",Result.LST)
            theta.line<-0
            theta.flag<-T
            while(theta.flag)
            {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
               {  theta.line<-theta.line+1
               } else
               {  theta.flag<-FALSE
               }  
            }
            temp<-NULL
            for(i in 1:theta.line)
               temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
            temp<-as.numeric(temp[temp!=""])
            THETA<-temp      
            seTHETA<-rep(NA,length(THETA))  
            if(length(theta.id)!=1)
            {  temp<-NULL
               for(i in 1:theta.line)
                  temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
               temp<-as.numeric(temp[temp!=""])
               seTHETA<-temp
            }       
            omega.id<-grep("OMEGA",Result.LST)
            omega.line<-0
            omega.flag<-T
            while(omega.flag)
            {  if(Result.LST[omega.id[1]+3+omega.line]!=" ")
               {  omega.line<-omega.line+1
               } else
               {  omega.flag<-FALSE
               }  
            }
            temp<-NULL
            for(i in 1:omega.line)
               temp<-c(temp,unlist(strsplit(Result.LST[omega.id[1]+2+i],split=" ")))
            temp<-temp[temp!=""]
            N.eta<-length(temp)
            OMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
            seOMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
            omega.name<-NULL
            id.current<-omega.id[1]+4+omega.line  
            id.current1<-omega.id[2]+4+omega.line        
            for(i in 1:N.eta)
            {  temp<-NULL;temp1<-NULL
               flag<-T
               while(flag)
               {  id.current<-id.current+1;id.current1<-id.current1+1
                  flag<-Result.LST[id.current]!=" "
                  if(flag)
                  {  temp<-c(temp,unlist(strsplit(Result.LST[id.current],split="+ ")))
                     temp1<-c(temp1,unlist(strsplit(Result.LST[id.current1],split="+ ")))            
                  }
               }
               temp<-temp[temp!=""];temp1<-temp1[temp1!=""]
               temp<-temp[temp!="+"] ;temp1<-temp1[temp1!="+"] 
               temp1[ temp1=="........."]<-NA         
               for(j in 1:N.eta)
               {  OMEGA[i,j]<-as.numeric(temp[j])
                  seOMEGA[i,j]<-as.numeric(temp1[j])
                  omega.name<-c(omega.name,paste("OMEGA(",i,"/",j,")",sep=""))
               }  
               id.current<-id.current+1 
               id.current1<-id.current1+1              
            } 
            sigma.id<-grep("SIGMA",Result.LST)
            temp<-unlist(strsplit(Result.LST[sigma.id[1]+6],split="  "))
            temp<-temp[temp!=""]; temp<-temp[temp!="+"]
            SIGMA<-as.numeric(temp)
            seSIGMA<-rep(NA,length(SIGMA))
            if(length(theta.id)!=1)
            {  temp<-unlist(strsplit(Result.LST[sigma.id[2]+6],split="  "))
               temp<-temp[temp!=""]; temp<-temp[temp!="+"]
               seSIGMA<-as.numeric(temp)      
            } 
            names.est<-c(paste("TH",1:length(THETA),sep=""),omega.name,paste("SIGMA",1:length(SIGMA)))      
            EST<-c(THETA,OMEGA,SIGMA)
            SE<-c(seTHETA,seOMEGA,seSIGMA)
            Prob<-k
            if(sum(!is.na(seTHETA))==0)
            {  COV.s<-"NONE"
               cond.num<-NA
            } else
            {  COV.s<-"OK"
            }
            Delta<-as.numeric(svalue(base.obj.t))-Obj
            temp<-c(Prob,Delta,Obj,Min,COV.s,EST,SE)
            random.table<-rbind(random.table,temp)
         }  
         dispose(RT.win)
         colnames(random.table)<-c("Prob","Delta","Obj","Min","Covariance",names.est,paste(names.est,":se",sep=""))
         random.table<-data.frame(random.table)
         sort.id<-sort.list(as.numeric(as.character(random.table$Delta[-1])),decreasing=T)
         random.T<-random.table[-1,]
         random.T<-rbind(random.table[1,],random.T[sort.id,])
         Quantile<-(0:iteration.N)/iteration.N
         random.T<-cbind(random.T[,1:3],Quantile,random.T[,-c(1:3)])
         Quantile<-Quantile[-1]
         data.delta<-as.numeric(as.character(random.table[1,2]))
         Delta<-as.numeric(as.character(random.T[-1,2]))
         plot(as.numeric(Delta),1-Quantile,type='n',xlim=c(min(Delta,0,data.delta),max(data.delta,Delta,0)*1.3),
              ylim=c(0,1),xlab="Delta : OBJ(Reference-test with permuted data)",
               main=paste(file.id,":",cov.T),               
                sub="Green line = OBJ(Reference-test with original data)",ylab="Quantile",axes=F,col="grey10")
         axis(1)
         for(i in 1:length(Delta))
            lines(c(0,Delta[i]),c(1-Quantile[i],1-Quantile[i]),col="grey")
         cutoff.id<-max(round(0.05*iteration.N),1)
         lines(c(0,Delta[cutoff.id]),c(1-Quantile[cutoff.id],1-Quantile[cutoff.id]),col=2,lwd=2)
         text(max(Delta),1-Quantile[cutoff.id],round(Delta[cutoff.id],3),cex=0.6)
         text(Delta[cutoff.id]/2,1-Quantile[cutoff.id],"5%",cex=0.6)
         data.quantile<-1-sum(Delta>data.delta)/iteration.N
         lines(c(0,data.delta),c(data.quantile,data.quantile),col=3,lwd=2)   
         text(data.delta*1.1,data.quantile,round(data.delta,3),cex=0.6)        
         RAN.result<<- random.T  
         tempR<<-random.table   
      }

      Rplot<-function()
      {  random.T<-RAN.result[-1,]
         data.delta<-as.numeric(as.character(RAN.result[1,2]))
         Delta<-as.numeric(as.character(random.T[,2]))
         iteration.N<-as.numeric(svalue(iteration.n))
         cov.T<-svalue(cov.t)
         Quantile<-(1:iteration.N)/iteration.N     
         plot(as.numeric(Delta),1-Quantile,type='n',xlim=c(min(Delta,0,data.delta),max(data.delta,Delta,0)*1.3),
              ylim=c(0,1),xlab="Delta : OBJ(Reference-test with permuted data)",
               main=paste(file.id,":",cov.T),               
                sub="Green line = OBJ(Reference-test with original data)",ylab="Quantile",axes=F,col="grey10")
         axis(1)
         for(i in 1:length(Delta))
            lines(c(0,Delta[i]),c(1-Quantile[i],1-Quantile[i]),col="grey")
         cutoff.id<-max(round(0.05*iteration.N),1)
         lines(c(0,Delta[cutoff.id]),c(1-Quantile[cutoff.id],1-Quantile[cutoff.id]),col=2,lwd=2)
         text(max(Delta),1-Quantile[cutoff.id],round(Delta[cutoff.id],3),cex=0.6)
         text(Delta[cutoff.id]/2,1-Quantile[cutoff.id],"5%",cex=0.6)
         data.quantile<-1-sum(Delta>data.delta)/iteration.N
         lines(c(0,data.delta),c(data.quantile,data.quantile),col=3,lwd=2)   
         text(data.delta*1.1,data.quantile,round(data.delta,3),cex=0.6)        
      }      
      timewin<-gwindow("Randomization test")
      Bgroup<-ggroup(cont=timewin,horizontal=TRUE)
      BBgroup<-ggroup(cont=Bgroup,horizontal=FALSE) 
      tmp<-gframe("Find run number subfolder",cont=BBgroup)
      control.t<-gedit("",width=50)
      button1<-gbutton("Open control file",handler=openControl)
      add(tmp,button1)
      add(tmp,control.t)
      tmp<-gframe("Find run number subfolder",cont=BBgroup)
      button2<-gbutton("Open data files",handler=opendata,width=20,height=10)
      data.t<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,data.t) 
      add(Bgroup,ggraphics())
   }   
   
   #### Predictive checks ######################################################
   NumericalCheck<-function(h,...)
   {  Quantile.tot<-Quantile.keep
      L1<-svalue(VarList.g1)
      L2<-svalue(VarList.g2)
      L3<-svalue(VarList.g3)
         
      DV.data<-read.csv(paste(VPC.RUN,".csv",sep=""),na.string=".")

      if(L1!="NONE")
      {  L1.d<-as.numeric(svalue(id.g1))
         if(L2!="NONE")
         {  L2.d<-as.numeric(svalue(id.g2))
            if(L3!="NONE")
            {  L3.d<-as.numeric(svalue(id.g3))              
               sel.id<-which(DV.data[,L1]==L1.d & DV.data[,L2]!=L2.d&DV.data[,L3]==L3.d)
            } else
            {  sel.id<-which(DV.data[,L1]==L1.d & DV.data[,L2]!=L2.d)
            }    
         } else if(L3!="NONE")
         {  L3.d<-as.numeric(svalue(id.g3))              
            sel.id<-which(DV.data[,L1]==L1.d & DV.data[,L3]==L3.d)
         } else
         {  sel.id<-which(DV.data[,L1]==L1.d)
         }  
      } else if(L3!="NONE")
      {  L3.d<-as.numeric(svalue(id.g3))              
         sel.id<-which(DV.data[,L3]==L3.d)
      } else
      {  sel.id<-1:nrow(DV.data)
      }  
      Dtemp.data<-DV.data[sel.id,]
  
      var.name<-tolower(colnames(Dtemp.data))
      time.id<-which(var.name=="time")
      DV.id<-which(var.name=="dv")
      if(sum(var.name=="mdv")==0)
      {  D.temp<-D.data
      } else
      {  D.temp<-D.data[(D.data[,which(var.name=="mdv")]==0),] 
      }
         
      time.list<-as.numeric(names(table(round(Dtemp.data[,time.id],3))))

      S1<-nrow(Dtemp.data)
      S2<-nrow(D.temp)
      S3<-VPC.N  
      Quantile.tot$TIME<-round(Quantile.tot$TIME,3)
      Dtemp.data$TIME<-round(Dtemp.data$TIME,3)
      NPC.tot<-NULL
      for(i in 1:length(time.list))
      {  temp.id<-which(Dtemp.data[,time.id]==time.list[i])
         temp.id1<-which(Quantile.tot$TIME==time.list[i])        
         S4<-sum(Dtemp.data[temp.id,DV.id]>Quantile.tot$Q50[temp.id1],na.rm=T)
         S5<-sum(Dtemp.data[temp.id,DV.id]<=Quantile.tot$Q50[temp.id1],na.rm=T)
         S6<-sum(Dtemp.data[temp.id,DV.id]>Quantile.tot$Q97.5[temp.id1],na.rm=T)
         S7<-sum(Dtemp.data[temp.id,DV.id]<Quantile.tot$Q02.5[temp.id1],na.rm=T)  
         S8<-sum(Dtemp.data[temp.id,DV.id]>Quantile.tot$Q95[temp.id1],na.rm=T)
         S9<-sum(Dtemp.data[temp.id,DV.id]<Quantile.tot$Q05[temp.id1],na.rm=T)
         S10<-sum(Dtemp.data[temp.id,DV.id]>Quantile.tot$Q90[temp.id1],na.rm=T)
         S11<-sum(Dtemp.data[temp.id,DV.id]<Quantile.tot$Q10[temp.id1],na.rm=T)
         NPC.tot<-rbind(NPC.tot,c(S4,S5,S6,S7,S8,S9,S10,S11))
      }
      NPC.sum<-apply(NPC.tot,2,sum)
      NPC.txt<-paste("---------------------------------------------\n")
      NPC.txt<-paste(NPC.txt,"Number of records        : ",nrow(Dtemp.data),"\n",sep="")
      NPC.txt<-paste(NPC.txt,"Total observations       : ",nrow(D.temp),"\n",sep="")
      NPC.txt<-paste(NPC.txt,"Number of iteration      : ",S3,"\n",sep="")
      NPC.txt<-paste(NPC.txt,"=============================================\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points above the medians : ",NPC.sum[1],
                            "(",round(NPC.sum[1]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points below the medians : ",NPC.sum[2],
                            "(",round(NPC.sum[2]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Ratio of points above to points below : ",
                                 round(NPC.sum[1]/NPC.sum[2],2),"\n",sep="")
      NPC.txt<-paste(NPC.txt,"=============================================\n",sep="")
      NPC.txt<-paste(NPC.txt,"***  95% Prediction interval ***\n",sep="")
      NPC.txt<-paste(NPC.txt,"---------------------------------------------\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points above 95% PI      : ",NPC.sum[3],
                          "(",round(NPC.sum[3]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points below 95% PI      : ",NPC.sum[4],
                          "(",round(NPC.sum[4]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Ratio of points above to points below : ",
                               round(NPC.sum[3]/NPC.sum[4],2),"\n",sep="")
      NPC.txt<-paste(NPC.txt,"---------------------------------------------\n",sep="")
      NPC.txt<-paste(NPC.txt,"***  90% Prediction interval ***\n",sep="")
      NPC.txt<-paste(NPC.txt,"---------------------------------------------\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points above 90% PI      : ",NPC.sum[5],
                          "(",round(NPC.sum[5]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points below 90% PI      : ",NPC.sum[6],
                          "(",round(NPC.sum[6]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Ratio of points above to points below : ",
                               round(NPC.sum[5]/NPC.sum[6],2),"\n",sep="")
      NPC.txt<-paste(NPC.txt,"---------------------------------------------\n",sep="")
      NPC.txt<-paste(NPC.txt,"***  80% Prediction interval ***\n",sep="")
      NPC.txt<-paste(NPC.txt,"---------------------------------------------\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points above 80% PI      : ",NPC.sum[7],
                          "(",round(NPC.sum[7]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Points below 80% PI      : ",NPC.sum[8],
                          "(",round(NPC.sum[8]/nrow(Dtemp.data)*100,1),"%)","\n",sep="")
      NPC.txt<-paste(NPC.txt,"Ratio of points above to points below : ",
                               round(NPC.sum[7]/NPC.sum[8],2),"\n",sep="")
      NPC.txt<-paste(NPC.txt,"---------------------------------------------\n",sep="")
      edit.win<-gwindow("Predictive checks")
      g<-ggroup(horizontal=FALSE,cont=edit.win)
      tmp<-gframe("",container=g)
      a<-gtext(NPC.txt,width=410,height=310,font.attr=c(sizes="large",family="monospace"))
      add(tmp,a)
      NPC<<-NPC.txt
   }
      
   NumericalCheck.save<-function(h,...)
   {  write.table(NPC,paste(gfile(text="Save predictive check result(name.txt)",
            type="save"),".txt",sep=""),row.names=F)
   }   

   updatePlot<-function(h,...)
   {  L1<-svalue(VarList.g1)
      L2<-svalue(VarList.g2)
      L3<-svalue(VarList.g3)
      Quantile.tot<-Quantile.keep

      DV.data<-read.csv(paste(VPC.RUN,".csv",sep=""),na.string=".")
      if(L1!="NONE")
      {  L1.d<-as.numeric(svalue(id.g1))
         if(L2!="NONE")
         {  L2.d<-as.numeric(svalue(id.g2))
            if(L3!="NONE")
            {  L3.d<-as.numeric(svalue(id.g3))              
               sel.id<-which(DV.data[,L1]==L1.d & DV.data[,L2]!=L2.d&DV.data[,L3]==L3.d)
            } else
            {  sel.id<-which(DV.data[,L1]==L1.d & DV.data[,L2]!=L2.d)
            }    
         } else if(L3!="NONE")
         {  L3.d<-as.numeric(svalue(id.g3))              
            sel.id<-which(DV.data[,L1]==L1.d & DV.data[,L3]==L3.d)
         } else
         {  sel.id<-which(DV.data[,L1]==L1.d)
         }  
      } else if(L3!="NONE")
      {  L3.d<-as.numeric(svalue(id.g3))              
         sel.id<-which(DV.data[,L3]==L3.d)
      } else
      {  sel.id<-1:nrow(DV.data)
      }  
      DV.data<-DV.data[sel.id,]
      CI.range<-svalue(CI.list)
      plot(DV.data$TIME,DV.data$DV,type='n',
                   xlab="TIME",ylab="DV",ylim=c(min(c(unlist(Quantile.tot[,3:11]),DV.data$DV),na.rm=T),max(c(unlist(Quantile.tot[,3:11]),DV.data$DV),na.rm=T)))
      points(DV.data$TIME,DV.data$DV)
      if(CI.range=="95%")
      {  lines(Quantile.tot$TIME,Quantile.tot$Q02.5,lty=2,col=2,lwd=2)
         lines(Quantile.tot$TIME,Quantile.tot$Q97.5,lty=2,col=2,lwd=2)
      } else if(CI.range=="90%")   
      {  lines(Quantile.tot$TIME,Quantile.tot$Q05,lty=2,col=2,lwd=2)
         lines(Quantile.tot$TIME,Quantile.tot$Q95,lty=2,col=2,lwd=2)
      } else if(CI.range=="80%")   
      {  lines(Quantile.tot$TIME,Quantile.tot$Q10,lty=2,col=2,lwd=2)
         lines(Quantile.tot$TIME,Quantile.tot$Q90,lty=2,col=2,lwd=2)
      }  
      VPC.data<<-Quantile.tot
   } 
      
   OpenPC<-function(h,...)
   {  file.PC<<-gfile(text="Open predictive check calculation",
          type="open",filter=list("csv files"=list(patterns=c("*.csv"))))
      svalue(edit14)<-file.PC       
      Quantile.keep<<-read.csv(file.PC)
   }     

   SavePC<-function(h,...)
   {  file.PC<<-gfile(text="Save Predictive check calculation as csv",
           type="save",filter=list("csv files"=list(patterns=c("*.csv"))))
      write.csv(Quantile.keep,paste(file.PC,".csv",sep=""),row.names=F)
      svalue(edit14)<-paste(file.PC,".csv",sep="")       
   }

   CalcNPC2<-function(h,...)
   {  VPC.N<-as.numeric(svalue(N.g1))
      temp.id<-which(TOT.RUN$data[,"ID"]==VPC.RUN)
      temp.id<-temp.id[length(temp.id)]
      current.dir<-TOT.RUN$data[temp.id,"path"]
 
      VPC.win<-gwindow("Predictive checks progress",width=300,height=50)
      VPC.progress<-gslider(from=0,to=100,by=1,value=0,cont=VPC.win)
      svalue(VPC.progress)<-0
      n.sim<-VPC.N
      n.data<-nrow(D.data)

      L1<-svalue(VarList.g1)
      L2<-svalue(VarList.g2)
      L3<-svalue(VarList.g3)

      if(L1!="NONE")
      {  L1.d<-as.numeric(svalue(id.g1))
         if(L2!="NONE")
         {  L2.d<-as.numeric(svalue(id.g2))
            if(L3!="NONE")
            {  L3.d<-as.numeric(svalue(id.g3))              
               sel.id<-which(D.data[,L1]==L1.d & D.data[,L2]!=L2.d&D.data[,L3]==L3.d)
            } else
            {  sel.id<-which(D.data[,L1]==L1.d & D.data[,L2]!=L2.d)
            }    
         } else if(L3!="NONE")
         {  L3.d<-as.numeric(svalue(id.g3))              
            sel.id<-which(D.data[,L1]==L1.d & D.data[,L3]==L3.d)
         } else
         {  sel.id<-which(D.data[,L1]==L1.d)
         }  
      } else if(L3!="NONE")
      {  L3.d<-as.numeric(svalue(id.g3))              
         sel.id<-which(D.data[,L3]==L3.d)
      } else
      {  sel.id<-1:nrow(D.data)
      }  
      TIME.id<-as.numeric(as.character(names(table(round(D.data[sel.id,"TIME"],3)))))
      print(length(TIME.id))
      Quantile.tot1<-NULL
      NN<-0
      for(i in 1:length(TIME.id))
      {  b<-round(i/length(TIME.id)*100)
         svalue(VPC.progress)<-b
         ID.list<-sel.id[which(round(D.data[sel.id,"TIME"],3)==TIME.id[i])]; NN<-NN+length(ID.list)
         con<-file("PC.SIM","r")
         TT<-list()
         temp<-scan(con,skip=1+ID.list[1],nlines=1,quiet=TRUE)
         if(length(ID.list)>1)
         {  for(ki in 2:length(ID.list))
               temp<-rbind(temp,scan(con,skip=ID.list[ki]-ID.list[ki-1]-1,nlines=1,quiet=TRUE))
         }
         for(k in 2:n.sim)
         {  temp<-rbind(temp,scan(con,skip=(n.data+1)-ID.list[length(ID.list)]+ID.list[1],nlines=1,quiet=TRUE))
            if(length(ID.list)>1)
            {  for(ki in 2:length(ID.list))
                  temp<-rbind(temp,scan(con,skip=ID.list[ki]-ID.list[ki-1]-1,nlines=1,quiet=TRUE))
            }
         }

         close(con)
         sort.DV<-sort(unlist(temp[,3]))
         n.DV<-length(sort.DV)
         Q02.5<-sort.DV[max(round(n.DV*0.025),1)]
         Q05<-sort.DV[max(round(n.DV*0.05),1)]
         Q10<-sort.DV[max(round(n.DV*0.10),1)]
         Q25<-sort.DV[max(round(n.DV*0.25),1)]
         Q50<-sort.DV[round(n.DV*0.50)]
         Q75<-sort.DV[round(n.DV*0.75)]
         Q90<-sort.DV[round(n.DV*0.90)]   
         Q95<-sort.DV[round(n.DV*0.95)]
         Q97.5<-sort.DV[round(n.DV*0.975)]
         N<-length(ID.list)
         q.temp<-c(TIME.id[i],N,Q02.5,Q05,Q10,Q25,Q50,Q75,Q90,Q95,Q97.5)
         Quantile.tot1<-rbind(Quantile.tot1,unlist(q.temp))            
      }

      dispose(VPC.win)
      colnames(Quantile.tot1)<-c("TIME","N","Q02.5","Q05","Q10","Q25","Q50","Q75","Q90","Q95","Q97.5")
      Quantile.keep<<-data.frame(Quantile.tot1)
      alarm()
   }

   openF<-function(h,...)
   {  VPC.dir<<-gfile("Select run number foler with PC.sim",type="selectdir")
      setwd(VPC.dir)
      svalue(dir.g1)<-VPC.dir
      temp<-strsplit(VPC.dir,split="\\\\")[[1]]
      VPC.RUN<<-temp[length(temp)]
      current.result<<-paste(VPC.RUN,".noh",sep="")
      D.data<-read.csv(paste(VPC.RUN,".csv",sep=""),na.string=".")
       
      Var.Name<-colnames(D.data)
      group1<<-ggroup(horizontal=TRUE,cont=group)   
      From.label.1<-glabel("")
      VarList.g1<<-gdroplist(c("MDV",Var.Name))
      VarList.g2<<-gdroplist(c("NONE",c("NONE",Var.Name)))
      id.g1<<-gedit("0",width=10)
      id.g2<<-gedit(" ",width=10)
      label.g1<<-glabel("Include")
      label.g2<<-glabel("Exclude")
      
      tmp<-gframe("Select data",container=group1)
      add(tmp,label.g1)
      add(tmp,VarList.g1)
      add(tmp,id.g1)
      add(tmp,label.g2)
      add(tmp,VarList.g2)
      add(tmp,id.g2)

      id.g3<<-gedit(" ",width=10)
      VarList.g3<<-gdroplist(c("NONE",c("NONE",Var.Name)))
      tmp<-gframe("Stratification",container=group)
      add(tmp,VarList.g3)
      add(tmp,id.g3)
   
      tmp<-gframe("Calculate predictive checks",cont=group)
      Button2<<-gbutton("OK",handler=CalcNPC2)
      add(tmp,Button2)
         
      tmp<-gframe("Save predictive checks calculation",cont=group)

      Button13<<-gbutton("OK",handler=SavePC)
      add(tmp,Button13)

      tmp<-gframe("Load predictive checks calculation",cont=group)

      Button14<<-gbutton("Open",handler=OpenPC)
      edit14<<-gedit("")
      add(tmp,Button14)
      add(tmp,edit14)
 
      CI.list<<-gdroplist(c("95%","90%","80%"))
      Button<<-gbutton("OK",handler=updatePlot)
      tmp<-gframe("PI",cont=group)
      add(tmp,CI.list)
      tmp<-gframe("Plot",cont=group)
      add(tmp,Button)  
                   
      tmp<-gframe("Summary",cont=group)
      Button3<<-gbutton("OK",handler=NumericalCheck)
      Button3.save<<-gbutton("save",handler=NumericalCheck.save)
      add(tmp,Button3)  
      add(tmp,Button3.save)  
      add(BigGroup,ggraphics())
   }

   VPC.try<-function()
   {  Current.CTL<-current.ctl
      temp<-strsplit(Current.CTL,split=" ")
      indicator<-NULL
      for(i in 1:length(temp))
         indicator<-rbind(indicator,temp[[i]][1])
      t.id<-which(indicator=="$EST" | indicator=="$ESTIMATION")
      VPC.CTL<-Current.CTL[1:(t.id[1]-1)]
      temp.ctl1<-paste("$SIMULATION (20030521) ONLYSIM SUBPROBLEMS=",VPC.N,sep="")
      input.id<-which(indicator=="$INPUT")
      input.CTL<-Current.CTL[input.id]
      temp.ctl2<-"$TABLE ID TIME DV NOAPPEND NOPRINT ONEHEADER FILE=PC.SIM"
      
      VPC.CTL<-c(VPC.CTL,temp.ctl1,temp.ctl2)
      write.table(VPC.CTL,"PC.CTL",quote=FALSE,row.names=FALSE,col.names=FALSE)
      VPC.command<-paste(Default.NMpath," PC.ctl PC.res")
      system(VPC.command,invisible=F,show.output.on.console=F)      
      VPC2GUI()
   }
 
   VPC2GUI<-function()
   {  VPC.N<<-as.numeric(svalue(VPC.input))    
      dispose(vpcwin)
      win<-gwindow("Predictive checks")
      BigGroup<<-ggroup(cont=win)
      group<<-ggroup(horizontal=FALSE,cont=BigGroup)     
      
      dir.g1<<-gedit("")
      button.g1<-gbutton("Open folder",handler=openF)

      tmp<-gframe("Run number folder with PC.sim",container=group)
      add(tmp,button.g1)
      add(tmp,dir.g1)
      N.g1<<-gedit(VPC.N)
      tmp<-gframe("Number of simulated sample",container=group)
      add(tmp,N.g1)
   }
    
   openControl<-function(h,...)
   {  control.file<-gfile(text="Open control file",type="open")
      current.ctl<<-readLines(control.file)
      svalue(control.t)<-control.file
      temp<-strsplit(control.file,split="\\\\")[[1]]
      VPC.RUN<-temp[length(temp)]
      VPC.dir<<-strsplit(control.file,split=VPC.RUN)[[1]]
      setwd(VPC.dir)
      VPC.RUN<<-strsplit(VPC.RUN,split="\\.")[[1]][1]
   }

   opendata<-function(h,...)
   {  data.file<<-gfile(text="Open data file",type="open")
      D.data<<-read.csv(data.file,na.string=".")
      svalue(data.t)<-data.file
   }
   
   set.VPCN<-function(h,...)
   {  VPC.N<<-as.numeric(svalue(VPC.input))
      VPC.try()
   }
   
   VPCwithFile<-function(h,...)
   {  VPC2GUI()
   }
   


   VPC1GUI<-function(h,...)
   {  vpcwin<<-gwindow("Predictive checks")
      BBgroup<-ggroup(cont=vpcwin,horizontal=FALSE)
 
      tmp<-gframe("Find run number subfolder",cont=BBgroup)
      control.t<<-gedit(" ",width=50)
      button1<-gbutton("Open control file",handler=openControl)
      add(tmp,button1)
      add(tmp,control.t)

      tmp<-gframe("Find run number subfolder ",cont=BBgroup)
      button2<-gbutton("Open data file for predictive checks",handler=opendata,width=20,height=10)
      data.t<<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,data.t)
 
      tmp=gframe("Number of simulations",container=BBgroup)
      VPC.label<-glabel("")
      VPC.input<<-gedit("1000",width=10)
      Button<-gbutton("Start predictive checks",handler=set.VPCN)
      add(tmp,VPC.label)
      add(tmp,VPC.input)
      tmp=gframe("",container=BBgroup)
      add(tmp,Button)
      ButtonF<-gbutton("Predictive checks with PC.sim ($TABLE ID TIME  FILE=PC.SIM)",handler=VPCwithFile)
      tmp=gframe("",container=BBgroup)
      add(tmp,ButtonF) 
   }
    
   #### Bootstrap ##############################################################
   show.BTsummary<-function()
   {  
      temp.id<-which(TOT.RUN$data[,"ID"]==Boot.RUN)
      temp.id<-temp.id[length(temp.id)]
      current.dir<-TOT.RUN$data[temp.id,"path"]
      D.LST<-readLines(paste(current.dir,"\\",Boot.RUN,".res",sep=""))
      D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])

      final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
      final.start.id<-final.start.id[length(final.start.id)]
      Result.LST<-D.LST[final.start.id:length(D.LST)]

      theta.id<-grep("THETA",Result.LST)
      theta.line<-0
      theta.flag<-T
      while(theta.flag)
      {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
         {  theta.line<-theta.line+1
         } else
         {  theta.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:theta.line)
         temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
      temp<-as.numeric(temp[temp!=""])
      THETA<-temp
      N.theta<-length(THETA)
      
      seTHETA<-rep(NA,length(THETA))  
      if(length(theta.id)!=1)
      {  temp<-NULL
         for(i in 1:theta.line)
            temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
         temp<-as.numeric(temp[temp!=""])
         seTHETA<-temp
      } 

      omega.id<-grep("OMEGA",Result.LST)
      omega.line<-0
      omega.flag<-T
      while(omega.flag)
      {  if(Result.LST[omega.id[1]+3+omega.line]!=" ")
         {  omega.line<-omega.line+1
         } else
         {  omega.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:omega.line)
         temp<-c(temp,unlist(strsplit(Result.LST[omega.id[1]+2+i],split=" ")))
      temp<-temp[temp!=""]

      N.eta<-length(temp)
      OMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      seOMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      omega.name<-NULL
      id.current<-omega.id[1]+4+omega.line  
      id.current1<-omega.id[2]+4+omega.line  
      
      for(i in 1:N.eta)
      {  temp<-NULL;temp1<-NULL
         flag<-T
         while(flag)
         {  id.current<-id.current+1;id.current1<-id.current1+1
            flag<-Result.LST[id.current]!=" "
            if(flag)
            {  temp<-c(temp,unlist(strsplit(Result.LST[id.current],split="+ ")))
               temp1<-c(temp1,unlist(strsplit(Result.LST[id.current1],split="+ ")))            
            }
         }
         temp<-temp[temp!=""];temp1<-temp1[temp1!=""]
         temp<-temp[temp!="+"] ;temp1<-temp1[temp1!="+"] 
         temp1[ temp1=="........."]<-NA 
        
         for(j in 1:N.eta)
         {  OMEGA[i,j]<-as.numeric(temp[j])
            seOMEGA[i,j]<-as.numeric(temp1[j])
            omega.name<-c(omega.name,paste("OMEGA(",i,"/",j,")",sep=""))
         }  
         id.current<-id.current+1 
         id.current1<-id.current1+1              
      }

      sigma.id<-grep("SIGMA",Result.LST)
      temp<-unlist(strsplit(Result.LST[sigma.id[1]+6],split="  "))
      temp<-temp[temp!=""]; temp<-temp[temp!="+"]
      SIGMA<-as.numeric(temp)
      N.eps<-length(SIGMA)

      s.Boot.total<-NULL
## All
      Boot.keep.t<-matrix(as.numeric(Boot.keep[,-c(1:4)]),nrow=nrow(Boot.keep))
      Boot.summary<-c(unlist(THETA),unlist(OMEGA),unlist(SIGMA))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) mean(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) sd(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.5,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.025,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.975,na.rm=T)}))      
      Boot.summary<-t(Boot.summary)
      Boot.summary<-cbind(c("All",rep(" ",nrow(Boot.summary)-1)),Boot.summary)
      colnames(Boot.summary)<-c("Condition","Estimates","Mean","SD","Median","2.5%","97.5%")
      o.id<-NULL
      for(i in 1:N.eta)
         for(j in i:N.eta)
            o.id<-rbind(o.id,c(i,j))

      select.id<-c(1:N.theta,N.theta+((o.id[,1]-1)*N.eta+o.id[,2]),
                   N.theta+N.eta*N.eta+(1:N.eps)*(1:N.eps))
      s.Boot.summary<-Boot.summary[select.id,]
      rownames(s.Boot.summary)<-c(paste("THETA(",1:N.theta,")",sep=""),
                                  paste("OMEGA(",o.id[,1],"/",o.id[,2],")",sep=""),  
                                  paste("SIGMA(",1:N.eps,"/",1:N.eps,")",sep=""))
      s.Boot.summary<-cbind(s.Boot.summary[,1],rownames(s.Boot.summary),s.Boot.summary[,-1])
      s.Boot.total<-rbind(s.Boot.total,s.Boot.summary)

## Minimization Successful

      Boot.keep.t<-matrix(as.numeric(Boot.keep[Boot.keep[,3]=="SUCCESSFUL",-c(1:4)]),nrow=sum(Boot.keep[,3]=="SUCCESSFUL"))
      Boot.summary<-c(unlist(THETA),unlist(OMEGA),unlist(SIGMA))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) mean(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) sd(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.5,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.025,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.975,na.rm=T)}))      
      Boot.summary<-t(Boot.summary)
      Boot.summary<-cbind(c("Min success",rep(" ",nrow(Boot.summary)-1)),Boot.summary)
      colnames(Boot.summary)<-c("Condition","Estimates","Mean","SD","Median","2.5%","97.5%")
      o.id<-NULL
      for(i in 1:N.eta)
         for(j in i:N.eta)
            o.id<-rbind(o.id,c(i,j))

      select.id<-c(1:N.theta,N.theta+((o.id[,1]-1)*N.eta+o.id[,2]),
                   N.theta+N.eta*N.eta+(1:N.eps)*(1:N.eps))
      s.Boot.summary<-Boot.summary[select.id,]
      rownames(s.Boot.summary)<-c(paste("THETA(",1:N.theta,")",sep=""),
                                  paste("OMEGA(",o.id[,1],"/",o.id[,2],")",sep=""),  
                                  paste("SIGMA(",1:N.eps,"/",1:N.eps,")",sep=""))
      s.Boot.summary<-cbind(s.Boot.summary[,1],rownames(s.Boot.summary),s.Boot.summary[,-1])
      s.Boot.total<-rbind(s.Boot.total,s.Boot.summary)

## Minimization Successful & covariance OK

      Boot.keep.t<-matrix(as.numeric(Boot.keep[Boot.keep[,3]=="SUCCESSFUL"&Boot.keep[,4]=="OK",-c(1:4)]),nrow=sum(Boot.keep[,3]=="SUCCESSFUL"&Boot.keep[,4]=="OK"))
      Boot.summary<-c(unlist(THETA),unlist(OMEGA),unlist(SIGMA))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) mean(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) sd(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.5,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.025,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.975,na.rm=T)}))      
      Boot.summary<-t(Boot.summary)
      Boot.summary<-cbind(c("COV OK",rep(" ",nrow(Boot.summary)-1)),Boot.summary)
      colnames(Boot.summary)<-c("Condition","Estimates","Mean","SD","Median","2.5%","97.5%")
      o.id<-NULL
      for(i in 1:N.eta)
         for(j in i:N.eta)
            o.id<-rbind(o.id,c(i,j))

      select.id<-c(1:N.theta,N.theta+((o.id[,1]-1)*N.eta+o.id[,2]),
                   N.theta+N.eta*N.eta+(1:N.eps)*(1:N.eps))
      s.Boot.summary<-Boot.summary[select.id,]
      rownames(s.Boot.summary)<-c(paste("THETA(",1:N.theta,")",sep=""),
                                  paste("OMEGA(",o.id[,1],"/",o.id[,2],")",sep=""),  
                                  paste("SIGMA(",1:N.eps,"/",1:N.eps,")",sep=""))
      s.Boot.summary<-cbind(s.Boot.summary[,1],rownames(s.Boot.summary),s.Boot.summary[,-1])
      s.Boot.total<-rbind(s.Boot.total,s.Boot.summary)
      colnames(s.Boot.total)[1:2]<-c("Condition","Parameter")
      win<-gwindow("Bootstrap Summary",width=600,height=300)
      gtable(s.Boot.total,cont=win)
      s.Boot.summary<<-s.Boot.total
  }


   show.BTsummary1<-function(Boot.keep.A,D.LST)
   {  D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
      final.start.id<-final.start.id[length(final.start.id)]
      Result.LST<-D.LST[final.start.id:length(D.LST)]

      theta.id<-grep("THETA",Result.LST)
      theta.line<-0
      theta.flag<-T
      while(theta.flag)
      {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
         {  theta.line<-theta.line+1
         } else
         {  theta.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:theta.line)
         temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
      temp<-as.numeric(temp[temp!=""])
      THETA<-temp
      N.theta<-length(THETA)

      seTHETA<-rep(NA,length(THETA))  
      if(length(theta.id)!=1)
      {  temp<-NULL
         for(i in 1:theta.line)
            temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
         temp<-as.numeric(temp[temp!=""])
         seTHETA<-temp
      } 

      omega.id<-grep("OMEGA",Result.LST)
      omega.line<-0
      omega.flag<-T
      while(omega.flag)
      {  if(Result.LST[omega.id[1]+3+omega.line]!=" ")
         {  omega.line<-omega.line+1
         } else
         {  omega.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:omega.line)
         temp<-c(temp,unlist(strsplit(Result.LST[omega.id[1]+2+i],split=" ")))
      temp<-temp[temp!=""]

      N.eta<-length(temp)
      OMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      seOMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      omega.name<-NULL
      id.current<-omega.id[1]+4+omega.line  
      id.current1<-omega.id[2]+4+omega.line  
      
      for(i in 1:N.eta)
      {  temp<-NULL;temp1<-NULL
         flag<-T
         while(flag)
         {  id.current<-id.current+1;id.current1<-id.current1+1
            flag<-Result.LST[id.current]!=" "
            if(flag)
            {  temp<-c(temp,unlist(strsplit(Result.LST[id.current],split="+ ")))
               temp1<-c(temp1,unlist(strsplit(Result.LST[id.current1],split="+ ")))            
            }
         }
         temp<-temp[temp!=""];temp1<-temp1[temp1!=""]
         temp<-temp[temp!="+"] ;temp1<-temp1[temp1!="+"] 
         temp1[ temp1=="........."]<-NA 
        
         for(j in 1:N.eta)
         {  OMEGA[i,j]<-as.numeric(temp[j])
            seOMEGA[i,j]<-as.numeric(temp1[j])
            omega.name<-c(omega.name,paste("OMEGA(",i,"/",j,")",sep=""))
         }  
         id.current<-id.current+1 
         id.current1<-id.current1+1              
      }

      sigma.id<-grep("SIGMA",Result.LST)
      temp<-unlist(strsplit(Result.LST[sigma.id[1]+6],split="  "))
      temp<-temp[temp!=""]; temp<-temp[temp!="+"]
      SIGMA<-as.numeric(temp)
      N.eps<-length(SIGMA)
      s.Boot.total<-NULL
      
## All
      Boot.keep.t<-matrix(as.numeric(unlist(Boot.keep.A[,-c(1:4)])),nrow=nrow(Boot.keep.A))
      Boot.summary<-c(unlist(THETA),unlist(OMEGA),unlist(SIGMA))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) mean(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) sd(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.5,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.025,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.975,na.rm=T)}))      
      Boot.summary<-t(Boot.summary)
      Boot.summary<-cbind(c("All",rep(" ",nrow(Boot.summary)-1)),Boot.summary)
      colnames(Boot.summary)<-c("Condition","Estimates","Mean","SD","Median","2.5%","97.5%")
      o.id<-NULL
      for(i in 1:N.eta)
         for(j in i:N.eta)
            o.id<-rbind(o.id,c(i,j))

      select.id<-c(1:N.theta,N.theta+((o.id[,1]-1)*N.eta+o.id[,2]),
                   N.theta+N.eta*N.eta+(1:N.eps)*(1:N.eps))
      s.Boot.summary<-Boot.summary[select.id,]
      rownames(s.Boot.summary)<-c(paste("THETA(",1:N.theta,")",sep=""),
                                  paste("OMEGA(",o.id[,1],"/",o.id[,2],")",sep=""),  
                                  paste("SIGMA(",1:N.eps,"/",1:N.eps,")",sep=""))
      s.Boot.summary<-cbind(s.Boot.summary[,1],rownames(s.Boot.summary),s.Boot.summary[,-1])
      s.Boot.total<-rbind(s.Boot.total,s.Boot.summary)

## Minimization Successful

      Boot.keep.t<-matrix(as.numeric(unlist(Boot.keep.A[Boot.keep.A[,3]=="SUCCESSFUL",-c(1:4)])),nrow=sum(Boot.keep.A[,3]=="SUCCESSFUL"))
      Boot.summary<-c(unlist(THETA),unlist(OMEGA),unlist(SIGMA))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) mean(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) sd(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.5,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.025,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.975,na.rm=T)}))      
      Boot.summary<-t(Boot.summary)
      Boot.summary<-cbind(c("Min success",rep(" ",nrow(Boot.summary)-1)),Boot.summary)
      colnames(Boot.summary)<-c("Condition","Estimates","Mean","SD","Median","2.5%","97.5%")
      o.id<-NULL
      for(i in 1:N.eta)
         for(j in i:N.eta)
            o.id<-rbind(o.id,c(i,j))

      select.id<-c(1:N.theta,N.theta+((o.id[,1]-1)*N.eta+o.id[,2]),
                   N.theta+N.eta*N.eta+(1:N.eps)*(1:N.eps))
      s.Boot.summary<-Boot.summary[select.id,]
      rownames(s.Boot.summary)<-c(paste("THETA(",1:N.theta,")",sep=""),
                                  paste("OMEGA(",o.id[,1],"/",o.id[,2],")",sep=""),  
                                  paste("SIGMA(",1:N.eps,"/",1:N.eps,")",sep=""))
      s.Boot.summary<-cbind(s.Boot.summary[,1],rownames(s.Boot.summary),s.Boot.summary[,-1])
      s.Boot.total<-rbind(s.Boot.total,s.Boot.summary)

## Minimization Successful & covariance OK

      Boot.keep.t<-matrix(as.numeric(unlist(Boot.keep.A[Boot.keep.A[,3]=="SUCCESSFUL"&Boot.keep.A[,4]=="OK",-c(1:4)])),
                                        nrow=sum(Boot.keep.A[,3]=="SUCCESSFUL"&Boot.keep.A[,4]=="OK"))
      Boot.summary<-c(unlist(THETA),unlist(OMEGA),unlist(SIGMA))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) mean(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x) sd(x,na.rm=T)))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.5,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.025,na.rm=T)}))
      Boot.summary<-rbind(Boot.summary,apply(Boot.keep.t,2,function(x){quantile(x,probs=0.975,na.rm=T)}))      
      Boot.summary<-t(Boot.summary)
      Boot.summary<-cbind(c("COV OK",rep(" ",nrow(Boot.summary)-1)),Boot.summary)
      colnames(Boot.summary)<-c("Condition","Estimates","Mean","SD","Median","2.5%","97.5%")
      o.id<-NULL
      for(i in 1:N.eta)
         for(j in i:N.eta)
            o.id<-rbind(o.id,c(i,j))

      select.id<-c(1:N.theta,N.theta+((o.id[,1]-1)*N.eta+o.id[,2]),
                   N.theta+N.eta*N.eta+(1:N.eps)*(1:N.eps))
      s.Boot.summary<-Boot.summary[select.id,]
      rownames(s.Boot.summary)<-c(paste("THETA(",1:N.theta,")",sep=""),
                                  paste("OMEGA(",o.id[,1],"/",o.id[,2],")",sep=""),  
                                  paste("SIGMA(",1:N.eps,"/",1:N.eps,")",sep=""))
      s.Boot.summary<-cbind(s.Boot.summary[,1],rownames(s.Boot.summary),s.Boot.summary[,-1])
      s.Boot.total<-rbind(s.Boot.total,s.Boot.summary)
      colnames(s.Boot.total)[1:2]<-c("Condition","Parameter")
      win<-gwindow("Bootstrap Summary",width=600,height=300)
      gtable(s.Boot.total,cont=win)
      s.Boot.summary<<-s.Boot.total
   }
      
   Boot.ctl<-function()
   {  Current.CTL<-current.ctl
      temp.CTL<-Current.CTL
      temp<-strsplit(temp.CTL,split=" ")
      indicator<-NULL
      for(i in 1:length(temp))
        indicator<-rbind(indicator,temp[[i]][1])
      n.CTL<-length(temp.CTL)

      line2<-c("  IF(ICALL.EQ.3) THEN;",    
               "     WRITE (52,*) OBJECT",
               "     WRITE (53,*) THETA",
               "     WRITE (54,*) OMEGA(BLOCK)",
               "     WRITE (55,*) SIGMA(BLOCK)",
               "     DO WHILE(DATA)",
               "       IF(NEWIND.LE.1) WRITE(56,*) ETA",
               "     ENDDO",
               "  ENDIF","",
               "\"LAST")
      id<-which(indicator=="$DATA")
      temp.CTL[id]<-"$DATA boot.csv"
      id<-which(indicator=="$TABLE")
      temp.CTL<-temp.CTL[1:(id[1]-1)]
      id<-which(indicator=="$THETA")
      temp.CTL<-c(temp.CTL[1:(id-1)],line2,
                   temp.CTL[id:length(temp.CTL)])
      Boot.CTL<<-temp.CTL
      write.table(Boot.CTL,"boot.ctl",quote=FALSE,row.names=FALSE,col.names=FALSE)
   }
   
   Boot.try<-function()
   {  var.name<-tolower(colnames(D.data))
      ID.id<-which(var.name=="x.id")
      set.seed(as.numeric(svalue(Seed.input)))
      ID.list<-as.numeric(names(table(D.data[,ID.id])))
      Boot.tot<-NULL
      win<-gwindow(paste("Bootstrap progress : B=",B,sep=""),width=300,height=50)
      Boot.progress<-gslider(from=0,to=B,by=1,value=0,cont=win)
 
      for(b in 1:B)
      {  svalue(Boot.progress)<-b
         Boot.sample.id<-sample(ID.list,length(ID.list),replace=T)
         Boot.sample<-NULL
         for(i in 1:length(Boot.sample.id))
         {  id.temp<-which(D.data[,ID.id]==Boot.sample.id[i])
            id.temp<-D.data[id.temp,]
            id.temp$X.ID<-i
            Boot.sample<-rbind(Boot.sample,id.temp)
         }
         write.table(Boot.sample,"boot.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,na=".",sep=",")
         Boot.command<-paste(Default.NMpath," Boot.CTL Boot.RES")
         system(Boot.command,invisible=F,show.output.on.console=F)      

         D.LST<-readLines("Boot.res")
         D.temp<-matrix(D.LST)
         indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
         term.id<-which(indicator=="0PROGRAM")
         
         if(length(term.id)==0)
         { min.id<- which(indicator=="0MINIMIZATION")
         min.id<-min.id[length(min.id)]
         Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]

         indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])
         temp<-strsplit(D.LST[length(D.LST)-4],split=" ")[[1]]
         temp<-as.numeric(temp[temp!=""&temp!="+"])
         cond.num<-round(max(temp)/min(temp),3)
         obj.id<-which(indicator==" #OBJV")
         if(length(obj.id)!=0)
         {  obj.id<-obj.id[length(obj.id)]
         } else
         {  obj.id<-9+which(indicator==" ********************                           MINIMUM VALUE OF OBJECTIVE FUNCTION                  ********************" )
         } 
         temp<-strsplit(D.LST[obj.id],split=" ")[[1]]
         temp<-temp[3:(length(temp)-3)]
         temp<-as.numeric(temp[temp!=""])
         Obj<-temp[!is.na(temp)]
         
         final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
         final.start.id<-final.start.id[length(final.start.id)]
         Result.LST<-D.LST[final.start.id:length(D.LST)]
         theta.id<-grep("THETA",Result.LST)
         theta.line<-0
         theta.flag<-T
         while(theta.flag)
         {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
            {  theta.line<-theta.line+1
            } else
            {  theta.flag<-FALSE
            }  
         }
         temp<-NULL
         for(i in 1:theta.line)
            temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
         temp<-as.numeric(temp[temp!=""])
         THETA<-temp
      
         seTHETA<-rep(NA,length(THETA))  
         if(length(theta.id)!=1)
         {  temp<-NULL
            for(i in 1:theta.line)
               temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
            temp<-as.numeric(temp[temp!=""])
            seTHETA<-temp
         }   

         omega.id<-grep("OMEGA",Result.LST)
         omega.line<-0
         omega.flag<-T
         while(omega.flag)
         {  if(Result.LST[omega.id[1]+3+omega.line]!=" ")
            {  omega.line<-omega.line+1
            } else
            {  omega.flag<-FALSE
            }  
         }
         temp<-NULL
         for(i in 1:omega.line)
            temp<-c(temp,unlist(strsplit(Result.LST[omega.id[1]+2+i],split=" ")))
         temp<-temp[temp!=""]

         N.eta<-length(temp)
         OMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
         seOMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
         omega.name<-NULL
         id.current<-omega.id[1]+4+omega.line  
         id.current1<-omega.id[2]+4+omega.line  
      
         for(i in 1:N.eta)
         {  temp<-NULL;temp1<-NULL
            flag<-T
            while(flag)
            {  id.current<-id.current+1;id.current1<-id.current1+1
               flag<-Result.LST[id.current]!=" "
               if(flag)
               {  temp<-c(temp,unlist(strsplit(Result.LST[id.current],split="+ ")))
                  temp1<-c(temp1,unlist(strsplit(Result.LST[id.current1],split="+ ")))            
               }
            }
            temp<-temp[temp!=""];temp1<-temp1[temp1!=""]
            temp<-temp[temp!="+"] ;temp1<-temp1[temp1!="+"] 
            temp1[ temp1=="........."]<-NA 
         
            for(j in 1:N.eta)
            {  OMEGA[i,j]<-as.numeric(temp[j])
               seOMEGA[i,j]<-as.numeric(temp1[j])
               omega.name<-c(omega.name,paste("OMEGA(",i,"/",j,")",sep=""))
            }  
            id.current<-id.current+1 
            id.current1<-id.current1+1              
         }      

         sigma.id<-grep("SIGMA",Result.LST)
         temp<-unlist(strsplit(Result.LST[sigma.id[1]+6],split="  "))
         temp<-temp[temp!=""]; temp<-temp[temp!="+"]
         SIGMA<-as.numeric(temp)
         seSIGMA<-rep(NA,length(SIGMA))
         if(length(theta.id)!=1)
         {  temp<-unlist(strsplit(Result.LST[sigma.id[2]+6],split="  "))
            temp<-temp[temp!=""]; temp<-temp[temp!="+"]
            seSIGMA<-as.numeric(temp)      
         } 

         names.est<-c(paste("TH",1:length(THETA),sep=""),omega.name,paste("SIGMA",1:length(SIGMA)))      
         EST<-c(THETA,OMEGA,SIGMA)
         SE<-c(seTHETA,seOMEGA,seSIGMA)
         if(sum(!is.na(seTHETA))==0)
         {  COV.s<-"NONE"
            cond.num<-NA
         } else
         {  COV.s<-"OK"
         }

         N.theta<<-length(THETA)
         N.eta<<-ncol(OMEGA)
         N.eps<<-length(SIGMA)
         tot<-c(b,Obj,Min,COV.s,EST,SE)
         } else
         {   tot<-c(b,NA,"TERMINATED","NONE",rep(NA,ncol(Boot.tot)-4))
         }
         Boot.tot<-rbind(Boot.tot,tot)
      colnames(Boot.tot)<-c("Prob","Obj","Min","Covariance",names.est,paste(names.est,":se",sep=""))
      Boot.keep<<-Boot.tot  
      }

      dispose(win)
      colnames(Boot.tot)<-c("Prob","Obj","Min","Covariance",names.est,paste(names.est,":se",sep=""))
      Boot.keep<<-Boot.tot  
      show.BTsummary()
   }


   Boot.B.init<-function()
   {  setB<-function(h,...)
      {  B<<-as.numeric(svalue(Boot.input))
         Boot.ctl()     
         Boot.try()
      }

      openControl<-function(h,...)
      {  control.file<-gfile(text="Open control file",type="open")
         current.ctl<<-readLines(control.file)
         svalue(control.t)<-control.file
         temp<-strsplit(control.file,split="\\\\")[[1]]
         Boot.RUN<-temp[length(temp)]
         setwd(strsplit(control.file,split=Boot.RUN)[[1]])
         Boot.RUN<<-strsplit(Boot.RUN,split="\\.")[[1]][1]
      }
 
      opendata<-function(h,...)
      {  data.file<<-gfile(text="Open data file",type="open")
         D.data<<-read.csv(data.file)
         svalue(data.t)<-data.file
      }
 
      save1<-function(h,...)
      {  write.csv(Boot.keep,paste(gfile(text="Save bootstrap raw data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }

      save2<-function(h,...)
      {  write.csv(s.Boot.summary,paste(gfile(text="Save bootstrap summary data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }
     
      save3<-function(h,...)
      {  boot.filename<-gfile("Open joined bootstrap raw data file",type="open")
         A<-read.csv(boot.filename)
         D.LST<-readLines(gfile("Open result file(runnumber.res)",type="open"))
         show.BTsummary1(A,D.LST)
      }
     
      bootwin<-gwindow("Bootstrap")
      BBgroup<-ggroup(cont=bootwin,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open control file",handler=openControl)
      add(tmp,button1)
      add(tmp,control.t)

      tmp<-gframe("",cont=BBgroup)
      button2<-gbutton("Open data files",handler=opendata,width=20,height=10)
      data.t<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,data.t)
  
      tmp=gframe("Number of bootstrap replicates / Seed number",container=BBgroup)
      Boot.label<-glabel("# of replicates")
      Boot.input<-gedit("2000",width=10)
      add(tmp,Boot.label)
      add(tmp,Boot.input)
      Seed.label<-glabel("Seed number")
      Seed.input<<-gedit("0",width=10)
      add(tmp,Seed.label)
      add(tmp,Seed.input)     
            
      Button<-gbutton("Start bootstrap",handler=setB)

      tmp=gframe("",container=BBgroup)
      add(tmp,Button)

      tmp<-gframe("",cont=BBgroup)
      Button1<-gbutton("Save bootstrap raw data as csv",handler=save1)
      add(tmp,Button1)
      tmp<-gframe("",cont=BBgroup)
      Button2<-gbutton("Save bootstrap summary data as csv",handler=save2)
      add(tmp,Button2)
      tmp<-gframe("",cont=BBgroup)
      Button3<-gbutton("Summary data from joined bootstrap raw data file (Prob,Obj,Min,COV,EST,SE)",handler=save3)
      add(tmp,Button3)   
      tmp<-gframe("",cont=BBgroup)
      Button4<-gbutton("Save bootstrap summary data from joined bootstrap raw data file",handler=save2)
      add(tmp,Button4)          
   }

   Boot<-function(h,...)
   {  Boot.B.init()
   }

        
                     
   #### Log-likelihood profiling ###############################################
   LLprofiling<-function(h,...)
   {  selRUNnum<-function(h,...)
      {  selectTH<-function(h,...)
         {  EST<-select.EST(D.LST)
            TH.cov<-select.COV(D.LST)           
            TH.est<-EST$TH
            TH.namelist<-names(TH.est)
            

            tot.n.param<- select.n.param(D.LST)
            TH.n<-length(TH.est)
            if(TH.n>10)
            {  TH.select<-c(svalue(THcheck1),svalue(THcheck2),svalue(THcheck3))
            } else
            {  TH.select<-svalue(THcheck)
            }

            TH.est.final<-TH.est[TH.select]
            TH.select.n<-which(TH.namelist==TH.select)
            TH.est.list<-TH.est.final*(1+c(-60,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,60)/100)
            if(as.numeric(svalue(range.E))>60)
              { temp.range<-as.numeric(svalue(range.E))
                temp.int<-as.numeric(svalue(Int.E))
                temp.R<-seq(from=(60+temp.int),to=temp.range,by=temp.int)/100
                TH.est.list<-sort(c((c(-temp.R,temp.R)+1)*TH.est.final,TH.est.list))
              }  
            LLT.FINAL<<-LL.prof(TH.select,TH.select.n,TH.est.list)  
         }

         saveLLT<-function(h,...)
         {   save.LLT<-as.matrix(LLT.FINAL[-1],ncol=1)
             colnames(save.LLT)<-LLT.FINAL[1]
             rownames(save.LLT)<-c("Estimate","2.5%","97.5%")
             save.LLT<-cbind(save.LLT,c(" ","smooth.spline in R","smooth.spline in R"))
             
             write.csv(save.LLT,paste(gfile(text="Save as csv",
             type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""))                  
         }
         saveLLT1<-function(h,...)
         {   write.csv(LLT.save,paste(gfile(text="Save as csv",
             type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""))                  
         }
         
         openLLT<-function(h,...)
         {   LLT.CI<<-read.csv(paste(gfile(text="Open 95% CI",type="open")))   
         }

         openLLT1<-function(h,...)
         {   LLT.tot<<-read.csv(paste(gfile(text="Plot plot-data", type="open")))
                     
         }  
         
         drawLLT<-function(h,...)
         {  plot(LLT.tot$est,LLT.tot$Delta.OBJ,type='l',xlab=colnames(LLT.CI)[2],ylab="Delta.OBJ")
            abline(h=3.84,col=2,lwd=2)
            text(min(LLT.tot$est,na.rm=T),3.84,"3.84",col=2)
            points(LLT.tot$est,LLT.tot$Delta.OBJ,pch=1)
            mid.point<-round((nrow(LLT.tot)+1)/2)
            last.point<-nrow(LLT.tot)            
            L.id<-LLT.tot$est[1:mid.point][max(which(LLT.tot$Delta.OBJ[1:mid.point]>3.84))]
            U.id<-LLT.tot$est[mid.point:last.point][min(which(LLT.tot$Delta.OBJ[mid.point:last.point]>3.84))]
            A<-predict(smooth.spline(LLT.tot$est,LLT.tot$Delta.OBJ),x=seq(L.id,U.id,length.out=100))
            LB<-max(A$x[1:50][A$y[1:50]>3.84])
            UB<-min(A$x[51:100][A$y[51:100]>3.84])
            abline(v=LB,lty=2,col=3)
            abline(v=UB,lty=2,col=3)
            text(LB,sort(LLT.tot$Delta.OBJ,decreasing=T)[2],as.character(round(LB,3)))
            text(UB,sort(LLT.tot$Delta.OBJ,decreasing=T)[2],as.character(round(UB,3)))
         }
            
         LL.prof<-function(EST.select,EST.select.n,EST.list)
         {  LLT.tot<-NULL
            win<-gwindow("Loglikelihood profiling progress ",width=400,height=50)
            LLT.progress<-gslider(from=0,to=length(EST.list),by=1,value=0,cont=win)
            for(i in 1:length(EST.list))
            {  svalue(LLT.progress)<-i
               orig.CTL<-readLines(paste(runnum.path,"\\",FILE.ID,".ctl",sep=""))
               mod.CTL<-addCTL.TH(orig.CTL,EST.select.n,EST.list[i])
               write.table(mod.CTL,"LLP.CTL",quote=FALSE,row.names=FALSE,col.names=FALSE)
               wam.command<-paste(Default.NMpath,"LLP.CTL LLP.RES")
               system(wam.command,invisible=F,show.output.on.console=F)
               LLS.LST<-readLines("LLP.RES")
               obj.v<-select.OBJ(LLS.LST)    
               print(c(EST.list[i],obj.v))         
               LLT.tot<-rbind(LLT.tot,c(EST.list[i],obj.v))
            }
            dispose(LLT.progress)
            colnames(LLT.tot)<-c("est","OBJ")
            LLT.tot<-data.frame(LLT.tot)
            LLT.tot$diffOBJ<-LLT.tot$OBJ-LLT.tot$OBJ[median(1:nrow(LLT.tot))]
            LLT.save<-LLT.tot
            mid.point<-round((nrow(LLT.save)+1)/2)
            last.point<-nrow(LLT.save)
            colnames(LLT.save)<-c("est","OBJ","Delta.OBJ")
            LLT.save<<-LLT.save          
            plot(LLT.tot$est,LLT.tot$diffOBJ,type='l',xlab=EST.select,ylab="Delta.OBJ")
            LLT.save<-LLT.tot
            mid.point<-round((nrow(LLT.save)+1)/2)
            last.point<-nrow(LLT.save)
            colnames(LLT.save)<-c("est","OBJ","Delta.OBJ")
            LLT.save<<-LLT.save
            abline(h=3.84,col=2,lwd=2)
            text(min(LLT.tot$est,na.rm=T),3.84,"3.84",col=2)
            points(LLT.tot$est,LLT.tot$diffOBJ,pch=1)
            L.id<-LLT.tot$est[1:mid.point][max(which(LLT.tot$diffOBJ[1:mid.point]>3.84))]
            U.id<-LLT.tot$est[mid.point:last.point][min(which(LLT.tot$diffOBJ[mid.point:last.point]>3.84))]
            A<-predict(smooth.spline(LLT.tot$est,LLT.tot$diffOBJ),x=seq(L.id,U.id,length.out=100))
            LB<-max(A$x[1:50][A$y[1:50]>3.84])
            UB<-min(A$x[51:100][A$y[51:100]>3.84])
            abline(v=LB,lty=2,col=3)
            abline(v=UB,lty=2,col=3)
            text(LB,sort(LLT.tot$diffOBJ,decreasing=T)[2],as.character(round(LB,3)))
            text(UB,sort(LLT.tot$diffOBJ,decreasing=T)[2],as.character(round(UB,3)))
            
            LLT.result<-c(EST.select,EST.list[8],LB,UB)
            names(LLT.result)<-c("","est","LB","UB")
            return(LLT.result)
         }
               
         file.id<-svalue(id.sel)
         FILE.ID<<-file.id
         runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==file.id),"path"]
         setwd(runnum.path)
         Data.temp<-read.csv(paste(runnum.path,"\\",file.id,".csv",sep=""))
         N<<-nrow(Data.temp)
         D.LST<<-readLines(paste(runnum.path,"\\",file.id,".res",sep=""))
         D.CTL<<-readLines(paste(runnum.path,"\\",file.id,".ctl",sep=""))
         
         EST<-select.EST(D.LST)
         n.theta<-length(EST$TH)
         TH.list<-paste("TH",1:n.theta,sep="")
         if(n.theta<10)
         {  TH.list<-paste("TH",1:n.theta,sep=" ")
         } else
         {  TH.list<-c(paste("TH",1:9,sep=" "),paste("TH",10:n.theta,sep=""))
         }         
         if(length(TH.list)>10)
         { 
           n.temp<-round(length(TH.list)/3)
           THcheck1<<-gcheckboxgroup(TH.list[1:n.temp])
           THcheck2<<-gcheckboxgroup(TH.list[(n.temp+1):(2*n.temp)])
           THcheck3<<-gcheckboxgroup(TH.list[(2*n.temp+1):length(TH.list)])
           Button1<-gbutton("OK",type="OK",handler=selectTH)
           tmp<-gframe(" Select one parameter ",container=group,horizontal=TRUE)
           add(tmp,THcheck1);add(tmp,THcheck2);add(tmp,THcheck3)
         } else  
         { THcheck<<-gcheckboxgroup(TH.list)
           Button1<-gbutton("OK",type="OK",handler=selectTH)
           tmp<-gframe(" Select one parameter ",container=group,horizontal=FALSE)
           add(tmp,THcheck)       
         }
         range.E<<-gedit("60")
         Int.E<<-gedit("10")         
         tmp<-gframe("Maximum percent increase",cont=group)
         add(tmp,range.E)
         tmp<-gframe("Percent interval",cont=group)
         add(tmp,Int.E)           
         tmp<-gframe("Apply log-likelihood profiling method",cont=group)
         add(tmp,Button1)          
         Button2<-gbutton("SAVE",type="SAVE",handler=saveLLT)      
         tmp<-gframe("Save 95% CI as csv",cont=group)
         add(tmp,Button2)   
         Button3<-gbutton("SAVE",type="SAVE",handler=saveLLT1)      
         tmp<-gframe("Save plot-data as csv",cont=group)
         add(tmp,Button3)   
         Button4<-gbutton("Open",type="OPEN",handler=openLLT)      
         tmp<-gframe("Open 95% CI",cont=group)
         add(tmp,Button4)   
         Button5<-gbutton("Open",type="OPEN",handler=openLLT1)      
         tmp<-gframe("Open plot-data",cont=group)
         add(tmp,Button5)   
         Button6<-gbutton("Plot",handler=drawLLT)      
         tmp<-gframe("Log-likelihood profiling plot",cont=group)
         add(tmp,Button6)        
        
         tmp<-gframe("Delta.OBJ = OBJ during profiling - OBJ of selected model",cont=group)
         tmp<-gframe("redline = delta. OBJ 3.84",cont=group)
    
         
      }
      Button<-gbutton("OK",handler=selRUNnum)
      runid.list<-unique(TOT.RUN$data[,"ID"])
      id.sel<-gdroplist(runid.list)

      win<-gwindow("Log-likelihood profiling method")
      BigGroup<<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe("Run number of selected model",container=group)
      add(tmp,id.sel)
      add(tmp,Button)
      add(BigGroup,ggraphics())
   }
   
   #### Wald approximation method ##############################################
   WaldApprox<-function(h,...)
   {  
      selRUNnum<-function(h,...)
      {  selectTH<-function(h,...)
         {  EST<-select.EST(D.LST)
            TH.cov<-select.COV(D.LST)           
            TH.est<-EST$TH
            TH.namelist<-names(TH.est)
            OBJ.full<<-select.OBJ(D.LST)

            tot.n.param<- select.n.param(D.LST)
            TH.n<-length(TH.est)
            if(TH.n>10)
            {  TH.select<-c(svalue(THcheck1),svalue(THcheck2),svalue(THcheck3))
            } else
            {  TH.select<-svalue(THcheck)
            }
            sel.id<-which(duplicated(c(TH.select,TH.namelist))[-c(1:length(TH.select))])          
            TH.cov.final<-as.matrix(TH.cov[sel.id,sel.id]   )
            TH.est.final<-TH.est[sel.id]
            WA.FINAL<<-wa.test(TH.est.final,TH.cov.final,tot.n.param,N)  
#            tmp<-gframe("WAM result",cont=BigGroup,width=600)
#            table1<-gtable(WA.FINAL)
#            add(tmp,table1) 
            gtable(WA.FINAL,cont=gwindow("Top 15 models based on the WAM algorithm",width=600,height=300)) 
         }
         
         integer.base.b <- function(x, b=2)
         {  xi <- as.integer(x)
	    if(any(is.na(xi) | ((x-xi)!=0)))
	    {	print(list(ERROR="x not integer", x=x))
	    }
	    N <- length(x)
	    xMax <- max(x)	
	    ndigits <- (floor(logb(xMax, base=2))+1)
	    Base.b <- array(NA, dim=c(N, ndigits))
	    for(i in 1:ndigits)
	    {  Base.b[, ndigits-i+1] <- (x %% b)
	       x <- (x %/% b)
	    }
	    if(N ==1) Base.b[1, ] else Base.b
         }
 
         saveWA<-function(h,...)
         {  write.csv(WA.FINAL,paste(gfile(text="Save as csv",
            type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)                 
         }
         
         openWA<-function(h,...)
         {  open.WA.t<-read.csv(paste(gfile(text="Open WAM result",
            type="open",filter=list("csv files"=list(patterns=c("*.csv"))))))   
            gtable(open.WA.t,cont=gwindow("Top 15 models based on the WAM algorithm",width=600,height=300))                               
            open.WA<-open.WA.t[which(!is.na(open.WA.t[,7])),]
            temp<-open.WA
            spearman.corr<-cor.test(temp$SBC.W,temp$SBC,method="spearman")$estimate
            pearson.corr<-cor.test(temp$SBC.W,temp$SBC,method="pearson")$estimate      
            plot(as.numeric(open.WA[1:min(nrow(open.WA),15),4]),
                 as.numeric(open.WA[1:min(nrow(open.WA),15),7]),xlab="Actual SBC",ylab="Approximate SBC",
                 xlim=c(0,max(as.numeric(open.WA[1:min(nrow(open.WA),15),4]))),
                 ylim=c(0,max(as.numeric(open.WA[1:min(nrow(open.WA),15),7]))),pch=16,
                 sub=paste("Pearson cor=",round(pearson.corr,3),
                           "  Spearman cor=",round(spearman.corr,3)))
            abline(a=0,b=1)   
            
                      
         }    
              
         wa.test<-function(EST,COV,p,n)
         {  n.param<-length(EST)
            TH.list<-names(EST)#paste("TH",1:n.theta,sep="")
            n.cont<-sum(choose(n.param,1:n.param))
            design.M<-as.matrix(integer.base.b(1:n.cont))
            Name.A<-NULL
            Lambda.A<-NULL
            SBC.A<-NULL
            for(i in 1:n.cont)
            {  sel.par<-as.logical(design.M[i,])
               q<-sum(design.M[i,])
               EST.t<-as.matrix(EST[sel.par])
               COV.t<-as.matrix(COV[sel.par,sel.par])
               if(det(COV.t)<0.00000000001)
                  warnings("Error due to singular covariance matrix")
               Lambda.t<-t(EST.t)%*%solve(COV.t)%*%EST.t
               SBC.t<- Lambda.t+(p-q)*log(n)
               Lambda.A<-c(Lambda.A,Lambda.t)
               SBC.A<-c(SBC.A,SBC.t)
               Name.t<-NULL
               for(j in 1:q)
                 Name.t<-ifelse(j!=1,paste(Name.t,TH.list[sel.par][j],sep=","),
                     paste(Name.t,TH.list[sel.par][j],sep=""))
               Name.A<-c(Name.A,Name.t)  
            }   
            sort.SBC<-sort.list(SBC.A)
            WA.tot.A<-cbind(Name.A[sort.SBC],round(Lambda.A[sort.SBC],3),round(SBC.A[sort.SBC],3))
            design.M<-design.M[sort.SBC,]
            Lambda.T<-NULL
            SBC.T<-NULL
            win<-gwindow(paste("WAM progress : top",min(nrow(WA.tot.A),15)," selections"),width=400,height=50)
            WAM.progress<-gslider(from=0,to=min(nrow(WA.tot.A),15),by=1,value=0,cont=win)
            WA.directory<-getwd()
            for(i in 1:min(nrow(WA.tot.A),15))
            {  setwd(WA.directory)
               svalue(WAM.progress)<-i
               th.list<-WA.tot.A[i,1]
               sel.list<-as.numeric(unlist(strsplit(strsplit(th.list,split="TH")[[1]],split=",")))
               orig.CTL<-readLines(paste(runnum.path,"\\",FILE.ID,".ctl",sep=""))
               mod.CTL<-addCTL.TH(orig.CTL,sel.list,0)
               data.file<-select.datafile.name(mod.CTL)                  
               runID<-paste(FILE.ID,"-WAM-",i,sep="")
               mod.CTL<-changeCTL.id(mod.CTL,FILE.ID,runID)
               mod.CTL.name<-paste(FILE.ID,"-WAM-",i,".ctl",sep="")
               write.table(mod.CTL,mod.CTL.name,quote=FALSE,row.names=FALSE,col.names=FALSE)              
               parent<-FILE.ID
               Description<-paste("WAM-",i,sep="")
               runNONMEM.addRUNtable(runID,mod.CTL.name,data.file,Description,parent)
               wam.LST<-readLines(paste(runID,".res",sep=""))
               n.param<-select.n.param(wam.LST)
               obj.v<-select.OBJ(wam.LST)-OBJ.full             
               sbc.t<-0.5*(obj.v+n.param*log(N))
               Lambda.T<-c(Lambda.T,obj.v)
               SBC.T<-c(SBC.T,sbc.t)
               dir.name<-getwd()
               ETA<-read.table(paste(dir.name,"\\",runID,".ETA",sep=""),skip=1,header=T)
               ETA<-ETA[,-1]
               SaveResult.RES(wam.LST,n.param,N,Description,ETA,runID,dir.name)
            }
#            for(i in 1:min(nrow(WA.tot.A),15))
#            {  svalue(WAM.progress)<-i
#               th.list<-WA.tot.A[i,1]
#               sel.list<-as.numeric(unlist(strsplit(strsplit(th.list,split="TH")[[1]],split=",")))
#               orig.CTL<-readLines(paste(runnum.path,"\\",FILE.ID,".ctl",sep=""))
#               mod.CTL<-addCTL.TH(orig.CTL,sel.list,0)
#               write.table(mod.CTL,"WAM.CTL",quote=FALSE,row.names=FALSE,col.names=FALSE)
#               wam.command<-paste(Default.NMpath,"WAM.CTL WAM.RES")
#               system(wam.command,invisible=F,show.output.on.console=F)
#               wam.LST<-readLines("WAM.RES")
#               n.param<-select.n.param(wam.LST)
#               obj.v<-select.OBJ(wam.LST)-OBJ.full             
#               sbc.t<-obj.v+n.param*log(N)
#               Lambda.T<-c(Lambda.T,obj.v)
#               SBC.T<-c(SBC.T,sbc.t)
#            }            
#fromPDxPOP
#  
#    lrt    <- t(theta2)%*%solve(c2)%*%theta2
#    s      <- p-length(idel)
#    sbc    <- -0.5*(lrt+s*log(n))

# obj <- nres[i,ncols] - objfull
#NLRT[i] <- obj 
#NSBC[i]  <- -0.5*(obj+nparams*log(n))          
            dispose(WAM.progress)
            temp.1<<-WA.tot.A
            temp.2<<-SBC.T
            open.WA<<-cbind(as.numeric(WA.tot.A[,3]),SBC.T)
            open.WA<-data.frame(open.WA)            
            Lambda.T<-c(round(Lambda.T,3),rep(NA,nrow(WA.tot.A)-length(Lambda.T)))
            SBC.T<-c(round(SBC.T,3),rep(NA,nrow(WA.tot.A)-length(SBC.T)))
            Rank.T<-c(rank(SBC.T),rep(NA,nrow(WA.tot.A)-length(SBC.T)))
            WA.tot.A<-cbind(1:length(Lambda.T),WA.tot.A,Rank.T,Lambda.T,SBC.T)
            colnames(WA.tot.A)<-c("WAM rank","Covariate parameters fixed to 0","Lambda.W","SBC.W"," Actual rank","Lambda","SBC")
            WA.tot.A<<-WA.tot.A
            spearman.corr<-cor.test(open.WA[,1],open.WA[,2],method="spearman")$estimate
            pearson.corr<-cor.test(open.WA[,1],open.WA[,2],method="pearson")$estimate      
            
            plot(as.numeric(WA.tot.A[1:min(nrow(WA.tot.A),15),4]),
                 as.numeric(WA.tot.A[1:min(nrow(WA.tot.A),15),7]),xlab="Actual SBC",ylab="Approximate SBC",
                 xlim=c(0,max(as.numeric(WA.tot.A[1:min(nrow(WA.tot.A),15),4]))),
                 ylim=c(0,max(as.numeric(WA.tot.A[1:min(nrow(WA.tot.A),15),7]))),pch=16,
                 sub=paste("Pearson cor=",round(pearson.corr,3),
                           "  Spearman cor=",round(spearman.corr,3)))               
            abline(a=0,b=1)
            return(WA.tot.A)
         }               
         file.id<-svalue(id.sel)
         FILE.ID<<-file.id
         runnum.path<-TOT.RUN$data[which(TOT.RUN$data[,"ID"]==file.id),"path"]
         setwd(runnum.path)
         Data.temp<-read.csv(paste(runnum.path,"\\",file.id,".csv",sep=""))
         N<<-nrow(Data.temp)
         D.LST<<-readLines(paste(runnum.path,"\\",file.id,".res",sep=""))
         EST<-select.EST(D.LST)
         n.theta<-length(EST$TH)
         TH.list<-paste("TH",1:n.theta,sep="")
         if(n.theta<10)
         {  TH.list<-paste("TH",1:n.theta,sep=" ")
         } else
         {  TH.list<-c(paste("TH",1:9,sep=" "),paste("TH",10:n.theta,sep=""))
         }         
         if(length(TH.list)>10)
         {  n.temp<-round(length(TH.list)/3)
            THcheck1<<-gcheckboxgroup(TH.list[1:n.temp])
            THcheck2<<-gcheckboxgroup(TH.list[(n.temp+1):(2*n.temp)])
            THcheck3<<-gcheckboxgroup(TH.list[(2*n.temp+1):length(TH.list)])
            Button1<-gbutton("OK",type="OK",handler=selectTH)
            tmp<-gframe(" Select covariate parameters ",container=group,horizontal=TRUE)
            add(tmp,THcheck1);add(tmp,THcheck2);add(tmp,THcheck3)
            tmp<-gframe("Apply Wald approximation method",cont=group)
            add(tmp,Button1)                    
         } else  
         {  THcheck<<-gcheckboxgroup(TH.list)
            Button1<-gbutton("OK",type="OK",handler=selectTH)
            tmp<-gframe(" Select covariate parameters ",container=group,horizontal=FALSE)
            add(tmp,THcheck)
            tmp<-gframe("Apply Wald approximation method",cont=group)
            add(tmp,Button1)         
         }
         Button2<-gbutton("SAVE",type="SAVE",handler=saveWA)         
         tmp<-gframe("Save as csv",cont=group)         
         add(tmp,Button2)  

         Button3<-gbutton("Open",type="OPEN",handler=openWA)         
         tmp<-gframe("Load WAM result",cont=group)         
         add(tmp,Button3)  
         
         tmp<-gframe("Solid line : line of identity",cont=group)         
      }
      Button<-gbutton("OK",handler=selRUNnum)
      runid.list<-unique(TOT.RUN$data[,"ID"])
      id.sel<-gdroplist(runid.list)

      win<-gwindow("Wald approximation method")
      BigGroup<<-ggroup(cont=win)
      group=ggroup(horizontal=FALSE,cont=BigGroup)
      tmp<-gframe("Run number of full model",container=group)
      add(tmp,id.sel)
      add(tmp,Button)

      add(BigGroup,ggraphics())      
   }
     
   #### External validation of TCI #############################################
   Performance.CCIP<-function(h,...)
   {  save.PE.pop<-function(h,...)
      {  PE.pop.name<-gfile(text="Save measures of performance - population as csv",type="save")
         write.csv(EX.pop,paste(PE.pop.name,".csv",sep=""))
      }

      save.PE.indiv<-function(h,...)
      {  PE.indiv.name<-gfile(text="Save measures of performance - individual as csv",type="save")
         write.csv(EX.Indiv,paste(PE.indiv.name,".csv",sep=""))
      }
      
      save.Raw.PE<-function(h,...)
      {  PE.Raw.name<-gfile(text="Save measures of performance - individual as csv",type="save")
         write.csv(Raw.PE,paste(PE.Raw.name,".csv",sep=""))
      }

      calc.PE<-function(h,...)
      {  ID.list<-NULL
         for(i in 1:4)
            ID.list<-c(ID.list,svalue(PE.list[[i]]))
         EX.data.temp<-EX.data.T[,ID.list]
         colnames(EX.data.temp)<-c("ID","TIME","Cm","Cp")
         EX.data<<-data.frame(EX.data.temp)
         print(dim(EX.data))
         Peformance.Error(EX.data)
         tmp<-gframe("Save external validation - Population",cont=Bgroup1)
         gbutton("OK",handler=save.PE.pop,cont=tmp)
         tmp<-gframe("Save external validation - Individual",cont=Bgroup1)
         gbutton("OK",handler=save.PE.indiv,cont=tmp)
         tmp<-gframe("Save raw data with PE and absolute PE(APE)",cont=Bgroup1)
         gbutton("OK",handler=save.Raw.PE,cont=tmp)         
      }
    
      select.PE<-function(h,...)
      {  noh.name<-gfile(text="data file",type="open")
         svalue(file.PE)<-noh.name
         EX.data.T<<-read.csv(noh.name,na.string=".")
         temp.list<-colnames(EX.data.T)  
         PEparam.input<-c("ID      ","TIME ","Observation    ","Prediction     ")
         PE.g<-list()
         PE.list<-list()
         for(i in 1:4)
         {  PE.g[[i]]<-ggroup(cont=Bgroup1)
            glabel(PEparam.input[i],cont=PE.g[[i]])
            temp<- gdroplist(temp.list,cont=PE.g[[i]])
            if(i==2)
               time.t<-gdroplist(c("sec","min","hour","day"),cont=PE.g[[i]])
            PE.list[[i]]<-temp
         }
         PE.list<<-PE.list
         time.t<<-time.t
         tmp<-gframe("Calculate performance errors",cont=Bgroup1)
         gbutton("OK",handler=calc.PE,cont=tmp)
         tmp<-gframe("Formulae",cont=Bgroup1,horizontal=FALSE)        
         glabel("Performance error(PE) = (observation-prediction)/prediction x 100(%)\n",cont=tmp)                
         glabel("Median performance error(MDPE) = median(PE)\n",cont=tmp)        
         glabel("Median absolute performance error(MDAPE) = median(|PE|)\n",cont=tmp)        
         glabel("Mean performance error(MPE) = mean(PE)\n",cont=tmp)        
         glabel("Mean absolute performance error(MAPE) = mean(|PE|)\n",cont=tmp)                   
         glabel("Divergence = the slope of |PE| ~ time(hour)\n",cont=tmp)        
         glabel("Wobble = median(|PE-individual MDPE|)\n",cont=tmp)              
         glabel("TS : Two-stage approach\n",cont=tmp)
         glabel("PD : Pooled data approach\n",cont=tmp)
         glabel("VW : Variance-weighted approach\n",cont=tmp)       
      }

      Peformance.Error<-function(Ex.data)
      {  ID.list<-unique(Ex.data$ID)
         Ex.data$PE<-ifelse(Ex.data$Cp!=0,(Ex.data$Cm-Ex.data$Cp)/Ex.data$Cp,0)*100
         Ex.data$APE<-abs(Ex.data$PE)
         Raw.PE<<-Ex.data
         Eff<-data.frame(N=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                E=c(1,0.743,0.838,0.697,0.776,0.679,0.743,0.669,0.723,0.663,
                0.709,0.659,0.699,0.656,0.692,0.653,0.686,0.651,0.681))
         EX.Indiv<-NULL
         EX.var<-NULL
         for(i in ID.list)
         {  sel.id<-which(Ex.data$ID==i)
            Ni<-length(sel.id)
            E<-ifelse(Ni<=20,Eff[Eff$N==Ni,2],2/pi)
            PE.temp<-Ex.data$PE[sel.id]
            TIME.temp<-Ex.data$TIME[sel.id]
            MDAPEi<-median(abs(PE.temp),na.rm=T)
            MAPEi<-mean(abs(PE.temp),na.rm=T)
            MDPEi<-median(PE.temp,na.rm=T)
            WOBBLEi<-median(abs(PE.temp-MDPEi),na.rm=T)
            RMSEi<-sqrt(mean(PE.temp^2,na.rm=T))
            MPEi<-mean(PE.temp,na.rm=T)            
            time.unit<-svalue(time.t)
            if(time.unit=="sec")
            {  TT<-TIME.temp/60/60
            } else if(time.unit=="min")
            {  TT<-TIME.temp/60
            } else if(time.unit=="day")
            {  TT<-TIME.temp*24
            } 
            if(Ni!=1)
            {  temp.lm<-summary(lm(abs(PE.temp)~TT))
               DIVERGENCEi<-temp.lm$coef[2,1]
               DIV.var<-temp.lm$coef[2,2]^2
               MDAPE.var<-var(abs(PE.temp),na.rm=T)/(Ni*E)  
               MDPE.var<-var(PE.temp,na.rm=T)/(Ni*E)                          
               WOB.var<-var(abs(PE.temp-MDPEi),na.rm=T)/(Ni*E)                            
            } else
            {  DIVERGENCEi<-0
               DIV.var<-0
               MDAPE.var<-0    
               MDPE.var<-0                          
               WOB.var<-0
            }   

            EX.Indiv<-rbind(EX.Indiv,c(MDAPEi,DIVERGENCEi,MDPEi,WOBBLEi,MAPEi,MPEi))
            EX.var<-rbind(EX.var,c(MDAPE.var,DIV.var,MDPE.var,WOB.var))
         }
         rownames(EX.Indiv)<-ID.list
         colnames(EX.Indiv)<-c("MDAPEi","DIVERGENCEi","MDPEi","WOBBLEi","MAPEi","MPEi")

         rownames(EX.var)<-ID.list
         colnames(EX.var)<-c("MDAPE.var","DIVERGENCE.var","MDPE.var","WOBBLE.var")

         EX.Indiv<-data.frame(EX.Indiv)

         N<-table(Ex.data$ID)
         M<-length(N)
         TS<-apply(EX.Indiv[,1:4],2,mean)
         TS.var<-apply(EX.var,2,mean)/(M^2)
         TS.025<-TS-1.96*sqrt(TS.var)
         TS.975<-TS+1.96*sqrt(TS.var)         
         PD<-apply(EX.Indiv[,1:4],2,function(x) {sum(N*x)/sum(N)})
         PD.var<-apply(EX.var[,1:4],2,function(x) {sum(N*N*x)/(sum(N)^2)})
         PD.025<-TS-1.96*sqrt(PD.var)
         PD.975<-TS+1.96*sqrt(PD.var)           
         Omega<-apply(EX.Indiv[,1:4],2,var)-apply(EX.var,2,mean)
         Omega[Omega<0]<-0
         O.S<-t(as.matrix(Omega)%*%matrix(rep(1,length(N)),nrow=1))+EX.var
         VW<-apply(EX.Indiv[,1:4]/O.S,2,sum)/apply(1/O.S,2,sum)
         VW.var<-apply(EX.var/(O.S*O.S),2,sum)/((apply(1/O.S,2,sum))^2)
         VW.025<-TS-1.96*sqrt(VW.var)
         VW.975<-TS+1.96*sqrt(VW.var)           

         EX.pop<-cbind(TS,sqrt(TS.var),TS.025,TS.975,
                       PD,sqrt(PD.var),PD.025,PD.975,
                       VW,sqrt(VW.var),VW.025,VW.975)
         row.names(EX.pop)<-c("MDAPE","DIVERGENCE","MDPE","WOBBLE")
         colnames(EX.pop)<-c("TS.est","TS.se","TS 2.5%","TS 97.5%",
                             "PD.est","PD.se","PD 2.5%","PD 97.5%",
                             "VW.est","VW.se","VW 2.5%","VW 97.5%")
         EX.pop<-t(EX.pop)
         EX.pop.temp<-cbind(rownames(EX.pop),round(EX.pop,5))

         g1<-gwindow("External validation of TCI - Population")
         gtable(EX.pop.temp,cont=g1)

         EX.Indiv.temp<-cbind(rownames(EX.Indiv),round(EX.Indiv,5))
         colnames(EX.Indiv.temp)<-c("ID",colnames(EX.Indiv))

         g2<-gwindow("External validation of TCI - Individual")
         gtable(EX.Indiv.temp,cont=g2)
         EX.Indiv<<-EX.Indiv
         EX.pop<<-EX.pop
      }

      PE.win<<-gwindow("External validation of TCI")
      BBgroup<-ggroup(cont=PE.win,horizontal=TRUE)
      Bgroup1<-ggroup(cont=BBgroup, horizontal=FALSE)
      tmp<-gframe("Open csv file with ID, TIME, observation, and prediction",cont=Bgroup1)
      file.PE<<-gedit(" ",cont=tmp)
      gbutton("Open",handler=select.PE,cont=tmp)
   }
   
   #### case deletion / jackkinfe  #############################################
   CDD<-function(h,...)
   {
      CDD.ctl<-function()
      {  Current.CTL<-current.ctl
         temp.CTL<-Current.CTL
         temp<-strsplit(temp.CTL,split=" ")
         indicator<-NULL
         for(i in 1:length(temp))
           indicator<-rbind(indicator,temp[[i]][1])
         n.CTL<-length(temp.CTL)

         id<-which(indicator=="$DATA")
         temp.CTL[id]<-"$DATA CDD.csv"
         CV.CTL<-temp.CTL
         write.table(CV.CTL,"CDD.ctl",quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
 

      CDDcalc<-function(h,...)
      {  
         CDD.ctl()
         var.name<-tolower(colnames(D.data))
         ID.id<-which(var.name=="x.id")
           
         ID.list<-unique(D.data[,ID.id])
         N<-length(ID.list)      
         win<-gwindow(paste("Case deletion progress : N=",N,sep=""),width=300,height=50)
         CDD.progress<-gslider(from=0,to=N,by=1,value=0,cont=win)
         D.LST1<<-readLines(paste(CDD.RUN,".RES",sep=""))
         EST<-select.EST(D.LST1)   
         COV<-select.COV(D.LST1)
         OBJ<-select.OBJ(D.LST1)
         CDD.tot<-NULL
         for(k in 1:N)
         {  svalue(CDD.progress)<-k
            RUN.sample<-NULL
            id.temp<-which(D.data[,ID.id]!=ID.list[k])
            RUN.sample<-D.data[id.temp,]
            write.table(RUN.sample,"CDD.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,na=".",sep=",")
            CDD.command<-paste(Default.NMpath," CDD.CTL CDD.RES")
            system(CDD.command,invisible=F,show.output.on.console=F)   
            
            D.LST<<-readLines("CDD.RES")
            RUN.EST<-select.EST(D.LST)   
            RUN.COV<-select.COV(D.LST)
            RUN.SE<-select.SE(D.LST)
            RUN.OBJ<-select.OBJ(D.LST)
            TH.temp<-as.matrix(RUN.EST$TH-EST$TH)
            COV.temp<-COV[1:length(TH.temp),1:length(TH.temp)]
            RUNCOV.temp<-RUN.COV[1:length(TH.temp),1:length(TH.temp)]
            CS<-t(TH.temp)%*%solve(COV.temp)%*%(TH.temp)
            CR<-det(COV.temp)/det(RUNCOV.temp) 
            AddInfo<-select.AddInfo(CDD.RUN)
            
            result<-c(k,AddInfo[c(5,3,4)],c(RUN.EST$TH),c(RUN.EST$OMEGA),c(RUN.EST$SIGMA),
                      c(RUN.SE$TH),c(RUN.SE$OMEGA),c(RUN.SE$SIGMA),CS,CR)   
            CDD.tot<-rbind(CDD.tot,result)  
         }
         omega.name.m<-outer(1:nrow(EST$OMEGA.m),1:nrow(EST$OMEGA.m),function(x,y) paste("OM",x,y,sep=""))
         colnames(CDD.tot)<-c("k","OBJ","MIN","COV", names(RUN.EST$TH),
                       omega.name.m[upper.tri(omega.name.m,diag=TRUE)],names(RUN.EST$SIGMA),
                       paste("SE-",c(names(RUN.EST$TH),omega.name.m[lower.tri(omega.name.m,diag=TRUE)],names(RUN.EST$SIGMA)),sep=""),
                       "CS","CR")
         CDD.tot<-data.frame(CDD.tot)
         CDD.tot<<-CDD.tot
         dispose(win)
         EST.tot<-c(EST$TH,EST$OMEGA,EST$SIGMA)
         EST.CDD<-CDD.tot[,(1:length(EST.tot))+4]
         n.TH<-length(RUN.EST$TH)
         bias<-(N-1)*(apply(EST.CDD,2,function(x) mean(as.numeric(as.character(x))))-EST.tot)
         Jack.est<-N*EST.tot-(N-1)/N*apply(EST.CDD,2,function(x) sum(as.numeric(as.character(x))))
         Jack.var<-sqrt(((N*EST.tot-Jack.est)^2)/(N-1)-2*(N*EST.tot-Jack.est)*apply(EST.CDD,2,function(x) sum(as.numeric(as.character(x))))/N+
                     (N-1)/N*apply(EST.CDD,2,function(x) sum(as.numeric(as.character(x))^2)))
         Jack.L<-EST.tot-qt(0.975,N-1)*sqrt(Jack.var)
         Jack.U<-EST.tot+qt(0.975,N-1)*sqrt(Jack.var)
         CDD.summary<-cbind(bias,EST.tot,Jack.est,Jack.var,Jack.L,Jack.U)   
         colnames(CDD.summary)<-c("Bias","Parameter est","Jackknife est.","Jackkinfe SE","Jackknife 2.5%","Jackknife 97.5%")     
         CDD.summary.P<<-data.frame(Parameter=rownames(CDD.summary),CDD.summary) 
         par(mfrow=c(2,2))
         plot(ID.list,as.numeric(as.character(CDD.tot$CS)),type='n',xlab="ID",ylab="Cook's distance")
         text(ID.list,as.numeric(as.character(CDD.tot$CS)),as.character(ID.list))
         plot(ID.list,as.numeric(as.character(CDD.tot$CR)),type='n',xlab="ID",ylab="Cov Ratio")
         text(ID.list,as.numeric(as.character(CDD.tot$CR)),as.character(ID.list))
         abline(h=1,col=2)
         plot(as.numeric(as.character(CDD.tot$CS)),as.numeric(as.character(CDD.tot$CR)),type='n',xlab="Cook's distance",ylab="Cov Ratio")
         text(as.numeric(as.character(CDD.tot$CS)),as.numeric(as.character(CDD.tot$CR)),as.character(ID.list))         
         abline(h=1,col=2)
         g1<-gwindow("Case deletion summary")         
         gtable(CDD.summary.P,cont=g1)         
      }

      openControl<-function(h,...)
      {  control.file<-gfile(text="Open control file(runnumber subfolder)",type="open")
         current.ctl<<-readLines(control.file)
         svalue(control.t)<-control.file
         temp<-strsplit(control.file,split="\\\\")[[1]]
         CDD.RUN<-temp[length(temp)]
         setwd(strsplit(control.file,split=CDD.RUN)[[1]])
         CDD.RUN<<-strsplit(CDD.RUN,split="\\.")[[1]][1]
      }
 
      opendata<-function(h,...)
      {  data.file<<-gfile(text="Open data file(runnumber subfolder)",type="open")
         D.data<<-read.csv(data.file,na.string=".")
         svalue(data.t)<-data.file
      }
 
      save1<-function(h,...)
      {  write.csv(CDD.tot,paste(gfile(text="Save case deletion raw data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }

      save2<-function(h,...)
      {  write.csv(CDD.summary.P,paste(gfile(text="Save case deletion summary data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }
     
      load1<-function(h,...)
      {  CDD.filename<-gfile("Open case deletion raw data file",type="open")
         CDD.tot<-read.csv(CDD.filename)
         CDD.res<-gfile("Open runnumber.RES file",type="open")
         D.LST1<-readLines(CDD.res)
         EST<-select.EST(D.LST1)   
         COV<-select.COV(D.LST1)
         OBJ<-select.OBJ(D.LST1)         
         EST.tot<-c(EST$TH,EST$OMEGA,EST$SIGMA)
         EST.CDD<-CDD.tot[,(1:length(EST.tot))+4]
         n.TH<-length(EST$TH)
         var.name<-tolower(colnames(D.data))
         ID.id<-which(var.name=="x.id")
           
         ID.list<-unique(D.data[,ID.id])
         N<-length(ID.list)      

         bias<-(N-1)*(apply(EST.CDD,2,function(x) mean(as.numeric(as.character(x))))-EST.tot)
         Jack.est<-N*EST.tot-(N-1)/N*apply(EST.CDD,2,function(x) sum(as.numeric(as.character(x))))
         Jack.var<-sqrt(((N*EST.tot-Jack.est)^2)/(N-1)-2*(N*EST.tot-Jack.est)*apply(EST.CDD,2,function(x) sum(as.numeric(as.character(x))))/N+
                     (N-1)/N*apply(EST.CDD,2,function(x) sum(as.numeric(as.character(x))^2)))
         Jack.L<-EST.tot-qt(0.975,N-1)*sqrt(Jack.var)
         Jack.U<-EST.tot+qt(0.975,N-1)*sqrt(Jack.var)
         CDD.summary<-cbind(bias,EST.tot,Jack.est,Jack.var,Jack.L,Jack.U)   
         colnames(CDD.summary)<-c("Bias","Parameter est","Jackknife est.","Jackkinfe SE","Jackknife 2.5%","Jackknife 97.5%")     
         CDD.summary.P<<-data.frame(Parameter=rownames(CDD.summary),CDD.summary) 
         par(mfrow=c(2,2))
         plot(ID.list,as.numeric(as.character(CDD.tot$CS)),type='n',xlab="ID",ylab="Cook's distance")
         text(ID.list,as.numeric(as.character(CDD.tot$CS)),as.character(ID.list))
         plot(ID.list,as.numeric(as.character(CDD.tot$CR)),type='n',xlab="ID",ylab="Cov Ratio")
         text(ID.list,as.numeric(as.character(CDD.tot$CR)),as.character(ID.list))
         abline(h=1,col=2)
         plot(as.numeric(as.character(CDD.tot$CS)),as.numeric(as.character(CDD.tot$CR)),type='n',xlab="Cook's distance",ylab="Cov Ratio")
         text(as.numeric(as.character(CDD.tot$CS)),as.numeric(as.character(CDD.tot$CR)),as.character(ID.list))         
         abline(h=1,col=2)
         g1<-gwindow("Case deletion summary")         
         gtable(CDD.summary.P,cont=g1)         
      }
     
      CDDwin<-gwindow("Case deletion diagnostics")
      BigGroup<-ggroup(cont=CDDwin,horizontal=TRUE)
      BBgroup<-ggroup(cont=BigGroup,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open control file(runnumber subfolder)",handler=openControl)
      add(tmp,button1)
      add(tmp,control.t)

      tmp<-gframe("",cont=BBgroup)
      button2<-gbutton("Open data files(runnumber subfolder)",handler=opendata,width=20,height=10)
      data.t<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,data.t)
            
      Button<-gbutton("Start case deletion",handler=CDDcalc)

      tmp=gframe("",container=BBgroup)
      add(tmp,Button)

      tmp<-gframe("",cont=BBgroup)
      Button1<-gbutton("Save case deletion raw data as csv",handler=save1)
      add(tmp,Button1)
      tmp<-gframe("",cont=BBgroup)
      Button2<-gbutton("Save case deletion summary data as csv",handler=save2)
      add(tmp,Button2)
      tmp<-gframe("",cont=BBgroup)
      Button4<-gbutton("Load case deletion raw data",handler=load1)
      add(tmp,Button4)
      
      tmp<-gframe("Formulae",cont=BBgroup,horizontal=FALSE)
      glabel("bias = (N-1)*(mean(TH(-k))-TH) \n",cont=tmp)      
      glabel("CS : Cook's distance(k) = sqrt((TH(-k)-TH)*COV(TH)^(-1)*(TH(-k)-TH)) : best=0 \n",cont=tmp)
      glabel("CR : Cov ratio(k) = det(COV(TH))/det(COV(TH(-k))) : best=1\n",cont=tmp)
      glabel("bias = (N-1)*(mean(TH(-k))-TH) \n",cont=tmp)      
      glabel("Jackknife est. = sum(N*TH-(N-1)*TH(-k))/N\n",cont=tmp)
      glabel("Jackknife SE = sqrt(sum((N*TH-(N-1)*TH(-k)-Jackknife est.)**2)/(N*N-N))\n",cont=tmp)
      glabel("Jackknife 95% CI lower bound = Jackknife est. - qt(0.975,N-1)*Jackknife SE\n",cont=tmp)
      glabel("Jackknife 95% CI upper bound = Jackknife est. + qt(0.975,N-1)*Jackknife SE\n",cont=tmp)
     
#      tmp<-gframe("",cont=BBgroup)
#      Button3<-gbutton("Summary data from joined bootstrap raw data file (Prob,Obj,Min,COV,EST,SE)",handler=save3)
#      add(tmp,Button3)   
     add(BigGroup,ggraphics())
     
   }
   #### cross validation  #############################################
   CV<-function(h,...)
   {
      CVRUN.ctl<-function()
      {  Current.CTL<-current.ctl
         temp.CTL<-Current.CTL
         temp<-strsplit(temp.CTL,split=" ")
         indicator<-NULL
         for(i in 1:length(temp))
           indicator<-rbind(indicator,temp[[i]][1])
         n.CTL<-length(temp.CTL)

         id<-which(indicator=="$DATA")
         temp.CTL[id]<-"$DATA RUN.csv"
         CV.CTL<-temp.CTL
         write.table(CV.CTL,"CV-RUN.ctl",quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
      
      CVTEST.ctl<-function()
      {  Current.CTL<-current.ctl
         temp.CTL<-Current.CTL
         temp<-strsplit(temp.CTL,split=" ")
         indicator<-NULL
         for(i in 1:length(temp))
           indicator<-rbind(indicator,temp[[i]][1])
         n.CTL<-length(temp.CTL)

         id<-which(indicator=="$DATA")
         temp.CTL[id]<-"$DATA TEST.csv"
         id1<-which(indicator=="$THETA")
         id2<-which(indicator=="$OMEGA")
         id3<-which(indicator=="$SIGMA")
         id4<-which(indicator=="$ESTIMATION")
         id5<-which(indicator=="$TABLE")
         if((id2-id1-1) > length(EST$TH))
         {  for(i in 1:length(EST$TH))
              temp.CTL[id1+i]<-paste(RUN.EST$TH[i], "  FIX",sep="")        
            if((id1+i+1)<=(id2[1]-1))
               temp.CTL[(id1+i+1):(id2[1]-1)]<-" "  
         } 
         
#         temp.T<-paste("$OMEGA BLOCK(",length(RUN.EST$TH),")",sep="")
#         for(i in 1:length(RUN.EST$OMEGA))
#           temp.T<-paste(temp.T,RUN.EST$OMEGA[i],sep=" ")
#         temp.T<-paste(temp.T," FIX",sep="")
#         temp.CTL[id2[1]]<-temp.T
#         for(i in (id2[1]+1):(id3-1))
#           temp.CTL[i]<-" "  
         OMEGA<-diag(RUN.EST$OMEGA.m)
         temp.T<-"$OMEGA"
         for(i in 1:length(OMEGA))
           temp.T<-paste(temp.T,OMEGA[i], " FIX ",sep=" ")
         temp.CTL[id2[1]]<-temp.T
         for(i in (id2[1]+1):(id3-1))
           temp.CTL[i]<-" "  

         temp.T<-"$SIGMA "
         for(i in 1:length(RUN.EST$SIGMA))
           temp.T<-paste(temp.T,RUN.EST$SIGMA[i],sep=" ")
         temp.T<-paste(temp.T," FIX",sep="")
         temp.CTL[id3]<-temp.T
         for(i in (id3+1):(id4[1]-1))
           temp.CTL[i]<-" "  
         temp.CTL[id4[1]]<-"$ESTIMATION MAXEVAL=0"
         for(i in (id4[1]+1):(id5[1]-1))
           temp.CTL[i]<-" "             
         temp.CTL[id5[1]]<-"$TABLE ID TIME DV IPRED IWRES FILE=CVTEST.LST NOPRINT ONEHEADER"
         temp.CTL<-temp.CTL[1:id5[1]]
         CV.CTL<-temp.CTL       
         write.table(CV.CTL,"CV-TEST.ctl",quote=FALSE,row.names=FALSE,col.names=FALSE)
      }


      CVcalc<-function(h,...)
      {  
         CVRUN.ctl()
         K<<-as.numeric(svalue(CV.input))
         var.name<-tolower(colnames(D.data))
         ID.id<-which(var.name=="x.id")
           
         ID.list<-unique(D.data[,ID.id])
         n<-length(ID.list)      
         if(K < n)
         {  n.rep<-floor(n/K)
            CV.list<-sample(c(rep(c(1:K),n.rep),sample(1:K,(length(ID.list)-n.rep*K),replace=F)))
         } else
         {  CV.list<-sample(1:n,replace=F)         
         }   
         CV.n<-max(CV.list)
         win<-gwindow(paste("Cross validation progress : K=",CV.n,sep=""),width=300,height=50)
         CV.progress<-gslider(from=0,to=CV.n,by=1,value=0,cont=win)
         D.LST<-readLines(paste(CV.RUN,".RES",sep=""))
         EST<<-select.EST(D.LST)   
         COV<-select.COV(D.LST)
         OBJ<-select.OBJ(D.LST)
         CV.tot<-NULL
         for(k in 1:CV.n)
         {  svalue(CV.progress)<-k
            RUN.sample.id<-which(CV.list!=k)
            TEST.sample.id<-which(CV.list==k)            
            RUN.sample<-NULL
            for(i in 1:length(RUN.sample.id))
            {  id.temp<-which(D.data[,ID.id]==RUN.sample.id[i])
               id.temp<-D.data[id.temp,]
               id.temp$X.ID<-i
               RUN.sample<-rbind(RUN.sample,id.temp)
            }
            TEST.sample<-NULL
            for(i in 1:length(TEST.sample.id))
            {  id.temp<-which(D.data[,ID.id]==TEST.sample.id[i])
               id.temp<-D.data[id.temp,]
               id.temp$X.ID<-i
               TEST.sample<-rbind(TEST.sample,id.temp)
            }
            write.table(TEST.sample,"TEST.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,na=".",sep=",")            
            write.table(RUN.sample,"RUN.csv",quote=FALSE,row.names=FALSE,col.names=FALSE,na=".",sep=",")
            CV.command<-paste(Default.NMpath," CV-RUN.CTL CV-RUN.RES")
            system(CV.command,invisible=F,show.output.on.console=F)   
            
            D.LST<<-readLines("CV-RUN.RES")
            D.temp<-readLines("CV-RUN.RES")
            RUN.EST<<-select.EST(D.temp)   
            RUN.COV<-select.COV(D.temp)
            RUN.SE<-select.SE(D.temp)            
            RUN.OBJ<-select.OBJ(D.temp)
            TH.temp<-as.matrix(RUN.EST$TH-EST$TH)
            COV.temp<-COV[1:length(TH.temp),1:length(TH.temp)]
            RUNCOV.temp<-RUN.COV[1:length(TH.temp),1:length(TH.temp)]
            CS<-t(TH.temp)%*%solve(COV.temp)%*%(TH.temp)
            CR<-det(COV.temp)/det(RUNCOV.temp) 
            file.id<-CV.RUN  
            AddInfo<-select.AddInfo(file.id)

            
# PE °è»ê
# MDWR / MDAWR / MAWR / MWR 
            CVTEST.ctl()
            CV.command<-paste(Default.NMpath," CV-TEST.CTL CV-TEST.RES")
            system(CV.command,invisible=F,show.output.on.console=F)   
            CVtest.result<-read.table("CVTEST.LST",skip=1,header=T)
            PE<-(CVtest.result$DV-CVtest.result$IPRED)/CVtest.result$IPRED*100
            PE<-PE[CVtest.result$IPRED!=0]
            MDAWR<-median(abs(PE),na.rm=T)
            MDWR<-median(PE,na.rm=T)
            MAWR<-mean(abs(PE),na.rm=T)
            MWR<-mean(PE,na.rm=T)      
            RMSWR<-sqrt(mean((PE)^2,na.rm=T))    

            result<-c(k,AddInfo[c(5,3,4)],c(RUN.EST$TH),c(RUN.EST$OMEGA),c(RUN.EST$SIGMA),
                      c(RUN.SE$TH),c(RUN.SE$OMEGA),c(RUN.SE$SIGMA),CS,CR,MDAWR,MDWR,MAWR,MWR,RMSWR)   
            CV.tot<-rbind(CV.tot,result)    
            CV.tot.t<<-CV.tot          
         }
         omega.name.m<-outer(1:nrow(EST$OMEGA.m),1:nrow(EST$OMEGA.m),function(x,y) paste("OM",x,y,sep=""))
         omega.name.m[lower.tri(omega.name.m,diag=TRUE)]
         colnames(CV.tot)<-c("k","OBJ","MIN","COV", names(RUN.EST$TH),
                       omega.name.m[upper.tri(omega.name.m,diag=TRUE)],names(RUN.EST$SIGMA),
                       paste("SE-",c(names(RUN.EST$TH),omega.name.m[lower.tri(omega.name.m,diag=TRUE)],names(RUN.EST$SIGMA)),sep=""),
                       "CS","CR","MDAWR","MDWR","MAWR","MWR","RMSWR") 
         CV.tot<-data.frame(CV.tot)
         CV.tot<<-CV.tot

         dispose(win)
         EST.tot<-c(EST$TH,EST$OMEGA,EST$SIGMA)
         EST.CV<-CV.tot[,(1:length(EST.tot))+4]
         n.TH<-length(RUN.EST$TH)
         N<-length(ID.list)
         bias<-(N-1)*(apply(EST.CV,2,function(x) mean(as.numeric(as.character(x))))-EST.tot)
#         Jack.est<-N*EST.tot-(N-1)/N*apply(EST.CV,2,function(x) sum(as.numeric(as.character(x))))
#         Jack.var<-((N*EST.tot-Jack.est)^2)/(N-1)-2*(N*EST.tot-Jack.est)*apply(EST.CV,2,function(x) sum(as.numeric(as.character(x))))/N+
#                     (N-1)/N*apply(EST.CV,2,function(x) sum(as.numeric(as.character(x))^2))
#         Jack.L<-EST.tot-qt(0.975,N-1)*sqrt(Jack.var)
#         Jack.U<-EST.tot+qt(0.975,N-1)*sqrt(Jack.var)
         MDAWR.t<-mean(CV.tot$MDAWR)
         MDWR.t<-mean(CV.tot$MDWR)
         MWR.t<-mean(CV.tot$MWR)
         RMSWR.t<-mean(CV.tot$RMSWR)
         n.last<-n.TH+length(RUN.EST$SIMGA)+length(omega.name.m[lower.tri(omega.name.m,diag=TRUE)])-1
         EST.CV.tot<-apply(matrix(as.numeric(as.character(unlist(
                            CV.tot[,5:(5+n.last)]))),nrow=nrow(CV.tot)),2,mean)
         names(EST.CV.tot)<-colnames(CV.tot)[5:(5+n.last)]
         CV.summary<-cbind(bias,EST.CV.tot,EST.tot)     
         colnames(CV.summary)<-c("Bias","CV parameter est.","Parameter est")     
         CV.summary.P<<-data.frame(parameter=rownames(CV.summary),CV.summary) 
         
         g1<-gwindow("Cross validation summary")         
         gtable(CV.summary.P,cont=g1)                 
      }
      
      

      openControl<-function(h,...)
      {  control.file<-gfile(text="Open control file(runnumber subfolder)",type="open")
         current.ctl<<-readLines(control.file)
         svalue(control.t)<-control.file
         temp<-strsplit(control.file,split="\\\\")[[1]]
         CV.RUN<-temp[length(temp)]
         setwd(strsplit(control.file,split=CV.RUN)[[1]])
         CV.RUN<<-strsplit(CV.RUN,split="\\.")[[1]][1]
      }
 
      opendata<-function(h,...)
      {  data.file<<-gfile(text="Open data file(runnumber subfolder)",type="open")
         D.data<<-read.csv(data.file,na.string=".")
         svalue(data.t)<-data.file
      }
 
      save1<-function(h,...)
      {  write.csv(CV.tot,paste(gfile(text="Save cross validation raw data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }

      save2<-function(h,...)
      {  write.csv(CV.summary.P,paste(gfile(text="Save cross validation summary data as csv",
              type="save",filter=list("csv files"=list(patterns=c("*.csv")))),".csv",sep=""),row.names=F)
      }
     
      load1<-function(h,...)
      {  CV.filename<-gfile("Open case deletion raw data file",type="open")
         CV.tot<-read.csv(CV.filename)
         CV.res<-gfile("Open runnumber.RES file",type="open")
         D.LST1<-readLines(CV.res)
         EST<-select.EST(D.LST1)   
         COV<-select.COV(D.LST1)
         OBJ<-select.OBJ(D.LST1)         
         EST.tot<-c(EST$TH,EST$OMEGA,EST$SIGMA)
         EST.CV<-CV.tot[,(1:length(EST.tot))+4]
         n.TH<-length(EST$TH)
         var.name<-tolower(colnames(D.data))
         ID.id<-which(var.name=="x.id")           
         ID.list<-unique(D.data[,ID.id])         
         N<-length(ID.list)
         bias<-(N-1)*(apply(EST.CV,2,function(x) mean(as.numeric(as.character(x))))-EST.tot)
#         Jack.est<-N*EST.tot-(N-1)/N*apply(EST.CV,2,function(x) sum(as.numeric(as.character(x))))
#         Jack.var<-((N*EST.tot-Jack.est)^2)/(N-1)-2*(N*EST.tot-Jack.est)*apply(EST.CV,2,function(x) sum(as.numeric(as.character(x))))/N+
#                     (N-1)/N*apply(EST.CV,2,function(x) sum(as.numeric(as.character(x))^2))
#         Jack.L<-EST.tot-qt(0.975,N-1)*sqrt(Jack.var)
#         Jack.U<-EST.tot+qt(0.975,N-1)*sqrt(Jack.var)
         MDAWR.t<-mean(CV.tot$MDAWR)
         MDWR.t<-mean(CV.tot$MDWR)
         MWR.t<-mean(CV.tot$MWR)
         RMSWR.t<-mean(CV.tot$RMSWR)
         n.last<-length(bias)-1
         EST.CV.tot<-apply(matrix(as.numeric(as.character(unlist(
                            CV.tot[,5:(5+n.last)]))),nrow=nrow(CV.tot)),2,mean)
         names(EST.CV.tot)<-colnames(CV.tot)[5:(5+n.last)]
         CV.summary<-cbind(bias,EST.CV.tot,EST.tot)     
         colnames(CV.summary)<-c("Bias","CV parameter est.","Parameter est")     
         CV.summary.P<<-data.frame(parameter=rownames(CV.summary),CV.summary) 
         
         g1<-gwindow("Cross validation summary")         
         gtable(CV.summary.P,cont=g1)                 
        }
     
      CVwin<-gwindow("Cross-validation")
      BBgroup<-ggroup(cont=CVwin,horizontal=FALSE)
 
      tmp<-gframe("",cont=BBgroup)
      control.t<-gedit(" ",width=50)
      button1<-gbutton("Open control file(runnumber subfolder)",handler=openControl)
      add(tmp,button1)
      add(tmp,control.t)

      tmp<-gframe("",cont=BBgroup)
      button2<-gbutton("Open data files(runnumber subfolder)",handler=opendata,width=20,height=10)
      data.t<-gedit(" ",width=50)
      add(tmp,button2)
      add(tmp,data.t)
  
      tmp=gframe("Number of subgroups",container=BBgroup)
      CV.label<-glabel("# of subgroups")
      CV.input<-gedit("10",width=10)
      add(tmp,CV.label)
      add(tmp,CV.input)
            
      Button<-gbutton("Start cross-validation",handler=CVcalc)

      tmp=gframe("",container=BBgroup)
      add(tmp,Button)

      tmp<-gframe("",cont=BBgroup)
      Button1<-gbutton("Save cross-validation raw data as csv",handler=save1)
      add(tmp,Button1)
      tmp<-gframe("",cont=BBgroup)
      Button2<-gbutton("Save cross-validation summary data as csv",handler=save2)
      add(tmp,Button2)
      tmp<-gframe("",cont=BBgroup)
      Button4<-gbutton("Load cross-validation raw data",handler=load1)
      add(tmp,Button4)
      
#      tmp<-gframe("",cont=BBgroup)
#      Button3<-gbutton("Summary data from joined bootstrap raw data file (Prob,Obj,Min,COV,EST,SE)",handler=save3)
#      add(tmp,Button3)   
   }  
   ############################################################################# 
   # NONMEM help
   #############################################################################
   #### Open default help ######################################################
   Help1<-function(h,...)
   {  browseURL(Default.Helppath)
   }
   
   #### Open alternative help ##################################################
   Help2<-function(h,...)
   {  browseURL(Alternative.Helppath)
   } 

   ############################################################################# 
   # Additional functions
   #############################################################################
 
   SaveResult.RES<-function(D.LST,param.num,data.n,Description,ETA,file.id,dir.name)
   {  RunID<-file.id
      n.lst<-length(D.LST)
      Date<-D.LST[n.lst-1]
      Time<-D.LST[n.lst]

      D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      min.id<- which(indicator=="0MINIMIZATION")
      min.id<-min.id[length(min.id)]
      Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]

      indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])

      cond.id<-grep(" EIGENVALUES ",D.LST)
      if(length(cond.id)!=0)
      {  
         flag<-T
         id.current<-cond.id+5
         cond.line<-0
         while(flag)
         {  ttemp<-D.LST[id.current]
            flag<-ttemp!=" "
            if(flag)
            {  cond.line<-cond.line+1
               id.current<-id.current+1
            }
         }
         temp<-NULL
         for(i in 1:cond.line)
            temp<-c(temp,strsplit(D.LST[cond.id+5+cond.line+i],split=" ")[[1]])
         temp<-as.numeric(temp[temp!=""&temp!="+"])
         cond.num<-round(max(temp)/min(temp),3)        
      } else
      {  cond.num<-NA
      }   
      
      obj.id<-which(indicator==" #OBJV")
      if(length(obj.id)!=0)
      { obj.id<-obj.id[length(obj.id)]
      } else
      { obj.id<-9+which(indicator==" ********************                           MINIMUM VALUE OF OBJECTIVE FUNCTION                  ********************" )
      } 
      temp<-strsplit(D.LST[obj.id],split=" ")[[1]]
      temp<-temp[3:(length(temp)-3)]
      temp<-as.numeric(temp[temp!=""])
      Obj<-temp[!is.na(temp)]
      AIC<-Obj+2*param.num
      AICc<-round(Obj+2*param.num+2*param.num*(param.num+1)/(data.n-param.num-1),3)
      SBC<-round(Obj+param.num*log(data.n),3)
      parent<-TOT.RUN$data[TOT.RUN$num,"parents"]

# choose THETA,seTHETA,OMEGA,seOMEGA,SIGMA,seSIGMA

      final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
      final.start.id<-final.start.id[length(final.start.id)]
      Result.LST<-D.LST[final.start.id:length(D.LST)]

      theta.id<-grep("THETA",Result.LST)
      theta.line<-0
      theta.flag<-TRUE
      while(theta.flag)
      {  if(Result.LST[theta.id[1]+3+theta.line]!=" ")
         {  theta.line<-theta.line+1
         } else
         {  theta.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:theta.line)
         temp<-c(temp,unlist(strsplit(Result.LST[theta.id[1]+3+theta.line+i],split=" ")))
      temp<-as.numeric(temp[temp!=""])
      THETA<-temp
      
      seTHETA<-rep(NA,length(THETA))  
      if(length(theta.id)!=1)
      {  temp<-NULL
         for(i in 1:theta.line)
            temp<-c(temp,unlist(strsplit(Result.LST[theta.id[2]+3+theta.line+i],split=" ")))
         temp<-as.numeric(temp[temp!=""])
         seTHETA<-temp
      } 

      omega.id<-grep("OMEGA",Result.LST)
      omega.line<-0
      omega.flag<-TRUE
      while(omega.flag)
      {  if(Result.LST[omega.id[1]+3+omega.line]!=" ")
         {  omega.line<-omega.line+1
         } else
         {  omega.flag<-FALSE
         }  
      }
      temp<-NULL
      for(i in 1:omega.line)
         temp<-c(temp,unlist(strsplit(Result.LST[omega.id[1]+2+i],split=" ")))
      temp<-temp[temp!=""]

      N.eta<-length(temp)
      OMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      seOMEGA<-matrix(NA,nrow=N.eta,ncol=N.eta)
      omega.name<-NULL
      id.current<-omega.id[1]+4+omega.line  
      id.current1<-omega.id[2]+4+omega.line  
      
      for(i in 1:N.eta)
      {  temp<-NULL;temp1<-NULL
         flag<-TRUE
         while(flag)
         {  id.current<-id.current+1;id.current1<-id.current1+1
            flag<-Result.LST[id.current]!=" "
            if(flag)
            {  temp<-c(temp,unlist(strsplit(Result.LST[id.current],split="+ ")))
               temp1<-c(temp1,unlist(strsplit(Result.LST[id.current1],split="+ ")))            
            }
         }
         temp<-temp[temp!=""];temp1<-temp1[temp1!=""]
         temp<-temp[temp!="+"] ;temp1<-temp1[temp1!="+"] 
         temp1[ temp1=="........."]<-NA 
         
         for(j in 1:N.eta)
         {  OMEGA[i,j]<-as.numeric(temp[j])
            seOMEGA[i,j]<-as.numeric(temp1[j])
            omega.name<-c(omega.name,paste("OMEGA(",i,"/",j,")",sep=""))
         }  
         id.current<-id.current+1 
         id.current1<-id.current1+1              
      }     

      sigma.id<-grep("SIGMA",Result.LST)
      temp<-unlist(strsplit(Result.LST[sigma.id[1]+6],split="  "))
      temp<-temp[temp!=""]; temp<-temp[temp!="+"]
      SIGMA<-as.numeric(temp)
      seSIGMA<-rep(NA,length(SIGMA))
      if(length(theta.id)!=1)
      {  temp<-unlist(strsplit(Result.LST[sigma.id[2]+6],split="  "))
         temp<-temp[temp!=""]; temp<-temp[temp!="+"]
         seSIGMA<-as.numeric(temp)      
      } 

      names.est<-c(paste("TH",1:length(THETA),sep=""),omega.name,paste("SIGMA",1:length(SIGMA)))      
      EST<-c(THETA,OMEGA,SIGMA)
      SE<-c(seTHETA,seOMEGA,seSIGMA)

      RSE<-round(SE/EST*100,4)
      Lower<-round(EST-1.96*SE,4)
      Upper<-round(EST+1.96*SE,4)

      if(sum(!is.na(seTHETA))==0)
      {  COV<-"NONE"
         cond.num<-NA
      } else
      {  COV<-"OK"
      }

      temp<-c(RunID,Date,Time,Min,COV,Obj,AIC,AICc,SBC,cond.num,parent,Description,param.num)

      se.ETA<-apply(ETA,2,sd) 
      shrinkage.ETA<-matrix(NA,ncol=N.eta,nrow=N.eta)
      diag(shrinkage.ETA)<-(1-se.ETA/sqrt(diag(OMEGA)))*100
      shrinkage<-c(rep("NA",length(THETA)),round(shrinkage.ETA,3),rep("NA",length(SIGMA)))

      CV.ETA<-matrix(NA,ncol=N.eta,nrow=N.eta)
      diag(CV.ETA)<-sqrt(diag(OMEGA))*100
      CV<-c(rep("NA",length(THETA)),round(CV.ETA,3),rep("NA",length(SIGMA)))
      
      tot.res<-cbind(names.est,EST,SE,RSE,Lower,Upper,shrinkage,CV)
      colnames(tot.res)<-c("Parameters","Estimates","SE","%RSE","Lower","Upper","%Shrinkage","%CV")
      tot.res[is.na(tot.res)| is.nan(tot.res) | tot.res=="NA"| tot.res=="NaN"]<-" "
      tot.res<-tot.res[-which(apply(tot.res,1,function(x) sum(x==" "))==7),]
      write.csv(tot.res,paste(dir.name,"\\",file.id,".sum",sep=""),quote=F)
      write.csv(tot.res,paste(dir.name,"\\",file.id,".sum.csv",sep=""),quote=F)      
   }
   
   runNONMEM.addRUNtable<-function(runID,file.ctl,data.file,Description,parent)
   {  current.dir<-getwd()
      print(data.file)
      new.dir<-paste(current.dir,runID,sep="/")
      dir.create(new.dir,showWarnings=F)
      file.copy(file.ctl,new.dir,overwrite=T)
      file.copy(data.file,new.dir,overwrite=T)
      file.copy(data.file,paste(new.dir,"/",runID,".csv",sep=""),overwrite=T)
      setwd(new.dir)
      NM.command<-paste(Default.NMpath, paste(runID,".ctl",sep=""),paste(runID,".res",sep=""))
      write.table(" ",paste(runID,".res",sep=""))
      system(NM.command,invisible=F,show.output.on.console=F)
      temp.LST<-readLines(paste(runID,".res",sep=""))
      AA<<-temp.LST
      data.n<- select.N(read.csv(paste(runID,".csv",sep="")))
      TOT.RUN$num<-TOT.RUN$num+1
      TOT.RUN$data<-rbind(TOT.RUN$data,c(runID,getwd(),parent))
      TOT.RUN<<-TOT.RUN
      addRunTable(temp.LST,Description,runID,dir.name,parent,data.n,runID)
   }
  
   select.datafile.name<-function(D.CTL)
   { data.id<-grep("\\$DATA",D.CTL)
     temp<-strsplit(D.CTL[data.id],split=" ")
     data.file<-temp[[1]][lapply(temp,function(x) grep(".csv",x))[[1]]]
     return(data.file)
   }  

   
   select.N<-function(data.d)
   { if(sum(colnames(data.d)=="MDV"))
     { data.d<-data.d[data.d$MDV==0,]
     }
     return(nrow(data.d))
   }
     
   addRunTable<-function(D.LST,Description,file.id,dir.name,parent,data.n,runID)
   {  n.lst<-length(D.LST)
      Date<-D.LST[n.lst-1]
      Time<-D.LST[n.lst]

      D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      min.id<- which(indicator=="0MINIMIZATION")
      min.id<-min.id[length(min.id)]
      Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]
      param.num<- select.n.param(D.LST)
      data.n<-select.N(read.csv(paste(file.id,".csv",sep="")))
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])
      cond.id<-grep(" EIGENVALUES ",D.LST)
      if(length(cond.id)!=0)
      {  
         flag<-T
         id.current<-cond.id+5
         cond.line<-0
         while(flag)
         {  ttemp<-D.LST[id.current]
            flag<-ttemp!=" "
            if(flag)
            {  cond.line<-cond.line+1
               id.current<-id.current+1
            }
         }
         temp<-NULL
         for(i in 1:cond.line)
            temp<-c(temp,strsplit(D.LST[cond.id+5+cond.line+i],split=" ")[[1]])
         temp<-as.numeric(temp[temp!=""&temp!="+"])
         cond.num<-round(max(temp)/min(temp),3)        
      } else
      {  cond.num<-NA
      }   
      
      Obj<- select.OBJ(D.LST)
      AIC<-Obj+2*param.num
      AICc<-round(Obj+2*param.num+2*param.num*(param.num+1)/(data.n-param.num-1),3)
      SBC<-round(Obj+param.num*log(data.n),3)
      COV.THETA<-select.COV(D.LST)    
      if(is.na(COV.THETA))
      {  COV<-"NONE"
         cond.num<-NA
      } else
      {  COV<-"OK"
      }  
      temp<-c(runID,Date,Time,Min,COV,Obj,AIC,AICc,SBC,cond.num,parent,Description,param.num)
      run.table[]<-rbind(run.table[],temp)
#      if(TOT.RUN$num>2)
#      {  run.table[]<-rbind(run.table[],temp)
#      } else
#      {  run.table[][TOT.RUN$num,]<-temp
#      }
      run.table<<-run.table
   }
   changeCTL.id<-function(D.CTL,old.ID,new.ID)
   {  
      oldID<-paste(old.ID,"\\.",sep="")
      newID<-paste(new.ID,"\\.",sep="")
      old.id<-grep(oldID,D.CTL)
      D.CTL[old.id]<-sub(oldID,newID,D.CTL[old.id])
      return(D.CTL)
   }   
   
   addCTL.TH<-function(D.CTL,fix.id,fix.value)
   {  TH.id<-grep("\\$THETA",D.CTL)
      OMEGA.id<-grep("\\$OMEGA",D.CTL)
      TH.n<-length(fix.id)
      D.CTL<-D.CTL[-(which(unlist(lapply(D.CTL[TH.id:(OMEGA.id-1)],
                  function(x) sum(strsplit(x,split=" ")[[1]]!="")))==0)+TH.id-1)]
      D.CTL[TH.id+fix.id]<-paste("   ",fix.value," FIX",sep=" ")
      return(D.CTL)
   }
   
   select.n.param<-function(D.LST)
   { EST<-select.EST(D.LST)
     ctl.start.id<-grep("\\$PROB",D.LST)
     ctl.end.id<-grep("\\$TABLE",D.LST)[1]
     n.fix<-length(grep("FIX",D.LST[ctl.start.id:ctl.end.id]))
     n.param<-length(EST$TH)+sum(EST$OMEGA!=0)+length(EST$SIGMA)-n.fix
     return(n.param)
   }  
   
   select.COV<-function(D.LST)
   {  cov.start.id<-grep("     COVARIANCE MATRIX OF ESTIMATE  ",D.LST)
      if(length(cov.start.id)!=0)
      {  cov.start.id<-cov.start.id[length(cov.start.id)]
         cov.last.id<-grep("      CORRELATION MATRIX OF ESTIMATE",D.LST)
         cov.last.id<-cov.last.id[length(cov.last.id)]
     
         Result.LST<-D.LST[(cov.start.id+4):(cov.last.id-5)]
         temp<-unlist(lapply(Result.LST,function(x) strsplit(x,split="   ")[[1]][1]))
         list1<-which(temp=="+")
         list2<-which(temp==" ")
         COV.name.temp<-unlist(strsplit(Result.LST[1:list2[1]],split="  "))
         COV.name.temp<-COV.name.temp[COV.name.temp!=""&COV.name.temp!=" "]
         n.COV<-length(COV.name.temp)

         sel.id<-NULL
         for(i in 1:length(list1))
         {  if(list2[i]<list1[i])
               list2<-list2[-i]
            sel.id<-c(sel.id,list1[i]:list2[i])        
         }
         COV.temp<-unlist(strsplit(Result.LST[sel.id],split=" "))
         COV.temp[COV.temp== "........."]<-"0"
         COV<-as.numeric(COV.temp[COV.temp!="+" & COV.temp!=""])
         COV.m<-matrix(0,nrow=n.COV,ncol=n.COV)
         k<-1
         for(i in 1:n.COV)
         {  for(j in 1:i)
            {  COV.m[i,j]<-COV[k]
               COV.m[j,i]<-COV[k]
               k<-k+1
            }
         }
         colnames(COV.m)<-COV.name.temp
         rownames(COV.m)<-COV.name.temp
      } else
      { COV.m<-NA
      }   
      return(COV.m)
   }
   
   select.OBJ<-function(D.LST)
   {     
      OBJ.id<-grep("#OBJV",D.LST)
      OBJ.id<-OBJ.id[length(OBJ.id)]
      temp<-strsplit(D.LST[OBJ.id],split="  ")[[1]]
      temp<-temp[-c(1,length(temp))]
      temp<-temp[temp!=""]
      OBJ.value<-as.numeric(temp)
      return(OBJ.value)  
   }

   select.SE<-function(D.LST)
   {  final.start.id<-grep("STANDARD ERROR OF ESTIMATE ",D.LST)
      final.start.id<-final.start.id[length(final.start.id)]
      Result.LST<-D.LST[final.start.id:length(D.LST)]
      TH.id<-grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS",Result.LST)[1]
      OMEGA.id<-grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS",Result.LST)[1]
      SIGMA.id<-grep("SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS ",Result.LST)[1]
      THname.id<-(which(Result.LST[TH.id:(TH.id+10)]==" ")[1]+TH.id-1)
      OMEGAname.id<-(which(Result.LST[OMEGA.id:(OMEGA.id+10)]==" ")[1]+OMEGA.id-1)
      SIGMAname.id<-(which(Result.LST[SIGMA.id:(SIGMA.id+10)]==" ")[1]+SIGMA.id-1)
      
      TH.name.temp<-unlist(strsplit(Result.LST[(TH.id+1):THname.id],split="  "))
      TH.name.temp<-TH.name.temp[TH.name.temp!=""&TH.name.temp!=" "]
      n.TH<-length(TH.name.temp)
      OMEGA.name.temp<-unlist(strsplit(Result.LST[(OMEGA.id+1):OMEGAname.id],split="  "))
      OMEGA.name.temp<-OMEGA.name.temp[OMEGA.name.temp!=""&OMEGA.name.temp!=" "]
      n.OMEGA<-length(OMEGA.name.temp)
      SIGMA.name.temp<-unlist(strsplit(Result.LST[(SIGMA.id+1):SIGMAname.id],split="  "))
      SIGMA.name.temp<-SIGMA.name.temp[SIGMA.name.temp!=""&SIGMA.name.temp!=" "]
      n.SIGMA<-length(SIGMA.name.temp)
      
      TH.temp<-unlist(strsplit(Result.LST[THname.id:(OMEGA.id-1)],split=" "))
      TH.temp[TH.temp== "........."]<-"0"      
      TH<-as.numeric(TH.temp[TH.temp!=""])
      names(TH)<-TH.name.temp

      OMEGA.temp<-strsplit(Result.LST[OMEGAname.id:(SIGMA.id-1)],split=" ")
      OMEGA.temp<-unlist(OMEGA.temp[unlist(lapply(OMEGA.temp,function(x){x[1]=="+"}))])
      OMEGA.temp[OMEGA.temp== "........."]<-"0"
      OMEGA<-as.numeric(OMEGA.temp[OMEGA.temp!="+" & OMEGA.temp!=""])
      OMEGA.m<-matrix(0,nrow=n.OMEGA,ncol=n.OMEGA)
      k<-1
      for(i in 1:n.OMEGA)
      {  for(j in 1:i)
         {  OMEGA.m[i,j]<-OMEGA[k]
            OMEGA.m[j,i]<-OMEGA[k]
            k<-k+1
         }
      }
      colnames(OMEGA.m)<-OMEGA.name.temp
      rownames(OMEGA.m)<-OMEGA.name.temp
      
      SIGMA.temp<-strsplit(Result.LST[SIGMAname.id:(which(Result.LST=="1")[1]-1)],split=" ")
      SIGMA.temp<-unlist(SIGMA.temp[unlist(lapply(SIGMA.temp,function(x){x[1]=="+"}))])
      SIGMA.temp[SIGMA.temp== "........."]<-"0"        
      SIGMA<-as.numeric(SIGMA.temp[SIGMA.temp!="+" & SIGMA.temp!=""])
      names(SIGMA)<-SIGMA.name.temp
      SE<-list(TH=TH,OMEGA=OMEGA,OMEGA.m=OMEGA.m,SIGMA=SIGMA)
      return(SE)  
   }  


   select.EST<-function(D.LST)
   {  final.start.id<-grep("FINAL PARAMETER ESTIMATE",D.LST)
      final.start.id<-final.start.id[length(final.start.id)]
      Result.LST<-D.LST[final.start.id:length(D.LST)]
      TH.id<-grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS",Result.LST)[1]
      OMEGA.id<-grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS",Result.LST)[1]
      SIGMA.id<-grep("SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS ",Result.LST)[1]
      THname.id<-(which(Result.LST[TH.id:(TH.id+10)]==" ")[1]+TH.id-1)
      OMEGAname.id<-(which(Result.LST[OMEGA.id:(OMEGA.id+10)]==" ")[1]+OMEGA.id-1)
      SIGMAname.id<-(which(Result.LST[SIGMA.id:(SIGMA.id+10)]==" ")[1]+SIGMA.id-1)
      
      TH.name.temp<-unlist(strsplit(Result.LST[(TH.id+1):THname.id],split="  "))
      TH.name.temp<-TH.name.temp[TH.name.temp!=""&TH.name.temp!=" "]
      n.TH<-length(TH.name.temp)
      if(n.TH<10)
      {   TH.name.temp<-paste("TH",1:n.TH,sep=" ")   
      } else
      {   TH.name.temp<-c(paste("TH",1:9,sep=" "),paste("TH",10:n.TH,sep=""))      
      }   
      OMEGA.name.temp<-unlist(strsplit(Result.LST[(OMEGA.id+1):OMEGAname.id],split="  "))
      OMEGA.name.temp<-OMEGA.name.temp[OMEGA.name.temp!=""&OMEGA.name.temp!=" "]
      n.OMEGA<-length(OMEGA.name.temp)
      SIGMA.name.temp<-unlist(strsplit(Result.LST[(SIGMA.id+1):SIGMAname.id],split="  "))
      SIGMA.name.temp<-SIGMA.name.temp[SIGMA.name.temp!=""&SIGMA.name.temp!=" "]
      n.SIGMA<-length(SIGMA.name.temp)
      
      TH.temp<-unlist(strsplit(Result.LST[THname.id:(OMEGA.id-1)],split=" "))
      TH<-as.numeric(TH.temp[TH.temp!=""])
      names(TH)<-TH.name.temp

      OMEGA.temp<-strsplit(Result.LST[OMEGAname.id:(SIGMA.id-1)],split=" ")
      OMEGA.temp<-unlist(OMEGA.temp[unlist(lapply(OMEGA.temp,function(x){x[1]=="+"}))])
      OMEGA<-as.numeric(OMEGA.temp[OMEGA.temp!="+" & OMEGA.temp!=""])
      OMEGA.m<-matrix(0,nrow=n.OMEGA,ncol=n.OMEGA)
      k<-1
      for(i in 1:n.OMEGA)
      {  for(j in 1:i)
         {  OMEGA.m[i,j]<-OMEGA[k]
            OMEGA.m[j,i]<-OMEGA[k]
            k<-k+1
         }
      }
      colnames(OMEGA.m)<-OMEGA.name.temp
      rownames(OMEGA.m)<-OMEGA.name.temp
      
      SIGMA.temp<-strsplit(Result.LST[SIGMAname.id:(which(Result.LST=="1")[1]-1)],split=" ")
      SIGMA.temp<-unlist(SIGMA.temp[unlist(lapply(SIGMA.temp,function(x){x[1]=="+"}))])
      SIGMA<-as.numeric(SIGMA.temp[SIGMA.temp!="+" & SIGMA.temp!=""])
      names(SIGMA)<-SIGMA.name.temp
      Est<-list(TH=TH,OMEGA=OMEGA,OMEGA.m=OMEGA.m,SIGMA=SIGMA)
      return(Est)  
   }  
   
   select.AddInfo<-function(file.id)
   {  n.lst<-length(D.LST)
      Date<-D.LST[n.lst-1]
      Time<-D.LST[n.lst]

      D.temp<-matrix(D.LST)
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=" ")[[1]][1])
      min.id<- which(indicator=="0MINIMIZATION")
      min.id<-min.id[length(min.id)]
      Min<-strsplit(D.LST[min.id],split=" ")[[1]][2]
      param.num<- select.n.param(D.LST)
      data.n<-select.N(read.csv(paste(file.id,".csv",sep="")))
      indicator<-apply(D.temp,1,function(x) strsplit(x,split=":")[[1]][1])
      cond.id<-grep(" EIGENVALUES ",D.LST)
      if(length(cond.id)!=0)
      {  
         flag<-T
         id.current<-cond.id+5
         cond.line<-0
         while(flag)
         {  ttemp<-D.LST[id.current]
            flag<-ttemp!=" "
            if(flag)
            {  cond.line<-cond.line+1
               id.current<-id.current+1
            }
         }
         temp<-NULL
         for(i in 1:cond.line)
            temp<-c(temp,strsplit(D.LST[cond.id+5+cond.line+i],split=" ")[[1]])
         temp<-as.numeric(temp[temp!=""&temp!="+"])
         cond.num<-round(max(temp)/min(temp),3)        
      } else
      {  cond.num<-NA
      }   
      
      Obj<- select.OBJ(D.LST)
      AIC<-Obj+2*param.num
      AICc<-round(Obj+2*param.num+2*param.num*(param.num+1)/(data.n-param.num-1),3)
      SBC<-round(Obj+param.num*log(data.n),3)
      COV.THETA<-select.COV(D.LST)    
      if(sum(is.na(COV.THETA))!=0)
      {  COV<-"NONE"
         cond.num<-NA
      } else
      {  COV<-"OK"
      }  
      temp<-c(Date,Time,Min,COV,Obj,AIC,AICc,SBC,cond.num,param.num)
      return(temp)   
   }
      
   ############################################################################# 
   # Program information
   #############################################################################
   PI<-function(h,...)
   {  gmessage("*** Program information ***
                \n\"fit4NM\" stands for \"Fit for NONMEM\".
                \nEun-Kyung Lee, Ph.D. (lee.eunk@ewha.ac.kr)
                \n         Assistant Professor
                \n         Department of Statistics
                \n         Ewha Womans University
                \n         Seoul, Korea.
                \nGyujeong Noh, M.D. & Ph.D. (gyujeong.noh@gmail.com)
                \n         Professor
                \n         Department of Clinical Pharmacology and Therapeutics
                \n         Department of Anesthesiology and Pain Medicine
                \n         Asan Medical Center
                \n         University of Ulsan College of Medicine
                \n         Seoul, Korea.
                \nVersion: 3.3.3 (Mar. 16, 2011)
                \nfit4NM 3.3.3 was developed in R 2.10.0 and tested in up to R 2.12.0.",cont=TRUE,width=600)   
   }

   ############################################################################# 
   # Exit
   #############################################################################

   ############################################################################# 
   # Main GUI
   #############################################################################
   TOT.temp<-list()
   TOT.temp$num<-0
   TOT.temp$data<-c(NULL,NULL,NULL)
   
   TOT.RUN<<-TOT.temp
   TOT.RESULT<<-list()

   NONMEM.win<<-gwindow("GUI for NONMEM",width=850,height=300)
   menu.list<-list(Configuration=
                       list('Notes before configuration'=list(handler=ConfigNotes),
                            'Set NONMEM path-default'=list(handler=NMpath1),
                            'Set NONMEM help-default'=list(handler=NMpathHelp1),
                            'Set NONMEM path-alternative'=list(handler=NMpath2),
                            'Set NONMEM help-alternative'=list(handler=NMpathHelp2),
                            'Set external editor'=list(handler=Editorpath),                            
                            'Save configuration(c:/fit4NM/.Rdata)'=list(handler=saveConfig)
                           ),
                   Data=
                       list(
                            'Data manipulation'=
                                 list('Calculate elapsed time'=list(handler=CalcTime),
                                      'Join data'=list(handler=DataJoinhandler),
                                      'Split data'=list(handler=DataSplit),                                     
                                      'Convert data'=list('Column to line'=list(handler=ColtoLine),
                                                          'Line to column'=list(handler=LinetoCol)),
                                      'Biosignal data preparation'
                                              =list('Notes before biosignal data preparation'=list(handler=BiosignalNotes),
                                                    'Selection'
                                                       =list('Batch process'=list(handler=BDS.batch),
                                                             'Individual process'=list(handler=BDS.indiv)),                                                         
                                                    'Central moving average'       
                                                       =list('Batch process'=list(handler=BDS.smooth.batch),
                                                             'Individual process'=list(handler=BDS.smooth)))                                                                                                                                                                                                                               
                                     ),
                            'NONMEM data'=
                                list('Notes before NONMEM data creation'=list(handler=dataNotes),
                                     'Create NONMEM data'=list(handler=DataPrep),
                                     'Create successive ID'=list(handler=data.ID)),
                            'Explore NONMEM data'=
                                 list('Select data file'=list(handler=OpenEDAData),
                                      'Summary'=list('Summary statistics-continuous'=list(handler=Summary.stat), 
                                                     'Summary statistics-continuous by ID'=list(handler=Summary.stat.ind),           
                                                     'Summary statistics-categorical-single level per person'=list(handler=Summary.cat),     
                                                     'Summary statistics-categorical-multiple levels per person'=list(handler=Summary.cat1),
                                                     'Summary statistics-categorical by ID'=list(handler=Summary.stat.cat.ind)),                                                     
                                      'Plot'=list('XY plot'=list(handler=XY.plot),
                                                  'DV vs TIME by ID'=list(handler=ID.plot),
                                                  'DV vs TIME by covariates'=list(handler=IDCOV.plot),
                                                  'Covariate vs covariate'=list(handler=COVvsCOV.plot))
                                     )
                           ),
                   'Control stream'=
                       list('Edit with default editor'=list(handler=EditEditor),
                            'Edit with external editor'=list(handler=ExternalEditor)
                           ),
                   'NONMEM run'=
                       list(                                    
                            'Run table'=list('Make run table from runnumber subfolders'=list(handler=AddRunTable),
                                             'Save run table as csv'=list(handler=saveRUNTABLE.handler),
                                             'Load run table'=list(handler=loadRUNTABLE.handler)),   
                            'Model tree'=list(handler=Tree.handler),
                            'Run'=list('Notes before run'=list(handler=BeforeRun),
                                        'From default editor'=list(handler=Editor),
                                        'From external editor'=list(handler=ExternalRun),
                                        'Direct run'=list(handler=DirectRun)),                                                                  
                            'Explore output'=
                                 list('Select output data' = list(handler=outputselect),
                                      'View run summary' = list(handler=showRES), 
                                      'Measures of performance' = list('Performance errors'=list(handler=Measure.Performance1),
                                                                       'Prediction probability'=list(handler=Measure.Performance2)),                                
                                      'Plot'=list('XY plot'=list(handler=postXY.plot),                                    
                                                  'PRED vs DV'=list(handler=DVvsPRED.plot),
                                                  'RES vs DV'=list(handler=DVvsRES.plot), 
                                                  'RES vs TIME'=list(handler=TIMEvsRES.plot),                  
                                                  'Predictions and DV vs TIME'=list(handler=TIMEvsDVandPRED.plot),
                                                  'Predictions and DV vs TIME by ID'=list(handler=TIMEvsDVandPREDID.plot),                                                 
                                                  'Covariate vs parameter'=list(handler=EBEvsCOV.plot)),
                                      'Covariate search'=list('Wald approximation method'=list(handler=WaldApprox),                                  
                                                  'Randomization test'=list(handler=RandomTest))                                                
                                      ),                                                 
                            'Explore simulated data'=list(handler=simulationD), 
                            'Convert PK parameters'=list(handler=PKparam.converter),                                                                                                                                                           
                            'Open R terminal for xpose'=list(handler=OpenXpose)
                           ),
                   'Model evaluations'=
                       list('Notes before model evaluations'=list(handler=VPCNote),                                                  
                            'Predictive checks'=list(handler=VPC1GUI),
                            'Bootstrap'=list(handler=Boot),                            
                            'Case deletion / Jackkinfe diagnostics'=  list(handler=CDD),   
                            'Log-likelihood profiling'=list(handler=LLprofiling),                              
                            'Cross validation'=  list(handler=CV),                                                                                     
                            'External validation of TCI'=  list(handler=Performance.CCIP)                         
                           ),   
                                               
                   'NONMEM help'=list('Open default help'=list(handler=Help1),
                                      'Open alternative help'=list(handler=Help2)),    
                   'Program information'=  list(handler=PI),                           
                   'Exit'=  list(handler=function(h,...) dispose(NONMEM.win))
                  ) 
   gmenu(menu.list,cont=NONMEM.win)
   table.name<-c("Run number","Date","Time","MIN","COV","OFV","AIC","AICc","SBC","Condtion number",
                "Parents","Model description","# of parameters") 
   nonmem.run<<-matrix("",ncol=length(table.name),nrow=2)
   colnames(nonmem.run)<-table.name
   run.table<<-gtable(nonmem.run,chosencol=length(table.name),cont=NONMEM.win)
}

       