#' @export
#' @import tcltk
#' @import tcltk2
#' @importFrom tm VCorpus VectorSource tm_map stripWhitespace removePunctuation DocumentTermMatrix
#' @importFrom wordcloud wordcloud
#' @importFrom RColorBrewer brewer.pal
#' @importFrom arules apriori
#' @import arulesViz
#' @importFrom ggplot2 ggplot aes ggtitle
#' @importFrom gridExtra tableGrob grid.arrange
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics abline barplot pie plot plot.new text title
#' @importFrom stats dist median na.omit quantile
#' @importFrom utils data read.table

telp=function()
{
  
  # Creating loading screen to run the package
  ttLoad=tktoplevel()
  tkwm.title(ttLoad,"TELP Package - Loading Package")
  labelText1Load=tclVar("Loading package... Please wait!")
  label1Load=tklabel(ttLoad,text=tclvalue(labelText1Load),width=25)
  tkconfigure(label1Load,textvariable=labelText1Load)
  labelText2Load=tclVar("Please press the 'Start' button to proceed.")
  label2Load=tklabel(ttLoad,text=tclvalue(labelText2Load),width=35)
  tkconfigure(label2Load,textvariable=labelText2Load)
  tkfocus(ttLoad)
  
  # Continuing with the loading of the package
  OnCarregar=function()
  {
    tkdestroy(ttLoad)
    # Loading standardized Database (TXT format)
    # Column 1: Identification code;  Column 2: First word;  Column 3: Second word;  Column 4: Third word;  Column 5: Fourth word;  Column 6: Fifth word;  Column 7: Meaning of the first word;  Other columns: Stratification Variables
    tt=tktoplevel()
    tkwm.title(tt,"TELP Package - The Free Evocation of Words Technique")
    labelText=tclVar("Choose the Database")
    label1=tklabel(tt,text=tclvalue(labelText),width=20)
    tkconfigure(label1,textvariable=labelText)
    nome=tclVar("")
    nomeaux=tclVar(0)
    tkfocus(tt)
    tt2=""
    tt3=""
    statistics=""
    dadosorig=""
    
    # Checking the status of the currently Database
    OnVerificar=function()
    {
      link=file.choose() 
      tclvalue(nome)=link
    }
    
    ver_but=tkbutton(tt,text="...",command=OnVerificar)
    tkgrid(label1,ver_but)
    
    # Performing Database reading
    OnLer=function()
    {
      if(tclvalue(nome)!="")
      {
        tclvalue(nomeaux)=1 
        tk2state.set(ver_but,state="disabled")
        tk2state.set(ler_but,state="disabled")
        tk2state.set(exemplo_but,state="disabled")
      }
      else
      {
        tkmessageBox(message=paste("No file selected."),icon="warning",type="ok",default="ok",title="Message")
      }
    }
    OnExemplificar=function()
    {
      if(exists(data('statistics',package='TELP',envir=environment()))==TRUE)
      {
        tclvalue(nomeaux)=3 
        tk2state.set(ver_but,state="disabled")
        tk2state.set(ler_but,state="disabled")
        tk2state.set(exemplo_but,state="disabled")
      }
      else
      {
        tkmessageBox(message=paste("An error occurred while reading the example database, please try again."),icon="warning",type="ok",default="ok",title="Message")
      }
    }
    
    ler_but=tkbutton(tt,text="Read file",command=OnLer)
    exemplo_but=tkbutton(tt,text="Use Example",command=OnExemplificar)    
    tkgrid(exemplo_but,ler_but,tklabel(tt,text=" "),tklabel(tt,text=" "))  
    tkwait.variable(nomeaux)
    
    if(tclvalue(nomeaux)==1)
    {
      dados=read.table(tclvalue(nome),header=T,sep="\t")
      rm(statistics)
    }
    if(tclvalue(nomeaux)==3)
    {
      data('statistics',package='TELP',envir=environment())
      dados=statistics
      rm(statistics)
    }
    
    # Checking for Database stratification
    labelTextest=tclVar("Stratification Variables")
    label1est=tklabel(tt,text=tclvalue(labelTextest),width=20)
    tkconfigure(label1est,textvariable=labelTextest)
    opcoes=c("General",if(ncol(dados)==8){colnames(dados[,7:8])[2]}else{colnames(dados[,8:ncol(dados)])})
    opcoes_estrat=tk2combobox(tt,values=opcoes)
    tkgrid(label1est,opcoes_estrat)
    estrat=tclVar("General")
    tkconfigure(opcoes_estrat,textvariable=estrat)
    nomeaux2=tclVar(0)
    
    # Graphic of sectors for variable that will be used in the Stratification
    OnGrafSetor=function()
    {
      if(tclvalue(estrat)=="General")
      {
        tkmessageBox(message=paste("You can not perform this function when there is no stratification selected."),icon="warning",type="ok",default="ok",title="Message")
      }
      else
      {
        if(tclvalue(nomeaux2)==2)
        {
          pie((table(as.factor(dadosorig[,(which(opcoes==tclvalue(estrat))+6)]))/sum(table(as.factor(dadosorig[,(which(opcoes==tclvalue(estrat))+6)]))))*100,main="Graphic of sectors for stratification variable")
        }
        else
        {
          pie((table(as.factor(dados[,(which(opcoes==tclvalue(estrat))+6)]))/sum(table(as.factor(dados[,(which(opcoes==tclvalue(estrat))+6)]))))*100,main="Graphic of sectors for stratification variable")
        }        
      }
    }
    
    # Stratifying the database using this variable 
    cats_new=tclVar("")
    OnEstratificar=function()
    { 
      if(tclvalue(estrat)!="General")
      {
        # Defining variable level for realization of stratification
        OnDefinirNivel=function()
        {
          tclvalue(cats_new)=tclvalue(cats)
          tk2state.set(opcoes_cats,state="disabled")
          tk2state.set(niv_but,state="disabled")            
          tclvalue(nomeaux2)=2 
        }
        tk2state.set(opcoes_estrat,state="disabled")
        tk2state.set(est_but,state="disabled")
        labelText2=tclVar("Level of stratification")
        label2=tklabel(tt,text=tclvalue(labelText2),width=20)
        tkconfigure(label2,textvariable=labelText2)
        opcoes2=levels(as.factor(dados[,(which(opcoes==tclvalue(estrat))+6)]))        
        opcoes_cats=tk2combobox(tt,values=opcoes2)
        tkgrid(label2,opcoes_cats)
        cats=tclVar(levels(as.factor(dados[,(which(opcoes==tclvalue(estrat))+6)]))[1])
        tkconfigure(opcoes_cats,textvariable=cats)
        tclvalue(cats_new)=tclvalue(cats)
        niv_but=tkbutton(tt,text="Choose level",command=OnDefinirNivel)
        tkgrid(niv_but)		
      }
      else
      {		
        tk2state.set(opcoes_estrat,state="disabled")
        tk2state.set(est_but,state="disabled")
        tk2state.set(grafset_but,state="disabled")
        tclvalue(nomeaux2)=1 
      }	
    }
    
    grafset_but=tkbutton(tt,text="Graphic of sectors",command=OnGrafSetor)
    est_but=tkbutton(tt,text="Stratify",command=OnEstratificar)
    tkgrid(grafset_but,est_but)
    tkwait.variable(nomeaux2)
    
    if(tclvalue(nomeaux2)==2)
    {
      dadosorig=dados
      dados=dados[as.factor(dados[,(which(opcoes==tclvalue(estrat))+6)])==tclvalue(cats_new),]
    }
    
    # Performing Outline of Analysis of The Free Evoction of Words Technique
    
    # Transforming Data-Table in Matrix
    matriz=as.matrix(dados[,c(2,3,4,5,6)])
    
    # Calculating frequency of words
    freq=table(matriz)
    
    # Storing words evoked
    palavras=rownames(freq)
    
    # Definition of OME levels
    ome=rep(0,length(palavras))
    for(i in 1:length(palavras))
    {
      ome[i]=round((sum(na.omit(as.integer(dados[,2]==palavras[i])))+2*sum(na.omit(as.integer(dados[,3]==palavras[i])))+3*sum(na.omit(as.integer(dados[,4]==palavras[i])))+4*sum(na.omit(as.integer(dados[,5]==palavras[i])))+5*sum(na.omit(as.integer(dados[,6]==palavras[i]))))/(freq[[i]]),digits=2)
    }
    
    # Engaging the words with their frequency and level of OME
    resini=cbind(palavras,freq,ome)
    rownames(resini)=rep("",nrow(resini))
    
    # Bar Chart of frequency for the words evoked position
    OnGrafBarra=function()
    {
      tabela_freq=rep(0,6)
      for(i in 1:nrow(dados))
      {
        if(is.na(dados[i,6]))
        {
          if(is.na(dados[i,5]))
          {
            if(is.na(dados[i,4]))
            {
              if(is.na(dados[i,3]))
              {
                if(is.na(dados[i,2]))
                {
                  tabela_freq[1]=tabela_freq[1]+1
                }
                else
                {
                  tabela_freq[2]=tabela_freq[2]+1
                }
              }
              else
              {
                tabela_freq[3]=tabela_freq[3]+1
              }
            }
            else
            {
              tabela_freq[4]=tabela_freq[4]+1
            }
          }
          else
          {
            tabela_freq[5]=tabela_freq[5]+1
          }
        }
        else
        {
          tabela_freq[6]=tabela_freq[6]+1
        }
      }
      barplot(tabela_freq,names.arg=c("0","1","2","3","4","5"),main="Frequency for the words evoked position")
    }
    
    # Outlining the wordsCloud by the frequency
    OnNuvemPalavra=function()
    {
      corpus=c(as.character(dados[,2]),as.character(dados[,3]),as.character(dados[,4]),as.character(dados[,5]),as.character(dados[,6]))
      corpus1=VCorpus(VectorSource(corpus))
      corpus2=tm_map(corpus1,stripWhitespace)
      corpus2=tm_map(corpus2,removePunctuation)
      wordcloud(corpus2,main="title",min.freq=0,random.order=F,colors=brewer.pal(5,"Dark2"))
    }
    
    grafbar_but=tkbutton(tt,text="Graphic of frequency for the words evoked position",command=OnGrafBarra)
    nuvem_but=tkbutton(tt,text="WordsCloud",command=OnNuvemPalavra)
    tkgrid(tklabel(tt,text=" "),tklabel(tt,text=" "),tklabel(tt,text=" "))
    tkgrid(grafbar_but,tklabel(tt,text=" "),nuvem_but)
    tkgrid(tklabel(tt,text=" "),tklabel(tt,text=" "),tklabel(tt,text=" "))
    labelText3=tclVar("    Cutoff frequency level:")
    label3=tklabel(tt,text=tclvalue(labelText3))
    tkconfigure(label3,textvariable=labelText3)
    lbl_cortefreq=tclVar("")
    resp_cortefreq=tk2spinbox(tt,from=0,to=max(freq),increment=1,textvariable=lbl_cortefreq)
    tkgrid(label3,resp_cortefreq,tklabel(tt,text="(Default value:)"),tklabel(tt,text=round(median(freq),digits=0)))
    labelText4=tclVar("    Cutoff OME level:")
    label4=tklabel(tt,text=tclvalue(labelText4))
    tkconfigure(label4,textvariable=labelText4)
    lbl_corteome=tclVar("")
    resp_corteome=tk2spinbox(tt,from=0,to=max(ome),increment=0.01,textvariable=lbl_corteome)
    tkgrid(label4,resp_corteome,tklabel(tt,text="(Default value:)"),tklabel(tt,text=round(median(ome),digits=2)))
    
    # Drawing the Chart Quadrants
    OnGrafQuad=function()
    {
      if((tclvalue(lbl_cortefreq)=="")||(as.integer(tclvalue(lbl_cortefreq))<=0)||(is.na(as.integer(tclvalue(lbl_cortefreq)))=="TRUE"))
      {
        cortefreq=round(median(as.integer(resini[,2])),digits=0)
      }
      else
      {
        cortefreq=as.integer(tclvalue(lbl_cortefreq))
      }
      if((tclvalue(lbl_corteome)=="")||(as.double(tclvalue(lbl_corteome))<=0)||(is.na(as.double(tclvalue(lbl_corteome)))=="TRUE"))
      {
        corteome=round(median(as.double(resini[,3])),digits=2)
      }
      else
      {
        corteome=as.double(tclvalue(lbl_corteome))
      }
      plot(as.integer(resini[,2]),as.double(resini[,3]),xlab="Word frequency",ylab="Average order of evocation (OME)",col=rainbow(nrow(resini)),pch=19)
      abline(v=cortefreq,h=corteome)
      text(x=quantile(as.integer(resini[,2]),0.75),y=quantile(as.double(resini[,3]),0.25),col="darkgray","1th Quadrant")
      text(x=quantile(as.integer(resini[,2]),0.75),y=quantile(as.double(resini[,3]),0.75),col="darkgray","2nd Quadrant")
      text(x=quantile(as.integer(resini[,2]),0.25),y=quantile(as.double(resini[,3]),0.25),col="darkgray","3rd Quadrant")
      text(x=quantile(as.integer(resini[,2]),0.25),y=quantile(as.double(resini[,3]),0.75),col="darkgray","4th Quadrant")
    }
    
    # Performing The Free Evocation of Words Technique
    nomeaux3=tclVar(0) 
    centro=tclArray()
    OnDetQuad=function()
    {
      tk2column(saida1,"delete","s1pal")
      tk2column(saida1,"delete","s1freq")
      tk2column(saida1,"delete","s1ome")
      tk2column(saida1,"add","s1pal",label="Word",width=35)
      tk2column(saida1,"add","s1freq",label="Frequency",width=10)
      tk2column(saida1,"add","s1ome",label="OME",width=10)
      tk2column(saida2,"delete","s2pal")
      tk2column(saida2,"delete","s2freq")
      tk2column(saida2,"delete","s2ome")
      tk2column(saida2,"add","s2pal",label="Word",width=35)
      tk2column(saida2,"add","s2freq",label="Frequency",width=10)
      tk2column(saida2,"add","s2ome",label="OME",width=10)
      tk2column(saida3,"delete","s3pal")
      tk2column(saida3,"delete","s3freq")
      tk2column(saida3,"delete","s3ome")
      tk2column(saida3,"add","s3pal",label="Word",width=35)
      tk2column(saida3,"add","s3freq",label="Frequency",width=10)
      tk2column(saida3,"add","s3ome",label="OME",width=10)
      tk2column(saida4,"delete","s4pal")
      tk2column(saida4,"delete","s4freq")
      tk2column(saida4,"delete","s4ome")
      tk2column(saida4,"add","s4pal",label="Word",width=35)
      tk2column(saida4,"add","s4freq",label="Frequency",width=10)
      tk2column(saida4,"add","s4ome",label="OME",width=10)
      if((tclvalue(lbl_cortefreq)=="")||(as.integer(tclvalue(lbl_cortefreq))<=0)||(is.na(as.integer(tclvalue(lbl_cortefreq)))=="TRUE"))
      {
        cortefreq=round(median(as.integer(resini[,2])),digits=0)
      }
      else
      {
        cortefreq=as.integer(tclvalue(lbl_cortefreq))
      }
      if((tclvalue(lbl_corteome)=="")||(as.double(tclvalue(lbl_corteome))<=0)||(is.na(as.double(tclvalue(lbl_corteome)))=="TRUE"))
      {
        corteome=round(median(as.double(resini[,3])),digits=2)
      }
      else
      {
        corteome=as.double(tclvalue(lbl_corteome))
      }
      QUAD1=resini[(as.integer(resini[,2])>=cortefreq&as.double(resini[,3])<corteome),]
      QUAD2=resini[(as.integer(resini[,2])>=cortefreq&as.double(resini[,3])>=corteome),]
      QUAD3=resini[(as.integer(resini[,2])<cortefreq&as.double(resini[,3])<corteome),]
      QUAD4=resini[(as.integer(resini[,2])<cortefreq&as.double(resini[,3])>=corteome),]
      tk2insert.multi(saida1,"end",as.matrix(QUAD1))
      tk2insert.multi(saida2,"end",as.matrix(QUAD2))
      tk2insert.multi(saida3,"end",as.matrix(QUAD3))
      tk2insert.multi(saida4,"end",as.matrix(QUAD4))
      tclvalue(nomeaux3)=1
    }
    
    grafquad_but=tkbutton(tt,text="Chart Quadrants",command=OnGrafQuad)
    quads_but=tkbutton(tt,text="Determination of the Four Quadrants",command=OnDetQuad)
    tkgrid(grafquad_but,tklabel(tt,text=" "),quads_but)
    tkgrid(tklabel(tt,text=" "),tklabel(tt,text="Figure of Four Quadrants:"))
    saida1=tk2mclistbox(tt,width=55,resizablecolumns=FALSE)
    tk2column(saida1,"add","s1pal",label="Word",width=35)
    tk2column(saida1,"add","s1freq",label="Frequency",width=10)
    tk2column(saida1,"add","s1ome",label="OME",width=10)
    saida2=tk2mclistbox(tt,width=55,resizablecolumns=FALSE)
    tk2column(saida2,"add","s2pal",label="Word",width=35)
    tk2column(saida2,"add","s2freq",label="Frequency",width=10)
    tk2column(saida2,"add","s2ome",label="OME",width=10)
    saida3=tk2mclistbox(tt,width=55,resizablecolumns=FALSE)
    tk2column(saida3,"add","s3pal",label="Word",width=35)
    tk2column(saida3,"add","s3freq",label="Frequency",width=10)
    tk2column(saida3,"add","s3ome",label="OME",width=10)
    saida4=tk2mclistbox(tt,width=55,resizablecolumns=FALSE)
    tk2column(saida4,"add","s4pal",label="Word",width=35)
    tk2column(saida4,"add","s4freq",label="Frequency",width=10)
    tk2column(saida4,"add","s4ome",label="OME",width=10)
    tkgrid(saida1,tklabel(tt,text=" "),saida2)
    tkgrid(saida3,tklabel(tt,text=" "),saida4)
    
    # Creating options for page changes
    nomeaux4=tclVar(0)
    nomeaux5=tclVar(0)
    id_sel_new=tclVar("")
    pal_sel_new=tclVar("")
    metod_sel_new=tclVar("")
    assoc_sel_new=tclVar("")
    OnMaisOpcoes=function()
    {
      if(exists("tt2")==TRUE)
      {
        rm(tt2)
      }
      tt2=tktoplevel()
      tkwm.title(tt2,"TELP Package - Additional Options Analysis")
      
      # Outlining the Meaning for First Words
      OnSelecionarPalavra=function()
      {
        OnSignificado=function()
        {
          tclvalue(pal_sel_new)=tclvalue(pal_sel)
          tclvalue(id_sel_new)=tclvalue(id_sel)    
          tk2column(saida5,"delete","s5sig")
          tk2column(saida5,"add","s5sig",label="Meaning",width=115)	
          signi=dados[as.integer(tclvalue(id_sel_new)),7]
          tk2insert.multi(saida5,"end",signi)
          tclvalue(nomeaux4)=1 
        }
        tclvalue(pal_sel_new)=tclvalue(pal_sel)
        tk2state.set(opcoes_pal,state="disabled")
        tk2state.set(pal_but,state="disabled")
        labelText6=tclVar("Choose the ID:")
        label6=tklabel(tt2,text=tclvalue(labelText6),width=20)
        tkconfigure(label6,textvariable=labelText6)
        opcoes4=na.omit(dados[dados[,2]==tclvalue(pal_sel_new),1])
        opcoes4=opcoes4[1:length(opcoes4)]
        opcoes_id=tk2combobox(tt2,values=opcoes4)
        sig_but=tkbutton(tt2,text="Choose the Meaning",command=OnSignificado)
        tkgrid(label6,opcoes_id,sig_but)
        id_sel=tclVar(opcoes4[1])        
        tkconfigure(opcoes_id,textvariable=id_sel)        
        tclvalue(id_sel_new)=tclvalue(id_sel)
        tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))
        tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text="Meaning of the selected word:"))
        saida5=tk2mclistbox(tt2,width=115,resizablecolumns=FALSE)
        tk2column(saida5,"add","s5sig",label="Meaning",width=115)
        tkgrid(tklabel(tt2,text=" "),saida5)
      }
      
      labelText5=tclVar("Choose the word:")
      label5=tklabel(tt2,text=tclvalue(labelText5),width=20)
      tkconfigure(label5,textvariable=labelText5)
      pal_signi=rownames(table(dados[,2]))
      opcoes3=pal_signi
      opcoes_pal=tk2combobox(tt2,values=opcoes3)
      pal_but=tkbutton(tt2,text="Set word",command=OnSelecionarPalavra)
      tkgrid(label5,opcoes_pal,pal_but)
      pal_sel=tclVar(palavras[1])      
      tkconfigure(opcoes_pal,textvariable=pal_sel)
      tclvalue(pal_sel_new)=tclvalue(pal_sel)
      tkwait.variable(nomeaux4)
      
      # Performing Associations Rules with words evoked
      labelText7=tclVar("Choose the word:")
      label7=tklabel(tt2,text=tclvalue(labelText7),width=20)
      tkconfigure(label7,textvariable=labelText7)
      pal_assoc=rownames(table(dados[,2]))
      opcoes5=pal_assoc
      opcoes_assoc=tk2combobox(tt2,values=opcoes5)
      labelText8=tclVar("Choose the method:")
      label8=tklabel(tt2,text=tclvalue(labelText8),width=20)
      tkconfigure(label8,textvariable=labelText8)
      opcoes6=c("Support","Confidence","Lift")
      opcoes_metod=tk2combobox(tt2,values=opcoes6)
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))
      tkgrid(label7,opcoes_assoc,label8,opcoes_metod)
      assoc_sel=tclVar(opcoes5[1])
      tkconfigure(opcoes_assoc,textvariable=assoc_sel)
      tclvalue(assoc_sel_new)=tclvalue(assoc_sel)
      metod_sel=tclVar(opcoes6[1])      
      tkconfigure(opcoes_metod,textvariable=metod_sel)
      tclvalue(metod_sel_new)=tclvalue(metod_sel)
      
      OnAssociar=function()
      {
        tclvalue(assoc_sel_new)=tclvalue(assoc_sel)
        tclvalue(metod_sel_new)=tclvalue(metod_sel)
        metodo=""
        if(tclvalue(metod_sel_new)=="Support")
        {
          metodo="support"
        }
        else if(tclvalue(metod_sel_new)=="Confidence")
        {
          metodo="confidence"
        }
        else
        {
          metodo="lift"
        }
        dadosmin=dados[,2:6]
        colnames(dadosmin)=c("WORD01","WORD02","WORD03","WORD04","WORD05")
        regras=apriori(dadosmin,parameter=list(supp=0.001,conf=0.8,maxlen=3),appearance=list(default="lhs",rhs=paste("WORD01=",tclvalue(assoc_sel_new),sep="")))
        options(digits=2)
        regras=sort(regras,by=metodo,decreasing=TRUE)	
        if(length(regras)>0)
        {
          plot(regras[1:if(length(regras)<5) { length(regras) } else { 5 }],method="graph")
        }
        else
        {
          tkmessageBox(message=paste("This word has no association, please select another."),icon="warning",type="ok",default="ok",title="Message")
        }
        rm(regras)
        tclvalue(assoc_sel_new)=tclvalue(assoc_sel)
        tclvalue(metod_sel_new)=tclvalue(metod_sel)
        tclvalue(nomeaux5)=1
      }
      
      assoc_but=tkbutton(tt2,text="The word associations",command=OnAssociar)
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))
      tkgrid(tklabel(tt2,text=" "),assoc_but)
      
      # Calculation of summary statistics Distances Social Representation for Central Nucleus
      OnMedidas=function()
      {
        tk2state.set(NC_but,state="disabled")
        if(tclvalue(nomeaux3)==1)
        {
          tk2column(saida6,"delete","s6min")
          tk2column(saida6,"delete","s6q1")
          tk2column(saida6,"delete","s6med")
          tk2column(saida6,"delete","s6mean")
          tk2column(saida6,"delete","s6q3")
          tk2column(saida6,"delete","s6max")
          tk2column(saida6,"add","s6min",label="Minimum",width=19)
          tk2column(saida6,"add","s6q1",label="1th Quartile",width=19)
          tk2column(saida6,"add","s6med",label="Median",width=19)
          tk2column(saida6,"add","s6mean",label="Mean",width=20)
          tk2column(saida6,"add","s6q3",label="3rd Quartile",width=19)
          tk2column(saida6,"add","s6max",label="Maximum",width=19)
          if((tclvalue(lbl_cortefreq)=="")||(as.integer(tclvalue(lbl_cortefreq))<=0)||(is.na(as.integer(tclvalue(lbl_cortefreq)))=="TRUE"))
          {
            cortefreq=round(median(as.integer(resini[,2])),digits=0)
          }
          else
          {
            cortefreq=as.integer(tclvalue(lbl_cortefreq))
          }
          if((tclvalue(lbl_corteome)=="")||(as.double(tclvalue(lbl_corteome))<=0)||(is.na(as.double(tclvalue(lbl_corteome)))=="TRUE"))
          {
            corteome=round(median(as.double(resini[,3])),digits=2)
          }
          else
          {
            corteome=as.double(tclvalue(lbl_corteome))
          }
          nc_pal=resini[(as.integer(resini[,2])>=cortefreq&as.double(resini[,3])<corteome),]
          nc_dist=rep(0,nrow(dados))
          for(i in 1:nrow(dados))
          {
            evocs=matrix("",max(length(nc_pal),5),2)
            evocs[,1]=c(nc_pal,rep("",max(length(nc_pal),5)-length(nc_pal)))
            evocs[,2]=c(t(dados[i,2:6]),rep("",max(length(nc_pal),5)-5))  
            colnames(evocs)=c("Central Nucleus","Data")
            evocs=as.data.frame(evocs)
            corp=VCorpus(VectorSource(evocs),readerControl=list(language="en"))
            dtm=DocumentTermMatrix(corp)
            mydtm=as.matrix(dtm)
            distancia=as.matrix(dist(mydtm))
            nc_dist[i]=distancia[2,1]  
          }	
          nc_infdist=round(summary(nc_dist),digits=2)
          tk2insert.multi(saida6,"end",nc_infdist)
        }
        else
        {
          tkmessageBox(message=paste("To perform this function its necessary to determine the Figure of Four Quadrants."),icon="warning",type="ok",default="ok",title="Message")
        }
      }
      
      NC_but=tkbutton(tt2,text="Central Nucleus Measures",command=OnMedidas)
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))
      tkgrid(tklabel(tt2,text=" "),NC_but)
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text="Central Nucleus Measures for the individuals:"))
      saida6=tk2mclistbox(tt2,width=115,resizablecolumns=FALSE)
      tk2column(saida6,"add","s6min",label="Minumum",width=19)
      tk2column(saida6,"add","s6q1",label="1th Quartile",width=19)
      tk2column(saida6,"add","s6med",label="Median",width=19)
      tk2column(saida6,"add","s6mean",label="Mean",width=20)
      tk2column(saida6,"add","s6q3",label="3rd Quartile",width=19)
      tk2column(saida6,"add","s6max",label="Maximum",width=19)
      tkgrid(tklabel(tt2,text=" "),saida6)
      OnVoltar=function()
      {      
        tkdestroy(tt2)
        tkfocus(tt)
      }	
      volta_but=tkbutton(tt2,text="Back",command=OnVoltar)
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))
      tkgrid(tklabel(tt2,text=" "),volta_but,tklabel(tt2,text=" "))
      tkgrid(tklabel(tt2,text=" "),tklabel(tt2,text=" "),tklabel(tt2,text=" "))    
    }  
    
    x=""
    OnExportar=function()
    {
      if(exists("tt3")==TRUE)
      {
        rm(tt3)
      }
      tt3=tktoplevel()
      tkwm.title(tt3,"TELP Package - Exporting Generated Results")
      labelTextExp=tclVar("Choose the Results:")
      labelExp=tklabel(tt3,text=tclvalue(labelTextExp),width=20)
      tkconfigure(labelExp,textvariable=labelTextExp)
      tkgrid(tklabel(tt3,text=" "),labelExp,tklabel(tt3,text=" "))
      tkgrid(tklabel(tt3,text=" "),tklabel(tt3,text=" "),tklabel(tt3,text=" "))
      check1=tkcheckbutton(tt3)
      check1_valor=tclVar("0")
      tkconfigure(check1,variable=check1_valor)
      tkgrid(check1,tklabel(tt3,text="Graphic of Sectors (Stratification Variable)"))
      if(tclvalue(estrat)=="General")
      {
        tk2state.set(check1,state="disabled")
      }
      else
      {
        tk2state.set(check1,state="normal")
      }
      check2=tkcheckbutton(tt3)
      check2_valor=tclVar("0")
      tkconfigure(check2,variable=check2_valor)
      tkgrid(check2,tklabel(tt3,text="Bar Chart of frequency of evocation"))
      check3=tkcheckbutton(tt3)
      check3_valor=tclVar("0")
      tkconfigure(check3,variable=check3_valor)
      tkgrid(check3,tklabel(tt3,text="WordsCloud"))
      check4=tkcheckbutton(tt3)
      check4_valor=tclVar("0")
      tkconfigure(check4,variable=check4_valor)
      tkgrid(check4,tklabel(tt3,text="Chart Quadrants"))
      check5=tkcheckbutton(tt3)
      check5_valor=tclVar("0")
      tkconfigure(check5,variable=check5_valor)
      tkgrid(check5,tklabel(tt3,text="Figure of Four Quadrants"))
      if(tclvalue(nomeaux3)==1)
      {
        tk2state.set(check5,state="normal")
      }
      else
      {
        tk2state.set(check5,state="disabled")
      }
      check6=tkcheckbutton(tt3)
      check6_valor=tclVar("0")
      tkconfigure(check6,variable=check6_valor)
      tkgrid(check6,tklabel(tt3,text="Meaning of the Selected Word"))
      if(tclvalue(nomeaux4)==1)
      {
        tk2state.set(check6,state="normal")
      }
      else
      {
        tk2state.set(check6,state="disabled")
      }
      check7=tkcheckbutton(tt3)
      check7_valor=tclVar("0")
      tkconfigure(check7,variable=check7_valor)
      tkgrid(check7,tklabel(tt3,text="Graphic for the Word Associations"))
      if(tclvalue(nomeaux5)==1)
      {
        tk2state.set(check7,state="normal")
      }
      else
      {
        tk2state.set(check7,state="disabled")
      }
      check8=tkcheckbutton(tt3)
      check8_valor=tclVar("0")
      tkconfigure(check8,variable=check8_valor)
      tkgrid(check8,tklabel(tt3,text="Central Nucleus Distance Measures"))
      if(tclvalue(nomeaux3)==1)
      {
        tk2state.set(check8,state="normal")
      }
      else
      {
        tk2state.set(check8,state="disabled")
      }
      OnGerar=function()
      {
        if(as.character(tclvalue(check1_valor))=="1"||as.character(tclvalue(check2_valor))=="1"||as.character(tclvalue(check3_valor))=="1"||as.character(tclvalue(check4_valor))=="1"||as.character(tclvalue(check5_valor))=="1"||as.character(tclvalue(check6_valor))=="1"||as.character(tclvalue(check7_valor))=="1"||as.character(tclvalue(check8_valor))=="1")
        {
          caminho=getwd()
          caminho=paste(caminho,"/telp_results.pdf",sep="")
          pdf(file=caminho)
          if(as.character(tclvalue(check1_valor))=="1")
          {
            # Including on PDF: Graphic of Sectors (Stratification Variable)
            pie((table(as.factor(dadosorig[,(which(opcoes==tclvalue(estrat))+6)]))/sum(table(as.factor(dadosorig[,(which(opcoes==tclvalue(estrat))+6)]))))*100,main="Graphic of sectors for stratification variable")
          }
          if(as.character(tclvalue(check2_valor))=="1")
          {
            # Including on PDF: Bar Chart of frequency of evocation
            tabela_freq=rep(0,6)
            for(i in 1:nrow(dados))
            {
              if(is.na(dados[i,6]))
              {
                if(is.na(dados[i,5]))
                {
                  if(is.na(dados[i,4]))
                  {
                    if(is.na(dados[i,3]))
                    {
                      if(is.na(dados[i,2]))
                      {
                        tabela_freq[1]=tabela_freq[1]+1
                      }
                      else
                      {
                        tabela_freq[2]=tabela_freq[2]+1
                      }
                    }
                    else
                    {
                      tabela_freq[3]=tabela_freq[3]+1
                    }
                  }
                  else
                  {
                    tabela_freq[4]=tabela_freq[4]+1
                  }
                }
                else
                {
                  tabela_freq[5]=tabela_freq[5]+1
                }
              }
              else
              {
                tabela_freq[6]=tabela_freq[6]+1
              }
            }
            barplot(tabela_freq,names.arg=c("0","1","2","3","4","5"),main="Frequency for the words evoked position")
          }
          if(as.character(tclvalue(check3_valor))=="1")
          {
            # Including on PDF: WordsCloud
            corpus=c(as.character(dados[,2]),as.character(dados[,3]),as.character(dados[,4]),as.character(dados[,5]),as.character(dados[,6]))
            corpus1=VCorpus(VectorSource(corpus))
            corpus2=tm_map(corpus1,stripWhitespace)
            corpus2=tm_map(corpus2,removePunctuation)
            wordcloud(corpus2,main="title",min.freq=0,random.order=F,colors=brewer.pal(5,"Dark2"))
          }
          if(as.character(tclvalue(check4_valor))=="1")
          {
            # Including on PDF: Chart Quadrants
            if((tclvalue(lbl_cortefreq)=="")||(as.integer(tclvalue(lbl_cortefreq))<=0)||(is.na(as.integer(tclvalue(lbl_cortefreq)))=="TRUE"))
            {
              cortefreq=round(median(as.integer(resini[,2])),digits=0)
            }
            else
            {
              cortefreq=as.integer(tclvalue(lbl_cortefreq))
            }
            if((tclvalue(lbl_corteome)=="")||(as.double(tclvalue(lbl_corteome))<=0)||(is.na(as.double(tclvalue(lbl_corteome)))=="TRUE"))
            {
              corteome=round(median(as.double(resini[,3])),digits=2)
            }
            else
            {
              corteome=as.double(tclvalue(lbl_corteome))
            }
            plot(as.integer(resini[,2]),as.double(resini[,3]),xlab="Word frequency",ylab="Average order of evocation (OME)",col=rainbow(nrow(resini)),pch=19)
            abline(v=cortefreq,h=corteome)
            text(x=quantile(as.integer(resini[,2]),0.75),y=quantile(as.double(resini[,3]),0.25),col="darkgray","1th Quadrant")
            text(x=quantile(as.integer(resini[,2]),0.75),y=quantile(as.double(resini[,3]),0.75),col="darkgray","2nd Quadrant")
            text(x=quantile(as.integer(resini[,2]),0.25),y=quantile(as.double(resini[,3]),0.25),col="darkgray","3rd Quadrant")
            text(x=quantile(as.integer(resini[,2]),0.25),y=quantile(as.double(resini[,3]),0.75),col="darkgray","4th Quadrant")
          }
          if(as.character(tclvalue(check5_valor))=="1")
          {
            # Including on PDF: Figure of Four Quadrants
            if((tclvalue(lbl_cortefreq)=="")||(as.integer(tclvalue(lbl_cortefreq))<=0)||(is.na(as.integer(tclvalue(lbl_cortefreq)))=="TRUE"))
            {
              cortefreq=round(median(as.integer(resini[,2])),digits=0)
            }
            else
            {
              cortefreq=as.integer(tclvalue(lbl_cortefreq))
            }
            if((tclvalue(lbl_corteome)=="")||(as.double(tclvalue(lbl_corteome))<=0)||(is.na(as.double(tclvalue(lbl_corteome)))=="TRUE"))
            {
              corteome=round(median(as.double(resini[,3])),digits=2)
            }
            else
            {
              corteome=as.double(tclvalue(lbl_corteome))
            }
            QUAD1=resini[(as.integer(resini[,2])>=cortefreq&as.double(resini[,3])<corteome),]
            QUAD2=resini[(as.integer(resini[,2])>=cortefreq&as.double(resini[,3])>=corteome),]
            QUAD3=resini[(as.integer(resini[,2])<cortefreq&as.double(resini[,3])<corteome),]
            QUAD4=resini[(as.integer(resini[,2])<cortefreq&as.double(resini[,3])>=corteome),]
            plot1_pdf=ggplot(data.frame(x=c(20, 150)),aes(x))
            ggtitle(paste("First Quadrant - Frequency >= ",cortefreq," & OME < ",corteome,sep=""))
            tbl1_pdf=tableGrob(QUAD1,rows=NULL)
            grid.arrange(tbl1_pdf,as.table=TRUE)
            plot2_pdf=ggplot(data.frame(x=c(20, 150)),aes(x))
            ggtitle(paste("Second Quadrant - Frequency >= ",cortefreq," & OME >= ",corteome,sep=""))
            tbl2_pdf=tableGrob(QUAD2,rows=NULL)
            grid.arrange(tbl2_pdf,as.table=TRUE)
            plot3_pdf=ggplot(data.frame(x=c(20, 150)),aes(x))
            ggtitle(paste("Third Quadrant - Frequency < ",cortefreq," & OME < ",corteome,sep=""))
            tbl3_pdf=tableGrob(QUAD3,rows=NULL)
            grid.arrange(tbl3_pdf,as.table=TRUE)
            plot4_pdf=ggplot(data.frame(x=c(20, 150)),aes(x))
            ggtitle(paste("Fourth Quadrant - Frequency < ",cortefreq," & OME >= ",corteome,sep=""))
            tbl4_pdf=tableGrob(QUAD4,rows=NULL)
            grid.arrange(tbl4_pdf,as.table=TRUE)
          }
          if(as.character(tclvalue(check6_valor))=="1")
          {
            # Including on PDF: Meaning of the Selected Word
            signi=dados[as.integer(tclvalue(id_sel_new)),7]
            plot5_pdf=ggplot(data.frame(x=c(20, 150)),aes(x))
            ggtitle(paste("Meaning of the Word: ",tclvalue(pal_sel_new)," from the individual ID: ",tclvalue(id_sel_new),sep=""))
            tbl5_pdf=tableGrob(signi,rows=NULL)
            grid.arrange(tbl5_pdf,as.table=TRUE)
          }
          if(as.character(tclvalue(check7_valor))=="1")
          {
            # Including on PDF: Graphic for the Word Associations
            metodo=""
            if(tclvalue(metod_sel_new)=="Support")
            {
              metodo="support"
            }
            else if(tclvalue(metod_sel_new)=="Confidence")
            {
              metodo="confidence"
            }
            else
            {
              metodo="lift"
            }
            dadosmin=dados[,2:6]
            colnames(dadosmin)=c("WORD01","WORD02","WORD03","WORD04","WORD05")
            regras=apriori(dadosmin,parameter=list(supp=0.001,conf=0.8,maxlen=3),appearance=list(default="lhs",rhs=paste("WORD01=",tclvalue(assoc_sel_new),sep="")))
            options(digits=2)
            regras=sort(regras,by=metodo,decreasing=TRUE)  
            if(length(regras)>0)
            {
              plot(regras[1:if(length(regras)<5) { length(regras) } else { 5 }],method="graph")
            }
            else
            {
              plot.new()
              title(main="Association chart was not generated")
              text(0.5,0.5,"Insufficient number of rules")            
            }
          }
          if(as.character(tclvalue(check8_valor))=="1")
          {
            # Including on PDF: Central Nucleus Distance Measures
            if((tclvalue(lbl_cortefreq)=="")||(as.integer(tclvalue(lbl_cortefreq))<=0)||(is.na(as.integer(tclvalue(lbl_cortefreq)))=="TRUE"))
            {
              cortefreq=round(median(as.integer(resini[,2])),digits=0)
            }
            else
            {
              cortefreq=as.integer(tclvalue(lbl_cortefreq))
            }
            if((tclvalue(lbl_corteome)=="")||(as.double(tclvalue(lbl_corteome))<=0)||(is.na(as.double(tclvalue(lbl_corteome)))=="TRUE"))
            {
              corteome=round(median(as.double(resini[,3])),digits=2)
            }
            else
            {
              corteome=as.double(tclvalue(lbl_corteome))
            }
            nc_pal=resini[(as.integer(resini[,2])>=cortefreq&as.double(resini[,3])<corteome),]
            nc_dist=rep(0,nrow(dados))
            for(i in 1:nrow(dados))
            {
              evocs=matrix("",max(length(nc_pal),5),2)
              evocs[,1]=c(nc_pal,rep("",max(length(nc_pal),5)-length(nc_pal)))
              evocs[,2]=c(t(dados[i,2:6]),rep("",max(length(nc_pal),5)-5))  
              colnames(evocs)=c("Central Nucleus","Data")
              evocs=as.data.frame(evocs)
              corp=VCorpus(VectorSource(evocs),readerControl=list(language="en"))
              dtm=DocumentTermMatrix(corp)
              mydtm=as.matrix(dtm)
              distancia=as.matrix(dist(mydtm))
              nc_dist[i]=distancia[2,1]  
            }  
            nc_infdist=round(summary(nc_dist),digits=2)
            plot6_pdf=ggplot(data.frame(x=c(20, 150)),aes(x))
            ggtitle("Summary Statistics Nucleus Central Distances for Individuals")
            tbl6_pdf=tableGrob(nc_infdist,rows=NULL)
            grid.arrange(tbl6_pdf,as.table=TRUE)
          }  	
          dev.off()
          tkmessageBox(message=paste("PDF successfully generated in the default Software R directory."),icon="warning",type="ok",default="ok",title="Message")
          tkdestroy(tt3)
          tkfocus(tt)
        }
        else
        {
          tkmessageBox(message=paste("Please select at least one result to export and generate the report."),icon="warning",type="ok",default="ok",title="Message")
        }
      }
      OnCancelar=function()
      {      
        tkdestroy(tt3)
        tkfocus(tt)
      }	
      gera_but=tkbutton(tt3,text="Generate PDF",command=OnGerar)
      canc_but=tkbutton(tt3,text="Cancel",command=OnCancelar)
      tkgrid(tklabel(tt3,text=" "),tklabel(tt3,text=" "),tklabel(tt3,text=" "))
      tkgrid(tklabel(tt3,text=" "),gera_but,tklabel(tt3,text=" "),canc_but,tklabel(tt3,text=" "))
      tkgrid(tklabel(tt3,text=" "),tklabel(tt3,text=" "),tklabel(tt3,text=" "))       
    }
    
    OnSair=function()
    {
      if((exists("tt2")==TRUE)&&(tt2!=""))
      {
        tkdestroy(tt2)
      }
      if((exists("tt3")==TRUE)&&(tt3!=""))
      {
        tkdestroy(tt3)
      }      
      tkdestroy(tt)      
    }
    
    mais_but=tkbutton(tt,text="More Options...",command=OnMaisOpcoes)
    expo_but=tkbutton(tt,text="Export Results",command=OnExportar)
    sair_but=tkbutton(tt,text="Exit",command=OnSair)
    tkgrid(tklabel(tt,text=" "),tklabel(tt,text=" "),tklabel(tt,text=" "))
    tkgrid(mais_but,expo_but,sair_but) 
    tkgrid(tklabel(tt,text=" "),tklabel(tt,text=" "),tklabel(tt,text=" "))  
  }
  
  # Stopping with the loading of the package
  OnDesistir=function()
  {
    tkdestroy(ttLoad)  
  }
  
  carrega_but=tkbutton(ttLoad,text="Start Package",command=OnCarregar)
  desiste_but=tkbutton(ttLoad,text="End Package",command=OnDesistir)  
  tkgrid(label1Load,label2Load)
  tkgrid(tklabel(ttLoad,text=" "))
  tkgrid(desiste_but,carrega_but)
  tkgrid(tklabel(ttLoad,text=" "))  
  
}