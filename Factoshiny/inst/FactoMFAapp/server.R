# server2. AFM

shinyServer(
  function(input, output) {
    values=reactive({
    nbfreq=0
    nbquali=0
    nbquanti=0
    quantisup=0
    listquali=c()
    gsup=c()
    dataselec=x[,input$variables1]
    groupe=length(input$variables1)
    nbgroupe=1
    if(input$typeG1=="quant"){
      nbquanti=nbquanti+1
      if(input$typeG12=="sup"){
        quantisup=quantisup+1
      }
      if(input$scale1=="sc"){
        typ="s"
      }
      if(input$scale1=="un"){
        typ="c"
      }
    }
    else if(input$typeG1=="qual"){
      typ="n"
      nbquali=nbquali+1
      listquali=c(listquali,input$variables1)
    }
    else if(input$typeG1=="freq"){
      typ="f"
      nbfreq=nbfreq+1
    }
    if(input$typeG12=="sup"){
      gsup=c(gsup,1)
    }
    types=typ
    nom=input$nameG1
    if(input$activeG2==TRUE && length(input$variables2)>0){
      dataselec=cbind(dataselec,x[,input$variables2])
      if(length(input$variables2)==1){
        colnames(dataselec)[dim(dataselec)[2]]=input$variables2
      }
      groupe=c(groupe,length(input$variables2))
      nbgroupe=2
      if(input$typeG2=="quant"){
        nbquanti=nbquanti+1
        if(input$typeG22=="sup"){
          quantisup=quantisup+1
        }
        if(input$scale2=="sc"){
          typ="s"
        }
        if(input$scale2=="un"){
          typ="c"
        }
      }
      else if(input$typeG2=="qual"){
        typ="n"
        nbquali=nbquali+1
        listquali=c(listquali,input$variables2)
      }
      else if(input$typeG2=="freq"){
        typ="f"
        nbfreq=nbfreq+1
      }
      if(input$typeG22=="sup"){
        gsup=c(gsup,2)
      }
      types=c(types,typ)
      nom=c(nom,input$nameG2)
      if(input$activeG3==TRUE && length(input$variables3)>0){
        dataselec=cbind(dataselec,x[,input$variables3])
        if(length(input$variables3)==1)
        {
          colnames(dataselec)[dim(dataselec)[2]]=input$variables3
        }
        groupe=c(groupe,as.numeric(length(input$variables3)))
        nbgroupe=3
        if(input$typeG3=="quant"){
          nbquanti=nbquanti+1
          if(input$typeG32=="sup"){
            quantisup=quantisup+1
          }
          if(input$scale3=="sc"){
            typ="s"
          }
          if(input$scale3=="un"){
            typ="c"
          }
        }
        else if(input$typeG3=="qual"){
          typ="n"
          nbquali=nbquali+1
          listquali=c(listquali,input$variables3)
        }
        else if(input$typeG3=="freq"){
          typ="f"
          nbfreq=nbfreq+1
        }
        if(input$typeG32=="sup"){
          gsup=c(gsup,3)
        }
        types=c(types,typ)
        nom=c(nom,input$nameG3)
        if(input$activeG4==TRUE && length(input$variables4)>0){
          dataselec=cbind(dataselec,x[,input$variables4])
          if(length(input$variables4)==1)
          {
            colnames(dataselec)[dim(dataselec)[2]]=input$variables4
          }
          groupe=c(groupe,length(input$variables4))
          nbgroupe=4
          if(input$typeG4=="quant"){
            nbquanti=nbquanti+1
            if(input$typeG42=="sup"){
              quantisup=quantisup+1
            }
            if(input$scale4=="sc"){
              typ="s"
            }
            if(input$scale4=="un"){
              typ="c"
            }
          }
          else if(input$typeG4=="qual"){
            typ="n"
            nbquali=nbquali+1
            listquali=c(listquali,input$variables4)
          }
          else if(input$typeG4=="freq"){
            typ="f"
            nbfreq=nbfreq+1
          }
          if(input$typeG42=="sup"){
            gsup=c(gsup,4)
          }
          types=c(types,typ)
          nom=c(nom,input$nameG4)
          if(input$activeG5==TRUE && length(input$variables5)>0){
            dataselec=cbind(dataselec,x[,input$variables5])
            if(length(input$variables5)==1)
            {
              colnames(dataselec)[dim(dataselec)[2]]=input$variables5
            }
            groupe=c(groupe,length(input$variables5))
            nbgroupe=5
            if(input$typeG5=="quant"){
              nbquanti=nbquanti+1
              if(input$typeG52=="sup"){
                quantisup=quantisup+1
              }
              if(input$scale5=="sc"){
                typ="s"
              }
              if(input$scale5=="un"){
                typ="c"
              }
            }
            else if(input$typeG5=="qual"){
              typ="n"
              nbquali=nbquali+1
              listquali=c(listquali,input$variables5)
            }
            else if(input$typeG5=="freq"){
              typ="f"
              nbfreq=nbfreq+1
            }
            if(input$typeG52=="sup"){
              gsup=c(gsup,5)
            }
            types=c(types,typ)
            nom=c(nom,input$nameG5)
            if(input$activeG6==TRUE && length(input$variables6)>0){
              dataselec=cbind(dataselec,x[,input$variables6])
              if(length(input$variables6)==1)
              {
                colnames(dataselec)[dim(dataselec)[2]]=input$variables6
              }
              groupe=c(groupe,length(input$variables6))
              nbgroupe=6
              if(input$typeG6=="quant"){
                nbquanti=nbquanti+1
                if(input$typeG62=="sup"){
                  quantisup=quantisup+1
                }
                if(input$scale6=="sc"){
                  typ="s"
                }
                if(input$scale6=="un"){
                  typ="c"
                }
              }
              else if(input$typeG6=="qual"){
                typ="n"
                nbquali=nbquali+1
                listquali=c(listquali,input$variables6)
              }
              else if(input$typeG6=="freq"){
                typ="f"
                nbfreq=nbfreq+1
              }
              if(input$typeG62=="sup"){
                gsup=c(gsup,6)
              }
              types=c(types,typ)
              nom=c(nom,input$nameG6)
              if(input$activeG7==TRUE && length(input$variables7)>0){
                dataselec=cbind(dataselec,x[,input$variables7])
                if(length(input$variables7)==1)
                {
                  colnames(dataselec)[dim(dataselec)[2]]=input$variables7
                }
                groupe=c(groupe,length(input$variables7))
                nbgroupe=7
                if(input$typeG7=="quant"){
                  nbquanti=nbquanti+1
                  if(input$typeG72=="sup"){
                    quantisup=quantisup+1
                  }
                  if(input$scale7=="sc"){
                    typ="s"
                  }
                  if(input$scale7=="un"){
                    typ="c"
                  }
                }
                else if(input$typeG7=="qual"){
                  typ="n"
                  nbquali=nbquali+1
                  listquali=c(listquali,input$variables7)
                }
                else if(input$typeG7=="freq"){
                  typ="f"
                  nbfreq=nbfreq+1
                }
                if(input$typeG72=="sup"){
                  gsup=c(gsup,7)
                }
                types=c(types,typ)
                nom=c(nom,input$nameG7)
                if(input$activeG8==TRUE && length(input$variables8)>0){
                  dataselec=cbind(dataselec,x[,input$variables8])
                  if(length(input$variables8)==1)
                  {
                    colnames(dataselec)[dim(dataselec)[2]]=input$variables8
                  }
                  groupe=c(groupe,length(input$variables8))
                  nbgroupe=8
                  if(input$typeG8=="quant"){
                    nbquanti=nbquanti+1
                    if(input$typeG82=="sup"){
                      quantisup=quantisup+1
                    }
                    if(input$scale8=="sc"){
                      typ="s"
                    }
                    if(input$scale8=="un"){
                      typ="c"
                    }
                  }
                  else if(input$typeG8=="qual"){
                    typ="n"
                    nbquali=nbquali+1
                    listquali=c(listquali,input$variables8)
                  }
                  else if(input$typeG8=="freq"){
                    typ="f"
                    nbfreq=nbfreq+1
                  }
                  if(input$typeG82=="sup"){
                    gsup=c(gsup,8)
                  }
                  types=c(types,typ)
                  nom=c(nom,input$nameG8)
                  if(input$activeG9==TRUE && length(input$variables9)>0){
                    dataselec=cbind(dataselec,x[,input$variables9])
                    if(length(input$variables9)==1)
                    {
                      colnames(dataselec)[dim(dataselec)[2]]=input$variables9
                    }
                    groupe=c(groupe,length(input$variables9))
                    nbgroupe=9
                    if(input$typeG9=="quant"){
                      nbquanti=nbquanti+1
                      if(input$typeG92=="sup"){
                        quantisup=quantisup+1
                      }
                      if(input$scale9=="sc"){
                        typ="s"
                      }
                      if(input$scale9=="un"){
                        typ="c"
                      }
                    }
                    else if(input$typeG9=="qual"){
                      typ="n"
                      nbquali=nbquali+1
                      listquali=c(listquali,input$variables9)
                    }
                    else if(input$typeG9=="freq"){
                      typ="f"
                      nbfreq=nbfreq+1
                    }
                    if(input$typeG92=="sup"){
                      gsup=c(gsup,9)
                    }
                    types=c(types,typ)
                    nom=c(nom,input$nameG9)
                    if(input$activeG10==TRUE && length(input$variables10)>0){
                      dataselec=cbind(dataselec,x[,input$variables10])
                      if(length(input$variables10)==1)
                      {
                        colnames(dataselec)[dim(dataselec)[2]]=input$variables10
                      }
                      groupe=c(groupe,length(input$variables10))
                      nbgroupe=10
                      if(input$typeG10=="quant"){
                        nbquanti=nbquanti+1
                        if(input$typeG102=="sup"){
                          quantisup=quantisup+1
                        }
                        if(input$scale10=="sc"){
                          typ="s"
                        }
                        if(input$scale10=="un"){
                          typ="c"
                        }
                      }
                      else if(input$typeG10=="qual"){
                        typ="n"
                        nbquali=nbquali+1
                        listquali=c(listquali,input$variables10)
                      }
                      else if(input$typeG10=="freq"){
                        typ="f"
                        nbfreq=nbfreq+1
                      }
                      if(input$typeG102=="sup"){
                        gsup=c(gsup,10)
                      }
                      types=c(types,typ)
                      nom=c(nom,input$nameG10)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if(length(gsup)==0){
      gsup=NULL
    }
    if(length(input$variables1)==1)
    {
      colnames(dataselec)[1]=input$variables1
    }
    list(res.MFA=(MFA(base=dataselec,group=groupe,type=types,name.group=nom,ncp=5,num.group.sup=gsup,graph=FALSE)),DATA=(dataselec),GROUPE=(groupe),NB=(nbgroupe),TYPE=(types),NBFREQ=(nbfreq),NBQUALI=(nbquali),NBQUANTI=(nbquanti),SUP=(gsup),QUANTISUP=(quantisup),LISTQUALI=(listquali))
    })
    
    error=function(){
      if(length(input$variables1)!=0 && length(input$variables2)!=0){
        etat="ok"
      }
      else{
       etat="not"
      }
      return(etat)
    }
    
    Plot1 <- function(){
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )      
      if(input$choixpartial==1){
        part=NULL
      }
      if(input$choixpartial==2){
        part="all"
      }
      if(input$choixpartial==3){
        part=input$indivpartiel
      }
      lapbar=TRUE
      if(input$choixpartial!=1 && input$partind==FALSE){
        lapbar=FALSE
      }
      if(input$choixpartial==1 && input$drawind=="a"){
        habi="group"
      }
      else if(input$choixpartial==1 && input$drawind=="b"){
        habi="ind"
      }
      else if((input$choixpartial==2 || input$choixpartial==3) && input$drawind=="a"){
        habi="ind"
      }
      else if((input$choixpartial==2 || input$choixpartial==3) && input$drawind=="b"){
        habi="group"
      }
      else if(input$drawind=="c"){
        habi=input$habiquali
      }
      invi="none"
      if(input$meanind1==FALSE){
        invi="ind" 
      }
      else if (input$qualind1==FALSE){
        invi="quali"
      }
      if (input$qualind1==FALSE && input$meanind1==FALSE){
        invi=c("ind","quali")
      }
      if(!(is.null(habi))){
      plot.MFA(values()$res.MFA,choix="ind",axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),title=input$title2,partial=part,lab.ind=input$meanind, lab.par=lapbar,lab.var=input$qualind, habillage=habi,invisible=invi)
      }
    }
    
    output$map <- renderPlot({
      p <- Plot1()
    })
    
    output$map2 <- renderPlot({
      
      if(input$colorgroup==TRUE){
        habi="group"
      }
      if(input$colorgroup==FALSE){
        habi="none"
      }
      if(input$selection=="no"){
        selec=NULL
      }
      if(input$selection=="contrib"){
        selec=paste("contrib ",input$slider2)
      }
      if(input$selection=="cos2"){
        if(input$slider3!=1){
          selec=paste("cos2 ",input$slider3)
        }
        else{
          selec="cos2 0.999"
        }
      }
      invi="none"
      plot.MFA(values()$res.MFA,choix="var",axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),title=input$title3,habillage=habi,select=selec,invisible=invi)
    })
    
    output$map22=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(values()$NBQUANTI!=0,"No quantitative group")
      )
      if(values()$NBQUANTI==0){
        return(p("No quantitative variable"))
      }
        else{
          plotOutput("map2", width = 500, height=500)
        }
    })
    
    Plot5 <- function(){
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      plot.MFA(values()$res.MFA,choix="group",title=input$title1,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)))
    }
    
    output$map5 <- renderPlot({
      p <- Plot5()
    })
    
    Plot4 <- function(){
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      if(input$coloraxe==TRUE){
        habi="group"
      }
      else{
        habi="none"
      }
      plot.MFA(values()$res.MFA,choix="axes",axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),title=input$title4,habillage=habi)
    }
    
    output$map4 <- renderPlot({
      p <- Plot4()
    })

    output$map6 <- renderPlot({
      if(input$affichcol==TRUE){
        col=TRUE
      }
      if(input$affichcol==FALSE){
        col=FALSE
      }
      plot.MFA(values()$res.MFA,choix="freq",axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),title=input$title5,lab.col=col)
    })
    
    output$map66=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      if(values()$NBFREQ ==0){
        return(p("No group of frequencies"))
      }
      else{
        return(plotOutput("map6", width = 500, height=500))
      }
    })

    
    output$drawindiv=renderUI({
      if(input$choixpartial==1){
        return(radioButtons("drawind","Select drawing",choices=list("No selection"= "a","By individual"="b","By categorical variable"="c"),inline=TRUE))
      }
      else{
        return(radioButtons("drawind","Select drawing",choices=list("By individual"= "a","By group"="b","By categorical variable"="c"),inline=TRUE))
      }
    })
    
    
    
    output$habillagequali=renderUI({
      if(input$activeG2==TRUE && length(input$variables2)>0){
        if(!(is.null(values()$LISTQUALI))){
          choix=values()$LISTQUALI
          if(length(choix)==1){
          return(selectInput("habiquali"," ",choices=choix))
          }
          else{
          num=c(1:length(choix))
          return(selectInput("habiquali"," ",choices=list(num=choix)))
          }
        }
      }
      else{
        p("No group of categorical variable")
      }
    })
    ###Recup codes
    observe({
      if(input$HCPCcode==0){
      }
      else {
        isolate({
          cat(codeGraph1(),sep="\n")
          cat(codeGraph2(),sep="\n")
          cat(codeGraph3(),sep="\n")
          cat(codeGraph4(),sep="\n")
          if(!is.null(code$freq)){
            cat(codeGraph5(),sep="\n")
          }
        })
      }
    })
    
    codeGraph1<-function(){
      if(input$choixpartial==1){
        part="NULL"
      }
      if(input$choixpartial==2){
        part="all"
      }
      if(input$choixpartial==3){
        part1=input$indivpartiel
        if(length(input$indivpartiel)==1){
          part=paste("'",part1,"'")
        }
        if(length(input$indivpartiel)>1){
          vec4=NULL
          vec4<-paste(vec4,"'",input$indivpartiel[1],"'",sep="")
          for (i in 2:(length(input$indivpartiel))){
            vec4<-paste(vec4,paste("'",input$indivpartiel[i],"'",sep=""),sep=",")
          }
          part=paste("c(",vec4,")",sep="")
        }
      }
      lapbar=TRUE
      if(input$choixpartial!=1 && input$partind==FALSE){
        lapbar=FALSE
      }
      habi="none"
      if(!(is.null(input$drawind))){
        if(input$choixpartial==1 && input$drawind=="a"){
          habi="group"
        }
        else if(input$choixpartial==1 && input$drawind=="b"){
          habi="ind"
        }
        else if((input$choixpartial==2 || input$choixpartial==3) && input$drawind=="a"){
          habi="ind"
        }
        else if((input$choixpartial==2 || input$choixpartial==3) && input$drawind=="b"){
          habi="group"
        }
        else if(input$drawind=="c"){
          habi=input$habiquali
        }
      }
      invi="none"
      if(input$meanind1==FALSE){
        invi="ind" 
      }
      else if (input$qualind1==FALSE){
        invi="quali"
      }
      if (input$qualind1==FALSE && input$meanind1==FALSE){
        invi=c("ind","quali")
      }
      Call1=as.name(paste("plot.MFA(res,choix='ind',axes=c(",input$nb1,",",input$nb2,"),partial=",part,",title='",input$title2,"',lab.ind=",input$meanind,",lab.par=",lapbar,",lab.var=",input$qualind, ",habillage='",habi,"',invisible='",invi,"')",sep=""))
      return(Call1)
    }
    
    codeGraph2<-function(){
      if(input$colorgroup==TRUE){
        habi="group"
      }
      if(input$colorgroup==FALSE){
        habi="none"
      }
      if(input$selection=="no"){
        selec="NULL"
      }
      if(input$selection=="contrib"){
        selec=paste("contrib ",input$slider2)
      }
      if(input$selection=="cos2"){
        if(input$slider3!=1){
          selec=paste("cos2 ",input$slider3)
        }
        else{
          selec="cos2 0.999"
        }
      }
      invi="none"
      Call2=paste("plot.MFA(res,choix='var',axes=c(",input$nb1,",",input$nb2,"),habillage='",habi,"',title='",input$title3,"',select=",selec,",invisible='",invi,"')",sep="")
      return(Call2)
    }
    
    codeGraph3<-function(){
      Call3=paste("plot.MFA(res,choix='group',title='",input$title1,"',axes=c(",input$nb1,",",input$nb2,"))",sep="")
      return(Call3)
    }
    
    codeGraph4<-function(){
      if(input$coloraxe==TRUE){
        habi="group"
      }
      else{
        habi="none"
      }
      Call4=paste("plot.MFA(res,choix='axes',title='",input$title4,"',axes=c(",input$nb1,",",input$nb2,"),habillage='",habi,"')",sep="")
      return(Call4)
    }
    
    codeGraph5<-function(){
      if(input$affichcol==TRUE){
        col=TRUE
      }
      if(input$affichcol==FALSE){
        col=FALSE
      }
      Call5=paste("plot.MFA(res,choix='freq',title='",input$title5,"',axes=c(",input$nb1,",",input$nb2,"),lab.col=",col,")",sep="")
      return(Call5)
    }
    
    
    ### Recuperation des parametres
    observe({
      if(input$Quit==0){
      }
      else{
        isolate({
          stopApp(returnValue=valeuretour())
        })
      }
    })
    
    valeuretour=function(){
      res=list()
      res$code=values()$res.MFA
      res$axe1=input$nb1
      res$axe2=input$nb2
      res$ind1=input$meanind1
      res$ind2=input$meanind
      res$ind3=input$qualind1
      res$ind4=input$qualind
      res$drawing=input$drawind
      res$drawing2=input$habiquali
      res$partial=input$choixpartial
      res$partial2=input$indivpartiel
      res$partial3=input$partind
      res$selectvar=input$selection
      sel=NULL
      if(input$selection=="contrib"){
        sel=input$slider2
      }
      if(input$selection=="cos2"){
        sel=input$slider3
      }
      res$selectvar2=sel
      res$hide=input$hides
      res$colorvar=input$colorgroup
      res$freq1=input$affichind
      res$freq2=input$affichcol
      res$partaxe=input$coloraxe
      res$nom=nomData
      res$code1=codeGraph1()
      res$code2=codeGraph2()
      res$code3=codeGraph3()
      res$code4=codeGraph4()
      res$code5=codeGraph5()
      res$title1=input$title1
      res$title2=input$title2
      res$title3=input$title3
      res$title4=input$title4
      res$title5=input$title5
      class(res)="MFAshiny"
      return(res)
    }
    
    output$slider1=renderUI({
      etat2=error()
      validate(
        need(etat2!="not"," ")
      )
      maxlength=dim(values()$res.MFA$quanti.var$coord)[1]
      if(input$selection=="contrib"){
        return(sliderInput("slider2","Number of the most contributive variables",min=1, max=maxlength, value=maxlength, step=1))
      }
      if(input$selection=="cos2"){
        return(sliderInput("slider3","Number of variables with highest cos2",min=0, max=maxlength, value=maxlength, step=1))
      }
    })
    
    output$hide2=renderUI({
      etat2=error()
      validate(
        need(etat2!="not"," ")
      )
      if(values()$QUANTISUP!=0){
        return(radioButtons("hides","Hide :",choices=list("Nothing"="non","Active variables"="act","Supplementary variables"="sup"),selected="non"))
      }
    })
    
    output$sorties=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$eig))
    })
    
    output$map3=renderPlot({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(barplot(values()$res.MFA$eig[,1],names.arg=rownames(values()$res.MFA$eig),las=2))
    })
    output$JDD=renderDataTable({
      cbind(Names=rownames(x),x)},
      options = list(    "orderClasses" = TRUE,
                         "responsive" = TRUE,
                         "pageLength" = 10))
    output$summary=renderPrint({
      summary(x)
    })
    output$summaryMFA=renderPrint({
      etat2=error()
      validate(
        need(etat2!="not"," ")
      )
      summary.MFA(values()$res.MFA)
    })  
    
      
    output$histo=renderPlot({
      par(mfrow=c(1,2))
      boxplot(x[,input$bam])
      plot(density(x[,input$bam]),main="",xlab="")
    })
    
    
    output$downloadData = downloadHandler(
      filename = function() { 
        paste('graph1','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot1()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData3 = downloadHandler(
      filename = function() { 
        paste('graph2','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot2()
        dev.off()
      },
      contentType='image/png')
    
    output$download3=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","")
      )
      if(values()$NBQUANTI==0){
        return()
      }
      else{
        return(downloadButton("downloadData3","Download as png"))
      }
    })
    
    output$downloadData11 = downloadHandler(
      filename = function() { 
        paste('graph3','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot5()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData12 = downloadHandler(
      filename = function() { 
        paste('graph3','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot5()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData13 = downloadHandler(
      filename = function() { 
        paste('graph3','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot5()
        dev.off()
      },
      contentType=NA)
    
    
    output$downloadData15 = downloadHandler(
      filename = function() { 
        paste('graph4','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot4()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData16 = downloadHandler(
      filename = function() { 
        paste('graph4','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot4()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData17 = downloadHandler(
      filename = function() { 
        paste('graph4','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot4()
        dev.off()
      },
      contentType=NA)
    
    
    output$downloadData19 = downloadHandler(
      filename = function() { 
        paste('graph5','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot6()
        dev.off()
      },
      contentType='image/png')
    
    output$download19=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","")
      )
      if(values()$NBFREQ==0){
        return()
      }
      else{
        return(downloadButton("downloadData19","Download as png"))
      }
    })
    
    output$downloadData20 = downloadHandler(
      filename = function() { 
        paste('graph5','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot6()
        dev.off()
      },
      contentType='image/jpg')
    
    output$download20=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","")
      )
      if(values()$NBFREQ==0){
        return()
      }
      else{
        return(downloadButton("downloadData20","Download as jpg"))
      }
    })
    
    output$downloadData21 = downloadHandler(
      filename = function() { 
        paste('graph5','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot6()
        dev.off()
      },
      contentType=NA)
    
    output$download21=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","")
      )
      if(values()$NBFREQ==0){
        return()
      }
      else{
        return(downloadButton("downloadData21","Download as pdf"))
      }
    })
    
    output$downloadData22 = downloadHandler(
      filename = function() { 
        paste('graph5','.emf', sep='') 
      },
      content = function(file) {
        emf(file)
        Plot6()
        dev.off()
      },
      contentType=NA)
    
    output$downloadData1 = downloadHandler(
      filename = function() { 
        paste('graph1','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot1()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData2 = downloadHandler(
      filename = function() { 
        paste('graph1','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot1()
        dev.off()
      },
      contentType=NA)
    
    output$downloadData4 = downloadHandler(
      filename = function() { 
        paste('graph2','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot2()
        dev.off()
      },
      contentType='image/jpg')
    
    output$download4=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","")
      )
      if(values()$NBQUANTI==0){
        return()
      }
      else{
        return(downloadButton("downloadData4","Download as jpg"))
      }
    })
    
    output$downloadData5 = downloadHandler(
      filename = function() { 
        paste('graph2','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot2()
        dev.off()
      },
      contentType=NA)
    
    output$download5=renderUI({
      etat2=error()
      validate(
        need(etat2!="not","")
      )
      if(values()$NBQUANTI==0){
        return()
      }
      else{
        return(downloadButton("downloadData5","Download as pdf"))
      }
    })
    
    output$downloadData6 = downloadHandler(
      filename = function() { 
        paste('graph2','.emf', sep='') 
      },
      content = function(file) {
        emf(file)
        Plot2()
        dev.off()
      },
      contentType=NA)
    
    output$downloadData7 = downloadHandler(
      filename = function() { 
        paste('graph1','.emf', sep='') 
      },
      content = function(file) {
        emf(file)
        Plot1()
        dev.off()
      },
      contentType=NA)
    
    ### Sorties
    
    output$sorties1=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$ind$coord))
    })
    
    output$sorties2=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$ind$contrib))
    })
    
    output$sorties3=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$ind$cos2))
    })
    
    output$sorties4=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$ind$within.inertia))
    })
    
    output$sorties5=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$ind$coord.partiel))
    })
    
    output$sorties6=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$ind$within.partial.inertia))
    })
    
    output$sorties11=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$quanti.var$coord))
    })
    
    output$sorties22=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$quanti.var$contrib))
    })
    
    output$sorties33=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$quanti.var$cos2))
    })
    
    output$sorties44=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$quanti.var$cor))
    })
    
    output$sorties12=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$partial.axes$coord))
    })
    
    output$sorties23=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$partial.axes$cor))
    })
    
    output$sorties34=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$partial.axes$contrib))
    })
    
    output$sorties45=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      return(as.data.frame(values()$res.MFA$partial.axes$cor.between))
    })
    
    
    
    output$sortiegroup=renderTable({
      etat2=error()
      validate(
        need(etat2!="not","Please select at least 2 groups")
      )
      write.infile(X=values()$res.MFA$group,file=paste(getwd(),"fichgroup.csv"),sep=";",nb.dec=5)
      baba=read.csv(paste(getwd(),"fichgroup.csv"),sep=";",header=FALSE)
      colnames(baba)=NULL
      file.remove(paste(getwd(),"fichgroup.csv"))
      baba
    },
    include.rownames=FALSE)
    
    
   ### Fonction permettant d'afficher la liste des variables disponibles, en fonction du type 
    
    output$listvarG1=renderUI({
      if(input$typeG1=="quant" || input$typeG1=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG1=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables1",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG2=renderUI({
      if(input$typeG2=="quant" || input$typeG2=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG2=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables2",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG3=renderUI({
      if(input$typeG3=="quant" || input$typeG3=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG3=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables3",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG4=renderUI({
      if(input$typeG4=="quant" || input$typeG4=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG4=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables4",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG5=renderUI({
      if(input$typeG5=="quant" || input$typeG5=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG5=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables5",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG6=renderUI({
      if(input$typeG6=="quant" || input$typeG6=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG6=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables6",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG7=renderUI({
      if(input$typeG7=="quant" || input$typeG7=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG7=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables7",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG8=renderUI({
      if(input$typeG8=="quant" || input$typeG8=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG8=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables8",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG9=renderUI({
      if(input$typeG9=="quant" || input$typeG9=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG9=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables9",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })
    
    output$listvarG10=renderUI({
      if(input$typeG10=="quant" || input$typeG10=="freq"){
        if(length(quanti)>1){
          choix=list(IdChoices=VariableChoices)
        }
        if(length(quanti)==1){
          choix=quanti
        }
      }
      if(input$typeG10=="qual"){
        if(length(quali)>1){
          choix=list(Idqualisup=QualiChoice)
        }
        if(length(quali)==1){
          choix=quali
        }
      }
      return(selectInput("variables10",label="Choose variables",
                         choices=choix,multiple=TRUE,selectize=TRUE))
    })

  }
)
      

