# server script for FAMD2
shinyServer(
  function(input, output) {
    values=reactive({
    if (input$selecactive=="Toutes"){
      data.selec=newdata
    }
    else{
      validate(
        need(length(input$supvar)>0 || length(input$supvar1)>0, "Please select at least one supplementary variables")
      )
      data.selec=newdata
    }
    choixsup=getactive()$sup
    if(length(input$indsup)==0){
      suple=NULL
    }
    else{
      suple=c()
      for (i in 1:length(nom)){
        if(nom[i]%in%input$indsup){
          suple=c(suple,i)
        }
      }
    }
    list(res.FAMD=(FAMD(data.selec,sup.var=choixsup,ind.sup=suple,graph=FALSE,ncp=5)),DATA=(data.selec),choixsuple=(suple),varsup=(choixsup))
    })
    
    Plot1 <- reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$select0=="cos2"){
        if(input$slider00!=1){
          selecindiv=paste("cos2 ",input$slider00)
        }
        else{
          selecindiv="cos2 0.999"
        }
        selecindivText=paste("'",selecindiv,"'",sep="")
      }
      if(input$select0=="NONE"){
        selecindiv=NULL
        selecindivText="NULL"
      }
      if(input$select0=="contrib"){
        selecindiv=paste("contrib ",input$slider4)
        selecindivText=paste("'",selecindiv,"'",sep="")
      }
      list(PLOT=(plot.FAMD(values()$res.FAMD,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="var",cex=input$cex2,cex.main=input$cex2,cex.axis=input$cex2,title=input$title2,select=selecindiv)),selecindivText=selecindivText)
    })
    
    output$map <- renderPlot({
      p <- Plot1()$PLOT
    })
    
    Plot2 <- reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      hab="none"
      if(length(QualiChoice)==0){
        hab="none"
      }
      else if(length(QualiChoice)==1 && input$habi==TRUE){
        if(input$habi==TRUE){
          hab=QualiChoice
          colquali="blue"
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      else if (length(QualiChoice)>1){
        if(input$habi==TRUE){
          if(is.null(input$habiller)){
            hab="none"
          }
          else{
          hab=input$habiller
          }
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      if(hab!="none"){
        hab=which(all==hab)
        hab=as.numeric(hab)
      }
      if(input$select=="cos2"){
        if(input$slider1!=1){
          selecindiv=paste("cos2 ",input$slider1)
        }
        else{
          selecindiv="cos2 0.999"
        }
        selecindivText=paste0("'",selecindiv,"'")
      }
      if(input$select=="NONE"){
        selecindiv=NULL
        selecindivText="NULL"
      }
      if(input$select=="contrib"){
        selecindiv=paste("contrib ",input$slider0)
        selecindivText=paste0("'",selecindiv,"'")
      }
      if(input$select=="Manuel"){
        selecindiv=c(input$indiv)
        selecindivText=paste0("'",selecindiv,"'")
      }
      list(PLOT=(plot.FAMD(values()$res.FAMD,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="ind",select=selecindiv,lab.var=input$labels,lab.ind=input$labels2,cex=input$cex,cex.main=input$cex,cex.axis=input$cex,habillage=hab,title=input$title1)),selecindivText=selecindivText,HABILLAGE=hab)
      
    })
    
    output$map2 <- renderPlot({
      p <- Plot2()$PLOT
    })
   
    
    Plot4 <- reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecti=="cos2"){
       if(input$slider000!=1){
        selecindiv=paste("cos2 ",input$slider000)
      }
      else{
       selecindiv="cos2 0.999"
      }
      selecindivText=paste("'",selecindiv,"'",sep="")
      }
      if(input$selecti=="NONE"){
       selecindiv=NULL
      selecindivText="NULL"
      }
      if(input$selecti=="contrib"){
       selecindiv=paste("contrib ",input$slider6)
      paste("'",selecindiv,"'",sep="")
      }
      list(PLOT=(plot.FAMD(values()$res.FAMD,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="quanti",cex=input$cex3,cex.main=input$cex3,cex.axis=input$cex3,title=input$title3,select=selecindiv)))      
    })
    
    output$map4 <- renderPlot({
      p <- Plot4()$PLOT
    })
    
    
    
    ### Boutton pour quitter l'application
    ### Recuperation parametres
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
      res$nomData=nomData
      res$data=values()$DATA
      res$b=input$supvar
      res$c=input$supvar1
      res$d=input$indsup
      res$e=input$nb1
      res$f=input$nb2
      habillage=NULL
      if(input$habi==TRUE){
        habillage=input$habiller
      }
      res$g=habillage
      if(input$select=="cos2"){
        selecindiv=input$slider1
      }
      if(input$select=="NONE"){
        selecindiv=NULL
      }
      if(input$select=="contrib"){
        selecindiv=input$slider0
      }
      if(input$select=="Manuel"){
        selecindiv=input$indiv
      }
      res$h=input$select
      res$i=selecindiv
      selecindiv2=NULL
      if(input$select0=="cos2"){
        selecindiv2=input$slider00
      }
      if(input$select0=="contrib"){
        selecindiv2=input$slider4
      }
      res$j=input$select0
      res$k=selecindiv2
      selecindiv3=NULL
      if(input$selecti=="cos2"){
        selecindiv3=input$slider000
      }
      if(input$selecti=="contrib"){
        selecindiv3=input$slider6
      }
      res$o=input$selecti
      res$p=selecindiv3
      res$l=input$cex
      res$m=input$cex2
      res$n=input$cex3
      res$code1=code()
      res$code2=codeGraphVar()
      res$code3=codeGraphInd()
      res$labind=input$labels2
      res$labvar=input$labels
      class(res) <- "FAMDshiny"
      return(res)
    }
    
    #### Fonction recuperation de code
    
    observe({
      if(input$FAMDcode==0){
      }
      else {
        isolate({
          if (length(input$habiller)==2 & input$habi==TRUE){
            cat(paste("newCol<-paste(",nomData,"[,'",input$habiller[1],"'],",nomData,"[,'",input$habiller[2],"'],","sep='/')",sep=""),sep="\n")
          }
          cat(code(),sep="\n")
          cat(codeGraphVar(),sep="\n")
          cat(codeGraphInd(),sep="\n")
        })
      }
    })
  
    
    
    codeGraphVar<-function(){
      
      if(length(input$slider4)==0){
        selection="NULL"
      }
      else{
        selection=Plot1()$selecindivText
      }
      Call1=paste("plot.FAMD(res.FAMD,axes=c(",input$nb1,",",input$nb2,"),choix='var',select=",selection,",cex=",input$cex2,",cex.main=",input$cex2,",cex.axis=",input$cex2,",unselect=0,col.quanti.sup='red')",sep="")
      return(Call1)
    }
    
    codeGraphInd<-function(){
      hab="none"
      if (length(input$habiller)<=1 & input$habi==TRUE || input$habi==FALSE){
        hab=paste(Plot2()$HABILLAGE,sep="")
      }
      
      Call2=paste("plot.FAMD(res.FAMD,","axes=c(",input$nb1,",",input$nb2,"),choix='ind',select=",Plot2()$selecindivText,",habillage=",hab,",cex=",input$cex,",cex.main=",input$cex,",cex.axis=",input$cex,",lab.var=",input$labels,",lab.ind=",input$labels2,",title='",input$title1,"')",sep="")
      return(Call2)
    }
    
    ##### Fin de la fonction recuperation du code
    
    
    output$out22=renderUI({
      choix=list("Summary of FAMD"="ACP","Eigenvalues"="eig","Results of the variables"="resvar","Results of the individuals"="resind")
      if(!is.null(values()$choixsuple)){
        choix=c(choix,"Results of the supplementary individuals"="supind")
      }
      if(!is.null(values()$varsup)){
        choix=c(choix,"Results of the supplementary variables"="varsup")
      }
      radioButtons("out","Which outputs do you want ?",
                   choices=choix,selected="ACP",inline=TRUE)
    })
    
    getactive=function(){
      if(input$selecactive=="choix"){
      sup=c()
      sup2=c()
      sup3=c()
      if(length(input$supvar)==0&&length(input$supvar1)==0){
        activevar=all
        supl=NULL
      }
      else if(length(input$supvar1)==0&&length(input$supvar)!=0){
        for (i in 1:length(all)){
          if(all[i]%in%input$supvar){
            sup=c(sup,i)
          }
        }
        activevar=all[-sup]
        supl=VariableChoices[sup]
        quanti=VariableChoices[-sup]
      }
      else if(length(input$supvar)==0&&length(input$supvar1)!=0){
        for (i in 1:length(all)){
          if(all[i]%in%input$supvar1){
            sup=c(sup,i)
          }
        }
        activevar=all[-sup]
        supl=QualiChoice[sup]
        quali=QualiChoice[-sup]
      }
      else if(length(input$supvar)!=0&&length(input$supvar1)!=0){
        for (i in 1:length(all)){
          if(all[i]%in%input$supvar1 || all[i]%in%input$supvar){
            sup=c(sup,i)
          }
        }
        activevar=all[-sup]
        supl=all[sup]
        for (i in 1:length(QualiChoice)){
          if(QualiChoice[i]%in%input$supvar1){
            sup2=c(sup2,i)
          }
        }
        for (i in 1:length(VariableChoices)){
          if(VariableChoices[i]%in%input$supvar){
            sup3=c(sup3,i)
          }
        }
      quanti=QualiChoice[-sup2]
      quali=VariableChoices[-sup3]
      }
      ind=NULL
      if(!is.null(supl)){
        for(i in 1:length(supl)){
          ind=c(ind,which(all==supl[i]))
        }
      }
      return(list(activevar=activevar,supl=supl,quanti=quanti,quali=quali,sup=sup))
    }
  }
    
    
    output$NB1=renderUI({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes" || length(getactive()$activevar)>5){
        return(selectInput("nb1", label = h6("x axis"), 
                    choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = axe1,width='80%'))
      }
      else{
        baba=c(1:length(getactive()$activevar))
        return(selectInput("nb1",label=h6("x axis"), choices=baba,selected=axe1,width='80%'))
      }
    })
    
    output$NB2=renderUI({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes" || length(getactive()$activevar)>5){
        return(selectInput("nb2", label = h6("y axis"), 
                           choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = axe2,width='80%'))
      }
      else{
        baba=c(1:length(getactive()$activevar))
        return(selectInput("nb2",label=h6("y axis"), choices=baba,selected=axe2,width='80%'))
      }
    })
    
    output$sorties=renderTable({
        return(as.data.frame(values()$res.FAMD$eig))
    })
    
    output$sorties12=renderTable({
        validate(
          need((length(input$supquali)>0 || input$supquali==TRUE), "No categorical variables selected")
        )
        validate(
          need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
        )
        return(as.data.frame(values()$res.FAMD$quali.sup$coord))
    })
    
    output$sorties13=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      validate(
        need((length(input$supquali)>0 || input$supquali==TRUE), "No categorical variables selected")
      )
      return(as.data.frame(values()$res.FAMD$quali.sup$v.test))
    })
    
    output$sorties2=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$var$coord))
    })
    
    output$sorties22=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$ind$coord))
    })
    
    output$sorties23=renderTable({
      validate(
        need(length(input$supvar)!=0, "No supplementary quantitative variables")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$var$coord.sup))
    })
    
    output$sorties32=renderTable({
      validate(
        need(length(input$supvar)!=0, "No supplementary quantitative variables")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$var$cos2.sup))
    })
    
    output$sorties36=renderTable({
      validate(
        need(length(input$indsup)!=0, "No supplementary individuals")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$ind.sup$coord))
    })
    
    output$sorties37=renderTable({
      validate(
        need(length(input$indsup)!=0, "No supplementary individuals")
      )
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$ind.sup$cos2))
    })
    
    
    output$sorties3=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$var$contrib))
    })
    
    output$sorties33=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$ind$contrib))
    })
    
    output$sorties4=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$var$cos2))
    })
    
    output$sorties44=renderTable({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.FAMD$ind$cos2))
    })
  
  output$sortieDimdesc3=renderTable({
    validate(
      need(length(getactive()$activevar)>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.FAMD)[[1]]$quanti))
  })
  
  output$sortieDimdesc4=renderTable({
    validate(
      need(length(getactive()$activevar)>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.FAMD)[[1]]$quali))
  })
  
  #DIM2
  
  output$sortieDimdesc33=renderTable({
    validate(
      need(length(getactive()$activevar)>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.FAMD)[[2]]$quanti))
  })
  
  output$sortieDimdesc44=renderTable({
    validate(
      need(length(getactive()$activevar)>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.FAMD)[[2]]$quali))
  })
  
  #DIM3
  
  output$sortieDimdesc333=renderTable({
    validate(
      need(length(getactive()$activevar)>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.FAMD)[[3]]$quanti))
  })
  
  output$sortieDimdesc444=renderTable({
    validate(
      need(length(getactive()$activevar)>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.FAMD)[[3]]$quali))
  })
    
    
    output$map3=renderPlot({
      return(barplot(values()$res.FAMD$eig[,1],names.arg=rownames(values()$res.FAMD$eig),las=2))
    })
    
    output$JDD=renderDataTable({
      cbind(Names=rownames(newdata),newdata)},
      options = list(    "orderClasses" = TRUE,
                         "responsive" = TRUE,
                         "pageLength" = 10))
  
    output$summary=renderPrint({
      summary(newdata)
    })
  
  
  code<-function(){
    vec=nomData
    part2=""
    if(!is.null(input$supvar)||!is.null(input$supvar1)){
      choixsup=getactive()$sup
      vect3=NULL
      vect3<-paste(vect3,choixsup[1],sep="")
      for(i in 2:length(choixsup)){
        vect3<-paste(vect3,paste(choixsup[i],sep=""),sep=",")
      }
      part2=paste(",sup.var=c(",vect3,"),")
    }
    part3=""
    if(!is.null(input$indsup)){
        suple=c()
        for (i in 1:length(nom)){
          if(nom[i]%in%input$indsup){
            suple=c(suple,i)
            
          }
        }
        vect4=NULL
        vect4<-paste(vect4,suple[1],sep="")
        if(length(suple)>1){
        for(i in 2:length(suple)){
          vect4<-paste(vect4,paste(suple[i],sep=""),sep=",")
        } 
        }
      if(part2!=""){
        part3=paste("ind.sup=c(",vect4,")")
      }
      else{
        part3=paste(",ind.sup=c(",vect4,")")
      }
    }
    Call1=as.name(paste("res.FAMD<-FAMD(",vec,part2,part3,",graph=FALSE,ncp=5)",sep=""))
    return(Call1)
  }
  
  # Attention, si le nombre d'individus passe en dessous de 10, bug
    output$summaryFAMD=renderPrint({
      validate(
        need(input$nbele!=0, "Please select at least one element")
      )
      a<-values()$res.FAMD
      a$call$call<-code()
      summary.FAMD(a,nbelements=input$nbele)
    })
  
    output$summary2=downloadHandler(filename = function() { 
      paste('summaryofFAMD','.txt', sep='') 
    },
    content = function(file) {
      summary.FAMD(values()$res.FAMD,nbelements=input$nbele,file=file)
    },
    contentType='text/csv')
  
    
    output$slider3=renderUI({
      validate(
        need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes"){
        maxvar=length(all)
      }
      if(input$selecactive=="choix"){
        maxvar=length(getactive()$activevar)
      }
      if(selection3=="contrib"){
        return(div(align="center",sliderInput("slider4",label="Number of the most contributive variables",
                                              min=1,max=maxvar,value=selection4,step=1)))  
      }
      else{
      return(div(align="center",sliderInput("slider4",label="Number of the most contributive variables",
                  min=1,max=maxvar,value=maxvar,step=1)))}
    })
  
  output$slider5=renderUI({
    validate(
      need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
    )
    if(input$selecactive=="Toutes"){
      maxvar=length(quanti)
    }
    if(input$selecactive=="choix"){
      maxvar=length(getactive()$quanti)
    }
    if(selection5=="contrib"){
      return(div(align="center",sliderInput("slider6",label="Number of the most contributive variables",
                                            min=1,max=maxvar,value=selection6,step=1)))  
    }
    else{
      return(div(align="center",sliderInput("slider6",label="Number of the most contributive variables",
                                            min=1,max=maxvar,value=maxvar,step=1)))}
  })
  
  output$slider7=renderUI({
    validate(
      need(length(getactive()$activevar)>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
    )
    if(is.null(input$indsup)){
      maxi=length(nom)
    }
    if(!is.null(input$indsup)){
      maxi=length(nom)-length(input$indsup)
    }
    if(selection=="contrib"){
      return(div(align="center",sliderInput("slider0", label = "Number of the most contributive individuals",
                                     min = 1, max = maxi, value =as.numeric(selection2),step=1)))
    }
    else{
      return(div(align="center",sliderInput("slider0", label = "Number of the most contributive individuals",
                                     min = 1, max = maxi, value =maxi,step=1))) 
    }
  })

    
    output$habillage2=renderUI({
      if(length(QualiChoice)==0){
        return(p("No categorical variable"))
      }
      if(length(QualiChoice)>1){
        if(is.null(habillageind)){
        num=c(1:length(QualiChoice))
        return(selectInput("habiller","", choices=list(num=QualiChoice),multiple=FALSE))
        }
        else{
          num=c(1:length(QualiChoice))
          return(selectInput("habiller","", choices=list(num=QualiChoice),multiple=FALSE,selected=habillageind))
        }
      }
      if(length(QualiChoice)==1){
        if(is.null(habillageind)){
          return(selectInput("habiller","", choices=QualiChoice,multiple=FALSE))
        }
        else{
          return(selectInput("habiller","", choices=QualiChoice,multiple=FALSE,selected=habillageind))
        }
      }
    })
      
    output$histo=renderPlot({
      if(input$bam%in%quanti){
      par(mfrow=c(1,2))
      boxplot(newdata[,input$bam])
      hist(newdata[,input$bam],main="",xlab="")
      }
      else{
        barplot(table(newdata[,input$bam]),cex.names=0.8)
      }
    })
    
    
    
    output$downloadData = downloadHandler(
      filename = function() { 
        paste('graph1','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot11()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData1 = downloadHandler(
      filename = function() { 
        paste('graph1','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot11()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData2 = downloadHandler(
      filename = function() { 
        paste('graph1','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot11()
        dev.off()
      },
      contentType=NA)
    
    output$downloadData3 = downloadHandler(
      filename = function() { 
        paste('graph2','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot22()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData4 = downloadHandler(
      filename = function() { 
        paste('graph1','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot22()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData5 = downloadHandler(
      filename = function() { 
        paste('graph1','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot22()
        dev.off()
      },
      contentType=NA)
  
  output$downloadData6 = downloadHandler(
    filename = function() { 
      paste('graph3','.jpg', sep='') 
    },
    content = function(file) {
      jpeg(file)
      Plot33()
      dev.off()
    },
    contentType='image/jpg')
  
  output$downloadData7 = downloadHandler(
    filename = function() { 
      paste('graph3','.png', sep='') 
    },
    content = function(file) {
      png(file)
      Plot33()
      dev.off()
    },
    contentType='image/png')
  
  output$downloadData8 = downloadHandler(
    filename = function() { 
      paste('graph3','.pdf', sep='') 
    },
    content = function(file) {
      pdf(file)
      Plot33()
      dev.off()
    },
    contentType=NA)
  
    
    Plot11=function(){
      if(input$select0=="cos2"){
        if(input$slider00!=1){
          selecindiv=paste("cos2 ",input$slider00)
        }
        else{
          selecindiv="cos2 0.999"
        }
        selecindivText=paste("'",selecindiv,"'",sep="")
      }
      if(input$select0=="NONE"){
        selecindiv=NULL
        selecindivText="NULL"
      }
      if(input$select0=="contrib"){
        selecindiv=paste("contrib ",input$slider4)
        selecindivText=paste("'",selecindiv,"'",sep="")
      }
      plot.FAMD(values()$res.FAMD,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="var",cex=input$cex2,cex.main=input$cex2,cex.axis=input$cex2,title=input$title2,select=selecindiv)
    }
    Plot22=function(){
      hab="none"
      if(length(QualiChoice)==0){
        hab="none"
      }
      else if(length(QualiChoice)==1 && input$habi==TRUE){
        if(input$habi==TRUE){
          hab=QualiChoice
          colquali="blue"
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      else if (length(QualiChoice)>1){
        if(input$habi==TRUE){
          if(is.null(input$habiller)){
            hab="none"
          }
          else{
            hab=input$habiller
          }
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      if(hab!="none"){
        hab=which(all==hab)
        hab=as.numeric(hab)
      }
      if(input$select=="cos2"){
        if(input$slider1!=1){
          selecindiv=paste("cos2 ",input$slider1)
        }
        else{
          selecindiv="cos2 0.999"
        }
        selecindivtext=paste0("'",selecindiv,"'")
      }
      if(input$select=="NONE"){
        selecindiv=NULL
      }
      if(input$select=="contrib"){
        selecindiv=paste("contrib ",input$slider0)
      }
      if(input$select=="Manuel"){
        selecindiv=c(input$indiv)
      }
      plot.FAMD(values()$res.FAMD,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="ind",select=selecindiv,lab.var=input$labels,lab.ind=input$labels2,cex=input$cex,cex.main=input$cex,cex.axis=input$cex,habillage=hab,title=input$title1)
    }
    Plot33=function(){
      if(input$selecti=="cos2"){
      if(input$slider000!=1){
        selecindiv=paste("cos2 ",input$slider000)
      }
      else{
        selecindiv="cos2 0.999"
      }
      selecindivText=paste("'",selecindiv,"'",sep="")
    }
    if(input$selecti=="NONE"){
      selecindiv=NULL
      selecindivText="NULL"
    }
    if(input$selecti=="contrib"){
      selecindiv=paste("contrib ",input$slider6)
      paste("'",selecindiv,"'",sep="")
    }
    plot.FAMD(values()$res.FAMD,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="quanti",cex=input$cex3,cex.main=input$cex3,cex.axis=input$cex3,title=input$title3,select=selecindiv)
    }
  }
)
      

