# server script for PCA2
shinyServer(
  function(input, output) {
    values=reactive({
    if (input$selecactive=="Toutes"){
      data.selec=newdata[,VariableChoices]
    }
    else{
      validate(
        need(length(input$supvar)>0, "Please select at least one supplementary variables")
      )
      data.selec=newdata[,c(getactive())]
    }
    if(length(QualiChoice)==0){
      choixquali=NULL
    }
    else if (length(QualiChoice)==1){
      if(input$supquali==FALSE){
        choixquali=NULL
      }
      else{
        data.selec=cbind(data.selec,newdata[,QualiChoice])
        colnames(data.selec)[dim(data.selec)[2]]=QualiChoice
        choixquali=length(data.selec)
      }
    }
    else{
      if(length(input$supquali)==0){
        choixquali=NULL
      }
      else{
        data.selec=cbind(data.selec,newdata[,input$supquali])
        if(length(input$supquali)==1){
          choixquali=length(data.selec)
          colnames(data.selec)[choixquali]=input$supquali
        }
        else{
          choixquali=seq((dim(data.selec)[2]-length(input$supquali)+1),dim(data.selec)[2])
          colnames(data.selec)[choixquali]=input$supquali
        }
      }
    }
    if(length(input$supvar)==0){
      choixquanti=NULL
    }
    else {
      data.selec=cbind(data.selec,newdata[,input$supvar])
      if(length(input$supvar)==1){
        choixquanti=length(data.selec)
        colnames(data.selec)[choixquanti]<-input$supvar
      }
      else{
        choixquanti=seq((dim(data.selec)[2]-length(input$supvar)+1),dim(data.selec)[2])
      }
    }
    if (length(input$habiller)==2 && input$habi==TRUE){
      data.selec <- data.frame(data.selec,newCol=paste(newdata[,input$habiller[1]],newdata[,input$habiller[2]],sep="/"))
      choixquali=c(choixquali,dim(data.selec)[2])
    }
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
    list(res.PCA=(PCA(data.selec,quali.sup=choixquali,quanti.sup=choixquanti,scale.unit=TRUE,graph=FALSE,ncp=5,ind.sup=suple)),DATA=(data.selec),choixquant=(choixquanti),choixqual=(choixquali),choixsuple=(suple))
    })
    
    Plot1 <- reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
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
      list(PLOT=(plot.PCA(values()$res.PCA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="var",select=selecindiv,unselect=0,col.quanti.sup="blue",cex=input$cex2,cex.main=input$cex2,cex.axis=input$cex2,title=input$title2)),SELECTION=(selecindiv),selecindivText=(selecindivText))
    })
    
    output$map <- renderPlot({
      p <- Plot1()$PLOT
    })
    
    Plot2 <- reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      validate(
        need(input$habiller == TRUE || input$habiller == FALSE || length(input$habiller)<=2,"Please select maximum 2 variables as habillage")
      )
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
        selecindivtext="NULL"
      }
      if(input$select=="contrib"){
        selecindiv=paste("contrib ",input$slider0)
        selecindivtext=paste0("'",selecindiv,"'")
      }
      if(input$select=="Manuel"){
        selecindiv=c(input$indiv)
      }
      if(input$supquali==FALSE || length(QualiChoice)==0 || length(input$supquali)==0 || input$habi==FALSE){
        hab="none"
        colquali="magenta"
      }
      else if(length(QualiChoice)==1 && input$supquali==TRUE){
        if(input$habi==TRUE){
          hab=QualiChoice
          colquali="blue"
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      else if (length(input$supquali)==1){
        if(input$habi==TRUE){
          hab=input$supquali
          colquali="blue"
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      if(length(input$supquali)>1){
        if(length(input$habiller)==0){
          hab="none"
          colquali="magenta"
        }
        if (length(input$habiller)==1 & input$habi==TRUE){
          hab=as.character(input$habiller)
          colquali="blue"
        }
        if (length(input$habiller)==2 & input$habi==TRUE){
          hab=dim(values()$DATA)[2]
          colquali="blue"
        }
      }
      
      if(input$select=="Manuel"){
        if(length(input$indiv)==0){
          selecindivtext="NULL"
        }
        if(length(input$indiv)>1){
          vec<-NULL
          vec<-paste(vec,"'",selecindiv[1],"'",sep="")
          
          for (i in 2:(length(selecindiv))){
            vec<-paste(vec,paste("'",selecindiv[i],"'",sep=""),sep=",")
          }
          selecindivtext<-paste("c(",vec,")",sep="")
        }
        else if (length(input$indiv)==1){
          selecindivtext=paste0("'",c(input$indiv),"'")
        }
      }
      list(PLOT=(plot.PCA(values()$res.PCA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="ind",cex=input$cex,cex.main=input$cex,cex.axis=input$cex,select=selecindiv,habillage=hab,col.quali=colquali,col.ind.sup="blue",title=input$title1)),SELECTION2=(selecindiv),SELECTION3=(selecindivtext),HABILLAGE=(hab),colquali=(colquali))
      
    })
   
    output$map2 <- renderPlot({
      p <- Plot2()$PLOT
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
      res$data=newdata
      res$a=values()$DATA
      if (length(QualiChoice)==1){
        if(input$supquali==FALSE){
          quali=NULL
        }
        else{
          quali=QualiChoice
        }
      }
      else{
        if(length(input$supquali)==0){
          quali=NULL
        }
        else{
          quali=input$supquali
        }
      }
      res$b=quali
      res$c=input$supvar
      res$d=input$indsup
      res$e=input$nb1
      res$f=input$nb2
      hab=NULL
      if(length(QualiChoice)==1 && input$supquali==TRUE){
        if(input$habi==TRUE){
          hab=QualiChoice
        }
      }
      else if (length(input$supquali)==1){
        if(input$habi==TRUE){
          hab=input$supquali
        }
      }
      if(length(input$supquali)>1){
        if (length(input$habiller)==1 & input$habi==TRUE){
          hab=as.character(input$habiller)
        }
        if (length(input$habiller)==2 & input$habi==TRUE){
          hab=input$habiller
        }
      }
      res$g=hab
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
      res$l=input$cex
      res$m=input$cex2
      res$code1=code()
      res$code2=codeGraphVar()
      res$code3=codeGraphInd()
      res$title1=input$title1
      res$title2=input$title2
      res$anafact=values()$res.PCA
      class(res) <- "PCAshiny"
      return(res)
    }
    
    #### Fonction recuperation de code
    
    observe({
      if(input$PCAcode==0){
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
  
    
    code<-function(){
      vecquant<-values()$choixquant
      choixqual<-values()$choixqual
      Datasel<-values()$DATA
      indsupl<-values()$choixsuple
      
      vec<-NULL
      for (i in 1:length(colnames(Datasel))){
        vec<-c(vec,colnames(Datasel)[i])
      }
      vec2<-NULL
      vec2<-paste(vec2,"'",vec[1],"'",sep="")
      for (i in 2:(length(vec))){
        vec2<-paste(vec2,paste("'",vec[i],"'",sep=""),sep=",")
      }
      vecfinal<-paste(nomData,"[,c(",vec2,")","]",sep="")
      
      vec4<-NULL
      vec4<-paste(vec4,vecquant[1],sep="")
      for (i in 2:(length(vecquant))){
        vec4<-paste(vec4,vecquant[i],sep=",")
      }
      vecquant1<-paste("c(",vec4,")",sep="")
      vecquant2<-vecquant
      
      vecqual<-choixqual
      vec5<-NULL
      vec5<-paste(vec5,vecqual[1],sep="")
      for (i in 2:(length(vecqual))){
        vec5<-paste(vec5,vecqual[i],sep=",")
      }
      vecqual1<-paste("c(",vec5,")",sep="")
      vecqual2<-vecqual
      
      vecind<-NULL
      vecind<-paste(vecind,indsupl[1],sep="")
      for (i in 2:(length(indsupl))){
        vecind<-paste(vecind,indsupl[i],sep=",")
      }
      vecind1<-paste("c(",vecind,")",sep="")
      vecind2<-indsupl
      vec<-vecfinal
      
      if(length(input$indsup)==0){
        indsupl<-"NULL"
      }
      else if(length(input$indsup)==1){
        indsupl<-vecind2
      }
      else if(length(input$indsup)>1){
        indsupl<-vecind1
      }
      
      
      if(length(input$supvar)>1){
        vecquant<-vecquant1
      }
      else if(length(input$supvar)==1){
        vecquant<-vecquant2
      }
      else if(length(input$supvar)==0){
        vecquant<-"NULL"
      }
      
      if (length(input$supquali)>1){ 
        vecqual<-vecqual1
      }
      if(length(QualiChoice)==1){
        if(input$supquali==TRUE){
          vecqual<-vecqual2 
        }
        else{
          vecqual<-"NULL"  
        }
      }
      
      else if(length(QualiChoice)>1){
        if(length(input$supquali)==1){
          vecqual<-vecqual2  
        }
        else if (length(input$supquali)>1){ 
          vecqual<-vecqual1
        }
        else if (length(input$supquali)==0){ 
          vecqual<-"NULL"
        }  
      }
      else if(length(QualiChoice)==0){
        vecqual<-"NULL"
      }
      Call1=as.name(paste("res.PCA<-PCA(",vec,",quali.sup=",vecqual,",","quanti.sup=",vecquant,",ind.sup=",indsupl,",scale.unit=TRUE,graph=FALSE,ncp=5)",sep=""))
      return(Call1)
    }
    
    
    codeGraphVar<-function(){
      
      if(length(input$slider4)==0){
        selection="NULL"
      }
      else{
        selection=Plot1()$selecindivText
      }
      Call1=paste("plot.PCA(res.PCA,axes=c(",input$nb1,",",input$nb2,"),choix='var',select=",selection,",cex=",input$cex2,",cex.main=",input$cex2,",cex.axis=",input$cex2,",title='",input$title2,"',unselect=0,col.quanti.sup='red')",sep="")
      return(Call1)
    }
    
    codeGraphInd<-function(){
      if (length(input$habiller)<=1 & input$habi==TRUE || input$habi==FALSE){
        hab=paste("'",Plot2()$HABILLAGE,"'",sep="")
      }
      else if (length(input$habiller)==2 & input$habi==TRUE){
        hab=Plot2()$HABILLAGE
      }
      
      Call2=paste("plot.PCA(res.PCA,","axes=c(",input$nb1,",",input$nb2,"),choix='ind',select=",Plot2()$SELECTION3,",habillage=",hab,",title='",input$title1,"',cex=",input$cex,",cex.main=",input$cex,",cex.axis=",input$cex,",col.quali='",Plot2()$colquali,"',col.ind.sup='blue')",sep="")
      return(Call2)
    }
    
    ##### Fin de la fonction recuperation du code
    
    
    output$out22=renderUI({
      choix=list("Summary of PCA"="ACP","Eigenvalues"="eig","Results of the variables"="resvar","Results of the individuals"="resind")
      if(!is.null(values()$choixsuple)){
        choix=c(choix,"Results of the supplementary individuals"="supind")
      }
      if(!is.null(values()$choixquant)){
        choix=c(choix,"Results of the supplementary variables"="varsup")
      }
      if(!is.null(values()$choixqual)){
        choix=c(choix,"Results of the categorical variables"="qualico")
      }
      radioButtons("out","Which outputs do you want ?",
                   choices=choix,selected="ACP",inline=TRUE)
    })
    
    getactive=function(){
      if(input$selecactive=="choix"){
      sup=c()
      if(length(input$supvar)==0){
        activevar=VariableChoices
      }
      else{
        for (i in 1:length(VariableChoices)){
          if(VariableChoices[i]%in%input$supvar){
            sup=c(sup,i)
          }
        }
        activevar=VariableChoices[-sup]
      }
      return(activevar)
    }
  }
    
    
    output$NB1=renderUI({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes" || length(getactive())>5){
        return(selectInput("nb1", label = h6("x axis"), 
                    choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = axe1,width='80%'))
      }
      else{
        baba=c(1:length(getactive()))
        return(selectInput("nb1",label=h6("x axis"), choices=baba,selected=axe1,width='80%'))
      }
    })
    
    output$NB2=renderUI({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes" || length(getactive())>5){
        return(selectInput("nb2", label = h6("y axis"), 
                           choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected = axe2,width='80%'))
      }
      else{
        baba=c(1:length(getactive()))
        return(selectInput("nb2",label=h6("y axis"), choices=baba,selected=axe2,width='80%'))
      }
    })
    
    output$sorties=renderTable({
        return(as.data.frame(values()$res.PCA$eig))
    })
    
    output$sorties12=renderTable({
        validate(
          need((length(input$supquali)>0 || input$supquali==TRUE), "No categorical variables selected")
        )
        validate(
          need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
        )
        return(as.data.frame(values()$res.PCA$quali.sup$coord))
    })
    
    output$sorties13=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      validate(
        need((length(input$supquali)>0 || input$supquali==TRUE), "No categorical variables selected")
      )
      return(as.data.frame(values()$res.PCA$quali.sup$v.test))
    })
    
    output$sorties2=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$var$coord))
    })
    
    output$sorties22=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$ind$coord))
    })
    
    output$sorties23=renderTable({
      validate(
        need(length(input$supvar)!=0, "No supplementary quantitative variables")
      )
      validate(
        need(length(input$getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$quanti.sup$coord))
    })
    
    output$sorties32=renderTable({
      validate(
        need(length(input$supvar)!=0, "No supplementary quantitative variables")
      )
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$quanti.sup$cor))
    })
    
    output$sorties36=renderTable({
      validate(
        need(length(input$indsup)!=0, "No supplementary individuals")
      )
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$ind.sup$coord))
    })
    
    output$sorties37=renderTable({
      validate(
        need(length(input$indsup)!=0, "No supplementary individuals")
      )
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$ind.sup$cos2))
    })
    
    
    output$sorties3=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$var$contrib))
    })
    
    output$sorties33=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$ind$contrib))
    })
    
    output$sorties4=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$var$cos2))
    })
    
    output$sorties44=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      return(as.data.frame(values()$res.PCA$ind$cos2))
    })
  
  output$sortieDimdesc3=renderTable({
    validate(
      need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.PCA)[[1]]$quanti))
  })
  
  output$sortieDimdesc4=renderTable({
    validate(
      need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.PCA)[[1]]$quali))
  })
  
  #DIM2
  
  output$sortieDimdesc33=renderTable({
    validate(
      need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.PCA)[[2]]$quanti))
  })
  
  output$sortieDimdesc44=renderTable({
    validate(
      need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.PCA)[[2]]$quali))
  })
  
  #DIM3
  
  output$sortieDimdesc333=renderTable({
    validate(
      need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.PCA)[[3]]$quanti))
  })
  
  output$sortieDimdesc444=renderTable({
    validate(
      need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
    )
    return(as.data.frame(dimdesc(values()$res.PCA)[[3]]$quali))
  })
    
    
    output$map3=renderPlot({
      return(barplot(values()$res.PCA$eig[,1],names.arg=rownames(values()$res.PCA$eig),las=2))
    })
    
    output$JDD=renderDataTable({
      cbind(Names=rownames(newdata),newdata)},
      options = list(    "orderClasses" = TRUE,
                         "responsive" = TRUE,
                         "pageLength" = 10))
  
    output$summary=renderPrint({
      summary(newdata)
    })
  
    output$summaryPCA=renderPrint({
      validate(
        need(input$nbele!=0, "Please select at least one element")
      )
      a<-values()$res.PCA
      a$call$call<-code()
      summary.PCA(a,nbelements=input$nbele)
    })
  
    output$summary2=downloadHandler(filename = function() { 
      paste('summaryofPCA','.txt', sep='') 
    },
    content = function(file) {
      summary.PCA(values()$res.PCA,nbelements=input$nbele,file=file)
    },
    contentType='text/csv')
  
    
    output$slider3=renderUI({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes"){
        maxvar=length(VariableChoices)
      }
      if(input$selecactive=="choix"){
        maxvar=length(getactive())
      }
      if(selection3=="contrib"){
        return(div(align="center",sliderInput("slider4",label="Number of the most contributive variables",
                                              min=1,max=maxvar,value=selection4,step=1)))  
      }
      else{
      return(div(align="center",sliderInput("slider4",label="Number of the most contributive variables",
                  min=1,max=maxvar,value=maxvar,step=1)))}
    })

    
    output$habillage2=renderUI({
      if(length(QualiChoice)==0 || input$supquali==FALSE || length(input$supquali)==0){
        return(p("No categorical variable"))
      }
      if(length(input$supquali)>1){
        if(is.null(habillageind)){
        num=c(1:length(input$supquali))
        return(selectInput("habiller","Select 1 or 2 variables", choices=list(num=input$supquali),multiple=TRUE))
        }
        else{
          num=c(1:length(input$supquali))
          return(selectInput("habiller","Select 1 or 2 variables", choices=list(num=input$supquali),multiple=TRUE,selected=habillageind))
        }
      }
    })
      
    output$histo=renderPlot({
      par(mfrow=c(1,2))
      boxplot(newdata[,input$bam])
      hist(newdata[,input$bam],main="",xlab="")
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
    
    Plot11=function(){
      if(input$select0=="cos2"){
        if(input$slider00!=1){
          selecindiv=paste("cos2 ",input$slider00)
        }
        else{
          selecindiv="cos2 0.999"
        }
      }
      if(input$select0=="NONE"){
        selecindiv=NULL
      }
      if(input$select0=="contrib"){
        selecindiv=paste("contrib ",input$slider4)
      }
      plot.PCA(values()$res.PCA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="var",select=selecindiv,unselect=0,col.quanti.sup="blue",cex=input$cex2,cex.main=input$cex2,cex.axis=input$cex2,title=input$title2)
    }
    Plot22=function(){
      if(input$select=="cos2"){
        if(input$slider1!=1){
          selecindiv=paste("cos2 ",input$slider1)
        }
        else{
          selecindiv="cos2 0.999"
        }
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
      if(input$supquali==FALSE || length(QualiChoice)==0 || length(input$supquali)==0 || input$habi==FALSE){
        hab="none"
        colquali="magenta"
      }
      else if(length(QualiChoice)==1 && input$supquali==TRUE){
        if(input$habi==TRUE){
          hab=QualiChoice
          colquali="blue"
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      else if (length(input$supquali)==1){
        if(input$habi==TRUE){
          hab=input$supquali
          colquali="blue"
        }
        else{
          hab="none"
          colquali="magenta"
        }
      }
      if(length(input$supquali)>1){
        if(length(input$habiller)==0){
          hab="none"
          colquali="magenta"
        }
        if (length(input$habiller)==1 & input$habi==TRUE){
          hab=as.character(input$habiller)
          colquali="blue"
        }
        if (length(input$habiller)==2 & input$habi==TRUE){
          hab=dim(values()$DATA)[2]
          colquali="blue"
        }
      }
      plot.PCA(values()$res.PCA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="ind",cex=input$cex,cex.main=input$cex,cex.axis=input$cex,select=selecindiv,habillage=hab,col.quali=colquali,col.ind.sup="blue",title=input$title1)
    }
  }
)
      

