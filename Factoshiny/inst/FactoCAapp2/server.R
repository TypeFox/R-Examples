# server scipt for CA2

shinyServer(
  function(input, output) {
    values=reactive({
      if (input$selecactive=="Toutes"){
        data.selec=newdata[,VariableChoices]
      }
      else{
        validate(
          need(length(getactive()!=0), "Please select at least one supplementary column")
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
          colnames(data.selec)[dim(data.selec)[2]]=input$supvar
          choixquanti=length(data.selec)
        }
        else{
          choixquanti=seq((dim(data.selec)[2]-length(input$supvar)+1),dim(data.selec)[2])
        }
      }
      if(length(input$rowsupl)!=0){
        indexes=c()
        for (i in 1:length(nom)){
          if(nom[i]%in%input$rowsupl){
            indexes=c(indexes,i)
          }
        }
      }
      else{
        indexes=NULL
      }
      indexes=c(indexes,rowna)
      choixquanti2=NULL
      if(length(withna)!=0){
        data.selec=cbind(data.selec,newdata[,withna])
        if(length(withna)==1){
          colnames(data.selec)[dim(data.selec)[2]]=withna
          if(is.null(choixquanti)){
            choixquanti2=length(data.selec)
          }
          else{
            choixquanti2=c(choixquanti,length(data.selec))
          }
        }
        else{
          if(is.null(choixquanti)){
            choixquanti2=seq((dim(data.selec)[2]-length(withna)+1),dim(data.selec)[2])
          }
          else{
            choixquanti2=c(choixquanti,seq((dim(data.selec)[2]-length(withna)+1),dim(data.selec)[2]))
          }
        }
      }
      else{
        choixquanti2=choixquanti
      }
      list(res.CA=(CA(data.selec,quali.sup=choixquali,col.sup=choixquanti2,row.sup=indexes,graph=FALSE,ncp=5)),DATA=(data.selec),CHOIXQUALI=(choixquali),CHOIXQUANTI=(choixquanti2),INDEXES=(indexes))
    })
    
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
      res$data=newdata
      res$nomData=nomData
      # a : colonnes supplementaires
      res$a=input$supvar
      # b : lignes supplementaires
      res$b=input$rowsupl
      # c : colonnes quali
      choixquali=NULL
      if (length(QualiChoice)==1){
        if(input$supquali==TRUE){
          choixquali=QualiChoice
        }
      }
      else{
        if(length(input$supquali)!=0){
          choixquali=input$supquali
        }
      }
      res$c=choixquali
      # d et e : axes
      res$d=input$nb1
      res$e=input$nb2
      # f : invisible points 
      invisi=NULL
      if(length(input$invis)!=0){
        invisi=input$invis
      }
      res$f=invisi
      res$type1=input$seleccol
      res$type2=input$selecrow
      res$selec1=NULL
      if(input$seleccol=="cos2"){
        res$selec1=input$slider3
      }
      if(input$seleccol=="contrib"){
        res$selec1=input$contrib1
      }
      res$selec2=NULL
      if(input$selecrow=="cos2"){
        res$selec2=input$slider4
      }
      if(input$seleccol=="contrib"){
        res$selec2=input$contrib2
      }
      res$taille=input$cex
      res$code1=Code()
      res$code2=CodeGraph()
      res$title1=input$title1
      res$anafact=values()$res.CA
      class(res) <- "CAshiny"
      return(res)
    }
    
    
    observe({
      if(input$CAcode==0){
      }
      else {
        isolate({
          cat(Code(),sep="\n")
          cat(CodeGraph(),sep="\n")
        })
      }
    })
    
    
    createVec=function(arg){
      vec<-NULL
      vec<-paste(vec,arg[1],sep="")
      for (i in 2:(length(arg))){
        vec<-paste(vec,arg[i],sep=",")
      }
      vec<-paste("c(",vec,")",sep="")
      return(vec)
    }
    
    
    Code=function(){
      
      vecquant<-values()$CHOIXQUANTI
      vecqual<-values()$CHOIXQUALI
      Datasel<-values()$DATA
      indexes<-values()$INDEXES
      
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
      
      vecquant1<-createVec(vecquant)
      vecquant2<-vecquant
      
      vecqual1<-createVec(vecqual)
      vecqual2<-vecqual
      
      indexes1<-createVec(indexes)
      indexes2<-indexes
      
      if(length(vecqual)==0){
        vecqual<-"NULL" 
      }
      else if(length(vecqual)==1){
        vecqual<-vecqual
      }
      else if(length(vecqual)>1){
        vecqual<-vecqual2
      }
      
      if(length(vecquant)==0){
        vecquant<-"NULL"  
      }
      else if(length(vecquant)==1){
        vecquant
      }
      else if(length(vecquant)>1){
        vecquant<-vecquant1
      }
      
      
      if(length(indexes)==0){
        indexes<-"NULL"  
      }
      else if(length(indexes)==1){
        indexes
      }
      else if(length(indexes)>1){
        indexes<-indexes1
      }
      Call1=as.name(paste("res.CA=(CA(",vecfinal,",quali.sup=",vecqual,",col.sup=",vecquant,",row.sup=",indexes,",graph=FALSE,ncp=5))",sep=""))
      return(Call1)
    }
    
    CodeGraph=function(){
      sel="NULL"
      if(input$seleccol=="cos2"){
        if(input$slider3!=1){
          sel=paste("cos2 ",input$slider3)
        }
        else{
          sel="cos2 0.999"
        }
      }
      if(input$seleccol=="contrib"){
        sel=paste("contrib ",input$contrib1)
      }
      sel2="NULL"
      if(input$selecrow=="cos2"){
        if(input$slider4!=1){
          sel2=paste("cos2 ",input$slider4)
        }
        else{
          sel2="cos2 0.999"
        }
      }
      if(input$selecrow=="contrib"){
        sel2=paste("contrib ",input$contrib2)
      }
      Call2=paste("plot.CA(res.CA,axes=c(",as.numeric(input$nb1),",",as.numeric(input$nb2),"),selectCol='",sel,"',selectRow='",sel2,"',unselect=0,col.sup='darkred',cex=",input$cex,",title='",input$title1,"',invisible='",Plot1()$invisiText,"')",sep="")
      return(Call2)
    }
    
    Plot1=reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      if(length(input$invis)==0){
        invisi="none"
        invisiText="'non'"
      }
      if(length(input$invis)!=0){
        invisi=input$invis
        invisiText=invisi
      }
      sel=NULL
      if(input$seleccol=="cos2"){
        if(input$slider3!=1){
          sel=paste("cos2 ",input$slider3)
        }
        else{
          sel="cos2 0.999"
        }
      }
      if(input$seleccol=="contrib"){
        sel=paste("contrib ",input$contrib1)
      }
      sel2=NULL
      if(input$selecrow=="cos2"){
        if(input$slider4!=1){
          sel2=paste("cos2 ",input$slider4)
        }
        else{
          sel2="cos2 0.999"
        }
      }
      if(input$selecrow=="contrib"){
        sel2=paste("contrib ",input$contrib2)
      }
      list(PLOT1=(plot.CA(values()$res.CA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),selectCol=sel,title=input$title1,selectRow=sel2,cex=input$cex,cex.main=input$cex,cex.axis=input$cex,unselect=0,col.col.sup="darkred",invisible=invisi)),invisiText=(invisiText))
    })
    
    output$map <- renderPlot({
      p <- Plot1()$PLOT1
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
    
    output$contribcol=renderUI({
      maxx=dim(values()$res.CA$col$coord)[1]
      if(selec1=="contrib"){
        return(sliderInput("contrib1",h6("Number of the most contributive active columns"),min=1,max=maxx,value=valueselec2,step=1))
      }
      else{
        return(sliderInput("contrib1",h6("Number of the most contributive active columns"),min=1,max=maxx,value=maxx,step=1))
      }
      
    })
    
    output$contribrow=renderUI({
      maxx=dim(values()$res.CA$row$coord)[1]
      if(selec2=="contrib"){
        return(sliderInput("contrib2",h6("Number of the most contributive active rows"),min=1,max=maxx,value=valueselec2,step=1))
      }
      else{
        return(sliderInput("contrib2",h6("Number of the most contributive active rows"),min=1,max=maxx,value=maxx,step=1))
      }
    })
    
    output$out22=renderUI({
      choix=list("Summary of CA"="CA","Eigenvalues"="eig","Coordinates for the columns"="var","Coordinates of the rows"="ind")
      if(!is.null(values()$INDEXES)){
        choix=c(choix,"Coordinates of the supplementary rows"="suprow")
      }
      if(!is.null(values()$CHOIXQUANTI)){
        choix=c(choix,"Coordinates of the supplementary columns"="supcol")
      }
      if(!is.null(values()$CHOIXQUALI)){
        choix=c(choix,"Coordinates of the categorical variables"="qualico")
      }
      radioButtons("out","Which outputs do you want ?",
                   choices=choix,selected="CA",inline=TRUE)
    })
    
    output$warn=renderPrint({
      if(length(withna)!=0){
        baba=paste(withna,collapse=", ")
        bibi=paste(nomrow,collapse=", ")
        a=paste0("Warning : ", baba, " have NA : they are considered as supplementary columns")
        b=paste0("Warning : ", bibi, " have NA : they are considered as supplementary rows")
        return(cat(a,b,sep="\n"))
      }
    })
    
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
      return(as.data.frame(values()$res.CA$eig))
    })
    
    output$sorties1=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Not enough active columns")
      )
      return(as.data.frame(values()$res.CA$col$coord))
    })
    
    output$sorties2=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Not enough active columns")
      )
      return(as.data.frame(values()$res.CA$col$cos2))
    })
    
    output$sorties3=renderTable({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Not enough active columns")
      )
      return(as.data.frame(values()$res.CA$col$contrib))
    })
    
    output$sorties4=renderTable({
      return(as.data.frame(values()$res.CA$row$coord))
    })
    
    output$sorties5=renderTable({
      return(as.data.frame(values()$res.CA$row$cos2))
    })
    
    output$sorties6=renderTable({
      return(as.data.frame(values()$res.CA$row$contrib))
    })
    
    output$sorties7=renderTable({
      validate(
        need((length(input$rowsupl)>0), "No supplementary rows selected")
      )
      return(as.data.frame(values()$res.CA$row.sup$coord))
    })
    
    output$sorties8=renderTable({
      return(as.data.frame(values()$res.CA$row.sup$cos2))
    })
    
    output$sorties9=renderTable({
      return(as.data.frame(values()$res.CA$col.sup$coord))
    })
    
    output$sorties10=renderTable({
      return(as.data.frame(values()$res.CA$col.sup$cos2))
    })
    
    output$sorties11=renderTable({
      validate(
        need((length(input$supquali)>0 || input$supquali==TRUE), "No categorical variables selected")
      )
      return(as.data.frame(values()$res.CA$quali.sup))
    })
    
    
    
    output$map3=renderPlot({
      return(barplot(values()$res.CA$eig[,1],names.arg=rownames(values()$res.CA$eig),las=2))
    })
    
    ### Fonction permettant l'affichage du JDD sous la forme d'un DataTable, qui permet la recherche de donnes. 
    output$JDD=renderDataTable({
      cbind(Names=rownames(newdata),newdata)},
      options = list(    "orderClasses" = TRUE,
                         "responsive" = TRUE,
                         "pageLength" = 10))
    
    ### Fonction permettant l'affichage du summary du JDD
    output$summary=renderPrint({
      summary(newdata)
    })
    
    
    ### Fonction permettant l'affichage du summary de la fonction CA sur le JDD
    output$summaryCA=renderPrint({
      a<-values()$res.CA
      a$call$call<-Code()
      summary.CA(a,nbelements=input$nbele)
    })
    
    output$summary2=downloadHandler(filename = function() { 
      paste('summaryofCA','.txt', sep='') 
    },
    content = function(file) {
      summary.CA(values()$res.CA,nbelements=input$nbele,file=file)
    },
    contentType='text/csv')
    
    
    ## Creation des fonctions permettant l'enregistrement des graphs sous les formats : png, jpeg, pdf et emf
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
    
    Plot11=function(){
      if(length(input$invis)==0){
        invisi="none"
      }
      if(length(input$invis)!=0){
        invisi=input$invis
      }
      sel=NULL
      if(input$seleccol=="cos2"){
        if(input$slider3!=1){
          sel=paste("cos2 ",input$slider3)
        }
        else{
          sel="cos2 0.999"
        }
      }
      if(input$seleccol=="contrib"){
        sel=paste("contrib ",input$contrib1)
      }
      sel2=NULL
      if(input$selecrow=="cos2"){
        if(input$slider4!=1){
          sel2=paste("cos2 ",input$slider4)
        }
        else{
          sel2="cos2 0.999"
        }
      }
      if(input$seleccol=="contrib"){
        sel2=paste("contrib ",input$contrib2)
      }
      plot.CA(values()$res.CA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),selectCol=sel,selectRow=sel2,cex=input$cex,cex.main=input$cex,cex.axis=input$cex,title=input$title1,unselect=0,col.col.sup="darkred",invisible=invisi)
    }
    
  }
)
