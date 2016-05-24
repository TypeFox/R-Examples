# server script for MCA2

shinyServer(
  function(input, output) {
    
    #Realisation de l'ACM    
    values=reactive({
      
      if (input$selecactive=="Toutes"){
        data.selec=newdata[,VariableChoices]
      }
      else{
        validate(
          need(getactive()!= "", "Please select active variables")
        )
        data.selec=newdata[,c(getactive())]
      }
      
      
      if(length(QuantiChoice)==0){
        choixquanti=NULL
      }
      else if (length(QuantiChoice)==1){
        if(input$supquanti==FALSE){
          choixquanti=NULL
        }
        else{
          data.selec=cbind(data.selec,newdata[,QuantiChoice])
          colnames(data.selec)[dim(data.selec)[2]]=QuantiChoice
          #Renomme les colonnes
          choixquanti=length(data.selec)
        }
      }
      #Si plusieurs quanti existent
      else{
        if(length(input$supquanti)==0){
          choixquanti=NULL
        }
        else{
          data.selec=cbind(data.selec,newdata[,input$supquanti])
          if(length(input$supquanti)==1){
            choixquanti=length(data.selec)
            colnames(data.selec)[choixquanti]=input$supquanti
          }
          else{
            choixquanti=seq((dim(data.selec)[2]-length(input$supquanti)+1),dim(data.selec)[2])
            colnames(data.selec)[choixquanti]=input$supquanti
          }
        }
      }
      if(length(input$supvar)==0){
        choixquali=NULL
      }
      else {
        data.selec=cbind(data.selec,newdata[,input$supvar])
        if(length(input$supvar)==1){
          choixquali=length(data.selec)
          #modif
          colnames(data.selec)[choixquali]=input$supvar
        }
        else{
          choixquali=seq((dim(data.selec)[2]-length(input$supvar)+1),dim(data.selec)[2])
        }
      }
      if (length(input$habiller)==2){
        data.selec <- data.frame(data.selec,newCol=paste(newdata[,input$habiller[1]],newdata[,input$habiller[2]],sep="/"))
        choixquali=c(choixquali,dim(data.selec)[2])
      }
      
      if(is.null(input$indsup)){
        indsuplem<-NULL
      }
      else {
        vec<-NULL
        for(i in 1:length(input$indsup)){
          vec<-c(vec,which(rownames(newdata)==input$indsup[i]))
        }
        indsuplem<-vec
      }
      list(res.MCA=(MCA(data.selec,quanti.sup=choixquanti,quali.sup=choixquali,ind.sup=indsuplem,graph=FALSE)),DATA=(data.selec),choixquant=(choixquanti),choixqual=(choixquali),indsup=(indsuplem))     
    })
    
    
    observe({
      if(input$MCAcode==0){
      }
      else {
        isolate({
          if (length(input$habiller)==2 & input$habi==TRUE){
            cat(paste("newCol=paste(",nomData,"['",input$habiller[1],"'],x[,'",input$habiller[2],"'],sep='/'))",sep=""),sep="\n")
          }
          cat(code(),sep="\n")
          cat(codeGraphVar(),sep="\n")
          cat(codeGraphInd(),sep="\n")
          
          if((length(values()$choixquant)!=0)){
            cat(codeGraphQuanti(),sep="\n") 
          }
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
    
    code<-function(){
      vecquant<-values()$choixquant
      choixqual<-values()$choixqual
      Datasel<-values()$DATA
      indsup<-values()$indsup
      
      
      vec<-NULL
      for (i in 1:length(colnames(Datasel))){
        vec<-c(vec,colnames(Datasel)[i])
      }
      vec<-paste0("'",vec,"'")
      vec<-createVec(vec)
      
      vecquant1<-createVec(vecquant)
      vecquant2<-vecquant
      
      vecqual<-choixqual
      vecqual1<-createVec(vecqual)
      vecqual2<-vecqual
      
      
      indsup1<-createVec(indsup)
      indsup2<-indsup
      
      if(length(input$supvar)>1){
        vecqual<-vecqual1
      }
      else if(length(input$supvar)==1){
        vecqual<-vecqual2
      }
      else if(length(input$supvar)==0){
        vecqual<-"NULL"
      }
      
      if(length(input$indsup)==0){
        indsuplem<-"NULL"
      }
      else if(length(input$indsup)==1){
        indsuplem<-indsup2
      }
      else if(length(input$indsup)>1){
        indsuplem<-indsup1
      }
      
      if(length(QuantiChoice)==0){
        vecquant<-"NULL"
      }
      
      else if(length(QuantiChoice)==1){
        if(input$supquanti==TRUE){
          vecquant<-vecquant2 
        }
        else{
          vecquant<-"NULL"  
        }
      }
      
      else if(length(QuantiChoice)>1){
        if(length(input$supquanti)==1){
          vecquant<-vecquant2  
        }
        else if (length(input$supquanti)>1){ 
          vecquant<-vecquant1
        }
        else if (length(input$supquanti)==0){ 
          vecquant<-"NULL"
        }  
      }
      Call1=as.name(paste("res.MCA<-MCA(",nomData,"[,",vec,"],quali.sup=",vecqual,",","quanti.sup=",vecquant,",ind.sup=",indsuplem,",graph=FALSE)",sep=""))  
      return(Call1)
    }
    
    
    codeGraphVar<-function(){
      Call2=paste("plot.MCA(res.MCA,choix='var',invisible=",Plot4()$invisible,",title='",input$title2,"',axes=c(",as.numeric(input$nb1),",",as.numeric(input$nb2),"))",sep="")  
      return(Call2)
    }
    
    codeGraphInd<-function(){
      Call3=paste("plot.MCA(res.MCA,choix='ind',invisible=",Plot1()$inv,",axes=c(",as.numeric(input$nb1),",",as.numeric(input$nb2),"),selectMod=",Plot1()$selm,",selec=",Plot1()$sel,",habillage=",Plot1()$hab,",title='",input$title1,"',col.quali='",Plot1()$colquali,"',col.ind.sup='",Plot1()$colindsup,"')",sep="")   
      return(Call3)
    }
    
    codeGraphQuanti<-function(){
      Call4=paste("plot.MCA(res.MCA,axes=c(",as.numeric(input$nb1),",",as.numeric(input$nb2),"),choix='quanti.sup',title='",input$title3,"')",sep="")
      return(Call4)
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
      res$nomData=nomData
      res$data=newdata
      res$a=values()$DATA#data of the factorial analysis
      class(res)<-"MCAshiny"#Class of the result
      
      #Supplementary quantitative variables selected
      if (length(QuantiChoice)==1){
        if(input$supquanti==FALSE){
          quanti=NULL
        }
        else{
          quanti=QuantiChoice
        }
      }
      else{
        if(length(input$supquanti)==0){
          quanti=NULL
        }
        else{
          quanti=input$supquanti
        }
      }
      res$b=quanti
      
      res$c=input$supvar#suplementary qualitative variables
      res$z=input$var_sup#1st graph multiple choice selected
      res$y=input$ind_var#2nd graph multiple choice selected
      res$lab=input$indvarpoint
      res$d=input$indsup#supplementary individuals selected
      
      res$e=input$nb1#axes selected
      res$f=input$nb2#
      
      #Selected habillage
      if(length(input$supvar)==0 || input$habi==FALSE){
        hab="none"
      }
      
      if(length(input$supvar)>1){
        if(length(input$habiller)==0){
          hab="none"
        }
        
        if (length(input$habiller)==1 & input$habi==TRUE){
          hab=as.character(input$habiller)
        }
        
        if (length(input$habiller)==2 & input$habi==TRUE){
          hab=dim(values()$DATA)[2]
        }
      }
      else if (length(input$supvar)==1){
        if(input$habi==TRUE){
          hab=values()$choixqual
        }
        else{
          hab="none"
        }
      } 
     
      res$g=hab
      
      #Selection for individuals
      if(input$select=="Manuel"){
        selecindiv=input$indiv 
      }
      else if(input$select=="cos2"){
        selecindiv=input$slider1
      }
      else if(input$select=="Contrib"){
        selecindiv=input$sliderContrib 
        }
      else if(input$select=="NONE"){
        selecindiv=NULL
      }
      res$h=input$select#Type of selections
      res$i=selecindiv#selection
    
    #Selection for modalities
    if(input$selectMod=="cos2"){
      selecMod=input$sliderCosMod
    }
    else if(input$selectMod=="Contrib"){
      selecMod=input$slider4
    }
    else if(input$selectMod=="NONE"){
      selecMod=NULL
    }
    res$j=input$selectMod
    res$k=selecMod
    res$code1=code()
    res$code2=codeGraphVar()
    res$code3=codeGraphInd()
    if((length(values()$choixquant)!=0)){
      res$code4=codeGraphQuanti() 
    }
    else{
      res$code4=NULL
    }
    res$title1=input$title1
    res$title2=input$title2
    res$title3=input$title3
    res$anafact=values()$res.MCA
    return(res)
    }
    
    #Getactive
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
    
    output$choixindvar=renderUI({
      choix=list("Individuals"="Ind","Modalities"="Mod")
      selec=list("Ind","Mod")
      if(!(is.null(input$indsup))){
        choix=c(choix,"Supplementary individuals"="Indsup")
        selec=c(selec,"Indsup")
      }
      if(!(is.null(input$supvar))){
        choix=c(choix,"Supplementary modalities"="Modsup")
        selec=c(selec,"Modsup")
      }
      div(align="center",checkboxGroupInput("ind_var","", choices=choix,
                                                   selected = indvar))
    })
    
    output$pointlabel=renderUI({
      validate(
        need(!is.null(input$ind_var),""))
      choix=list()
      selec=c()
      reponse=input$ind_var
      if("Ind"%in% reponse){
        choix=c(choix,"Individuals"="Ind")
        selec=c(selec,"Ind")
      }
      if("Mod" %in% reponse){
        choix=c(choix,"Modalities"="Mod")
        selec=c(selec,"Mod")
      }
      if("Indsup" %in% reponse){
        choix=c(choix,"Supplementary individuals"="Indsup")
        selec=c(selec,"Indusp")
      }
      if("Modsup"%in% reponse){
        choix=c(choix,"Supplementary modalities"="Modsup")
        selec=c(selec,"Modsup")
      }
      div(align="center",checkboxGroupInput("indvarpoint","",choices=choix,selected=labvar))
    })
    
    output$out22=renderUI({
      choix=list("Summary of MCA"="MCA","Eigenvalues"="eig","Results of the variables"="resvar","Results of the individuals"="resind")
      if(!is.null(values()$indsup)){
        choix=c(choix,"Results of the supplementary individuals"="Isup")
      }
      if(!is.null(values()$choixquant)){
        choix=c(choix,"Results of the supplementary quantitative variables"="quantico")
      }
      if(!is.null(values()$choixqual)){
        choix=c(choix,"Results of the supplementary categorical variables"="varsup")
      }
      radioButtons("out","Which outputs do you want ?",
                   choices=choix,selected="MCA",inline=TRUE)
    })
    
    
    #Getinv
    getinv=function(){
      
      inv<-c()
      if(!("Ind"%in%input$ind_var)){
        inv<-c(inv,"ind")
      }
      
      if(!("Mod"%in%input$ind_var)){
        inv<-c(inv,"var")
      }
      if(!(is.null(values()$choixqual))){
      if(!("Modsup"%in%input$ind_var)){
        inv<-c(inv,"quali.sup")
      }
      }
      if(!(is.null(values()$indsup))){
      if(!("Indsup"%in%input$ind_var)){
        inv<-c(inv,"ind.sup")
      }
      }
      
      vecinv<-NULL
      vecinv<-paste("'",vecinv,inv[1],"'",sep="")
      for (i in 2:(length(inv))){
        vecinv<-paste(vecinv,paste("'",inv[i],"'",sep=""),sep=",")
      }
      
      if(length(inv)>1){
        vecinv<-paste("c(",vecinv,")",sep="")
      }
      else if(length(inv)==1){
        vecinv<-paste("'",inv,"'",sep="")
      }
      else if(length(inv)==0){
        vecinv<-"NULL"
      }
      
      list(inv=(inv),vecinv=(vecinv))
    }
    
    getinv2=function(){
      inv<-c()
      if(!("suplquali"%in%input$var_sup)){
        inv<-c(inv,"quali.sup")
      }
      
      if(!("suplquanti"%in%input$var_sup)){
        inv<-c(inv,"quanti.sup")
      }
      
      if(!("act"%in%input$var_sup)){
        inv<-c(inv,"var")
      }
      
      vecinv<-NULL
      vecinv<-paste("'",vecinv,inv[1],"'",sep="")
      for (i in 2:(length(inv))){
        vecinv<-paste(vecinv,paste("'",inv[i],"'",sep=""),sep=",")
      }
      
      if(length(inv)>1){
        vecinv<-paste("c(",vecinv,")",sep="")
      }
      else if(length(inv)==1){
        vecinv<-paste("'",inv,"'",sep="")
      }
      else if(length(inv)==0){
        vecinv<-"NULL"
      }
      
      list(inv=(inv),vecinv=(vecinv))
    }
    
    #getlab=function(){
     # labels=input$indvarpoint
    #  if(!is.null(labels)){
     #   ok=TRUE
    #  }
     # else{
    #    ok=FALSE
     # }
    #  hide=c("ind","var","ind.sup","var.sup")
    #  if(!"Ind"%in%labels){
     #   hide=hide[-which(hide=="ind")]
      #}
      #if(!"Mod"%in% labels){
       # hide=hide[-which(hide=="var")]
      #}
      #if(!"Indsup"%in%labels){
      #  hide=hide[-which(hide=="ind.sup")]
      #}
    #  if(!"Modsup"%in%labels){
     #   hide=hide[-which(hide=="var.sup")]
    #  }
     # if(length(hide)==0&& ok==FALSE){
    #    hide=c("ind","var")
    #  }
     # if(length(hide)==0&& ok==TRUE){
    #    hide="none"
    #  }
     # return(hide)
    #}
    
    #GRAPHIQUE 3: Variables
    
    Plot4=reactive({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      inv=getinv2()$inv
      invtext=getinv2()$vecinv
      list(PLOT4=(plot.MCA(values()$res.MCA,choix="var",invisible=inv,title=input$title2,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)))),invisible=(invtext))    
    })
    
    output$map4 <- renderPlot({
      p <- Plot4()$PLOT4
    })
    
    
    #GRAPIQUE 1   
    
    Plot1=reactive({
      
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      )
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      
      validate(
        need(length(input$ind_var)>=1,"Please select the object you want to plot: Individuals, variables or both")
      )
      validate(
        need(input$habiller == TRUE || input$habiller == FALSE || length(input$habiller)<=2,"Please select maximum 2 variables as habillage")
      )
      
      
      #Selection des individus
      if(input$select=="Manuel"){
        selecindiv=c(input$indiv) 
        selecindivText=createVec(selecindiv)
      }
      else if(input$select=="cos2"){
        if(input$slider1!=1){
          selecindiv=paste("cos2",input$slider1)
        }
        else{
          selecindiv="cos2 0.999"
        }
        selecindivText=paste0("'",selecindiv,"'")
      }
      else if(input$select=="Contrib"){
        selecindiv=paste("contrib ",input$sliderContrib) 
        selecindivText=paste0("'",selecindiv,"'")
      }
      else if(input$select=="NONE"){
        selecindiv=NULL
        selecindivText="NULL"
      }
      
      #Selection des modalites
      
      if(input$selectMod=="cos2"){
        if(input$sliderCosMod!=1){
          selecMod=paste("cos2",input$sliderCosMod)
        }
        else{
          selecMod="cos2 0.999"
        }
        selecModText=paste0("'",selecMod,"'")
      }
      else if(input$selectMod=="Contrib"){
        selecMod=paste("contrib ",input$slider4)
        selecModText=paste0("'",selecMod,"'")
      }
      else if(input$selectMod=="NONE"){
        selecMod=NULL
        selecModText="NULL"
      }
      
      
      if(length(input$supvar)==0 || input$habi==FALSE){
        hab="none"
        habText<-"'none'"
        colquali="magenta"
      }
      
      if(length(quali)>1){
        if(length(input$habiller)==0){
          hab="none"
          habText<-"'none'"
          colquali="magenta"
        }
        
        if (length(input$habiller)==1 & input$habi==TRUE){
          hab=as.character(input$habiller)
          habText<-paste("'",input$habiller,"'",sep="")
          colquali="blue"
        }
        
        if (length(input$habiller)==2 & input$habi==TRUE){
          hab=dim(values()$DATA)[2]
          habText<-hab
          colquali="blue"
        }
      }
      ###
      else if (length(input$supvar)==1){
        if(input$habi==TRUE){
          hab=values()$choixqual
          habText<-hab
          colquali="blue"
        }
        else{
          hab="none"
          habText<-"'non'"
          colquali="magenta"
        }
      }
      
      ###
      
      validate(
        need(length(input$ind_var)!="","Please select which object you would like to print")
      )
      choixText="NULL"
      inv<-getinv()$inv
      invText<-getinv()$vecinv
      sel<-selecindiv
      selm<-selecMod
      colindsup<-"green"
     # hide=getlab()
      list(PLOT1=(plot.MCA(values()$res.MCA,choix="ind",invisible=inv,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),selectMod=selm,selec=sel,habillage=hab,col.quali=colquali,col.ind.sup=colindsup,title=input$title1)),choix=(choixText),inv=(invText),selm=(selecModText),sel=(selecindivText),hab=(habText),colquali=(colquali),colindsup=(colindsup))  
    })
    
    output$map <- renderPlot({
      p <- Plot1()$PLOT1
    })
    
    #GRAPHIQUE 2
    
    Plot2=function(){
      plot.MCA(values()$res.MCA,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),choix="quanti.sup",title=input$title3) 
    }
    
    output$map2 <- renderPlot({
      p=Plot2()
    })  
    
    output$map22=renderUI({
      validate(
        need(input$nb1 != input$nb2, "Please select two different dimensions")
      ) 
      validate(
        need(input$selecactive=="Toutes" || length(getactive())>2,"Please select quantitative variable and more active variables")
      )
      
      if(length(values()$choixquant)==0){
        return(p("No quantitative variables"))
      }
      else{
        plotOutput("map2", width = 500, height=500)
      }
    })
    ####
    
    
    output$choixchange=renderUI({
      if(length(values()$choixquant)==0){
        return(radioButtons("MCAgraph",h6("Which graph do you want to modify ?"),
                            choices=list("Individual and categories"="ind","Variables"="var"),inline=TRUE))
      }
      else{
        return(radioButtons("MCAgraph",h6("Which graph do you want to modify ?"),
                            choices=list("Individual and categories"="ind","Variables"="var","Quantitative variables"="quant"),inline=TRUE))
      }
    })
    
    
    output$habillage2=renderUI({
      #if(length(input$supvar)==0){
       # return(p("No supplementary categorical variable"))
      #}
      #if(length(input$supvar)>1){
        if(is.null(habillageind)){
        num=c(1:length(quali))
        return(selectInput("habiller","Select 1 or 2 variables", choices=list(num=quali),multiple=TRUE))
      }
      else{
        num=c(1:length(quali))
        return(selectInput("habiller","Select 1 or 2 variables", choices=list(num=quali),multiple=TRUE,selected=habillageind))
      }
      #}
    }) 
    
    #CALCUL DE LA CONTRIBUTION DES MODALITES
    
    output$slider3=renderUI({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      maxvar=dim(values()$res.MCA$var$coord)[1]

      if(selection3=="Contrib"){return(sliderInput("slider4",label="Contribution",
                                                   min=1,max=maxvar,value=as.numeric(selection4),step=1)) }
      else{
        return(sliderInput("slider4",label="Contribution",
                           min=1,max=maxvar,value=maxvar,step=1))
      }
    })
    
    ###
    
    
    ###
    #SUMMARY
    
    output$summary=renderPrint({
      summary(newdata)
    })
    
    
    #Histogramme des valeurs propres
    output$map3=renderPlot({
      return(barplot(values()$res.MCA$eig[,1],names.arg=rownames(values()$res.MCA$eig),las=2,density=TRUE))
    })
    
    #Histogramme du summary
    output$histo=renderPlot({
      barplot(prop.table(table(newdata[,input$bam]))*100)
    })
    
    #Summary de l'ACM
    
    
    output$summaryMCA=renderPrint({
      validate(
        need(input$nbele!=0, "Please select at least one element")
      )
      a<-values()$res.MCA  
      a$call$call<-code()
      
      summary.MCA(a,nbelements=input$nbele)
    })
    
    output$summary2=downloadHandler(filename = function() { 
      paste('summaryofMCA','.txt', sep='') 
    },
    content = function(file) {
      summary.MCA(values()$res.MCA,nbelements=input$nbele,file=file)
    },
    contentType='text/csv')
    
    #autre
    
    output$sorties=renderTable({
      return(as.data.frame(values()$res.MCA$eig))
    })
    
    output$sorties2=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(values()$res.MCA$var$coord))
    })
    
    output$sorties3=renderTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(values()$res.MCA$var$contrib))
    })
    
    output$sorties4=renderTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(values()$res.MCA$var$cos2))
    })
    
    output$sorties22=renderDataTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      tab<-as.data.frame(values()$res.MCA$ind$coord)
      tab<-round(tab, 3)
      tab<-cbind(Names=rownames(tab),tab)
      return(tab)
    })
    
    output$sorties33=renderDataTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      tab1<-as.data.frame(values()$res.MCA$ind$contrib)
      tab1<-round(tab1,3)
      tab1<-cbind(Names=rownames(tab1),tab1)
      return(tab1)
    })
    
    output$sorties44=renderDataTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      tab2<-as.data.frame(values()$res.MCA$ind$cos2)
      tab2<-round(tab2,3)
      tab2<-cbind(Names=rownames(tab2),tab2)
      return(tab2)
    })
    
    output$sorties23=renderTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$supvar)!=0, "No supplementary categorical variables")
      )
      return(as.data.frame(values()$res.MCA$quali.sup$coord))
    })
    
    output$sorties232=renderTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$supvar)!=0, "No supplementary categorical variables")
      )
      return(as.data.frame(values()$res.MCA$quali.sup$cos2))
    })
    
    output$sorties233=renderTable({
      
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$supvar)!=0, "No supplementary categorical variables")
      )
      return(as.data.frame(values()$res.MCA$quali.sup$v.test))
    })
    
    output$sorties43=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$supquanti)!=0 || input$supquanti==TRUE, "No supplementary quantitative variables")
      )
      return(as.data.frame(values()$res.MCA$quanti.sup$coord))
    })
    
    output$sortiesIsupC=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$indsup)!=0,"No supplementary individuals")
      )
      return(as.data.frame(values()$res.MCA$ind.sup$coord))
    })
    
    output$sortiesIsupCos=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$indsup)!=0,"No supplementary individuals")
      )
      return(as.data.frame(values()$res.MCA$ind.sup$cos2))
    })
    
    #DIM1
    
    output$sortieDimdesc=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[1]]$category))
    })
    
    output$sortieDimdesc2=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[1]]$quali))
    })
    output$sortieDimdesc3=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(input$supquanti)>0,"No quantitative variable"))
      validate(
        need(length(dimdesc(values()$res.MCA)[[1]]$quanti)!=0,"")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[1]]$quanti))
    })
    
    #DIM2
    output$sortieDimdesc00=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[2]]$category))
    })
    output$sortieDimdesc22=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[2]]$quali))
    })
    output$sortieDimdesc33=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(dimdesc(values()$res.MCA)[[2]]$quanti)!=0,"")
      )
      validate(
        need(length(input$supquanti)>0,"No quantitative variable"))
      return(as.data.frame(dimdesc(values()$res.MCA)[[2]]$quanti))
    })
    
    #DIM3
    output$sortieDimdesc000=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[3]]$category))
    })
    output$sortieDimdesc222=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      return(as.data.frame(dimdesc(values()$res.MCA)[[3]]$quali))
    })
    output$sortieDimdesc333=renderTable({
      validate(
        need(length(getactive())>2 || input$selecactive=="Toutes","Please select more variable")
      )
      validate(
        need(length(dimdesc(values()$res.MCA)[[3]]$quanti)!=0,"")
      )
      validate(
        need(length(input$supquanti)>0,"No quantitative variable"))
      return(as.data.frame(dimdesc(values()$res.MCA)[[3]]$quanti))
    })
    
    #Le JDDONNEES
    output$JDD=renderDataTable({
      cbind(Names=rownames(newdata),newdata)},
      
      options = list(    "orderClasses" = TRUE,
                         "responsive" = TRUE,
                         "pageLength" = 10))
    
    
    
    ####
    
    
    output$downloadData0 = downloadHandler(
      filename = function() { 
        paste('graph4','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot44()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData10 = downloadHandler(
      filename = function() { 
        paste('graph4','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot44()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData20 = downloadHandler(
      filename = function() { 
        paste('graph4','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot44()
        dev.off()
      },
      contentType=NA)
    
    ####
    
    
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
    
    
    output$download3 = renderUI({
      if(length(values()$choixquant)==0){
        return()
      }
      else{
        return(downloadButton("downloadData3","Download as png"))
      }
    })
    
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
    
    output$download4 = renderUI({
      if(length(values()$choixquant)==0){
        return()
      }
      else{
        return(downloadButton("downloadData4","Download as jpg"))
      }
    })
    
    output$downloadData4 = downloadHandler(
      filename = function() { 
        paste('graph1','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot2()
        dev.off()
      },
      contentType='image/jpg')
    
    
    output$download5 = renderUI({
      if(length(values()$choixquant)==0){
        return()
      }
      else{
        return(downloadButton("downloadData5","Download as pdf"))
      }
    })
    
    output$downloadData5 = downloadHandler(
      filename = function() { 
        paste('graph1','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot2()
        dev.off()
      },
      contentType=NA)    
    
    ####AXES
    
    output$NB1=renderUI({
      validate(
        need(length(getactive())>1 || input$selecactive=="Toutes","Please select at least one supplementary variables")
      )
      if(input$selecactive=="Toutes" || length(getactive())>5){
        return(selectInput("nb1", label = h6("x axis"), 
                           choices = list("1" = 1, "2" = 2, "3" = 3,"4"= 4,"5" =5), selected =axe1,width='80%'))
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
    
    Plot11=function(){
      if(input$select=="Manuel"){
      selecindiv=c(input$indiv) 
    }
    else if(input$select=="cos2"){
      if(input$slider1!=1){
        selecindiv=paste("cos2",input$slider1)
      }
      else{
        selecindiv="cos2 0.999"
      }
    }
    else if(input$select=="Contrib"){
      selecindiv=paste("contrib ",input$sliderContrib) 
    }
    else if(input$select=="NONE"){
      selecindiv=NULL
    }
    
    if(input$selectMod=="cos2"){
      if(input$sliderCosMod!=1){
        selecMod=paste("cos2",input$sliderCosMod)
      }
      else{
       selecMod="cos2 0.999" 
      }
    }
    else if(input$selectMod=="Contrib"){
      selecMod=paste("contrib ",input$slider4)
    }
    else if(input$selectMod=="NONE"){
      selecMod=NULL
    }
    
    
    if(length(input$supvar)==0 || input$habi==FALSE){
      hab="none"
      habText<-"'none'"
      colquali="magenta"
    }
    
    if(length(input$supvar)>1){
      if(length(input$habiller)==0){
        hab="none"
        habText<-"'none'"
        colquali="magenta"
      }
      
      if (length(input$habiller)==1 & input$habi==TRUE){
        hab=as.character(input$habiller)
        habText<-paste("'",input$habiller,"'",sep="")
        colquali="blue"
      }
      
      if (length(input$habiller)==2 & input$habi==TRUE){
        hab=dim(values()$DATA)[2]
        habText<-hab
        colquali="blue"
      }
    }
    ###
    else if (length(input$supvar)==1){
      if(input$habi==TRUE){
        hab=values()$choixqual
        habText<-hab
        colquali="blue"
      }
      else{
        hab="none"
        habText<-"'non'"
        colquali="magenta"
      }
    }
    choixText="NULL"
    inv<-getinv()$inv
    invText<-getinv()$vecinv
    sel<-selecindiv
    selm<-selecMod
    colindsup<-"green"
    plot.MCA(values()$res.MCA,choix="ind",title=input$title1,invisible=inv,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),selectMod=selm,selec=sel,habillage=hab,col.quali=colquali,col.ind.sup=colindsup)}
    Plot44=function(){inv=getinv2()$inv
                      plot.MCA(values()$res.MCA,choix="var",title=input$title2,invisible=inv,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)))}
    
  }
)