# server script for HCPC for dataframe
shinyServer(
  function(input, output) {
    
    values=reactive({
        data.selec=x[,VariableChoices]
        choixquali=NULL
        if(length(quali)!=0){
          data.selec=cbind(data.selec,x[,QualiChoice])
          choixquali=seq((length(data.selec)-length(quali)),length(data.selec))
        }
        
      list(res.PCA=(PCA(data.selec,quali.sup=choixquali,scale.unit=FALSE,graph=FALSE,ncp=Inf)),DATA=(data.selec))
      })
    
    res.HCPCdef=reactive({
      (HCPC(values()$res.PCA,nb.clust=-1,graph=FALSE)$call$t$nb.clust)
    })
    
    res.HCPC=reactive({
    (HCPC(values()$res.PCA,nb.clust=input$clust,consol=input$consoli,graph=FALSE,metric=input$metric))
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
      res$nomData=nomData
      res$data=x
      res$classx<-"data.frame"
      class(res) <- c("HCPCshiny")
      res$clust=input$clust
      res$consoli=input$consoli
      res$metric=input$metric
      res$drawtree=input$drawtree 
      res$nom3D=input$nom3D
      res$center=input$center
      res$num=input$num
      res$nb1=as.numeric(input$nb1)
      res$nb2=as.numeric(input$nb2)
      res$code1=Code()
      res$code2=Plot1Code()
      res$code3=Plot2Code()
      res$code4=Plot3Code()
      res$title1=input$title1
      res$title2=input$title2
      res$title3=input$title3
      return(res)
    }
    
    observe({
      if(input$HCPCcode==0){
      }
      else {
        isolate({
          cat(Code(),sep="\n")
          cat(Plot1Code(),sep="\n")
          cat(Plot2Code(),sep="\n")
          cat(Plot3Code(),sep="\n")
        })
      }
    })
    
    Code <- function(){
      Call1=as.name(paste("res.HCPC<-HCPC(",nomData,",nb.clust=",input$clust,",consol=",input$consoli,",graph=FALSE,metric='",input$metric,"')",sep="")) 
      return(Call1)
    }
    
    Plot1Code <- function(){
      Call2=paste("plot.HCPC(res.HCPC,choice='map',draw.tree=",input$drawtree,",title='",input$title2,"',axes=c(",as.numeric(input$nb1),",",as.numeric(input$nb2),"))",sep="") 
      return(Call2)
    }
    
    Plot2Code <- function(){
      Call3=paste("plot.HCPC(res.HCPC,choice='3D.map',ind.names=",input$nom3D,",centers.plot=",input$center,",title='",input$title,"',angle=",input$num,",axes=c(",as.numeric(input$nb1),",",as.numeric(input$nb2),"))",sep="") 
      return(Call3)
    }
    
    Plot3Code <- function(){
      Call4=paste("plot.HCPC(res.HCPC,choice='tree',title='",input$title3,"')",sep="")
      return(Call4)
    }
    
    getactive=function(){
      if(input$quantisup==TRUE){
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
    
    Plot1 <- function(){
      if(is.null(input$clust)){
        return()
      }
      else{
      return(plot.HCPC(res.HCPC(),choice="map",draw.tree=input$drawtree,title=input$title2,axes=c(as.numeric(input$nb1),as.numeric(input$nb2))))
      }
    }
    output$map <- renderPlot({
      p <- Plot1()
    })
    
    
    Plot2 <- function(){
      if(is.null(input$clust)){
        return()
      }
      else{
      return(plot.HCPC(res.HCPC(),choice="3D.map",ind.names=input$nom3D,title=input$title1,centers.plot=input$center,axes=c(as.numeric(input$nb1),as.numeric(input$nb2)),angle=input$num))
      }
    }
    
    output$map2 <- renderPlot({
      p <- Plot2()
    })
    
    output$map4=renderPlot({
      if(is.null(input$clust)){
        return()
      }
      else{
      return(plot.HCPC(res.HCPC(),choice="tree",title=input$title3))
      }
    })
    
    output$sorties=renderTable({
      if(input$out=="axe"){
        return(as.data.frame(res.HCPC()$desc.axes))
      }
      if(input$out=="para"){
        return(as.data.frame(res.HCPC()$ind.desc))
      }
    })
    

    output$clusters=renderUI({
      choix=res.HCPCdef()
      if(is.data.frame(x)==TRUE){
        if(nbindiv<=11){
          sliderInput("clust","Number of clusters",min=2,max=(nbindiv-1),value=choix,step=1)
        }
        else{
          sliderInput("clust","Number of clusters",min=2,max=10,value=choix,step=1)
        }
      }
      else{
      if(nbindiv<=11){
        sliderInput("clust","Number of clusters",min=2,max=(nbindiv-1),value=clustdf,step=1)
      }
      else{
        sliderInput("clust","Number of clusters",min=2,max=10,value=clustdf,step=1)
      }
      }
    })
    
    output$JDD=renderDataTable({
      cbind(Names=rownames(x),x)},
      options = list(    "orderClasses" = TRUE,
                         "responsive" = TRUE,
                         "pageLength" = 10))

    
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
    
    output$downloadData6 = downloadHandler(
      filename = function() { 
        paste('graph3','.png', sep='') 
      },
      content = function(file) {
        png(file)
        Plot4()
        dev.off()
      },
      contentType='image/png')
    
    output$downloadData7 = downloadHandler(
      filename = function() { 
        paste('graph3','.jpg', sep='') 
      },
      content = function(file) {
        jpeg(file)
        Plot4()
        dev.off()
      },
      contentType='image/jpg')
    
    output$downloadData8 = downloadHandler(
      filename = function() { 
        paste('graph3','.pdf', sep='') 
      },
      content = function(file) {
        pdf(file)
        Plot4()
        dev.off()
      },
      contentType=NA)
    


    
    ### Fonction permettant d'afficher la description des classes par les variables
    output$descript=renderTable({
      write.infile(X=res.HCPC()$desc.var$quanti,file=paste(getwd(),"essai.csv"),sep=";")
      baba=read.csv(paste(getwd(),"essai.csv"),sep=";",header=FALSE)
      colnames(baba)=NULL
      b=which(baba[,1]=="format non affichable")
      file.remove(paste(getwd(),"essai.csv")) 
      baba
    },
    include.rownames=FALSE)
    
    ### Fonction permettant d'afficher les parangons des classes
    output$parangons=renderTable({
      bibi=list()
      for (i in 1:input$clust){
        bibi[[i]]=rbind(colnames(res.HCPC()$desc.ind$para[[i]]),res.HCPC()$desc.ind$para[[i]])
      }
      write.infile(X=bibi,file=paste(getwd(),"essai3.csv"),sep=";",nb.dec=8)
      baba=read.csv(paste(getwd(),"essai3.csv"),sep=";",header=FALSE)
      colnames(baba)=NULL
      file.remove(paste(getwd(),"essai3.csv"))
      baba
    },
    include.rownames=FALSE)
    
    ### Fonction permettant d'afficher la description des classes par les axes 
    output$axes=renderTable({
      write.infile(X=res.HCPC()$desc.axes$quanti,file=paste(getwd(),"essai2.csv"),sep=";",nb.dec=8)
      baba=read.csv(paste(getwd(),"essai2.csv"),sep=";",header=FALSE)
      colnames(baba)=NULL
      file.remove(paste(getwd(),"essai2.csv"))
      baba
    },
    include.rownames=FALSE)
    
    
  }
)
      

