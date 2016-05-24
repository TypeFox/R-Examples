
shinyServer(function(input, output , session) {

  output$PIPEImage <- renderImage(list(src="PIPE.jpg") , deleteFile=FALSE)
  
  
  ndimA <- reactive({ input$dimA })
  ndimB <- reactive({ input$dimB })
  targ <- reactive({ input$theta  })
  
  prior.ss <- reactive({ m <- matrix(NA , nrow=ndimA() , ncol=ndimB() )
                         for(i in 1:ndimA() ) {
                         for(j in 1:ndimB() ) {
                            m[i,j] <- input[[paste("prior",i,j,sep="_")]]
                         } 
                         }
                         m
                         })
  
  pi <- reactive({ m <- matrix(NA , nrow=ndimA() , ncol=ndimB() )
  for(i in 1:ndimA() ) {
    for(j in 1:ndimB() ) {
      m[i,j] <- input[[paste("true",i,j,sep="_")]]
    } 
  }
  m
  })
  
  dat <- reactive({
    cohortsize <- input$cohortsize
    dose1s <- sapply(1:20 , function(i) input[[paste("dose1",i,sep="_")]])
    dose2s <- sapply(1:20 , function(i) input[[paste("dose2",i,sep="_")]])    
    dlts <- sapply(1:20 , function(i) input[[paste("DLT",i,sep="_")]])      
    dat <- data.frame(doseA=dose1s , doseB=dose2s , tox=dlts)
    dat <- dat[apply(!is.na(dat),1,all) , ]
    longdat <- data.frame(doseA=numeric(0) , doseB=numeric(0) , tox=numeric(0) )
    if(nrow(dat)>0) {
      longdat <- data.frame(doseA = rep(dat$doseA,each=cohortsize) , doseB=rep(dat$doseB,each=cohortsize))
      longdat$tox <- unlist(lapply(dat$tox , function(n) c(rep(1,n) , rep(0,cohortsize-n)) ))
      longdat$patient <- 1:nrow(longdat)
    }  
    longdat
  })
  
 
  des <- reactive({
    if(all(!is.na( prior.ss() ))) {
        des <- pipe.design(
        N = nrow( dat() ) + input$cohortsize ,
        c = input$cohortsize , 
        theta = input$theta , 
        prior.med = prior.ss() ,
        prior.ss = matrix( 1/(ndimA()*ndimB()*input$priorss) , ndimA() , ndimB() ), 
        strategy = input$strategy , 
        admis = input$admis , 
        constraint = input$constraint , 
        epsilon = input$epsilon ,
        mode = "sim",
        alternate = input$alternate ,
        data = dat(),
        seed=input$seed
        ) 
    }    
    else des <- NULL
    des
  })
  
simfn<-reactive({
    input$simulate
  
    isolate({
        return(
        
        pipe.design(
        N = input$N ,
        S = input$S ,
        c = input$cohortsize ,
        theta = input$theta ,
        pi = pi(),
        prior.med = prior.ss() ,
        prior.ss = matrix( 1/(ndimA()*ndimB()*input$priorss) , ndimA() , ndimB() ),
        strategy = input$strategy ,
        admis = input$admis ,
        constraint = input$constraint ,
        epsilon = input$epsilon ,
        mode = "sim",
        alternate = input$alternate ,
        data = dat(),
        seed=input$seed
      )
      ) 
    })
})
  
observe({
   if(input$update){
      isolate({
        n<-length(des()$rec.i)
        updateNumericInput( session , paste("dose1",n,sep="_") , value=des()$rec.i[n] )
        updateNumericInput( session , paste("dose2",n,sep="_") , value=des()$rec.j[n] )    
      })
  }
})
  
 observe({ 
   if(input$interpolate) {
    isolate(
      if(all(!is.na(prior.ss()[ c(1,ndimA()) , c(1,ndimB())]))) {
             m <- prior.ss()     
             for( j in 2:(ndimB()-1) ) {
               m[1,j] <- approx(x=c(1,ndimB()) , y=m[1,c(1,ndimB())] , xout=j)$y
               updateNumericInput( session , paste("prior",1,j,sep="_") , value=m[1,j] ) 
               m[ndimA(),j] <- approx(x=c(1,ndimB()) , y=m[ndimA(),c(1,ndimB())] , xout=j)$y
               updateNumericInput( session , paste("prior",ndimA(),j,sep="_") , value=m[ndimA(),j] ) 
             }
             for( j in 1:ndimB() ) {
               for(i in 2:(ndimA()-1) ) {
                 m[i,j] <- approx(x=c(1,ndimA()) , y=m[c(1,ndimA()),j] ,xout=i)$y
                 updateNumericInput( session , paste("prior",i,j,sep="_") , value=m[i,j] )         
             }}
             
             
            }
      ) 
    }
 })

observe({
  if(input$flat.prior) {
    isolate( {
      updateNumericInput( session , "prior_1_1" , value=targ()-0.03)
      updateNumericInput( session , paste("prior",1,ndimB(),sep="_") , value=targ()-0.02)    
      updateNumericInput( session , paste("prior",ndimA(),1,sep="_") , value=targ()-0.02)    
      updateNumericInput( session , paste("prior",ndimA(),ndimB(),sep="_") , value=targ())    
     } )    
  }
})


observe({ 
  if(input$reset) {
 #   isolate( {   
      for(i in 1:15) {
        updateNumericInput( session , paste("dose1",i,sep="_") , value=NA ) 
        updateNumericInput( session , paste("dose2",i,sep="_") , value=NA ) 
        updateNumericInput( session , paste("DLT",i,sep="_") , value=NA )       
      }
 #  } ) 
  }
})


observe({
  if(input$useprior) {
    isolate( {
      for(i in 1:ndimA() ) {
        for(j in 1:ndimB()){
          updateNumericInput( session , paste("true",i,j,sep="_"),value=input[[paste("prior",i,j,sep="_")]])
        }
      }
    }  )    
  }
})


 output$histplot <- renderPlot({
   if (!is.null( des() )) {  pipe.design:::plothists( des() )  }
  })    

 output$segplot <- renderPlot({
   if (!is.null( des() )) {
     pipe.design:::plotsegs(des())
  }   
})
 
 output$simplot <- renderPlot({
   if(input$simulate!=0){
     withProgress(message = 'Calculations in progress', value = 0, {plot(simfn())})
   }
 }) 
 
 
output$table<-renderTable({
  if(input$simulate!=0){
    exp.table<-matrix(print(simfn())$exp.table,nrow=1,dimnames=list("Experimentation percentages by true toxicity",names(print(simfn())$exp.table)))
    rec.table<-matrix(print(simfn())$rec.table,nrow=1,dimnames=list("Recommendation percentages by true toxicity",names(print(simfn())$rec.table)))
    tab<-rbind(exp.table,rec.table)
    xtable::xtable(tab)
  }
})
 
output$n <-reactive({
  n<-1
  for(i in 1:20){
    suppressWarnings(n<-n+all(!is.na(input[[paste("dose1",i,sep="_")]]),!is.na(input[[paste("dose2",i,sep="_")]])))
  }
  n
})
outputOptions(output, 'n', suspendWhenHidden=FALSE)
  
#   output$prior.ss <- renderTable({
#     print( rbind(  des()$admis.list[[length( des()$admis.list  )]]  , des()$dom.list[[length( des()$dom.list  )]]) )
#     if(nrow( dat() )>0) print( dat() )
#     print(prior.ss() )
# })    
  
 #make dynamic min and max Dose A
 output$doseA_1 <- renderUI({
    numericInput("dose1_1" , value=NA , min=1, max=input$dimA, label="Dose A")
 })
 output$doseA_2 <- renderUI({
   numericInput("dose1_2" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_3 <- renderUI({
   numericInput("dose1_3" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_4 <- renderUI({
   numericInput("dose1_4" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_5 <- renderUI({
   numericInput("dose1_5" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_6 <- renderUI({
   numericInput("dose1_6" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_7 <- renderUI({
   numericInput("dose1_7" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_8 <- renderUI({
   numericInput("dose1_8" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_9 <- renderUI({
   numericInput("dose1_9" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_10 <- renderUI({
   numericInput("dose1_10" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_11 <- renderUI({
   numericInput("dose1_11" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_12 <- renderUI({
   numericInput("dose1_12" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_13 <- renderUI({
   numericInput("dose1_13" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_14 <- renderUI({
   numericInput("dose1_14" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_15 <- renderUI({
   numericInput("dose1_15" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_16 <- renderUI({
   numericInput("dose1_16" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_17 <- renderUI({
   numericInput("dose1_17" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_18 <- renderUI({
   numericInput("dose1_18" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_19 <- renderUI({
   numericInput("dose1_19" , value=NA , min=1, max=input$dimA, label="")
 })
 output$doseA_20 <- renderUI({
   numericInput("dose1_20" , value=NA , min=1, max=input$dimA, label="")
 })
 
 #make dynamic min and max Dose B
 output$doseB_1 <- renderUI({
   numericInput("dose2_1" , value=NA , min=1, max=input$dimB, label="Dose B")
 })
 output$doseB_2 <- renderUI({
   numericInput("dose2_2" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_3 <- renderUI({
   numericInput("dose2_3" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_4 <- renderUI({
   numericInput("dose2_4" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_5 <- renderUI({
   numericInput("dose2_5" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_6 <- renderUI({
   numericInput("dose2_6" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_7 <- renderUI({
   numericInput("dose2_7" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_8 <- renderUI({
   numericInput("dose2_8" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_9 <- renderUI({
   numericInput("dose2_9" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_10 <- renderUI({
   numericInput("dose2_10" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_11 <- renderUI({
   numericInput("dose2_11" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_12 <- renderUI({
   numericInput("dose2_12" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_13 <- renderUI({
   numericInput("dose2_13" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_14 <- renderUI({
   numericInput("dose2_14" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_15 <- renderUI({
   numericInput("dose2_15" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_16 <- renderUI({
   numericInput("dose2_16" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_17 <- renderUI({
   numericInput("dose2_17" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_18 <- renderUI({
   numericInput("dose2_18" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_19 <- renderUI({
   numericInput("dose2_19" , value=NA , min=1, max=input$dimB, label="")
 })
 output$doseB_20 <- renderUI({
   numericInput("dose2_20" , value=NA , min=1, max=input$dimB, label="")
 })
 
 
 #make dynamic min and max DLTs
 output$dlt1 <- renderUI({
   numericInput("DLT_1" , value=NA , min=0, max=input$cohortsize,step=1, label="DLTs")
 })
 output$dlt2 <- renderUI({
   numericInput("DLT_2" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt3 <- renderUI({
   numericInput("DLT_3" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt4 <- renderUI({
   numericInput("DLT_4" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt5 <- renderUI({
   numericInput("DLT_5" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt6 <- renderUI({
   numericInput("DLT_6" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt7 <- renderUI({
   numericInput("DLT_7" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt8 <- renderUI({
   numericInput("DLT_8" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt9 <- renderUI({
   numericInput("DLT_9" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt10 <- renderUI({
   numericInput("DLT_10" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt11 <- renderUI({
   numericInput("DLT_11" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt12 <- renderUI({
   numericInput("DLT_12" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt13 <- renderUI({
   numericInput("DLT_13" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt14 <- renderUI({
   numericInput("DLT_14" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt15 <- renderUI({
   numericInput("DLT_15" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt16 <- renderUI({
   numericInput("DLT_16" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt17 <- renderUI({
   numericInput("DLT_17" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt18 <- renderUI({
   numericInput("DLT_18" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt19 <- renderUI({
   numericInput("DLT_19" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 output$dlt20 <- renderUI({
   numericInput("DLT_20" , value=NA , min=0, max=input$cohortsize,step=1, label="")
 })
 
 
})

