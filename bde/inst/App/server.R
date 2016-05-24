library(bde)

message_plot<-function(mssg,size=10){
  ggplot(data.frame(x=0,y=0),aes(x=x,y=y)) + 
    annotate("text",label=mssg,x=0,y=0,size=size) + 
    theme(rect=element_blank(),line=element_blank(),axis.text=element_blank(),axis.title=element_blank())
}

shinyServer(function(input, output) {
  #**********************************
  # SESSION VARIABLES  -------------- 
  #**********************************
  values<-reactiveValues()
  values$dataPath<-""
  values$data<-NULL
  values$bandwidth<-NULL
  values$m<-NULL
  values$min_value<-NULL
  values$max_value<-NULL
  values$update_count<-0
  
  #*********************************
  # GET DATA FUNCTION -------------- 
  #*********************************
  
  getData<-reactive({
    data<-values$data
    if(!is.null(input$file) && values$dataPath!=input$file$datapath){
      data<-unlist(read.csv(file=input$file$datapath,header=F))
      values$data<-data
      values$dataPath=input$file$datapath
      values$bandwidth<-signif(length(data)^{-2/5},4)
      values$m<-max(1,round(length(data)^{2/5}))
      values$min_value<-signif(min(0,data),4)
      values$max_value<-signif(max(1,data),4)
    }
    return(data)
  })
  
  
  #***************************
  # GET DENSITY -------------- 
  #***************************
  
  getDensity<-reactive({
    sample<-getData()
    app_dens<-NULL
    ##The second argument is to avoid running this before the options are displayed
    if (!is.null(sample)){
      values$update_count<-input$update
      b<-input$bandwidth
      opt<-switch(input$method,
                    "betakernel" = {
                      l<-list(modified=input$beta_modified,normalization=input$beta_normalization,mbc=input$beta_mbc)
                      if(!is.null(input$beta_c)) l$c=input$beta_c
                      l
                    },
                    "vitale" = {
                      l<-list(biasreduced=input$vitale_br, M=input$vitale_M)
                      b<-1/input$m
                      l
                    },
                    "boundarykernel"={
                      l<-list(mu=as.numeric(input$boundary_mu),nonegative=input$boundary_nonegative,method=input$boundary_kern)
                      l
                    })
        m<-min(input$low_bound, min(sample))
        M<-max(input$up_bound,max(sample))
        x<-seq(m,M,(M-m)/input$num_points)
        dens<-bde(dataPoints=unlist(sample),dataPointsCache=x,estimator=input$method,b=b,options=opt,lower.limit=input$low_bound,upper.limit=input$up_bound)
        app_dens<-switch(input$kaki_method,
                         "none"={
                           a<-dens
                           a
                         },
                         "b1"={
                           a<-bde(dataPoints=unlist(sample),dataPointsCache=x,estimator="kakizawa",options=list(method="b1",c=input$kaki_c),lower.limit=input$low_bound,upper.limit=input$up_bound)
                           a
                         },
                         "b2"={
                           a<-bde(dataPoints=unlist(sample),dataPointsCache=x,estimator="kakizawa",options=list(method="b2"),lower.limit=input$low_bound,upper.limit=input$up_bound)
                           a
                         },
                         "b3"={
                           a<-bde(dataPoints=unlist(sample),dataPointsCache=x,estimator="kakizawa",options=list(method="b3"),lower.limit=input$low_bound,upper.limit=input$up_bound)
                           a
                         })
      
    }
    
    return(app_dens)
  })
  
  
  #****************************
  # PLOT DENSITY -------------- 
  #****************************
  
  getDensityPlot<-reactive({
    dens<-getDensity()
    if (!is.null(dens)){
      if (input$points_alpha==0){
        plot<-gplot(dens)
      }else{
        a<-as.numeric(input$points_alpha)/100
        plot<-gplot(dens,includePoints=T,alpha=a)
      }
    }else{
      plot<-message_plot("There is no data loaded.\nClick on the \'Choose File\'  button and select\na csv file containing the sample points",size=7)
    }
    plot
  })
  
  output$plot_density<-renderPlot({
    plot<-getDensityPlot()
    print(plot)
  })
  

#*********************************
# SAVE PLOT DENSITY -------------- 
#*********************************

output$density_plot_save<-downloadHandler(filename=paste("plot.",input$density_plot_format,sep=""),
                                     content=function(file){
                                       switch(input$density_plot_format,
                                              "pdf" = {
                                                pdf(file=file,
                                                    width=input$density_plot_width/input$density_plot_resolution,
                                                    height=input$density_plot_height/input$density_plot_resolution)
                                                },
                                              "eps" = {
                                                postscript(file=file,
                                                           width=input$density_plot_width/input$density_plot_resolution,
                                                           height=input$density_plot_height/input$density_plot_resolution)
                                                },
                                              "png" = {
                                                png(file=file,
                                                    width=input$density_plot_width,
                                                    res=input$density_plot_resolution,
                                                    height=input$density_plot_height)
                                                },
                                              "jpg" = {
                                                jpg(file=file,
                                                    width=input$density_plot_width,
                                                    res=input$density_plot_resolution,
                                                    height=input$density_plot_height)
                                                }
                                       )
                                       plot<-getDensityPlot()
                                       print(plot)
                                       dev.off()           
                                     },
                                     contentType=paste('.',input$density_plot_format,sep=''))

#*********************************
# PLOT DISTRIBUTION -------------- 
#*********************************

getDistributionPlot<-reactive({
  dens<-getDensity()
  if (!is.null(dens)){
    x<-getdataPointsCache(dens)
    cdf<-distribution(dens,x=x,discreteApproximation=T,scaled=F)
    df<-data.frame(X=x,CDF=cdf)
    plot<-ggplot(df,aes(x=X,y=CDF)) + geom_line(size=1.1)
  }else{
    plot<-message_plot("There is no data loaded.\nClick on the \'Load data\'  button and select\na csv file containing the sample points",size=7)
  }
  plot
})

output$plot_distribution<-renderPlot({
  plot<-getDistributionPlot()
  print(plot)
})


#**************************************
# SAVE PLOT DISTRIBUTION -------------- 
#**************************************

output$distribution_plot_save<-downloadHandler(filename=paste("plot.",input$distribution_plot_format,sep=""),
                                          content=function(file){
                                            switch(input$distribution_plot_format,
                                                   "pdf" = {
                                                     pdf(file=file,
                                                         width=input$distribution_plot_width/input$distribution_plot_resolution,
                                                         height=input$distribution_plot_height/input$distribution_plot_resolution)
                                                   },
                                                   "eps" = {
                                                     postscript(file=file,
                                                                width=input$distribution_plot_width/input$distribution_plot_resolution,
                                                                height=input$distribution_plot_height/input$distribution_plot_resolution)
                                                   },
                                                   "png" = {
                                                     png(file=file,
                                                         width=input$distribution_plot_width,
                                                         res=input$distribution_plot_resolution,
                                                         height=input$distribution_plot_height)
                                                   },
                                                   "jpg" = {
                                                     jpg(file=file,
                                                         width=input$distribution_plot_width,
                                                         res=input$distribution_plot_resolution,
                                                         height=input$distribution_plot_height)
                                                   }
                                            )
                                            plot<-getDistributionPlot()
                                            print(plot)
                                            dev.off()           
                                          },
                                          contentType=paste('.',input$distribution_plot_format,sep=''))



#***************************
# DATA TABLE  -------------- 
#***************************

output$samplePoints<-renderDataTable({
  data<-getData()
  df<-NULL
  if (!is.null(data)){
    df<-data.frame(Sample=unlist(data))
  }
  df
})


  #****************************
  # MAIN OPTIONS  ------------- 
  #****************************
  
  output$main_opts<-renderUI({
    data<-getData()
    if (is.null(data)){
      list()
    }else{
      m<-values$min_value
      M<-values$max_value
      list(h4("Main options"),
           selectInput(inputId="method",label="Method",
                       choices=list("Beta kernel"="betakernel",
                                    "Vitale"="vitale",
                                    "Boundary kernel"="boundarykernel")),
           conditionalPanel(condition="input.method==\'vitale\'",
                            h5("Approximation degree"),
                            numericInput(inputId="m",label="",value=values$m,min=1,max=10*values$m,step=1)),
           conditionalPanel(condition="input.method!=\'vitale\'",
                            h5("Bandwidth"),
                            numericInput(inputId="bandwidth",label="",value=values$bandwidth,min=values$bandwidth/10,max=10*values$bandwidth,step=values$bandwidth/10)),
           h5("Lower bound"),
           numericInput(inputId="low_bound",label="",value=m,step = (max(data)-min(data))/200,max=min(data)),
           h5("Upper bound"),
           numericInput(inputId="up_bound",label="",value=M,step = (max(data)-min(data))/200,min=max(data)),
           h5("Number of evaluation points"),
           numericInput(inputId="num_points",label="",min=10,max=1000,value=200),
           h5("Transparency of the sample points"),
           numericInput(inputId="points_alpha",label="",min=0,max=100,value=0))
    }
  })

  
  #********************************
  # ADITIONAL OPTIONS ------------- 
  #********************************
  
  output$additional_opts<-renderUI({
    if (!is.null(input$method)){
      opt<-switch(input$method,
             "betakernel"={
               elems<-list(
                 selectInput(inputId="beta_normalization",label="Normalization",
                                                  choices=list("None" = "none",
                                                               "Density-wise" = "densitywise",
                                                               "Kernel-wise" = "kernelwise")),
                 checkboxInput(inputId="beta_modified",label="Modified beta kernel"),
                 selectInput(inputId="beta_mbc",label="Multiplicative bias correction",
                                                  choices=list("None"="none",
                                                               "JNL" = "jnl",
                                                               "TS" = "ts")),
                 conditionalPanel(condition="input.beta_mbc==\'ts\'",
                                  h5("C parameter"),
                                  numericInput(inputId="beta_c",label="",
                                               value=0.5,min=0,max=1,step=0.01)))
               elems
           },
           "vitale"={
             theom<-round(1/(2*input$bandwidth)) 
             elems<-list(
               checkboxInput(inputId="vitale_br",label="Bias reduction",value=T),
               conditionalPanel(condition="input.vitale_br",
                                numericInput(inputId="vitale_M","M",value=max(1,theom),min=max(1,theom/10),max=theom*10,step=1)))
             elems
           },
           "boundarykernel"={
             elems<-list(
               selectInput(inputId="boundary_kern",label="Boundary kernel",
                                                choices=list("Muller '94" = "muller94",
                                                             "Muller '91" = "muller91",
                                                             "Normalized" = "normalized",
                                                             "None"="none")),
               conditionalPanel(condition="input.boundary_kern==\'muller94\' | input.boundary_kern==\'muller91\'",
                                                      checkboxInput(inputId="boundary_nonegative",label="Nonegativity correction")),
               selectInput(inputId="boundary_mu",label="Type of kernel",
                           choices=list("Rectangular"="0",
                                        "Epanechnikov" = "1",
                                        "Bicubic" = "2",
                                        "Tricubic" = "3")))
             elems
           },
           "kakizawa"={
             elems<-list()
             elems
           })
      list(h4("Additional options"),
           opt,
           br(),hr(),br(),
           h4("Kakizawa's approximation"),
           selectInput(inputId="kaki_method",label="Method",
                                          choices=list("None"  = "none",
                                                       "B1" = "b1",
                                                       "B2" = "b2",
                                                       "B3" = "b3")),
           numericInput(inputId="kaki_c",label="C parameter",value=0.5,min=0,max=1,step=0.001))
     }
  })
})
