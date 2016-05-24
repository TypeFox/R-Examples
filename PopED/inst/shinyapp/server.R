library(PopED)
#source("models/MODEL.2.PK.one.comp.oral.R")

#sfg <- function(){}

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  updateDesign <- reactive({
    xt_txt <- input$xt
    xt <-  eval(parse(text=paste("c(",xt_txt,")")))
    return(list(xt=xt))
  })
  
  updateModel <- reactive({
    struct_model <- input$struct_model
    ruv_model <- input$ruv_model
    bsv_model <- input$bsv_model
    
    sfg <- build_sfg(model=input$struct_model,etas=input$bsv_model)
    #environment(eval(parse(text=input$struct_model)))
    #parent.env(environment())

    #browser()
    
    nbpop <- find.largest.index(func.str=sfg,lab="bpop") 
    #bpop_vals=c(CL=0.15, V=8, KA=1.0, Favail=1)
    #bpop_vals=c(CL=1, V=1, KA=1, Favail=1)
    #bpop_vals <- rep(1,nbpop)
    nb <- find.largest.index(func.str=sfg,lab="b")    
    
    
    
    if(input$struct_model=="ff.PK.1.comp.oral.sd.CL"){ 
      bpop_vals=c(CL=0.15, V=8, KA=1.0, Favail=1) 
      notfixed_bpop=c(1,1,1,0)
      d_vals=c(CL=0.07, V=0.02, KA=0.6) 
      sigma_vals=c(0.1,0.1)
      groupsize=32
      #xt=c( 0.5,1,2,6,24,36,72,120),
      minxt=0
      maxxt=120
      a=70
    }
    return(list(bpop=bpop_vals,d=d_vals,sigma=sigma_vals,
                notfixed_bpop=notfixed_bpop,
                sfg=sfg))
  })
  
  # Return the formula text for printing as a caption
  #output$caption <- renderText({
  #  "Model predictions"
  #})
  
  #get_parameters <- 
  #codetools::findGlobals(ff.PK.1.comp.oral.sd.CL,merge=F)
  
  # Generate a plot of the requested variable against mpg and only 
  # include outliers if requested
  output$modelPlot <- renderPlot({
    #     if(input$model=="one.comp") {
    #       source("Model.2.PK.one.comp.oral.R")
    #       poped.db <- poped.db.1
    #       facet_scales="fixed"
    #       #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    #     }
    #     if(input$model=="db.1") {
    #       source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.model.POPED.R")
    #       source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/one.comp.emax.design.POPED.R") 
    #       poped.db <- poped.db.2
    #       facet_scales="free"
    #       #print(plot_model_prediction(poped.db.2,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups,facet_scales="free"))
    #     }
    #     if(input$model=="db.2"){
    #       source("/Users/ahooker/Documents/_PROJECTS/PopED_in_R/poped_r/models/warfarin.model.design.all_in_one.POPED.R") # 4-group, add+prop, as pfim
    #       poped.db <- create.poped.database(warfarin.design.1.red.input())
    #       poped.db$design$a <- rbind(50,60,70,80)
    #       facet_scales="fixed"
    #       #print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #     }
    model <- updateModel()
    design <- updateDesign()
    
    #     poped.db <- create.poped.database(ff_file=input$struct_model,
    #                                       fError_file=input$ruv_model,
    #                                       fg_file="sfg",
    #                                       groupsize=32,
    #                                       sigma=model$sigma,
    #                                       #bpop=model$bpop,  
    #                                       bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
    #                                       d=model$d, 
    #                                       xt=design$xt,
    #                                       maxxt=120, 
    #                                       minxt=0, 
    #                                       a=70) 
    
    
    ff <- function(model_switch,xt,parameters,poped.db){
      ##-- Model: One comp first order absorption
      with(as.list(parameters),{
        y=xt
        y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
        return(list(y=y,poped.db=poped.db))
      })
    }
    
    sfg <- function(x,a,bpop,b,bocc){
      ## -- parameter definition function 
      parameters=c(CL=bpop[1]*exp(b[1]),
                   V=bpop[2]*exp(b[2]),
                   KA=bpop[3]*exp(b[3]),
                   Favail=bpop[4],
                   DOSE=a[1])
      return(parameters) 
    }
    
    feps <- function(model_switch,xt,parameters,epsi,poped.db){
      ## -- Residual Error function
      ## -- Proportional 
      returnArgs <- ff(model_switch,xt,parameters,poped.db) 
      y <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      y = y*(1+epsi[,1])
      
      return(list(y=y,poped.db=poped.db)) 
    }
    
    ## -- Define initial design  and design space
    poped.db <- create.poped.database(ff_file=input$struct_model,
                                      #ff_file="ff",
                                      fg_fun=model$sfg,
                                      #fError_file="feps",
                                      fError_file=input$ruv_model,
                                      #bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                      bpop=model$bpop,  
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(CL=0.07, V=0.02, KA=0.6), 
                                      sigma=c(0.01,0.1),
                                      groupsize=32,
                                      xt=design$xt,
                                      minxt=0,
                                      maxxt=120,
                                      a=70)
    
    #print(plot_model_prediction(poped.db))
    
    print(plot_model_prediction(poped.db,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #print(plot_model_prediction(poped.db.1,IPRED=input$IPRED,DV=input$DV,separate.groups=input$separate.groups))
    #print(plot_model_prediction(poped.db.2,IPRED=TRUE,DV=TRUE))
  })
})