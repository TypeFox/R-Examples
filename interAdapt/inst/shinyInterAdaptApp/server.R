# Notes:
# Avoid unicode
# The central function of this script is the regen() function, which determines, based on the current parameters, whether or not table_constructor() from the "[...]Design.R" file should be called to update the results in the main panel


  # ______                         _     _      
  # | ___ \                       | |   | |     
  # | |_/ / __ ___  __ _ _ __ ___ | |__ | | ___ 
  # |  __/ '__/ _ \/ _` | '_ ` _ \| '_ \| |/ _ \
  # | |  | | |  __/ (_| | | | | | | |_) | |  __/
  # \_|  |_|  \___|\__,_|_| |_| |_|_.__/|_|\___|
                                              
                                              


#Are we on the shiny server?
onRStudioServer <- 'onRStudio.txt' %in% dir()
# The local and RStudio versions of the shiny app must differ in a few ways.
# The onRStudioServer variable lets us copy paste from the local version to the RStudio
# version without having to adjust files by hand.


###########
#Functions

subH0 <- function(x){ #make a function that does the same thing as "strong()" but for <sub></sub>
#finds all H0C and H01 terms and subs them
  x <- strong(x)
  x <- gsub("<strong>", "", x)
  x <- gsub("</strong>", "", x)
  x <- gsub("H0C", "H<sub>0C</sub>", x)
  x <- gsub("H01", "H<sub>01</sub>", x)
  x <- gsub("H02", "H<sub>02</sub>", x)
  return(x)
}

#To be used in xtable function
subH01sanitize<-function(x){
  x <- gsub("H0C", "H<sub>0C</sub>", x)
  x <- gsub("H01", "H<sub>01</sub>", x)
  return(x)
}
options(xtable.sanitize.text.function=subH01sanitize)


#To print2log stores errors that we can read in a log later, useful for checking sessions on RStudio (glimmer or spark)
#first override last session (start new log)
#In final release, set "print2R" to FALSE
cat(file='session_log.txt',paste(Sys.time(),'\n \n')) 
print2log<-function(x,logFileName='session_log.txt',print2R=FALSE){ #takes a string as input
  if(print2R) print(x)
  cat(file=logFileName,paste(x,'\n'),append=TRUE)
}
###############



print2log("source'ing code...")

#Load initial inputs
#Source code must be sourced inside the server input function for the
#functions to look in the local user's environment 
#for the relevant variables. If defined here, the functions
#look in the global env.
#See scoping information for shiny apps
st<-read.csv(file= "sliderTable.csv",header=TRUE,as.is=TRUE)
bt<-read.csv(file= "boxTable.csv",header=TRUE,as.is=TRUE)
print2log("found code locally...")

#If online, make sure the max time limit is not infinite
#Also be sure to load kniter packages (which are imported in the local package version).
# on the local version, it's OK for the time limit input max to be inf
# but on the RStudio server this would crash things.
if(onRStudioServer){
  time_limit_ind<-which(bt[,1]=='time_limit')
  bt[time_limit_ind,'max']<- min(90,bt[time_limit_ind,'max'])
  library(knitr)
  library(knitcitations)
}

print2log("...supplementary files found and loaded...")






  #  _____                          _____           _      
  # /  ___|                        /  __ \         | |     
  # \ `--.  ___ _ ____   _____ _ __| /  \/ ___   __| | ___ 
  #  `--. \/ _ \ '__\ \ / / _ \ '__| |    / _ \ / _` |/ _ \
  # /\__/ /  __/ |   \ V /  __/ |  | \__/\ (_) | (_| |  __/
  # \____/ \___|_|    \_/ \___|_|   \____/\___/ \__,_|\___|
                                                         
                                                         




shinyServer(function(input, output) {


  ##########
  # (_)     (_) | (_)     | (_)        
  #  _ _ __  _| |_ _  __ _| |_ _______ 
  # | | '_ \| | __| |/ _` | | |_  / _ \
  # | | | | | | |_| | (_| | | |/ /  __/
  # |_|_| |_|_|\__|_|\__,_|_|_/___\___|
                                     
  # Initialize some static & reactive variables
  ##########

  if(onRStudioServer){
    #make sure file size hasn't blown up.
    #if not, write the current user time.
    if(file.info("user_log.txt")$size < 1000000)
      cat(paste(date(),'\n'),file='user_log.txt',append=TRUE)
  }

  #Functions must be defined in local env. as they call user specific objects
  source("Adaptive_Group_Sequential_Design.R", local=TRUE) #need local=TRUE for the functions to be loaded to the user-specific env.

  #####
  #process initial inputs from CSVs
  allVarNames<-c(st[,'inputId'],bt[,'inputId'])
  allVarLabels<-c(st[,'label'],bt[,'label'])
  lastAllVars<-rep(0,length(allVarNames)) #for use later on, for shiny to tell when inputs have changed or not.
  names(lastAllVars)<-allVarNames

  for(i in 1:dim(st)[1]){
    assign(st[i,'inputId'], st[i,'value'])
    lastAllVars[st[i,'inputId']] <- st[i,'value']
  }
  for(i in 1:dim(bt)[1]){
    assign(bt[i,'inputId'], bt[i,'value']) 
    lastAllVars[bt[i,'inputId']] <- bt[i,'value']
  }

    
  # If the default inputs to the files have not changed since last time interAdapt was run, then we don't have to redo the initial calculations.
  # table1 stores the results needed to display performance of each trial.
  # The following code answers the question: do we need to regenerate table1? If we need to update it, do so & save the results
  stillNeedTable1<-TRUE
  try({
  load('last_default_inputs.RData') #This may fail the first time the app is run.
    if(all(bt==lastBt)&all(st==lastSt)){ #lastBt and lastSt are from the last time we generated table1, and are contained in last_default_inputsR.Data. If we're a match, then:
        load('last_default_table1_&_assigned_vars.RData')
        stillNeedTable1<-FALSE
        print2log("loaded table1...")
    }
  })
  if(stillNeedTable1){ 
    
    #these variables are all assigned via <<- in table_constructor(). 
    #we define them here so that their value is stored in the local env. when
    #we call table_constructor, next.

    futility_boundaries_standard_design_H0C<-
    futility_boundaries_standard_design_H01<-
    H0C_efficacy_boundary_proportionality_constant_adaptive_design<-
    H01_efficacy_boundary_proportionality_constant_adaptive_design<-
    H0C_efficacy_boundary_proportionality_constant_standard_design<-
    H01_efficacy_boundary_proportionality_constant_standard_design<-
    subpop_1_efficacy_boundaries_adaptive_design<-
    subpopulation_2_stopping_boundaries_adaptive_design<-
    subpop_1_futility_boundaries_adaptive_design<-
    risk_difference_list<-NULL

    table1<- table_constructor()
    lastBt<-bt
    lastSt<-st
    save(list=c('table1',
        'futility_boundaries_standard_design_H0C',
        'futility_boundaries_standard_design_H01',
        'H0C_efficacy_boundary_proportionality_constant_adaptive_design',
        'H01_efficacy_boundary_proportionality_constant_adaptive_design',
        'H0C_efficacy_boundary_proportionality_constant_standard_design',
        'H01_efficacy_boundary_proportionality_constant_standard_design',
        'subpop_1_efficacy_boundaries_adaptive_design',
        'subpopulation_2_stopping_boundaries_adaptive_design',
        'subpop_1_futility_boundaries_adaptive_design',
        'risk_difference_list'),
        file='last_default_table1_&_assigned_vars.RData')

    save(list=c('lastBt','lastSt'),file='last_default_inputs.RData')
    print2log("built table1...")
    stillNeedTable1<-FALSE
  }
  ######

  
  ########
  #For tracking how much we need to update things
  lastApplyValue <- 0 # need to put -1 here if we don't load table 1 beforehand
  totalCalls<-0 #total number of times regen() has been called so far
  #for use in uploading files:
  uploadCsvTicker<-0
  uploadDatasetTicker<-0
  inCsvValues<-NULL
  inDatasetValues<-NULL

  #current value of all the data, need to store this if we want to change the sliders to all be animated in interactive mode, but not change their values.
  #all use a comparison against static lastAllVars
  #Also do error checks for invalid inputs
  allVars<-reactive({
    x<-c()
    for(i in 1:length(allVarNames)) x[allVarNames[i]]<- input[[allVarNames[i] ]]

    #Check to make sure all box inputs are within the required ranges. 
    #Don't need to do this for sliders, since you can't set them past the min/max
    minMaxErrs<-rep('',dim(bt)[1]) #vector to store error/warnings messages. The ith entry is '' if allVars()[i] is valid
    for(i in 1:dim(bt)[1]){
      nameInd<- i+dim(st)[1]
      minMaxErrs_ind<-FALSE
      if( x[allVarNames[nameInd]]>bt[i,'max'])  {
        x[allVarNames[nameInd]]<-bt[i,'max']
        minMaxErrs_ind<-TRUE
      }
      if( x[allVarNames[nameInd]]<bt[i,'min'])  {
        x[allVarNames[nameInd]]<-bt[i,'min']
        minMaxErrs_ind<-TRUE
      }
      if(minMaxErrs_ind) minMaxErrs[i]<- paste0('Warning: the variable "',bt[i,'label'], '" is outside the allowed range, and has been set to ',x[allVarNames[nameInd]],'. ')
    }
    output$warn3<-renderText({paste(minMaxErrs,collapse='')})

    #Other error checks on inputs.
    warn2<-""
    if(x['total_number_stages']<x['last_stage_subpop_2_enrolled_adaptive_design']){
        x['last_stage_subpop_2_enrolled_adaptive_design']<-x['total_number_stages']
        warn2<-paste("Warning: the last stage sub population 2 is enrolled must be less than the total number of stages. Here the last stage in which sub population 2 is enrolled is set to",x['total_number_stages'],"the total number of stages")
    }
    output$warn2<-renderText({warn2})

    #Done! Send back the error-checked list of inputs
    x
  })


  params<-reactive({ input$Parameters1 + input$Parameters2 }) #Tracks the number of button presses for advanced paremeter panel (input$Parameters2), one for basic parameter panel (intput$Parameters1).

  #When advanced parameters are visible to user (which_params!=1), we force batch mode.
  effectivelyBatch<- reactive({input$Batch == "1" | input$Which_params != "1"})

  # STOP alert
  output$stop <- renderText({
	x <- input$stopButton
	if(x > 0) {
		stopApp(x)
		x <- "((( application stopped )))"
	}
	else
		x <- ""
	x
  })

  #                          (_)                
  # __      ____ _ _ __ _ __  _ _ __   __ _ ___ 
  # \ \ /\ / / _` | '__| '_ \| | '_ \ / _` / __|
  #  \ V  V / (_| | |  | | | | | | | | (_| \__ \
  #   \_/\_/ \__,_|_|  |_| |_|_|_| |_|\__, |___/
  #                                    __/ |    
  #                                   |___/     
  output$warn1<-renderText({
    x<-""
    if(input$Which_params!="1" & input$Batch=="2")x<-"Note: interactive mode not enabled for advanced parameters, defaulting to batch mode."
    x    
  })
#See also the allVars() for other warnings
  


  #  _                     _  ______      _        
  # | |                   | | |  _  \    | |       
  # | |     ___   __ _  __| | | | | |__ _| |_ __ _ 
  # | |    / _ \ / _` |/ _` | | | | / _` | __/ _` |
  # | |___| (_) | (_| | (_| | | |/ / (_| | || (_| |
  # \_____/\___/ \__,_|\__,_| |___/ \__,_|\__\__,_|

  #below are 2 reactive chunks to feed to the dynamicBoxes & dynamicSliders



  csvUpload<-reactive({
    upFile <- input$uploadCsvInput
    x<-NULL
    if (!is.null(upFile)){
      x<-c(read.csv(file=upFile$datapath, row.names=1, header=FALSE))[[1]]
      names(x)<-allVarNames
      uploadCsvTicker<<-0
      print2log('resetting upload csv all inputs')
    }
    inCsvValues<<-x
  })





  datasetVarNames<-c('p1_user_defined','p10_user_defined','p11_user_defined','p20_user_defined','p21_user_defined')

  datasetUpload<-reactive({
    upFile <- input$uploadDataTable
    x<-NULL
    if (!is.null(upFile)){
      tmp <- read.csv(file=upFile$datapath, header=TRUE)
      S <- tmp$S        # subpopulation, 1 or 2
      A <- tmp$A        # study arm, 0 or 1
      Y <- tmp$Y        # outcome, 0 or 1
      x['p1_user_defined'] <- mean(S==1)
      x['p10_user_defined'] <- mean(Y*(S==1)*(A==0))/mean((S==1)*(A==0))
      x['p11_user_defined'] <- mean(Y*(S==1)*(A==1))/mean((S==1)*(A==1))
      x['p20_user_defined'] <- mean(Y*(S==2)*(A==0))/mean((S==2)*(A==0))
      x['p21_user_defined'] <- mean(Y*(S==2)*(A==1))/mean((S==2)*(A==1))
      uploadDatasetTicker<<-0
      print2log('new Data calculated from')
      print2log(x)
      print2log('resetting upload dataset')      
    }
    inDatasetValues<<-x
  })






#      _ _     _                          _                        
#     | (_)   | |                 ___    | |                       
#  ___| |_  __| | ___ _ __ ___   ( _ )   | |__   _____  _____  ___ 
# / __| | |/ _` |/ _ \ '__/ __|  / _ \/\ | '_ \ / _ \ \/ / _ \/ __|
# \__ \ | | (_| |  __/ |  \__ \ | (_>  < | |_) | (_) >  <  __/\__ \
# |___/_|_|\__,_|\___|_|  |___/  \___/\/ |_.__/ \___/_/\_\___||___/

#Code to create dynamic inputs, which are rendered by ui.R

  #NOTE - June 26 2013: When sliders update, regen() thinks it needs to be called again because sliders have updated values and you're now in interactive mode.
        #solution -- added lastAllVars (nonreactive) variable to cancel this out.
  #explanation of uploadCsvTicker: Only if it's the first time since uploading do we change the sliders to the uploaded values from the full csv list. If it's not the first time, we use the current value of the variable, taken from allVars.
  #for uploadDatasetTicker, it's the same basic idea, but we only check it for variables in the vector 'datasetVarNames'.
  sliders <- reactive({ #NOTE: if you upload the same file again it won't update because nothing is techincally new
    print2log('Slider inputs')
    print2log(uploadCsvTicker)
    print2log(uploadDatasetTicker)
    csvUpload() #reactive input to uploaded file; sets uploadCsvTicker to zero when it's activated
    datasetUpload() #reactive input to uploaded file; sets uploadDatasetTicker to zero when it's activated
    labelSliderList<-list()
    animate<-FALSE
    if(input$Batch=='2' & input$Which_params=='1' ) animate<-TRUE # reactive input
    for(i in 1:dim(st)[1]){
      #each of these cases overrides the previous one
      #case1: upfile=null  & uploadCsvTicker=zero  (sliders haven't changed yet)
      value_i<-st[i,'value']
      #case2: something has been uploaded in csv input file, but hasn't been incorporated yet (uploadCsvTicker==0)
      if(!is.null(inCsvValues)) value_i<-inCsvValues[st[i,'inputId']]
      #case3: uploaded csv has been used already, so ignore uploaded input csv data.
      if(uploadCsvTicker>0) isolate(value_i<-allVars()[st[i,'inputId']])
      #case4: user just uploaded a dataset that we parsed to get slider values
        #(only do this for sliders that are in the datasetVarNames vector)
      if(uploadDatasetTicker==0 & (!is.null(inDatasetValues)) & st[i,'inputId'] %in% datasetVarNames)
        value_i<-inDatasetValues[st[i,'inputId']]
      #end of cases

      labelListi<-subH0(st[i,'label']) #Labels are stored as non slider text objects, so that we can apply subscript styling
      sliderListi<-sliderInput(inputId=st[i,'inputId'], label='', min=st[i,'min'], max=st[i,'max'], value=value_i, step=st[i,'step'], animate=animate)
      ind<-length(labelSliderList) #starts at 0, grows on each repeat of this loop
      labelSliderList[[ind+1]]<-labelListi #alternate labels inbetween sliders
      labelSliderList[[ind+2]]<-sliderListi
    }
    uploadCsvTicker<<-uploadCsvTicker+1
    uploadDatasetTicker<<-uploadDatasetTicker+1
    print2log('............ sliders updating')
    labelSliderList
  })
  output$fullSliders<-renderUI({sliders()})


  boxes <- reactive({ #NOTE: if you upload the same file again it won't update because nothing's techincally new
    print2log('Box inputs')
    csvUpload()
    labelBoxList<-list()
    for(i in 1:dim(bt)[1]){
      value_i<-bt[i,'value']
      if( (!is.null(inCsvValues)) ){ value_i<-inCsvValues[bt[i,'inputId']]}
      boxLabeli<-subH0(bt[i,'label'])
      boxListi<-numericInput(inputId=bt[i,'inputId'], label='', min=bt[i,'min'], max=bt[i,'max'], value=value_i, step=bt[i,'step'])
      ind<-length(labelBoxList)           
      #add extra text between boxes:
      #'Lower bound for...'
      if(grepl('Lowest value to plot',bt[i,'label'])){
        labelBoxList[[ind+1]]<-strong("For use in Plots of Power vs. Average Treatment Effect for Subpopulation 2:")
        labelBoxList[[ind+2]]<-br()
        ind<-length(labelBoxList)
      }
      #text to follow "Delta"
      if(grepl('Delta',bt[i,'label'])){
        labelBoxList[[ind+1]]<-strong("Applies to all Parameters:")
        labelBoxList[[ind+2]]<-br()
        ind<-length(labelBoxList)
      }
      labelBoxList[[ind+1]]<-boxLabeli
      labelBoxList[[ind+2]]<-boxListi
      labelBoxList[[ind+3]]<-br()      
    }
    labelBoxList
  })
  output$fullBoxes<-renderUI({boxes()})








# ______ _____ _____  _____ _   _ 
# | ___ \  ___|  __ \|  ___| \ | |
# | |_/ / |__ | |  \/| |__ |  \| |
# |    /|  __|| | __ |  __|| . ` |
# | |\ \| |___| |_\ \| |___| |\  |
# \_| \_\____/ \____/\____/\_| \_/

  # In interactive mode, we re-export the parameters and rebuild table1 (see Design.R file) every time.  In batch, we rebuild table1 only on the first call for a given push of the Apply button.
  # applyValue is usally fed in as params(), the sum of the two apply buttons on the basic and advanced parameter input panels shown to the user.



  regen <- reactive({
    applyValue<-params()
    totalCalls<<-totalCalls+1
    print2log(paste('total regen calls =',totalCalls))
    #ESCAPE SCENARIOS
    #escape if it's called when dynamic sliders/buttons are still loading (just when app starts up)
    if(is.null(input$Batch) || is.null(applyValue))
      {print2log('regen null batch or apply -> out') ; return()}
    for(name in allVarNames){
      isolate(if(is.null(input[[name]])) {print2log('null regen input -> out');return()})
    }
    #in batch mode: if buttons have not been pressed since last time, 
  	if (  effectivelyBatch() && (applyValue <= lastApplyValue) )
  	    return()
    #In batch or interactive, check for no change -- especially useful for slider asthetic changes.
        #no need to isolate this next line -- if we're not in interative mode we would have already exited by now.
    if(all(allVars()==lastAllVars)) {
      print2log('no change')
      lastApplyValue <<- applyValue #fixes bug where if you hit apply with no changes, the next slider change will interactively call table_constructor()
      return()
    }


    #If we haven't escaped, continue...


    #if Batch==1, we don't want things updating automatically, so we use isolate()
    #JUNE 26 2013 - This must be done in conjunction with isolating the is.null tests at the start of the function.
    #need to use allVars() so that we head off some potential errors
    assignAllVars<-function(){
      for(i in 1:length(allVarNames))
        assign(allVarNames[i], allVars()[allVarNames[i]], inherits=TRUE)
    }
    print2log("assigning Variables ...")
  	if (effectivelyBatch()){ isolate(assignAllVars() )
    } else { assignAllVars() }

    print2log("making table1 ...")
    table1 <<- table_constructor()
    print2log("Done \n")

  	if (effectivelyBatch() ) lastApplyValue <<- applyValue
    lastAllVars<<-allVars()
})





  # ______ _       _       
  # | ___ \ |     | |      
  # | |_/ / | ___ | |_ ___ 
  # |  __/| |/ _ \| __/ __|
  # | |   | | (_) | |_\__ \
  # \_|   |_|\___/ \__|___/
                         


  #############
  #Performance plots

  output$power_curve_plot <- renderPlot({
	regen() 
  print2log('power plot')
	power_curve_plot()
  })

  output$expected_sample_size_plot <- renderPlot({
	regen()
  print2log('sample plot')
	expected_sample_size_plot()
  })

  output$expected_duration_plot <- renderPlot({
	regen()
  print2log('duration plot')
	expected_duration_plot()
  })


  #############
  #Boundary Plots
  output$standard_H0C_boundary_plot <-renderPlot({
    regen()
    print2log('H0C Boundary Plot')
    boundary_standard_H0C_plot()
  })

  output$standard_H01_boundary_plot <-renderPlot({
    regen()
    print2log('H01 Boundary Plot')
    boundary_standard_H01_plot()
  })

  output$adapt_boundary_plot <-renderPlot({
    regen()
    print2log('H0C Boundary Plot')
    boundary_adapt_plot()
  })





  #  _____     _     _           
  # |_   _|   | |   | |          
  #   | | __ _| |__ | | ___  ___ 
  #   | |/ _` | '_ \| |/ _ \/ __|
  #   | | (_| | |_) | |  __/\__ \
  #   \_/\__,_|_.__/|_|\___||___/
  #####


  # In our code, the input to xtable has both a dynamic data value, and a dynamic digits arguemnt.
  # We want to pass both of these reactively to xtable, but shiny's renderTable is set up to just take the data, and feed that data to xtable's "x" argument.
  # To adjust for this, we create a new xtable.list function, create an environment to enclose a copy of renderTable, and add our new xtable function to that env. Then, our new function can intercept our new renderTable function's calls to xtable.list. (our_renderTable will look for xtable.list in the environment we specify).
  # The result is that we pass a list to our_renderTable, our_renderTable passes that list to xtable.list, and that last call is intercepted by our_table.list within our new env.
  our_xtable.list <- function(x, ...){
    # our_xtable.list is a custom xtable function to be applied to be applied to objects of type 'list'.
    # input to our_xtable.list is a list of length 3.
    xtable::xtable(x[[1]], digits=x$digits, caption=x$caption, ...) 
  }
  env_for_our_renderTable <- new.env(parent = environment(shiny::renderTable)) #create a child of the shiny package env.
  our_renderTable <- shiny::renderTable 
  environment(our_renderTable) <- env_for_our_renderTable #change the environment of renderTable to our new custom env. Since this env. is a child of renderTable's natural env., any shiny objects we haven't defined will still be found by renderTable.
  env_for_our_renderTable$xtable.list <- our_xtable.list #Bind our new xtable func within our new env.This intercepts calls to xtable (applied to list objects), and redirect them towards our_table.list function.
  #This code is adapted from:
  #http://stackoverflow.com/questions/8204008/redirect-intercept-function-calls-within-a-package-function


  #copies of all the tables have to be made to put in different panels of the shiny app (hence "[...]table.2")

  output$adaptive_design_sample_sizes_and_boundaries_table.2 <-
  output$adaptive_design_sample_sizes_and_boundaries_table <- our_renderTable(
    {
    	regen()
      print2log('adaptive design table')
    	adaptive_design_sample_sizes_and_boundaries_table()
    })



  output$standard_H0C_design_sample_sizes_and_boundaries_table.2 <-
  output$standard_H0C_design_sample_sizes_and_boundaries_table <- our_renderTable({
    	regen()
      print2log('H0C standard trial table')
    	standard_H0C_design_sample_sizes_and_boundaries_table()
    })



  output$standard_H01_design_sample_sizes_and_boundaries_table.2 <-
  output$standard_H01_design_sample_sizes_and_boundaries_table <-our_renderTable(
    {
      regen()
      print2log('H01 standard trial table')
      standard_H01_design_sample_sizes_and_boundaries_table()
    })


  output$performance_table <- our_renderTable(expr=
    {
    	regen()
      print2log('performance table')
    	transpose_performance_table(performance_table())
    },include.colnames=FALSE)






  #                            _       _        
  #                           | |     | |       
  #  ___  __ ___   _____    __| | __ _| |_ __ _ 
  # / __|/ _` \ \ / / _ \  / _` |/ _` | __/ _` |
  # \__ \ (_| |\ V /  __/ | (_| | (_| | || (_| |
  # |___/\__,_| \_/ \___|  \__,_|\__,_|\__\__,_|

  #saving current parameters to a csv
  output$downloadInputs <- downloadHandler(
    filename =  paste0('inputs_',gsub('/','-',format(Sys.time(), "%D")),'.csv'),
    contentType =  'text/csv',
    content = function(filename) {
      inputCsv<-rep(NA,length=length(allVarNames))
      for(i in 1:length(allVarNames)) inputCsv[i]<- input[[ allVarNames[i] ]]
      write.table(inputCsv, filename, row.names=allVarLabels, col.names=FALSE, sep=',')
    }
  )

  #generate a knitr report
  output$knitr <- downloadHandler(
    filename =  'report.html',
    contentType =  'text/html',
    content = function(filename) {
      if (file.exists('knitr_report.html')) file.remove('knitr_report.html')
      if (file.exists('knitr_report.md')) file.remove('knitr_report.md')
      htmlKnitted<-knit2html('knitr_report.Rmd',quiet=TRUE) #"plain" version, without knitrBootstrap
      x<-readLines(con=htmlKnitted) #"plain" version, without knitrBootstrap
      #library(knitrBootstrap)
      #knit_bootstrap('knitr_report.Rmd') #fancy knitrBootstrap version
      #x<-readLines(con='knitr_report.html')#fancy knitrBootstrap version
      writeLines(x,con=filename)
      # file.rename('knitr_report.html', filename)

    }
  )

  #functions for downloading table output.
  roundTable<-function(tab,digits){
    newTab<-array(0,dim=dim(tab))
    for(i in 1:dim(tab)[1]){
      for(j in 1:dim(tab)[2]){
        newTab[i,j]<-round(tab[i,j],digits=digits[i,j])
      }
    }
    rownames(newTab)<-rownames(tab)
    return(newTab)
  }
  designTable2csv<-function(t1,filename){
      K<-dim(t1[[1]])[2]
      designsCsv<-rbind( 
        'labeltext'=rep(NA,K),
        'Stage'=1:K,
        roundTable(tab=t1[[1]],digits=t1[[2]][,-1])
      )
      rownames(designsCsv)[1]<-t1[[3]]
      write.table(designsCsv, filename, row.names=TRUE, col.names=FALSE, sep=',')

  }

  # Generate csv tables
  # It could be good to double check these output files below if we change the decision rules at some point, to make sure the tables still come out right.
  
  output$downloadDesignAD.1<-
  output$downloadDesignAD.2 <- downloadHandler(
    filename =  paste0('DesignAD_',gsub('/','-',format(Sys.time(), "%D")),'.csv'),
    contentType =  'text/csv',
    content = function(filename) {
      t1<-adaptive_design_sample_sizes_and_boundaries_table()
      designTable2csv(t1,filename)
    }
  )
  output$downloadDesignSC.1<-
  output$downloadDesignSC.2<- downloadHandler(
    filename =  paste0('DesignSC_',gsub('/','-',format(Sys.time(), "%D")),'.csv'),
    contentType =  'text/csv',
    content = function(filename) {
      t1<-standard_H0C_design_sample_sizes_and_boundaries_table()
      designTable2csv(t1,filename)
    }
  )
  output$downloadDesignSS.1<-
  output$downloadDesignSS.2 <- downloadHandler(
    filename =  paste0('DesignSS_',gsub('/','-',format(Sys.time(), "%D")),'.csv'),
    contentType =  'text/csv',
    content = function(filename) {
      t1<-standard_H01_design_sample_sizes_and_boundaries_table()
      designTable2csv(t1,filename)
    }
  )

  output$downloadPerformance.1<-
  output$downloadPerformance.2<-
  output$downloadPerformance.3<-
  output$downloadPerformance.4<- downloadHandler(
    filename =  paste0('Performance_',gsub('/','-',format(Sys.time(), "%D")),'.csv'),
    contentType =  'text/csv',
    content = function(filename) {
      t1<-transpose_performance_table(performance_table())
      perfTab<-t1[[1]]
      for(row in 1:dim(perfTab)[1]) perfTab[row,]<-round(t1[[1]][row,],digits=t1[[2]][row,2])
      perfCsv<-rbind('labeltext'=NA,perfTab)
      rownames(perfCsv)[1]<-t1[[3]]
      write.table(perfCsv, filename, row.names=TRUE, col.names=FALSE, sep=',')
    }
  )




})
