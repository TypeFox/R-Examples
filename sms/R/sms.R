
##' Generate small area population microdata from census and survey datasets. 
##' Fit the survey data to census area descriptions and export the population of small areas (microdata).
##' 
##' Generate small area population microdata from census and panel datasets. 
##' Fit the survey data to census area descriptions and export the popultion of small areas.
##' 
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @name sms-package
##' @docType package
##' @references Dimitris Kavroudakis D (2015). \strong{sms: An R Package for the Construction of Microdata for 
##' Geographical Analysis.} \emph{Journal of Statistical Software}, \strong{68}(2), pp. 1-23. \url{http://10.18637/jss.v068.i02}
##' @title Spatial Microsimulation Library
NULL


# ---------------------   Data   ------------------------------------------------

#' A survey dataset of 200 individuals
#' 
#' A sample survey dataset containing binary (0 or 1) information about 200 individuals. 
#' Those individuals will be used to populate the simulated areas. 
#'    The variables in the dataset are as follows:
#' \itemize{
#'   \item pid: The unique indentifier of the individual  
#'   \item female: Binary value of the sex of the individual. 1-Female, 0-Male 
#'   \item agemature: Binary value indicating if the individual belongs to the mature age group. 
#'        0-No, 1-Yes
#'   \item car_owner: Binary value indicating if the individual owns a car. 0-No, 1-Yes
#'   \item house_owner: Binary value indicating if the individual owns a house. 0-No, 1-Yes
#'   \item working: Binary value indicating if the individual is working. 0-No, 1-Yes
#' }
#' @docType data
#' @keywords datasets
#' @name survey
#' @usage data(survey)
#' @format A data frame with 200 rows and 7 variables
NULL


#' A census dataset of 10 areas
#' 
#' A sample census dataset containing descriptive information about 10 geographical areas. 
#' The variables in the dataset are as follows:
#' \itemize{
#'   \item areaid: The unique indentifier of the area  
#'   \item population: The number of indivisuals in the area. 
#'   \item he: Number of individuals in the area, with at least Higher Education degree
#'   \item females: Number of female individuals in the area
#' }
#' @docType data
#' @keywords datasets
#' @name census
#' @usage data(census)
#' @format A data frame with 10 rows and 4 variables
NULL

#====================== Class ========================

##' A microsimulation object
##' 
##' It holds all microsimulation details and objects such as data, results etc.
##' @param census: A census data.frame where each row contains census information about 
##' a geographical area
##' @param panel: A data.frame containing the individual based records from a panel survey. 
##' Those data will be fitted to small area contrains and will populate each vrtual area.
##' @param lexicon: A data.frame containing the association of columns between census data and 
##' panel data. Each row contain a conection between census and panel data.frame.
##' @param resuls: A list of results from the fitting process.
##' @param iterations: The number of itertions until th end of the fitting process.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @name microsimulation-class
##' @exportClass microsimulation
setClass("microsimulation",
         representation(
           census="data.frame",panel="data.frame",lexicon="data.frame",results="ANY",
           iterations="numeric"),
         prototype(
           census=data.frame(),panel=data.frame(),lexicon=data.frame(),results=list(),iterations=0)
)

##' getInfo Generic
##' 
##' getInfo Generic
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @param object A microsimulation object to get its information.
##' @aliases getInfo,microsimulation-generic
setGeneric("getInfo", function(object) {
  standardGeneric("getInfo")
  #cat("Generic getInfo \n")
})

##' getInfo Method
##' 
##' Get information from a microsimulation object
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @param object A microsimulation object to get its information.
##' @exportMethod getInfo
##' @aliases getInfo,microsimulation-method
setMethod("getInfo", signature(object = "microsimulation"), function(object) {
  cat(paste0("A microsimulation object with ",length(object@panel)," panel records and ",
             length(object@census)," census areas\n"))
})



##' getTAEs Generic
##' 
##' Get the TAE from a microsimulation object.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @param object A microsimulation object to get its information.
##' @aliases getTAEs,microsimulation-generic
setGeneric("getTAEs", function(object) {
  standardGeneric("getTAEs")
})

##' getTAEs Method
##' 
##' getTAEs Method
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @param object A microsimulation object to get its information.
##' @return taes A list of numbers indicating the Total Absolute Error of the fitting process 
##' for each of the census areas.
##' @exportMethod getTAEs
##' @aliases getTAEs,microsimulation-method
setMethod("getTAEs", signature(object = "microsimulation"), function(object) {
  taes=sapply(object@results, '[[', "tae")
  return(taes)
})

#====================== Simulation Methods ========================


##' mysetSeed
##' 
##' mysetSeed 
##' @title mysetSeed
##' @param inseed A number to set as a random seed.
##' @export
##' @examples library(sms)
##' sms::mysetSeed(1900)
mysetSeed=function(inseed){
  set.seed(inseed)
}

##' Select n random rows from a dataframe
##'
##' Select n random rows from a dataframe
##' @title random_panel_selection
##' @param indf The initial dataframe from wich a selection will be made.
##' @param n The number of random rows
##' @return a selection of rows as a dataframe
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @export
##' @examples library(sms)
##' data(survey) #load the data
##' data(census)
##'     
##' some.individuals=random_panel_selection(survey,4)
##' print(some.individuals)     # Print the selection of individuals
random_panel_selection=function(indf,n){
  #print(paste0("Selecting ", n, " random individuals from the survey data."))
  
  return(indf[sample(nrow(indf), n, replace=TRUE), ]) 
}


##' Calculate the error of a selection.
##'
##' Calculates the Total Absolute Error (TAE) of a selection for a census area.
##' @title Calculate error of a selection
##' @param selection A population selection, to evaluate its error
##' @param area_census An area from census (a row)
##' @param lexicon A data.frame with details about data connections
##' @return TAE Total Absolute Error of this selection against the census description of this area.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com} 
##' @export
##' @examples library(sms)
##' data(survey) #load the data
##' data(census)
##' in.lexicon=createLexicon() # Create a data lexicon for holding the associated column names.
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' #Select the first area from the census table
##' this_area=as.data.frame(census[1,]) 
##' 
##' #make a random selection of individuals for this area.
##' selection=random_panel_selection( survey, this_area$population ) 
##' 
##' #evaluate the Total Absolute Error (TAE) for this selection
##' error=calculate_error( selection, this_area, in.lexicon ) 
##' print( error )    # print the error of the selection
calculate_error=function(selection,area_census,lexicon){
  check_lexicon(lexicon)  #checks
  nvar=length(lexicon)
  myNewColumnNames=c()
  for (i in 1:nvar){  # for every variable in the lexicon
    column_name=colnames( lexicon )[i]# con_01, con_02
    myNewColumnNames=c(myNewColumnNames,as.character(lexicon[column_name]["census_row",]))
  }
  area.errors=c(areaid="")  #Add one column
  area.errors[myNewColumnNames]="" #add the other columns: cars, mature
  this_area_id=area_census$areaid
  area.error=data.frame(areaid=this_area_id)  # prepare the Dataframe with the errors for this area.
  tae=0
  for (i in 1:nvar){# for every variable in the lexicon
    column_name=colnames( lexicon )[i]# con_01, con_02
    #print(area_census)
    at.census= area_census[[ as.character(lexicon[column_name]["census_row",]) ]]# cars 5, mature 6
    panel.column=as.character(lexicon[column_name]["survey_row",])
    at.selection= sum(selection[[panel.column]])   
    this.var.error=abs(at.census - at.selection)
    tae=tae+this.var.error
    area.errors[panel.column]=this.var.error
  }  
  return(tae)
}


##' Make a single selection of individual records for a census area.
##'
##' Select a number of individual records from panel dataset, 
##' to represent a census description of an area.
##' @title selection_for_area
##' @param inpanel The panel dataset
##' @param area_census A census area
##' @param inlexicon A data lexicon showing the variable associations.
##' @return list A list of results (#areaid, #selection, #error)
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @export
##' @examples library(sms)
##' data(survey) #load the data
##' data(census)
##' in.lexicon=createLexicon() # Create a data lexicon for holding the associated column names.
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' # Select the first area from the census table
##' this_area=as.data.frame(census[1,]) 
##' 
##' #make a representation for this area.
##' sel=selection_for_area(survey, this_area, in.lexicon) 
##' 
##' print(sel) #print the representation
selection_for_area=function(inpanel, area_census, inlexicon){
  this_area_id=area_census$areaid
  this_area_pop=area_census$population
  selection=random_panel_selection(inpanel,this_area_pop)
  error=sms::calculate_error(selection,area_census,inlexicon)
  return(list(areaid=this_area_id,selection=selection,error=abs(error)))
}

##' Plot the selection process of an area from a microsimulation object.
##'
##' Plot errors during selection process for an area.
##' @title Plot selection results
##' @param insms The input results
##' @param number the number of the area to plot
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' data(survey) #load the data
##' data(census)
##' in.lexicon=createLexicon() # Create a data lexicon for holding the associated column names.
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' ansms = new("microsimulation", census=census, panel=survey, lexicon=in.lexicon, iterations=5)
##' sa = run_parallel_SA(ansms, inseed=1900)
##' plotTries( sa, 1 )
##' @export
plotTries <- function(insms, number){
  if (!is(insms, "microsimulation")){
    stop("You gave me a non-microsimulation object to plot.\n\t\t
         Please give me a proper microsimulation object.")
  }
  if ((length(insms@results)<1)){
    stop("There are no results to plot. Run the simulation first and then plot the results.")
  }
  #print("Inresult")
  #print(nrow(inresult$selection))
  inresult=insms@results[[number]]
  states=inresult$error_states
  tries=inresult$tries
  lim=range(states,tries)
  plot(1:length(states),ylim=lim,states, type="b", pch=20, xlab="Iterations", 
       ylab="Total Absolute Error")
  points(1:length(tries), tries, col="red")
  title(main = list(paste("Area",inresult$areaid), cex=0.9, font=4))
  mtext(cex=0.8,paste("Improvements:",length(unique(states))-1," ","TAE:",inresult$tae, 
                      "\nPopulation:",nrow(inresult$selection) )
  )#the error improvents towards 0 error.
}

##' Find the best selection of individual records for a census area.
##'
##' Calculate the best area representation, after a series of selection tries.
##' @title find_best_selection
##' @param area A census area
##' @param insms A microsimulation object which holds the data and details of 
##' the simulation such as iterations, lexicon.
##' @param inseed test
##' @return list A list with results (#areaid, #selection, #tae, #tries, #error_states).
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @export
##' @examples library(sms)
##' data(survey) #load the data
##' data(census)
##' in.lexicon=createLexicon() # Create a data lexicon for holding the associated column names.
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' this_area=as.data.frame(census[1,]) #Select the first area from the census table
##' insms= new("microsimulation",census=census,panel=survey, lexicon=in.lexicon, iterations=10)
##' best=find_best_selection(this_area, insms)
##' print(best)
find_best_selection<- function(area,insms, inseed=-1){
  area_census=as.data.frame(area)
  it=insms@iterations
  bhps=insms@panel
  my.lexicon=insms@lexicon
  old_error=NULL
  errors_list=c()
  current_error_list=c()
  result_selection=NULL
  
  if (inseed>0){
    set.seed(inseed)
  }
  
  
  for (i in 1:it){
    current_selection=sms::selection_for_area(bhps,area_census, my.lexicon)
    if (i==1){# At first time
      result_selection=current_selection #Just in case we got the best selection.
      old_error=abs(current_selection$error) #Keep the initial error
    }
    # Append the current selection error.
    current_error_list=c(current_error_list,current_selection$error)
    
    if (old_error > abs(current_selection$error) ) {#If this new error is better
      old_error = abs(current_selection$error) #keep the error
      result_selection=current_selection #keep the selection
    }
    errors_list=c(errors_list,old_error) #append the error to the errors_list
  }
  return(list(areaid=result_selection$areaid,
              selection=result_selection$selection,
              tae=result_selection$error,
              tries = current_error_list, 
              error_states = errors_list)
  )
}



##' Run a simulation in serial mode with Hill Climbing
##'
##' Run a simulation in serial mode with Hill Climbing
##' @title run_parallel_HC
##' @param insms A microsimulation object which holds the data and details 
##' of the simulation such as iterations, lexicon.
##' @param inseed A number to be used for random seed.
##' @return msm_results An object with the results of the simulation, for each area.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' data(survey) #load the data
##' data(census)
##' in.lexicon=createLexicon() # Create a data lexicon for holding the associated column names.
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' insms= new("microsimulation",census=census,panel=survey, lexicon=in.lexicon, iterations=10)
##' re=run_parallel_HC(insms, inseed=1900)
##' print(re)
##' 
##' @export
#run_parallel <- function(census,iterations,bhps, the.lexicon){
run_parallel_HC <- function(insms, inseed=-1){
  #options(cores=2)
  #library(parallel)
  #library(doParallel)
  #library(foreach)
  #library(doSNOW)
  #library(iterators)
  cores=2
  cl <- parallel::makePSOCKcluster(cores ) 
  doParallel::registerDoParallel(cl)
  #cores=parallel::detectCores()
  #cl<-parallel::makeCluster(cores,type="PSOCK")
  #doSNOW::registerDoSNOW(cl)
  #doParallel::registerDoParallel(cl)
  myFunctions=c("find_best_selection",  "calculate_error", 
                "selection_for_area", "random_panel_selection")
  
  msm_results=
    foreach::foreach(i=iterators::iter(insms@census, by='row'), .export=myFunctions, combine=c) %dopar% {
      
      sms::find_best_selection(i,insms, inseed)
    }
  parallel::stopCluster(cl)
  insms@results= msm_results
  i <- NULL 
  rm(i)
  return(insms)
}



##' Run a simulation in serial mode
##'
##' Run a simulation in serial mode.
##' @title Run_serial
##' @param insms A microsimulation object which holds the data and details 
##' of the simulation such as iterations, lexicon.
##' @return msm_results An object with the results of the simulation, for each area.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' data(survey)
##' data(census)
##' in.lexicon=createLexicon()
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' insms= new("microsimulation",census=census, panel=survey, lexicon=in.lexicon, iterations=5)
##' results= run_serial( insms)
##' print(results)
##' @export
run_serial <- function(insms){
  #library(parallel)
  #library(foreach)
  #library(doSNOW)
  #library(doParallel)
  myFunctions=c("find_best_selection_SA", "calculate_error", "selection_for_area", 
                "random_panel_selection")
  msm_results=
    foreach::foreach(i=iterators::iter(insms@census, by='row'), .export=myFunctions, combine=c) %do% {
      sms::find_best_selection(i,insms)
    }
  i <- NULL 
  rm(i)
  return(msm_results)
}

##' Run a simulation in parallel mode with Simulated Annealing
##'
##' @title find_best_selection_SA
##' @param area_census A census dataset consisting of various areas rows.
##' @param insms A microsimulation object which holds the data and details 
##' of the simulation such as iterations, lexicon.
##' @param inseed A number to be used for random seed.
##' @return msm_results An object with the results of the simulation, of this area.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' data(survey)
##' data(census)
##' in.lexicon=createLexicon()
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' this_area=as.data.frame(census[1,]) #Select the first area from the census table
##' insms= new("microsimulation",census=census, panel=survey, lexicon=in.lexicon, iterations=5)
##' myselection= find_best_selection_SA( this_area, insms, inseed=1900)
##' print(myselection)
##' @export
find_best_selection_SA <- function (area_census, insms, inseed=-1  ) {
  iterations=insms@iterations
  bhps=insms@panel
  in.lexicon=insms@lexicon
  #area_census=as.data.frame(area_census)
  old_error=NULL
  errors_list=c()
  current_error_list=c()
  result_selection=NULL
  tolerance=iterations/2
  
  if (inseed>0){
    set.seed(inseed)
  }
  
  
  for (i in 1:iterations){
    current_selection=selection_for_area(bhps,area_census, in.lexicon)
    if (i==1){# At first time
      old_error=abs(current_selection$error) #Keep the initial error
      result_selection=current_selection #Just in case we got the best selection.
    }
    
    # Append the current selection error.
    current_error_list=c(current_error_list,current_selection$error)
    
    if (old_error > abs(current_selection$error) ) {#If this new error is better
      old_error = abs(current_selection$error) #keep the error
      result_selection = current_selection #keep the selection
    }
    else if( tolerance > 0 ){
      print("before dice")
      dice=sample(c(F,T),1)
      
      #roll a boolean dice
      #dice=sample(c(F,T),1, prob=c((iterations/2) , (tolerance/(iterations/2)) )) 
      print("after dice")
      if(dice){ #if the dice gives TRUE
        old_error = abs(current_selection$error) #keep the error
        result_selection = current_selection #keep the selection
      }
    }
    tolerance=tolerance-1
    errors_list=c(errors_list,old_error) #append the error to the errors_list
  }
  return(list(areaid=result_selection$areaid,selection=result_selection$selection,
              tae=result_selection$error,
              tries = current_error_list,error_states = errors_list))
}


##' Run a simulation in parallel mode with Simulated Annealing
##'
##' @title run_parallel_SA
##' @param insms A microsimulation object which holds the data and details 
##' of the simulation such as iterations, lexicon.
##' @param inseed A random number to be used for random seed.
##' @return msm_results An object with the results of the simulation, for each area.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' data(survey)
##' data(census)
##' in.lexicon=createLexicon()
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' 
##' insms= new("microsimulation",census=census, panel=survey, lexicon=in.lexicon, iterations=5)
##' results= run_parallel_SA(insms, inseed=1900)
##' print(results)
##' @export
run_parallel_SA <- function(insms, inseed=-1){
  census=insms@census 
  iterations=insms@iterations 
  bhps=insms@panel 
  in.lexicon=insms@lexicon 
  #options(cores=2)
  #library(parallel)
  #library(foreach)
  #library(doParallel)
  
  #library(Rmpi)
  #library(doSNOW)
  #cores=parallel::detectCores()
  cores=2
  cl <- parallel::makePSOCKcluster(cores ) 
  doParallel::registerDoParallel(cl)
  #cl<-parallel::makeCluster(parallel::detectCores(),type="SOCK")
  #registerDoSNOW(cl)
  
  doParallel::registerDoParallel(cl)
  myFunctions=c("find_best_selection_SA", "calculate_error", "selection_for_area", 
                "random_panel_selection")
  msm_results=
    foreach::foreach(i=iterators::iter(census, by='row'), .export=myFunctions, combine=c) %dopar% {
      
        find_best_selection_SA(i,insms, inseed)
      
      
    }
  parallel::stopCluster(cl)
  i <- NULL 
  rm(i) 
  insms@results= msm_results
  return(insms)
}









##' Create a data lexicon for holding the associated column names
##'
##' @title createLexicon
##' @return dataLexicon A data.frame holding the associated column names.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @export
##' @examples library(sms)
##' data(survey)
##' data(census)
##' in.lexicon=createLexicon()
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' print(in.lexicon)
createLexicon=function(){
  test=data.frame(con_1=c(NA,NA))
  row.names(test)=c('census_row', 'survey_row')
  test$con_1=NULL
  return(test)
}

##' Create a data lexicon for holding the associated column names
##'
##' @title addDataAssociation
##' @export
##' @param indf A data Lexicon (data.frame) created from the function:  \code{\link{createLexicon}}
##' @param data_names A vector vith two elements. The first element should be the name of 
##'     the \code{\link{census}} data column, and the second element should be the name of 
##'     the \code{\link{survey}} data column
##' @return indf The imported data lexicon with one extra column.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' data(survey)
##' data(census)
##' in.lexicon=createLexicon()
##' in.lexicon=addDataAssociation(in.lexicon, c("he","he"))
##' in.lexicon=addDataAssociation(in.lexicon, c("females","female"))
##' print(in.lexicon)
addDataAssociation=function(indf, data_names ){
  #checks=checkIfNamesInDataColumns(data_names)
  #if (checks==1){
    indf[[paste0("con_",ncol(indf) +1)]] = data_names
    #row.names(indf)=c('census_row', 'survey_row')
    return(indf)
  #} else if (checks==0) {
   # print(paste0("Cannot add this data association: ",paste(data_names, collapse=', ')))
   # stop()
  #}
}

##' Check the integrisy of the data Lexicon
##' 
##' @title checkIfNamesInDataColumns
##' @param names A vector with names to check if they exist as column names 
##'   in the data (census and survey)
##' @param incensus The census data
##' @param insurvey The survey data
##' @return anumber If both names are valid then it return '1' else 
##'     if the names are not valid data column names, it returns '0'.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
checkIfNamesInDataColumns=function(names, incensus, insurvey){
  census.name=names[1]
  survey.name=names[2]
  if (!(census.name %in% names(incensus))){
      cat(paste0("\n'",census.name,"' is not in the names of the census datset: ",paste(names(incensus), collapse=', '),"\n"))
      return(0)
      stop()
  }
  if (!(survey.name %in% names(insurvey))){
    cat(paste0("\n'",survey.name,"' is not in the names of the survey datset: ",paste(names(insurvey), collapse=', '),"\n"))
    return(0)
    stop()
  }
  return(1)
}

##' Check the lexicon data.frame
##'
##' @title check_lexicon
##' @export
##' @param inlex A data.frame which will be used a data lexicon 
##'   for listing the associated data columns.
##' @author Dimitris Kavroudakis \email{dimitris123@@gmail.com}
##' @examples library(sms)
##' df=createLexicon()
##' df=addDataAssociation(df, c("ena","duo"))
##' check_lexicon(df) 
check_lexicon=function(inlex){
  this_rows_names=row.names(inlex)
  correct_row_names=c('census_row', 'survey_row')
  
  #cat("\n----  Checking Lexicon  ----\n")
  
  if ( !identical(this_rows_names,correct_row_names) ){
    cat("\nNot correct row.names in your lexicon.
        Please rename the row.names of your lexicon as folows:
        row.names(mylexicon)= c('census_row', 'survey_row')\n\n") 
  return(-1)
  }
  else{ 
    #print("Data columns associations in the Lexicon are OK.")
    return(1)}
  
  #Check if he first row values exist in the column names of the census dataset
  
  #Check if the second row values exist in the column names of the survey dataset
  
  #if ( this_rows_names %%in%%  ){
  #  cat("\nNot correct row.names in your lexicon.
  #Please rename the row.names of your lexicon as folows:
  #row.names(mylexicon)= c('census_row', 'survey_row')\n\n") }
  
  #cat("\n----  End Checking Lexicon  ----\n\n")
}





