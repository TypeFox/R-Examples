##########################
##########################
##above all else code 
##(hidden here because 
## first .r code file 
## alphabetically)
##########################
##########################

#example
#utils::globalVariables(c("a", "othervar"))

globalVariables(c("pems.scheme", "ref.unit.conversions", "ref.chem", "ref.diesel", "ref.petrol",
                  "panel.surfaceSmooth"))

#panel.surfaceSmooth only needed until new version of loa gets on CRAN...


##########################
##########################
##emission calculations
##########################
##########################

#kr

#description
##########################
#functions to calculate emissions


#includes 
##########################
#calcEm
#calcEMHoribaPitot
#


#to do
##########################

#comments
##########################
#



##########################
##########################
##calcVSP
##########################
##########################

#kr 23/01/2012 v 0.0.6

#what it does
##########################
#calculates Emissions

#most urgent
######################################
#confirm inputs to calcEmHoribaPitot
##what exhaust parameters are needed? 
##what delay offset is being applied to conc.co2
##round up/round down?
#

#urgent
##########################
#need to tidy the this.call
#handling give current
#not parent call!!!
##########################
#need to tidy calcEmHoribaPitot
#error catchers to add
#tempGet to be tidied 
#fix to be tidied
#




#to do
##########################
#make test more robust?

#comments
##########################
#



################################
################################
##calcEm
################################
################################

calcEm <- function(conc = NULL, time = NULL, 
                    calc.method = calcEmHoribaPitot, analyte = NULL,
                    ..., data = NULL, fun.name = "calcEm", force = FALSE, 
                    this.call = NULL, hijack= FALSE){
  
    #setup
#temp fix
#think about using listUpdate in loa
    if(is.null(this.call)) 
        this.call <- match.call() 
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 
    if(!hijack){   
        conc <- checkInput(conc, data=data, if.missing = "stop", fun.name = fun.name)  
        time <- checkInput(time, data=data, if.missing = "return", fun.name = fun.name)
    }

    temp <- attr(conc, "name")
    if(!force){
        if(length(grep("^conc.", temp))<1) 
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = paste("'", temp, "' should be concentration, \n\tdenoted conc.[analyte]", sep=""), 
                           suggest = "select suitable input or force if sure ?calcEm", if.warning = NULL, 
                           fun.name = fun.name)
    }
    temp <- gsub("^conc.", "", temp)
    if(is.null(analyte))
        analyte <- temp
    if(!force){
        if(temp != analyte)
            checkIfMissing(if.missing = settings$if.missing, 
                           reply =  "Input type does not match assigned 'analyte'", 
                           suggest = "select suitable input or force if sure ?calcEm", if.warning = NULL, 
                           fun.name = fun.name)
    }


    if(is.function(calc.method))
        return(calc.method(conc = conc, data = data, analyte = analyte,
                    ..., this.call = this.call, force = force, hijack= TRUE))

    checkIfMissing(if.missing = settings$if.missing, 
                   reply = "could not run calc.method!", 
                           suggest = "check ?calcEm if reason unclear", if.warning = "returning NULL", 
                           fun.name = fun.name)
    return(NULL)    
}





##################################
##################################
##calcEmHoribaPitot
##################################
##################################



#various error catchers needed


calcEmHoribaPitot <- function(conc = NULL, time = "local.time",
                    exflow = "exh.flow.rate", extemp = "exh.temp", 
                    express = "exh.press", analyte = NULL, delay = NULL, mm = NULL,  
                    ..., force = force, data = NULL, fun.name = "calcEmHoribaPitot", 
                    this.call = NULL, hijack= FALSE){
  
    #setup
#temp fix
#think about using listUpdate in loa
    if(is.null(this.call)) 
        this.call <- match.call() 
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 
    if(!hijack){   
        conc <- checkInput(conc, data=data, if.missing = "return")  
    }

    time <- checkInput(time, data=data, if.missing = "return")
    exflow <- checkInput(exflow, data=data, if.missing = "return") 
    extemp <- checkInput(extemp, data=data, if.missing = "return") 
    express <- checkInput(express, data=data, if.missing = "return")

#    if(!force & !is.null(time)){
        #check that time is equidistant
        

#    }


#could pull this out as function for next few?
#note some check ref.chem, some don't

tempGet <- function(..., id=NULL, data=NULL, default=NULL){

           extra.args <- list(...)
           if(id %in% names(extra.args)) return(extra.args[[id]])
           if(!is.null(data)){
              if(isPEMS(data))  
                  if(id %in% names(pemsConstants(data)))
                      return(pemsConstants(data)[[id]])
              if(id %in% names(data)) return(data[[id]])
           }
           if(!is.null(default)){
              if(id %in% names(default)) return(default[[id]])
           }
           return(NULL)
} 

##########################################
#need to rethink mm.hc options
#maybe extend this to calcVSP etc?
##########################################

     if(is.null(mm))
         mm <- tempGet(..., id=paste("mm.", analyte[1], sep=""), data=data, default=ref.chem)
     if(is.null(mm) & analyte[1]=="hc"){
         mm.h <- tempGet(..., id="mm.h", data=data, default=ref.chem)
         alpha.exhaust.hc  <- tempGet(..., id="alpha.exhaust.hc", data=data, default=ref.chem)
         mm.c <- tempGet(..., id="mm.c", data=data, default=ref.chem)
         temp <- mm.h * alpha.exhaust.hc + mm.c
         thc.c6 <- tempGet(..., id="thc.c6", data=data, default=ref.chem)
         mm <- temp * thc.c6

#################################
#error catchers for missing cases
#################################

     }
     if(is.null(delay))
         delay <- tempGet(..., id=paste("delay.", analyte[1], sep=""), data=data)

###################################
#to do
###################################

    #units of all inputs
    conc <- convertUnits(conc, to = "vol%", hijack = TRUE, 
                         if.missing = settings$if.missing, 
                         unit.convesions = settings$unit.conversions)
    exflow <- convertUnits(exflow, to = "L/min", hijack = TRUE, 
                           if.missing = settings$if.missing, 
                           unit.convesions = settings$unit.conversions)

    if(delay[1]>=1){
        conc <- c(conc[(floor(delay[1])+1): length(conc)], rep(NA, floor(delay[1])))
    }

    em <- conc * mm * exflow * (1/60) * (1/100) * (1/22.415)  * (273.15/293.15)

    #my structure
    em <- makePEMSElement(em, name=paste("em.", analyte[1], sep=""), units="g/s")

    
########################
#temp fix
########################
    #class(data) <- old.class

    #make output
    calcPack(output = em, data = data, settings = settings, 
             fun.name = fun.name, this.call = this.call) 

}
