##########################
##########################
##common calculations
##########################
##########################

#kr

#description
##########################
#functions to do common calculates


#includes 
##########################
#calcDistance


#to do
##########################

#comments
##########################
#template to think about




##########################
##########################
##calcDistance
##########################
##########################

#kr 23/01/2012 v 0.0.6

#what it does
##########################
#calculates distance travelled 
#since previous measurement


#to do
##########################
#make test more robust?

#URGENT


#comments
##########################
#


#############################
#############################
##calcChecks, calcPack
#############################
#############################

#front end management

#changed data = data to data = null

calcChecks <- function(fun.name = "calcChecks", ..., data = NULL,
                   if.missing = c("stop", "warning", "return"), 
                   output = c("input", "data.frame", "pems", "special"),
                   unit.conversions = NULL, overwrite = FALSE){

    #output handling
    output <- checkOption(output[1], formals(setUnits)$output, 
                       "output", "allowed outputs", 
                       fun.name = fun.name)

    if(output == "special"){
        output <- if(is.null(data))
                      "input" else if(comment(isPEMS(data)) == "other")
                                       "input" else comment(isPEMS(data))
    }

    #if.missing handling
    if.missing <- checkOption(if.missing[1], formals(setUnits)$if.missing, 
                              "if.missing", "allowed if.missings", 
                              fun.name = fun.name)

    list(output = output, if.missing = if.missing, overwrite = overwrite, 
         unit.conversions = unit.conversions)

}


calcPack <- function(output = NULL, data = NULL, settings = NULL, 
               fun.name = "calcPack", this.call = NULL){

    #make output
    output <- checkOutput(input = output, data = data, if.missing = settings$if.missing, 
                          fun.name = fun.name, output = settings$output, 
                          overwrite = settings$overwrite)
    if(isPEMS(output)){
        old.class <- class(output)
        class(output) <- "not.pems"
        output$history[[length(output$history)+1]] <- this.call
        class(output) <- old.class 
    }


    output
    
}





#############################
#############################
##calcDistance
#############################
#############################


calcDistance <- function(speed = NULL, time = NULL, data = NULL,
                     ..., fun.name = "calcDistance", hijack= FALSE){
    
    #setup
    this.call <- match.call()

    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #speed
    if(!hijack)   
        speed <- checkInput(speed, data=data)
    if(is.null(speed))
         checkIfMissing(settings$if.missing, reply = "argument 'speed' not supplied or null")
    speed <- convertUnits(speed, to = "m/s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #time
    if(!hijack)   
        time <- checkInput(time, data=data)
    if(is.null(time))
         checkIfMissing(settings$if.missing, reply = "argument 'time' not supplied or null")
    time <- convertUnits(time, to = "s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #my assumption
    #first unit resolution is average of rest
    #rest are time.now - time.last

    temp <- diff(time)
    temp <- c(mean(temp, na.rm=TRUE), temp)

    #my calculation
    distance <- speed * temp

    #my structure
    distance <- makePEMSElement(distance, name="distance", units="m")

    calcPack(output = distance, data = data, settings = settings, 
               fun.name = fun.name, this.call = this.call) 
   
}








#############################
#############################
##calcSpeed
#############################
#############################


calcSpeed <- function(distance = NULL, time = NULL, data = NULL,
                     ..., fun.name = "calcSpeed", hijack= FALSE){
    
    #setup
    this.call <- match.call()

    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #distance
    if(!hijack)   
        distance <- checkInput(distance, data=data)
    if(is.null(distance))
         checkIfMissing(settings$if.missing, reply = "argument 'distance' not supplied or null")
    distance <- convertUnits(speed, to = "m", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #time
    if(!hijack)   
        time <- checkInput(time, data=data)
    if(is.null(time))
         checkIfMissing(settings$if.missing, reply = "argument 'time' not supplied or null")
    time <- convertUnits(time, to = "s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #my assumption
    #first unit resolution is average of rest
    #rest are time.now - time.last

    temp <- diff(time)
    temp <- c(mean(temp, na.rm=TRUE), temp)

    #my calculation
    speed <- distance / temp

    #my structure
    speed <- makePEMSElement(speed, name="speed", units="m/s")

    calcPack(output = distance, data = data, settings = settings, 
               fun.name = fun.name, this.call = this.call) 
   
}






###########################
###########################
##calcAccel
###########################
###########################


calcAcceleration <- function(...) calcAccel(...)


calcAccel <- function(speed = NULL, time = NULL, data = NULL, 
                     ...,fun.name = "calcAccel", hijack= FALSE){
    
    #setup
    this.call <- match.call()
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)


    #speed
    if(!hijack)   
        speed <- checkInput(speed, data=data)
    if(is.null(speed))
         checkIfMissing(settings$if.missing, reply = "argument 'speed' not supplied or null")
    speed <- convertUnits(speed, to = "m/s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #time
    if(!hijack)   
        time <- checkInput(time, data=data)
    if(is.null(time))
         checkIfMissing(settings$if.missing, reply = "argument 'time' not supplied or null")
    time <- convertUnits(time, to = "s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #my assumption
    #first d.speed/d.time is 0
    #rest are ...now - ....last

    d.speed <- diff(speed)
    d.time <- diff(time)

    #my calculation
    accel <- c(0, d.speed / d.time)

    #my structure
    accel <- makePEMSElement(accel, name="accel", units="m/s/s")

    #make output
    calcPack(output = accel, data = data, settings = settings, 
               fun.name = fun.name, this.call = this.call) 
    
}







############################
############################
##calcJerk
############################
############################



calcJerk <- function(accel = NULL, time = NULL, data = NULL, 
                     ...,fun.name = "calcJerk", hijack= FALSE){
    
    #setup
    this.call <- match.call()
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)


    #accel
    if(!hijack)   
        accel <- checkInput(accel, data=data)
    if(is.null(accel))
         checkIfMissing(settings$if.missing, reply = "argument 'accel' not supplied or null")
    accel <- convertUnits(accel, to = "m/s/s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #time
    if(!hijack)   
        time <- checkInput(time, data=data)
    if(is.null(time))
         checkIfMissing(settings$if.missing, reply = "argument 'time' not supplied or null")
    time <- convertUnits(time, to = "s", fun.name = fun.name, 
                       unit.conversions = settings$unit.conversions, hijack = TRUE)

    #my assumption
    #first d.accel/d.time is 0
    #rest are ...now - ....last

    d.accel <- diff(accel)
    d.time <- diff(time)

    #my calculation
    jerk <- c(0, d.accel / d.time)

    #my units
    attr(jerk, "name") <- "jerk"
    attr(jerk, "units") <- "m/s/s/s"
    #my structure
    jerk <- makePEMSElement(jerk, name="jerk", units="m/s/s/s")

    #make output
    calcPack(output = jerk, data = data, settings = settings, 
               fun.name = fun.name, this.call = this.call) 
    
}
