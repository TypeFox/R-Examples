##########################
##########################
##VSP calculations
##########################
##########################

#kr

#description
##########################
#functions to calculate VSP


#includes 
##########################
#calcVSP


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
#calculates VSP

#urgent
##########################
#need to tify the this.call
#handling give current
#not parent call!!!
#


#to do
##########################
#make test more robust?

#comments
##########################
#



calcVSP <- function(speed = NULL, accel = NULL, slope = NULL, 
                    time = NULL, distance = NULL, data = NULL,
                    calc.method = calcVSPJimenezPalaciosCMEM,
                    ..., fun.name = "calcVSP", this.call = NULL, 
                    hijack= FALSE){
  
    #setup
#temp fix
#think about using listUpdate in loa
    if(is.null(this.call)) 
        this.call <- match.call() 
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 
    if(!hijack){   
        speed <- checkInput(speed, data=data, if.missing = "return")  
        accel <- checkInput(accel, data=data, if.missing = "return")
        slope <- checkInput(slope, data=data, if.missing = "return")
        time <- checkInput(time, data=data, if.missing = "return")
        distance <- checkInput(distance, data=data, if.missing = "return")
    }

    if(is.null(speed) & is.null(accel) & is.null(time) &is.null(distance))
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "want speed and accel but insufficient inputs\n\t can make do with time and distance and work up", 
                           suggest = "add something I can work with to call", if.warning = NULL, 
                           fun.name = fun.name)
        
    if(is.null(speed)){
        if(is.null(time) | is.null(distance)){
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "want speed but insufficient inputs\n\t can make do with time and distance and work up", 
                           suggest = "add speed or time and distance to call", if.warning = NULL, 
                           fun.name = fun.name)
        } else {
            speed <- calcSpeed(distance = distance, time = time, if.missing = settings$if.missing, 
                               unit.conversions= settings$unit.conversions, hijack = TRUE)
        }
    }


    if(is.null(accel)){
        if(is.null(time) | is.null(speed)){
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "want accel but insufficient inputs\n\t can make do with time and distance or time and speed", 
                           suggest = "add speed and time or distance and time to call", if.warning = NULL, 
                           fun.name = fun.name)
        } else {
            accel <- calcAccel(speed = speed, time = time, if.missing = settings$if.missing, 
                               unit.conversions= settings$unit.conversions, hijack = TRUE)
        }
    }


###################################
#rearrange this bit so 
#it can pack here?
#not sure????
###################################


    if(is.function(calc.method))
        return(calc.method(speed = speed, accel = accel, slope = slope, data = data, 
                    ..., fun.name = fun.name, this.call = this.call, hijack= TRUE))

    #not good
    checkIfMissing(if.missing = settings$if.missing, 
                   reply = "could not run calc.method!", 
                           suggest = "check ?calcVSP if reason unclear", if.warning = "returning NULL", 
                           fun.name = fun.name)

    return(NULL)    
}



calcVSPJimenezPalaciosCMEM <- function(speed = NULL, accel = NULL, 
                    slope = NULL, m = NULL, a = NULL, b = NULL, 
                    c = NULL, g = NULL, ..., data = NULL,  
                    fun.name = "calcVSPJimenezPalaciosCMEM", 
                    this.call = NULL, hijack= FALSE){
  
    #setup
#temp fix
    if(is.null(this.call)) 
        this.call <- match.call() 

    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 
    if(!hijack){   
        speed <- checkInput(speed, data=data, if.missing = settings$if.missing)  
        accel <- checkInput(accel, data=data, if.missing = settings$if.missing)
        slope <- checkInput(slope, data=data, if.missing = "return")
    }

    if(is.null(speed) | is.null(accel))
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "Need speed and accel", 
                           suggest = "add speed and accel to see or see ?calcVSP", if.warning = NULL, 
                           fun.name = fun.name)
        
    if(is.null(slope)){
            checkIfMissing(if.missing = "warning", 
                           reply = "slope not supplied\n\t assuming 0", 
                           suggest = "add slope to call and rerun if required", if.warning = NULL, 
                           fun.name = fun.name)
            slope <- 0
    }
    
    #want specific units
    speed <- convertUnits(speed, to = "m/s", hijack = TRUE, unit.conversions = settings$unit.conversions, 
                          if.missing = settings$if.missing, fun.name = fun.name)  
    accel <- convertUnits(accel, to = "m/s/s", hijack = TRUE, unit.conversions = settings$unit.conversions, 
                          if.missing = settings$if.missing, fun.name = fun.name)  

#slope to sort out

    #make data always pems
    if(!isPEMS(data)) data <- makePEMS(data)

###################
#temp fix
###################
    old.class <- class(data)
    class(data) <- "not.pems"
       
    if(is.null(m)){        
        m <- if(is.null(data$constants$vsp.m))
                 1.5 else data$constants$vsp.m
    }

#note no special bus handling

    if(is.null(a)){        
        a <- if(is.null(data$constants$vsp.a)){
                    if(m < 3.855) 1.1
                       else if(m < 6.35) (0.0996 * m) / 2204.6 
                           else if(m < 14.968) (0.0875 * m) / 2204.6
                               else (0.0661 * m) / 2204.6 
                 } else data$constants$vsp.a
    }
    
    if(is.null(b)){        
        b <- if(is.null(data$constants$vsp.b))
                 0.132 else data$constants$vsp.b
    }

    if(is.null(c)){        
        c <- if(is.null(data$constants$vsp.c)){
                    if(m < 3.855) 0.000302
                       else if(m < 6.35) (1.47 + (5.2e-05 * m)) / 2205 
                           else if(m < 14.968) (1.93 + (5.90e-5 * m)) / 2205
                               else (2.89 + (4.21e-5 * m)) / 2205
                 } else data$constants$vsp.c
    }

    if(is.null(g)){        
        g <- if(is.null(data$constants$vsp.g))
                 9.81 else data$constants$vsp.g
    }


    vsp <- speed * (a * accel + (g * slope) + b) + (c * speed^3)

        this.vsp.units <- "kW/metric Ton"

    #my units
    vsp <- makePEMSElement(vsp, name="vsp", units="kW/metric ton")


########################
#temp fix
########################
    class(data) <- old.class

    #make output
    calcPack(output = vsp, data = data, settings = settings, 
             fun.name = fun.name, this.call = this.call) 

}








