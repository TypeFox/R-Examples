##########################
##########################
##summary reports
##########################
##########################

#kr

#description
##########################
#functions to generate summary reports
#for pems objects, etc.


#includes 
##########################
#summaryReport
#


#to do
##########################
#

#comments
##########################
#


##########################
##########################
##summaryReport
##########################
##########################

#kr 23/01/2012 v 0.0.6

#what it does
##########################
#generates a summary of journey metrics


#to do
##########################
#make test more robust?

#document other report

#Average speed of entire drive cycle	km/h
#Average running speed	km/h
#Distance travelled	km
#Average accel	m/s2
#Average decel	m/s2
#Time percentage of idling	%
#Time percentage of accel: accel>0.1m/s2	%
#Time percentage of decel: accel <0.1m/s2)	%
#Time percentage in extremely low speed bracket [0-10]	%
#Average speed in extremely low speed bracket [0-10]	km/h
#Time percentage in very low speed bracket [10-20]	%
#Average speed in very low speed bracket [10-20]	km/h
#Time percentage in low speed bracket [20-50]	%
#Average speed in low speed bracket [20-50]	km/h
#Time percentage in moderate speed bracket [50-80]	%
#Average speed in moderate speed bracket [50-80]	km/h





#comments
##########################
#


#text for help would be very like 
#basic calcs


summaryReport <- function(speed = NULL, time = NULL, accel = NULL,  
                    distance = NULL, data = NULL, ..., 
                    lod.speed = 0.1, lod.accel = 0.1,
                    fun.name = "summaryReport", hijack= FALSE){
  
    #setup
    this.call <- match.call()
    
    #run checks

    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 

    if(!hijack){   
        speed <- checkInput(speed, data=data, if.missing = "return")  
        accel <- checkInput(accel, data=data, if.missing = "return")
        time <- checkInput(time, data=data, if.missing = "return")
        distance <- checkInput(distance, data=data, if.missing = "return")
    }

#######################
#suggestion 
#######################

   #could have lod test here for those not null?
   #then subsequent calcs, e.g. calcAccel(speed) would use this?


    if(!is.null(speed) & !is.null(lod.speed))
        speed[speed<lod.speed] <- 0

    #moved from later
    if(length(lod.accel)==1) 
        lod.accel <- c(-lod.accel, lod.accel)

    if(!is.null(accel) & !is.null(lod.accel))
        accel[accel<max(lod.accel) & accel>min(lod.accel)] <- 0



#######################

########################
#suggestion
########################

    #get what we need if incomplete
    #could make this a function 
    #tiggered by any(c(is.null(...), is.null(...),...))?

    if(is.null(speed) & is.null(accel) & is.null(time) &is.null(distance))
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "want time, speed and accel but insufficient inputs\n\t can make do with time and distance and work up", 
                           suggest = "add something I can work with to call", if.warning = NULL, 
                           fun.name = fun.name)

    if(is.null(distance)){
        if(is.null(time) | is.null(speed)){
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "want distance but insufficient inputs\n\t can make do with time and speed and work up", 
                           suggest = "add distance or time and speed to call", if.warning = NULL, 
                           fun.name = fun.name)
        } else {
            distance <- calcDistance(speed = speed, time = time, if.missing = settings$if.missing, 
                                     unit.convesions = settings$unit.conversions, hijack = TRUE)
        }
    }        

    if(is.null(speed)){
        if(is.null(time) | is.null(distance)){
            checkIfMissing(if.missing = settings$if.missing, 
                           reply = "want speed but insufficient inputs\n\t can make do with time and distance and work up", 
                           suggest = "add speed or time and distance to call", if.warning = NULL, 
                           fun.name = fun.name)
        } else {
            speed <- calcSpeed(distance = distance, time = time, if.missing = settings$if.missing, 
                               unit.convesions = settings$unit.conversions, hijack = TRUE)
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
                               unit.convesions = settings$unit.conversions, hijack = TRUE)
        }
    }

################################

    


    #make sure distance, speed and accel are in all right units for next bit

    #note: 
    #need this to be separate to bit before

    speed <- convertUnits(speed, to = "km/h", hijack = TRUE, if.missing = settings$if.missing, unit.convesions = settings$unit.conversions)
    distance <- convertUnits(distance, to = "km", hijack = TRUE, if.missing = settings$if.missing, unit.convesions = settings$unit.conversions)
    accel <- convertUnits(accel, to = "m/s/s", hijack = TRUE, if.missing = settings$if.missing, unit.convesions = settings$unit.conversions)

    #calculate time interval

    #note:
    #is this longwinded or safer than using time?
    
    if(is.null(time)){    
        d.time <- (distance / speed) * 3660 #interval in seconds
        temp <- mean(d.time, na.rm = TRUE)
        d.time <- ifelse(is.na(d.time), temp, d.time)
    } else {
        d.time <- diff(time)
        temp <- mean(d.time, na.rm = TRUE)
        d.time <- c(temp, ifelse(is.na(d.time), temp, d.time))
    }

    #calculations

#now earlier 

#    if(length(lod.accel)==1) 
#        lod.accel <- c(-lod.accel, lod.accel)

    distance.travelled.km <- sum(distance, na.rm = TRUE)
    time.total.s <- sum(d.time, na.rm = TRUE)

    avg.speed.km.h <- mean(speed, na.rm=TRUE)
    temp <- subset(speed, speed > lod.speed)
    avg.running.speed.km.h <- mean(temp, na.rm=TRUE)

    temp <- subset(d.time, speed < lod.speed)    
    time.idle.s <- sum(temp, na.rm = TRUE)    
    time.idle.pc <- (time.idle.s / time.total.s) * 100 

    temp <- subset(accel, accel > max(lod.accel))    
    avg.accel.m.s.s <- mean(temp, na.rm = TRUE)

    temp <- subset(d.time, accel > max(lod.accel))    
    time.accel.s <- sum(temp, na.rm = TRUE)    
    time.accel.pc <- (time.accel.s / time.total.s) * 100 
            
    temp <- subset(accel, accel < min(lod.accel))    
    avg.decel.m.s.s <- mean(temp, na.rm = TRUE)

    temp <- subset(d.time, accel < min(lod.accel))    
    time.decel.s <- sum(temp, na.rm = TRUE)    
    time.decel.pc <- (time.decel.s / time.total.s) * 100 

    #output

    data.frame(
               distance.travelled.km = distance.travelled.km,
               time.total.s = time.total.s,
               avg.speed.km.h = avg.speed.km.h, 
               avg.running.speed.km.h = avg.running.speed.km.h,     
               time.idle.s = time.idle.s,    
               time.idle.pc = time.idle.pc, 
               avg.accel.m.s.s = avg.accel.m.s.s,
               time.accel.s = time.accel.s,  
               time.accel.pc = time.accel.pc,            
               avg.decel.m.s.s = avg.decel.m.s.s, 
               time.decel.s = time.decel.s,    
               time.decel.pc = time.decel.pc 
    )
   
}



























