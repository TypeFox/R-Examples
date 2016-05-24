nicetime = function(seconds){
    lapseconds = round(seconds)
    seconds = lapseconds%%60
    minutes = ((lapseconds - seconds)/60)%%60
    hours = ((lapseconds - minutes*60 - seconds)/3600)%%24
    days = ((lapseconds - hours*3600 - minutes*60 - seconds)/86400)
    lapline={}
    if(days!=0){
        if(days==1){lapline=paste(days,"d, ",sep="")   # day
        }else{lapline=paste(days,"d, ",sep="")         # days
        }
    }
    if(hours!=0 | days!=0){
        if(hours==1){lapline=paste(lapline,hours,"h, ",sep="")     # hour
        }else{lapline=paste(lapline,hours,"h, ",sep="")            # hours
        }
    }
    if(minutes!=0 | hours!=0 | days!=0){
        if(minutes==1){lapline=paste(lapline,minutes,"m, ",sep="")     # minute
        }else{lapline=paste(lapline,minutes,"m, ",sep="")              # minutes
        }
    }
    if(seconds==1){lapline=paste(lapline,seconds,"s",sep="")           # second
    }else{lapline=paste(lapline,seconds,"s",sep="")                    # seconds
    }
    return(lapline)
}

