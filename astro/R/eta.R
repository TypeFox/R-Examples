eta = function(i, total, start){
    #start = proc.time()[3]
    c = (total-i)+1
    if(c!=total){
        timenow = proc.time()[3]
        timeper = (timenow-start)/i
        timeleft = c*timeper*1.1 # 1.1 uncertainty factor
        eta = nicetime(timeleft)
    }else{
        eta = ""
    }
    return(eta)
}

