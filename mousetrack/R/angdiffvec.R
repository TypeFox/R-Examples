.packageName <- 'mousetrack'

angdiffvec <- function(xinfo,yinfo){
    
    allang = vector()

    for (j in 2:length(xinfo)){
        
        x1diff = xinfo[j] - xinfo[j-1]
        y1diff = yinfo[j] - yinfo[j-1]
        
        x2base = 0
        y2base = abs(y1diff)
        
        ang =  atan2(x1diff*y2base - y1diff*x2base,
                     x1diff*x2base + y1diff*y2base)*(180/pi)
        
        allang = c(allang , ang)
        
    }

    return(allang)

}
