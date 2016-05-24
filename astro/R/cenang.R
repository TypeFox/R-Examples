cenang = function(a1, d1, a2, d2, units = "deg", method = "vincenty"){
    
    # http://en.wikipedia.org/wiki/Great-circle_distance#Formulas
    
    # conversion to radians?
    if(units == "deg"){
        a1 = a1 * (pi/180)
        a2 = a2 * (pi/180)
        d1 = d1 * (pi/180)
        d2 = d2 * (pi/180)
    }
    
    # calculation
    if(method == "sloc" | method == "s"){
        
        # spherical law of consines - good for large angles
        theta = acos( (sin(d1)*sin(d2)) + (cos(d1)*cos(d2)*cos(a1-a2)) )
        
    }else if(method == "haversine" | method == "h"){
        
        # haversine - good for most angles
        theta = 2 * asin( sqrt( (sin((d1-d2)/2)^2) + (cos(d1)*cos(d2)*((sin((a1-a2)/2))^2)) ) )
        
    }else{
        
        # vincenty - good for all angles
        numerator = sqrt( (( cos(d2)*sin(a1-a2) )^2) + (( (cos(d1)*sin(d2)) - (sin(d1)*cos(d2)*cos(a1-a2)) )^2) )
        denominator = ( sin(d1)*sin(d2) ) + ( cos(d1)*cos(d2)*cos(a1-a2) )
        theta = atan2(numerator,denominator)
    
    }
    
    # unit conversion?
    if(units == "deg"){
        theta = theta * (180/pi)
    }
    
    # return results
    return(theta)
    
}

