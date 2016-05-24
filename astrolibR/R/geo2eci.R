#source('ct2lst.R')

 geo2eci = function(incoord,jdtim){
        if(!is.matrix(incoord)) incoord = matrix(incoord,1)
        re=6378.137     # earth's equatorial radius, in km
        lat = incoord[,1]*pi/180.
        lon = incoord[,2]*pi/180.
        alt = incoord[,3]
        jdtime= jdtim
        
        gst = ct2lst(0,0,jdtime)
        angle_sid=gst*2.*pi/24.       # sidereal angle
        theta=lon+angle_sid             # azimuth
        r=(alt+re)*cos(lat)
        x=r*cos(theta)
        y=r*sin(theta)
        z=(alt+re)*sin(lat)
                
        return(c(x,y,z))
}
