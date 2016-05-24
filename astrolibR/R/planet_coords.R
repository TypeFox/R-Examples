# source('helio.R')
# source('euler.R')
# source('juldate.R')

planet_coords = function( date, planet=planet, jd = FALSE) {

    radeg = 180/pi
    c = 2.99792458e5

    if(jd){
        jj = date
        if(length(jj)>0 && length(planet)>1 )
            stop('A planet name must be supplied for vector dates')
    }
    else {
        jj = juldate(date)
        jj = jj + 2400000
    }
    if(!missing(planet) ){
        planetlist = c('mercury','venus','mars',
        'jupiter','saturn','uranus','neptune','pluto')
        index = which(planetlist==tolower(planet))
        if(length(index)==0 ) stop(paste('unrecognized planet:' + planet))
        afterEarth = (index>=4 )
        index[afterEarth] = index[afterEarth] + 1
    }
    else {
        index = c(1,2,4,5,6,7,8,9)
    }

    tmp = helio(jj,index,radian=TRUE)
    rad = tmp$hrad
    lon = tmp$hlong
    lat = tmp$hlat

    tmp = helio(jj,3,radian=TRUE)
    rade = as.numeric(tmp$hrad)
    lone = as.numeric(tmp$hlong)
    late = as.numeric(tmp$hlat)

    browser("End")
    x = rad * cos(lat) * cos(lon) - rade * cos(late) * cos(lone)
    y = rad * cos(lat) * sin(lon) - rade * cos(late) * sin(lone)
    z = rad * sin(lat)            - rade * sin(late)
    lambda = atan2(y,x) * radeg
    beta   = atan2(z,sqrt(x*x + y*y)) * radeg

    tmp = euler(lambda, beta, 4)
    return(list(ra=tmp$ao, dec=tmp$bo))
}
