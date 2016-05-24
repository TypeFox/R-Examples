## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

limit_value <- function(value,maxvalue,minvalue){
#Make sure that a value is set to be within a border (maxvalue and
#minvalue)

if((value>maxvalue)){
    returnvalue=maxvalue
} else {
    if((value<minvalue)){
        returnvalue=minvalue
    } else {
        returnvalue=value
    }
}
return( returnvalue) 
}
