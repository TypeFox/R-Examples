stereo <- 
function(left, right){
    if(missing(right))
        return(left@stereo)
    ln <- deparse(substitute(left))
    rn <- deparse(substitute(right))
    if(missing(right)) 
        right <- left
    else 
        equalWave(left, right)
    if(left@stereo)
        stop(ln, " already is a stereo 'Wave' object")
    if(right@stereo)
        stop(rn, " already is a stereo 'Wave' object")
    if(length(right@left) != length(left@left))
        stop("Channel length of ", ln, " and ", rn, " differ")
    left@stereo <- TRUE
    left@right <- right@left
    return(left)
}
