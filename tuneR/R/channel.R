channel <- 
function(object, which = c("both", "left", "right", "mirror")){
    if(!is(object, "Wave")) 
        stop("Object not of class 'Wave'")
    validObject(object)
    which <- match.arg(which)
    return(
        switch(which,
            both = object,
            left = {
                object@stereo <- FALSE
                object@right <- numeric(0)
                object
            },
            right = {
                if(!object@stereo)
                    stop("Object is mono and does not contain a right channel")
                object@left <- object@right
                object@stereo <- FALSE
                object@right <- numeric(0)
                object
            },
            mirror = {
                if(!object@stereo)
                    stop("Object must be stereo in order to mirror")
                temp <- object@right
                object@right <- object@left
                object@left <- temp
                object
            }
        )
    )
}
