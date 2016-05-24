"vec2pwm" <-
function(vec, as.list=FALSE) {
    if(as.list) {
       z <- list(BETA0 = vec[1], BETA1 = vec[2], BETA2 = vec[3],
                 BETA3 = vec[4], BETA4 = vec[5], source="vec2pwm")
    } else {
       z <- list(betas=vec, source="vec2pwm")
    }
    return(z)
}

