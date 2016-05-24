"Ops.spray" <- function (e1, e2 = NULL) 
{
    oddfunc <- function(...){stop("odd---neither argument has class spray?")}
    unary <- nargs() == 1
    lclass <- nchar(.Method[1]) > 0
    rclass <- !unary && (nchar(.Method[2]) > 0)
    
    if (!is.element(.Generic, c("+", "-", "*", "/", "^", "=="))){
        stop("operator '", .Generic, "' is not implemented for sprays")
    }

    if(unary){
        if (.Generic == "+") {
            return(e1)
        } else if (.Generic == "-") {
            return(spray_negative(e1))
        } else {
            stop("Unary operator '", .Generic, "' is not implemented for sprays")
        }
    }
    
    if (.Generic == "*") {
        if (lclass && rclass) {
            return(spray_times_spray(e1, e2))
        } else if (lclass) {
            return(spray_times_scalar(e1, e2))
        } else if (rclass) {
            return(spray_times_scalar(e2, e1))
        } else {
            oddfunc()
        }
    } else if (.Generic == "+") {
        if (lclass && rclass) {
            return(spray_plus_spray(e1, e2))  # S1+S2
        } else if (lclass) {
            return(spray_plus_scalar(e1, e2)) # S+x
        } else if (rclass) {
            return(spray_plus_scalar(e2, e1)) # x+S
        } else {
            oddfunc()
        }
    } else if (.Generic == "-") {
        if (lclass && rclass) {
            return(spray_plus_spray(e1, spray_negative(e2)))  # S1-S2
        } else if (lclass) {
            return(spray_plus_scalar(e1, -e2))                # S-x
        } else if (rclass) {
            return(spray_plus_scalar(spray_negative(e2), e1)) # x-S
        } else {
            oddfunc()
        }
    } else if (.Generic == "^") {
        if(lclass && !rclass){
            return(spray_power_scalar(e1,e2)) # S^n
        } else {
            stop("Generic '^' not implemented in this case")
        }
    } else if (.Generic == "==") {
        if(lclass && rclass){
            return(spray_eq_spray(e1,e2))
        } else {
            stop("Generic '==' only compares two spray objects with one another")
        }          
    } else if (.Generic == "/") {
        if(lclass && !rclass){
            return(spray_times_scalar(e1,1/e2))
          } else {
        stop("don't use '/', use ooom() instead")
    }
  }
}

spray_negative <- function(S){
    if(is.zero(S)){
        return(S)
    } else {
        return(spray(index(S),-value(S)))
    }
}

spray_times_spray <- function(S1,S2){
    if(is.zero(S1) || is.zero(S2)){return(zero)}
    stopifnot(arity(S1) == arity(S2))
    spraymaker(spray_mult(
        index(S1),value(S1),
        index(S2),value(S2)
        ))
}

spray_times_scalar <- function(S,x){
    stopifnot(length(x)==1)
    return(spray(index(S), x*value(S)))
}

spray_plus_spray <- function(S1,S2){
    if(is.zero(S1)){
        return(S2)
    } else if(is.zero(S2)){
        return(S1)
      }

    stopifnot(arity(S1)==arity(S2))
    return(spraymaker(spray_add(
    index(S1),value(S1),
    index(S2),value(S2)
        )))
}

spray_plus_scalar <- function(S,x){
    spray_plus_spray(S, x*one(S))
}

spray_power_scalar <- function(S,n){
    stopifnot(n==round(n))
    if(n<0){
        stop("use ooom() for negative powers")
    } else if(n==0){
        return(constant(S))
    } else if(n==1){
        return(S)
    } else {
        return(spray_times_spray(S,Recall(S,n-1)))   # the meat
    }
}

`spray_eq_spray` <- function(S1,S2){
    if(arity(S1) != arity(S2)){
        return(FALSE)
    } else if(length(S1) != length(S2)){
        return(FALSE)
    } else {
        return(spray_equality(index(S1),value(S1),index(S2),value(S2)))
    }
}
