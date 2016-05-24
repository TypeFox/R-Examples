"are.pargld.valid" <-
function(para,verbose=FALSE,nowarn=FALSE) {
    if(! is.gld(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    La2 <- para$para[2]
    La3 <- para$para[3]
    La4 <- para$para[4]

    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)

    if(verbose == TRUE) cat("Checking by Theorem 1.3.33\n")

    # The regional tests below are applicable if and only if
    # the following condition is met--this ensures that Lambda2
    # is suitable for the GLD.
    for(F in seq(0,1,by=0.00001)) {
      tmp <- La2 * (La3*F^(La3-1) +
                    La4*(1-F)^(La4-1))
      if(is.nan(tmp)) next; # figure that 0^(-1) is seen on end points
      if(tmp < 0) {
        warning("Parameters are invalid by Theorem 1.3.33 of Karian and Dudewicz (2000)")
        options(op)
        return(FALSE)
      }
    }

    if(verbose == TRUE) cat("Checking by region\n")
    ratio <- -1/La3
    # See Theorem 1.3.33 of Karian and Dudewicz and figure1.3-1
    if(La3 <= -1 && La4 >=  1) {         # REGION 1
       warning("Parameters are valid but L-moments are not")
       GO <- FALSE # ordinary L-moments do not exist in REGION 1
       #return(TRUE)  # However the GLD is valid
    }
    else if(La3 >=  1 && La4 <= -1) {    # REGION 2
       warning("Parameters are valid but L-moments are not")
       GO <- FALSE # ordinary L-moments do not exist in REGION 2
       #return(TRUE) # However the GLD is valid
    }
    else if(La3 >= 0 && La4 >= 0) {      # REGION 3
       GO <- TRUE
    }
    else if(La3 <= 0 && La4 <= 0) {      # REGION 4
       GO <- TRUE
    }
    else if((La4 <= ratio && La4 >= -1) && La3 >= 1) {  # REGION 6
       GO <- TRUE
    }
    else if(La4 >= ratio && (La3 >= -1 && La3 <= 0)) { # REGION 5
       GO <- TRUE
    }
    if(GO) {
      options(op)
      return(TRUE)
    }
    warning("Parameters are invalid by figure 1.3-1 of Karian and Dudewicz (2000)")
    options(op)
    return(FALSE)
}
