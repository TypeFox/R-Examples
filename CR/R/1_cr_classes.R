setClass("CureRate",
    representation(
        cureobs  = "numeric",
        medobs   = "numeric",
        curerx   = "numeric",
        medrx    = "numeric",
        actime   = "vector",
        futime   = "vector",
        info     = "vector",
        crits    = "vector",
        alpha    = "numeric",
        rho      = "numeric",

        acrate   = "numeric",
        probrx   = "numeric",
        numreps  = "integer", 
       
        numobs   = "matrix",
        timept   = "array",
        deaths   = "array",
        testname = "character",
        power    = "array",
        beta     = "matrix",

        indac    =  "integer",
        indfu    =  "integer",
        printflag = "integer"
    ) 
)

