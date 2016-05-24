setMethod("show","CureRate",function(object){
### 1 - power table only (DEFAULT); 
### 2 - (everything) all tables, inclusing the power table; 
### 3 - required tables and power table.

    numac=length(object@actime)
    numfu=length(object@futime)
    numinfo=length(object@info)
    if(object@printflag == 2L) {object@indac=1:numac;object@indfu=1:numfu;}

    cat("Input:\n")
    cat(" numreps =",object@numreps,"\n cureobs =",object@cureobs,
        "\n curerx =",object@curerx,"\n medobs =",object@medobs,
        "\n medrx =",object@medrx, "\n acrate =",object@acrate,
        "\n probrx =",object@probrx, 
        "\n actime =", format(object@actime,digit=1,nsmall=1),
        "\n futime =", format(object@futime,digit=1,nsmall=1),
        "\n info =", format(object@info,digit=1,nsmall=2),
        "\n crits =", format(object@crits,digit=1,nsmall=4),
        "\n alpha =", object@alpha,"\n\n")
    

    if(object@printflag == 2L || object@printflag == 3L){
        if(object@printflag == 2L) cat("Full Results:\n")
        if(object@printflag == 3L) cat("Results:\n")
        for(i in object@indac){
            for(j in object@indfu){
                cat(" Total Eligible Patients: ",object@numobs[i,j],"\n")
                cat(" Accrual: ",object@actime[i],"   Followup: ",object@futime[j],"\n")
                cat(" Null Hypothesis: Cure Rate ",object@cureobs,"   Median Survival ",object@medobs,"\n")
                cat(" Alternative    : Cure Rate ",object@curerx,"   Median Survival ",object@medrx,"\n")    
                cat("",object@testname,"\n")

                outTable=matrix(cbind(object@timept[i,j,],object@info,
                    object@deaths[i,j,],object@crits,
                    object@power[i,j,]), 
                    nrow = length(object@info), ncol=5, byrow=FALSE,
                    dimnames=list(rep(c(""),length(object@info)),
                                  c("Real Time","Info Time","Deaths","Boundary","Rejection Prob")))
                print(outTable)
                cat("\n")
            }
        }
    }

    if(object@printflag == 1L){
        cat("Results:\n")
        cat("",object@testname,"\n\n")
    }

    if(object@printflag == 1L || object@printflag == 2L  || object@printflag == 3L){
        cat(" Power with information times {",format(object@info,digit=1,nsmall=2),"}\n")  
        cat("             Follow-up durations: \n Accrual       ")
        for(j in 1:numfu)
            cat(format(object@futime[j],digit=1, nsmall=1),"        ")
        cat("\n durations: \n")
        for(i in 1:numac){
            cat("  ", format(object@actime[i],digit=1, nsmall=1),"")
            for(j in 1:numfu) 
                cat("      ",format(object@beta[i,j],digit=3, nsmall=3))
            cat("\n")
        }
    }

}
)
