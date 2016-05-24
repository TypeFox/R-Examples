# TODO: Generates 34 testfunctions, one for each testcase taken from the Roche Diagnostics method comparison
#       intranet module testcase collection. Results obtained from the approximative Passing-Bablok implementation
#       (PaBaLarge) are compared to results of the exact implementation (PaBa).
# 
# Author: schueta6
###############################################################################



# one could increase the number of bins (NBins) to minimize the differences, here we use the default setting of NBins=1e06

PaBaLargePrecision <- 1e-04                                                       

cat("\n\n********************************************\nmcPaBaLarge.R method comparison test cases\n********************************************")

load(".\\TestCaseCollection\\testcases.RData", .GlobalEnv)

TCnames <- names(testcases)

# 

genericPaBaLargeTest <- function(Data, Name, Exception=FALSE)
{
    X <- Data[,1]
    Y <- Data[,2]
    
    if(Exception)
    {
        cat("\n\nTestcase:", Name, "\n")
        NData <- nrow(Data)
        cat("\nN(Data)         =", NData)
        NUData <- nrow(unique(Data))
        cat("\nN(unique(Data)) =", NUData)
        cat("\nN(Ties)         =", NData-NUData)
        cat("\n\n")
        checkException(mcreg(X, Y, method.reg="PaBaLarge", method.ci="analytical",NBins=1e06)@para)
        checkException(mcreg(X, Y, method.reg="PaBa",      method.ci="analytical")@para)
    }
    else
    {
        cat("\n\nTestcase: ", Name, "\n")
        NData <- nrow(Data)
        cat("\nN(Data)         =", NData)
        NUData <- nrow(unique(Data))
        cat("\nN(unique(Data)) =", NUData)
        cat("\nN(Ties)         =", NData-NUData)
        cat("\n\n")
        resPaBaL <- mcreg(X, Y, method.reg="PaBaLarge", method.ci="analytical", NBins=1e06)@para
        resExact <- mcreg(X, Y, method.reg="PaBa",      method.ci="analytical")@para
        
        checkEquals(resPaBaL, resExact, tolerance=PaBaLargePrecision)
    }    
}

# imitate call-by-value argument passing

cloneLocalArgs <- function(Data, Name, Exception)
{
    cloneData <- Data                           # generate local copies
    cloneName <- Name
    cloneExpt <- Exception
    locFunc <- function(){genericPaBaLargeTest(cloneData, cloneName, cloneExpt)}
    return(locFunc) 
}

# generate test-function for each dataset of the testcase collection

for( i in 1:length(testcases))
{
    Fname <- paste("test.PaBaLargeTestcase_", TCnames[i], sep="")

    assign( Fname, cloneLocalArgs(testcases[[i]], TCnames[i], TCnames[i] %in% c("part_1_dataset_2", "part_1_dataset_12")) )
}

