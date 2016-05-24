### Unit tests of function pgig

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.pgig <- function()
{
    ## Purpose: Level 1 test of pgig
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: David Scott, Date:  8 Jan 2012

    ## compare with figures supplied by Slevinsky and Safouhi

    testData <- matrix(c(0.1, 1, 1, 0,
                         0.2, 1, 1, 0,
                         0.5, 1, 1, 0,
                         2, 2, 1, 1,
                         2, 2, 1, -1,
                         3, 2, 2, 2,
                         3, 2, 2, -2,
                         5, 3, 3, 4,
                         6, 3, 4, 5,
                         20, 1, 1, 0),
                       byrow = TRUE, ncol = 4)
    testResults <- c(0.9986941930676488,
                     0.9726685983580131,
                     0.7797339136097861,
                     0.6181518477661696,
                     0.1159938316382597,
                     0.3060636385450593,
                     0.1469821513077778*10^(-2),
                     0.8868482756840722*10^(-1),
                     0.1199275449037456*10^(-1),
                     0.4824320137970867*10^(-5))
    numTests <- NROW(testData)

    results <- numeric(numTests)

    ## compare with pgig using defaults
    for(i in 1:numTests){
        results[i] <- pgig(testData[i, 1], param = testData[i, 2:4],
                           lower.tail = FALSE)
    }

    results - testResults
    maxDiff <- max(abs(results - testResults))

    checkTrue(maxDiff < 10^(-11), msg =
              paste("pgig with defaults. maxDiff =", maxDiff))

    ## compare with pgig using machine tolerance
    for(i in 1:numTests){
        results[i] <- pgig(testData[i, 1], param = testData[i, 2:4],
                           lower.tail = FALSE, ibfTol = .Machine$double.eps)
    }

    results - testResults
    maxDiff <- max(abs(results - testResults))

    checkTrue(maxDiff < 10^(-13), msg =
              paste("pgig with machine tolerance. maxDiff =", maxDiff))

    ## compare with results obtained from numerical integration
    intResults <- numeric(numTests)
    for(i in 1:numTests){
        intResults[i] <- integrate(dgig, lower = testData[i, 1], upper = Inf,
                                   param = testData[i, 2:4])$value
    }

    intResults - testResults
    maxDiff <- max(abs(testResults - intResults))

    checkTrue(maxDiff < 10^(-9), msg =
              paste("test results compare with integrate. maxDiff =", maxDiff))

    results - intResults
    maxDiff <- max(abs(results - intResults))

    checkTrue(maxDiff < 10^(-9), msg =
              paste("pgig compared with integrate. maxDiff =", maxDiff))
}
