## Tests if the calibration via "cal..." returns a threshold that
## gives the desired property by checking the correspoinding property
## of the chart. Currently only consider Shewhart and CUSUM charts.

context("Consistency of Calibration")

test_that("Calibration not possible CUSUM",{
    chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1))
    q <- getq(chart,property="calhitprob",params=list(target=0.9,nsteps=100))
    P <- list(mu=-2,sd=1)
    xi <- list(mu=0,sd=1)
    expect_error(q$q(P,xi),"not possible")
})



testofcalibration <- function(chart,X,type="ARL",params,allowerror=FALSE){
    testcal <- function(P,xi){
        qcal <- getq(chart,property=paste("cal",type,sep=""),params=params)
        tryCatch({
            threshold <- qcal$trafo(qcal$q(P=P,xi=xi))
            propraw=getq(chart,property=type,params=c(threshold=threshold,params))
            propraw$trafo(propraw$q(P=P,xi=xi))},
            error=function(e) {
                expect_match(e$message,"not possible")
                if (allowerror){
                    return(params$target)
                }else{
                    stop(e);
                }
            }
                 )
    }
    P <- chart@model$Pofdata(X)
    xi <- chart@model$xiofP(P)
    res <- max(abs((c(testcal(P,xi),
                      replicate(10,{
                          Presample <- chart@model$Pofdata(chart@model$resample(P))
                          testcal(Presample,xi)
                      }))-params$target)/params$target))
    id <- paste(class(chart)[1],type,paste(params,collapse=" "),"Maximum relative error:",format(res))
    expect_true(res<1e-3,id)
}

test_that("Calibration CUSUM ARL",{
    skip_on_cran()
    for (samplesize in c(5,100,10000)){
        X <- rnorm(samplesize)
        chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1))
        testofcalibration(params=list(target=100),type="ARL",chart=chart,X=X,allowerror=TRUE)
        testofcalibration(params=list(target=1000),type="ARL",chart=chart,X=X,allowerror=TRUE)
    }
}
)

test_that("Calibration CUSUM hitprob",{
    skip_on_cran()
    set.seed(123)
    for (samplesize in c(5,100,10000)){
        X <- rnorm(samplesize)
        chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1))
        testofcalibration(params=list(target=0.9,nsteps=1e2),chart=chart,X=X,type="hitprob",allowerror=TRUE)
        testofcalibration(params=list(target=0.1,nsteps=1e2),chart=chart,X=X,type="hitprob",allowerror=TRUE)
        testofcalibration(params=list(target=0.01,nsteps=1e2),chart=chart,X=X,type="hitprob",allowerror=TRUE)
        testofcalibration(params=list(target=0.001,nsteps=1e2),chart=chart,X=X,type="hitprob",allowerror=TRUE)
        testofcalibration(params=list(target=0.01,nsteps=1000),chart=chart,X=X,type="hitprob",allowerror=TRUE)
    }
}
)

test_that("Calibration Shewhart ARL",{
    skip_on_cran()
    for (samplesize in c(5,100,10000)){
        X <- rnorm(samplesize)
        chartShew <- new("SPCShew",model=SPCModelNormal())
        testofcalibration(params=list(target=100),chart=chartShew,X=X,type="ARL")
        testofcalibration(params=list(target=1000),chart=chartShew,X=X,type="ARL")
    }
}
)

test_that("Calibration Shewhart Hitprob",{
    skip_on_cran()
    for (samplesize in c(5,100,10000)){
        X <- rnorm(samplesize)
        chartShew <- new("SPCShew",model=SPCModelNormal())
        testofcalibration(params=list(target=0.01,nsteps=1e2),chart=chartShew,X=X,type="hitprob")
        testofcalibration(params=list(target=0.01,nsteps=1000),chart=chartShew,X=X,type="hitprob")
    }
}
)
