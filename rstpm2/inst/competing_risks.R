## Stata data for comparison with stpm2cif
require(foreign)
require(rstpm2)

if (FALSE) {
    system("R CMD INSTALL ~/src/R/rstpm2")
    try(detach("package:rstpm2",unload=TRUE))
    library(rstpm2)
}

hivsi <- read.dta("http://www.stata-press.com/data/r12/hiv_si.dta")
if (test <- FALSE)
    hivsi <- transform(hivsi, status=ifelse(status=="event free",0,ifelse(status=="AIDS",1,2)))
crprep <- function(data,status,causes=NULL,cause="cause",event="event",refLevel=1) {
    if (is.null(causes)) {
        causes <- levels(as.factor(data[[status]]))
        warning(sprintf("Automagical reference level assumed to be '%s'.\n",causes[refLevel]))
        causes <- causes[-refLevel]
    }
    newdata <- do.call("rbind",lapply(causes,function(causek) {
        newdata <- data
        newdata[[cause]] <- causek
        newdata
    }))
    newdata[[cause]] <- factor(newdata[[cause]],levels=causes)
    if (is.factor(data[[status]])) 
        newdata[[event]] <- as.character(newdata[[status]]) == as.character(newdata[[cause]])
    else
        newdata[[event]] <- newdata[[status]] == newdata[[cause]]
    newdata
}
if (!test) {
    hivsi <- crprep(hivsi,"status")
    hivsi <- crprep(hivsi,"status",c("AIDS","SI")) ## new variables: cause, event
} else
    hivsi <- crprep(hivsi,"status",1:2)
fit <- stpm2(Surv(time,event) ~ 1, data=hivsi, # coxph.strata=as.factor(cause) # FAILS?
             tvc.formula=~as.factor(cause):nsx(log(time),df=3),
             baseoff=TRUE)
if (!test) {
    plot(fit,newdata=data.frame(cause="AIDS"))
    plot(fit,newdata=data.frame(cause="SI"),add=TRUE,ci=F,rug=F,line.col="blue")
} else {
    plot(fit,newdata=data.frame(cause=1))
    plot(fit,newdata=data.frame(cause=2),add=TRUE,ci=F,rug=F,line.col="blue")
}

stpm2cif <- function(obj,x,causeVar,causeVal) {
    causes <- levels(as.factor(obj@data[[causeVar]]))
    stpm2if <- function(x) {
        ## setup up newdata
        newdata <- data.frame(x=x)
        names(newdata) <- obj@timeVar
        ## calculate cause-specific survival
        Sk <- lapply(causes, function(cause) {
            newdata2 <- newdata
            newdata2[[causeVar]] <- cause
            predict(obj, newdata=newdata2)
        })
        S <- apply(do.call("cbind",Sk),1,prod)
        ## calculate cause-specific hazard
        newdata[[causeVar]] <- causeVal
        S*predict(obj,newdata=newdata,type="hazard")
    }
    sapply(x, function(xi)
           ifelse(xi==0,0,
                  integrate(stpm2if,0,xi)$value))
}
if (!test) {
    aidscif <- stpm2cif(fit,0:15,"cause","AIDS")
    sicif <- stpm2cif(fit,0:15,"cause","SI")
} else {
    aidscif <- stpm2cif(fit,0:15,"cause",1)
    sicif <- stpm2cif(fit,0:15,"cause",2)
}
matplot(0:15,cbind(aidscif,sicif),type="l",xlab="Time from diagnosis (years)",ylab="Cumulative probability")



## stpm2cif <- function(obj,x,causek,cause1,cause2,cause3=NULL,cause4=NULL) {
##     ## setup causes
##     causes <- c(cause1,cause2)
##     if (!is.null(cause3)) {
##         causes <- c(causes,cause3)
##     }
##     if (!is.null(cause4)) {
##         causes <- c(causes,cause4)
##     }
##     stpm2if <- function(x) {
##         ## setup up newdata
##         newdata <- data.frame(x=x)
##         names(newdata) <- obj@timeVar
##         for (cause in causes) newdata[[cause]] <- 0
##         ## calculate cause-specific survival
##         Sk <- lapply(causes, function(cause) {
##             newdata2 <- newdata
##             newdata2[[cause]] <- 1;
##             predict(obj, newdata=newdata2)
##         })
##         S <- apply(do.call("cbind",Sk),1,prod)
##         ## calculate cause-specific hazard
##         newdata[[causek]] <- 1
##         S*predict(obj,newdata=newdata,type="hazard")
##     }
##     sapply(x, function(xi)
##            ifelse(xi==0,0,
##                   integrate(stpm2if,0,xi)$value))
## }
## aidscif <- stpm2cif(fit,0:15,"aids","aids","si")
## sicif <- stpm2cif(fit,0:15,"si","aids","si")
## matplot(0:15,cbind(aidscif,sicif),type="l")

