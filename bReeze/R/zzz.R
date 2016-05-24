.onAttach <- 
function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(" ")
    packageStartupMessage(paste("This is", pkgname, ver))
    packageStartupMessage(" ")
    packageStartupMessage("Type changes(\"bReeze\") to see changes/bug fixes, help(bReeze) for documentation")
    packageStartupMessage("or citation(\"bReeze\") for how to cite bReeze.")
    packageStartupMessage(" ")
}


# show NEWS file
changes <- 
function(pkg="bReeze") {
    if(pkg=="bReeze") file.show(file.path(system.file(package="bReeze"), "NEWS"))
}


# short name wrapper functions
avail <- 
function(mast, v.set, dir.set, subset, digits=1, print=TRUE) availability(mast, v.set, dir.set, subset, digits, print)

cln <- 
function(mast, set, v.avg.min=0.4, v.avg.max=50, dir.clean=TRUE, turb.clean=4, icing=FALSE, rep=NULL, n.rep=5) clean(mast, set, v.avg.min, v.avg.max, dir.clean, turb.clean, icing, rep, n.rep)

en <- 
function(wb, rho=1.225, bins=c(5,10,15,20), digits=0, print=TRUE) energy(wb, rho, bins, digits, print)

freq <- 
function(mast, v.set, dir.set, num.sectors=12, bins=c(5,10,15,20), subset, digits=3, print=TRUE) frequency(mast, v.set, dir.set, num.sectors, bins, subset, digits, print)

ts <- 
function(timestamp, pattern, tz) timestamp(timestamp, pattern, tz)

map <- 
function(mast, type=c("satellite", "terrain", "hybrid", "roadmap"), zoom, label, ...) map.plot(mast, type, zoom, label, ...)

ms <- 
function(mast, set, signal="v.avg", fun=c("mean", "median", "min", "max", "sd"), subset, digits=3, print=TRUE) month.stats(mast, set, signal, fun, subset, digits, print)

day <- 
function(mast, set, dir.set=set, signal, num.sectors=NULL, subset, ...) day.plot(mast, set, dir.set, signal, num.sectors, subset, ...)

pol <- 
function(mast, v.set=1, dir.set=1, subset, ...) polar.plot(mast, v.set, dir.set, subset, ...)

iec <- 
function(mast, set, subset, ...) turb.iec.plot(mast, set, subset, ...)

pro <- 
function(mast, v.set, dir.set, num.sectors=12, method=c("hellman", "loglm", "fixed"), alpha=NULL, subset, digits=3, print=TRUE) profile(mast, v.set, dir.set, num.sectors, method, alpha, subset, digits, print)

turb <- 
function(mast, turb.set, dir.set, num.sectors=12, bins=c(5,10,15,20), subset, digits=3, print=TRUE) turbulence(mast, turb.set, dir.set, num.sectors, bins, subset, digits, print)

uc <- 
function(aep, uc.values, uc.names, prob=seq(5,95,5), digits=c(0,0), print=TRUE) uncertainty(aep, uc.values, uc.names, prob, digits, print)

wb <- 
function(mast, v.set, dir.set, num.sectors=12, subset, digits=3, print=TRUE) weibull(mast, v.set, dir.set, num.sectors, subset, digits, print)
