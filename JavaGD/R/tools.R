.javaGD.get.size <- function(devNr=dev.cur()) .Call(javaGDgetSize, devNr - 1)

.javaGD.copy.device <- function(devNr=dev.cur(), device=pdf, width.diff=0, height.diff=0, ...) {
    s<-.javaGD.get.size(devNr)
    pd<-dev.cur()
    dev.set(devNr)
    dev.copy(device, width=par()$din[1]/0.72+width.diff, height=par()$din[2]/0.72+height.diff, ...)
    dev.off()
    dev.set(pd)
    invisible(devNr)
}

.javaGD.version <- function() {
	v <- .Call(javaGDversion)
	list(major=v[1]%/%65536, minor=(v[1]%/%256)%%256, patch=(v[1]%%256), numeric=v[1])
}

.javaGD.set.display.parameters <- function(dpiX=100, dpiY=100, aspect=1)
	invisible(.Call(javaGDsetDisplayParam,c(dpiX, dpiY, aspect)))

.javaGD.get.display.parameters <- function()
    .Call(javaGDgetDisplayParam)

.javaGD.set.class.path <- function(cp) {
    cp <- as.character(cp)
    if (length(cp) < 1) stop("Invalid class path")
    invisible(.Call(setJavaGDClassPath, cp))
}

.javaGD.get.class.path <- function()
    .Call(getJavaGDClassPath)

.javaGD.resize <- function(devNr = dev.cur())
    invisible(.Call(javaGDresizeCall, devNr))
