.image <- function(device) {
    a <- .Call("get_img_backplane", device, PACKAGE="Cairo")
    names(a) <- c('ref', 'info')
    a$width <- a[[2]][1]
    a$height <- a[[2]][2]
    a$format <- c("ARGB","RGB","A8","A1","dep","RGB16")[a[[2]][3]+1]
    class(a) <- "CairoImageRef"
    a
}

.ptr.to.raw <- function(ptr, begin, length)
    .Call("ptr_to_raw", ptr, begin, length, PACKAGE="Cairo")

.raw.to.ptr <- function(ptr, offset=0, raw, begin=0, length=length(raw))
    invisible(.Call("raw_to_ptr", ptr, offset, raw, begin, length, PACKAGE="Cairo"))

