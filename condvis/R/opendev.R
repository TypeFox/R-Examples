# open a device suitable for interactivity with getGraphicsEvent
opendev <-
function (width = 7, height = 7)
{
    orig <- options("device")
    if (identical(.Platform$OS.type, "windows")){
        options(device = "windows")
    } else {
        options(device = "X11")
        if (identical(version$os, "linux-gnu")){
            X11.options(type = "Xlib")
        }
    }
    dev.new(width = width, height = height)
    options(orig)
}
