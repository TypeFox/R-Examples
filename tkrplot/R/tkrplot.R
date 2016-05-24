if (Sys.info()["sysname"] == "Windows") {
    .my.tkdev <- function(hscale=1, vscale=1)
        win.metafile(width=4*hscale,height=4*vscale, restoreConsole=FALSE)
} else if (exists("X11", env=.GlobalEnv)) {
    .my.tkdev <- function(hscale=1, vscale=1)
        X11("XImage", 480*hscale, 480*vscale)
} else stop("tkrplot only supports Windows and X11")

.make.tkindex <-
    local({
        .My.Tk.index <- 0
        function() {
            .My.Tk.index <<- .My.Tk.index + 1
            .My.Tk.index
        }
    })

tkrplot <- function(parent, fun, hscale=1, vscale=1) {
    image <- paste("Rplot", .make.tkindex(), sep="")
    .my.tkdev(hscale, vscale)
    try(fun())
    .Tcl(paste("image create Rplot", image))
    lab<-tklabel(parent,image=image) #**** use try, delete image on failure
    tkbind(lab,"<Destroy>", function() .Tcl(paste("image delete", image)))
    lab$image <- image
    lab$fun <- fun
    lab$hscale <- hscale
    lab$vscale <- vscale
    lab
}

tkrreplot <- function(lab, fun = lab$fun,
                      hscale=lab$hscale, vscale=lab$vscale) {
    .my.tkdev(hscale, vscale)
    try(fun())
    .Tcl(paste("image create Rplot", lab$image))
}

#**** this is unspeakably crude
tkpersp <- function(x,y,z, theta = 30,phi = 30,expand = 0.5, r = sqrt(3), ...) {
    base<-tktoplevel()

    draw <- function() {
        par(bg = "white")
        try(persp(x, y, z, theta = theta, phi = phi, expand = expand, r = r,
                  ...))
    }

    img<-tkrplot(base, draw)

    getTclVar <- function(name)
        .Tcl(paste("set", name))
    setTclVar <- function(name, value)
        .Tcl(paste("set ", name, " {", value, "}", sep = ""))

    make.scale <- function(parent, from, to, resolution, title, command,
                           initial=from, showvalue =FALSE, orient="horiz") {
        var <- paste("rgv_", .make.tkindex(), sep="")
        setTclVar(var, initial)
        fun <- function(...) command(as.numeric(getTclVar(var)))

        frame <-tkframe(parent, relief="groove", borderwidth=2)
        tkpack(tklabel (frame, text=title))
        tkpack(tkscale(frame, command=fun, from=from, to=to,
                       showvalue=showvalue, variable=var,
                       resolution=resolution, orient=orient))
        tkbind(frame,"<Destroy>", function() .Tcl(paste("unset", var)))
        frame
    }

    frame <- tkframe(base)

    s.theta <- make.scale(frame, from=0, to=360, resolution=5,
                title="Theta", initial=theta, showvalue=TRUE,
                command=function(x) {
                    if (x != theta) {
                        theta <<- x
                        tkrreplot(img)
                    }
                })

    s.phi <- make.scale(frame, from=0, to=360, resolution=5,
                title="Phi", initial=phi, showvalue=TRUE,
                command=function(x) {
                    if (x != phi) {
                        phi <<- x
                        tkrreplot(img)
                    }
                })

    s.expand <- make.scale(frame, from=0.05, to=1, resolution=0.05,
                title="Expand", initial=expand, showvalue=TRUE,
                command=function(x) {
                    if (x != expand) {
                        expand <<- x
                        tkrreplot(img)
                    }
                })

    s.r <- make.scale(frame, from=0.05, to=3, resolution=0.05,
                title="R", initial=r, showvalue=TRUE,
                command=function(x) {
                    if (x != r) {
                        r <<- x
                        tkrreplot(img)
                    }
                })

    tkpack(s.theta, s.phi, s.expand, s.r, side="left", anchor="n")
    tkpack(img,frame)
}

.onLoad <- function(libname, pkgname) {
    chname <- "tkrplot"
    file.ext <- .Platform$dynlib.ext
    dlname <- paste(chname, file.ext, sep = "")
    if (is.character(.Platform$r_arch) && .Platform$r_arch != "")
        path <- file.path("libs", .Platform$r_arc, dlname)
    else path <- file.path("libs", dlname)
    file <- system.file(path, package = pkgname, lib.loc = libname)[1]
    tryCatch(tcl("load", file, "Rplot"),
             error = function(e)
             warning("loading Rplot failed", call. = FALSE))
}
