r.dist <- function(rho, r, n){

rd <- function(rho, r, n){
((n-2)*gamma(n-1)*((1-rho^2)^((n-1)/2))*((1-r^2)^((n-4)/2)))/
(sqrt(2*pi)*gamma(n-0.5)*((1-rho*r)^(n-(3/2)))) *
(1 + 1/4 * (((rho * r) + 1)/((2 * n) - 1))  + 9/16 * ((((rho * r) + 1)^2)/((2 * n) - 1)* (2 * n + 1)))  
}
r.unst <- rd(rho, r, n)
integrand <-function(x)rd(rho = rho, n = n, r = x)
den <- integrate(integrand, -1, 1)
std.r <- r.unst/den$value
std.r
}


#----------------------------- tcltk implementation ---------------------------------#


see.r.dist.tck <-function(){
plot_r <- function(rho, r, n, ...){
vals <- r.dist(rho = rho, r = seq(-1, 1, .001) , n = n)
dev.hold()
plot(seq(-1, 1, .001), vals, type = "l", ylab = expression(paste(italic(f),"(",italic(r),")", sep = "")), xlab = expression(italic(r)), ...)
dev.flush()
}


    if (!exists("slider.env")) 
        slider.env <- NULL
    suppressWarnings(rm(slider.env))
    slider.env <<- new.env()
    n <- 5
    rho <- 0
    assign("n", tclVar(n), envir = slider.env)
    assign("rho", tclVar(rho), envir = slider.env)
    xmin <- -1
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 1
    assign("xmax", tclVar(xmax), envir = slider.env)
    
    norm.refresh <- function(...) {
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        rho <- as.numeric(evalq(tclvalue(rho), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
  
        plot_r(n = n, rho = rho, xlim = c(xmin, xmax))
    }
   
   
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.geometry(m, "+600+4")
    tkwm.title(m, "Pearson product moment correlation distribution")
    tkpack(tklabel(m, text = "Sampling distribution of r"))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 5, 
        to = 50, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C1", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -.99, 
        to = .99, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = rho), envir = slider.env)
 
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir = slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir = slider.env)
 
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
