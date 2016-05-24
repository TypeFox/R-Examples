see.ancova.tck <- function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check

    sigma <- 1
    assign("sigma", tclVar(sigma), envir = slider.env)
    n <- 10
    assign("n", tclVar(n), envir = slider.env)
    con.power <- 0.4
    assign("con.power", tclVar(con.power), envir = slider.env)
   
       
       see.anc <- function(mu1 = 10, mu2 = 12, mu3 = 16, sigma = 2, beta = 1, n = 10, var.exp = .4){ 
               
        Y1 <- rnorm(n, mu1, sigma)
        Y2 <- rnorm(n, mu2, sigma)
        Y3 <- rnorm(n, mu3, sigma)
        
        mu <- (mu1+mu2+mu3)/3
        alpha1 <- rep(mu1 - mu, n)
        alpha2 <- rep(mu2 - mu, n)
        alpha3 <- rep(mu3 - mu, n)
        alpha <- c(alpha1, alpha2, alpha3)
        
        dev.hold()
        layout(matrix(c(1,2,3,4,4,4),2,3, byrow = TRUE))
        par(mar = c(4,0,0,.5), oma=c(0,4.5,0.5,0.5))
        plot(1:n,Y1, xaxt = "n", ylab = "", xlab = "", ylim = c(5, 20), cex.lab = 1.3, cex.axis = 1.3, cex = 1.3); abline(h = mean(Y1), lty = 1); legend("top", legend = "Treatment 1", bty = "n", cex = 1.1) 
        mtext(side=2, line = 3, expression(italic(Y)))
        axis(1,at=seq(1,n), labels=FALSE)
        plot(1:n,Y2, xaxt = "n", yaxt = "n", ylab = "", xlab = "", pch = 19, ylim = c(5, 20), cex.lab = 1.3, cex.axis = 1.3, cex = 1.3); abline(h = mean(Y2), lty = 2); legend("top", legend = "Treatment 2", bty = "n", cex = 1.1); mtext(side = 1, line = 1.8, "Experimental unit", cex = .9) 
        axis(1,at=seq(1,n), labels=FALSE)
        plot(1:n,Y3, xaxt = "n", yaxt = "n", ylab = "", xlab = "", pch = 2, ylim = c(5, 20), cex.lab = 1.3, cex.axis = 1.3, cex = 1.3); abline(h = mean(Y3), lty = 3); legend("top", legend = "Treatment 3", bty = "n", cex = 1.1)  
        axis(1,at=seq(1,n), labels=FALSE)
        
        Xij <- runif(n * 3, 0, 10) 
        Yij <- mu + alpha + beta * (Xij - (mean(Xij))) + rnorm(n * 3, 0, sigma * (1 - var.exp))
        
        par(mar = c(4.6,0,0,0.5))
        plot(Xij - mean(Xij), Yij, pch = c(rep(1,n), rep(19, n), rep(2, n)), xlab = expression(italic(X-bar(X))), ylab = "", cex.lab = 1.4, cex.axis = 1.3, cex = 1.3)
        mtext(side=2, line = 3, expression(italic(Y)))
        abline(mu + alpha1, beta, lty = 1)
        abline(mu + alpha2, beta, lty = 2)
        abline(mu + alpha3, beta, lty = 3)
        dev.flush()
        }
        
norm.refresh <- function(...){
     sigma <- as.numeric(evalq(tclvalue(sigma), envir = slider.env))
     n <- as.numeric(evalq(tclvalue(n), envir = slider.env))   
     con.power <- as.numeric(evalq(tclvalue(con.power), envir = slider.env))
     
     see.anc(sigma = sigma, n = n, var.exp = con.power)
     }
     
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing ANCOVA")
    tkpack(tklabel(m, text = "      Visualizing ANCOVA      "))
    tkwm.geometry(m, "+0+0")
    
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C3", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "left")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.1, 
        to = 2, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = sigma), envir = slider.env)
  
   tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "left")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 2, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env) 
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Concomitant \n Exp. power", font = c("Helvetica", "9", 
        "normal"), width = "12"), side = "left")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 1, orient = "horiz", resolution = .1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = con.power), envir = slider.env)  
   
   
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}

                                                                                                                    



