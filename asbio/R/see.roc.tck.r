see.roc.tck <- function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    mu1 <- 0.4
    assign("mu1", tclVar(mu1), envir = slider.env)
    mu2 <- 0.6
    assign("mu2", tclVar(mu2), envir = slider.env)
    c <- 0.5
    assign("c", tclVar(c), envir = slider.env)
      
    dev.new(width=10,height=5)
    
    xlim <- c(0, 1.8)
    ylim <- c(0, dnorm(0, 0, 0.07))
    
    norm.refresh <- function(...) {
        dev.hold()
        
        mu1 <- as.numeric(evalq(tclvalue(mu1), envir = slider.env))
        mu2 <- as.numeric(evalq(tclvalue(mu2), envir = slider.env))
        c <- as.numeric(evalq(tclvalue(c), envir = slider.env))
        sigma = 0.08
        layout(matrix(c(1,1,1,1,1,2,2), 1, 7))
        par(mar=c(5,4.1,3,1))
        
        
        x<-NULL; rm(x)
        
        cols <- c(rgb(red=1,blue=1,green=1,alpha=.6), rgb(red=0,blue=0,green=0,alpha=.6), rgb(red=.5,blue=.5,green=0.5,alpha=.6))
        curve(dnorm(x,mu1,sigma),from=0,to=1,ylab = "",ylim=ylim,xlab = "",xlim=xlim, main = "Dichotomous populations",cex.lab=1.4,cex.axis=1.3, cex.main=1.4,xaxt="n",yaxt="n")
        axis(side = 1, at = c(0, .2, .4, .6, .8,1))
        mtext(side = 1, at = 0.5, line = 3, expression(italic(c)))
        mtext(side = 2, line = 1.7, "Density")
        legend("topleft", pch=22, pt.cex=1.7, pt.bg=c(cols[1],cols[2]), legend = c(expression(paste(italic(Y)," = 0")), expression(paste(italic(Y)," = 1"))), bty = "n")
                
        xx <- seq(-0.5, 1.5, length = 500)
        yy1 <- dnorm(xx, mu1, sigma)
        yy2 <- dnorm(xx, mu2, sigma)
                       
        polygon(c(xx[xx<=1.5],1.5),c(yy2[xx<=1.5],yy2[xx==-0.5]),col=cols[2])
        polygon(c(xx[xx<=1.5],1.5),c(yy1[xx<=1.5],yy1[xx==-0.5]),col=cols[1])
        
       
        ys <- c(ylim[2]*.6, ylim[2]*.8, ylim[2])
        xs <- c(.6*xlim[2], .8*xlim[2], xlim[2])
       
        rect(xs[1],ys[2],xs[2],ys[3], col=cols[1])   #topleft
        rect(xs[2],ys[2],xs[3],ys[3], col=cols[3])   #topright
        rect(xs[1],ys[1],xs[2],ys[2], col=cols[3])   #bottomleft
        rect(xs[2],ys[1],xs[3],ys[2], col=cols[2])   #bottomright
       
        xv <- c(rnorm(200, mu1, sigma), rnorm(200, mu2, sigma))
        yv <- c(rep(0, 200), rep(1, 200))
        
        glmm <- suppressWarnings(glm(yv ~ xv, family = "binomial"))
        abline(v = c, lwd = 2)
        
        t <- table(ifelse(fitted(glmm) > c, 1, 0), yv)
        TP <- round(t[,2][2]/sum(t[,2]),2)
        FN <- round(1 - t[,1][1]/sum(t[,1]),2)
        TN <- round(t[,1][1]/sum(t[,1]),2)
        FP <- round(1 - t[,2][2]/sum(t[,2]),2)
        
        text(mean(c(xs[1],xs[2])), .97*ys[3], "True negative rate")
        text(mean(c(xs[2],xs[3])), .97*ys[3], "False positive rate")
        text(mean(c(xs[1],xs[2])), .97*ys[2], "False negative rate")
        text(mean(c(xs[2],xs[3])), .97*ys[2], "True positive rate")
        
        text(mean(c(xs[1],xs[2])), .9*ys[3], bquote(.(TN)))
        text(mean(c(xs[2],xs[3])), .9*ys[3], bquote(.(FP)))
        text(mean(c(xs[1],xs[2])), .86*ys[2], bquote(.(FN)))
        text(mean(c(xs[2],xs[3])), .86*ys[2], bquote(.(TP)))
        
        int <- seq(0.01, 0.99, by =.01)
        sens <- rep(NA, length(int))
        omspec <- rep(NA, length(int))
        
        
        
        for(i in 1:length(int)){
        t <- table(ifelse(fitted(glmm) > int[i], 1, 0), yv)
        if(ncol(t)==1){
            if(colnames(t)=="0") t <- cbind(t, c(0,0))
            if(colnames(t)=="1") t <- cbind(c(0,0), t)
                      }
        if(nrow(t)==1){
            if(rownames(t)=="0") t <- rbind(t, c(0,0))
            if(rownames(t)=="1") t <- rbind(c(0,0), t)
                      }
                                    
        sens[i] <- t[,2][2]/sum(t[,2])
        omspec[i] <- 1 - t[,1][1]/sum(t[,1])
        }
        
        plot(c(1,omspec,0), c(1,sens,0), xlab = "False positive rate", ylab = "True positive rate", type = "l", xlim=c(0,1), ylim=c(0,1), lwd = 2.5, main = "ROC curve",cex.lab=1.4,cex.main=1.4,cex.axis=1.3)
        abline(0, 1, lty = 2)
        
        dev.flush()
    }
    
    
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing ROC curves")
    tkpack(tklabel(m, text = "      Visualizing ROC curves      "))
    tkwm.geometry(m, "+0+0")
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Y = 0", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.3, 
        to = 0.48, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu1), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Y = 1", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.54, 
        to = 0.7, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu2), envir = slider.env)
    
 tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "c", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.2, 
        to = 0.95, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = c), envir = slider.env)
    
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
                                                                     