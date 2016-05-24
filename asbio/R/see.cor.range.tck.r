see.cor.range.tck <-function(sd = .5){
draw.cor.range <- function(y1, y2, a = 0, b = 2, col.arg = TRUE){
dev.hold()
plot(y1, y2, xlab=expression(italic(Y)[1]),ylab=expression(italic(Y)[2]), type = "n")
if(col.arg==TRUE) col = rgb(red=0.4,blue=0.4,green=0.5,alpha=.9)
if(col.arg==FALSE) col = rgb(red=0.5,blue=0.5,green=0.5)

ny1 <- c(y1[y1 >= a & y1 <= b])
ny2 <- c(y2[y1 >= a & y1 <= b])
rect(a,-10, b, 10, col=col)
points(y1, y2)
points(ny1, ny2, pch = 19)

legend("topleft", bty = "n", title = expression(paste("     ",italic(r),"        ",italic(r[k]), "       ",italic(r[s]))), legend = c(paste(round(cor(y1, y2), 2),round(cor(y1, y2, method = "k"), 2), round(cor(y1, y2, method = "s"), 2)),
paste(round(cor(ny1, ny2), 2),round(cor(ny1, ny2, method = "k"),2), round(cor(ny1, ny2, method = "s"),2))), pch = c(1, 19))
dev.flush()
}

y1 <- runif(50,-3,3)
y2 <- y1
for(i in 1 : length(y1)){
y2[i] <- y1[i] + rnorm(1,sd = sd)}


    if (!exists("slider.env"))
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    a <- 0
    b <- 1
    assign("a", tclVar(a), envir= slider.env)
    assign("b", tclVar(b), envir= slider.env)
    color<-tclVar(1)
    norm.refresh <- function(...) {
        a <- as.numeric(evalq(tclvalue(a), envir= slider.env))
        b <- as.numeric(evalq(tclvalue(b), envir= slider.env))
        color <- as.logical(tclObj(color))
        draw.cor.range(y1 = y1, y2 = y2, as.numeric(a), as.numeric(b), col.arg = color)
        }
        
        tclServiceMode(TRUE)
    m <- tktoplevel()
    col.cbut <- tkcheckbutton(m, text="Color", variable=color)
    tkwm.title(m, "The effect of range on correlation")
    tkpack(tklabel(m, text = "      Effect of range on correlation      "))
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "lower range", font = c("Helvetica",
        "9"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -4,
        to = 0, orient = "horiz", resolution = 0.1, showvalue = TRUE),
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = a), envir= slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "upper range", font = c("Helvetica",
        "9"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0,
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE),
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = b), envir= slider.env)
    
    
    tkpack(col.cbut)
        }