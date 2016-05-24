
library(gridGraphics)

text1 <- function() {
    plot(-1:1, -1:1, type = "n", xlab = "Re", ylab = "Im")
    K <- 16; text(exp(1i * 2 * pi * (1:K) / K), col = 2)
}

text2 <- function() {
    ## The following two examples use latin1 characters: these may not
    ## appear correctly (or be omitted entirely).
    plot(1:10, 1:10, main = "text(...) examples\n~~~~~~~~~~~~~~",
         sub = "R is GNU ©, but not ® ...")
    mtext("«Latin-1 accented chars»: éè øØ å<Å æ<Æ", side = 3)
    points(c(6,2), c(2,1), pch = 3, cex = 4, col = "red")
    text(6, 2, "the text is CENTERED around (x,y) = (6,2) by default",
         cex = .8)
    text(2, 1, "or Left/Bottom - JUSTIFIED at (2,1) by 'adj = c(0,0)'",
         adj = c(0,0))
    text(4, 9, expression(hat(beta) == (X^t * X)^{-1} * X^t * y))
    text(4, 8.4, "expression(hat(beta) == (X^t * X)^{-1} * X^t * y)",
         cex = .75)
    text(4, 7, expression(bar(x) == sum(frac(x[i], n), i==1, n)))
    
    ## Two more latin1 examples
    text(5, 10.2,
         "Le français, c'est façile: Règles, Liberté, Egalité, Fraternité...")
    text(5, 9.8,
         "Jetz no chli züritüütsch: (noch ein bißchen Zürcher deutsch)")
}

# Test 'family'
text3 <- function() {
    plot(1:3, type="n")
    families <- c("sans", "serif", "mono")
    for (i in 1:3)
        text(i, i, "test", family=families[i])
}

# Test 'vfont'
text4 <- function() {
    plot(1, type="n")
    text(1, 1, "test", vfont=c("serif", "plain"))
}

plotdiff(expression(text1()), "text-1")
plotdiff(expression(text2()), "text-2")
plotdiff(expression(text3()), "text-3")
plotdiff(expression(text4()), "text-4")

plotdiffResult()

