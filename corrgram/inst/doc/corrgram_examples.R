## ----setup, echo=FALSE, results="hide"--------------------------------------------------
library("knitr")
opts_chunk$set(fig.align="center", fig.width=7, fig.height=7)
options(width=90)

## ---------------------------------------------------------------------------------------
library("corrgram")

## ----fig2-------------------------------------------------------------------------------
vars2 <- c("Assists","Atbat","Errors","Hits","Homer","logSal",
           "Putouts","RBI","Runs","Walks","Years")
corrgram(baseball[,vars2], order=TRUE,
         main="Baseball data PC2/PC1 order",
         lower.panel=panel.shade, upper.panel=panel.pie,
         diag.panel=panel.minmax, text.panel=panel.txt)

## ----fig3-------------------------------------------------------------------------------
baseball.cor <- cor(baseball[,vars2], use='pair')
baseball.eig <- eigen(baseball.cor)$vectors[,1:2]
e1 <- baseball.eig[,1]
e2 <- baseball.eig[,2]
plot(e1,e2,col='white', xlim=range(e1,e2), ylim=range(e1,e2))
text(e1,e2, rownames(baseball.cor), cex=1)
arrows(0, 0, e1, e2, cex=0.5, col="red", length=0.1)

## ----fig4-------------------------------------------------------------------------------
corrgram(baseball[,vars2], main="Baseball data (alphabetic order)")

corrgram(baseball[,vars2], order=TRUE,
         main="Baseball data (PC order)",
         panel=panel.shade, text.panel=panel.txt)

## ----fig5-------------------------------------------------------------------------------
corrgram(baseball, order=TRUE, main="Baseball data (PC order)")

## ----fig6-------------------------------------------------------------------------------
corrgram(auto, order=TRUE, main="Auto data (PC order)")

## ----fig7-------------------------------------------------------------------------------
rinv <- function(r){
  # r is a correlation matrix
  # calculate r inverse and scale to correlation matrix
  # Derived from Michael Friendly's SAS code

  ri <- solve(r)
  s <- diag(ri)
  s <- diag(sqrt(1/s))
  ri <- s %*% ri %*% s
  n <- nrow(ri)
  ri <- ri * (2*rep(1,n) - matrix(1, n, n))
  diag(ri) <- 1  # Should already be 1, but could be 1 + epsilon
  colnames(ri) <- rownames(ri) <- rownames(r)
  return(ri)
}

vars7 <- c("Years", "logSal", "Homer", "Putouts", "RBI", "Walks",
           "Runs", "Hits", "Atbat", "Errors", "Assists")
cb <- cor(baseball[,vars7], use="pair")
corrgram(-rinv(cb), main=expression(paste("Baseball data ", R^-1)))


## ----fig8-------------------------------------------------------------------------------
require(Matrix) # For block diagonal function

partial <- function(r, xvar){
  # r is a correlation matrix
  # Calculate partial correlation of y|x
  yvar <- setdiff(colnames(r), xvar)
  ri <- r[yvar,yvar] - r[yvar,xvar] %*% solve(r[xvar,xvar]) %*% r[xvar,yvar]
  s <- diag(ri)
  s <- diag(sqrt(1/s))
  ri <- s %*% ri %*% s
  ri <- as.matrix(bdiag(ri, r[xvar, xvar]))
  diag(ri) <- 1  # Should already be 1, but could be 1 + epsilon
  colnames(ri) <- rownames(ri) <- c(yvar, xvar)
  return(ri)
}

vars8a <- c("Gratio", "Rep78", "Rep77", "Hroom", "Trunk", "Rseat",
            "Length", "Weight", "Displa", "Turn")
vars8b <- c("MPG", "Price")
vars8 <- c(vars8a, vars8b)
auto.cor <- cor(auto[, vars8], use="pair")
auto.par <- partial(auto.cor, vars8b)
corrgram(auto.par, lower.panel=panel.pie, upper.panel=panel.pie,
         main="Auto data, partialing out Price,MPG")

## ---------------------------------------------------------------------------------------
corrgram(baseball[,vars2], order=TRUE,
         main="Baseball correlation ellipses",
         panel=panel.ellipse,
         text.panel=panel.txt, diag.panel=panel.minmax)

## ----finish, echo=FALSE, results="asis"-------------------------------------------------
toLatex(sessionInfo(), locale=FALSE)

