`Round.modelShow` <-
function(x, digits=6) {
 cat("Note: Function Valid Only For modelShow Class\n")
 cat("      For Other Classes Use round(x, digits) Instead.\n\n")
 cond     <- NULL
 for (i in 1:dim(x)[2]) if (is.numeric(x[,i]) == TRUE) cond <- c(cond, i)
 x[,cond] <- round(x[,cond], digits)
 print(x)
 invisible(x)
 }
# Round(essai,2)