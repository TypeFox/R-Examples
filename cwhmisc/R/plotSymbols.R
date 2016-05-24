plotSymbols <- function(interactive=FALSE) {
  interactive <- interactive && interactive()
  op <- options(warn=-2) 
  i <- 0:255
  j <- i
  j[c(27:32,129:256)] <- 47 # "."
  ncol <-16
  top <- 3 + 2*interactive
  opar <- par(cex.axis=0.7, mar=c(3,3,top,3)+0.1)
  on.exit(par(opar))

  plot(i%%ncol,1+i%/%ncol, pch=j-1, xlim=c(0,ncol-1), xlab="", ylab="", axes=FALSE)
  axis(1, at=0:15)
  axis(2, at=1:16, labels=0:15*16, las=2)
  axis(3, at=0:15)
  axis(4, at =1:16, labels=0:15*16+15, las=2)
  if (interactive) {
    title(main="Click on a symbol to add it to the resulting data frame.\n (Double) Click in margin to quit!", cex.main=0.8, line=3.5)
  }
  if (interactive) {
    df <- list()
    usr <- par("usr")
    ready <- FALSE
    while (!ready) {
      click <- locator(n=1)
      x <- click$x
      y <- click$y - 1
      ready <- (x < -0.25 || x > 16.25 || y < -0.25 || y > 16.25)
      if (!ready) {
        x <- round(x)
        y <- round(y)
        z <- 16*y + x
        ch  <- intToASCII(z)
        dec <- as.character(z) 
        hex <- intToHex(z)
        oct <- intToOct(z)
        spc <- paste(rep("0", 2-nchar(hex)), collapse="")
        hex <- paste(spc, hex, sep="")
        spc <- paste(rep("0", 3-nchar(oct)), collapse="")
        oct <- paste(spc, oct, sep="")
        df$ch  <- c(df$ch , ch )
        df$dec <- c(df$dec, dec)
        df$hex <- c(df$hex, hex)
        df$oct <- c(df$oct, oct)
        if (z == 0) ch <- " "
        spc <- paste(rep(" ", 3-nchar(dec)), collapse="")
        dec <- paste(spc, dec, sep="")
        cat("Selected ASCII character '", ch, "' ", dec, " 0x", hex, " \\", oct, "\n",sep="")
      }
    }
    return(df)
  }
  options(op)
  invisible()
}  ## plotSymbols()

availColors <- function (indx = 0:6)
{
    for (ii in unique(indx)) {
        is <- 100 * ii + 1:100
        if (min(is) > length(colors())) {
            cat("Maximum value of arg is", floor(length(colors())/100),
                "\n")
            return(NULL)
        }
        foo <- matrix(colors()[is], nrow = 10)
        par(mar = c(3, 3, 0.25, 0.25))
        plot(1:10, 1:10, type = "n", yaxt = "n", xlab = "", ylab = "")
        axis(2, at = 1:10, labels = 10:1)
        for (j in 1:10) {
            for (i in 1:10) {
                points(j, 11 - i, col = foo[i, j], pch = 16,
                  cex = 4)
                text(j, 11 - i - 0.3, foo[i, j], cex = 0.8)
            }
        }
        if (length(indx) > 1 & ii < max(indx))
            readline(paste("Currently showing group", ii, "  CR to continue "))
    }
    invisible(foo)
}  ## availColors()

plotSymbolsFonts <- function (fn=1) {
    i <- 0:255
    ncol <- 16
    opar <- par(cex.axis = 0.7, mar = c(3, 3, 3, 3) + 0.1)
    plot(i%%ncol, 1 + i%/%ncol, pch=i, font=fn, xlab = "", ylab = "", 
        axes = FALSE)
    axis(1, at = 0:15)
    axis(2, at = 1:16, labels = 0:15 * 16, las = 2)
    axis(3, at = 0:15)
    axis(4, at = 1:16, labels = 0:15 * 16 + 15, las = 2)
    par(opar)
}  ## plotSymbolsFonts()
