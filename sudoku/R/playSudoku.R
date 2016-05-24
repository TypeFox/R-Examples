## Original by Greg Snow <Greg.Snow@intermountainmail.org>

playSudoku <- function(z=NULL, hist.len=100, solve=TRUE,
                        display=c("guess","windows","tk"),
                        hscale=1.25, vscale=1.25, ...) {

  dsp <- substring(match.arg(display), 1,1)
  if (dsp=="g") dsp <- switch(.Platform$OS.type, windows="w", "t")
  if (dsp=="t" && !require(tkrplot)) stop("'tkrplot' package needed\n")
  
  if (identical(z,0)) {z <- matrix(0, 9,9); solve <- FALSE}
  if (is.null(z))      z <- generateSudoku(...)
  if (length(z)==1)    z <- readSudoku(z)
  if (solve) {cat("Solving..."); zz <- solveSudoku(z, print.it=FALSE); cat("done!\n")}
  cols <- ifelse(z, "blue","black")

  hst <- list(z)                   # Keep a history of z's to length "hist.len"
  ah <- function(newz) {hst <<- c(hst, list(newz))
                        if (length(hst) > hist.len) hst <<- hst[-1]}
    
  cusr <- cplt <- rep(0+NA, 4)
  replot <- function() {
    par(mar=c(0,0,0,0), bg="white")
    plot(0.5:9.5, 0.5:9.5, type="n", axes=FALSE, xlab="", ylab="")
    cusr <<- par("usr"); cplt <<- par("plt")
    segments(0.5:9.5, rep(0.5,10), 0.5:9.5, rep(9.5,10), col="grey")
    segments(rep(0.5,10), 0.5:9.5, rep(9.5,10), 0.5:9.5, col="grey")
    segments(c(0,3,6,9)+0.5, rep(0.5,4), c(0,3,6,9)+0.5, rep(9.5,4), lwd=3)
    segments(rep(0.5,4), c(0,3,6,9)+0.5, rep(9.5,4), c(0,3,6,9)+0.5, lwd=3)
    for (i in 1:9) for (j in 1:9) if (z[i,j]) {
      if (cols[i,j]=="red") text(j, 10-i, "X", col="pink", cex=3)
      text(j, 10-i, z[i,j], col=cols[i,j], font=ifelse(cols[i,j]=="blue",2,1),
           cex=ifelse(cols[i,j]=="blue", 2.0, 1.8))
    }
  }

  if (dsp=="t") {
    tt <- tktoplevel()
    tkwm.title(tt,"Sudoku")
    img <- tkrplot(tt, replot, hscale=hscale, vscale=vscale)
    txt <- tktext(tt, bg="white", font="courier")
    scr <- tkscrollbar(tt, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    tkconfigure(txt, yscrollcommand=function(...)tkset(scr,...))
    tkpack(img, side='top')
    tkpack(txt, side="left", fill="both", expand=TRUE)
    tkpack(scr, side="right", fill="y")
    iw <- as.numeric(tcl('image','width', tkcget(img,'-image')))
    ih <- as.numeric(tcl('image','height',tkcget(img,'-image')))
  }

  showz <- function() switch(dsp, w=replot(), t=tkrreplot(img))
  showz()
  
  cc <- function(x, y) {           # Convert mouse position to cell coordinates
    if (dsp=="t") {x <- (as.double(x)-1)/iw;  y <- 1 - (as.double(y)-1)/ih}
    px <- (x-cplt[1])/(cplt[2]-cplt[1])
    py <- (y-cplt[3])/(cplt[4]-cplt[3])
    ux <- px*(cusr[2]-cusr[1])+cusr[1]
    uy <- py*(cusr[4]-cusr[3])+cusr[3]
    c(10-round(uy), round(ux))
  }
  
  help.txt <- paste(" ?     -- this help",
                    "1-9   -- insert digit",
                    "0,' ' -- clear cell",
                    "r     -- replot the puzzle",
                    "q     -- quit",
                    "h     -- hint/help",
                    "c     -- correct wrong entries (show in red)",
                    "u     -- undo last entry",
                    "s     -- show number in cell",
                    "a     -- show all (solve the puzzle)",
                    "\n", sep="\n")
  type <- function(s) switch(dsp, w=cat(s),
                                  t={tkinsert(txt,'end',s); tksee(txt,'end')})
  ij <- c(5,5)                                                # Initial "point"
  mm.w <- function(buttons, x, y) {ij <<- cc(x,y); return()}
  mm.t <- function(x, y)          {ij <<- cc(x,y); return()}

  kb <- function(A) {
    i <- ij[1];  j <- ij[2]
    z[cols=="red"] <<- 0;  cols[cols=="red"] <<- "black"
    key <- switch(A, " "="0", "/"="?", tolower(A))
    if (key=="q") switch(dsp, t=tkdestroy(tt), w=return(1))
    if (key %in% c(0:9,"h","s") && (i < 1 || i > 9 || j < 1 || j > 9))
      {type("Must be over puzzle cell\n"); return()}
    if (key %in% c("c","s","a") && !solve)
      {type("Solution not available\n"); return()}
    if (key %in% c(0:9,"c","s","a")) ah(z)
    if (key %in% 0:9) {z[i,j] <<- as.double(key);  cols[i,j] <<- "black"}
    if (key=="?") type(help.txt)
    if (key=="h") type(hintSudoku(z, i,j))
    if (key=="c") {cols[z != 0 & z != zz] <<- "red"
                   if (!any(cols=="red")) {type("All Correct\n"); return()}}
    if (key=="u") {h <- length(hst); z <<- hst[[h]]; if (h>1) hst <<- hst[-h]}
    if (key=="s") {z[i,j] <<- zz[i,j];  cols[i,j] <<- "green3"}
    if (key=="a") {cols[z != zz] <<- "green3";  z <<- zz}
    if (key %in% c(0:9,"r","c","u","s","a")) showz()
    if (solve && all(z==zz)) type("You got it!\n")
    return()
  }    

  kb("?")
  if (solve && is.null(zz)) {type("Puzzle not solvable.\n"); solve <- FALSE}
  switch(dsp, w=getGraphicsEvent("Ready!", onMouseMove=mm.w, onKeybd=kb),
              t={tkbind(img,'<Motion>',mm.t); tkbind(tt,'<Key>',kb);
                 tkwait.window(tt)})
  return(invisible(z))
}
