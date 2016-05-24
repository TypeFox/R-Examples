`PCC` <-
function(theta=0, S=0, C=0, D=0, s=1/1.702, b=seq(-5,5,length=300), c=0, d=1,
         groups=TRUE, ID="ID", main="Person Characteristic Curve",
         xlab="Item Difficulty Parameter (b)", ylab="P(x = 1)", type=c("g","a") ) {
   if (length(b) <= 1) stop("The length of the item parameters vector must be greater than 1")
   if (length(s) != 1 && length(s) != length(b)) stop("Item parameters not well defined")
   if (length(c) != 1 && length(c) != length(b)) stop("Item parameters not well defined")
   if (length(d) != 1 && length(d) != length(b)) stop("Item parameters not well defined")
   rep            <- max(length(theta),length(S),length(C),length(D))
   if (is.null(ID)) ID <- "NULL"
   subMain <- switch(ID,
    ID    = "\n (Subject ID)",
    ALL   = "\n (THETA, S, C, D)",
    THETA = "\n (THETA)",
    NULL  = "")
   main <- paste(main,subMain)

   prob           <- data.frame(matrix(NA, ncol=10, nrow=length(b)*rep))
   colnames(prob) <- c("ID","b","s","c","d","THETA","S","C","D","p")

   prob$b <- rep(b,rep)
   if (length(s) == 1) s1 <- rep(s,length(b))
   if (length(c) == 1) c1 <- rep(c,length(b))
   if (length(d) == 1) d1 <- rep(d,length(b))
   prob$s     <- rep(s1,rep)
   prob$c     <- rep(c1,rep)
   prob$d     <- rep(d1,rep)

   if (length(theta) != 1 && length(theta) != rep)
    stop("Length of person parameters vectors not well defined")
   if (length(S) != 1 && length(S) != rep)
    stop("Length of person parameters vectors not well defined")
   if (length(C) != 1 && length(C) != rep)
    stop("Length of person parameters vectors not well defined")
   if (length(D) != 1 && length(D) != rep)
    stop("Length of person parameters vectors not well defined")
   prob$THETA <- rep(theta,length(b))
   prob$S     <- rep(S,length(b))
   prob$C     <- rep(C,length(b))
   prob$D     <- rep(D,length(b))
   prob$ID    <- rep((1:rep),length(b))

   prob$p <- pm4pl(theta=prob$THETA,S=prob$S,C=prob$C,D=prob$D,s=prob$s,
                   b=prob$b,c=prob$c,d=prob$d)
   #prob   <- data.frame(round(prob,2))

   if (groups == TRUE) {
    res <- xyplot(p ~ b, groups=ID, data=prob, type=type,
                  main=main, ylab=ylab, xlab=xlab, lty=1:rep, pch=1:rep,
                  col.line="black", col.symbol="black")
    yCoord <- unlist(aggregate(prob$p, by=list(factor(prob$ID)), mean)$x)
    
    label <- switch(ID,
     ID    = paste(factor(1:rep), sep=""),
     ALL   = paste("(",factor(theta),", ",factor(S),", ",factor(C),", ",factor(D),")",sep=""),
     THETA = paste("Theta=",factor(theta), sep=""),
     NULL  = "")

    res <- update(res,
                  panel = function(...) {
                  panel.xyplot(...)
                  panel.text(x=theta, y=yCoord, cex=0.6,
                             labels=label)
                  })
    }
   if (groups == FALSE)
    res <- xyplot(p ~ b | factor(ID), data=prob, type=type,
                  main=main, ylab=ylab, xlab=xlab, lty=1:rep, pch=1:rep,
                  col.line="black", col.symbol="black")

   res$condlevels <- switch(ID,
    ALL = list(paste("(",levels(factor(theta)),", ",
                        levels(factor(S)),", ",
                        levels(factor(C)),", ",
                        levels(factor(D)),
                        ")",sep="")),
    ID    = list(levels(factor(1:rep))),
    THETA = list(levels(factor(theta))),
    NULL  = list(levels(factor(1:rep))))
   res$y.limits   <- c(-0.1,1.1)
   invisible(list(graphic=res, probability=prob))
   }

