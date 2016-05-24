################################################################################
# print.summary.mjca(): Printing summarized mjca objects (ca package 0.70)
################################################################################
print.summary.mjca <- function(x, ...){
  object <- x
  if (!is.null(object$scree)){
    cat("\n")
   # init:
    nchars <- 25
    Dim    <- object$scree[,1]
    ev     <- object$scree[,2]
    rev    <- object$scree[,3]
    crev   <- object$scree[,4]
    Value  <- ev[Dim]
    EV     <- rev[Dim]
    CUMEV  <- crev[Dim]
    sev <- object$sev
    if (length(rev)>1) {
      st <- round(nchars * rev/sum(rev), 0)
      } else {
      st <- nchars
      }
    scree <- character(length(Dim))
    for (q in Dim) {
      if (!is.na(st[q])){
        s1 <- paste(rep("*", st[q]), collapse = "")
        s2 <- paste(rep(" ", nchars - st[q]), collapse = "")
        scree[q] <- paste(" ", s1, s2, sep = "")
        } else {
        scree[q] <- " "
        }
      }
    temp0 <- c(" "," ------"," "," "," ")
    temp1 <- c("Total:", sev, "", "", "")
    valuesum <- round(object$tin, 6)
    gluezero <- function(item, dig = 8, point = NA){
      item0 <- paste(item, paste(rep(0, dig), collapse=""), sep = "")
      item1 <- strsplit(item0,"", fixed = TRUE)
      pastebit <- function(x, digits = dig, poin = point){
        if(!is.na(poin)) { 
          x[poin] <- "."
          }
        paste(x[1:dig], collapse = "")
        }
      unlist(lapply(item1, pastebit))
      }
    remzero <- function(x, doub = FALSE){
      x0 <- strsplit(x, "", fixed = TRUE)
      pastebit2 <- function(x, doubl = doub){
        if (doubl){
           if (x[1]==0 & x[2]==0){
             x[1] <- " "
             x[2] <- " "
             }
          }
        if (x[1]==0) x[1] <- " "
        paste(x, collapse = "")
        }
      unlist(lapply(x0, pastebit2))
      }
    if (!is.na(EV[1])){
      EV.1 <- floor(log(EV, base = 10))
      EV.1[EV.1 < 0] <- 0
      EV.2 <- as.character(EV.1)
      EV.2[EV.1 == 1] <- ""
      EV.2[EV.1 == 0] <- "0"
      EV <- remzero(gluezero(paste(EV.2, EV, sep = ""), 4, 3))
      EV.sp <- paste(rep(" ", ifelse(max(EV.1==2), 0, 1)), collapse = "", sep = "")
      EV    <- paste(EV.sp, EV, sep = "")
      }
   # add leading space:
    if (!is.na(CUMEV[1])){
      CUMEV.1 <- floor(log(CUMEV, base = 10))
      CUMEV.1[CUMEV.1 < 0] <- 0
      CUMEV.2 <- as.character(CUMEV.1)
      CUMEV.2[CUMEV.1 == 2] <- ""
      CUMEV.2[CUMEV.1 == 1] <- "0"
      CUMEV.2[CUMEV.1 == 0] <- "00"
      CUMEV <- remzero(gluezero(paste(CUMEV.2, CUMEV, sep = ""), 5, 4), doub = TRUE)
      }
    scree.out <- data.frame(Dim   = c(Dim, "", "Total:"), 
                            Value = c(gluezero(as.character(Value)), "--------", gluezero(as.character(valuesum), 8, 2)), 
                            EV    = c(EV, "-----", ifelse(!is.na(sev), gluezero(sev,5,4), "")), 
                            CUMEV = c(CUMEV, "", ""), 
                            scree = c(scree, "", ""))
    colnames(scree.out) <- c("dim", "value", "  %", "cum%", " scree plot")
    if (is.na(object$JCA.nit[1])){
      cat("Principal inertias (eigenvalues):\n\n")
      scree.out <- as.matrix(scree.out)
      rownames(scree.out) <- rep("", nrow(scree.out))
      print(as.matrix(scree.out), quote = FALSE)
      cat("\n")
      } else {
      cat("Principal inertias (eigenvalues):\n\n")
      scree.out <- as.matrix(scree.out[,1:2])
#      dimnames(scree.out)[[1]] <- rep("", length(dimnames(scree.out)[[1]]))
      rownames(scree.out) <- rep("", nrow(scree.out))
      print(as.matrix(scree.out), quote = FALSE)
      cat(paste("\n Diagonal inertia discounted from eigenvalues: ", round(object$JCA.ind, 7), sep = ""))
      cat(paste("\n Percentage explained by JCA in ", object$JCA.nd, " dimensions: ", sev, "%", sep = ""))
      cat("\n (Eigenvalues are not nested)")
      cat(paste("\n [Iterations in JCA: ", object$JCA.nit, " , epsilon = ", round(object$JCA.eps, 7), "]\n\n", sep = ""))
      }
    }

 # print row/column summary:
  if (!is.null(object$rows)){
    r.out   <- object$rows
    n1      <- dim(r.out)[1]
    n2      <- dim(r.out)[2]
    r.names <- dimnames(r.out)[[2]]
    r.dummy <- rep("|", n1)
    r.new   <- cbind(r.dummy, r.out[,1], r.dummy, r.out[,2:4])
    r.nn    <- c("", r.names[1], "", r.names[2:4])
    for (q in 1:((n2 - 4) / 3)){
      r.new <- cbind(r.new, r.dummy, r.out[,(5 + (q - 1) * 3):(5 + q * 3 - 1)])
      r.nn  <- c(r.nn, "", r.names[(5 + (q - 1) * 3):(5 + q * 3 - 1)])
      }
    r.new <- cbind(r.new, r.dummy)
    r.nn  <- c(r.nn, "")
    colnames(r.new) <- r.nn
    rownames(r.new) <- 1:n1
    # print rows
    cat("\nRows:\n")
    print(as.matrix(r.new), quote = FALSE, right = TRUE)
    }
 ### COLUMNS:
  if (!is.null(object$columns)){
    c.out   <- object$columns
    n1      <- dim(c.out)[1]
    n2      <- dim(c.out)[2]
    c.names <- dimnames(c.out)[[2]]
    c.dummy <- rep("|", n1)
    if (is.na(object$JCA.nit[1])){
      c.new   <- cbind(c.dummy, c.out[,1], c.dummy, c.out[,2:4])
      c.nn    <- c("", c.names[1], "", c.names[2:4])
      for (q in 1:((n2 - 4) / 3)){
        c.new <- cbind(c.new, c.dummy, c.out[,(5 + (q - 1) * 3):(5 + q * 3 - 1)])
        c.nn <- c(c.nn, "", c.names[(5 + (q - 1) * 3):(5 + q * 3 - 1)])
        }
      } else { #JCA BELOW:
      c.new <- cbind(c.dummy, c.out[,1], c.dummy, c.out[,2:3], c.dummy, 
                     c.out[,4:(4+object$JCA.nd-1)], c.dummy, c.out[,(n2-1):n2])
      c.nn  <- c("", c.names[1], "", c.names[2:3], "", 
                 c.names[4:(4+object$JCA.nd-1)], "", c.names[(n2-1):n2])
      }
    c.new <- cbind(c.new, c.dummy)
    c.nn  <- c(c.nn, "")
    colnames(c.new) <- c.nn
    rownames(c.new) <- 1:n1
   ### PRINT COLUMNS:
    cat("\nColumns:\n")
    print(as.matrix(c.new), quote = FALSE, right = TRUE)
    cat("\n")
    }
  }
################################################################################
