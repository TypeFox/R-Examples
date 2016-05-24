################################################################################
# print.summary.ca(): Printing summarized ca objects (ca package 0.70)
################################################################################
print.summary.ca <- function(x, ...){
  object <- x
  r.out  <- object$rows
  c.out  <- object$columns
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
    if (length(rev)>1) {
      st <- round(nchars * rev/sum(rev), 0)
      } else {
      st <- nchars
      }
    scree <- character(length(Dim))
    for (q in Dim) {
      s1 <- paste(rep("*", st[q]), collapse = "")
      s2 <- paste(rep(" ", nchars - st[q]), collapse = "")
      scree[q] <- paste(" ", s1, s2, sep = "")
      }
    temp0  <- c(" "," ------"," "," "," ")
    temp1  <- c("Total:", sum(EV), "", "", "")
    Value0 <- round(Value, 6)
    Value1 <- round(sum(Value), 6)
    EV1    <- round(EV, 1)
    EV2    <- round(sum(EV), 1)
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
    EV.1 <- floor(log(EV1, base = 10))
    EV.1[EV.1 < 0] <- 0
    EV.2 <- as.character(EV.1)
    EV.2[EV.1 == 1] <- ""
    EV.2[EV.1 == 0] <- "0"
    EV1     <- remzero(gluezero(paste(EV.2, EV1, sep = ""), 4, 3))
    EV.sp   <- paste(rep(" ", ifelse(max(EV.1==2), 0, 1)), collapse = "", sep = "")
    EV1     <- paste(EV.sp, EV1, sep = "")
    CUMEV1  <- round(CUMEV, 1)
    CUMEV.1 <- floor(log(CUMEV1, base = 10))
    CUMEV.1[CUMEV.1 < 0] <- 0
    CUMEV.2 <- as.character(CUMEV.1)
    CUMEV.2[CUMEV.1 == 2] <- ""
    CUMEV.2[CUMEV.1 == 1] <- "0"
    CUMEV.2[CUMEV.1 == 0] <- "00"
    CUMEV1 <- remzero(gluezero(paste(CUMEV.2, CUMEV1, sep = ""), 5, 4), 
                      doub = TRUE)
    scree.out <- data.frame(Dim   = c(Dim, "", "Total:"), 
                            Value = c(gluezero(Value0), "--------", gluezero(Value1)), 
                            EV    = c(EV1, "-----", gluezero(EV2, 5, 4)), 
                            CUMEV = c(CUMEV1, "", ""), 
                            scree = c(scree, "", ""))
    colnames(scree.out) <- c("dim", "value", "  %", "cum%", " scree plot")
    cat("Principal inertias (eigenvalues):\n\n")
    scree.out <- as.matrix(scree.out)
    rownames(scree.out) <- rep("", nrow(scree.out))
    print(scree.out, quote = FALSE)
    cat("\n")
    }
 # print row/column summary:
  if (!is.null(object$rows)){
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
    cat("\nRows:\n")
    print(as.matrix(r.new), quote = FALSE, right = TRUE)
    }
  if (!is.null(object$columns)){
    n1 <- dim(c.out)[1]
    n2 <- dim(c.out)[2]
    c.names <- dimnames(c.out)[[2]]
    c.dummy <- rep("|", n1)
    c.new   <- cbind(c.dummy, c.out[,1], c.dummy, c.out[,2:4])
    c.nn    <- c("", c.names[1], "", c.names[2:4])
    for (q in 1:((n2 - 4) / 3)){
      c.new <- cbind(c.new, c.dummy, c.out[,(5 + (q - 1) * 3):(5 + q * 3 - 1)])
      c.nn  <- c(c.nn, "", c.names[(5 + (q - 1) * 3):(5 + q * 3 - 1)])
      }
    c.new <- cbind(c.new, c.dummy)
    c.nn  <- c(c.nn, "")
    colnames(c.new) <- c.nn
    rownames(c.new) <- 1:n1
    cat("\nColumns:\n")
    print(as.matrix(c.new), quote = FALSE, right = TRUE)
    }
  }
################################################################################
