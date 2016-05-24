make.dummy <- function(fact) {
   if (!is.factor(fact)) {
      stop("make.dummy needs a factor")
   }
   dum <- diag(nlevels(fact))[fact, ]
   colnames(dum) <- levels(fact)
   names(dum) <- rownames(fact)
   dum
}

dummy2factor <- function(x) {
   if (min(x) != 0 & max(x) != 1)
      stop("dummy2factor needs a matrix of dummy variables")
   ord <- apply(x, 1, order)[ncol(x), ]
   nms <- colnames(x)
   if (is.null(nms)) {
      nms <- paste("Var", as.character(1:ncol(x)), sep="")
   }
   as.factor(nms[ord])
}

sp.summ <- function(y, n.cut=c(5, 10, 20)) {   
   mx <- apply(y, 2, max, na.rm=TRUE)
   n.occur <- apply(y>0, 2, sum, na.rm=TRUE)
   n2 <- Hill.N2(y)
   mm <- matrix(NA, nrow=ncol(y), ncol=length(n.cut))
   for (i in 1:length(n.cut)) {
      mm[, i] <- apply(y>n.cut[i], 2, sum)
   }
   colnames(mm) <- paste("N", sprintf("%03d", n.cut))
   data.frame(N.occur=n.occur, N2 = n2, Max.abun = mx, mm)
}

site.summ <- function(y, max.cut=c(2, 5, 10, 20)) {   
   tot <- rowSums(y, na.rm=TRUE)
   n.taxa <- apply(y>0, 1, sum, na.rm=TRUE)
   mx <- apply(y, 1, max, na.rm=TRUE)
   n2 <- Hill.N2(y, margin=1)
   mm <- matrix(NA, nrow=nrow(y), ncol=length(max.cut))
   for (i in 1:length(max.cut)) {
      mm[, i] <- apply(y>max.cut[i], 1, sum)
   }
   colnames(mm) <- paste("M.", sprintf("%03d", max.cut), sep="")
   data.frame(N.taxa=n.taxa, N2 = n2, Max=mx, Total=tot, mm)
}

write.list.Excel <- function(x, fName) {
  if(.Platform$OS.type == "windows" & .Machine$sizeof.pointer < 5) {
     if (file.exists(fName)) 
         if (!file.remove(fName))
            stop("Could not remove existing file - is it open?")
     if (! ("list" %in% class(x)))
        stop("object should be a list")
     if (requireNamespace("RODBC", quietly=TRUE)==FALSE) {
        stop("This function requires package RODBC")
     }
     on.exit(RODBC::odbcCloseAll())
     #     fp <- RODBC:::full.path(fName)
#     con <- paste("Driver={Microsoft Excel Driver (*.xls)};DriverId=790;Dbq=", fp, ";DefaultDir=", dirname(fp), ";", sep = "")
#     con = paste(con, "ReadOnly=False", sep = ";")
#     channel <- odbcDriverConnect(con, tabQuote = c("[", "]"))
     channel <- RODBC::odbcConnectExcel(fName, readOnly=FALSE)
     if (channel == -1)
        stop(paste("Could not open file ", fName, "for writing. Is it already open?", sep=""))
     nms <- names(x)
     for (i in 1:length(x)) {
        if (class(x[[i]]) == "data.frame") {
           colnames(x[[i]]) <- gsub("\\.", "_", colnames(x[[i]]))
           RODBC::sqlSave(channel, x[[i]], nms[i], colnames=FALSE, rownames=TRUE)
        }
     }    
     RODBC::odbcCloseAll()
  } else {
    stop("This function is only available on 32 bit Windows. See package WriteXLS for a Linux / MacOS X alternative.")
  }
}

Hill.N2 <- function(df, margin=2) {
   if (margin == 2) {
     yk <- colSums(df)
     N2 <- 1/colSums(sweep(df, margin, yk, "/")^2, na.rm=TRUE)
   } else {
      if (margin == 1) {
         yk <- rowSums(df)
        N2 <- 1/rowSums(sweep(df, margin, yk, "/")^2, na.rm=TRUE)
      } else {
        stop("Margin out of bounds in Hill.N2")
      }
   }
   N2
}

dot <- function (x, ...) {
   UseMethod("dot")
}

dot.default <- function (x, ...) {
   stop(paste("No appropriate method for dot for object of class", class(x)[1]))
}

dot.data.frame <- function(x, head = 3, tail=1, dotrows=2, ...)
{
# From code posted to R help request on 8/2/07
   x <- as.data.frame(x)
   x <- format(rbind(head(x,head + dotrows), tail(x,tail)))
   if(dotrows>0)
   {
      x[(head + 1):(head + dotrows),] <- "."
      for(i in 1:dotrows){ rownames(x)[head+i]<-paste(".", substring("           ", 1, i-1))}
   }
   x
}


