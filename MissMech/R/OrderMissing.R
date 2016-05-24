OrderMissing <- function(data, del.lesscases = 0)
{
# case order has the order of the original data
  y <- data
  if (is.data.frame(y)) {
  y <- as.matrix(y)
  } 
  if(!is.matrix(y))
  {
    cat("Warning data is not a matrix or data frame")
    stop("")
  }
  if(length(y)==0)
 {
   cat("Warning: data is empty")
   return
 }
  names <- colnames(y)
  n <- nrow(y)
  pp <- ncol(y)
  yfinal <- c()
  patused <- c()
  patcnt <- c()
  caseorder <- c()
  removedcases <- c()
  ordertemp <- c(1:n)
  ntemp <- n
  ptemp <- pp
  done <- FALSE
  yatone <- FALSE
  while(!done)
  {
    pattemp <- is.na(y[1, ])
    indin <- c()
    indout <- c()
    done <- TRUE
    for(i in 1:ntemp)
    {
      if(all(is.na(y[i, ]) == pattemp))
      {
        indout <- c(indout, i)
      }  else {
        indin <- c(indin, i)
        done <- FALSE
      }
    }
    if(length(indin) == 1) yatone = TRUE
    yfinal <- rbind(yfinal, y[indout, ])
    y <- y[indin, ]
    caseorder <- c(caseorder, ordertemp[indout])
    ordertemp <- ordertemp[indin]
    patcnt <- c(patcnt, length(indout))
    patused <- rbind(patused, pattemp)
    if(yatone)
    {
      pattemp <- is.na(y)
      yfinal <- rbind(yfinal, matrix(y,ncol=pp))
      y <- c()
      indin <- c()
      indout <- c(1)
      caseorder <- c(caseorder, ordertemp[indout])
      ordertemp <- ordertemp[indin]
      patcnt <- c(patcnt, length(indout))
      patused <- rbind(patused, pattemp)
      done <- TRUE
    }
    
    if(!done) ntemp <- nrow(y)
  }
  #yfinal <- rbind(yfinal, y)
  caseorder <- c(caseorder, ordertemp)
  patused <- ifelse(patused, NA, 1)
  rownames(patused) <- NULL
  colnames(patused) <- names
  spatcnt <- cumsum(patcnt)
  dataorder <- list(data = yfinal, patused = patused, patcnt = patcnt, 
                   spatcnt = spatcnt, g = length(patcnt), 
                   caseorder = caseorder, removedcases = removedcases)
  dataorder$call <- match.call()
  class(dataorder) <- "orderpattern"
  if(del.lesscases > 0)
  {
    dataorder <- DelLessData(dataorder, del.lesscases)
  }   
  dataorder$patused <- matrix(dataorder$patused, ncol = pp)
  colnames(dataorder$patused) <- names
  dataorder
}
#Order <- function(x, ...) UseMethod("Order")
#Order.default <- function(x, ...) {
# temp <- OrderMissing(x)
# temp$call <- match.call()
# class(temp) <- "ordered"
# temp
#}
print.orderpattern <- function(x, ...) {
 cat("Call:\n")
 print(x$call)
 cat("\nNumber of Ptterns: ", x$g, "\n")
 cat("\nPttern used:\n")
 ni <- x$patcnt
 disp.patt <- cbind(x$patused, ni)
 colnames(disp.patt)[ncol(disp.patt)] <- "Number of cases"
 rownames(disp.patt) <- rownames(disp.patt, do.NULL = FALSE, prefix = "group.")
 print(disp.patt, print.gap = 3) 
}

summary.orderpattern <- function(object, ...) {
  summary(object$data)
}
