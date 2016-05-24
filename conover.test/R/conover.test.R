# version 1.1.1 January 6, 2016 by alexis.dinno@pdx.edu
# perform Conover-Iman test of multiple comparisons using rank sums

p.adjustment.methods <- c("none","bonferroni","sidak","holm","hs","hochberg","bh","by")

conover.test <- function(x=NA, g=NA, method=p.adjustment.methods, kw=TRUE, label=TRUE, wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05) {

  # FUNCTIONS

  # kwallis.test: a custom Kruskal Wallis test function to support conover.test
  # Note: does not currently play nicely with missing values.
  kwallis.test <- function(x=NA,g=NA,pass=0) {
    # FUNCTIONS
    # tiedranks: enumerates tied values from a vector of ranks
    # ranks: a vector of rank values
    # returns: ties, a vector of values that are tied
    tiedranks <- function(ranks) {
      ranks <- sort(ranks)
      #enumerate tied values
      ties <- c()
      for (i in 2:length(ranks)) {
        if (ranks[i-1] == ranks[i]) {
          if (length(ties) > 0) {
            if (ranks[i-1] != tail(ties,n=1)) {
              ties <- c(ties,ranks[i-1])
              }
            }
           else {
            ties <- c(ranks[i-1])
            }
          }
        }
      return(ties)
      }

    if (pass == 1) {
#      x <- x[[1]]
      }
  
    #set up data by lists
    if (is.list(x)) {
      N <- 0
      for (i in 1:length(x)) {
        N <- N + length(x[[i]])
        }
      Data <- matrix(NA,N,4)
      for (i in 1:N) {
        Data[i,1] <- i
        }
      obs <- c()
      group <- c()
      for (i in 1:length(x)) {
        obs <- c(obs,x[[i]])
        group <- c(group,rep(i,length(x[[i]])))
        }
      Data[,2] <- obs
      if (length(g) > 1) {
#        Data[,3] <- as.integer(factor(g))
        Data[,3] <- g
        }
       else {
        Data[,3] <- group        
        }
      Data[,4] <- rank(Data[,2], ties.method="average", na.last=NA)
      }
  
    #set up data by groups
    if (!is.list(x)) {
#      g <- as.integer(factor(g))
      N <- length(x)
      Data <- matrix(NA,length(x),4)
      Data[,1] <- 1:length(x)
      Data[,2] <- x
      Data[,3] <- g
      Data[,4] <- rank(Data[,2], ties.method="average", na.last=NA)
      }
    k <- length(unique(Data[,3]))
  
    #calculate ties adjustment
    ranks <- Data[,4]
    ranks <- ranks[!is.na(ranks)]
    ties <- tiedranks(ranks)
    r <- length(ties)
    tiesadj <- 0
    if (r > 0) {
      for (s in 1:r) {
        tau <- sum(ranks==ties[s])
        tiesadj <- tiesadj + (tau^{3} - tau)
        }
      }
    tiesadj <- 1-(tiesadj/((N^3) - N))
  
    #calculate H
    ranksum <- 0
    for (i in 1:k) {
      ranksum <- ranksum + ((sum(Data[,4][Data[,3]==i]))^2)/(sum(Data[,3]==i))
      }
    H  <- ((12/(N*(N+1)))*ranksum - 3*(N+1))/tiesadj
    df <- k-1
    p  <- pchisq(H,k-1,lower.tail=FALSE)
  
    #present output
    output <- paste("Kruskal-Wallis chi-squared = ",round(H,digits=4),", df = ",df,", p-value = ",round(p,digits=2),"\n\n" ,sep="")
  
    invisible(list(output=output,H=H,df=df,p=p,N=N,Data=Data,data.name=data))
    }

  # all.integers robustly tests whether all elements of a vector are integers
  all.integers <- function(x) {
    for (i in length(x)) {
      if (is.na(x[i]) | is.list(x[i]) | length(x[i]) > 1 | !is.numeric(x[i])) {
        return(FALSE)
        }
      if (x[i]%%1!=0) {
        return(FALSE)
        }
      }
    return(TRUE)
    }

  # tpad: returns the pad string, n times, defaulting to spaces
  # n: a number of replications; pad: string to replicate
  tpad <- function(n=1, pad=" ") {
    if (n == 0) {
      return("")
      }
    return( paste0(rep(x=pad,times=n),collapse="") )
    }

  # tformat: formats t values for display in table
  # t: a real t-value
  # returns: a formatted string
  tformat <- function(t) {
    if (t < 0) {
      sign <- "-"
      }
     else {
       sign <- " "
       }
     leftspaces <- max(2,floor(log10(abs(t)))+2)
     leftdigits <- floor(abs(t))
     rightspaces <- 8 - leftspaces
     rightdigits <- substr(paste0(abs(t) - floor(abs(t)),"00000000"),3,rightspaces+2)
     return(paste0(sign,leftdigits,".",rightdigits))
    }

  # pformat: formats p values for display in table
  # p: a real p-value
  # returns: a formatted string
  pformat <- function(p) {
     if (p == 1) {
       return("1.0000")
       }
      else {
       rightspaces <- 4
       rightdigits <- substr(paste0(sprintf("%1.4f",p),"000000000"),3,rightspaces+2)
       return(paste0("0.",rightdigits))
       }
    }

  # centertext: centers a string within a specific width
  centertext <- function(text,width=80,lmargin=2,rmargin=2) {
    textwidth <- nchar(text)
    if (textwidth <= width-lmargin-rmargin) {
      text <- substr(text,1,width-lmargin-rmargin)
      }
    buff <- (width-lmargin-rmargin-textwidth)
    if (buff%%2 == 0) {
      return(paste(paste(rep(" ",buff/2),collapse=""),text,paste(rep(" ",buff/2),collapse=""),"\n",collapse=""))
      }
     else {
      return(paste(paste(rep(" ",1+(buff-1)/2),collapse=""),text,paste(rep(" ",(buff-1)/2),collapse=""),"\n",collapse=""))
      }
    }

  # conovertestheader: displays the Conover-Iman test table headers.
  conovertestheader <- function(groupvar,colstart,colstop,rmc) {
    if (rmc==FALSE) {
      cat("Col Mean-|\nRow Mean |")
      }
     else {
      cat("Row Mean-|\nCol Mean |")
      }
    groupvalues  <- levels(factor(groupvar))
    for (col in colstart:colstop) {
      vallab  <- substr(groupvalues[col],1,8)
      pad     <- 8-nchar(vallab)
      colhead <- paste0(tpad(n=pad),vallab)
      cat(paste0("   ",substr(colhead,1,8)),sep="")
      }
    cat("\n")  
    separatorlength = 10 + 11*(colstop-colstart)+1
    cat(tpad(n=9,pad="-"),"+",tpad(n=separatorlength,pad="-"),"\n",sep="")
    }

  # conovertestttable displays the Conover-Iman test t values
  conovertestttable <- function(groupvar, index, Tvalues, colstart, colstop) {
    groupvalues  <- levels(factor(groupvar))

    # Row headers
    vallab  <- substr(groupvalues[index],1,8)
    pad     <- 8-nchar(vallab)
    rowhead <- paste0(tpad(n=pad),vallab)
    cat(paste0(rowhead," |"))

    # Table t entries
    for (i in colstart:colstop) {
      t <- Tvalues[(((index-2)*(index-1))/2) + i]
      cat(paste0("  ",tformat(t)))
      }
    cat("\n")
    }

  # conovertestptable: displays the Conover-Iman test p values
  conovertestptable <- function(index, P, colstart, colstop, Reject,last) {
    # Blank row header
    cat("         |")

    # Table p entries
    for (i in colstart:colstop) {
      p <- P[(((index-2)*(index-1))/2) + i]
      if ( Reject[(((index-2)*(index-1))/2) + i] == 0) {
        cat(paste0("     ",pformat(p)),sep="")
        }
       else {
        cat(paste0("    ",pformat(p)),"*",sep="")
        }
      }
    cat("\n")

    # Close out with another blank row header
    if (last == 0) {
      cat("         |\n") 
      }
    }

  # VALIDATIONS & PREPARATIONS
  # names for output
  xname <- paste(if(is.name(substitute(x))) {deparse(substitute(x))} else {"x"})
  gname <- paste(if(is.name(substitute(g))) {deparse(substitute(g))} else {"group"})
  xaslist <- is.list(x)
  tempg <- 1:length(x)
  
  # casewise deletion of missing data
  if (length(g) > 1) {
    if (label==TRUE) {
      if (is.factor(g)) {
        glevels <- levels(g)
        }
       else {
        glevels <- unique(g)
        }
      }
     else {
      glevels <- names(table(addNA(g, ifany = TRUE)))
      }
    Data <- data.frame(cbind(x,g))
    Data <- Data[!is.na(unlist(Data$x)),]
    Data <- Data[!is.na(Data$g),]
    x <- as.numeric(factor(Data$x))
    g <- factor(Data$g)
    levels(g) <- glevels
    }
   else {
     g <- c()
     for (i in 1:length(x)) {
       for (j in 1:length(x[[i]])) {
         g <- c(g,i)
         }
       }
    x <- as.numeric(unlist(x)[!is.na(unlist(x))])
    }
  # validate method
  if (length(method) > 1) {
    method <- "none"
    }
    
  if (tolower(method) != "none" & tolower(method) != "bonferroni" & tolower(method) != "sidak" & tolower(method) != "holm" & tolower(method) != "hs" & tolower(method) != "hochberg" & tolower(method) != "bh" & tolower(method) != "by") {
    stop("method must be one of: none, bonferroni, sidak, hs, bh, or by")
    }

  if (tolower(method)=="none") {
    Name <- "(No adjustment)"
    }
  if (tolower(method)=="bonferroni") {
    Name <- "(Bonferroni)"
    }
  if (tolower(method)=="sidak") {
    Name <- "(\u0160id\u00E1k)"
    }
  if (tolower(method)=="holm") {
    Name <- "(Holm)"
    }
  if (tolower(method)=="hs") {
    Name <- "(Holm-\u0160id\u00E1k)"
    }
  if (tolower(method)=="hochberg") {
    Name <- "(Hochberg)"
    }
  if (tolower(method)=="bh") {
    Name <- "(Benjamini-Hochberg)"
    }
  if (tolower(method)=="by") {
    Name <- "(Benjamini-Yekuteili)"
    }

  # validate that x is longer than 1
  if (length(unlist(x))==1) {
    stop("too few observations in x.")
    }
  # validate that x is a numeric vector of data values, or a list of numeric 
  # data vectors, and is not NA
  if (TRUE %in% is.na(unlist(x)) | !is.vector(unlist(x), mode="numeric") ) {
    stop("x must contain a numeric vector of data values, or a list of numeric data vectors.")
    }
  # validate that g is not missing if x is a vector
  if (!xaslist & TRUE %in% is.na(g) ) {
    stop("when specifying x as a vector, you must include g.")
    }
  # validate that g is not NA
  if (length(g) > 1 & TRUE %in% is.na(g)) {
    stop("g must have no missing values.")
    }
  # validate that g is factor or vector.
  if (length(g) > 1 & ( !is.factor(g) & !is.vector(g) ) ) {
    
    stop("g must be a factor, character vector, or integer vector.")
    }
  # validate that g is a vector of mode = character or mode = integer.
  if (length(g) > 1 & is.vector(g) ) {
    if ( !is.vector(g, mode="character") & !all.integers(g) ) {
      stop("g must be a factor, character vector, or integer vector.")
      }
    }

  # CALCULATIONS
  out <- NULL
  if (xaslist & length(g)==1) {
    kwallis.test(x, 1:length(x), pass=1) -> out
    }
  if (length(g)>1 & xaslist) {
    kwallis.test(x, g, pass=1) -> out
    }
  if (length(g)>1 & !xaslist) {
    kwallis.test(x, g, pass=0) -> out
    }

  if (kw==TRUE) {
    cat("  Kruskal-Wallis rank sum test\n\ndata: ",xname," and ",gname,"\n",sep="")
    cat(out$output)
    }
    chi2    <- out$H
    k       <- out$df + 1
    N       <- out$N
    Data    <- out$Data
    m       <- k*(k-1)/2
    Tvalues <- rep(0,m)
    P       <- rep(0,m)
    nu      <- N - k

    # calculate pooled variance to be used in generating t statistics
    ranks   <- Data[,4]
    varA  <- 1/(N-1)
    varB  <- 0
    for (i in 1:N) {
      varB <- varB + ranks[i]^2
      } 
    varB <- varB - (N*(((N+1)^2)/4))
    S2 <- varA*varB

    # Calculate t test statistics
    for (i in 2:k) {
      for (j in 1:(i-1)) {
        ni <- sum(Data[,3]==i)
        nj <- sum(Data[,3]==j)
        meanranki <- mean(Data[,4][Data[,3]==i])
        meanrankj <- mean(Data[,4][Data[,3]==j])
        if (rmc == FALSE) {
          t <- (meanrankj-meanranki)/(sqrt(S2*((N-1-chi2)/(N-k)))*sqrt(1/nj + 1/ni))
          }
        if (rmc == TRUE) {
          t <- (meanranki-meanrankj)/(sqrt(S2*((N-1-chi2)/(N-k)))*sqrt(1/nj + 1/ni))
          }
        index <- ((i-2)*(i-1)/2) + j
        Tvalues[index] <- t
        }
      }

  # Calculate p-values for t statistics, and adjust as needed
  P <- pt(q=abs(Tvalues), df=nu, lower.tail=FALSE)
  
  #calculate adjusted p-values based on method argument
  Reject <- rep(0,m)
  # No adjustment for multiple comparisons
  if (tolower(c(method))=="none") {
    P.adjust <- P
    }

  # Control FWER using (Dunn's) Bonferroni
  if (tolower(c(method))=="bonferroni") {
    P.adjust <- pmin(1,P*m)
    }

  # Control FWER using Šidák
  if (tolower(c(method))=="sidak") {
    P.adjust <- pmin(1,1 - (1-P)^m)
    }

  # Control FWER using Holm(-Bonferroni)
  if (tolower(c(method))=="holm") {
    Psort <- matrix(c(P,1:m,rep(0,m)),3,m,byrow=TRUE)
    Psort <- Psort[,order(Psort[1,])]
    for (i in 1:m) {
      adjust <- m+1-i
      Psort[1,i] <- pmin(1,Psort[1,i]*adjust)
      if (i==1) {
        Psort[3,i] <- Psort[1,i] <= alpha/2
        }
       else {
         Psort[3,i] <- ((Psort[1,i] <= alpha/2) & Psort[3,i-1] != 0)
         }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }

  # Control FWER using Holm-Šidák
  if (tolower(c(method))=="hs") {
    Psort <- matrix(c(P,1:m,rep(0,m)),3,m,byrow=TRUE)
    Psort <- Psort[,order(Psort[1,])]
    for (i in 1:m) {
      adjust <- m+1-i
      Psort[1,i] <- pmin(1,(1 - ((1 - Psort[1,i])^adjust)))
      if (i==1) {
        Psort[3,i] <- Psort[1,i] <= alpha/2
        }
       else {
         Psort[3,i] <- ((Psort[1,i] <= alpha/2) & Psort[3,i-1] != 0)
         }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }

  # Control FWER using Hochberg
  if (tolower(c(method))=="hochberg") {
    Psort <- matrix(c(P,1:m,rep(0,m)),3,m,byrow=TRUE)
    Psort <- Psort[,order(Psort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      adjust <- i
      Psort[1,i] <- min(1,Psort[1,i]*adjust)
      if (i==1) {
        Psort[3,i] <- Psort[1,i] <= alpha/2
        }
       else {
         Psort[3,i] <- ((Psort[1,i] <= alpha/2) | Psort[3,i-1] == 1)
         }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }

  # Control FDR using Benjamini-Hochberg
  if (tolower(c(method))=="bh") {
    Psort <- matrix(c(P,1:m,rep(0,m)),3,m,byrow=TRUE)
    Psort <- Psort[,order(Psort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      adjust <- (m/(m+1-i))
      Psort[1,i] <- min(1,Psort[1,i]*adjust)
      if (i==1) {
        Psort[3,i] <- Psort[1,i] <= alpha/2
        }
       else {
         Psort[3,i] <- ((Psort[1,i] <= alpha/2) | Psort[3,i-1] == 1)
         }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }

  # Control FDR using Benjamini-Yekuteili
  if (tolower(c(method))=="by") {
    Psort <- matrix(c(P,1:m,rep(0,m)),3,m,byrow=TRUE)
    Psort <- Psort[,order(Psort[1,], decreasing=TRUE)]
    for (i in 1:m) {
      adjust <- (m/(m+1-i))*sum(1/(1:m))
      Psort[1,i] <- min(1,Psort[1,i]*adjust)
      if (i==1) {
        Psort[3,i] <- Psort[1,i] <= (alpha/2) # reverse sorted, so m-i+1, rather than i
        }
       else {
         Psort[3,i] <- ((Psort[1,i] <= alpha/2) | Psort[3,i-1] == 1)
         }
      }
    Psort <- Psort[,order(Psort[2,])]
    P.adjust <- Psort[1,]
    Reject <- Psort[3,]
    }

  # OUTPUT
  cat("\n")
  if (table==TRUE | list==TRUE) {
    if ((TRUE %in% is.na(g))) {
      title <- paste("Comparison of ",xname," across ",k," groups",sep="")
      }
     else {
      title <- paste("Comparison of ",xname," by ",gname,sep="")
      }
    cat(centertext(title))
    cat(centertext(Name))
    }

  if (table==TRUE) {
    # Need to determine how many tables (reps) to output
    reps      <- floor((k-1)/6)
    laststart <- k - (reps*6) + 1
    kminusone <- k - 1
    if (label==FALSE) {
      g <- as.numeric(g)
      }
    if (length(g)==1) {
      g <- 1:k
      }
  
    # Replication loop for >7 groups, no wrap
    if (wrap==FALSE) {
      if (k > 7) {
        for (rep in 1:reps) {
          colstart <- (6*rep)-5
          colstop  <- 6*rep
          conovertestheader(g,colstart,colstop,rmc)
          # Table body
          for (i in (colstart+1):k) {
            colstop <- min(i-1,6*rep)
            conovertestttable(g,i,Tvalues,colstart,colstop)
            if (i < k) {
              conovertestptable(i,P.adjust,colstart,colstop,Reject,0)
              }
             else {
              conovertestptable(i,P.adjust,colstart,colstop,Reject,1)
              }
            }
          }
        # End of table
        if (laststart < k) {
          conovertestheader(g,laststart,kminusone,rmc)
          # Table body
          for (i in (laststart+1):k) {
            conovertestttable(g,i,Tvalues,laststart,kminusone) 
            if (i < k) {
              conovertestptable(i,P.adjust,laststart,kminusone,Reject,0)
              }
             else {
              conovertestptable(i,P.adjust,laststart,kminusone,Reject,1)
              }
            }
          }
        }
  
      # Replication loop for <=7 groups
      if (k <= 7) {
        conovertestheader(g,1,kminusone,rmc)
        # Table body
        for (i in 2:k) {
          colstop <- i-1
          conovertestttable(g,i,Tvalues,1,colstop) 
          if (i < k) {
            conovertestptable(i,P.adjust,1,colstop,Reject,0)
            }
           else {
            conovertestptable(i,P.adjust,1,colstop,Reject,1)
            }
          }
        }
      }
  
    # Replication loop for >7 groups, with wrap
    if (wrap==TRUE) {
      conovertestheader(g,1,kminusone,rmc)
      # Table body
      for (i in 2:k) {
        colstop <- i-1
        conovertestttable(g,i,Tvalues,1,colstop) 
        if (i < k) {
          conovertestptable(i,P.adjust,1,colstop,Reject,0)
          }
         else {
          conovertestptable(i,P.adjust,1,colstop,Reject,1)
          }
        }
      }
      
    cat("\n")
    }

  # Output pairwise comparisons as list if requested.
  if (list==TRUE) {
    groupvalues  <- levels(factor(g))
    stringlength <- 2*max(nchar(groupvalues)) + 4
    
    # Output list header
    cat("\nList of pairwise comparisons: t (p-value)")
    cat("\n---------------------------------------------------\n")
    index <- 0
    for (i in 2:k) {
      for (j in 1:(i-1)) {
        index <- index + 1
          buffer <- stringlength - (nchar(groupvalues[i]) + nchar(groupvalues[j]) + 4) - 2
        if (rmc==FALSE) {
          cat(groupvalues[j]," - ",groupvalues[i],paste(rep(" ",buffer))," :\t",tformat(Tvalues[index]), " (",pformat(P.adjust[index]),")\n", sep="")
          }
        if (rmc==TRUE) {
          cat(groupvalues[i]," - ",groupvalues[j],paste(rep(" ",buffer))," :\t",tformat(Tvalues[index]), " (",pformat(P.adjust[index]),")\n", sep="")
          }
        }
      }
    cat("\n")
    }

  # Create comparisons variable for returned values (whether the list option
  # is TRUE or FALSE
  comparisons <- rep(NA,(k*(k-1)/2))
  groupvalues  <- levels(factor(g))
  index <- 0
  for (i in 2:k) {
    for (j in 1:(i-1)) {
      index <- index + 1
      if (rmc==FALSE) {
        comparisons[index] <- paste0(groupvalues[j]," - ",groupvalues[i])
        }
      if (rmc==TRUE) {
        comparisons[index] <- paste0(groupvalues[i]," - ",groupvalues[j])
        }
      }
    }
   
  invisible(list(chi2=chi2,T=Tvalues,P=P,P.adjusted=P.adjust,comparisons=comparisons))
  }