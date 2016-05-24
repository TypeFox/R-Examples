# MedOr package for R (http://www.R-project.org)
# Copyright (C) 2012 Adriano Polpo, Carlos A. de B. Pereira.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


################################################################################
# Evaluate the confidence interval for
# the population median.
#
# Obs. This code works fine for "small" sample size
#      (n <= 1000), for large samples (n > 1000) will work too,
#      however it may be needed more computer memory and time
#      to evaluate the confidence interval.
conf.interval <- function(x,alpha=0.95,verbose=TRUE) {
  # X is a data vector.
  # alpha is the confidence level for the interval.
  # verbose: print the results.
  Call <- match.call()

  # Conditions to x
  if (!is.numeric(x))
    stop("x is not numeric")
  if (!is.vector(x))
    stop("x is not a vector")
  x <- sort(x)
  m <- length(x)
  if (m <= 2)
    stop("Sample size is too small (n <= 2)")
  # Conditions to alpha
  if ((!is.numeric(alpha)) && (length(alpha) == 1))
    stop("alpha is not numeric or have length bigger than 1")
  if ((alpha <= 0) || (alpha >= 1))
    stop("alpha must be a value in (0,1)")
  q <- 0.5
  # Conditions to q
  if ((!is.numeric(q)) && (length(q) == 1))
    stop("q is not numeric or have length bigger than 1")
  if ((q <= 0) || (q >= 1))
    stop("q must be a value in (0,1)")
  # Conditions to verbose
  if ((verbose != TRUE) && (verbose != FALSE))
    stop("verbose must be TRUE or FALSE")

  if (verbose) {
    if (m > 10^3)
      cat("Evaluating confidence interval. It may take some minutes...\n")
  }
  initial.time <- proc.time()[3]

  aux <- matrix(NA,choose(m,2),3)
  z <- 1-pbinom(seq(0,m-1,1),m,q)

  k <- 0
  i <- 1
  it <- 1
  repeat {
    j <- m
    repeat {
      k <- k+1
      aux[k,] <- c(i,j,z[i]-z[j])
      if (aux[k,3] < alpha)
        break
      if (j <= i+1)
        break
      j <- j-1

      if (verbose) {
        run.time <- (proc.time()[3]-initial.time)/60
        if (run.time > it) {
          it <- it+1
          cat("Time elapsed: ", sprintf("%.0f",floor(run.time))," min, so far...\n",sep="")
        }
      }
      next
    }
    if ((aux[k,3] < alpha) && (aux[k,2] == m))
      break
    if (i >= m-1)
      break
    i <- i+1
    next
  }
  aux <- aux[!is.na(aux[,1]),]

  if ((is.matrix(aux)) && (all(!is.na(aux))))  {
    aux <- aux[order(aux[,3]),]
    i <- which(aux[,3] >= alpha)[1]
    if (aux[(i-2),3] == aux[(i-1),3]) {
      if ((x[aux[(i-2),2]]-x[aux[(i-2),1]]) > (x[aux[(i-1),2]]-x[aux[(i-1),1]])) {
        result1 <- c(x[aux[(i-1),1]],x[aux[(i-1),2]],aux[(i-1),3])
      } else {
        result1 <- c(x[aux[(i-2),1]],x[aux[(i-2),2]],aux[(i-2),3])
      }      
    } else {
      result1 <- c(x[aux[(i-1),1]],x[aux[(i-1),2]],aux[(i-1),3])
    }
    if (aux[i,3] == aux[(i+1),3]) {
      if ((x[aux[i,2]]-x[aux[i,1]]) > (x[aux[(i+1),2]]-x[aux[(i+1),1]])) {
        result2 <- c(x[aux[(i+1),1]],x[aux[(i+1),2]],aux[(i+1),3])
      } else {
        result2 <- c(x[aux[i,1]],x[aux[i,2]],aux[i,3])
      }      
    } else {
      result2 <- c(x[aux[i,1]],x[aux[i,2]],aux[i,3])
    }

    cint1 <- result1[1:2]
    attr(cint1,"conf.level") <- result1[3]
    cint2 <- result2[1:2]
    attr(cint2,"conf.level") <- result2[3]
    run.time <- (proc.time()[3]-initial.time)/60
    out <- list(cint1=cint1, cint2=cint2, alpha=alpha, run.time=run.time)
  } else {
    if (all(!is.na(aux)))  {
      cint1 <- c(x[aux[1]],x[aux[2]])
      attr(cint1,"conf.level") <- aux[3]
      run.time <- (proc.time()[3]-initial.time)/60
      out <- list(cint1=cint1, cint2=NULL, alpha=alpha, run.time=run.time)

      msg <- paste("\n","  It was not possible to find a confidence interval \n",
                   "  with significance level bigger than or equal to alpha.\n",
                   "  Probably reason: sample size is small.",sep="")
      warning(msg)
    }
  }
  out$call <- Call

  class(out) <- "conf.interval"
  return(out)
}

print.conf.interval <- function(x, ...) {
  print(x$call)
  cat("\n",sep="")
  cat("-----------------------------------------------------\n",sep="")
  cat("  Interval for population median.\n",sep="")
  cat("\n",sep="")
  cat("  Interval with confidence level lower than alpha:\n",sep="")
  cat("     confidence level: ",sprintf("%.4f",attr(x$cint1,"conf.level")),"\n",sep="")
  cat("  confidence interval: (",sprintf("%.4f",x$cint1[1]),", ",sprintf("%.4f",x$cint1[2]),")\n",sep="")
  cat("\n",sep="")
  cat("    distance to alpha: ",sprintf("%.4f",abs(x$alpha-attr(x$cint1,"conf.level"))),"\n",sep="")
  cat("   interval amplitude: ",sprintf("%.4f",x$cint1[2]-x$cint1[1]),"\n",sep="")
  if (!is.null(x$cint2)) {
    cat("\n",sep="")
    cat("\n",sep="")
    cat("  Interval with confidence level bigger than alpha:\n",sep="")
    cat("     confidence level: ",sprintf("%.4f",attr(x$cint2,"conf.level")),"\n",sep="")
    cat("  confidence interval: (",sprintf("%.4f",x$cint2[1]),", ",sprintf("%.4f",x$cint2[2]),")\n",sep="")
    cat("\n",sep="")
    cat("    distance to alpha: ",sprintf("%.4f",abs(x$alpha-attr(x$cint2,"conf.level"))),"\n",sep="")
    cat("   interval amplitude: ",sprintf("%.4f",x$cint2[2]-x$cint2[1]),"\n",sep="")
  }
  cat("\n",sep="")
  cat("     Total time spent: ", sprintf("%.2f",x$run.time)," min \n",sep="")
  cat("-----------------------------------------------------\n",sep="")
  cat("\n",sep="")
}


################################################################################
# Evaluate the confidence statement for
# the ordered population median.
conf.statement <- function(data,verbose=TRUE) {
  #    data: is a list with all groups to be test,
  #          where the first element of the list is considered
  #          as the lower population quantile, the second elemente of 
  #          the list is considered as the second lower population 
  #          quantile, and so...
  # verbose: if TRUE, display results on screen.
  Call <- match.call()

  # Conditions to data
  if (!is.list(data))
    stop("data is not a list")
  if (length(data) < 2)
    stop("data must have a minimum of 2 groups")
  for (i in 1:length(data)) {
    if (!is.vector(data[[i]]))
      stop(paste("data[[",i,"]] is not a vector",sep=""))
    if (!is.numeric(data[[i]]))
      stop(paste("data[[",i,"]] is not numeric",sep=""))
    if (length(data[[i]]) <= 2)
      stop(paste("Sample size of data[[",i,"]] is too small (n <= 2)",sep=""))
  }
  q <- 0.5
  # Conditions to q
  if ((!is.numeric(q)) && (length(q) == 1))
    stop("q is not numeric or have length bigger than 1")
  if ((q <= 0) || (q >= 1))
    stop("q must be a value in (0,1)")
  # Conditions to verbose
  if ((verbose != TRUE) && (verbose != FALSE))
    stop("verbose must be TRUE or FALSE")

  m <- rep(NA,length(data))
  for (i in 1:length(data)) {
    data[[i]] <- sort(data[[i]])
    m[i] <- length(data[[i]])
  }

  if (verbose) {
    if (prod(m) > 10^5)
      cat("Evaluating confidence statement. It may take some minutes...\n")
  }
  initial.time <- proc.time()[3]

  pr <- matrix(NA,max(m),length(data))
  pr[1:m[1],1] <- pbinom(seq(0,m[1]-1,1),m[1],q)
  if (length(pr[1,]) > 2) {
    for (i in 2:(length(data)-1)) {
      pr[1:m[i],i] <- 1-pbinom(seq(0,m[i]-1,1),m[i],q)
    }
  }
  pr[1:m[length(m)],length(pr[1,])] <- 1-pbinom(seq(0,m[i]-1,1),m[i],q)

  fb <- function(x,v) {
    aux <- which(v > x)
    if (length(aux) != 0) {
      return(aux[1])
    } else {
      return(NA)
    }
  }

  e1 <- new.env()
  assign("aux", NULL, envir=e1)
  assign("r1", rep(0,2*(length(data)-1)+1), envir=e1)
  assign("r2", rep(0,2*(length(data)-1)+1), envir=e1)
  assign("k", 0, envir=e1)
  assign("pos", rep(0,2*(length(data)-1)), envir=e1)
  assign("it", 1, envir=e1)

  recursive <- function(a,j) {
    if (j < length(m)) {
      for (i in a:m[j]) {
        a1 <- fb(data[[j]][i], data[[j+1]])
        if (verbose) {
          run.time <- (proc.time()[3]-initial.time)/60
          if (run.time > e1$it) {
            e1$it <- e1$it+1
            cat("Time elapsed: ", sprintf("%.0f",floor(run.time))," min, so far...\n",sep="")
          }
        }
        if (!is.na(a1)) {
          eval(parse(text=paste("e1$aux$a",j," <- c(",i,",",a1,")",sep="")))
          if ((j == length(m)-1) && (length(e1$aux) == length(m)-1)) {
            e1$k <- e1$k+1
            aux <- data[[1]][e1$aux[[1]][1]]
            aux2 <- pr[e1$aux[[1]][1],1]
            aux3 <- e1$aux[[1]][1]
            if (length(m) > 2) {
              for (i1 in 2:(length(m)-1)) {
                aux <- c(aux,data[[i1]][e1$aux[[i1-1]][2]],data[[i1]][e1$aux[[i1]][1]])
                aux2 <- aux2*(pr[e1$aux[[i1-1]][2],i1]-pr[e1$aux[[i1]][1],i1])
                aux3 <- c(aux3,e1$aux[[i1-1]][2],e1$aux[[i1]][1])
              }
            }
            aux <- c(aux,data[[length(m)]][e1$aux[[length(m)-1]][2]])
            aux2 <- aux2*pr[e1$aux[[length(m)-1]][2],length(pr[1,])]
            aux3 <- c(aux3,e1$aux[[length(m)-1]][2])

            e1$r1 <- c(aux,aux2)
            if ((all(!is.na(e1$r1))) && (e1$r1[length(e1$r1)] > e1$r2[length(e1$r2)])) {
              e1$r2 <- e1$r1
              e1$pos <- aux3
            }
            rm(aux)
          } else {
            b <- recursive(a1+1,j+1)
          }
        }
      }
    } 
  }
  recursive(1,1)

  out <- NULL
  out$call <- Call
  out$statement.level <- NA

  result <- e1$r2 
  if (result[length(result)] >= 0) {
    out$statement.level  <- result[length(result)]

    i <- 1
    eval(parse(text=paste("out$stat.order.",i," <- ",e1$pos[i],sep="")))
    eval(parse(text=paste("out$conf.statement.",i," <- ",result[i],sep="")))
    if (length(m) >= 3) {
      for (i in 2:(length(m)-1)) {
        eval(parse(text=paste("out$stat.order.",i," <- c(",e1$pos[2*(i-1)],",",e1$pos[2*(i-1)+1],")",sep="")))
        eval(parse(text=paste("out$conf.statement.",i," <- c(",result[2*(i-1)],",",result[2*(i-1)+1],")",sep="")))
      }
    }
    i <- length(result)-1
    eval(parse(text=paste("out$stat.order.",length(m)," <- ",e1$pos[i],sep="")))
    eval(parse(text=paste("out$conf.statement.",length(m)," <- ",result[i],sep="")))

    out$total.groups <- length(m)
    out$run.time <- (proc.time()[3]-initial.time)/60
  }
  class(out) <- "conf.statement"
  return(out)
}

print.conf.statement <- function(x, ...) {
  print(x$call)

  if (is.na(x$statement.level)) {
    cat("It was not possible to find a confidence statement.\n")
  } else {
    aux1 <- 0
      aux2 <- 0
    for (i in 1:x$total.groups) {
    eval(parse(text=paste("aux1 <- max(0,x$stat.order.",i,")",sep="")))
      eval(parse(text=paste("aux2 <- max(0,x$conf.statement.",i,")",sep="")))
    }
    nmc <- nchar(aux1)
    nmr <- nchar(trunc(aux2))
  
    cat("-----------------------------------------------------------\n")
    cat("Confidence Statement of ordered population medians:\n",sep="")
    cat("\n")
    cat(rep(" ",nmc+nmr-1),"                M1 < ",sep="")
    cat(rep(" ",nmr-nchar(trunc(x$conf.statement.1))),sprintf("%.4f",x$conf.statement.1),
        " = X1(",sprintf("%.0f",x$stat.order.1),")","\n",sep="")
    if (x$total.groups >= 3) {
      for (i in 2:(x$total.groups-1)) {
        eval(parse(text=paste("aux1 <- x$stat.order.",i,sep="")))
        eval(parse(text=paste("aux2 <- x$conf.statement.",i,sep="")))
  
        cat(rep(" ",nmc-nchar(aux1[1])),"X",sprintf("%.0f",i),"(",sprintf("%.0f",aux1[1]),") = ",
            rep(" ",nmr-nchar(trunc(aux2[1]))),sprintf("%.4f",aux2[1]),sep="")
        cat(" < M",sprintf("%.0f",i)," < ",sep="")
        cat(rep(" ",nmr-nchar(trunc(aux2[2]))),sprintf("%.4f",aux2[2]),
            " = X",sprintf("%.0f",i),"(",sprintf("%.0f",aux1[2]),")","\n",sep="")
      }
    }
    eval(parse(text=paste("aux1 <- x$stat.order.",x$total.groups,sep="")))
    eval(parse(text=paste("aux2 <- x$conf.statement.",x$total.groups,sep="")))
    cat(rep(" ",nmc-nchar(aux1)),"X",sprintf("%.0f",x$total.groups),
        "(",sprintf("%.0f",aux1),") = ",
        sprintf("%.4f",aux2),sep="")
    cat(" < M",sprintf("%.0f",x$total.groups),"\n",sep="")
    cat("\n")
    cat("Confidence statement: ",sprintf("%.4f",x$statement.level),"\n",sep="")
    cat("\n")
    cat("Total time spent: ", sprintf("%.2f",x$run.time)," min \n",sep="")
    cat("-----------------------------------------------------------\n")
    cat("\n")
  }
}


