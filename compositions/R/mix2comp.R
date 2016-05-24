#  ==========================================================================#  transformations from mixtures to  'compositions' classes 'aplus', 'acomp', 'rcomp', 'rplus'
#  ========================================================================# transform mixture to aplus
mix.2aplus <- function (X)
{
    Y <- X$mat
     class(Y) <- "aplus"
    Y
}

# transform mixture to acomp
mix.2acomp <- function (X)
 {
    Y <- X$mat
     class(Y) <- "acomp"
    Y
 }


# transform mixture to rcomp
mix.2rcomp <- function (X)
{
    Y <- X$mat
     class(Y) <- "rcomp"
    Y
}


# transform mixture to rplus
mix.2rplus <- function (X)
{
    Y <- X$mat
     class(Y) <- "rplus"
    Y
}


# transform mixture to rmult
mix.2rmult <- function (X)
{
    Y <- X$mat
     class(Y) <- "rmult"
    Y
}


mix.Read <- function(file, eps=1e-6){
  mix.Check <- function(m, eps=1e-6){
    m$sum <- NA
    if(any(m$mat < 0)) m$sta <- -2              # matrix contains negative elements
    else{ s <- drop(m$mat %*% rep(1, ncol(m$mat)))
          if( any(s<=0)) m$sta <- -1          #  zero sum row exists
          else{ if(any(m$mat <eps)) {m$sta <-  0      #  matrix contains zero elements
                                     if(((max(s)-min(s))/max(s)) < eps ){m$sum <- mean(s)}  # with constant row sum
                                   }
          else { if(((max(s)-min(s))/max(s)) < eps ){m$sum <- mean(s)
                                                     if (abs(m$sum-1) < eps) m$sta <- 3          # normalized
                                                     else m$sta <- 2       # mixture with constant row sum
                                                   } else m$sta <- 1 }        #  matrix contains positive element, rows with different row sum(s)
              }
        }
    m
  }

  con <- file(file,"r")
  t <- readLines(con,n=1)
  m <- as.matrix(read.table(con, header = TRUE, sep = "", quote = "\"'",
       dec = ".", row.names = 1, as.is = FALSE, na.strings = "NA",
       colClasses = NA, nrows = -1, skip = 0, check.names = TRUE,
       strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#"))
  close(con)
  mix.Check(structure(list(tit=t,sum=NA,sta=NA,mat=m),class=c("mixture")),eps)
}

