"contrMatRatio" <-
function(n, type="Tukey", base=1)
{

# check:

type<-match.arg(arg=type, choices=c("Dunnett", "Tukey", "Sequen", "AVE", "GrandMean", "Changepoint", "Marcus", "McDermott", "Williams", "UmbrellaWilliams"))

if (length(n) < 2) 
 {stop("less than 2 groups")}

if(!is.numeric(n))
 {stop("n must be a numeric vector")}
  
if (any(n < 2)) 
 {stop("less than 2 observations in at least one group")}

k <- length(n)

if (base < 1 || base > k) 
 {stop("base is not between 1 and ", k)}

# define values 

numC <- c()
denC <- c()
rnames <- c()

if (!is.null(names(n))) 
        varnames <- names(n)
    else varnames <- 1:length(n)

kindx <- 1:k

if(type=="Dunnett")

{
 for (i in kindx[-base])
  {numC <- rbind(numC, as.numeric(kindx == i) )}

   denC <- matrix(0, ncol=k, nrow=k-1 )
   denC[,base]<-1

rnames <- paste(varnames[kindx[-base]], "/", varnames[base], 
            sep = "")
}

 if(type=="Tukey")
{
  for (i in 1:(k - 1))
   {
    for (j in (i + 1):k)
     {
      numC <- rbind(numC, as.numeric(kindx == j))
      denC <- rbind(denC, as.numeric(kindx == i))

      rnames <- c(rnames, paste(varnames[j], "/", varnames[i], sep = ""))    
     }
   }
}

if(type=="Sequen") 
{
 for (i in 2:k)
  {
   numC <- rbind(numC, as.numeric(kindx == i) )
   denC <- rbind(denC, as.numeric(kindx == i - 1) )

   rnames <- c(rnames, paste(varnames[i], "/", varnames[i - 1], sep = ""))
  }
}

if(type=="Williams")
 {
 if(k==2)
  {
   numC<-matrix(c(0,1), ncol=2)
   denC<-matrix(c(1,0), ncol=2)
   rnames<-"C1"
  }
  else{
   for (i in 1:(k - 2)) {
   help <- c(0, rep(0, k - i - 1), n[(k - i + 1):k]/sum(n[(k - i + 1):k]))
   numC <- rbind(numC, help)

   denC <- rbind(denC, c(1, rep(0, times=k-1)))
  }
   help <- c(0, n[2:k]/sum(n[2:k]))
   numC <- rbind(numC, help)
   denC <- rbind(denC, c(1, rep(0, times=k-1)))

   rnames <- c(rnames, paste("C", 1:nrow(numC), sep = ""))
 }
}  

if(type=="UmbrellaWilliams")
{
denC <- matrix( c( rep(1,(k)*(k-1)/2 ), rep(0,(k)*((k-1)^2)/2)), ncol=k, byrow=FALSE )

numC <- c()

 for(j in 1:(k-1))
   {
    for(i in 1:(k - j))
      {
       helper <- c(0, rep(0, k - i - j),
              n[((k - i + 1):k)-(j-1)]/sum(n[((k - i + 1):k)-(j-1)]), rep(0, j-1))
       numC <- rbind(numC, helper)
      }
   }
rnames <- paste("C", 1:nrow(numC), sep = "")
}

if(type=="Changepoint") 
 {
  for (i in 1:(k - 1))
   {
    helpnum <- c(rep(0,times=i), n[(i + 1):k]/sum(n[(i + 1):k]))
    helpden <- c(n[1:i]/sum(n[1:i]), rep(0, times=k-i) )

    numC <- rbind(numC, helpnum)
    denC <- rbind(denC, helpden)
        }
    rnames <- c(rnames, paste("C", 1:nrow(numC), sep = ""))
}

if(type=="AVE")
{

 if(k==2)
  {
   numC<-matrix(c(0,1), ncol=2)
   denC<-matrix(c(1,0), ncol=2)
   rnames<-"C1"
  }
 else{
   helpnum <- c(1, rep(0, times=k-1))
   helpden <- c(0, n[2:k]/sum(n[2:k]))
   numC <- rbind(numC, helpnum)
   denC <- rbind(denC, helpden)


        for (i in 2:(k - 1)) {
            x <- sum(n[1:(i - 1)]) + sum(n[(i + 1):k])
            helpnum <- c( as.numeric(kindx == i))
            helpden <- c(n[1:(i - 1)]/x, 0, n[(i + 1):k]/x)
            numC <- rbind(numC, helpnum )
            denC <- rbind(denC, helpden)
        }

        helpnum <- c(rep(0, times=k-1), 1)
        helpden <- c(n[1:(k - 1)]/sum(n[1:(k - 1)]), 0)

        numC <- rbind(numC, helpnum)
        denC <- rbind(denC, helpden)

        rnames <- paste("C", 1:nrow(numC), sep = "")
  }
}

if(type=="GrandMean")
{
 numC <- diag(1,nrow=k)
 denC <- matrix(rep(n/sum(n), times=k), byrow=TRUE, ncol=k)
 rnames <- paste(varnames, "Grand mean", sep="/")
}


if(type=="McDermott")
{
 if(k==2)
  {
   numC<-matrix(c(0,1), ncol=2)
   denC<-matrix(c(1,0), ncol=2)
   rnames<-"C1"
  }
 else{
  for (i in 1:(k - 2))
   {
    helpnum <- c(rep(0, times=i), 1, rep(0, times = k - i - 1))
    helpden <- c(n[1:i]/sum(n[1:i]),  rep(0, times= k - i ))

    numC <- rbind(numC, helpnum)
    denC <- rbind(denC, helpden)
   }

   helpnum <- c(rep(0, times=k-1), 1)
   helpden <- c(n[1:(k - 1)]/sum(n[1:(k - 1)]), 0)

   numC <- rbind(numC, helpnum)
   denC <- rbind(denC, helpden)

   rnames <- c(rnames, paste("C", 1:nrow(numC), sep = ""))
  }
}

if(type=="Marcus")
{
  cm1 <- matrix(0, nrow = k - 1, ncol = k)
  cm2 <- cm1
  for (i in 1:(k - 1))
   {
    cm1[i, (i + 1):k] <- n[(i + 1):k]/sum(n[(i + 1):k])
    cm2[i, 1:i] <- n[1:i]/sum(n[1:i])
   }

row <- k * (k - 1)/2
index <- 1
  for (i in 1:(k - 1))
   {
    for (j in 1:i) 
     {
      helpnum <- cm1[i, ]
      helpden <- cm2[j, ]
      numC <- rbind(numC, helpnum)
      denC <- rbind(denC, helpden)
      index <- index + 1
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(numC), sep = ""))
}

  rownames(numC)<-rnames
  rownames(denC)<-rnames
  colnames(numC)<-varnames
  colnames(denC)<-varnames

  out<-list(numC=numC, denC=denC, rnames=rnames)
  attr(out, which="type")<-type


return(out)
}

