Biograph.msm <-
function(Bdata) 
{# remove records where start > end
  if (is.null(attr(Bdata,"param"))) print ("Biograph.msm: Parameters missing. Run Parameters . . . . ",quote=FALSE)
# require (reshape)
z <- check.par(Bdata)
namstates <-  attr(Bdata,"param")$namstates
ncovariates  <- locpath(Bdata) - 7
Bdata2 <- subset (Bdata,Bdata$start <= Bdata$end)
attr(Bdata2,"format.date") <- format.in <- attr(Bdata,"format.date")
attr(Bdata2,"param") <- attr(Bdata,"param")
tmat <- attr(Bdata,"param")$tmat
z <- check.par (Bdata2)
locpat <- locpath(Bdata2)
z88 <- ifelse (format.in=="CMC",12,1)
print (". . . . .  running reshape  . . . . . ",quote=FALSE)
nn <- ncol(Bdata2)

Bdata2$start2 <- Bdata2$start
zx <- reshape(Bdata2, idvar = "ID", varying = list(c(3, 
        (locpat + 1):nn, 4)), v.names = "date", direction = "long",drop = NULL)
print (" . . .  . Adjusting long format  . . . . ",quote=FALSE)
zx2 <- zx[!is.na(zx$date),]
D <- zx2[do.call(order,list(zx2$ID,zx2$date)),] # sort by 2 variables
D$age <- (D$date-D$born)/z88
D$ns <- nchar(D$path)
D$time <- ifelse (D$time==max(D$time),D$ns+1,D$time)   # time = line number of episode in trajectory (first= from birth; last = open to censored)
D$OR <- ifelse (D$date==D$start2,"#",ifelse (D$time >D$ns, substr(D$path,D$ns,D$ns), substr(D$path,(D$time-1),(D$time-1))))  
D$DES <-   ifelse (D$time > D$ns, "cens",substr(D$path,(D$time),(D$time)))
nnz <- ncovariates+9+3  # = DES
tmat <- attr(Bdata2,"param")$tmat
D$trans <- apply(D,1,function (x) {ifelse (x[nnz]=="cens",grep(x[nnz-1],namstates),tmat[grep(x[nnz-1],namstates),grep(x[nnz],namstates)])})
D$firstobs <- ifelse (is.na(D$trans),1,0)
D$stateN <- ifelse (D$DES=="cens",D$OR,D$DES)
D$state <- match (D$stateN,namstates)
attr(D,"param") <- attr(Bdata,"param")
attr(D,"format.date") <- attr(Bdata,"format.date")
return (D)
}
