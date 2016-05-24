"smdenreg" <-
function (x, verbose = FALSE, bandwidth = -1,
    maxkuipnr = 19, asympbounds = FALSE, squeezing.factor=0.9,firstlambda=10,smqeps=1/length(x),fsign=double(0),gensign=TRUE,...)
{
    nsamp <- length(x)
    lambda <- rep(firstlambda,nsamp-2)
    if (asympbounds || nsamp > max(kuipdiffbounds.x)) 
        currbounds <- kuipdiffbounds[length(kuipdiffbounds.x), ] * sqrt(max(kuipdiffbounds.x))/sqrt(nsamp)
    else {
        currbounds <- double(maxkuipnr)
        for (i in 1:maxkuipnr) currbounds[i] <- approx(kuipdiffbounds.x, 
            kuipdiffbounds[, i], nsamp, rule = 2)$y
    }
    if (maxkuipnr > dim(kuipdiffbounds)[2]) 
        stop("maxkuipnr is too large")
    dataemp <- c(0, rep(1/(nsamp - 1), nsamp - 1))
    x <- sort(x)
    fdist.y <- c(seq(0, 1, length = nsamp))
    if(gensign)
      fsign <- gendensgn(x)
    repeat {
        tmp <- smqden(x,lambda=lambda,eps=smqeps,fsign=fsign)
        x.string <- tmp$y
        lastunif <- c(0,cumsum((x[-1]-x[-nsamp])*x.string))
        if (verbose) {
            par(mfrow = c(3, 1))
            hist(x, 40, prob = TRUE)
            lines(tmp,col="red")
            plot(x[-c(1,nsamp)],lambda,ty="l")
        }
        if (bandwidth > 0) 
            break
        diff <- cumsum((lastunif[-1]-lastunif[-nsamp]) -1/(nsamp-1))

        tmpkkuip <- kkuip(diff, maxkuipnr)
        currkkuip <- tmpkkuip$met
        kuipinds <- c(currkkuip[1], currkkuip[-1] - currkkuip[-maxkuipnr]) > currbounds + 1e-08
        if (sum(kuipinds) == 0) 
            break

        lambda <- lambda * squeezing.factor

        if (verbose) {
            print("Press Enter")
            readline()
        }
    }
    list(x = tmp$x, y = x.string, nmax = findmod(x.string)$mod,trans = lastunif)
}

"smqreg" <-
function (y, thr.const = 2.5, verbose = FALSE,  
    bandwidth = -1, sigma = -1, localsqueezing = TRUE, squeezing.factor = 0.5, DYADIC=TRUE, 
    firstlambda=100,smqeps=1/length(y),fsign=double(0),gensign=TRUE,tolerance=1e-12,...) 
{
    n <- length(y)
    lambda <- rep(firstlambda,n-1)
    sigma <- mad((y[-1] - y[-n])/sqrt(2))
    if (verbose)
      print(c("sigma is ", sigma))

    if(gensign)
      fsign <- gensign(y,thr.const=thr.const,extrema.mean=TRUE,sigma=sigma,localsqueezing=localsqueezing,squeezing.factor=squeezing.factor)

    repeat {
        f <- smqnew(y=y,lambda=lambda,eps=smqeps,fsign=fsign,tolerance=tolerance)

        if (bandwidth < 0) {
            residuals <- y - f
            residuals <- residuals - mean(residuals)
            if(DYADIC) 
              residuals.wr <- multiwdwr(residuals, sqrt(thr.const * log(n)) * sigma)
            else
              residuals.wr <- nondymwdr(residuals, sqrt(thr.const * log(n)) * sigma)
        }
        if (verbose) {
            par(mfrow = c(2, 1))
            plot(y, col = "grey",...)
            lines(f, col = "red")
            lines(residuals.wr, type = "l", col = "green")
            lines(fsign-2,col="blue")
            plot(lambda,ty="b",...)
            print("Press Enter")
            dum <- readline()
        }
        if (bandwidth > 0) 
            break
        ind <- (abs(residuals.wr) > 1e-10)
        ind2 <- ind[-1] | ind[-n]
        if (sum(ind)==0)
            break
        if (localsqueezing) 
            lambda[ind2] <- lambda[ind2] * squeezing.factor
        else lambda <- lambda * squeezing.factor
        if(min(lambda)<1e-08)
          {
          par(mfrow = c(2, 1))
          plot(y, col = "grey",...)
          lines(f, col = "red")
          lines(residuals.wr, type = "l", col = "green")
          lines(fsign-2,col="blue")
          plot(lambda,ty="b",...)

          print("ERROR, This should not happen")
          break
          }
    }
    list(y = f, nmax = findmod(f)$mod,sigma = sigma)
}

"gensign" <-
function(y,...)
{
n <- length(y)
tmp <- pmreg(y,...)
fsign <- rep(-2,n-1)

kx <- tmp$knotsind[-length(tmp$knotsind)]
ky <- tmp$y[kx]

klen <- length(kx)

if(klen==1)
  {
  print("SPECIAL CASE! NOT YET DONE!")
  stop()
  }

for(i in 1:(kx[2]-1))
  fsign[i] <- (ky[2] > ky[1])

for(j in 2:(klen-1))
  {
  if((ky[j-1]<ky[j])&&(ky[j]<ky[j+1]))
    for(i in kx[j]:(kx[j+1]-1))
      fsign[i] <- 1
  else
  if((ky[j-1]>ky[j])&&(ky[j]>ky[j+1]))
    for(i in kx[j]:(kx[j+1]-1))
      fsign[i] <- 0
  else
  if((ky[j-1]<ky[j])&&(ky[j]>ky[j+1]))
    {
    for(i in kx[j]:floor(0.5*(kx[j]+kx[j+1])))
      fsign[i] <- 1
    for(i in floor(0.5*(kx[j]+kx[j+1])):(kx[j+1]-1))
      fsign[i] <- 0
    }  
  else
    {
    for(i in kx[j]:floor(0.5*(kx[j]+kx[j+1])))
      fsign[i] <- 0
    for(i in floor(0.5*(kx[j]+kx[j+1])):(kx[j+1]-1))
      fsign[i] <- 1
    }  
  }

for(i in kx[klen]:(n-1))
  fsign[i] <- (ky[klen] > ky[klen-1])

fsign

}

"gendensgn" <-
function(x,...)
{
n <- length(x)
tmp <- pmden(x,...)
fsign <- rep(-2,n-2)

kx <- tmp$ind[-length(tmp$ind)]
ky <- tmp$y[kx]

klen <- length(kx)

if(klen==1)
  {
  print("SPECIAL CASE! NOT YET DONE!")
  stop()
  }

for(i in 1:(kx[2]-1))
  fsign[i] <- (ky[2] > ky[1])

for(j in 2:(klen-1))
  {
  if((ky[j-1]<ky[j])&&(ky[j]<ky[j+1]))
    for(i in kx[j]:(kx[j+1]-1))
      fsign[i] <- 1
  else
  if((ky[j-1]>ky[j])&&(ky[j]>ky[j+1]))
    for(i in kx[j]:(kx[j+1]-1))
      fsign[i] <- 0
  else
  if((ky[j-1]<ky[j])&&(ky[j]>ky[j+1]))
    {
    for(i in kx[j]:floor(0.5*(kx[j]+kx[j+1])))
      fsign[i] <- 1
    for(i in floor(0.5*(kx[j]+kx[j+1])):(kx[j+1]-1))
      fsign[i] <- 0
    }  
  else
    {
    for(i in kx[j]:floor(0.5*(kx[j]+kx[j+1])))
      fsign[i] <- 0
    for(i in floor(0.5*(kx[j]+kx[j+1])):(kx[j+1]-1))
      fsign[i] <- 1
    }  
  }

for(i in kx[klen]:(n-2))
  fsign[i] <- (ky[klen] > ky[klen-1])

fsign

}

"smqden" <-
function(x,lambda,eps=1,fsign=double(0))
{
n <- length(x)
if(length(eps)==1)
  eps <- rep(eps,n-1)
xeval <- 0.5*(x[-1]+x[-n])
if(length(fsign)==0)
  fsign=rep(-1,n-1)
tmp <- .C("smqden",x=as.double(x),xeval=as.double(xeval),f=double(n-1),as.integer(n),as.double(lambda),as.double(eps),as.integer(fsign),PACKAGE="ftnonpar")


list(x=xeval,y=tmp$f)
}

"smqnew" <-
function(y,lambda,eps=1,fsign,tolerance=1e-12)
{
n <- length(y)
if(length(fsign)==0)
  fsign<-rep(-1,n-1)
if(length(eps)==1)
  eps <- rep(eps,n-1)
tmp <- .C("smqnew",y=as.double(y),f=double(n),as.integer(n),as.double(lambda),as.double(eps),as.integer(fsign),as.double(tolerance),PACKAGE="ftnonpar")


tmp$f
}

"nondymwdr" <-
function (y, thresh, firstwidth = 1) 
{
    .C("nondymwdwr", y = as.double(y), as.integer(length(y)), 
        as.double(thresh), as.integer(firstwidth),PACKAGE="ftnonpar")$y
}

"findmod" <-
function (y,tol=1e-05)
{
n <- length(y)
x <- 1:n
ind <- abs(y[-1]-y[-n]) < tol
if(sum(ind)>0)
  {
  y <- y[- (1:(n-1))[ind]]
  x <- x[- (1:(n-1))[ind]]
  n <- length(y)
  }
ind <- (y[-c(1,2)]<y[-c(1,n)])&(y[-c(1,n)]>y[-c(n-1,n)])|(y[-c(1,2)]>y[-c(1,n)])&(y[-c(1,n)]<y[-c(n-1,n)])
list(modality=sum(ind),x=x[ind])
}

