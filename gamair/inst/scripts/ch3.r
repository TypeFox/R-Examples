# R code for chapter 3 of Wood (2006) "GAMs: An Introduction with R"

## 3.2.1 Representing a smooth function: Regression splines
#### Using the cubic spline basis
par(mfrow=c(1,2))
size <- c(1.42,1.58,1.78,1.99,1.99,1.99,2.13,2.13,2.13,
          2.32,2.32,2.32,2.32,2.32,2.43,2.43,2.78,2.98,2.98)
wear <- c(4.0,4.2,2.5,2.6,2.8,2.4,3.2,2.4,2.6,4.8,2.9,
          3.8,3.0,2.7,3.1,3.3,3.0,2.8,1.7)
x <- size-min(size); x <- x/max(x)
plot(x,wear,xlab="Scaled engine size",ylab="Wear index")

rk<-function(x,z) # R(x,z) for cubic spline on [0,1]
{ ((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4-
  ((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2+7/240)/24
}

spl.X<-function(x,xk)
# set up model matrix for cubic penalized regression spline
{ q<-length(xk)+2  # number of parameters
  n<-length(x)     # number of data
  X<-matrix(1,n,q) # initialized model matrix
  X[,2]<-x         # set second column to x
  X[,3:q]<-outer(x,xk,FUN=rk) # and remaining to R(x,xk)
  X
}

xk<-1:4/5      # choose some knots
X<-spl.X(x,xk) # generate model matrix
mod.1<-lm(wear~X-1) # fit model
xp<-0:100/100    # x values for prediction
Xp<-spl.X(xp,xk) # prediction matrix
lines(xp,Xp%*%coef(mod.1)) # plot fitted spline

## 3.2.2 Controlling the degree of smoothing with penalized regression splines

spl.S<-function(xk)
# set up the penalized regression spline penalty matrix,
# given knot sequence xk
{ q<-length(xk)+2;S<-matrix(0,q,q) # initialize matrix to 0
  S[3:q,3:q]<-outer(xk,xk,FUN=rk)  # fill in non-zero part
  S
}

mat.sqrt<-function(S) # A simple matrix square root
{ d<-eigen(S,symmetric=TRUE)
  rS<-d$vectors%*%diag(d$values^0.5)%*%t(d$vectors)
}

prs.fit<-function(y,x,xk,lambda)
# function to fit penalized regression spline to x,y data,
# with knots xk, given smoothing parameter, lambda.
{ q<-length(xk)+2           # dimension of basis
  n<-length(x)              # number of data
  # create augmented model matrix ....
  Xa <- rbind(spl.X(x,xk),mat.sqrt(spl.S(xk))*sqrt(lambda))
  y[(n+1):(n+q)]<-0  # augment the data vector
  lm(y~Xa-1) # fit and return penalized regression spline
}

xk<-1:7/8  # choose some knots
mod.2<-prs.fit(wear,x,xk,0.0001) # fit pen. reg. spline
Xp<-spl.X(xp,xk) # matrix to map params to fitted values at xp
plot(x,wear);lines(xp,Xp%*%coef(mod.2)) # plot data & spl. fit

## 3.2.3 Choosing the smoothing parameter, \lambda: Cross validation

lambda<-1e-8;n<-length(wear);V<-rep(0,60)
for (i in 1:60)            # loop through smoothing parameters
{ mod<-prs.fit(wear,x,xk,lambda)     # fit model, given lambda
  trA<-sum(influence(mod)$hat[1:n])  # find tr(A)
  rss<-sum((wear-fitted(mod)[1:n])^2)# residual sum of squares
  V[i]<-n*rss/(n-trA)^2              # obtain GCV score
  lambda<-lambda*1.5                 # increase lambda
}
plot(1:60,V,type="l",main="GCV score",xlab="i")  # plot score

i<-(1:60)[V==min(V)]                 # extract index of min(V)
mod.3<-prs.fit(wear,x,xk,1.5^(i-1)*1e-8)  # fit optimal model
Xp<-spl.X(xp,xk)                          # .... and plot it
plot(x,wear);lines(xp,Xp%*%coef(mod.3))

## 3.3.1 Penalized regression spline representation of an additive model

am.setup<-function(x,z,q=10)
# Get X, S_1 and S_2 for a simple 2 term AM
{ # choose knots ...
  xk <- quantile(unique(x),1:(q-2)/(q-1))
  zk <- quantile(unique(z),1:(q-2)/(q-1))
  # get penalty matrices ...
  S <- list()
  S[[1]] <- S[[2]] <- matrix(0,2*q-1,2*q-1)
  S[[1]][2:q,2:q] <- spl.S(xk)[-1,-1]
  S[[2]][(q+1):(2*q-1),(q+1):(2*q-1)] <- spl.S(zk)[-1,-1]
  # get model matrix ...
  n<-length(x)
  X<-matrix(1,n,2*q-1)
  X[,2:q]<-spl.X(x,xk)[,-1]           # 1st smooth
  X[,(q+1):(2*q-1)]<-spl.X(z,zk)[,-1] # 2nd smooth
  list(X=X,S=S)
}

fit.am<-function(y,X,S,sp)
# function to fit simple 2 term additive model
{ # get sqrt of total penalty matrix ...
  rS <- mat.sqrt(sp[1]*S[[1]]+sp[2]*S[[2]])
  q.tot <- ncol(X)                # number of params
  n <- nrow(X)                    # number of data
  X1 <- rbind(X,rS)               # augmented X
  y1 <- c(y,rep(0,q.tot))         # augment data
  b<-lm(y1~X1-1)                  # fit model
  trA<-sum(influence(b)$hat[1:n]) # tr(A)
  norm<-sum((y-fitted(b)[1:n])^2) # RSS
  list(model=b,gcv=norm*n/(n-trA)^2,sp=sp)
}

data(trees)
rg <- range(trees$Girth)
trees$Girth <- (trees$Girth - rg[1])/(rg[2]-rg[1])
rh <- range(trees$Height)
trees$Height <- (trees$Height - rh[1])/(rh[2]-rh[1])

am0 <- am.setup(trees$Girth,trees$Height)

sp<-c(0,0)   # initialize smoothing parameter (s.p.) array
for (i in 1:30) for (j in 1:30)      # loop over s.p. grid
{ sp[1]<-1e-5*2^(i-1);sp[2]<-1e-5*2^(j-1) # s.p.s
  b<-fit.am(trees$Volume,am0$X,am0$S,sp)  # fit using s.p.s
  if (i+j==2) best<-b else                # store 1st model
  if (b$gcv<best$gcv) best<-b             # store best model
}
best$sp   # GCV best smoothing parameters found

# plot fitted against data ...
plot(trees$Volume,fitted(best$model)[1:31],
     xlab="Fitted Volume",ylab="Actual Volume")
# evaluate and plot f_1 against Girth ...
b<-best$model
b$coefficients[1]<-0     # zero the intercept
b$coefficients[11:19]<-0 # zero the second smooth coefs
f0<-predict(b)           # predict f_1 only, at data values
plot(trees$Girth,f0[1:31],xlab="Scaled Girth",
     ylab=expression(hat(f[1])))


## 3.4 Generalized Additive Models

fit.gamG<-function(y,X,S,sp)
# function to fit simple 2 term generalized additive model
# Gamma errors and log link
{ # get sqrt of combined penalty matrix ...
  rS <- mat.sqrt(sp[1]*S[[1]]+sp[2]*S[[2]])
  q <- ncol(X)   # number of params
  n <- nrow(X)   # number of data
  X1 <- rbind(X,rS)     # augmented model matrix
  eta <- log(y)         # initialize linear predictor
  norm <- 0;old.norm <- 1   # initialize convergence control
  while (abs(norm-old.norm)>1e-4*norm)  # repeat un-converged
  { mu <- exp(eta)          # fitted values
    z <- (y-mu)/mu + eta    # pseudodata (recall w_i=1, here)
    z[(n+1):(n+q)] <- 0     # augmented pseudodata
    m <- lm(z~X1-1)         # fit penalized working model
    b <- m$coefficients     # current parameter estimates
    eta <- (X1%*%b)[1:n]    # 'linear predictor'
    trA <- sum(influence(m)$hat[1:n]) # tr(A)
    old.norm <- norm        # store for convergence test
    norm <- sum((z-fitted(m))[1:n]^2) # RSS of working model
  }
  list(model=m,gcv=norm*n/(n-trA)^2,sp=sp)
}
