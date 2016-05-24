## taken from numerics/Density/explicit.densities.r
## updated on January 9th 2003

f <- function(d, rowwidth, colwidth, ax, bx, i, j)
  {

    ## calculates the density of the distance (not its sqaure)
    ## one point is in rectangle of sides (rw by cw)
    ## rectangle is at location (i=1,j=1)
    ## rectangles have rs and cs seperations between rows and columns
    ## second point is in rectangle (i,j)

    ## depending on the values of i and j there are three different algorithms
    ## to compute the density of uu, the first is a formula calculated by Ghosh
    ## the other two take his technique and apply it more generally

    ## all three densities fcase1, fcase2 and fcase3 are for the square of the distance
    ## that is why we pass uu^2 and multiply the result by 2*uu

    ## only pass vectors uu with values on the postive real axis

    ## pass the dimensions of the plot and the coordinates of the lower left hand corner of the second plot, see report - July 20th 2003

    ## TP for temporary parameter
    uu <- d
    TPa <- rowwidth
    rw <- rowwidth
    TPb <- colwidth
    cw <- colwidth
    TPr <- ax
    TPs <- bx

    if(i != 1 && j != 1) {
      ret <- fcase3(uu^2,TPa=TPa,TPb=TPb,TPr=TPr,TPs=TPs)*2*uu
    } else if(j==1 && i!=1) {
      ret <- fcase2(uu^2,TPa=TPa,TPb=TPb,TPr=TPr)*2*uu
    } else if(j!=1 && i==1) {
      ret <- fcase2(uu^2,TPa=cw,TPb=rw,TPr=bx)*2*uu
    } else if(i==1 && j==1) {
      ret <- fcase1(uu^2,TPa=TPa,TPb=TPb)*2*uu
    }
    
    ret
  }

h <- function(uVect,alpha)
  {
    uTEMP <- uVect[uVect<=alpha^2]
    ret <- rep(alpha*pi/2,length(uTEMP))
    uTEMP <- uVect[uVect>alpha^2]
    ret <- c(ret,2*sqrt(uTEMP-alpha^2)-alpha*asin(1-2*alpha^2/uTEMP))
    ret
  }

## case 1

fcase1 <- function(uVect,TPa,TPb)
  {
    ## density for distance squared
    ## Ghosh, Bull Calcutta
    ## assumes a <= b

    lim1 <- 0
    lim2 <- min(TPa^2,TPb^2)
    lim3 <- max(TPa^2,TPb^2)
    lim4 <- TPa^2+TPb^2

    ret <- NULL
    uTEMP <- uVect[uVect < lim1]
    ret <- c(ret,rep(0,length(uTEMP)))
    uTEMP <- uVect[uVect >= lim1 & uVect < lim2]
    ret <- c(ret,TPa*TPb*pi-2*sqrt(uTEMP)*(TPa+TPb)+uTEMP)
    uTEMP <- uVect[uVect >= lim2 & uVect < lim3]
    ret <- c(ret,TPa*TPb*pi/2-lim2-2*max(TPa,TPb)*sqrt(uTEMP)+max(TPa,TPb)*h(uTEMP,min(TPa,TPb)))
    uTEMP <- uVect[uVect >= lim3 & uVect < lim4]
    ret <- c(ret,TPa*h(uTEMP,TPb)+TPb*h(uTEMP,TPa)-lim2-lim3-uTEMP)
    uTEMP <- uVect[uVect >=lim4]
    ret <- c(ret,rep(0,length(uTEMP)))
    
    ret <- ret/(lim2*lim3)
    ret
  }

## case 2
f2.1 <- function(u,TPa,TPb,TPr)
  {
    -(TPr-TPa)*TPb*pi/2-(TPr-TPa)^2-u+2*(TPr-TPa)*sqrt(u)+TPb*h(u,TPr-TPa)
  }

f2.2 <- function(u,TPa,TPb,TPr)
  {
    ((TPr+TPa)*TPb*pi/2 + TPr^2 + 2*TPr*TPa-TPa^2) + u - 2*(TPr+TPa)*sqrt(u) - 2*TPb*h(u,TPr)+TPb*h(u,TPr-TPa)
  }

f2.3 <- function(u,TPa,TPb,TPr)
  {
    -2*TPa^2-2*TPb*h(u,TPr)+TPb*h(u,TPr-TPa)+TPb*h(u,TPr+TPa)
  }

f2.4 <- function(u,TPa,TPb,TPr)
  {
    (TPr^2-2*TPr*TPa-TPa^2+TPb^2)+u-(TPr-TPa)*h(u,TPb)-2*TPb*h(u,TPr)+TPb*h(u,TPr+TPa)
  }

f2.5 <- function(u,TPa,TPb,TPr)
  {
    -(TPr+TPa)^2-TPb^2-u+(TPr+TPa)*h(u,TPb)+TPb*h(u,TPr+TPa)
  }

f2.6 <- function(u,TPa,TPb,TPr)
  {
    (TPb^2+(TPr+TPa)*TPb*pi/2+2*TPr^2)+2*u-2*(TPr+TPa)*sqrt(u)-2*TPb*h(u,TPr)-(TPr-TPa)*h(u,TPb)
  }

f2.7 <- function(u,TPa,TPb,TPr)
  {
    ((TPr+TPa)*TPb*pi/2-TPb^2)-2*(TPr+TPa)*sqrt(u)+(TPr+TPa)*h(u,TPb)
  }

f2.8 <- function(u,TPa,TPb,TPr)
  {
    (-(TPr-TPa)*TPb*pi/2+TPb^2)+2*(TPr-TPa)*sqrt(u)-(TPr-TPa)*h(u,TPb)
  }

fcase2 <- function(uVect,TPa,TPb,TPr)
  {
    ret <- NULL

    ## in the notes these are u1, u2, u3, u4, u5, u6
    ## but this will be confusing with uVect
    lim1 <- (TPr-TPa)^2
    lim2 <- TPr^2
    lim3 <- (TPr+TPa)^2
    lim4 <- (TPr-TPa)^2+TPb^2
    lim5 <- TPr^2+TPb^2
    lim6 <- (TPr+TPa)^2+TPb^2
    
    if(lim3<=lim4)
      {
        uTEMP <- uVect[uVect < lim1]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim1 & uVect < lim2]
        ret <- c(ret,f2.1(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim2 & uVect < lim3]
        ret <- c(ret,f2.2(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim3 & uVect < lim4]
        ret <- c(ret,f2.3(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim4 & uVect < lim5]
        ret <- c(ret,f2.4(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim5 & uVect < lim6]
        ret <- c(ret,f2.5(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim6]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim4<lim3 && lim3<=lim5)
      {
        uTEMP <- uVect[uVect < lim1]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim1 & uVect < lim2]
        ret <- c(ret,f2.1(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim2 & uVect < lim4]
        ret <- c(ret,f2.2(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim4 & uVect < lim3]
        ret <- c(ret,f2.6(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim3 & uVect < lim5]
        ret <- c(ret,f2.4(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim5 & uVect < lim6]
        ret <- c(ret,f2.5(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim6]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim2<=lim4 && lim5<lim3)
      {
        uTEMP <- uVect[uVect < lim1]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim1 & uVect < lim2]
        ret <- c(ret,f2.1(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim2 & uVect < lim4]
        ret <- c(ret,f2.2(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim4 & uVect < lim5]
        ret <- c(ret,f2.6(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim5 & uVect < lim3]
        ret <- c(ret,f2.7(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim3 & uVect < lim6]
        ret <- c(ret,f2.5(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim6]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim4 < lim2)
      {
        uTEMP <- uVect[uVect < lim1]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim1 & uVect < lim4]
        ret <- c(ret,f2.1(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim4 & uVect < lim2]
        ret <- c(ret,f2.8(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim2 & uVect < lim5]
        ret <- c(ret,f2.6(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim5 & uVect < lim3]
        ret <- c(ret,f2.7(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim3 & uVect < lim6]
        ret <- c(ret,f2.5(uTEMP,TPa,TPb,TPr))
        uTEMP <- uVect[uVect >= lim6]
        ret <- c(ret,rep(0,length(uTEMP)))
      }
    
    ret/(2*TPa^2*TPb^2)
  }

## CASE 3

## PARALLELOGRAM 1

f11 <- function(u,TPa,TPb,TPr,TPs)
  {
    ((TPr-TPa)^2+(TPs-TPb)^2)+u-2*(TPr-TPa)*sqrt(u-(TPs-TPb)^2)-2*(TPs-TPb)*sqrt(u-(TPr-TPa)^2)+(TPr-TPa)*(TPs-TPb)*(asin(1-2*(TPr-TPa)^2/u) + asin(1-2*(TPs-TPb)^2/u))
  }

f12 <- function(u,TPa,TPb,TPr,TPs)
  {
    TPa^2+2*(TPs-TPb)*(sqrt(u-TPr^2)-sqrt(u-(TPr-TPa)^2))-(TPr-TPa)*(TPs-TPb)*(asin(1-2*TPr^2/u)-asin(1-2*(TPr-TPa)^2/u))
  }

f13 <- function(u,TPa,TPb,TPr,TPs)
  {
    -(TPr*(TPr-2*TPa)+TPs*(TPs-2*TPb))-u+2*(TPr-TPa)*sqrt(u-TPs^2)+2*(TPs-TPb)*sqrt(u-TPr^2)-(TPr-TPa)*(TPs-TPb)*(asin(1-2*TPr^2/u)+asin(1-2*TPs^2/u))
  }

f14 <- function(u,TPa,TPb,TPr,TPs)
  {
    TPb^2-2*(TPr-TPa)*(sqrt(u-(TPs-TPb)^2)-sqrt(u-TPs^2))+(TPr-TPa)*(TPs-TPb)*(asin(1-2*(TPs-TPb)^2/u)-asin(1-2*TPs^2/u))
  }

f1 <- function(uVect,TPa,TPb,TPr,TPs)
  {
    lim <- c((TPr-TPa)^2+(TPs-TPb)^2,TPr^2+(TPs-TPb)^2,(TPr-TPa)^2+TPs^2,TPr^2+TPs^2)
    ret <- NULL
    if(lim[2]<=lim[3])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[2]]
        ret <- c(ret,f11(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[3]]
        ret <- c(ret,f12(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[4]]
        ret <- c(ret,f13(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim[3] < lim[2])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[3]]
        ret <- c(ret,f11(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[2]]
        ret <- c(ret,f14(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[4]]
        ret <- c(ret,f13(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }
    
    ret/(4*TPa^2*TPb^2)
  }
    
## PARALLELOGRAM 2

f21 <- function(u,TPa,TPb,TPr,TPs)
  {
    (-TPr^2-2*TPr*TPa-(TPs-TPb)^2)-u+2*(TPr+TPa)*sqrt(u-(TPs-TPb)^2)+2*(TPs-TPb)*sqrt(u-TPr^2)-(TPr+TPa)*(TPs-TPb)*(asin(1-2*(TPs-TPb)^2/u)+asin(1-2*TPr^2/u))
  }

f22 <- function(u,TPa,TPb,TPr,TPs)
  {
    TPa^2 - 2*(TPs-TPb)*(sqrt(u-(TPr+TPa)^2)-sqrt(u-TPr^2))+(TPr+TPa)*(TPs-TPb)*(asin(1-2*(TPr+TPa)^2/u)-asin(1-2*TPr^2/u))
  }

f23 <- function(u,TPa,TPb,TPr,TPs)
  {
    ((TPr+TPa)^2+TPs^2-2*TPs*TPb)+u-2*(TPr+TPa)*sqrt(u-TPs^2)-2*(TPs-TPb)*sqrt(u-(TPr+TPa)^2)+(TPr+TPa)*(TPs-TPb)*(asin(1-2*(TPr+TPa)^2/u)+asin(1-2*TPs^2/u))
  }

f24 <- function(u,TPa,TPb,TPr,TPs)
  {
    -TPb^2+2*(TPr+TPa)*(sqrt(u-(TPs-TPb)^2)-sqrt(u-TPs^2))-(TPr+TPa)*(TPs-TPb)*(asin(1-2*(TPs-TPb)^2/u)-asin(1-2*TPs^2/u))
  }

f2 <- function(uVect,TPa,TPb,TPr,TPs)
  {
    lim <- c((TPr)^2+(TPs-TPb)^2,(TPr+TPa)^2+(TPs-TPb)^2,(TPr)^2+(TPs)^2,(TPr+TPa)^2+(TPs)^2)
    ret <- NULL
    if(lim[2]<=lim[3])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[2]]
        ret <- c(ret,f21(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[3]]
        ret <- c(ret,f22(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[4]]
        ret <- c(ret,f23(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim[3] < lim[2])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[3]]
        ret <- c(ret,f21(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[2]]
        ret <- c(ret,f24(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[4]]
        ret <- c(ret,f23(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }
    
    ret/(4*TPa^2*TPb^2)
  }

## PARALLELOGRAM 3

f31 <- function(u,TPa,TPb,TPr,TPs)
  {
    (TPr*(TPr+2*TPa)+TPs*(TPs+2*TPb))+u-2*(TPr+TPa)*sqrt(u-TPs^2)-2*(TPs+TPb)*sqrt(u-TPr^2)+(TPr+TPa)*(TPs+TPb)*(asin(1-2*TPs^2/u)+asin(1-2*TPr^2/u))
  }

f32 <- function(u,TPa,TPb,TPr,TPs)
  {
    -TPa^2+2*(TPs+TPb)*(sqrt(u-(TPr+TPa)^2)-sqrt(u-TPr^2))-(TPr+TPa)*(TPs+TPb)*(asin(1-2*(TPr+TPa)^2/u)-asin(1-2*TPr^2/u))
  }

f33 <- function(u,TPa,TPb,TPr,TPs)
  {
    (-(TPr+TPa)^2-(TPs+TPb)^2)-u+2*(TPr+TPa)*sqrt(u-(TPs+TPb)^2)+2*(TPs+TPb)*sqrt(u-(TPr+TPa)^2)-(TPr+TPa)*(TPs+TPb)*(asin(1-2*(TPr+TPa)^2/u)+asin(1-2*(TPs+TPb)^2/u))
  }

f34 <- function(u,TPa,TPb,TPr,TPs)
  {
    -TPb^2+2*(TPr+TPa)*(sqrt(u-(TPs+TPb)^2)-sqrt(u-TPs^2))-(TPs+TPb)*(TPr+TPa)*(asin(1-2*(TPs+TPb)^2/u)-asin(1-2*TPs^2/u))
  }

f3 <- function(uVect,TPa,TPb,TPr,TPs)
  {
    lim <- c((TPr)^2+(TPs)^2,(TPr+TPa)^2+(TPs)^2,(TPr)^2+(TPs+TPb)^2,(TPr+TPa)^2+(TPs+TPb)^2)
    ret <- NULL
    if(lim[2]<=lim[3])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[2]]
        ret <- c(ret,f31(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[3]]
        ret <- c(ret,f32(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[4]]
        ret <- c(ret,f33(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim[3] < lim[2])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[3]]
        ret <- c(ret,f31(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[2]]
        ret <- c(ret,f34(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[4]]
        ret <- c(ret,f33(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }
    
    ret/(4*TPa^2*TPb^2)
  }

## PARALLELOGRAM 4

f41 <- function(u,TPa,TPb,TPr,TPs)
  {
    -(TPr-TPa)^2-TPs^2-2*TPb*TPs-u+2*(TPr-TPa)*sqrt(u-TPs^2)+2*(TPs+TPb)*sqrt(u-(TPr-TPa)^2)-(TPr-TPa)*(TPs+TPb)*(asin(1-2*TPs^2/u)+asin(1-2*(TPr-TPa)^2/u))
  }

f42 <- function(u,TPa,TPb,TPr,TPs)
  {
    -TPa^2-2*(TPs+TPb)*(sqrt(u-TPr^2)-sqrt(u-(TPr-TPa)^2))+(TPr-TPa)*(TPs+TPb)*(asin(1-2*TPr^2/u)-asin(1-2*(TPr-TPa)^2/u))
  }

f43 <- function(u,TPa,TPb,TPr,TPs)
  {
    (TPs+TPb)^2+TPr^2-2*TPa*TPr+u-2*(TPr-TPa)*sqrt(u-(TPs+TPb)^2)-2*(TPs+TPb)*sqrt(u-TPr^2)+(TPr-TPa)*(TPs+TPb)*(asin(1-2*TPr^2/u)+asin(1-2*(TPs+TPb)^2/u))
  }

f44 <- function(u,TPa,TPb,TPr,TPs)
  {
    TPb^2+2*(TPr-TPa)*(sqrt(u-TPs^2)-sqrt(u-(TPs+TPb)^2))-(TPs+TPb)*(TPr-TPa)*(asin(1-2*TPs^2/u)-asin(1-2*(TPs+TPb)^2/u))
  }

f4 <- function(uVect,TPa,TPb,TPr,TPs)
  {
    lim <- c((TPr-TPa)^2+(TPs)^2,(TPr)^2+(TPs)^2,(TPr-TPa)^2+(TPs+TPb)^2,(TPr)^2+(TPs+TPb)^2)
    ret <- NULL
    if(lim[2]<=lim[3])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[2]]
        ret <- c(ret,f41(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[3]]
        ret <- c(ret,f42(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[4]]
        ret <- c(ret,f43(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }

    if(lim[3] < lim[2])
      {
        uTEMP <- uVect[uVect < lim[1]]
        ret <- c(ret,rep(0,length(uTEMP)))
        uTEMP <- uVect[uVect >= lim[1] & uVect < lim[3]]
        ret <- c(ret,f41(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[3] & uVect < lim[2]]
        ret <- c(ret,f44(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >= lim[2] & uVect < lim[4]]
        ret <- c(ret,f43(uTEMP,TPa,TPb,TPr,TPs))
        uTEMP <- uVect[uVect >=lim[4]]
        ret <- c(ret,rep(0,length(uTEMP)))
      }
    
    ret/(4*TPa^2*TPb^2)
  }

fcase3 <- function(u,TPa,TPb,TPr,TPs) f1(u,TPa,TPb,TPr,TPs)+f2(u,TPa,TPb,TPr,TPs)+f3(u,TPa,TPb,TPr,TPs)+f4(u,TPa,TPb,TPr,TPs)
