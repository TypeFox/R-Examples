## ----setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
library("mbbefd")
library("knitcitations")
cleanbib()
options("citation_format" = "pandoc")
bib <- read.bibtex("mbbefd.bib")
library(lattice)
my.settings <- canonical.theme(color=FALSE)
my.settings[['fontsize']] = list(text = 8, points = 4)
my.settings[['strip.background']]$col <- "darkgrey"
my.settings[['strip.border']]$col<- "black"  

## ---- fig.height=2.5, fig.width=6, echo=FALSE, warning=FALSE-------------
library(lattice)
n <- 100
Loss <- seq(1, 150, length=n)
mu <- 4.13
sigma <- 0.29
CDF <- plnorm(Loss, mu, sigma)
Density <- dlnorm(Loss, mu, sigma)
Survival <- 1 - CDF
dat <- data.frame(Loss=rep(Loss, 3),
                  Value=c(Density, CDF, Survival),
                  Type=gl(3, n,
                         labels=c("PDF: f(x)", 
                                  "CDF: F(x)", 
                                  "SF: S(x)=1-F(x)"),
                         ordered=TRUE))
xyplot(Value ~ Loss | Type, data=dat, ylab="", 
       layout=c(3,1), as.table=TRUE, t="l",
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2), 
       scales=list(relation="free", alternating=1))

## ----LEV, fig.height=2.5, fig.width=5, echo=FALSE------------------------
alpha <- 100
xyplot(Value ~ Loss | Type, 
       data=subset(dat, !Type %in% "PDF: f(x)"), 
       ylab="", 
       panel=function(x,y,...){
         x1 <- c(x[x<=alpha])
         if(panel.number()<2){
           panel.polygon(c(x1, rev(x1)),
                          rev(c(rep(1, length(x[x<=alpha])),
                           rev(y[x<=alpha]))),
                         col="skyblue", border=NA)
         }else{
           panel.polygon(c(x1, rev(x1)),
                         c(rep(0, length(x[x<=alpha])),
                           rev(y[x<=alpha])),
                        col="skyblue", border=NA)
         }
         panel.xyplot(x,y,...)
         #panel.text(x=20, y=0.4, label=paste("LEV[X]", cex=2)
       },
       layout=c(2,1), as.table=TRUE, t="l",
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2), 
       scales=list(relation="free", alternating=1))

## ----LEV2, fig.height=2.5, fig.width=5, echo=FALSE-----------------------
alpha1 <- 80
alpha2 <- 100
xyplot(Value ~ Loss | Type, 
       data=subset(dat, !Type %in% "PDF: f(x)"), 
       ylab="", 
       panel=function(x,y,...){
         x1 <- c(x[x >= alpha1 & x<=alpha2])
         if(panel.number()<2){
           panel.polygon(
             c(x1, rev(x1)),
             rev(c(rep(1, length(x[x >= alpha1 & x<=alpha2])),
               rev(y[x >= alpha1 & x<=alpha2]))),
              col="skyblue", border=NA)
         }else{
           panel.polygon(
             c(x1, rev(x1)),
             c(rep(0, length(x[x >= alpha1 & x<=alpha2])),
               rev(y[x >= alpha1 & x<=alpha2])),
             col="skyblue", border=NA)
         }
         panel.xyplot(x,y,...)
       },
       layout=c(2,1), as.table=TRUE, t="l",
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2), 
       scales=list(relation="free", alternating=1))

## ------------------------------------------------------------------------
S <- function(x){ 1 - plnorm(x, mu, sigma) }
(lyr <- integrate(S, 0, 100)$value - integrate(S, 0, 80)$value)

## ------------------------------------------------------------------------
(ILF <- integrate(S, 0, 100)$value / integrate(S, 0, 80)$value)

## ----ExposureCurve, dev.args=list(pointsize=8), fig.height=3, fig.width=3----
MPL <- 200
ExpectedLoss <- 65
Deductible <- seq(0, MPL, 1)
G <- sapply(Deductible, function(x){
  LEVx <- integrate(S, 0, x)$value
  LEVx/ExpectedLoss
})
plot(Deductible/MPL, G, 
     t="l", bty="n", lwd=1.5,
     main="Exposure Curve",
     xlab="Deductible as % of MPL",
     ylab="% of expected loss paid by insured",
     xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, lty=2)

## ----simufit, echo=TRUE, fig.height=5, fig.width=5, warning=FALSE--------
library(mbbefd)
set.seed(123456)
x <- c(rbeta(50, 3, 1/2), rmbbefd(50, 1/2, 1/10))
f1 <- fitDR(x, "mbbefd", method="mle")
summary(f1)
b1 <- bootDR(f1, niter=20)
summary(b1)
par(mfrow=c(1, 2))
denscomp(f1, demp=TRUE)
plot(b1, enhance=TRUE)

## ----simufit2, echo=TRUE, fig.height=5, fig.width=10, warning=FALSE------
f2 <- fitDR(x, "oibeta", method="mle")
f3 <- fitDR(x, "oiunif", method="mle")
gofstat(list(f1, f2, f3))
par(mfrow=c(1, 2))
cdfcomp(list(f1, f2, f3), leg=c("mbbefd", "oibeta", "oiunif"))
denscomp(list(f1, f2, f3), leg=c("mbbefd", "oibeta", "oiunif"), ylim=c(0,4), xleg="topleft")

## ----simufit3, echo=TRUE, fig.height=5, fig.width=5, warning=FALSE-------
eccomp(list(f1, f2, f3), leg=c("mbbefd", "oibeta", "oiunif"), do.points=FALSE)

## ----SwissReExample, message=FALSE, echo=FALSE, warning=FALSE------------
Client <- scan(textConnection(
'150 75 250 200 400 325 600 500 800 700 1000 900 1250 1125 1500 1375 1750 1625 2000 1875 2500 2250 3000 2750 4000 3500 5500 4750 9000 7250 12500 10750 18000 15250 24000 21000 36000 30000 48000 42000 72000 60000 90000 81000'))
MaxMPL <- Client[seq(from=1, to=length(Client),by = 2)]
MeanMPLGrossLoss <- Client[seq(from=2, to=length(Client),by = 2)]
GrossPremium <- scan(textConnection(
  '33434 14568 6324 4584 3341 1405 1169 683 613 554 700 552 1194 1490 4177 3527 3249 2712 2588 1988 657 1918'
))
ExposureCurve <- c(rep(1.5,3), rep(2.0,3), rep(3.0, 4), rep(4.0, 12))
ClientData <- data.frame(MaxMPL, MeanMPLGrossLoss, 
                         GrossPremium, ExposureCurve)
ClientData2 <- ClientData
names(ClientData2) <- c("Max MPL '000", "Mean MPL Gross Loss '000",
                        "Gross Premium '000", "Exposure Curve Parameter c")
library(pander)
panderOptions('big.mark', ',')
panderOptions('table.split.cells', 8)
pander(ClientData2, justify = rep('right', ncol(ClientData2)))
## Example layer
D <- 1.5e3*457/550
L <- 5e3*457/550
MPLoss <- 3.5e6
retainedDed <- D/MPLoss*1000
retainedLoss <- ecMBBEFD(retainedDed, b=swissRe(4)['b'], g=swissRe(4)['g'])
prem <- 1194 * (1 - retainedLoss)

## ----SwissRe2------------------------------------------------------------
XL_ROL <- function(Deductible, Limit, MaxMPL, 
                   MeanMPL, GrossPremium, C, ULR){
  DedPerMPL <- Deductible / ifelse(MaxMPL < Deductible, Deductible, 
                                   ifelse(MaxMPL<Limit, MaxMPL, Limit))
  LossShare <- 1 - apply(cbind(DedPerMPL, C), 1, 
                         function(x){ 
                           ecMBBEFD(x[1], 
                                          b=swissRe(x[2])['b'], 
                                          g=swissRe(x[2])['g'])
                           })
  NetPremium <- GrossPremium * ifelse(MaxMPL < Limit, 1, 
                                      Limit/MaxMPL)
  XL_Premium <- NetPremium * LossShare
  # Rate on Line
  sum(XL_Premium)/sum(NetPremium)*ULR
}

rol <- XL_ROL(Deductible = D,
              Limit = L,
              MaxMPL = ClientData$MaxMPL,
              MeanMPL = ClientData$MeanMPLGrossLoss,
              GrossPremium = ClientData$GrossPremium,
              C=ClientData$ExposureCurve, ULR=0.55)

