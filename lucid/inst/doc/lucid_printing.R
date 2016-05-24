## ----SETUP, echo=FALSE, results="hide"--------------------------------------------------
library("knitr")
opts_chunk$set(fig.align="center", fig.width=7, fig.height=7)
options(width=90)

## ----echo=FALSE-------------------------------------------------------------------------
df1 <- data.frame(effect=c(-13.5, 4.5,  24.5, 6.927792e-14, -1.75,
                    16.5, 113.5000))
rownames(df1) <- c("A","B","C","C1","C2","D","(Intercept)")
print(df1)

## ----message=FALSE----------------------------------------------------------------------
require("lucid")
options(digits=7) # knitr defaults to 4, R console uses 7
lucid(df1)

## ---------------------------------------------------------------------------------------
print(antibiotic)

## ---------------------------------------------------------------------------------------
lucid(antibiotic)

## ----echo=FALSE, message=FALSE, fig.height=4, fig.width=6-------------------------------
require(lattice)
anti=antibiotic # make a copy of the data to reverse the levels
anti$bacteria <- factor(anti$bacteria, levels=rev(anti$bacteria))

cust <- myyscale.component <- function(...) {  #Custom y-scale function component
  ans <- yscale.components.default(...)
  ans$right <- ans$left
  foo <- ans$right$labels$at
  ans$right$labels$labels <- rev(c(" 870    ","   1    ","   0.001","   0.005"," 100    ",
                               "850    "," 800    ","   3    "," 850    ","   1    ",
                               " 10    ","   0.007", "  0.03 ","   1    ",
                               "  0.001","   0.005"))
  return(ans)
}

dotplot(bacteria~ -log10(penicillin), anti,
        cex=1, xlim=c(-4,4), #xlab="variance component (log10 scale)",
        scales=list(x=list(at= c(-2,0,2),
                      labels=c('100','1','.01')),
          y=list(relation="free", fontfamily='mono')), # 'free' required for 2nd axis
        yscale.components=cust,
        #this creates more space on the right hand side of the plot
        par.settings=list(layout.widths=list(left.padding=10,right.padding=10))
        )


## ----message=FALSE----------------------------------------------------------------------
require("nlme")
data(Rail)
mn <- lme(travel~1, random=~1|Rail, data=Rail)
vc(mn)

require("lme4")
m4 <- lmer(travel~1 + (1|Rail), data=Rail)
vc(m4)

# require("asreml")
# ma <- asreml(travel~1, random=~Rail, data=Rail)
# vc(ma)
##         effect component std.error z.ratio constr
##  Rail!Rail.var    615.3      392.6     1.6    pos
##     R!variance     16.17       6.6     2.4    pos

## ----message=FALSE----------------------------------------------------------------------
require("nlme")
data(Rail)
require("rjags")
m5 <-
"model {
for(i in 1:nobs){
  travel[i] ~ dnorm(mu + theta[Rail[i]], tau)
}
for(j in 1:6) {
  theta[j] ~ dnorm(0, tau.theta)
}
mu ~ dnorm(50, 0.0001) # Overall mean. dgamma() 
tau ~ dgamma(1, .001)
tau.theta ~ dgamma(1, .001)
residual <- 1/sqrt(tau)
sigma.rail <- 1/sqrt(tau.theta)
}"
jdat <- list(nobs=nrow(Rail), travel=Rail$travel, Rail=Rail$Rail)
jinit <- list(mu=50, tau=1, tau.theta=1)
tc5 <- textConnection(m5)
j5 <- jags.model(tc5, data=jdat, inits=jinit, n.chains=2, quiet=TRUE)
close(tc5)
c5 <- coda.samples(j5, c("mu","theta", "residual", "sigma.rail"), 
                   n.iter=100000, thin=5)

## ---------------------------------------------------------------------------------------
vc(c5)

## ---------------------------------------------------------------------------------------
m4
ranef(m4)

## ----echo=FALSE-------------------------------------------------------------------------
# Results are from lme4_1.1-7, as.data.frame(VarCorr(m2b))
d1 <- structure(list(grp = c("new.gen", "one", "one.1", "one.2", "one.3",
"one.4", "one.5", "one.6", "one.7", "one.8", "one.9", "one.10",
"one.11", "one.12", "one.13", "Residual"), var1 = c("(Intercept)",
"r1:c3", "r1:c2", "r1:c1", "c8", "c6", "c4", "c3", "c2", "c1",
"r10", "r8", "r4", "r2", "r1", NA), var2 = c(NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_, NA_character_),
    vcov = c(2869.44692139271, 5531.57239635089, 58225.767835444,
    128004.156092455, 6455.74953933247, 1399.72937329085, 1791.65071661348,
    2548.88470543732, 5941.79076230161, 0, 1132.95013713932,
    1355.22907294114, 2268.72957045473, 241.789424531994, 9199.9021721834,
    4412.1096176349), sdcor = c(53.5672187199663, 74.3745413185916,
    241.300161283502, 357.7766846686, 80.3476791160297, 37.4129572914365,
    42.3278952537623, 50.4864804223598, 77.0830121511972, 0,
    33.6593246684974, 36.8134360382339, 47.6311827530529, 15.5495795612613,
    95.9161205021523, 66.4237127661116)), .Names = c("grp", "var1",
    "var2", "vcov", "sdcor"), row.names = c(NA, -16L), class = "data.frame")

d2 <- structure(list(grp = c("new.gen", "one", "one.1", "one.2", "one.3",
"one.4", "one.5", "one.6", "one.7", "one.8", "one.9", "one.10",
"one.11", "one.12", "one.13", "Residual"), var1 = c("(Intercept)",
"r1:c3", "r1:c2", "r1:c1", "c8", "c6", "c4", "c3", "c2", "c1",
"r10", "r8", "r4", "r2", "r1", NA), var2 = c(NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_, NA_character_),
    vcov = c(3228.41890564251, 7688.13916836557, 69747.5508913552,
    107427.043198654, 6787.00354507896, 1636.12771714548, 12268.4603217744,
    2686.30159414561, 7644.72994565782, 0.00122514315152732,
    1975.50482871438, 1241.42852423718, 2811.24084391436, 928.227473340838,
    10363.5849610346, 4126.83169047631), sdcor = c(56.8191772700249,
    87.6820344675326, 264.097616216722, 327.760649252856, 82.3832722406616,
    40.4490756031022, 110.763081944186, 51.8295436420735, 87.4341463368736,
    0.0350020449620779, 44.4466514904597, 35.23391156595, 53.02113582256,
    30.4668257838069, 101.801694293536, 64.2404210017051)), .Names = c("grp",
"var1", "var2", "vcov", "sdcor"), row.names = c(NA, -16L), class = "data.frame")
out <- cbind(d1[, c(2,4,5)], sep="   ",d2[,4:5])
names(out) <- c('term','vcov-bo','sdcor-bo','sep','vcov-ne','sdcor-ne')

## ---------------------------------------------------------------------------------------
print(out)

## ---------------------------------------------------------------------------------------
lucid(out, dig=4)

## ---------------------------------------------------------------------------------------
noquote(lucid(as.matrix(head(mtcars)),2))

## ----finish, echo=FALSE, results="asis"-------------------------------------------------
toLatex(sessionInfo(), locale=FALSE)

