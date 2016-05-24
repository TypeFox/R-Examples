## distribution.name
## mean0
## mean1
## xbar
## sd
## df
## n
## alpha.right
## alpha.left
##
## stderr
##
## zc.right
## xbarc.right
## zc.left
## xbarc.left
## sided
## xbar.left
## xbar.right
## xbar.otherside
## pvalue.left
## pvalue.right
## pvalue
## power.left
## power.right
## power.total


NTplotTable <- function(
distribution.name=distribution.name,
type=             type,
mean0=            mean0,
mean1=            mean1,
xbar=             xbar,
## sd=               sd,
df=               df,
n=                n,
alpha.right=      alpha.right,
alpha.left=       alpha.left,
stderr=           stderr,
sigma.p1=         sigma.p1,
zc.right=         zc.right,
xbarc.right=      xbarc.right,
zc.left=          zc.left,
xbarc.left=       xbarc.left,
sided=            sided,
xbar.left=        xbar.left,
xbar.right=       xbar.right,
xbar.otherside=   xbar.otherside,
pvalue.left=      pvalue.left,
pvalue.right=     pvalue.right,
pvalue=           pvalue,
power.left=       power.left,
power.right=      power.right,
power.total=      power.total,
conf.left=        conf.left,
conf.right=       conf.right,
conf=             conf,
mean0.alist=      mean0.alist,
mean1.alist=      mean1.alist
) {

normalTable <- as.table(
  matrix(NA, 9, 12,
         dimnames=list(
           c(
             "xscale","xbarscale","zscale","z1scale",
             "alpha","beta","power","p","conf"),
           c("mean0","mean1","xbar","xbar.otherside","xbar.left","xbar.right",
             "xbarc.left","xbarc.right","Prob","sigma","n","df"))))



   normalTable[ "xscale" , "mean0"          ] <- mean0
   normalTable[ "xscale" , "mean1"          ] <- mean1
   normalTable[ "xscale" , "xbar"           ] <- xbar
## normalTable[ "xscale" , "xbar.otherside" ] <-
## normalTable[ "xscale" , "xbar.left"      ] <-
## normalTable[ "xscale" , "xbar.right"     ] <-
## normalTable[ "xscale" , "xbarc.left"     ] <-
## normalTable[ "xscale" , "xbarc.right"    ] <-
## normalTable[ "xscale" , "Prob"           ] <-
##   normalTable[ "xscale" , "sigma"          ] <- sd
   normalTable[ "xscale" , "sigma"          ] <- stderr*sqrt(n)
   normalTable[ "xscale" , "n"              ] <- n
   normalTable[ "xscale" , "df"             ] <- df

   normalTable[ "xbarscale" , "mean0"          ] <- mean0
   normalTable[ "xbarscale" , "mean1"          ] <- mean1
   normalTable[ "xbarscale" , "xbar"           ] <- xbar
   normalTable[ "xbarscale" , "xbar.otherside" ] <- xbar.otherside
   normalTable[ "xbarscale" , "xbar.left"      ] <- xbar.left
   normalTable[ "xbarscale" , "xbar.right"     ] <- xbar.right
   normalTable[ "xbarscale" , "xbarc.left"     ] <- xbarc.left
   normalTable[ "xbarscale" , "xbarc.right"    ] <- xbarc.right
## normalTable[ "xbarscale" , "Prob"           ] <-
   normalTable[ "xbarscale" , "sigma"          ] <- stderr
   normalTable[ "xbarscale" , "n"              ] <- 1
   normalTable[ "xbarscale" , "df"             ] <- df

   normalTable[ "zscale" , "mean0"          ] <- 0
   normalTable[ "zscale" , "mean1"          ] <- (mean1 - mean0)/stderr
   normalTable[ "zscale" , "xbar"           ] <- (xbar - mean0)/stderr
   normalTable[ "zscale" , "xbar.otherside" ] <- (xbar.otherside - mean0)/stderr
   normalTable[ "zscale" , "xbar.left"      ] <- (xbar.left - mean0)/stderr
   normalTable[ "zscale" , "xbar.right"     ] <- (xbar.right - mean0)/stderr
   normalTable[ "zscale" , "xbarc.left"     ] <- (xbarc.left - mean0)/stderr
   normalTable[ "zscale" , "xbarc.right"    ] <- (xbarc.right - mean0)/stderr
## normalTable[ "zscale" , "Prob"           ] <-
   normalTable[ "zscale" , "sigma"          ] <- stderr
   normalTable[ "zscale" , "n"              ] <- 1
   normalTable[ "zscale" , "df"             ] <- df

   normalTable[ "z1scale" , "mean0"          ] <- (mean0 - mean1)/sigma.p1
   normalTable[ "z1scale" , "mean1"          ] <- (mean1 - mean1)/sigma.p1
   normalTable[ "z1scale" , "xbar"           ] <- (xbar - mean1)/sigma.p1
   normalTable[ "z1scale" , "xbar.otherside" ] <- (xbar.otherside - mean1)/sigma.p1
   normalTable[ "z1scale" , "xbar.left"      ] <- (xbar.left - mean1)/sigma.p1
   normalTable[ "z1scale" , "xbar.right"     ] <- (xbar.right - mean1)/sigma.p1
   normalTable[ "z1scale" , "xbarc.left"     ] <- (xbarc.left - mean1)/sigma.p1
   normalTable[ "z1scale" , "xbarc.right"    ] <- (xbarc.right - mean1)/sigma.p1
## normalTable[ "z1scale" , "Prob"           ] <-
   normalTable[ "z1scale" , "sigma"          ] <- sigma.p1
   normalTable[ "z1scale" , "n"              ] <- 1
   normalTable[ "z1scale" , "df"             ] <- df

## normalTable[ "alpha" , "mean0"          ] <-
## normalTable[ "alpha" , "mean1"          ] <-
## normalTable[ "alpha" , "xbar"           ] <-
## normalTable[ "alpha" , "xbar.otherside" ] <-
## normalTable[ "alpha" , "xbar.left"      ] <-
## normalTable[ "alpha" , "xbar.right"     ] <-
   normalTable[ "alpha" , "xbarc.left"     ] <- if (sided != "right") alpha.left else NA
   normalTable[ "alpha" , "xbarc.right"    ] <- if (sided != "left") alpha.right else NA
   normalTable[ "alpha" , "Prob"           ] <- if (type=="hypothesis") sum(alpha.left, alpha.right, na.rm=TRUE) else NA
## normalTable[ "alpha" , "sigma"          ] <-
## normalTable[ "alpha" , "n"              ] <-
## normalTable[ "alpha" , "df"             ] <-

## normalTable[ "power" , "mean0"          ] <-
## normalTable[ "power" , "mean1"          ] <-
## normalTable[ "power" , "xbar"           ] <-
## normalTable[ "power" , "xbar.otherside" ] <-
## normalTable[ "power" , "xbar.left"      ] <-
## normalTable[ "power" , "xbar.right"     ] <-
   normalTable[ "power" , "xbarc.left"     ] <- if (!is.na(mean1)) power.left else NA
   normalTable[ "power" , "xbarc.right"    ] <- if (!is.na(mean1)) power.right else NA
   normalTable[ "power" , "Prob"           ] <- if (!is.na(mean1)) sum(power.left, power.right, na.rm=TRUE) else NA
## normalTable[ "power" , "sigma"          ] <-
## normalTable[ "power" , "n"              ] <-
## normalTable[ "power" , "df"             ] <-

## normalTable[ "beta" , "mean0"          ] <-
## normalTable[ "beta" , "mean1"          ] <-
## normalTable[ "beta" , "xbar"           ] <-
## normalTable[ "beta" , "xbar.otherside" ] <-
## normalTable[ "beta" , "xbar.left"      ] <-
## normalTable[ "beta" , "xbar.right"     ] <-
   normalTable[ "beta" , "xbarc.left"     ] <- if (!is.na(mean1)) 1 - power.left else NA
   normalTable[ "beta" , "xbarc.right"    ] <- if (!is.na(mean1)) 1 - power.right else NA
   normalTable[ "beta" , "Prob"           ] <- if (!is.na(mean1)) 1 - sum(power.left, power.right, na.rm=TRUE) else NA
## normalTable[ "beta" , "sigma"          ] <-
## normalTable[ "beta" , "n"              ] <-
## normalTable[ "beta" , "df"             ] <-

## normalTable[ "conf" , "mean0"          ] <-
## normalTable[ "conf" , "mean1"          ] <-
## normalTable[ "conf" , "xbar"           ] <-
## normalTable[ "conf" , "xbar.otherside" ] <-
## normalTable[ "conf" , "xbar.left"      ] <-
## normalTable[ "conf" , "xbar.right"     ] <-
   normalTable[ "conf" , "xbarc.left"     ] <- if (type=="confidence") conf.left else NA
   normalTable[ "conf" , "xbarc.right"    ] <- if (type=="confidence") conf.right else NA
   normalTable[ "conf" , "Prob"           ] <- if (type=="confidence") conf else NA
## normalTable[ "conf" , "sigma"          ] <-
## normalTable[ "conf" , "n"              ] <-
## normalTable[ "conf" , "df"             ] <-

## normalTable[ "p" , "mean0"          ] <-
## normalTable[ "p" , "mean1"          ] <-
## normalTable[ "p" , "xbar"           ] <-
## normalTable[ "p" , "xbar.otherside" ] <-
   normalTable[ "p" , "xbar.left"      ] <- if (type=="hypothesis" && !is.na(xbar)) pvalue.left else NA
   normalTable[ "p" , "xbar.right"     ] <- if (type=="hypothesis" && !is.na(xbar)) pvalue.right else NA
## normalTable[ "p" , "xbarc.left"     ] <-
## normalTable[ "p" , "xbarc.right"    ] <-
   normalTable[ "p" , "Prob"           ] <- if (type=="hypothesis" && !is.na(xbar)) sum(pvalue.left, pvalue.right, na.rm=TRUE) else NA
## normalTable[ "p" , "sigma"          ] <-
## normalTable[ "p" , "n"              ] <-
## normalTable[ "p" , "df"             ] <-


## > dimnames(normalTable)
## [[1]]
## [1] "xscale"    "xbarscale" "zscale"    "z1scale"
## [5] "alpha"     "beta"      "power"     "p"
## [9] "conf"

## [[2]]
##  [1] "mean0"          "mean1"          "xbar"
##  [4] "xbar.otherside" "xbar.left"      "xbar.right"
##  [7] "xbarc.left"     "xbarc.right"    "Prob"
## [10] "sigma"          "n"              "df"

dimnames(normalTable) <- list(expression(x, bar(x), z, z[1], alpha, beta, power, p, "Confidence"),
                              c(as.expression(mean0.alist[[1]]), as.expression(mean1.alist[[1]]), ## mu[0], mu[1],
                              expression(w["obs"], ##==bar(x)["obs"],
                                         w["other"], ##==bar(x)["other"],
                                         w["left"], ##==bar(x)["left"],
                                         w["right"], ##==bar(x)["right"],
                                         w[crit.L], ##==bar(x)[crit.L],
                                         w[crit.R], ##==bar(x)[crit.R],
                                         "Probability", sigma, n, df)))
if (distribution.name=="t")
  dimnames(normalTable)[[1]][3:4] <- c(expression(t), expression(t[1]))

column.sequence <- if (sided=="both")
                     c("w[crit.L]", "Probability", "w[crit.R]",
                       'w["left"]', 'w["right"]')
                   else
                     "Probability"
prob <- normalTable[c("p","alpha","power","beta","Confidence"),
                    column.sequence,
                    drop=FALSE]
if (type == "confidence") {a
  dimnames(normalTable)[[2]][7:8] <- expression(w[LCL], w[UCL]) ##mu[LCL], mu[UCL])
}
## recover()
prob <- prob[!is.na(prob[,"Probability"]), , drop=FALSE]
if (ncol(prob) >= 3) dimnames(prob)[[2]][1:3] <- c("Left", "Combined", "Right")
if (type == "hypothesis" && sided == "both") {
  if (!is.na(mean1)) prob["beta", c("Left","Right")] <- NA
  if ("p" %in% dimnames(prob)[[1]])
    prob["p", c("Left","Right")] <- prob["p", c('w["left"]', 'w["right"]')]
  prob <- prob[,1:3, drop=FALSE]
}
if (type == "confidence" && sided == "both") {

  prob <- prob["Confidence",1:3, drop=FALSE]
  dimnames(prob)[[2]][2] <- "Confidence"
  dimnames(prob)[[1]] <- "Probability"
}

scales <- normalTable[c(2, 3, if(!is.na(mean1)) 4),
                      c(if (type == "hypothesis") 1,
                        if (!is.na(mean1) && type == "hypothesis") 2,
                        if (!is.na(xbar))  3,
                        if (!is.na(xbar) && sided == "both" && type == "hypothesis") 4,
                        ## if (!is.na(xbar) && sided == "both" && type == "hypothesis") 5,
                        ## if (!is.na(xbar) && sided == "both" && type == "hypothesis") 6,
                        if (sided != "right") 7,
                        if (sided != "left")  8)]



## if (FALSE) {
## prob <- t(normalTable[c("alpha","beta","power","p","conf"), "Prob", drop=FALSE])
## row.names(prob) <- "Probability"
## prob <- prob[, !is.na(prob), drop=FALSE]

## scales <- normalTable[c("xbarscale","zscale", if(!is.na(mean1)) "z1scale"),
##                       c(if (type == "hypothesis") "mean0",
##                         if (!is.na(mean1) && type == "hypothesis") "mean1",
##                         if (!is.na(xbar))  "xbar",
##                         if (!is.na(xbar) && sided == "both" && type == "hypothesis")  "xbar.otherside",
##                         ## if (!is.na(xbar) && sided == "both" && type == "hypothesis") "xbar.left",
##                         ## if (!is.na(xbar) && sided == "both" && type == "hypothesis") "xbar.right",
##                         if (sided != "right") "xbarc.left",
##                         if (sided != "left")  "xbarc.right")]

## if (distribution.name=="t")
##   dimnames(scales)[[1]] <- c(expression(bar(x)), expression(t), expression(t[1]))
## else
##   dimnames(scales)[[1]] <- c(expression(bar(x)), expression(z), expression(z[1]))

## if (type == "hypothesis") {
##   dimnames(scales)[[2]] <- c(expression(mu[0]), expression(mu[1]),
##                                 expression(bar(x)), expression(bar(x)[other]),
##                                 expression(bar(x)[crit.L]), expression(bar(x)[crit.R]))
## }

## if (type == "confidence") {
##   dimnames(normalTable)[[2]][ dimnames(normalTable)[[2]] %in% c("xbarc.left", "xbarc.right") ] <- c("LCL", "UCL")
##   dimnames(scales)[[2]][ dimnames(scales)[[2]] %in% c("xbar", "xbarc.left", "xbarc.right")] <-
##     c(expression(bar(x), mu[LCL], mu[UCL]))
##   dimnames(pro2b)[[2]] <- "Confidence"
## }
## }

list(normalTable=normalTable, prob=prob, scales=scales)
}
