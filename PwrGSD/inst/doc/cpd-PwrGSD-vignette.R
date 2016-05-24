### R code from vignette source 'cpd-PwrGSD-vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PwrGSDpd
###################################################
options(keep.source = TRUE, width = 70)
PwrGSDpd <- packageDescription("PwrGSD")


###################################################
### code chunk number 2: outcome-space
###################################################
  tlook <- c(7.14, 8.14, 9.14, 10.14, 10.64, 11.15, 12.14, 13.14,
             14.14, 15.14, 16.14, 17.14, 18.14, 19.14, 20.14)
  t0 <- 0:19      
  h0 <- c(rep(3.73e-04, 2), rep(7.45e-04, 3), rep(1.49e-03, 15))

  rhaz <-c(1, 0.9125, 0.8688, 0.7814, 0.6941, 0.6943, 0.6072, 0.5202,  
           0.4332, 0.652, 0.6524, 0.6527, 0.653, 0.6534, 0.6537, 
           0.6541, 0.6544, 0.6547, 0.6551, 0.6554)
  hc <- c(rep(1.05e-02, 2), rep(2.09e-02, 3), rep(4.19e-02, 15))
  hd1B <- c(0.1109, 0.1381, 0.1485, 0.1637, 0.2446, 0.2497, 0)


###################################################
### code chunk number 3: test-example
###################################################
  library(PwrGSD)
  test.example <- 
    PwrGSD(EfficacyBoundary=LanDemets(alpha=0.05, spending= ObrienFleming),
           FutilityBoundary=LanDemets(alpha=0.1,spending=ObrienFleming),
           RR.Futility = 0.82, sided="1<",method="A",accru =7.73, accrat=9818.65,
           tlook =tlook, tcut0 =t0, h0=h0, tcut1=t0, rhaz=rhaz, 
           tcutc0=t0, hc0=hc, tcutc1=t0, hc1=hc, 
           tcutd0B =c(0, 13), hd0B =c(0.04777, 0),
           tcutd1B =0:6, hd1B =hd1B,
           noncompliance =crossover, gradual =TRUE,
           WtFun =c("FH", "SFH", "Ramp"),
           ppar =c(0, 1, 0, 1, 10, 10))


###################################################
### code chunk number 4: max-effect
###################################################
  max.effect <- 0.80 + 0.05*(0:8)
  n.me <- length(max.effect)


###################################################
### code chunk number 5: censamt-settings
###################################################
  cens.amt <- 0.75 + 0.25*(0:2)
  n.ca <- length(cens.amt)


###################################################
### code chunk number 6: boundary-settings
###################################################
  Eff.bound.choice <- 1:2
  ebc.nms <- c("LanDemets(alpha=0.05, spending=ObrienFleming)",
               "LanDemets(alpha=0.05, spending=Pow(1))")
  n.ec <- length(Eff.bound.choice)


###################################################
### code chunk number 7: create-descr
###################################################
  descr <- as.data.frame(
              cbind(Eff.bound.choice=rep(Eff.bound.choice, each=n.ca*n.me),
                    cens.amt=rep(rep(cens.amt, each=n.me), n.ec),
                    max.effect=rep(max.effect, n.ec*n.ca)))

  descr$Eff.bound.choice <- ebc.nms[descr[["Eff.bound.choice"]]]


###################################################
### code chunk number 8: create-cpdPwrGSD-obj
###################################################
  test.example.set <- cpd.PwrGSD(descr)


###################################################
### code chunk number 9: loop-over-elements
###################################################
n.descr <- nrow(descr)

for(k in 1:n.descr){

  test.example.set$Elements[[k]]$call <- test.example$call

  test.example.set$Elements[[k]]$call$EfficacyBoundary <-
    parse(text=as.character(descr[k,"Eff.bound.choice"]))[[1]]

  test.example.set$Elements[[k]]$call$rhaz <-
    exp(descr[k,"max.effect"] * log(rhaz))
 
  test.example.set$Elements[[k]]$call$hc0 <- descr[k, "cens.amt"] * hc
  test.example.set$Elements[[k]]$call$hc1 <- descr[k, "cens.amt"] * hc

  test.example.set$Elements[[k]] <- update(test.example.set$Elements[[k]])
}


###################################################
### code chunk number 10: elements
###################################################
  test.example.subset <- 
      Elements(test.example.set, 
               subset=(substring(Eff.bound.choice, 32, 34)=="Obr" &   
                                 max.effect >= 1))


###################################################
### code chunk number 11: fig1LanDemets
###################################################
  plot(test.example.set, formula = ~ max.effect | stat * cens.amt,
       subset=(substring(Eff.bound.choice, 32, 34)=="Obr"))


###################################################
### code chunk number 12: fig2StochCurt
###################################################
  plot(test.example.set, formula = ~ max.effect | stat * cens.amt,
       subset=(substring(Eff.bound.choice, 32, 34)=="Pow"))


###################################################
### code chunk number 13: fig1-OBrien
###################################################
  plot(test.example.set, formula = ~ max.effect | stat * cens.amt,
       subset=(substring(Eff.bound.choice, 32, 34)=="Obr"))


###################################################
### code chunk number 14: fig2-Pow
###################################################
  plot(test.example.set, formula = ~ max.effect | stat * cens.amt,
       subset=(substring(Eff.bound.choice, 32, 34)=="Pow"))


