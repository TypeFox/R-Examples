### R code from vignette source 'gRain-intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: gRain-intro.Rnw:25-29
###################################################
require( gRain )
prettyVersion <- packageDescription("gRain")$Version
prettyDate <- format(Sys.Date())
dir.create( "figures" )


###################################################
### code chunk number 2: gRain-intro.Rnw:107-120
###################################################
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)
s    <- cptable(~smoke, values=c(5,5), levels=yn)
l.s  <- cptable(~lung|smoke, values=c(1,9,1,99), levels=yn)
b.s  <- cptable(~bronc|smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
x.e  <- cptable(~xray|either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)
plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
plist
net1 <- grain(plist)
net1


###################################################
### code chunk number 3: LS
###################################################
require(Rgraphviz)
plot(net1)


###################################################
### code chunk number 4: gRain-intro.Rnw:200-202
###################################################
library(gRain)
options("prompt"="> ","width"=85)


###################################################
### code chunk number 5: gRain-intro.Rnw:244-253
###################################################
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)
s    <- cptable(~smoke, values=c(5,5), levels=yn)
l.s  <- cptable(~lung|smoke, values=c(1,9,1,99), levels=yn)
b.s  <- cptable(~bronc|smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
x.e  <- cptable(~xray|either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)


###################################################
### code chunk number 6: gRain-intro.Rnw:259-265
###################################################
plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
plist
plist$tub
plist$either ## Notice: a logical node
net1 <- grain(plist)
net1


###################################################
### code chunk number 7: gRain-intro.Rnw:279-280
###################################################
querygrain(net1, nodes=c("lung","bronc"), type="marginal")


###################################################
### code chunk number 8: gRain-intro.Rnw:286-287
###################################################
querygrain(net1,nodes=c("lung","bronc"), type="joint")


###################################################
### code chunk number 9: gRain-intro.Rnw:292-295
###################################################
net12  <- setEvidence(net1, evidence=list(asia="yes", dysp="yes"))
net12  <- setEvidence(net1,
                      nodes=c("asia", "dysp"), states=c("yes", "yes"))


###################################################
### code chunk number 10: gRain-intro.Rnw:300-301
###################################################
pEvidence( net12 )


###################################################
### code chunk number 11: gRain-intro.Rnw:308-310
###################################################
querygrain( net12, nodes=c("lung","bronc") )
querygrain( net12, nodes=c("lung","bronc"), type="joint" )


###################################################
### code chunk number 12: gRain-intro.Rnw:320-322
###################################################
net13 <- setEvidence(net1,nodes=c("either", "tub"),
                     states=c("no","yes"))


###################################################
### code chunk number 13: gRain-intro.Rnw:327-328
###################################################
pEvidence( net13 )


###################################################
### code chunk number 14: gRain-intro.Rnw:334-335
###################################################
querygrain( net13, nodes=c("lung","bronc"), type="joint" )


###################################################
### code chunk number 15: gRain-intro.Rnw:342-344
###################################################
tt <- querygrain( net1, type="joint")
sum(tt==0)/length(tt)


###################################################
### code chunk number 16: gRain-intro.Rnw:349-350
###################################################
sum(tableSlice(tt, c("either","tub"), c("no","yes")))


###################################################
### code chunk number 17: gRain-intro.Rnw:357-365
###################################################
yn <- c("yes","no")
eps <- 1e-100
a    <- cptable(~a,   values=c(1,eps),levels=yn)
b.a  <- cptable(~b+a, values=c(1,eps,eps,1),levels=yn)
c.b  <- cptable(~c+b, values=c(1,eps,eps,1),levels=yn)
plist <- compileCPT(list(a, b.a, c.b))
bn   <- grain(plist)
( tt   <- querygrain(bn, type="joint") )


###################################################
### code chunk number 18: gRain-intro.Rnw:369-370
###################################################
querygrain(setEvidence(bn, nodes=c("a","c"), state=c("no", "yes")))


###################################################
### code chunk number 19: gRain-intro.Rnw:376-384
###################################################
eps <- 1e-200
a    <- cptable(~a,   values=c(1,eps),levels=yn)
b.a  <- cptable(~b+a, values=c(1,eps,eps,1),levels=yn)
c.b  <- cptable(~c+b, values=c(1,eps,eps,1),levels=yn)
plist <- compileCPT(list(a, b.a, c.b))
bn   <- grain(plist)
( tt   <- querygrain(bn, type="joint") )
querygrain(setEvidence(bn, nodes=c("a","c"), state=c("no", "yes")))


###################################################
### code chunk number 20: gRain-intro.Rnw:446-455
###################################################
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)

( plist1 <- compileCPT( list( a, t.a ) ) )
plist1[[1]]
plist1[[2]]
( chest1 <- grain(plist1) )
querygrain( chest1 )


###################################################
### code chunk number 21: gRain-intro.Rnw:469-472
###################################################
setFinding(  chest1, nodes="asia", states="yes")
setEvidence( chest1, nodes="asia", states="yes")
setEvidence( chest1, evidence=list(asia="yes"))


###################################################
### code chunk number 22: gRain-intro.Rnw:476-477
###################################################
querygrain( setEvidence( chest1, evidence=list(asia="yes")) )


###################################################
### code chunk number 23: gRain-intro.Rnw:499-501
###################################################
g.a <- parray(c("guess.asia", "asia"), levels=list(yn, yn),
              values=c(.8,.2, .1,.9))


###################################################
### code chunk number 24: gRain-intro.Rnw:510-513
###################################################
( plist2 <- compileCPT( list( a, t.a, g.a ) ) )
( chest2 <- grain(plist2) )
querygrain( chest2 )


###################################################
### code chunk number 25: gRain-intro.Rnw:521-522
###################################################
querygrain( setEvidence( chest2, evidence=list(guess.asia="yes")) )


###################################################
### code chunk number 26: gRain-intro.Rnw:532-533
###################################################
querygrain( setEvidence( chest1, evidence=list(asia=c(.8, .1))) )


###################################################
### code chunk number 27: gRain-intro.Rnw:539-540
###################################################
querygrain( setEvidence( chest1, evidence=list(asia=c(1, 0))) )


###################################################
### code chunk number 28: gRain-intro.Rnw:579-582
###################################################
dG  <- dag(~A:B)
uG  <- ug(~A:B)
par(mfrow=c(1,2)); plot( dG ); plot( uG )


###################################################
### code chunk number 29: gRain-intro.Rnw:590-592
###################################################
dat <-as.table(parray(c("A","B"), levels=c(2,2), values=c(0,0,2,3)))
class( dat )


###################################################
### code chunk number 30: gRain-intro.Rnw:598-600
###################################################
gr.dG <- compile( grain( dG, dat ) )
gr.uG <- compile( grain( uG, dat ) )


###################################################
### code chunk number 31: gRain-intro.Rnw:614-616
###################################################
extractCPT( dat, dG )
c( extractPOT( dat, uG ) )


###################################################
### code chunk number 32: gRain-intro.Rnw:635-638
###################################################
p.A.g.B <- tableDiv(dat, tableMargin(dat, "B"))
p.B <- tableMargin(dat, "B")/sum(dat)
p.AB <- tableMult( p.A.g.B, p.B)


###################################################
### code chunk number 33: gRain-intro.Rnw:657-659
###################################################
e <- 1e-2
(dat.e <- dat + e)


###################################################
### code chunk number 34: gRain-intro.Rnw:663-666
###################################################
pe.A.g.B <- tableDiv(dat.e, tableMargin(dat, "B"))
pe.B <- tableMargin(dat.e, "B")/sum(dat.e)
pe.AB  <- tableMult( pe.A.g.B, pe.B )


###################################################
### code chunk number 35: gRain-intro.Rnw:672-673
###################################################
dat.e / sum( dat.e )


###################################################
### code chunk number 36: gRain-intro.Rnw:683-684
###################################################
gr.dG <- compile( grain( dG, dat, smooth=e ) )


###################################################
### code chunk number 37: gRain-intro.Rnw:689-690
###################################################
extractCPT( dat, dG, smooth=e)


###################################################
### code chunk number 38: gRain-intro.Rnw:695-697
###################################################
querygrain( gr.dG )
querygrain( gr.uG )


###################################################
### code chunk number 39: gRain-intro.Rnw:702-704
###################################################
querygrain(setFinding(gr.dG, nodes="B", states="B1"))
querygrain(setFinding(gr.uG, nodes="B", states="B1"))


###################################################
### code chunk number 40: gRain-intro.Rnw:712-713
###################################################
gr.uG <- compile( grain( uG, dat, smooth=e) )


###################################################
### code chunk number 41: gRain-intro.Rnw:717-718
###################################################
c( extractPOT( dat, uG, smooth=e ) )


###################################################
### code chunk number 42: gRain-intro.Rnw:724-726
###################################################
querygrain( gr.uG )
querygrain( gr.dG )


###################################################
### code chunk number 43: gRain-intro.Rnw:731-733
###################################################
querygrain( setFinding(gr.uG, nodes="B", states="B1") )
querygrain( setFinding(gr.dG, nodes="B", states="B1") )


