##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 2   ##
## model formulations       ##
##############################

opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 2  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 2.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################

########################################
# Fig 2.1. conceptual model scheme
########################################

par(mar=c(0,0,0,0))
par(mfrow=c(2,1))
names <- c("STATE 1","STATE 2","STATE 3")
M <- matrix(nrow=3,ncol=3,byrow=TRUE,data=c(
#   1 2 3
    0,1,0,     # 1
    0,0,0,     # 2
    1,1,0  ))  # 3
pp<-plotmat(M,pos=c(1,2),curve=0,name=names,lwd=1,cex.txt=0.0,
            box.lwd=2,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.5)

# extra arrows: flow 5 to Detritus and flow 2 to detritus
state1   <-pp$comp[names=="STATE 1",]
state2   <-pp$comp[names=="STATE 2",]
state3   <-pp$comp[names=="STATE 3",]

# forcing
m1 <- 0.5*(state1+state2)
m2 <- c(0.25,0.85)
segments (m1[1],m1[2],m2[1],m2[2],lwd=1,lty=2)
text(m2[1]-0.01,m2[2]+0.03,"forcing function",adj=c(0.5,0.5))

# output variable
m1 <- state1
m1[1] <-m1[1]+ pp$radii[1,1]
m2 <- m1
m2[1]<-m2[1]+0.25
segments (m1[1],m1[2],m2[1],m2[2],lwd=1)
textellipse(mid=m2,radx=pp$radii[1,1],rady=pp$radii[1,2],lwd=1,shadow.size=0,lab="output variable")

m1 <- state3
m1[1] <-m1[1]+ pp$radii[3,1]
m2<-m1
m2[1]<- m2[1]+0.1
straightarrow(from=m1,to=m2,lwd=1,arr.pos=1,arr.type="triangle")

m1 <- 0.5*(state2+state3)
m1[2]<-m1[2]-0.05
text(m1[1],m1[2],"flow,interaction",adj=c(0.5,0.5))
subtitle()

########################################
# Fig. 2.2 Mass balance
########################################

par(mfrow=c(2,1))
par(mar=c(0,0,0,0))
openplotmat()
names <- c("SV1","SV2")
pos <- coordinates(2)
pp <- straightarrow (pos[1,]+c(0,0.4), pos[1,],arr.pos=0.4)
text(pp[1]+0.035,pp[2],"F0")
pp<-bentarrow(pos[2,],pos[2,]+c(0.2,-0.2), arr.pos=1)
text(pp[1],pp[2]-0.1,"F2")

M <- matrix(nrow=2,ncol=2,byrow=TRUE,data=c(
#   1 2 
    0,0,    # 1
    1,0 ))  # 2

pp<-plotmat(M,pos=pos,curve=0,name=names,lwd=1,prefix="F",dtext=0.6,
            box.lwd=2,box.type="square",box.prop=0.5,arr.lwd=2,
            arr.pos=0.5,add=TRUE)
rect(0.1,0.3,0.9,0.7)
text(0.9,0.31,adj=c(1,0),"Model domain", font=3)

text(0.1,0.1,"External world", adj=c(0,0), font=3)
box(col="grey")
subtitle()

########################################
# Fig 2.3. NPZZDD model
########################################
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
names <- c("PHYTO",expression(NH[4]^"+"),"ZOO","DETRITUS","BotDET","FISH")
M <- matrix(nrow=6,ncol=6,byrow=TRUE,data=c(
#   p n z  d  b  f
    0,1,0, 0, 0, 0, #p
    0,0,4, 10,11,0, #n
    2,0,0, 0, 0, 0, #z
    8,0,13,0, 0, 12,#d
    9,0,0, 7, 0, 0, #b
    0,0,5, 0, 0, 0  #f
    ))
pp<-plotmat(M,pos=c(1,2,1,2),curve=0,name=names,lwd=1,my=0.,cex.txt=0.8,
            box.lwd=2,box.size=0.08,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.4,shadow.size=0.01,prefix="f")

# extra arrows: flow 5 to Detritus and flow 2 to detritus
phyto   <-pp$comp[1,]
zoo     <-pp$comp[3,]
nh3     <-pp$comp[2,]
detritus<-pp$comp[4,]
fish    <-pp$comp[6,]

# flow6->detritus
m2 <- 0.5*(zoo+fish)
m1 <- detritus
m1[1]<-m1[1]+ pp$radii[4,1]
mid<-straightarrow (from=m2,to=m1,arr.type="triangle",arr.pos=0.4,lwd=1)
text(mid[1],mid[2]+0.03,"f6",cex=0.8)

# flow3->detritus
m2 <- 0.5*(zoo+phyto)
m1 <- detritus
m1[1] <-m1[1]+ pp$radii[3,1]*0.2
m1[2]<-m1[2] + pp$radii[3,2]
mid<-straightarrow (from=m2,to=m1,arr.type="triangle",arr.pos=0.3,lwd=1)
text(mid[1]-0.01,mid[2]+0.03,"f3",cex=0.8)

# solar radiation
m1 <- 0.5*(nh3+phyto)
m2 <- c(0.25,0.9)
segments (m1[1],m1[2],m2[1],m2[2],lwd=1,lty=2)
text(m2[1]-0.01,m2[2]+0.03,"solar radiation",adj=c(0.5,0.5))

# chlorophyll
m1 <- phyto
m1[1] <-m1[1]+ pp$radii[1,1]
m2 <- m1
m2[1]<-m2[1]+0.25
segments (m1[1],m1[2],m2[1],m2[2],lwd=1)
textellipse(mid=m2,radx=pp$radii[1,1],rady=pp$radii[1,2],lwd=1,shadow.size=0,lab="Chlorophyll")
subtitle()

########################################
# Fig 2.4. PZZDD model (C-units)
########################################

par(mar=c(0,0,0,0))
names <- c("PHYTOC","ZOOC","DETRITUSC","BotDETC","FISHC")
M <- matrix(nrow=5,ncol=5,byrow=TRUE,data=c(
#   p z  d  b  f
    0,0, 0, 0, 0, #p
    2,0, 0, 0, 0, #z
    8,13,0, 0, 12, #d
    9,0, 7, 0, 0, #b
    0,5, 0, 0, 0  #f
    ))
pp<-plotmat(M,pos=c(1,3,1),curve=0,name=names,lwd=1,my=-0.05,cex.txt=0.8,
            box.lwd=2,box.size=0.08,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.4,shadow.size=0.01,prefix="f")

# extra arrows: flow 5 to Detritus and flow 2 to detritus
phyto   <-pp$comp[1,]
zoo     <-pp$comp[2,]
detritus<-pp$comp[3,]
botdet  <-pp$comp[4,]
fish    <-pp$comp[5,]

# flow6->detritus
m2 <- 0.5*(zoo+fish)
m1 <- detritus
m1[1]<-m1[1]- pp$radii[3,1]
m1[2]<-m1[2]- pp$radii[3,2]
mid<-straightarrow (from=m2,to=m1,arr.type="triangle",arr.pos=0.4,lwd=1)
text(mid[1],mid[2]+0.03,"f6",cex=0.8)

# flow3->detritus
m2 <- 0.5*(zoo+phyto)
m1 <- detritus
m1[1] <-m1[1]- pp$radii[3,1]
m1[2]<-m1[2] + pp$radii[3,2]
mid<-straightarrow (from=m2,to=m1,arr.type="triangle",arr.pos=0.3,lwd=1)
text(mid[1],mid[2]-0.03,"f3",cex=0.8)

# f1: production
m1 <- phyto+c(0,0.15)
m2 <- phyto+ c(0,pp$radii[1,2])
mid<-straightarrow (from=m1,to=m2,arr.type="triangle",arr.pos=0.5,lwd=1)
text(mid[1]+0.03,mid[2],"f1",cex=0.8)

# respiration flows
mid <- bentarrow(from=zoo-c(pp$radii[1,1],0), to=zoo-c(0.125,0.075),
                  arr.type="triangle", lwd=1,arr.pos=1)
text(mid[1]+0.02,mid[2],"f4",cex=0.8)

mid <- bentarrow(from=detritus+c(pp$radii[1,1],-0.03),
                 to=detritus+c(0.125,-0.075), lwd=1,arr.pos=1,arr.type="triangle")
text(mid[1]+0.02,mid[2],"f10",cex=0.8)

mid <- bentarrow(from=botdet+c(pp$radii[4,1],0), to=botdet+c(0.125,-0.075),
                 lwd=1,arr.pos=1,arr.type="triangle")
text(mid[1]-0.02,mid[2],"f11",cex=0.8)

mid <- bentarrow(from=fish+c(pp$radii[5,1],0), to=fish+c(0.125,-0.075),
                 lwd=1,arr.pos=1,arr.type="triangle")
text(mid[1]+0.02,mid[2],"f14",cex=0.8)

# chlorophyll
m1 <- phyto
m1[1] <-m1[1]+ pp$radii[1,1]
m2 <- m1
m2[1]<-m2[1]+0.25
segments (m1[1],m1[2],m2[1],m2[2],lwd=1)
textellipse(mid=m2,radx=pp$radii[1,1],rady=pp$radii[1,2],lwd=1,shadow.size=0,lab="Chlorophyll")
subtitle()

########################################
# Types of interactions
########################################

basicplot <- function(lab1,lab2,mx=0.0,sx = 0.2,sy = 0.1,bcol="white",...)
{
openplotmat()
elpos  <-coordinates (c(1,1),mx=mx)

straightarrow (to=elpos[2,],from=elpos[1,],arr.pos=0.75,... )
textrect(elpos[1,],sx,sy,lab=lab1,shadow.size=0.02)
textrect(elpos[2,],sx,sy,lab=lab2,shadow.size=0.02,box.col=bcol)   
box(col="grey")
return(elpos)
}

########################################
# Figure 2.5. order of chemical reactions
########################################

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))

A <- matrix(nrow=2,ncol=2,0);A[2,1]<-1
pp <-plotmat(A,curve=0,pos=c(1,1),box.type="rect",name=c("A","B"),cex.txt=0)
text(0.6,0.5,expression(k[1]*"[A]"),font=3)
box(col="grey")
title("reaction order 1")
writelabel("A",line=0,at=-0.05)

# order 2
openplotmat()
elpos<-coordinates (c(2,1))
tt<-treearrow(from=elpos[1:2,],to=elpos[3,],arr.side=2)
for ( i in 1:3) textrect (elpos[i,],0.1,0.1,lab=LETTERS[i],cex=1.5)
text(tt[1]+0.05,tt[2]+0.05,expression(k[2]*"[A]"*"[B]"),font=3,adj=0)
box(col="grey")
title("reaction order 2")
writelabel("B",line=0,at=-0.05)

# order 3
openplotmat()
elpos<-coordinates (c(3,1))
tt<-treearrow(from=elpos[1:3,],to=elpos[4,],arr.side=2,path="H")
lab <- c("A","A","B","C")
for ( i in 1:4) textrect (elpos[i,],0.1,0.1,lab=lab[i],cex=1.5)
text(tt[1]+0.05,tt[2]+0.05,expression(k[3]*"[A]"*"[A]"*"[B]"),font=3,adj=0)
box(col="grey")
title("reaction order 3")
writelabel("C",line=0,at=-0.05)
subtitle()

########################################
# Fig. 2.6 Chemical reactions
########################################

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
# reversible reaction
# reversible reaction
openplotmat()
elpos<-coordinates (c(3,1))
treearrow(from=elpos[1:3,],to=elpos[4,],arr.side=2,path="H")
treearrow(from=elpos[4,],to=elpos[1:3,],arr.side=2,path="H")
labs <- c("C","D","D","E")
text(0.55,0.4,expression(k[1]),font=3,adj=0,cex=0.8)
text(0.55,0.6,expression(k[2]),font=3,adj=0,cex=0.8)
for ( i in 1:4) textrect (elpos[i,],0.1,0.1,lab=labs[i],cex=1.5)
box(col="grey")
title("reversible reaction")
writelabel("A",line=0,at=-0.05)

# enzymatic reaction
openplotmat()
elpos<-coordinates (c(3,2,3))
elpos <- elpos[-c(5,6),]
elpos[4,1]<-0.3333
elpos[6,1]<-0.7

treearrow(from=elpos[1:2,],to=elpos[4,],arr.side=2,path="H")
treearrow(to=elpos[1:2,],from=elpos[4,],arr.side=2,path="H")
treearrow(from=elpos[3:4,],to=elpos[5:6,],arr.side=2,path="H")
labs <- c("E","D","F","I","E","G")
for ( i in 1:6) textrect (elpos[i,],0.075,0.07,lab=labs[i],cex=1.5)
text(0.35,0.6,expression(k[1]),font=3,adj=0,cex=0.8)
text(0.52,0.7,expression(k[2]),font=3,adj=0,cex=0.8)
text(0.72,0.3,expression(k[3]),font=3,adj=0,cex=0.8)
box(col="grey")
title("enzymatic reaction")
writelabel("B",line=0,at=-0.05)
subtitle()

########################################
# Fig 2.7. DIN-Phytoplankton
########################################

par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
# 1-st order reaction
openplotmat()
elpos<-coordinates (c(1,1))
tt<-straightarrow(from=elpos[1,],to=elpos[2,])
lab<-c(expression(NH[4]^"+"),"PHYTO")
for ( i in 1:2) textrect (elpos[i,],0.15,0.1,lab=lab[i],cex=1.5)
box(col="grey")
writelabel("A")

# 2nd order reaction between NH4 and phyto
openplotmat()
elpos<-coordinates (c(2,1))
tt<-treearrow(from=elpos[1:2,],to=elpos[3,],arr.side=2)
lab<-c(expression(NH[4]^"+"),"PHYTO","PHYTO")
for ( i in 1:3) textrect (elpos[i,],0.15,0.1,lab=lab[i],cex=1.5)
box(col="grey")
writelabel("B")
subtitle()

########################################
# Fig 2.8. Consumer-resource interactions
########################################


par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
basicplot("RESOURCE","WORKER",bcol=grey(0.95))
writelabel("A")

basicplot("PREY","PREDATOR",bcol=grey(0.95))
writelabel("B")

basicplot("NUTRIENT","ALGAE",bcol=grey(0.95))
writelabel("C")

basicplot("HOST","PARASITE",bcol=grey(0.95))
writelabel("D")
subtitle()

########################################
# Fig 2.9. Biochemical transformations
########################################

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))

bb<-basicplot("RESOURCE","SINK",mx=-0.05,sx=0.25,sy=0.08)
mm <- colMeans(bb)
textrect(mm+c(0.1,+0.04),0.25,0.08,lab="WORKER",shadow.size=0.02,box.col=grey(0.95))
writelabel("A",line=0,at=-0.05)

bb<-basicplot("SEMILABILE DOC","LABILE DOC",mx=-0.1,sx=0.35,sy=0.08)
mm <- colMeans(bb)
textrect(mm+c(0.1,+0.04),0.35,0.08,lab="BACTERIA",shadow.size=0.02,box.col=grey(0.95))
title("DOC hydrolysis")
writelabel("B",line=0,at=-0.05)

bb<-basicplot(expression(CO[2]),"Carbohydrates",mx=-0.1,sx=0.3,sy=0.08)
mm <- colMeans(bb)
textrect(mm+c(0.1,+0.04),0.35,0.08,lab=c("Photosynthetic","proteins"),shadow.size=0.02,box.col=grey(0.95))
title("Photosynthesis")
writelabel("C",line=0,at=-0.05)

bb<-basicplot("BACTERIA","DOC",mx=-0.05,sx=0.25,sy=0.08)
mm <- colMeans(bb)
textrect(mm+c(0.1,+0.04),0.25,0.08,lab="VIRUSES",shadow.size=0.02,box.col=grey(0.95))
title("Bacterial lysis")
writelabel("D",line=0,at=-0.05)
subtitle()

########################################
# Figure 2.10. Functional responses
########################################

# define functions
fr1 <- function(k ,resource)     pmin(1,resource/k)                  # functional response 1
fr2 <- function(ks,resource)     resource/(resource+ks)              # functional response 2
fr3 <- function(ks,pow,resource) resource^pow/(resource^pow+ks^pow)  # functional response 3

# parameters
k   <- 2
ks  <- 2
pow <- 3

# graphs
par(mar=c(5.1,4.1,4.1,2.1))
par (mfrow=c(2,2))     # figures in one row, three columns,

# draw the responses, resource levels (x) varying from 0 to 10
# functional response I - no axis labels and scales
curve(expr=fr1(k,x)    ,from=0,to=10, ylim=c(0,1),
      ylab="",xlab="",lwd=2,col="darkgrey",xaxt="n",yaxt="n")
axis(2)
mtext(side=3,line=1,"Type I",cex=1.2)         # text in upper margin
segments(k,0,k,1,lty=2)
segments(k,1,0,1,lty=2)
mtext(side=1,at=k,"k",cex=0.8)
text(5,0.5,expression(over(R,k)))

# functional response 2
curve(expr=fr2(ks=1,x)    ,from=0,to=10, ylim=c(0,1),
      ylab="",xlab="",lwd=2,col="darkgrey",xaxt="n",yaxt="n")
mtext(side=3,line=1,"Type II",cex=1.2)
text(5,0.5,expression(over(R,ks+R)))
segments(ks,0,ks,0.5,lty=2)
segments(ks,0.5,0,0.5,lty=2)
mtext(side=1,at=ks,"ks",cex=0.8)

# functional response 3
curve(expr=fr3(ks,pow,x),from=0,to=10, ylim=c(0,1),
      ylab="",xlab="",lwd=2,col="darkgrey",xaxt="n",yaxt="n")
mtext(side=3,line=1,"Type III",cex=1.2)
text(5,0.5,expression(over(R^p,ks^p+R^p)))
segments(ks,0,ks,0.5,lty=2)
segments(ks,0.5,0,0.5,lty=2)
mtext(side=1,at=ks,"ks",cex=0.8)

mtext(side=3,line=-2,outer=TRUE,"Functional responses",cex=1.5)    # text in upper margin
mtext(side=1,line=-2,outer=TRUE,"Resource, R",cex=1.2)
subtitle()


########################################
# Figure 2.11. Inhibition
########################################

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))

# schema
A <- matrix(nrow=2,ncol=2,0);A[2,1]<-1
pp <-plotmat(A,curve=0,pos=c(1,1),box.type="rect",name=c("RESOURCE","WORKER"),
        cex.txt=0,mx=-0.2,arr.pos=0.75,box.size=0.14,box.prop=0.75,
        box.col=c("white",grey(0.95)))
text(pp$arr$ArrowX,0.5,"X")
segments(pp$arr$ArrowX,0.5,0.7,0.5)
textrect(c(0.7,0.5),0.14,0.14*0.75,lab="INHIBITOR",lwd=2)
box(col="grey")
writelabel("A",line=0,at=-0.05)


# denitrification
openplotmat()
elpos<-coordinates (c(2,1))
treearrow(from=elpos[1:2,],to=elpos[3,],arr.side=2,path="H",line.pos=0.3,arr.pos=0.65)
labs <- c("OM",expression(NO[3]),expression(CO[2]~","~N[2]~",..."))

for ( i in 1:2) textrect (elpos[i,],0.14,0.1,lab=labs[i],cex=1.0,lwd=2)
textrect (elpos[3,],0.14,0.1,lab=labs[3],cex=1.0,lwd=2)
text(0.5,0.5,"X")
segments(0.5,0.5,0.9,0.5)
textrect(c(0.85,0.5),0.14,0.075,lab=expression(O[2]),lwd=2)
box(col="grey")
title("Denitrification")
writelabel("B",line=0,at=-0.05)


# Nitrate uptake
A <- matrix(nrow=2,ncol=2,0);A[2,1]<-1
pp <-plotmat(A,curve=0,pos=c(1,1),box.type="rect",name=c(expression(NO[3]),"ALGAE"),
        cex.txt=0,mx=-0.2,arr.pos=0.75,box.size=0.14,box.prop=0.75,
        box.col=c("white",grey(0.95)))
text(pp$arr$ArrowX,0.5,"X")
segments(pp$arr$ArrowX,0.5,0.7,0.5)
textrect(c(0.7,0.5),0.14,0.14*0.75,lab=expression(NH[4]^"+"),lwd=2)
box(col="grey")
title("Nitrate uptake")
writelabel("C",line=0,at=-0.05)
subtitle()

########################################
# fig. 2.12. Inhibition formulations
########################################
par(mfrow=c(2,2))
par(mar=c(4,4,3,1))
# 1-Monod like inhibition
kin <- 1
curve(expr=kin/(kin+x)  ,from=0,to=10,
      ylab="",xlab="",lwd=2,col="darkgrey",xaxt="n",yaxt="n")
axis(2,cex.axis=0.7)
text(5,0.5,expression(over(kin,kin+I)),cex=1)         # text in upper margin
segments(kin,0,kin,0.5,lty=2)
segments(kin,0.5,0,0.5,lty=2)
mtext(side=1,at=kin,"kin",cex=1)
writelabel("A")

# 2. exponential inhibition
curve(expr=exp(-0.5*x),from=0,to=10,
      ylab="",xlab="",lwd=2,col="darkgrey",xaxt="n",yaxt="n")
text(5,0.5,expression(e^{-c*I}))

mtext(side=3,line=1,outer=TRUE,"Inhibition",cex=1.25)    # text in upper margin
mtext(side=1,line=-1,outer=TRUE,"Inhibitory substance, I",cex=1.15)
writelabel("B")
subtitle()

########################################
# fig. 2.13 Coupled reactions
########################################
par(mfrow=c(2,2))

sx = 0.2
sy = 0.05
openplotmat()
elpos  <-coordinates (c(1,1,1),my=0.17)

straightarrow (to=elpos[3,],from=elpos[2,],arr.pos=0.72,lwd=2 )
straightarrow (to=c(0.5,0.1),from=elpos[3,],arr.pos=1,lwd=2)
textrect(elpos[2,],sx,sy,lab="PREY",shadow.size=0.02)
textrect(elpos[3,],sx,sy,lab="PREDATOR",shadow.size=0.02)
Arrows(0.75,0.97,elpos[2,1],0.72,code=2,arr.type="triangle",arr.adj=1,lwd=2)
Arrows(elpos[1,1],0.55,0.75,0.4,code=2,arr.type="triangle",arr.adj=1,lwd=2)
text(0.6,0.93,"mu",cex=1.3)
text(0.8,0.45,expression(1-gamma),cex=1.3,adj=0)
text(0.45,0.45,expression(gamma),cex=1.3,adj=0.5)
text(0.5,0.05,"r",cex=1.3,adj=0.5)
title("Source-sink coupling")
box(col="grey")
writelabel("A",line=0,at=-0.05)

openplotmat()
elpos<-coordinates (c(2,4))
elpos[1:2,2]<-0.85
elpos[3:6,2]<-0.15
tt<-treearrow(from=elpos[1:2,],to=elpos[3:6,],arr.side=2,path="H",line.pos=0.5,arr.pos=0.5)
labs <- c("OM",expression(O[2]),expression(CO[2]),expression(NH[3]),expression(H[3]~PO[4]),expression(H[2]~O))
for ( i in 1:6) textrect (elpos[i,],0.07,0.06,lab=labs[i],cex=1.0,lwd=2)
box(col="grey")
coeff<-c("1","106","106","16","1","106")
x2<-elpos[,1]
y2<-c(rep(0.75,2),rep(0.275,4))
for ( i in 1:6) textempty (c(x2[i],y2[i]),lab=coeff[i],font=2)

title("Stoichiometric coupling")
writelabel("B",line=0,at=-0.05)
subtitle()

########################################
# Figure 2.14. Relationships between flows
########################################
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
openplotmat()
straightarrow(from= c(0.1,0.5), to = c(0.8,0.5),arr.pos=0.85)
Arrows(0.45,0.5,0.63,0.85,arr.adj=1,lwd=2)
Arrows(0.35,0.5,0.47,0.17,arr.adj=1,lwd=2)
segmentarrow(c(0.8,0.5),c(0.7,0.92),dd=-0.15,arr.pos=0.5)

textrect(c(0.1,0.5),0.1,0.1,lab="PREY")
textrect(c(0.8,0.5),0.1,0.1,lab="PREDATOR")
textrect(c(0.7,0.925),0.07,0.07,lab=expression(CO[2]))
textrect(c(0.55,0.1),0.08,0.07,lab="DETRITUS")


text(0.22,0.55,"Ingestion",adj=0,font=3)
text(0.67,0.55,"Growth",adj=1,font=3)
textempty(c(0.9,0.955),c("Basal respiration"),font=3)
textempty(c(0.4,0.35),lab="Defecation",font=3)
textempty(c(0.55,0.65),lab="Activity respiration",font=3)
subtitle()

########################################
# Figure 2.15. Carrying capacity model
########################################
# 1. functions   #

# population density as a function of time and parameters
Nt <- function (r,K,N0,t) { K/(1+(K-N0)/N0*exp(-r*t))}

# the rate of change as a function of population density and parameters
dN <- function (r,K,N)    { r*N*(1-N/K)}

# 2. parameters  #

r <- 0.1   # intrinsic rate of increase /day
K <- 10    # density

N0a <- 0.1 # initial density, case a
N0b <- 15  # initial density, case b

# 3. graphs      #

par(mfrow=c(2,2),mar=c(4,4,3,1))

curve (dN(r,K,x),from=0,to=20,lwd=2,xlab="density, N",ylab="rate of change, dN/dt")
legend ("topright",expression(frac(dN,dt)==r*N*(1-frac(N,K))),bty="n")
writelabel("A",at=-0.13,line=0)

curve (Nt(r,K,N0a,x),from=0,to=100,lwd=2,xlab="time",ylab="density, N",ylim=c(0,15))
curve (Nt(r,K,N0b,x),add=TRUE,lwd=2,lty=2)

legend ("bottomright",c("K=10","r=0.1"))
#legend ("topright",expression(N[t]==frac(K,1+frac(K-N[0],N[0])*e^-rt)),bty="n")
writelabel("B",at=-0.13,line=0)
subtitle()

########################################
# Figure 2.16. Simplifications
########################################

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))

names <- c("", "Prey","Pred")
M <- matrix(nrow=3,ncol=3,byrow=TRUE,data=c(0,0,0,1,0,0,0,1,0))

pp<-plotmat(M,pos=c(1,1,1),curve=0,name=names,lwd=1,cex.txt=0.,mx=-0.2,
            box.lwd=2,box.size=0.12,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.6)
writelabel("A",line=0,at=-0.05)

selfarrow(pp$comp[3,]+c(0.45,0.05), lwd=1, arr.pos=0.25, arr.len=0,
          arr.type="triangle",,path="R", curve=c(0.1,0.1))
textrect(mid=pp$comp[3,]+c(0.45,0),radx=0.12,rady=0.06,lab="Pred",lwd=2)

box(col="grey")

mtext(side=3,line=1,cex=1.2,"Carrying capacity")


names <- c("A","B","C","D")
M <- matrix(nrow=4,ncol=4,byrow=TRUE,data=0)
M[2,1]<-1
M[3,2]<-1
M[4,3]<-1

pp<-plotmat(M,pos=c(1,1,1,1),curve=0,name=names,lwd=1,cex.txt=0.,mx=-0.2,
            box.lwd=2,box.size=0.12,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.6)
writelabel("B",line=0,at=-0.05)
p2 <-plotmat(M[1:2,1:2],pos=pp$comp[1:2,]+c(0.4,0.4,0,0),curve=0,name=names[1:2],lwd=1,
            cex.txt=0.,box.lwd=2,box.size=0.12,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.6,add=TRUE)
box(col="grey")
xy1<-p2$comp[2,1]+p2$radii[2,1]
x2 <-xy1[1]+0.05
segments(xy1[1],0.6,x2,0.6)     
Arrows(x2,0.6,x2,0.4,arr.type="triangle")
text(x2,0.345,"Mortality",font=3)
text(x2,0.295 ,"flux",font=3)

mtext(side=3,line=1,cex=1.2,"Closure")

# organic matter CLOSURE
names <- c("DETRITUS","DOC","BACTERIA","CO2")
M <- matrix(nrow=4,ncol=4,byrow=TRUE,data=0)
M[2,1]<-1
M[3,2]<-1
M[4,3]<-1

pp<-plotmat(M,pos=c(1,1,1,1),curve=0,name=names,lwd=1,cex.txt=0.,mx=-0.2,
            box.lwd=2,box.size=0.15,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.6)
writelabel("C",line=0,at=-0.05)
p2 <-plotmat(M[1:2,1:2],pos=pp$comp[c(1,4),]+c(0.4,0.4,0,0),curve=0,name=names[c(1,4)],
            lwd=1,cex.txt=0.,box.lwd=2,box.size=0.15,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.6,add=TRUE)

text(0.72,0.5,"Mineralisation",font=3,adj=0)
text(0.72,0.435 ,"flux",font=3,adj=0)
box(col="grey")
mtext(side=3,line=1,cex=1.2,"Shunt")

subtitle()

########################################
# Fig. 2.17 Temperature factor
########################################
par(mar=c(5.1,4.1,4.1,2.1))
# the temperature function:
par(mfrow=c(1,1))
Tt    <- seq(3,31, by= 4)
temp  <- seq(0,30, by=0.1)
y     <- exp((temp-9.2)/10*log(2))
plot(0,0,xlab="Temperature, dgC",ylab="Temperature factor (-)",
     ylim=c(0,4.5),xlim=c(0,30),type="n",main="Temperature acclimation")
polygon(c(temp,30,0),c(y,0,0),col=grey(0.95),border=NA)
for (i in Tt) curve(12.5*exp((i-10)/10*log(2))*dnorm(x,i,5),0,30,lwd=1,add=TRUE)
curve(exp((x-9.2)/10*log(2)),lwd=3,add=TRUE,0,30)
legend("topleft",c("individual species response","ecosystem response"),lty=1,lwd=c(1,2))
subtitle()

########################################
# Fig. 2.18 Light forcing function
########################################

par(oma = c(0,0,3,0))
par(mfrow=c(2,2))
# radiation as a sine with random component, period 1 year, paramters per day

meanRad   <-      530       #    muEinst/m2/s
ampRad    <-      0.85      #    -
phaseRad  <-      -1.4      #    -
periodRad <-      365       #    hours

tefoldRanRad <-   1.61      #    -         e-folding decay time scale of random light comp
ampRanRad    <-   520       #    W/m2      Amplitude of stochastic light component
tseq         <- 1:365
SolSine      <- meanRad *(1+ampRad*sin(2*pi*tseq/periodRad+phaseRad))
RadRan       <- 0
for ( i in 1:length(tseq))
 {
  RadRan <- RadRan*(1-1/tefoldRanRad)+ SolSine[i]/(meanRad+ampRad)*ampRanRad*(runif(1)-0.5)
  SolSine[i] <- SolSine[i]+RadRan}

plot(tseq,SolSine,xlab="day",ylab="muEinst/m2/s",main="Solar radiation",type="l")
writelabel("A")
# vertical extinction of light
x <- seq(0,100,0.1)
light <- 1200*exp(-0.075*x)
plot(light,x,ylim=c(100,0),ylab="depth,m",xlab="muEinst/m2/s",type="l",
     lwd=2,main="Light extinction in water")

text(200,50,expression (I==I[0]*e^{-lambda*I}))
writelabel("B")

# light inhibition function

iopt  <- 200
beta  <- 0.005
pmax <- 2
curve(expr=pmax*(2*(1+beta)*x/iopt)/((x/iopt)^2+2*beta*x/iopt+1),-50,500,xlab="muEinst/m2/s",ylab="/d",
      main="Photosynthesis versus light",lwd=2)
text(50,0.6,"quasi-linear",adj=0,font=3,cex=1.2)
text(130,1.75,"saturation",adj=0,font=3,cex=1.2)
text(365,1.9,"inhibition",adj=0,font=3,cex=1.2)
writelabel("C")

# Chlorophyll maximum
par(mar=c(1.5,4.1,1.5,2.1))
x <- seq(0,100,0.1)
light <- 1*exp(-0.075*x)
chl   <- dnorm(x,40,8)
chl   <- chl/(max(chl)*1.1)
ksnit <- 10
nit <- (x-40)/(x-40+ksnit)
nit[x<40] <- 0
plot(light,x,ylim=c(100,0),ylab="depth,m",xaxt="n",xlab="",type="l",lwd=2,axes=FALSE)
axis(side=2)
polygon(c(0,light),c(0,x),col="lightgrey")
polygon(c(nit,0,0),c(x,100,0),col="grey")
polygon(chl,x,col="darkgrey")
lines(light,x)
lines(nit,x,lwd=2,lty=2)
text(0.5,5,"PAR",font=2)
text(0.5,40,"Chlorophyll",font=2)
text(0.5,90,"Nitrate",font=2)
writelabel("D")
mtext(outer=TRUE,line=1,cex=1.3,"Forcing function light")
subtitle()

########################################
# Fig. 2.19 Light limitation functions
########################################

col   <- "black"          ## color
lty   <- 1:10             ## line type
par(mar=c(5.1,4.1,4.1,2.1))
par(oma = c(0,0,3,0),mfrow=c(2,2))   ## outer margin, multiple figures 2 rows, 2 cols

# STEELE function
#----------------
# exprerssion to write the formula in the title
exp   <- expression (frac(I ,iopt) %*% e^(1-I/iopt))

# parameter
iopt  <- 150
curve(x/iopt*exp(1-x/iopt),0,500,lwd=2,xlab=expression("I, muEinst"~ m^{-2}~s^{-1}),ylab="-",main="1-par" )
legend("bottomright","Steele",lty=1,lwd=2)
writelabel("A")

# EVANS function
#----------------
ks<-50
exp <- expression (frac(I ,sqrt(I^2 + iopt^2)))
curve(x/ sqrt(x*x+iopt*iopt),0,500,lwd=2,xlab=expression("I, muEinst"~ m^{-2}~s^{-1}),ylab="-",main="1-par" )
curve(x/(ks+x),add=TRUE,lwd=2,lty=2)
legend("bottomright",c("Evans","Monod"),lty=1:2,lwd=2)
writelabel("B")

# PLATT, Eelers-Peeters function
#----------------

alpha     <-0.008
beta      <- 0.005
pmax      <- 2

curve(pmax*(1-exp(-x*alpha/pmax)) * exp(-beta*x/pmax),0,500,lwd=2,
     xlab=expression("I, muEinst"~ m^{-2}~s^{-1}),ylab="-",main="2-par" ,ylim=c(0,1))
curve((2*(1+beta)*x/iopt)/((x/iopt)^2+2*beta*x/iopt+1),add=TRUE,lwd=2,lty=2)
legend ("bottomright",c("Platt","Eilers-Peeters"),lty=1:2,lwd=2)
writelabel("C")

mtext(outer=TRUE,line=1,cex=1.3,"Light limitation functions")
subtitle()

#----------------------------#
#      Example figures       #
#----------------------------#

########################################
# Fig. 2.20 NPZD model
########################################
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
names <- c("PHYTO","DIN","ZOO","DETRITUS")
M <- matrix(nrow=4,ncol=4,byrow=TRUE,data=c(
#   p n z  d
    0,1,0, 0, #p
    0,0,4, 6,#n
    2,0,0, 0, #z
    0,0,5,0  #d
    ))

pp<-plotmat(M,pos=c(1,2,1),curve=0,name=names,lwd=1,my=0.0,cex.txt=0.8,prefix="f",
            box.lwd=2,box.size=0.08,box.type="square",box.prop=0.5,
            arr.type="triangle",arr.pos=0.6,shadow.size=0.01,main="NPZD")

# extra arrows: flow 5 to Detritus and flow 2 to detritus
phyto   <-pp$comp[names=="PHYTO"]
zoo     <-pp$comp[names=="ZOO"]
nh3     <-pp$comp[names=="DIN"]
detritus<-pp$comp[names=="DETRITUS"]

# flow2->detritus
m2 <- 0.5*(zoo+phyto)
m1 <- detritus
m1[1] <-m1[1]+ pp$radii[3,1]*0.2
m1[2] <-m1[2] + pp$radii[3,2]
mid<-straightarrow (to=m1,from=m2,arr.type="triangle",arr.pos=0.7,lwd=1)
text(mid[1]-0.01,mid[2]+0.03,"f3",cex=0.8)

# solar radiation
m1 <- 0.5*(nh3+phyto)
m2 <- c(0.25,0.8)
segments (m1[1],m1[2],m2[1],m2[2],lwd=1,lty=2)
text(m2[1]-0.01,m2[2]+0.03,"solar radiation",adj=c(0.5,0.5))

# chlorophyll
m1 <- phyto
m1[1] <-m1[1]+ pp$radii[1,1]
m2 <- m1
m2[1]<-m2[1]+0.25
segments (m1[1],m1[2],m2[1],m2[2],lwd=1)
textellipse(m2,pp$radii[1,1],pp$radii[1,2],lwd=1,shadow.size=0,lab="Chlorophyll")
subtitle()

################################################################################
# Fig. 2.21. NPZD model solution
################################################################################

#----------------------#
# the model equations: #
#----------------------#

NPZD<-function(t,state,parameters)
 {
 with(as.list(c(state,parameters)),{  # unpack the state variables, parameters
    # Light, a sine function; 50% of light is PAR
    PAR            <- 0.5*(540+440*sin(2*pi*t/365-1.4))
    din            <- max(0,DIN)      # to avoid errors when DIN becomes slightly negative..
    Nuptake        <- maxUptake * PAR/(PAR+ksPAR) * din/(din+ksDIN)*PHYTO

    Grazing        <- maxGrazing* PHYTO/(PHYTO+ksGrazing)*ZOO
    Faeces         <- pFaeces * Grazing
    Excretion      <- excretionRate * ZOO
    Mortality      <- mortalityRate * ZOO * ZOO
    Mineralisation <- mineralisationRate * DETRITUS
    Chlorophyll    <- chlNratio * PHYTO
    TotalN         <- PHYTO + ZOO + DETRITUS + DIN

    dPHYTO    <- Nuptake - Grazing
    dZOO      <- Grazing - Faeces - Excretion - Mortality
    dDETRITUS <- Mortality - Mineralisation + Faeces
    dDIN      <- Mineralisation + Excretion - Nuptake

    # the output, packed as a list
    list(c(dPHYTO,dZOO,dDETRITUS,dDIN),                          # the rate of change
        c(Chlorophyll = Chlorophyll, PAR=PAR, TotalN= TotalN))   # the ordinary output variables
    })
  }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

parameters<-c(maxUptake          =1.0,       # /day
              ksPAR              =140,       # muEinst/m2/s
              ksDIN              =0.5,       # mmolN/m3
              maxGrazing         =1.0,       # /day
              ksGrazing          =1.0,       # mmolN/m3
              pFaeces            =0.3,       # -
              excretionRate      =0.1,       # /day
              mortalityRate      =0.4,       # /(mmolN/m3)/day
              mineralisationRate =0.1,       # /day
              chlNratio          =1)         # mgChl/mmolN

#-------------------------#
# the initial conditions: #
#-------------------------#

state     <-c(PHYTO   =1,                    # state variable initial conditions, units mmolN/m3
              ZOO     =0.1,
              DETRITUS=5.0,
              DIN     =5.0)

#----------------------#
# RUNNING the model:   #
#----------------------#
#  2 steps
#  step 1  : run the model for 365 days, no intermediate output required

times     <-c(0,365)
out       <-as.data.frame(ode(state,times,NPZD,parameters)) # ode is integrator

#  step 2  : the model is reinitialised with the final conditions of previous run

num       <- length(out$PHYTO)   # last element
state     <- c(PHYTO   =out$PHYTO[num],
               ZOO     =out$ZOO[num],
               DETRITUS=out$DETRITUS[num],
               DIN     =out$DIN[num])

# new run; this time the output is kept: time: from 1 to 365 days, outputsteps of days
times     <-seq(0,730,by=1)
out       <-as.data.frame(ode(state,times,NPZD,parameters))

#------------------------#
# PLOTTING model output: #
#------------------------#
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(2,2), oma=c(0,0,3,0))   # set number of plots (mfrow) and margin size (oma)

plot (times,out$PAR        ,type="l",lwd=2,main="PAR"        ,xlab="time, hours",ylab="muEinst/m2/s")
writelabel("A")
plot (times,out$Chlorophyll,type="l",lwd=2,main="Chlorophyll",xlab="time, hours",ylab="mug/l")
writelabel("B")
plot (times,out$ZOO        ,type="l",lwd=2,main="Zooplankton",xlab="time, hours",ylab="mmolN/m3")
writelabel("C")
plot (times,out$DIN        ,type="l",lwd=2,main="DIN"        ,xlab="time, hours",ylab="mmolN/m3")
writelabel("D")

mtext(outer=TRUE,side=3,"NPZD model",cex=1.5)

subtitle()
########################################
# Fig. 2.22 AQUAPHY
########################################
par(mfrow=c(1,1))
par(mar=c(1,1,1,1),oma=c(0,0,0,0))
openplotmat(main="AQUAPHY")
textrect(c(0.45,0.5),0.35,rady=0.3,lab="",shadow.size=0.02,cex=1.5)
segments(0.1,0.5,0.8,0.5)
segments(0.45,0.5,0.45,0.8)

Arrows(0.53,0.75,0.37,0.75,arr.length= 0.5,arr.width=0.2,arr.type="triangle",arr.adj=1)
Arrows(0.37,0.67,0.53,0.67,arr.length= 0.5,arr.width=0.2,arr.type="triangle",arr.adj=1)
text(0.55,0.75,"f3")
text(0.35,0.67,"f4")

text(0.13,0.65,adj=0,"RESERVE",cex=1.15)
text(0.13,0.6,adj=0,"Carbohydrates",cex=1.15)

text(0.77,0.65,adj=1,"LMW",cex=1.15)
text(0.77,0.6,adj=1,"Carbohydrates",cex=1.15)

text(0.45,0.35,adj=0.5,"Biosynthetic and Photosynthetic",cex=1.15)
text(0.45,0.3,adj=0.5,"PROTEINS",cex=1.15)

Arrows(0.65,0.55,0.65,0.4,arr.length= 0.5,arr.width=0.22,arr.type="triangle",arr.adj=1)
text(0.68,0.4,"f5")

textrect(c(0.93,0.525),0.05,rady=0.05,lab="DIN",shadow.size=0.015,cex=1.5)
segments(0.88,0.52,0.65,0.52)

Arrows(0.65,0.9,0.65,0.7,arr.length= 0.5,arr.width=0.25,,arr.type="triangle",arr.adj=1)
text(0.63,0.88,"f1")

Arrows(0.65,0.85,0.75,0.85,arr.length= 0.5,arr.width=0.25,arr.type="triangle",arr.adj=1)
text(0.77,0.85,"f2")

Arrows(0.8,0.75,0.9,0.75,arr.length= 0.5,arr.width=0.25,arr.type="triangle",arr.adj=1)
text(0.92,0.75,"f6")

segments(0.1,0.5,0.05,0.5)
bentarrow(c(0.1,0.5),c(0.02,0.35),arr.length= 0.5,arr.width=0.25,arr.pos=1,
          arr.type="triangle",arr.adj=1)
text(0.025,0.32,"f7")

text(0.11,0.51,adj=c(0,0),"C",cex=1.3,col="darkgrey")
text(0.46,0.51,adj=c(0,0),"C",cex=1.3,col="darkgrey")

text(0.11,0.21,adj=c(0,0),"C,N,Chl",cex=1.3,col="darkgrey")
subtitle()

################################################################################
# Fig. 2.23. AQUAPHY model solution
################################################################################

#----------------------#
# the model equations: #
#----------------------#

model<-function(t,state,parameters)
  {
  with(as.list(c(state,parameters)),{  # unpack the state variables, parameters

    # PAR, on-off function depending on the hour within a day
    hourofday       <- t%%24
    PAR <- ifelse (hourofday  < dayLength, parMean , 0)

    ## the output variables
    PhytoC           <- PROTEIN + RESERVE + LMW       # all components contain carbon
    PhytoN           <- PROTEIN * rNCProtein          # only proteins contain nitrogen
    NCratio          <- PhytoN / PhytoC
    Chlorophyll      <- PhytoN * rChlN
    TotalN           <- PhytoN + DIN
    ChlCratio        <- Chlorophyll / PhytoC

    ## the rates, in mmol/hr
    PartLMW          <- LMW / PhytoC
    Limfac           <- max(0,min(1,(maxpLMW -PartLMW)/(maxpLMW-minpLMW)))
    PhotoSynthesis   <- maxPhotoSynt*Limfac*(1-exp(alpha*PAR/maxPhotoSynt)) * PROTEIN
    Exudation        <- pExudation * PhotoSynthesis
    MonodQuotum      <- max(0,LMW / PROTEIN - minQuotum)
    ProteinSynthesis <- maxProteinSynt*MonodQuotum * DIN / (DIN+ksDIN)      * PROTEIN
    Storage          <- maxStorage    *MonodQuotum                          * PROTEIN
    Respiration      <- respirationRate * LMW + pResp*ProteinSynthesis
    Catabolism       <- catabolismRate  * RESERVE

    ## the rates of change of state variables; includes dilution effects (last term)
    dLMW     <- ( PhotoSynthesis + Catabolism
                - Exudation - Storage  - Respiration - ProteinSynthesis
                - dilutionRate * LMW)

    dRESERVE <-  Storage - Catabolism          - dilutionRate * RESERVE

    dPROTEIN <-  ProteinSynthesis              - dilutionRate * PROTEIN

    dDIN     <- -ProteinSynthesis * rNCProtein - dilutionRate * (DIN - inputDIN)


    ## the output, as a list
    list(c(dDIN,dPROTEIN,dRESERVE,dLMW),              ## the rate of change of state variables
           c(PAR               = PAR,                 ## the ordinary variables
             TotalN            = TotalN,
             PhotoSynthesis    = PhotoSynthesis,
             NCratio           = NCratio,
             ChlCratio         = ChlCratio,
             Chlorophyll       = Chlorophyll))
    })
 }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

parameters<-c(maxPhotoSynt   =0.125,      #molC/molC/hr      Maximal protein C-specific rate of photsynthesis at 20 dg
              rMortPHY       =0.001,      #/hr               Mortality rate of Phytoplankton (lysis and zooplankton grazing)
              alpha          =-0.125/150, #muEinst/m2/s/hr   Light dependency factor
              pExudation     =0.0,        #-                 Part of photosynthesis that is exudated
              maxProteinSynt =0.136,      #molC/molC/hr      Maximal Biosynthetic C-specific N-uptake rate
              ksDIN          =1.0,        #mmolN/m3          Half-saturation ct of N uptake Phytoplankton
              minpLMW        =0.05,       #molC/molC         Minimum metabolite/totalC ratio in algae
              maxpLMW        =0.15,       #molC/molC         Maximum metabolite/totalC ratio in algae
              minQuotum      =0.075,      #molC/molC         Minimum metabolite/Protein ratio for synthesis
              maxStorage     =0.23,       #/h                Maximum storage rate for Phytoplankton
              respirationRate=0.0001,     #/h                Respiration rate of LMW
              pResp          =0.4,        #-                 Part of protein synthesis that is respired (cost of biosynthesis)
              catabolismRate =0.06,       #/h                Catabolism rate of Phytoplankton reserves
              dilutionRate   =0.01,       #/h                dilution rate in chemostat
              rNCProtein     =0.2,        #molN/molC         Nitrogen/carbon ratio of proteins
              inputDIN       =10.0,       #mmolN/m3          DIN in inflowing water
              rChlN          =1,          #gChl/molN         Chl to nitrogen ratio
              parMean        =250.,       #mumolPhot/m2/s    PAR during the light phase
              dayLength      =15.         #hours             Length of illuminated period
              )

#-------------------------#
# the initial conditions: #
#-------------------------#

# assume the amount of reserves = 50% amount of proteins
# 10% LMW

state     <-c(DIN     =6.,     #mmolN/m3
              PROTEIN =20.0,   #mmolC/m3
              RESERVE =5.0,    #mmolC/m3
              LMW     =1.0)    #mmolC/m3

#----------------------#
# RUNNING the model:   #
#----------------------#

times <-seq(0,24*10,1)
out   <-as.data.frame(ode(state,times,model,parameters))


#------------------------#
# PLOTTING model output: #
#------------------------#
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(2,2), oma=c(0,0,3,0))         # set number of plots (mfrow) and margin size (oma)
col <- grey(0.9)
ii <- 1:length(out$PAR)                   # output over entire period

plot (times[ii],out$Chlorophyll[ii],type="l",main="Chlorophyll",xlab="time, hours",ylab="mug/l")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$Chlorophyll[ii]  ,lwd=2 )


plot (times[ii],out$DIN[ii]        ,type="l",main="DIN"        ,xlab="time, hours",ylab="mmolN/m3")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$DIN[ii]  ,lwd=2 )


plot (times[ii],out$NCratio[ii]    ,type="n",main="NCratio"    ,xlab="time, hours",ylab="molN/molC")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$NCratio[ii]  ,lwd=2 )


plot (times[ii],out$PhotoSynthesis[ii],type="l",main="PhotoSynthesis",xlab="time, hours",ylab="mmolC/m3/hr")
polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
lines (times[ii],out$PhotoSynthesis[ii]  ,lwd=2 )

mtext(outer=TRUE,side=3,"AQUAPHY",cex=1.5)
subtitle()

###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################

########################################
# Fig. 2.24 CARRYING CAPACITY
########################################

par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
curve (expr= r*x*(1-x/K),from=0,to=20,lwd=2,
       xlab="density, N",ylab="rate of change, dN/dt")
legend ("topright",c("K=10","r=0.1"))
subtitle()

########################################
# Fig. 2.25 MONOD FUNCTION, SEVERAL ks
########################################

par(mar=c(5.1,4.1,4.1,2.1))
food <- seq(0,30,by=0.1)
ks   <- seq (1,10, by=2)

foodfun  <- outer(food,ks,function(x,y) x /(y +x))

matplot(x=food,foodfun,type="l",lty=1:10, col=1,
        xlab="food",ylab="-",
        main= expression (frac(food ,food+ks)))

legend ("bottomright", as.character(ks), title="ks=",
        col=1,lty=1:10)
subtitle()

###############################################################################
####======================================================================#####
####                                Projects                              #####
####======================================================================#####
###############################################################################

########################################
# Fig. 2.26 PHYTODIN
########################################

par(mar=c(0,0,0,0))
emptyplot()
rect(0.1,0.1,0.9,0.9)
rect(0.11,0.11,0.89,0.7,border=NA,col="lightgrey")
elpos <- matrix(nrow=2,data=c(0.25,0.75,0.4,0.4))
#straightarrow (from=elpos[1,],to=elpos[2,],arr.pos=0.5,lwd=2)
textrect(elpos[1,],0.1,0.1,lab="PHYTO",shadow.size=0.02)
textrect(elpos[2,],0.1,0.1,lab="DIN",shadow.size=0.02)
subtitle()

########################################
# Fig. 2.27. Detritus degradation
########################################


par(mar=c(0,0,0,0))
A <- matrix(nrow=4,ncol=4,data=0)
b<-plotmat(A,pos=c(1,2,1),box.type="rect",box.size=0.12,box.prop=0.7,name=c("BACT","POC", "LMWC", "HMWC"))
subtitle()

########################################
# Fig. 2.28 Autocatalytic degradation
########################################


par(mar=c(0,0,0,0))
openplotmat()
elpos<-coordinates (c(2,3))
treearrow(from=elpos[1:2,],to=elpos[3:5,],arr.side=2)
lab<-c("A","B","B","B","C")
for ( i in 1:5) textrect (elpos[i,],0.1,0.1,lab=lab[i],cex=1.5)
box(col="grey")

subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

