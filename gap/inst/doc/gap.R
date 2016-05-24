### R code from vignette source 'gap.Rnw'

###################################################
### code chunk number 1: gap.Rnw:172-176
###################################################
library(gap)
search()
lsf.str("package:gap")
data(package="gap")$results


###################################################
### code chunk number 2: gap.Rnw:263-270
###################################################
# pedigree diagram
data(lukas,package="gap")
library(kinship2)
ped <- with(lukas,pedigree(id,father,mother,sex))
pdf("figures/lukas.pdf",height=14,width=15)
plot(ped)
dev.off()


###################################################
### code chunk number 3: gap.Rnw:280-306
###################################################
# unordered individuals
library(gap)
gk1 <- kin.morgan(lukas)
write.table(gk1$kin.matrix,"results/gap_1.txt",quote=FALSE)

library(kinship2)
kk1 <- kinship(lukas[,1],lukas[,2],lukas[,3])
write.table(kk1,"results/kinship_1.txt",quote=FALSE)

d <- gk1$kin.matrix-kk1
sum(abs(d))

# order individuals so that parents precede their children
library(pedigree)
op <- orderPed(lukas)
olukas <- lukas[order(op),]
gk2 <- kin.morgan(olukas)

write.table(olukas,"olukas.csv",quote=FALSE)
write.table(gk2$kin.matrix,"results/gap_2.txt",quote=FALSE)

kk2 <- kinship(olukas[,1],olukas[,2],olukas[,3])
write.table(kk2,"results/kinship_2.txt",quote=FALSE)

z <- gk2$kin.matrix-kk2
sum(abs(z))


###################################################
### code chunk number 4: gap.Rnw:314-352
###################################################
library(gap)
models <- matrix(c(
         4.0, 0.01,
         4.0, 0.10,
         4.0, 0.50, 
         4.0, 0.80,
         2.0, 0.01,
         2.0, 0.10,
         2.0, 0.50,
         2.0, 0.80,
         1.5, 0.01,    
         1.5, 0.10,
         1.5, 0.50,
         1.5, 0.80), ncol=2, byrow=TRUE)
outfile <- "fbsize.txt"
cat("gamma","p","Y","N_asp","P_A","H1","N_tdt","H2","N_asp/tdt","L_o","L_s\n",
    file=outfile,sep="\t")
for(i in 1:12) {
    g <- models[i,1]
    p <- models[i,2]
    z <- fbsize(g,p)
    cat(z$gamma,z$p,z$y,z$n1,z$pA,z$h1,z$n2,z$h2,z$n3,z$lambdao,z$lambdas,
        file=outfile,append=TRUE,sep="\t")
    cat("\n",file=outfile,append=TRUE)
}
table1 <- read.table(outfile,header=TRUE,sep="\t")
nc <- c(4,7,9)
table1[,nc] <- ceiling(table1[,nc])
dc <- c(3,5,6,8,10,11)
table1[,dc] <- round(table1[,dc],2)
unlink(outfile)
# APOE-4, Scott WK, Pericak-Vance, MA & Haines JL
# Genetic analysis of complex diseases 1327
g <- 4.5
p <- 0.15
cat("\nAlzheimer's:\n\n")
fbsize(g,p)
table1


###################################################
### code chunk number 5: gap.Rnw:357-385
###################################################
library(gap)
kp <- c(0.01,0.05,0.10,0.2)
models <- matrix(c(
          4.0, 0.01,
          4.0, 0.10,
          4.0, 0.50, 
          4.0, 0.80,
          2.0, 0.01,
          2.0, 0.10,
          2.0, 0.50,
          2.0, 0.80,
          1.5, 0.01,    
          1.5, 0.10,
          1.5, 0.50,
          1.5, 0.80), ncol=2, byrow=TRUE)
outfile <- "pbsize.txt"
cat("gamma","p","p1","p5","p10","p20\n",sep="\t",file=outfile)
for(i in 1:dim(models)[1])
{
   g <- models[i,1]
   p <- models[i,2]
   n <- vector()
   for(k in kp) n <- c(n,ceiling(pbsize(k,g,p)))
   cat(models[i,1:2],n,sep="\t",file=outfile,append=TRUE)
   cat("\n",file=outfile,append=TRUE)
} 
table5 <- read.table(outfile,header=TRUE,sep="\t")
table5


###################################################
### code chunk number 6: gap.Rnw:390-435
###################################################
library(gap)
# ARIC study
outfile <- "aric.txt"
n <- 15792
pD <- 0.03
p1 <- 0.25
alpha <- 0.05
theta <- c(1.35,1.40,1.45)
beta1 <- 0.8
s_nb <- c(1463,722,468)
cat("n","pD","p1","hr","q","power","ssize\n",file=outfile,sep="\t")
for(i in 1:3)
{
  q <- s_nb[i]/n
  power <- ccsize(n,q,pD,p1,alpha,log(theta[i]))
  ssize <- ccsize(n,q,pD,p1,alpha,log(theta[i]),beta1)
  cat(n,"\t",pD,"\t",p1,"\t",theta[i],"\t",q,"\t",signif(power,3),"\t",ssize,"\n",
      file=outfile,append=TRUE)
}
read.table(outfile,header=TRUE,sep="\t")
unlink(outfile)
# EPIC study
outfile <- "epic.txt"
n <- 25000
alpha <- 0.00000005
power <- 0.8
s_pD <- c(0.3,0.2,0.1,0.05)
s_p1 <- seq(0.1,0.5,by=0.1)
s_hr <- seq(1.1,1.4,by=0.1)
cat("n","pD","p1","hr","alpha","ssize\n",file=outfile,sep="\t")
# direct calculation
for(pD in s_pD)
{
   for(p1 in s_p1)
   {
      for(hr in s_hr)
      {
         ssize <- ccsize(n,q,pD,p1,alpha,log(hr),power)
         if (ssize>0) cat(n,"\t",pD,"\t",p1,"\t",hr,"\t",alpha,"\t",ssize,"\n",
                          file=outfile,append=TRUE)
      }
   }
}
read.table(outfile,header=TRUE,sep="\t")
unlink(outfile)


###################################################
### code chunk number 7: gap.Rnw:443-451
###################################################
library(gap)
pdf("figures/qqunif.pdf",height=10,width=10)
u_obs <- runif(1000)
r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
u_exp <- r$y
hits <- u_exp >= 2.30103
points(r$x[hits],u_exp[hits],pch=21,bg="green")
dev.off()


###################################################
### code chunk number 8: gap.Rnw:459-475
###################################################
library(gap)
load("4w.rda")
ord <- with(d,order(chr,pos))
d <- d[ord,]
pdf("figures/4w.pdf",height=9,width=10)
oldpar <- par()
par(cex=0.6)
colors <- c(rep(c("blue","red"),15),"red")
mhtplot(d,control=mht.control(colors=colors,gap=1000),pch=19,srt=0)
axis(2,cex.axis=2)
suggestiveline <- -log10(3.60036E-05)
genomewideline <- -log10(1.8E-06)
abline(h=suggestiveline, col="blue")
abline(h=genomewideline, col="green")
abline(h=0)
dev.off()


###################################################
### code chunk number 9: gap.Rnw:481-497
###################################################
library(gap)
png("figures/mhtplot.png",height=10,width=16,units="cm",res=300)
data <- with(mhtdata,cbind(chr,pos,p))
glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3","PPP1R3B",
           "RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
hdata <- subset(mhtdata,gene%in%glist)[c("chr","pos","p","gene")]
color <- rep(c("lightgray","gray"),11)
glen <- length(glist)
hcolor <- rep("red",glen)  
par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
ops <- mht.control(colors=color,yline=1.5,xline=3)
hops <- hmht.control(data=hdata,colors=hcolor)
mhtplot(data,ops,hops,pch=19)
axis(2,pos=2,at=1:16)
title("Manhattan plot with genes highlighted",cex.main=1.8)
dev.off()


###################################################
### code chunk number 10: gap.Rnw:508-513
###################################################
library(gap)
pdf("figures/asplot.pdf",height=14,width=14)
asplot(CDKNlocus, CDKNmap, CDKNgenes, best.pval=5.4e-8, sf=c(3,6))
title("CDKN2A/CDKN2B Region")
dev.off()


###################################################
### code chunk number 11: gap.Rnw:523-533
###################################################
library(gap)
pdf("figures/ESplot.pdf",height=10,width=10)
options(stringsAsFactors=FALSE)
testdata <- data.frame(models=c("Basic model","Adjusted","Moderately adjusted",
                       "Heavily adjusted","Other"),
OR = c(4.5,3.5,2.5,1.5,1),
SElogOR = c(0.2,0.1,0.5,0.5,0.2))
ESplot(testdata,v=1)
title("This is a fictitious plot")
dev.off()


