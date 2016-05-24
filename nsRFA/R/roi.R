roi <- function (p.ungauged,p.gauged,cod.p,x=NULL,cod=NULL) {

 parametri <- rbind(p.ungauged,p.gauged)
 n <- dim(parametri)[1]
 k <- dim(parametri)[2]
 matrice.distanze <- dist(parametri)
 distanze <- as.matrix(matrice.distanze)[-1,1]
 ordini <- order(distanze)
 if(!((is.null(x))&(is.null(cod)))) {
  dum.reg <- factor(cod,levels=cod.p)
  x <- x[!is.na(dum.reg)]
  cod <- cod[!is.na(dum.reg)]
  n.reg <- tapply(x,cod,length)
  l.reg <- t(sapply(split(x,cod),Lmoments))
  tabella <- data.frame(cod.p,cbind(distanze,n.reg,l.reg))
  names(tabella)[1:3] <- c("cod","dist","n")
  output <- tabella[ordini,]
 }
 else {
  tabella <- data.frame(cod.p,distanze)
  names(tabella) <- c("cod","dist")
  output <- tabella[ordini,]
 }
 return(output)
}


# -------------------------------------------------------------------- #

roi.hom <- function (p.ungauged,p.gauged,cod.p,x,cod,test="HW",limit=2,Nsim=500,index=2) {

 parametri <- rbind(p.ungauged,p.gauged)
 n <- dim(parametri)[1]
 k <- dim(parametri)[2]
 matrice.distanze <- dist(parametri)
 distanze <- as.matrix(matrice.distanze)[-1,1]
 names(distanze) <- cod.p
 distanze.ordinate <- sort(distanze)

 i=20
 t=FALSE
 if (test=="HW") {
  while((t==FALSE)&&(i>1)) {
   bacini <- names(distanze.ordinate)[1:i]
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   H1 <- HW.tests(x.reg,cod.reg,Nsim)[1]
   if (H1<=limit) t=TRUE
   i=i-1
  }
 }
 else if (test=="AD") {
  while((t==FALSE)&&(i>1)) {
   bacini <- names(distanze.ordinate)[1:i]
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   P.AD <- ADbootstrap.test(x.reg,cod.reg,Nsim,index)[2]
   if (P.AD<=limit) t=TRUE
   i=i-1
  }
 }
 if (i <= 1) regione <- NULL
 else regione <- bacini
 return(regione)
}


# -------------------------------------------------------------------- #

roi.st.year <- function (p.ungauged,p.gauged,cod.p,x,cod,test="HW",station.year=500,Nsim=500,index=2) {

 parametri <- rbind(p.ungauged,p.gauged)
 n <- dim(parametri)[1]
 k <- dim(parametri)[2]
 matrice.distanze <- dist(parametri)
 distanze <- as.matrix(matrice.distanze)[-1,1]
 names(distanze) <- cod.p
 distanze.ordinate <- sort(distanze)

 ni <- tapply(x,cod,length)
 sum.ni <- cumsum(ni[names(distanze.ordinate)])
 bacini <- names(sum.ni)[1:(sum(sum.ni<station.year)+1)]

 if (test=="HW") {
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   H.HW <- HW.tests(x.reg,cod.reg,Nsim)
   homtest <- H.HW
 }
 else if (test=="AD") {
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   P.AD <- ADbootstrap.test(x.reg,cod.reg,Nsim,index)
   homtest <- P.AD
 }
 else if (test=="HW and AD") {
   dum.reg <- factor(cod,levels=bacini)
   x.reg <- x[!is.na(dum.reg)]
   cod.reg <- cod[!is.na(dum.reg)]
   H.HW <- HW.tests(x.reg,cod.reg,Nsim)
   P.AD <- ADbootstrap.test(x.reg,cod.reg,Nsim,index)
   homtest <- c(H.HW,P.AD)
 }
 li <- t(sapply(split(x.reg,cod.reg),Lmoments))[bacini,]
 ni <- tapply(x.reg,cod.reg,length)[bacini]
 Si <- 1 - c(0,cumsum(ni)[-length(ni)])/sum(ni) #1 - cumsum(ni)/sum(ni)
 wi <- ni*Si
 di <- distanze[bacini]
 tabella <- data.frame(cbind(li,ni,Si,wi,di))
 Rli <- apply(li*matrix(wi, nrow=dim(li)[1], ncol=dim(li)[2]),2,sum)/sum(wi)

 roiout <- list(region=bacini,test=homtest,char.sites=tabella,regLmom=Rli)

 return(roiout)
}
