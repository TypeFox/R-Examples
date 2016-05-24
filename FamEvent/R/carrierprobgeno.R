 # carrier probability

carrierprobgeno <-function(data, method="data", mode="dominant", q=0.02){
  carrp <- data$mgene
  id.na<-data$indID[is.na(data$mgene)]
  mut.ca <- data$relation[data$mgene==1 & !is.na(data$mgene)]
  
  cfam.id <- data$famID[data$proband==1 & data$mgene==1]
  nfam.id <- data$famID[data$proband==1 & data$mgene==0]
  i.cfam <- is.element(data$famID,cfam.id)
  i.nfam <- is.element(data$famID,nfam.id)
  
  if(method=="data"){
    for(g in 0:3){
      for(s in c(0,1)){
        # carrier families
        carrp[i.cfam & is.na(data$mgene) & data$generation==g & data$gender==s] <- mean(data$mgene[i.cfam & !is.na(data$mgene) & data$generation==g & data$gender==s])
        # non-carrier famiiies
        carrp[i.nfam & is.na(data$mgene) & data$generation==g & data$gender==s] <- mean(data$mgene[i.nfam & !is.na(data$mgene) & data$generation==g & data$gender==s])
      }
    }
  }
  else if (method=="mendelian"){
  	
  	if(is.null(q)) q <- ifelse(mode=="recessive", sqrt(mean(data$mgene[data$generation==0], na.rm=T)), 
  	                           1-sqrt(1-mean(data$mgene[data$generation==0], na.rm=T)) )
  	else if(q>1 | q<0) stop("The allele frequency (q) should lie between 0 and 1.")

cat("Estimate allele frequency = ", q, "\n")
  	#G1=mutation status of proband
  	#G2=mutation status of sib
  	#G3=mutation status of child
  	#G4=mutation status of parent
  	#G5=mutation status of sib's child

if(mode=="dominant"){
# dominant  	
p1 <- (1+q-q^2)/(2-q) # P(G3=1|G1=1) = P(G4=1|G1=1)
p2 <- (4+5*q-6*q^2+q^3)/(8-4*q) # P(G2=1|G1=1)
p4 <- p1*p2+q*(1-p2) #P(G5=1|G1=1) = P(G5=1|G2=1)P(G2=1|G1=1) + P(G5=1|G2=0)P(G2=0|G1=1) 
p3 <- q-q^2/4 #P(G2=1|G1=0)
p0 <- 1-q+q^2/4 # P(G2=0|G1=0)
p5 <- p1*p3+q*(1-q)  #P(G5=1|G1=0) = P(G5=1|G2=1)P(G2=1|G1=0) + P(G5=1|G2=0)P(G2=0|G1=0)
}
else if(mode=="recessive"){

# recessive
p1 <- q # P(G3=1|G1=1) = P(G4=1|G1=1)
p2 <- (1+q)^2/4 # P(G2=1|G1=1)
p4 <- p1*p2+(1-q)*(1-p2) #P(G5=1|G1=1) = P(G5=1|G2=1)P(G2=1|G1=1) + P(G5=1|G2=0)P(G2=0|G1=1) 
p3 <- (3-2*q-q^2)/4 #P(G2=1|G1=0)
p0 <- (1+q)^2/4
p5 <- p1*p3+(1-q)*p0  #P(G5=1|G1=0) = P(G5=1|G2=1)P(G2=1|G1=0) + P(G5=1|G2=0)P(G2=0|G1=0)
 }
else stop("Unrecognized inheritance mode")


for(i in id.na){

i.rel <-data$relation[data$indID==i]
i.fam <- data$famID[data$indID==i]

faid <- data$fatherID[data$indID==i]
moid <- data$motherID[data$indID==i]

fag <- ifelse(faid==0, NA, data$mgene[data$indID==faid])
mog <- ifelse(moid==0, NA, data$mgene[data$indID==moid])

#if( is.parent(i, faid, moid))  p <- p1
#else if( is.sibling(faid, moid)) p <- p2
#else if( is.niece(i.rel, mut.ca)) p <- p4
#else if( is.grandson(i.rel, mut.ca)) p <- p6
#else if( is.cousin(i.rel, mut.ca)) p <- p6*(1-(1-q)^2) + p7*(1-q)^2
#else if(i.rel=="NBS") p <- q*(2-q)

# 1=self, 2=sib, 3=child, 4=parent, 5=sibchild, 6=hus,7=sibspouse

if(is.element(data$famID[data$indID==i], cfam.id)){ # carrier families (proband is carrier)
	if( i.rel==3|i.rel==4)  p <- p1
	else if( i.rel==2) p <- p2
	else if( i.rel==5){
 		if(sum(c(fag, mog), na.rm=TRUE)>0) p <- p1 # parent is a carrier
 		else if(is.na(fag) & is.na(mog)) p <- p4 # P(G5=1|G1=1) both parents' genotypes are unknown
 		else p <- q # both parents are not carriers.
		} 
	else if(i.rel==6|i.rel==7) p <- q*(2-q)
	carrp[data$indID==i] <- p
	}
else{ # non-carrier families (proband is not a carrier)
	if( i.rel==3)  p <- q #P(G3=1|G1=0)
	else if (i.rel==4 ) {
		fam1 <- data[data$famID==i.fam,]
		if(sum(fam1$mgene[fam1$generation==2], na.rm=TRUE)>0) p <- p1
		else p <- q
	}
	else if( i.rel==2) p <- p3
	else if( i.rel==5){
 		if(sum(c(fag, mog), na.rm=TRUE)>0) p <- p1 # parent is a carrier
 		else if(is.na(fag)&is.na(mog)) p <- p5 # P(G5=1|G1=0) both parents' genotypes are unknown
 		else p <- q # both parents are not carriers.
		} 
	else if(i.rel==6|i.rel==7) p <- q*(2-q)
	carrp[data$indID==i] <- p
	}

	
	}
} 
  else stop("Unrecognized method")

data$carrp <- carrp
return(data)
}