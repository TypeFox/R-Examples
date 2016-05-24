######### Hotelling's T^2 statistic using Shrinkage covariance matrix estimates ###############
######### x:data matrix; Row: sample; Column: variables(genes)                  ###############
######### y: vector defining two-group of the samples                           ###############


library(corpcor)

Hott2 <- function(x, y, var.equal=TRUE){
    if(!is.null(y)){
         cl <- as.factor(y)
         lab <- levels(cl)
         ind1 <- which(cl==lab[1])
         ind2 <- which(cl==lab[2])
     }
     if(is.null(ind1) | is.null(ind2))stop("Error: classes 1 and 2 are undefined.")

    data1 <- x[ind1,]
    data2 <- x[ind2,]
    nn <- dim(data1)[1]
    nd <- dim(data2)[1]
    SN <- cov.shrink(data1,verbose=FALSE, lambda.var=0)
    SD <- cov.shrink(data2,verbose=FALSE, lambda.var=0)
    xbar1 <- colMeans(data1)
    xbar2 <- colMeans(data2)
    xdiff <- xbar1 - xbar2

    if (var.equal){
       S <- ((nd - 1) * SD + (nn - 1) * SN)/(nd + nn - 2)
       t2 <- ((nd * nn)/(nd + nn)) * (xdiff %*% solve(S) %*% xdiff)
    }
    else{
       S <- SN/nn + SD/nd
       t2 <- xdiff %*% solve(S) %*% xdiff
    }

    return(as.numeric(t2))
}

######### Ordinary least square(OLS) statistic                                  ###############
######### x:data matrix; Row: sample; Column: variables(genes)                  ###############
######### y: vector defining two-group of the samples  
Tols <- function (x, y ) {
if(!is.null(y)){
         cl <- as.factor(y)
         lab <- levels(cl)
         ind1 <- which(cl==lab[1])
         ind2 <- which(cl==lab[2])
     }
     if(is.null(ind1) | is.null(ind2))stop("Error: classes 1 and 2 are undefined.")
    data1 <- x[ind1,]
    data2 <- x[ind2,]
    nn <- dim(data1)[1]
    nd <- dim(data2)[1]
    #-------------O'Brien's OLS----------------

  olsf <- function(x){
    t.test(x[ind1],x[ind2],var.equal=TRUE)$statistic
  }
 
  t.emat <- apply(x,2,olsf)
  cor.emat <- cor(x)
  tols <- sum(t.emat)/sqrt(sum(cor.emat))
   

    return(abs(as.numeric(tols)))
}




######### Wilks Lambda for n-group comparisons                           ###############
######### Y:data matrix; Row: sample; Column: variables(genes)               ###############
######### Class: vector defining the clinical outcome of the samples         ###############
######### type: type of contrast;                                            ###############
######### base: which group is considered the baseline group for Dunnett contrasts ###############
library(MASS)
library(multcomp)
"design.matrix" <-
function(factors)
{
### Keep the initial order
  n<-length(factors)
  fac<-factor(factors)
  lev<-levels(fac)
  l<-length(lev)
  if(l<2)
    stop("Should have at least two groups")
  X<-matrix(0,n,l)
  for(i in 1:l)
    X[factors==lev[i],i]<-1
  X
}



library(MASS)
ma.estimate <- function (Y, X) {
    ginv(t(X) %*% X) %*% t(X) %*% Y
}



Wilksn <- function(Y, class,type = c("Tukey", "Dunnett", "Sequence"), base=1){
   X <- design.matrix(class)
   hatB <- ma.estimate(Y, X)
   type <- match.arg(type)
   lab <- unique(class)

   E=0
   for(i in 1:length(lab)){
   ind0 <- which(class==lab[i])
   data0<- Y[ind0,]
   n0   <- length(ind0)
   S0   <- cov.shrink(data0,verbose=FALSE)
   E    <- E+(n0-1)*S0
   }

   k<- length(lab)
    names(Y) <- paste("T", 1:k, sep="")
    CM1 <- c()
    CM <- c()
    rnames <- c()
    if (!is.null(names(Y)))
        varnames <- names(Y)
   else varnames <- 1:length(Y)
    kindx <- 1:k
    for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
                CM1 <- rbind(CM1, as.numeric(kindx == j) - as.numeric(kindx == i))
                }}
   L <- CM1
   H <- t(hatB) %*% t(L) %*% ginv(L %*% ginv(t(X) %*% X) %*% t(L)) %*% L %*% hatB
   temp <- solve(E) %*%H
   stat <- 1/det(temp+diag(1,dim(temp)[1],dim(temp)[1]))            
    switch(type, Dunnett = {
        for (i in kindx[-base]) CM <- rbind(CM, as.numeric(kindx == i) - as.numeric(kindx == base))
        rnames <- paste(varnames[kindx[-base]], "-", varnames[base])
        },
                  Tukey = {
        for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
                CM <- rbind(CM, as.numeric(kindx == j) - as.numeric(kindx == i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i]))
            }}}, 
                  Sequence = {
        for (i in 2:k) {
            CM <- rbind(CM, as.numeric(kindx == i) - as.numeric(kindx == i - 1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i - 1]))
        }})


   L <- CM
   for (i in 1:nrow(L)){
   D  <- matrix(0,nrow(L),ncol(L))
   D[i,]=1
   L0  <- D*L
   H0 <- t(hatB) %*% t(L0) %*% ginv(L0 %*% ginv(t(X) %*% X) %*% t(L0)) %*% L0 %*% hatB
   temp0 <- solve(E) %*%H0
   stat0 <- 1/det(temp0+diag(1,dim(temp0)[1],dim(temp0)[1]))
   stat=c(stat, stat0)

   }
   names(stat)=c("ANOVA", rnames)
   return(stat)
}


library(MASS)
library(multcomp)
library(corpcor)

MAVTn <- function(DATA, GS, MCP=1 , alpha=0.01, nbPerm=5000){


 # DATA : expression data with rows=genes, columns=samples
 #        Note that the first row is the group of sample
 #
 # GS : gene sets
 #      -> a data matrix with rows=genes,
 #                        columns= gene sets,
 #                        GS[i,j]=1 if gene i in gene set j
 #                        GS[i,j]=0 otherwise
 #
 #
 # alpha: signinificant level
 #
 #
 # MCP: type of multiple comparison methods;
 #      Dunnett = 1, Tuckey = 2, Sequential Group = 3 
 # 
     cl <- as.numeric(DATA[1,])
     DATA <- DATA[-1,]
     genes <- rownames(DATA)
     k <- length(levels(factor(cl)))
     n.Samples  <- ncol(DATA)
     n.GeneSets <- ncol(GS)
     GeneSets.sizes <- apply(GS,2, function(z) sum(z==1))
     
     if (k>=3) {
     base <- 1
     if (MCP==1){type <- "Dunnett"}
     if (MCP==2){type <- "Tukey"}
     if (MCP==3){type <- "Sequence"}
     if ((1*(MCP==1)+1*(MCP==2)+1*(MCP==3))==0) stop("Error: MCP must be 1, 2 or 3")
     } 

     # observed statitic for each gene set
     if (k==2) {
      stat.ols.obs <- apply(GS, 2, function(z) Tols(t(DATA[which(z==1),]), cl))
      stat.hott.obs <- apply(GS, 2, function(z) Hott2(t(DATA[which(z==1),]), cl))
      }
     if (k>=3) stat.obs <- apply(GS, 2, function(z) Wilksn(t(DATA[which(z==1),]), cl, type , base))

     # stats obtained on 'permuted' data

     stat.ols.permut <- matrix(NA,nbPerm,n.GeneSets)
     stat.hott.permut <- matrix(NA,nbPerm,n.GeneSets)

     if (k>=3){
     stat.temp.permut <- matrix(NA,nbPerm*nrow(stat.obs),n.GeneSets)
     }
      for(i in 1:nbPerm) {
         ind <- sample(n.Samples)
         p.data <- DATA[,ind]
         if (k==2){
           stat.ols.permut[i,] <- apply(GS, 2, function(z) Tols(t(p.data[which(z==1),]), cl))
           stat.hott.permut[i,] <- apply(GS, 2, function(z) Hott2(t(p.data[which(z==1),]), cl))
           }

         if (k>=3){
         k.p <- apply(GS, 2, function(z) Wilksn(t(p.data[which(z==1),]), cl, type , base))
           stat.temp.permut[c(seq(1,nrow(stat.obs)))+nrow(stat.obs)*(i-1),] <- k.p
           }
         }



     if (k==2) {
        GeneSets.ols.pval <- apply(t(stat.ols.permut) >= stat.ols.obs, 1, sum)/nbPerm
        GeneSets.hott.pval <- apply(t(stat.hott.permut) >= stat.hott.obs, 1, sum)/nbPerm
     s.pvalue <- apply(DATA,1, function(z) unlist(summary(aov(z~as.factor(cl))))["Pr(>F)1"])
     nb.ols.sign <- which(GeneSets.ols.pval<=alpha)
     sg.ols.pvalue <- apply(as.matrix(GS[,nb.ols.sign]),2,function(z) c(list("GS size"=sum(z==1),"p-value"=round(s.pvalue[z==1],4))))
     nb.hott.sign <- which(GeneSets.hott.pval<=alpha)
     sg.hott.pvalue <- apply(as.matrix(GS[,nb.hott.sign]),2,function(z) c(list("GS size"=sum(z==1),"p-value"=round(s.pvalue[z==1],4))))
     fwe.ols.pvalue.permut <- matrix(NA,n.GeneSets,nbPerm)
     fwe.hott.pvalue.permut <- matrix(NA,n.GeneSets,nbPerm)
     
     for (i in 1:nbPerm){
         ind.p <- sample(nbPerm) 
         fwe.ols.pvalue.permut[,i] <- apply(t(stat.ols.permut)[,which(ind.p<=ceiling(nbPerm/10))] <= stat.ols.obs, 1, sum)/(nbPerm/10)
         fwe.hott.pvalue.permut[,i] <- apply(t(stat.hott.permut)[,which(ind.p<=ceiling(nbPerm/10))] <= stat.hott.obs, 1, sum)/(nbPerm/10)
        }
     q.ols.ind <- rank(GeneSets.ols.pval,ties.method ="first")
     q.hott.ind <- rank(GeneSets.hott.pval,ties.method ="first") 
     q.ols.val <- minp(fwe.ols.pvalue.permut,q.ols.ind,n.GeneSets,nbPerm)
     q.hott.val <- minp(fwe.hott.pvalue.permut,q.hott.ind,n.GeneSets,nbPerm)

     
     
#------------significant marker--------------------------------

     ols.star <- c()
     hott.star <- c()
     for (i in 1:n.GeneSets){
     if (GeneSets.ols.pval[i]<=alpha) ols.star[i] <- "*"
     else ols.star[i] <- " "

     if (GeneSets.hott.pval[i]<=alpha) hott.star[i] <-"*"
     else hott.star[i] <- " "
     }
     ols.sign <- data.frame(ols.star)
     colnames(ols.sign) <- " "
     hott.sign <- data.frame(hott.star)
     colnames(hott.sign) <- " "    
#-------------adjust p-value-----------------------------------
     ad.ols.pvalue <- round(p.adjust(GeneSets.ols.pval,method="BH"),4)
     fwe.ols.pvalue <- apply(q.ols.val<=GeneSets.ols.pval,1,sum)/nbPerm
     ad.hott.pvalue <- round(p.adjust(GeneSets.hott.pval,method="BH"),4)
     fwe.hott.pvalue <- apply(q.hott.val<=GeneSets.hott.pval,1,sum)/nbPerm
        res <- as.data.frame(cbind("GS size"              = GeneSets.sizes,
                                   "OLS p-value"          = GeneSets.ols.pval,
                                       ols.sign,
                                   "OLS adjusted p-value (FDR)"   = ad.ols.pvalue,
                                   "OLS adjusted p-value (FWE)"   = fwe.ols.pvalue,
                                   "T-square p-value"     = GeneSets.hott.pval,
                                        hott.sign,
                                   "T-square adjusted p-value (FDR)"= ad.hott.pvalue,
                                   "T-square adjusted p-value (FWE)"= fwe.hott.pvalue
                                        ))
     result <- list("p value"=res,"Singinificant gene set (one-sided)"=sg.ols.pvalue,"Singinificant gene set (two-sided)"=sg.hott.pvalue)
 #----------------gsaplot---------------------------------------   
     p.sort.ols <- sort(GeneSets.ols.pval)
     p.sort.hott <- sort(GeneSets.hott.pval)
     plot(p.sort.ols,type="l",xaxs = "i",  
        xlab = "Ranked p-value", ylab = "p-value", main = "")
     par(new=TRUE)
     plot(p.sort.hott,type="l",axes = FALSE, col="red", xaxs = "i",  
        xlab = "", ylab = "", main = "")    
     #lines( par()$usr[1:2], par()$usr[3:4] , lty=2)
     abline(a=0,b=1/n.GeneSets,lty=2)
     leg.names <- c("OLS",expression("Hotelling's "*T^2))
     legend(locator(1),leg.names,lty=c(1,1),col=c(1,2))
     }
     
     
     
     if (k>=3){
         stat.temp.pval <- matrix(NA,nrow(stat.obs),n.GeneSets)
         for (i in 1:nrow(stat.obs)){
         t.pval=apply(t(stat.temp.permut[seq(i,nrow(stat.temp.permut),nrow(stat.obs)),])<= stat.obs[i,],1,sum)/nbPerm
           stat.temp.pval[i,] <- t.pval
           }
           fwe.ad.pval.permut <- matrix(NA,n.GeneSets,nbPerm)
           stat.permut.obs <- t(stat.temp.permut[seq(1,nrow(stat.temp.permut),nrow(stat.obs)),])
     for (i in 1:nbPerm){
         ind.p <- sample(nbPerm) 
         fwe.ad.pval.permut[,i] <- apply(stat.permut.obs[,which(ind.p<=ceiling(nbPerm/10))] <= stat.obs[1,], 1, sum)/(nbPerm/10)
        }
     q.mon.ind <- rank(t(stat.temp.pval)[,1],ties.method ="first") 
     q.ad.val <- minp(fwe.ad.pval.permut,q.mon.ind,n.GeneSets,nbPerm)
           
     r1<-(1:(k*(k-1)/2))
     r2<-NULL
     r3<-NULL
     for(i in 1:k){
     r2<-c(r2,rep(i,k-i))
     if (i+1<=k)
     r3<-c(r3,((i+1):k))
     }
     names=paste(rownames(k.p)," p-value",sep="")
     dimnames(stat.temp.pval) <- list(names)
     s.pvalue <- apply(DATA,1, function(z) unlist(summary(aov(z~as.factor(cl))))["Pr(>F)1"])
     nb.sign <- which(stat.temp.pval[1,]<=alpha)
     sg.pvalue <- apply(as.matrix(GS[,nb.sign]),2,function(z) c(list("GS size"=sum(z==1),"p-value"=round(s.pvalue[z==1],4))))
     
     star <- c()
     for (i in 1:n.GeneSets){
     if (stat.temp.pval[1,i]<=alpha) star[i] <- "*"
     else star[i] <- " "
     }
     sign <- data.frame(star)
     ad.pval <- round(p.adjust(t(stat.temp.pval)[,1],method="BH"),4)
     fwe.ad.pval <- 1- apply(q.ad.val<=t(stat.temp.pval)[,1],1,sum)/nbPerm
     colnames(sign) <- " "
     res <- as.data.frame(cbind("GS size"              = GeneSets.sizes,
                                "MANOVA p-value"   =t(stat.temp.pval)[,1],
                                sign,
                                "adjusted p-value (FDR)"   = ad.pval,
                                "adjusted p-value (FWE)"   = fwe.ad.pval,
                                t(stat.temp.pval)[,-1]     ))
   
     rownames(res)=colnames(GS)
     result <- list("p value"=res, "singinificant gene set"=sg.pvalue)
     p.sort.manova <- sort(t(stat.temp.pval)[,1])

     plot(p.sort.manova,type="l",xaxs = "i", 
        xlab = "Ranked p-value", ylab = "MANOVA p-value", main = "")
     #lines( par()$usr[1:2], par()$usr[3:4] , lty=2)
     abline(0,1/n.GeneSets,lty=2)
     }


     
     
     
   
    return(result)
}

#------------gstplot-------------------------
GSTplot <- function (data,gs,geneset.name=NULL,alpha=0.01){
     ccl <- as.numeric(data[1,])
     k <- length(levels(factor(ccl)))
     group <- c(1,2)
if (k==2){

if (is.null(geneset.name)||sum(1*(colnames(gs)==geneset.name))==0) {
     stop("Gene set name is unmatched", 
     call. = FALSE)
    }
    GS=gs[,which(colnames(gs)==geneset.name)]
    
    if (!is.null(data)){
         cl <- as.factor(group)
         lab <- levels(cl)
         ind1 <- which(data[1,]==lab[1])
         ind2 <- which(data[1,]==lab[2])
     }
    g=as.matrix(2:dim(data)[[1]])
    t.stat.two <- apply(g, 1, function(z) t.test(data[z,ind2],data[z,ind1],)$statistic )
    names(t.stat.two)=row.names(data[-1,])
    t.stat.sorted <- sort(t.stat.two)
    t.stat.order <- order(t.stat.two)
    gene.names.sorted <- names(t.stat.sorted)
    
aa=rt(length(t.stat.sorted),dim(data)[[2]]-2)   
m <- dim(data[-1,])[[1]]
ings <- which(GS==1)
stopa1 <- numeric(length(ings))
for (i in 1:length(ings)){
stopa1[i] <- which(t.stat.order==ings[i])
}
stopa <- sort(stopa1)
stopb <- stopa
    for (i in 2:(length(stopa)-1)){
    if (stopa[i]<=m/2) {
        if  (stopa[i]-max(stopa[i-1],stopb[i-1])<130)
        stopb[i] <- max(stopa[i],stopb[i-1])+130
        else stopb[i] <- stopa[i]
        }
    if (stopa[i]>m/2)   {
        if  (min(stopa[i+1],(stopb[i+1]-260)-stopa[i])<130) 
        stopb[i] <- min(stopa[i],stopb[i+1]-260)-130
        else stopb[i] <- stopa[i]
        }
    }
stopb[length(stopa)] <- stopa[length(stopa)]-130    
    nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), c(4, 4), 
        c(1, 3), TRUE)
    par(adj = 0.5, plt = c(0.2, 0.9, 0, 0.8))
    plot(1:m, rep(0,m), ylim = c(0.04,1), xaxs = "i", axes = FALSE, 
        xlab ="", ylab ="", main = "", col = 0)
        
    for (i in 1:length(stopa))  {
text(stopb[i]-200,0.3,gene.names.sorted[stopa[i]],srt = 90, cex = 0.5, font = 1, adj=-1, pos=4)     
lines(c(stopa[i],stopa[i],stopb[i]),c(0,0.15,0.25),type="l",col=1)
}

    par(adj = 0.5, plt = c(0.2, 0.9, 0.65, 1))
    plot(c(0, m), c(-6, 6), xlim = c(1, dim(data[-1])[[1]]), xaxs = "i", yaxs = "i", 
        xlab = "Ranked t-statistics", ylab = "Observed t-statistics", main = "", col = 0)
        text(m/25, 3.87, geneset.name, cex = 0.8, font = 2, adj = 0)
        polygon(c(0, m*alpha, m*alpha, 0), 
            c(-6, -6, 10, 10), border = FALSE, col = 8)
        polygon(c( m*(1-alpha),m, m,m*(1-alpha)), 
            c(-6, -6, 10, 10), border = FALSE, col = 8)

        lines(sort(t.stat.sorted),type="l",col=2)
        for (i in 1:length(stopa)){
        lines(c(stopa[i],stopa[i]),c(0,sort(t.stat.two[which(GS==1)])[i]),type="l",col=2)
        }
        
        abline(h=0,lty=2)
        #lines( par()$usr[1:2], par()$usr[3:4] , lty=2) 
        
par(new=TRUE)
Fnt <- ecdf(sort(t.stat.two[which(GS==1)]))
ecdf.y <- Fnt(knots(Fnt))
plot(stopa,ecdf.y, type="s",xlim = c(1, dim(data[-1])[[1]]),xaxs = "i",axes=FALSE,xlab="",ylab="")

}


#---------------plot for MANOVA---------------------------------------

if (k>=3)
{
if (is.null(geneset.name)||sum(1*(colnames(gs)==geneset.name))==0) {
     stop("Gene set name is unmatched", 
     call. = FALSE)
    }
    GS=gs[,which(colnames(gs)==geneset.name)]
    cl <- as.numeric(data[1,])
    df1 <- length(levels(factor(cl)))-1
    df2 <- length(cl)-length(levels(factor(cl)))    
    g=as.matrix(2:dim(data)[[1]])
    F.sata <- apply(g, 1, function(z) unlist(summary(aov(data[z,]~as.factor(cl))))["F value1"] )
    names(F.sata)=row.names(data[-1,])
    t.stat.sorted <- sort(F.sata)
    t.stat.order <- order(F.sata)
    gene.names.sorted <- names(t.stat.sorted)
    lables.sig <- qf(1-alpha,df1,df2)
    lowb <- (dim(data)[1]-1)*(1-alpha)
aa=rt(length(t.stat.sorted),dim(data)[[2]]-2)   
m <- dim(data[-1,])[[1]]
maxf <- max(F.sata)
ings <- which(GS==1)
stopa1 <- numeric(length(ings))
for (i in 1:length(ings)){
stopa1[i] <- which(t.stat.order==ings[i])
}
stopa <- sort(stopa1)
stopb.big <- sort(which(stopa>=lowb),decreasing=TRUE)
stopb.i <- max(stopa)-0.01*m
if (length(stopb.big)>1) {
    for (i in 2:length(stopb.big)){
        if  (stopb.i[i-1]-stopa[stopb.big[i]]>m*0.013) 
        stopb.i <- c(stopb.i,stopa[stopb.big[i]])
        else stopb.i <- c(stopb.i,stopb.i[i-1]-m*0.013)
    }
    stopb <- sort(stopb.i)}
    
     nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), c(4, 4), 
        c(1, 3), TRUE)
    par(adj = 0.5, plt = c(0.2, 0.9, 0, 0.8))
    plot(1:m, rep(0,m), ylim = c(0.04,1), xaxs = "i", axes = FALSE, 
        xlab ="", ylab ="", main = "", col = 0)
        
    for (i in 1:length(stopa))  {
    
lines(c(stopa[i],stopa[i]),c(0,0.15),type="l",col=1)
}
for (i in 1:length(stopb.i))    {

lines(c(sort(stopa,decreasing=TRUE)[i],stopb.i[i]),c(0.15,0.25),type="l",col=1)
}


name.sig <- which(stopa>m*(1-alpha))
for (i in 1:length(stopb.i)){
text(stopb.i[i]-m*0.026,0.27,gene.names.sorted[sort(stopa,decreasing=TRUE)[i]],srt = 90, cex = 0.5, font = 2, adj=-1, pos=4)    
}

    par(adj = 0.5, plt = c(0.2, 0.9, 0.65, 1))
    plot(c(0, m), c(0, maxf), xlim = c(1, m), xaxs = "i", yaxs = "i", 
        xlab = "Ranked F-statistics", ylab = "Observed F-statistics", main = "", col = 0)
        text(m/25, maxf-5, geneset.name, cex = 0.8, font = 2, adj = 0)
        polygon(c( lowb,m, m,lowb), 
            c(0, 0, maxf, maxf), border = FALSE, col = 8)

        lines(sort(t.stat.sorted),type="l",col=2)
        for (i in 1:length(stopa)){
        lines(c(stopa[i],stopa[i]),c(0,sort(F.sata[which(GS==1)])[i]),type="l",col=2)
        }
        
        abline(h=0,lty=2)
        #lines( par()$usr[1:2], par()$usr[3:4] , lty=2) 
        
par(new=TRUE)
Fnt <- ecdf(sort(F.sata[which(GS==1)]))
ecdf.y <- Fnt(knots(Fnt))
plot(c(0,stopa),c(0,ecdf.y), type="s",xlim = c(1, m),xaxs = "i",axes=FALSE,xlab="",ylab="")

}
}

#--------------minp----------------
minp <- function (p,rank,n.GeneSets,nbPerm){
     q.val <- matrix(NA,n.GeneSets,nbPerm)
     q.val.ra <- matrix(NA,n.GeneSets,nbPerm)
     q.val.ra[n.GeneSets,] <- p[which.max(rank),]
     for (j in (n.GeneSets-1):1){
         q.val.ra[j,] <- pmin(q.val.ra[(j+1),],p[which(rank==j),])
        }
     for (i in n.GeneSets:1) {
         q.val[which(rank==i),] <- q.val.ra[i,] 
        }
     return(q.val)
    }
    
#library(MCMCpack)
library(randomForest)
library(foreach)
library(corpcor)
#library(doParallel)


############### Random Forest ################
MAVTp <- function(DATA, GS, nbPerm=5000, numoftree=500,type=c("cont","cate"),impt=TRUE){

     cl <- as.numeric(DATA[1,])
     DATA_1 <- DATA[-1,]
     p.value=rep(0,dim(GS)[2])
     temp <- "importance measure"
     for (i in 1:dim(GS)[2]){
     ind_gs <- which(GS[,i]==1)
     DATA <- DATA_1[ind_gs,]

if (type=="cont"){
################  continuous group  ################


group <- as.numeric(cl)
nb.Samples  <- ncol(DATA)
test.rf <- randomForest(t(DATA),group, ntree=numoftree, importance=TRUE, proximity=FALSE, mtry=5)
value=round(test.rf$mse[numoftree],4)



nullV<-foreach(j=1:nbPerm, .combine=cbind) %do%{
ind=sample(1:length(group));
testp.rf <- randomForest(t(DATA)[ind,], group, ntree=numoftree, importance=FALSE, proximity=FALSE,mtry=5);
round(testp.rf$mse[numoftree],4)
}

Pvalue.rf=sum(nullV<=value)/nbPerm
#RF_pvalue=Pvalue.rf
p.value[i]<- Pvalue.rf
  temp <- list(temp, Imp=importance(test.rf))
}
################  categorical data  ################
if (type=="cate"){

group <- as.factor(cl)
  test.rf <- randomForest(t(DATA),group, ntree=numoftree, importance=TRUE, proximity=FALSE)
  value=round(test.rf$err.rate[numoftree,1],4)

  nullV<-foreach(j=1:nbPerm, .combine=cbind) %do%{
    ind=sample(1:length(group))
    testp.rf <- randomForest(t(DATA)[ind,], group, ntree=numoftree, importance=FALSE, proximity=FALSE);
    round(testp.rf$err.rate[numoftree,1],4)
  }

  Pvalue=sum(nullV<=value)/length(nullV)
  #RF_pvalue=Pvalue
  p.value[i]<- Pvalue
  temp <- list(temp, Imp=importance(test.rf)[,3:4])
}



}


if (impt==TRUE) result <- list("p.value"=p.value, "importance"=temp)
else result <- list(p.value)
return(result)
}

