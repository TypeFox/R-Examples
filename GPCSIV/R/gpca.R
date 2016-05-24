gpca <-
function (xmin = list, xmax = list, reduire = 0, nomVar = NULL, 
    axes = c(1, 2), axes2=c(1,2,3),nomInd = NULL, legend = NULL, xlim = NULL, 
    ylim = NULL, nametable = NULL, plot3d.table=NULL) 
{
Centrage<-function(x,reduire=0){

x1<-scale(x,center=TRUE, scale=FALSE);

x2<-scale(x,center=TRUE, scale=TRUE);

if (reduire==1) { y=x1}
else 
if (reduire==0)  {y=x2};
 return(as.matrix(y));
}
    Build3d<-function(dimens,GenPC,p,axes2,n,nomInd,axes,lab.x,lab.y,nametable)
{        

dev.new()
k=dimens;
       
        PCinterval1 <- fintab(as.data.frame(GenPC[[k]]))
         tempGenPC<-GenPC[[k]]         

regad<-ucad(p)
p2=nrow(regad)

for(ii in 1:p2)
{
colnames( PCinterval1)[ii]<-paste(regad[ii,1])
}

  united=seq(1,2*p,2)
 res<-matrix(0,p,2)
 for(i in 1:p)
 {
 res[i,]<-c(united[i], united[i]+1)
 }

 zPC1<-vector('list',p)
 for(i in 1:p)
 {
 zPC1[[i]]<-PCinterval1[,res[i,]]
 }

 pc1<-zPC1[[axes2[1]]]
 pc2<-zPC1[[axes2[2]]]
 pc3<-zPC1[[axes2[3]]]


        xmin1 <- min(pc1[, 1])
        xmin2 <- min(pc1[, 2])
        xxmin <- min(xmin1, xmin2)
        xmax1 <- max(pc1[, 1])
        xmax2 <- max(pc1[, 2])
        xxmax <- max(xmax1, xmax2)
        ymin1 <- min(pc2[, 1])
        ymin2 <- min(pc2[, 2])
        yymin <- min(ymin1, ymin2)
        ymax1 <- max(pc2[, 1])
        ymax2 <- max(pc2[, 2])
        yymax <- max(ymax1, ymax2)
        zmin1 <- min(pc3[, 1])
        zmin2 <- min(pc3[, 2])
        zzmin <- min(zmin1, zmin2)
        zmax1 <- max(pc3[, 1])
        zmax2 <- max(pc3[, 2])
        zzmax <- max(zmax1, zmax2)
        picto <- scatterplot3d(c(xxmin, xxmax), c(yymin, yymax), c(zzmin, 
            zzmax), xlab = paste('PC',axes2[1],sep='.'), ylab =paste('PC',axes2[2],sep='.'), 
            zlab = paste('PC',axes2[3],sep='.'), type = "n",main='3d-Projection of table onto  principal axes of the compromise' ,sub=nametable[k])
cube <- rbind(c(1, 1, 1), c(2, 1, 1), c(2, 1, 2), c(1, 
            1, 2), c(1, 1, 1), c(1, 2, 1), c(1, 2, 2), c(2, 2, 
            2), c(2, 2, 1), c(1, 2, 1), c(1, 2, 2), c(1, 1, 2), 
            c(2, 1, 2), c(2, 2, 2), c(2, 2, 1), c(2, 1, 1))

for (i in 1:n) {
            vec.x <- pc1[i, cube[, 1]]
            vec.y <- pc2[i, cube[, 2]]
            vec.z <- pc3[i, cube[, 3]]
            picto$points3d(vec.x, vec.y, vec.z, type = "l", lty = 1, 
                col = i)
}


ltext <- nomInd
            textPoints <- cbind((pc1[,1] + pc1[,2])/2, (pc2[,1] + pc2[,2])/2, (pc3[,1] + pc3[,2])/2)
            text(picto$xyz.convert(textPoints), labels = ltext,col=1:n)
      

dev.new() 

xlimbis=range(tempGenPC)
ylimbis=range(tempGenPC)

 n1 <- 1
        n2 <- n
        m1 <- n2 + 1
        m2 <- 2 * n2
        plot(tempGenPC[, axes[1]], tempGenPC[, axes[2]], xlab = lab.x, 
            ylab = lab.y, asp = 1, main = "Projection of table onto axes of the compromise.", xlim = xlimbis, 
            ylim = xlimbis,sub=nametable[k])
        text(tempGenPC[1:n2, axes[1]], tempGenPC[1:n2, axes[2]], 
            labels = nomInd, pos = 4,col=1:n)
        abline(h = 0, col = "blue")
        abline(v = 0, col = "blue")
        rect(tempGenPC[n1:n2, axes[1]],tempGenPC[n1:n2, axes[2]], 
            tempGenPC[m1:m2, axes[1]], tempGenPC[m1:m2, axes[2]], 
            border = 1:n)
        
}

#####  fonction restoratrice des min max
fintab<-function(gen1)
{
k=nrow(gen1)
k3=k/2
k1=2*ncol(gen1)
k2=k1/2
t5=vector('list',k2)
t1=gen1[1:(k/2),]
t1=data.frame(t1)
t2=gen1[((k/2)+1):k,]
t2=data.frame(t2)
t3=cbind(t1,t2)
t3=data.frame(t3)
for (i in 1:k2)
t5[[i]]=cbind(t3[,i],t3[,i+k2])

t6=vector('list',k2)
for (i in 1:k2)
#t6[[i]]= matrix(nrow=k3,ncol=2)
t6[[i]]= matrix(0,nrow=k3,ncol=2)

for (j in 1:k3)
{
for (i in 1:k2)
t6[[i]][j,]= t(as.data.frame(range((t5[[i]])[j,])))
}

t7=t6[[1]]
{
for (s in 2:k2)
t7=cbind(t7,t6[[s]])
}
t7=data.frame(t7)
return(t7)
}

ucad<-function(p)
{
res<-vector('list',p)
auxi<-NULL
for(i in 1:p)
{
res[[i]]<-(paste(c('PCMin','PCmax'),i, sep='.')  )
auxi<-c(auxi,res[[i]])

}

auxi<-(as.data.frame(auxi))

auxi
}

    xmin <- as.list(xmin)
    d <- length(xmin)
    for (i in 1:d) {
        xmin[[i]] <- as.data.frame(xmin[[i]])
    }
    for (i in 1:d) {
        xmax[[i]] <- as.data.frame(xmax[[i]])
    }

    if (d > 1) {
        xmax <- as.list(xmax)
        for (i in 1:d) {
            xmin[[i]] <- Centrage(xmin[[i]], reduire = reduire)
        }
        for (i in 1:d) {
            xmax[[i]] <- Centrage(xmax[[i]], reduire = reduire)
        }
        cpmin <- xmin
        for (i in 1:d) {
            cpmin[[i]] <- as.matrix(cpmin[[i]])
        }
        for (i in 2:d) {
            cpmin[[i]] <- cpmin[[i - 1]] + cpmin[[i]]
        }
        resmin <- cpmin[[d]]/d
        cpmax <- xmax
        for (i in 1:d) {
            cpmax[[i]] <- as.matrix(cpmax[[i]])
        }
        for (i in 2:d) {
            cpmax[[i]] <- cpmax[[i - 1]] + cpmax[[i]]
        }
        resmax <- cpmax[[d]]/d
        X <- vector("list", d)
        for (i in 1:d) {
            X[[i]] <- as.matrix(rbind(xmin[[i]], xmax[[i]]))
        }
        resmin <- as.matrix(resmin)
        resmax <- as.matrix(resmax)
        cpres <- rbind(resmin, resmax)
        n <- nrow(xmin[[1]])
        p <- ncol(cpres)
        V <- (1/n) * t(cpres) %*% cpres
        Vect <- svd(cpres, nu = n, nv = p)$v
        colnames(Vect)<-paste('Eigen Vector',1:p,sep='.')
        Val <- svd(V, nu = p)$d
        Pval <- Val/(sum(Val))
        Pval <- 100 * Pval
        Pval <- as.matrix(Pval)
        rownames(Pval) <- paste("Component n.", 1:p, "")
        colnames(Pval) <- c("% of variability")
        pval2 <- as.data.frame(matrix(NA, length(Val), 3))
        rownames(pval2) <- paste("comp", 1:length(Val))
        colnames(pval2) <- c("eigenvalue", "percentage of variance", 
            "cumulative percentage of variance")
        pval2[, "eigenvalue"] <- Val
        pval2[, "percentage of variance"] <- (Val/sum(Val)) * 
            100
        pval2[, "cumulative percentage of variance"] <- cumsum(pval2[, 
            "percentage of variance"])
        cp1 <- round(Pval[axes[1]], digits = 2)
        cp2 <- round(Pval[axes[2]], digits = 2)
        lab.x <- paste("component n.", axes[1], " (", cp1, "%)", 
            sep = "")
        lab.y <- paste("component n.", axes[2], " (", cp2, "%)", 
            sep = "")
        PC <- cpres %*% Vect
        Correl <- cor(cpres, PC)
        if (is.null(nomVar)) {
            nomVar <- paste("Variable ", 1:p, sep = "")
        }

         dev.new()
         n1 <- 1
        n2 <- nrow(xmin[[1]])
        m1 <- n2 + 1
        m2 <- 2 * n2
        lab.x <- paste(" Dim n.", axes[1], " (", cp1, "%)", sep = "")
        lab.y <- paste(" Dim n.", axes[2], " (", cp2, "%)", sep = "")
        if (is.null(nomInd)) {
            nomInd <- paste("Row ", 1:n2, sep = "")
        }
        if (is.null(xlim)) {
            xlim <- range(PC)
        }
        if (is.null(ylim)) {
            ylim <- range(PC)
        }
        if (is.null(legend)) {
            legend[1] <- max(xlim) + 0.03
            legend[2] <- max(ylim) + 0.03
        }
        xlegend <- legend[1]
        ylegend <- legend[2]
        plot(PC[, axes[1]], PC[, axes[2]], xlab = lab.x, 
            ylab = lab.y, asp = 1, main = "Compromise Factorial Map", xlim = xlim, 
            ylim = ylim)
        text(PC[1:n2, axes[1]], PC[1:n2, axes[2]], 
            labels = nomInd, pos = 4,col=1:n)
        abline(h = 0, col = "blue")
        abline(v = 0, col = "blue")
        rect(PC[n1:n2, axes[1]],PC[n1:n2, axes[2]], 
            PC[m1:m2, axes[1]], PC[m1:m2, axes[2]], 
            border = 1:n)
       
         dev.new()

       # PCinterval0 <- vector("list", d)
        PCinterval0 <- fintab(as.data.frame(PC))
 
regad<-ucad(p)
p2=nrow(regad)

for(ii in 1:p2)
{
colnames( PCinterval0)[ii]<-paste(regad[ii,1])
}

  united=seq(1,2*p,2)
 res<-matrix(0,p,2)
 for(i in 1:p)
 {
 res[i,]<-c(united[i], united[i]+1)
 }



 zPC0<-vector('list',p)
 for(i in 1:p)
 {
 zPC0[[i]]<-PCinterval0[,res[i,]]
 }




 pc1<-zPC0[[axes2[1]]]
 pc2<-zPC0[[axes2[2]]]
 pc3<-zPC0[[axes2[3]]]




        xmin1 <- min(pc1[, 1])
        xmin2 <- min(pc1[, 2])
        xxmin <- min(xmin1, xmin2)
        xmax1 <- max(pc1[, 1])
        xmax2 <- max(pc1[, 2])
        xxmax <- max(xmax1, xmax2)
        ymin1 <- min(pc2[, 1])
        ymin2 <- min(pc2[, 2])
        yymin <- min(ymin1, ymin2)
        ymax1 <- max(pc2[, 1])
        ymax2 <- max(pc2[, 2])
        yymax <- max(ymax1, ymax2)
        zmin1 <- min(pc3[, 1])
        zmin2 <- min(pc3[, 2])
        zzmin <- min(zmin1, zmin2)
        zmax1 <- max(pc3[, 1])
        zmax2 <- max(pc3[, 2])
        zzmax <- max(zmax1, zmax2)
        picto <- scatterplot3d(c(xxmin, xxmax), c(yymin, yymax), c(zzmin, 
            zzmax), xlab = paste('PC',axes2[1],sep='.'), ylab =paste('PC',axes2[2],sep='.'), 
            zlab = paste('PC',axes2[3],sep='.'), type = "n",main='Scatterplot3d of the compromise')

cube <- rbind(c(1, 1, 1), c(2, 1, 1), c(2, 1, 2), c(1, 
            1, 2), c(1, 1, 1), c(1, 2, 1), c(1, 2, 2), c(2, 2, 
            2), c(2, 2, 1), c(1, 2, 1), c(1, 2, 2), c(1, 1, 2), 
            c(2, 1, 2), c(2, 2, 2), c(2, 2, 1), c(2, 1, 1))
for (i in 1:n) {
            vec.x <- pc1[i, cube[, 1]]
            vec.y <- pc2[i, cube[, 2]]
            vec.z <- pc3[i, cube[, 3]]
            picto$points3d(vec.x, vec.y, vec.z, type = "l", lty = 1, 
                col = i)
}



ltext <- nomInd
            textPoints <- cbind((pc1[,1] + pc1[,2])/2, (pc2[,1] + pc2[,2])/2, (pc3[,1] + pc3[,2])/2)
            text(picto$xyz.convert(textPoints), labels = ltext,col=1:n)
 
         dev.new()

        rownames(Correl) <- nomVar
        plot(Correl[, axes[1]], Correl[, axes[2]], xlab = lab.x, 
            ylab = lab.y, xlim = c(-1.5, 1.5), ylim = c(-1.5, 
                1.5), asp = 1, main = "Correlation Map")
        text(Correl[, axes[1]], Correl[, axes[2]], labels = rownames(Correl), 
            pos = 4)
        segments(0, 0, Correl[, axes[1]], Correl[, axes[2]], 
            col = 1:p)
        abline(h = 0, col = "blue")
        abline(v = 0, col = "blue")
        GenPC <- vector("list", d)
        PCinterval <- vector("list", d)
        MeanTable1 <- NULL
        MeanTable2 <- NULL
        MeanTable <- NULL
        for (i in 1:d) {
            GenPC[[i]] <- as.matrix(X[[i]] %*% Vect)
            MeanTable1[[i]] <- (t(as.matrix(apply(xmin[[i]], 
                2, mean)))) %*% as.matrix(Vect)
            MeanTable2[[i]] <- (t(as.matrix(apply(xmax[[i]], 
                2, mean)))) %*% as.matrix(Vect)
            MeanTable[[i]] <- rbind(MeanTable1[[i]], MeanTable2[[i]])
            rownames(MeanTable[[i]]) <- paste("Dim", c(axes[1], 
                axes[2]))
            colnames(MeanTable[[i]]) <- nomVar
        }
        for (i in 1:d) {
            PCinterval[[i]] <- fintab(as.data.frame(GenPC[[i]]))
##
regad<-ucad(p)
p2=nrow(regad)
for(ii in 1:p2)
{
colnames( PCinterval[[i]])[ii]<-paste(regad[ii,1])
}
##
        }
        dev.new()
        n1 <- 1
        n2 <- nrow(xmin[[1]])
        m1 <- n2 + 1
        m2 <- 2 * n2
        lab.x <- paste(" Dim n.", axes[1], " (", cp1, "%)", sep = "")
        lab.y <- paste(" Dim n.", axes[2], " (", cp2, "%)", sep = "")
        if (is.null(nomInd)) {
            nomInd <- paste("Row ", 1:n2, sep = "")
        }
        if (is.null(xlim)) {
            xlim <- range(GenPC)
        }
        if (is.null(ylim)) {
            ylim <- range(GenPC)
        }
        if (is.null(legend)) {
            legend[1] <- max(xlim) + 0.03
            legend[2] <- max(ylim) + 0.03
        }
        xlegend <- legend[1]
        ylegend <- legend[2]
        plot(GenPC[[1]][, axes[1]], GenPC[[1]][, axes[2]], xlab = lab.x, 
            ylab = lab.y, asp = 1, main = "Factorial Map", xlim = xlim, 
            ylim = ylim)
        text(GenPC[[1]][1:n2, axes[1]], GenPC[[1]][1:n2, axes[2]], 
            labels = nomInd, pos = 4)
        abline(h = 0, col = "blue")
        abline(v = 0, col = "blue")
        rect(GenPC[[1]][n1:n2, axes[1]], GenPC[[1]][n1:n2, axes[2]], 
            GenPC[[1]][m1:m2, axes[1]], GenPC[[1]][m1:m2, axes[2]], 
            border = 1)
        for (i in 2:d) {
            points(GenPC[[i]][, axes[1]], GenPC[[i]][, axes[2]], 
                col = i, type = "p")
            text(GenPC[[i]][1:n2, axes[1]], GenPC[[i]][1:n2, 
                axes[2]], col = i, labels = nomInd)
            rect(GenPC[[i]][n1:n2, axes[1]], GenPC[[i]][n1:n2, 
                axes[2]], GenPC[[i]][m1:m2, axes[1]], GenPC[[i]][m1:m2, 
                axes[2]], border = i)
        }
        if (is.null(nametable)) {
            nametable <- paste("table n.", 1:d)
        }
        legend(xlegend, ylegend, nametable, col = c(1:d), 
            text.col = 1:d, lty = 1:d, pch = 1:d, merge = TRUE, 
            bg = "gray90")
        MeanRange = range(MeanTable)
        minRange = min(MeanRange)
        maxRange = max(MeanRange)
        dev.new()
        nametable2<-nametable
        if (is.null(nametable2)) {
            nametable2 <- paste("table n.", 1:d)
        }
        plot(MeanTable1[[1]][, axes[1]], MeanTable1[[1]][, axes[2]], 
            xlim = c(minRange, maxRange), ylim = c(minRange, 
                maxRange), xlab = lab.x, ylab = lab.y)
        points(MeanTable2[[1]][, axes[1]], MeanTable2[[1]][, 
            axes[2]])
        rect(MeanTable1[[1]][, axes[1]], MeanTable1[[1]][, axes[2]], 
            MeanTable2[[1]][, axes[1]], MeanTable2[[1]][, axes[2]], 
            border = 1)
        text(MeanTable1[[1]][, axes[1]], MeanTable1[[1]][, axes[2]], 
            col = 1, labels = nametable2[1])
        for (i in 2:d) {
            points(MeanTable1[[i]][, axes[1]], MeanTable1[[i]][, 
                axes[2]], col = i)
            points(MeanTable2[[i]][, axes[1]], MeanTable2[[i]][, 
                axes[2]], col = i)
            text(MeanTable2[[i]][, axes[1]], MeanTable2[[i]][, 
                axes[2]], col = i, labels = nametable2[i])
            rect(MeanTable1[[i]][, axes[1]], MeanTable1[[i]][, 
                axes[2]], MeanTable2[[i]][, axes[1]], MeanTable2[[i]][, 
                axes[2]], border = i)
        }
        legend(xlegend, ylegend, nametable2, col = c(1:d), text.col = 1:d, 
            lty = 1:d, pch = 1:d, merge = TRUE, bg = "gray90")
        colnames(Correl) <- paste("Component", 1:p, sep = " ")
        colnames(PC) <- paste("Component", 1:p, sep = " ")
        for (i in 1:d) {
            rownames(GenPC[[i]]) <- c(nomInd, paste(1:n))
            colnames(GenPC[[i]]) <- nomVar
        }
        write.table(cbind(nomVar, Correl), file = "CorrelInter.txt", 
            sep = " ", row.names = FALSE)
        bsun <- ".txt"
        for (i in 1:d) {
            write.table(PCinterval[i], paste(bsun, i, sep = "."), 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
        }

      
#######################################################################################
####

if (length(plot3d.table)>0){
for(m in 1:length(plot3d.table))
{ 
Build3d(dimens=plot3d.table[m],GenPC=GenPC,p=p,axes2=axes2,n=n,nomInd=nomInd,axes=axes,lab.x,lab.y,nametable=nametable)
}

}
#####################################################################################

        return(list(Correl = Correl, 
            Vect = Vect, Pval = pval2, PCinterval = PCinterval,PCCompromise=PCinterval0))
    }
    else {
        zmax <- as.list(xmax)
        zmin <- as.list(xmin)
        hmax <- zmax[[1]]
        hmin <- zmin[[1]]
        hmax <- as.matrix(hmax)
        hmin <- as.matrix(hmin)
        vmin <- Centrage(hmin, reduire = reduire)
        vmax <- Centrage(hmax, reduire = reduire)
        X <- as.matrix(rbind(vmin, vmax))
        cpres <- X
        n <- nrow(xmin[[1]])
        p <- ncol(cpres)
        V <- (1/n) * t(cpres) %*% cpres
        Vect <- svd(V, nu = n, nv = p)$v
        colnames(Vect)<-paste('Eigen Vector',1:p,sep='.')
        Val <- svd(V, nu = p)$d
        Pval <- Val/(sum(Val))
        Pval <- 100 * Pval
        Pval <- as.matrix(Pval)
        rownames(Pval) <- paste("component n.", 1:p, "")
        colnames(Pval) <- c("% of variability")
        pval2 <- as.data.frame(matrix(NA, length(Val), 3))
        rownames(pval2) <- paste("comp", 1:length(Val))
        colnames(pval2) <- c("eigenvalue", "percentage of variance", 
            "cumulative percentage of variance")
        pval2[, "eigenvalue"] <- Val
        pval2[, "percentage of variance"] <- (Val/sum(Val)) * 
            100
        pval2[, "cumulative percentage of variance"] <- cumsum(pval2[, 
            "percentage of variance"])
        cp1 <- round(Pval[axes[1]], digits = 2)
        cp2 <- round(Pval[axes[2]], digits = 2)
        lab.x <- paste("component n.", axes[1], " (", cp1, "%)", 
            sep = "")
        lab.y <- paste("component n.", axes[2], " (", cp2, "%)", 
            sep = "")
        PC <- cpres %*% Vect
        Correl <- cor(cpres, PC)
        if (is.null(nomVar)) {
            nomVar <- paste("Variable ", 1:p, sep = "")
        }
        rownames(Correl) <- nomVar
        plot(Correl[, axes[1]], Correl[, axes[2]], xlab = lab.x, 
            ylab = lab.y, xlim = c(-1.5, 1.5), ylim = c(-1.5, 
                1.5), asp = 1, main = "Correlation Map")
        text(Correl[, axes[1]], Correl[, axes[2]], labels = rownames(Correl), 
            pos = 4)
        segments(0, 0, Correl[, axes[1]], Correl[, axes[2]], 
            col = 1:p)
        abline(h = 0, col = "blue")
        abline(v = 0, col = "blue")
        GenPC <- X %*% Vect
        PCinterval <- vector("list", d)
        PCinterval <- fintab(as.data.frame(GenPC))
#
regad<-ucad(p)
p2=nrow(regad)

for(ii in 1:p2)
{
colnames( PCinterval)[ii]<-paste(regad[ii,1])
}
#
        dev.new()
        n1 <- 1
        n2 <- nrow(xmin[[1]])
        m1 <- n2 + 1
        m2 <- 2 * n2
        lab.x <- paste(" Dim n.", axes[1], " (", cp1, "%)", sep = "")
        lab.y <- paste(" Dim n.", axes[2], " (", cp2, "%)", sep = "")
        if (is.null(nomInd)) {
            nomInd <- paste("Row ", 1:n2, sep = "")
        }
        if (is.null(xlim)) {
            xlim <- range(GenPC)
        }
        if (is.null(ylim)) {
            ylim <- range(GenPC)
        }
        if (is.null(legend)) {
            legend[1] <- max(xlim)
            legend[2] <- max(ylim) + 0.1
        }
        xlegend <- legend[1]
        ylegend <- legend[2] + 0.1
        plot(GenPC[, axes[1]], GenPC[, axes[2]], xlab = lab.x, 
            ylab = lab.y, asp = 1, main = "Factorial Map", xlim = xlim, 
            ylim = ylim)
        text(GenPC[1:n2, axes[1]], GenPC[1:n2, axes[2]], labels = nomInd, 
            pos = 4,col=1:n)
        abline(h = 0, col = "blue")
        abline(v = 0, col = "blue")
        rect(GenPC[n1:n2, axes[1]], GenPC[n1:n2, axes[2]], GenPC[m1:m2, 
            axes[1]], GenPC[m1:m2, axes[2]], border = 1:n)
        if (is.null(nametable)) {
            nametable <- paste("table n.")
        }
        #legend(xlegend, ylegend, nametable, col = c(1:d), 
         #   text.col = 1:d, lty = 1:d, pch = 1:d, merge = TRUE, 
          #  bg = "gray90")
        colnames(Correl) <- paste("Component", 1:p, sep = " ")
        colnames(PC) <- paste("Component", 1:p, sep = " ")
        write.table(cbind(nomVar, Correl), file = "CorrelInter.txt", 
            sep = " ", row.names = FALSE)
        bsun <- ".txt"
        write.table(PCinterval, paste(bsun, sep = "."), quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

#### graphique 3d des 3 premieres composantes

united=seq(1,2*p,2)
res<-matrix(0,p,2)
for(i in 1:p)
{
res[i,]<-c(united[i], united[i]+1)
}


zPC<-vector('list',p)
for(i in 1:p)
{
zPC[[i]]<-PCinterval[,res[i,]]
}


pc1<-zPC[[axes2[1]]]
pc2<-zPC[[axes2[2]]]
pc3<-zPC[[axes2[3]]]

dev.new()
 xmin1 <- min(pc1[, 1])
        xmin2 <- min(pc1[, 2])
        xmin <- min(xmin1, xmin2)
        xmax1 <- max(pc1[, 1])
        xmax2 <- max(pc1[, 2])
        xmax <- max(xmax1, xmax2)
        ymin1 <- min(pc2[, 1])
        ymin2 <- min(pc2[, 2])
        ymin <- min(ymin1, ymin2)
        ymax1 <- max(pc2[, 1])
        ymax2 <- max(pc2[, 2])
        ymax <- max(ymax1, ymax2)
        zmin1 <- min(pc3[, 1])
        zmin2 <- min(pc3[, 2])
        zmin <- min(zmin1, zmin2)
        zmax1 <- max(pc3[, 1])
        zmax2 <- max(pc3[, 2])
        zmax <- max(zmax1, zmax2)
        picto <- scatterplot3d(c(xmin, xmax), c(ymin, ymax), c(zmin, 
            zmax), xlab = paste('PC',axes2[1],sep='.'), ylab =paste('PC',axes2[2],sep='.'), 
            zlab = paste('PC',axes2[3],sep='.'), type = "n",main='Scatterplot3d')

cube <- rbind(c(1, 1, 1), c(2, 1, 1), c(2, 1, 2), c(1, 
            1, 2), c(1, 1, 1), c(1, 2, 1), c(1, 2, 2), c(2, 2, 
            2), c(2, 2, 1), c(1, 2, 1), c(1, 2, 2), c(1, 1, 2), 
            c(2, 1, 2), c(2, 2, 2), c(2, 2, 1), c(2, 1, 1))
for (i in 1:n) {
            vec.x <- pc1[i, cube[, 1]]
            vec.y <- pc2[i, cube[, 2]]
            vec.z <- pc3[i, cube[, 3]]
            picto$points3d(vec.x, vec.y, vec.z, type = "l", lty = 1, 
                col = i)
        }


ltext <- nomInd
            textPoints <- cbind((pc1[,1] + pc1[,2])/2, (pc2[,1] + pc2[,2])/2, (pc3[,1] + pc3[,2])/2)
            text(picto$xyz.convert(textPoints), labels = ltext,col=1:n)
       
        return(list(Correl = Correl, 
            PCinterval = PCinterval, Pval = pval2, Vect=Vect))
    }
}
