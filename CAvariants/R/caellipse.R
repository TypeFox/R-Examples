caellipse <-
function (N, a1=1,a2=2,alpha = 0.05, cols = c(2, 4), M = 2, cex = .8, cex.lab = 0.5, mar = c(5, 4, 4, 2) + 0.1, 
    prop=.8,Imass,Jmass,a,b,g,fr,dmu,inertiapc,plottype="biplot",biptype="row",pos=2,arrow=T,length=0,graphy=T) 
{
    I <- nrow(Imass)
    J <- nrow(Jmass)
n<-sum(N)
rowgroup <- list(1:I,rep(1,I))
rowgrlab <- list(1,"","*","blue","T")
colgroup <- list(1:J,rep(1,J))
colgrlab <- list(1,"","+","red","T")
    Inames <- dimnames(fr)[1]
    Jnames <- dimnames(g)[1]
    dIh <- solve(Imass^0.5)
    dJh <- solve(Jmass^0.5)
    t.inertia <- sum(dmu)
dmu<-sqrt(dmu)
#browser()
#if ((catype=="SOCA")|(catype=="DOCA")){
 #a <- dIh %*% a
 #b <- dJh %*% b}
chisq.val <- qchisq(1 - alpha, df = (I - 1) * (J - 1))
    hlax1.row <- vector(mode = "numeric", length = I)
    hlax2.row <- vector(mode = "numeric", length = I)
    hlax1.col <- vector(mode = "numeric", length = J)
    hlax2.col <- vector(mode = "numeric", length = J)
#browser()
    if (M > 2) {
        for (i in 1:I) {
            hlax1.row[i] <- dmu[1,1] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Imass[i,i] - sum(a[i, 3:M]^2))))
            hlax2.row[i] <- dmu[2,2] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Imass[i,i] - sum(a[i, 3:M]^2))))
        }
        for (j in 1:J) {
            hlax1.col[j] <- dmu[1,1] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
            hlax2.col[j] <- dmu[2,2] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
        }
    }
    else {
        for (i in 1:I) {
            hlax1.row[i] <- dmu[1,1] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/(Imass)[i,i])))
            hlax2.row[i] <- dmu[2,2] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/(Imass)[i,i])))
        }
        for (j in 1:J) {
            hlax1.col[j] <- dmu[1,1] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/(Jmass)[j,j])))
            hlax2.col[j] <- dmu[2,2] * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/(Jmass)[j,j])))
        }
    }
if (graphy==TRUE){
picsize<-c(range(fr[,c(a1,a2)], g[,c(a1,a2)])/prop)
    plot.new()
    par(pty = "s", mar = mar)
#browser()
#plot(fr[, a1], fr[, a2], xlim=picsize, ylim=picsize, 
plot(0,0,pch=" ", xlim=picsize, ylim=picsize, 
 xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
 ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  asp=1,col=rowgrlab[[4]][as.integer(rowgroup[[2]])],cex=cex)
abline(h = 0, v = 0) 
text(fr[, a1], fr[, a2], labels = Inames[[1]], pos= pos, col = cols[1], cex = cex)
points(fr[, a1], fr[, a2], pch="*",col=cols[1])
    text(g[, a1], g[, a2], labels = Jnames[[1]], pos= pos, col = cols[2], 
        cex = cex)
points(g[, a1], g[, a2], pch="+",col=cols[2])
       title(main = paste(100 * (1 - alpha), "% Confidence Ellipses"))
if (plottype=="biplot") {
if ((biptype=="column")|(biptype=="col")|(biptype=="c")) {
#----------------------------------------------arrow on column principal coords
nv <- rep(0, nrow(g))
#vec <- g[, c(a1,a2)]
if(arrow) {
arrows(nv, nv, g[,a1], g[, a2], length = length)
}
 for (j in 1:J) {
        ellipse(hlax1.col[j], hlax2.col[j], xc = g[j, a1], yc = g[j, 
            a2], col = cols[2])
    }
} #end col
if ((biptype=="row")|(biptype=="r")|(biptype=="rows")){
#----------------------------------------------arrow on row principal coords
nv <- rep(0, nrow(fr))
#vec <- fr[, c(a1, a2)]
if(arrow) {
arrows(nv, nv, fr[, a1], fr[,a2], length = length)
}

 for (i in 1:I) {
        ellipse(hlax1.row[i], hlax2.row[i], xc = fr[i, a1], yc = fr[i, 
            a2], col = cols[1])
    }
} #end row
}#end biplot

else{ # when classical 

 for (j in 1:J) {
        ellipse(hlax1.col[j], hlax2.col[j], xc = g[j, a1], yc = g[j, 
            a2], col = cols[2])
    }
 
for (i in 1:I) {
        ellipse(hlax1.row[i], hlax2.row[i], xc = fr[i, a1], yc = fr[i, 
           a2], col = cols[1])
    }
}#end else
}#end graphy
    eccentricity <- sqrt(1 - (dmu[2, 2]/dmu[1, 1])^2)
    area.row <- vector(mode = "numeric", length = I)
    area.col <- vector(mode = "numeric", length = J)
    for (i in 1:I) {
        area.row[i] <- 3.14159 * hlax1.row[i] * hlax2.row[i]
    }
    for (j in 1:J) {
        area.col[j] <- 3.14159 * hlax1.col[j] * hlax2.col[j]
    }
    pvalrow <- vector(mode = "numeric", length = I)
    pvalcol <- vector(mode = "numeric", length = J)
    for (i in 1:I) {
        if (M > 2) {
            pvalrow[i] <- 1- pchisq(t.inertia *n* ((1/Imass[i,i] - sum(a[i, 
                3:M]^2))^(-1)) * (fr[i, 1]^2/dmu[1, 1]^2 + fr[i, 
                2]^2/dmu[2, 2]^2), df = (I - 1) * (J - 1))
        }
        else {
            pvalrow[i] <- 1-pchisq(t.inertia*n * (Imass[i,i]) * 
                (fr[i, 1]^2/dmu[1, 1]^2 + fr[i, 2]^2/dmu[2, 2]^2), 
                df = (I - 1) * (J - 1))
        }
    }
    for (j in 1:J) {
        if (M > 2) {
            pvalcol[j] <-  1-pchisq(t.inertia *n* ((1/Jmass[j,j] - 
                sum(b[j, 3:M]^2))^(-1)) * (g[j, 1]^2/dmu[1, 1]^2 + 
                g[j, 2]^2/dmu[2, 2]^2), df = (I - 1) * (J - 1))
        }
        else {
            pvalcol[j] <- 1- pchisq(t.inertia *n* (Jmass[j,j]) * 
                (g[j, 1]^2/dmu[1, 1]^2 + g[j, 2]^2/dmu[2, 2]^2), 
                df = (I - 1) * (J - 1))
        }
    }
    summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
    row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
    dimnames(row.summ) <- list(paste(Inames[[1]]), paste(summ.name))
    col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
    dimnames(col.summ) <- list(paste(Jnames[[1]]), paste(summ.name))
#browser() 
 invisible(list(eccentricity = eccentricity, row.summ = row.summ, 
        col.summ = col.summ))
}
