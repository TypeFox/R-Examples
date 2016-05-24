roc.plot.ade <-
function( pred, event, group=NULL, data=NULL, vnames=NULL, main='', xlab='1-Specificity', ylab='Sensitivity', digits=3, pdigs=4, lty=1, lwd=2, col=NULL, tcol=NULL,  bgcol=NULL, wall=0, test=FALSE, CC=TRUE, auc=TRUE, diag=TRUE, spec=FALSE){
if(any(par('mfg')!=c(1,1,1,1)) & any(par('mai') < c(1.02, 0.82, 0.82, 0.42))){
maidiff<-rep(0, 4)
norm<-c(1.02, 0.82, 0.82, 0.42)
maidiff[par('mai')<norm]<-  norm[par('mai')<norm] - par('mai')[par('mai')<norm]
par(mai=par('mai')+maidiff)
}
oldpar<-par(no.readonly =TRUE)
oldpar<-oldpar[-which(names(oldpar)%in%c('usr',  'plt',   'pin', 'fin', 'fig', 'mfg', 'mfcol', 'mfrow', 'omd', 'omi', 'oma'))]
on.exit(par(oldpar))
drawauc<-auc


#####################################
# without Data.Frame or with
if((is.list(pred) | is.matrix(pred) | is.numeric(pred)) & (is.numeric(event) | is.factor(event)) ){
data<-NULL
if(is.matrix(pred)) prednames<-colnames(pred)
if(is.list(pred))   prednames<-paste('Predictor_', 1:length(pred), sep='')
if(!is.matrix(pred) & is.numeric(pred))   prednames<-gsub('[(]{0}[A-Za-z0-9]*[$]', '', deparse(substitute(pred)))
if(is.null(prednames) & is.matrix(pred))  prednames<-paste('Predictor_', 1:dim(pred)[2], sep='')
if(is.list(pred)){
pred[[length(pred)+1]]<-event
if(!is.null(group)) pred[[length(pred)+1]]<-group
data<-as.data.frame(pred)
if(!is.null(group))  names(data)<-c(prednames, 'a.event', 'a.group')
if( is.null(group))  names(data)<-c(prednames, 'a.event')
}

if(is.matrix(pred) | is.numeric(pred)){
pred<-cbind(pred, event)
if(!is.null(group)) pred<-cbind(pred, group)
data<-as.data.frame(pred)
if(!is.null(group)) names(data)<-c(prednames, 'a.event', 'a.group')
if( is.null(group)) names(data)<-c(prednames, 'a.event')
data$a.group<-group
}


if(!is.null(group)) group<-'a.group'
pred<-prednames
event<-'a.event'
}
#####################################



a.N<-length(pred)
if (length(pred)>1) group=NULL
if (is.null(vnames)) vnames <- group
if (is.null(vnames)) vnames <- pred
if (length(pred)>1 & !(length(vnames)==length(pred)))  vnames <- pred
if(!is.null(data)){ if(!is.data.frame(data))  stop("(data) must be a data.frame!") }




# Remove Missings
if(CC){
for (j in 1:length(pred)) {
data<- subset(data, !is.na(eval(parse(text=paste(pred[j], sep='')))))
}
}
data<- subset(data, !is.na(eval(parse(text=paste(event, sep='')))))

if(!is.null(group)) data<- subset(data, !is.na(eval(parse(text=paste(group, sep='')))))
if(!is.null(group)){ g<-eval(parse(text=paste("data$",group, sep='')))
if(!is.factor(g))   g <- as.factor(g)
a.N<-nlevels(g)
}


# Test

#______________________________________________________________________________#
#------------------------------------------------------------------------------#
################################################################################
a.roc.test.ade<- function (y, x, L = NULL)
{
    trapezarea <- function(x, y) {
        if (x[1] > x[length(x)]) {
            x <- rev(x)
            y <- rev(y)
        }
        if (length(x) != length(y))
            stop("length x must equal length y")
        if (length(unique(x)) < 2)
            return(NA)
        ya <- approx(x, y, 0, ties = max, rule = 2)$y
        yb <- approx(x, y, 1, ties = max, rule = 2)$y
        x <- c(0, x, 1)
        y <- c(ya, y, yb)
        h <- diff(x)
        lx <- length(x)
        area <- 0.5 * sum(h * (y[-1] + y[-lx]))
        area
    }
    th <- NULL
    sens <- spec <- list(rep(NULL, length(x)))
    for (j in 1:length(x)) {
        DD <- table(-x[[j]], y)
        sens[[j]] <- c(0, cumsum(DD[, 2])/sum(DD[, 2]))
        spec[[j]] <- c(1, 1 - cumsum(DD[, 1])/sum(DD[, 1]))
        th[j] <- trapezarea(1 - spec[[j]], sens[[j]])
    }
    if (!is.null(names(x))) {
        names(sens) <- names(x)
        names(spec) <- names(x)
    }
    else {
        names(sens) <- paste("Test", LETTERS[1:length(x)])
        names(spec) <- paste("Test", LETTERS[1:length(x)])
    }
    if (!is.null(L)) {
        V10 <- matrix(NA, nrow = length(y[y == 1]), ncol = length(x))
        V01 <- matrix(NA, nrow = length(y[y == 0]), ncol = length(x))
        for (j in 1:length(x)) {
            x.s <- split(x[[j]], y)
            for (i in 1:length(x.s$"1")) V10[i, j] <- (length(x.s$"0"[x.s$"0" <
                x.s$"1"[i]]) + 0.5 * length(x.s$"0"[x.s$"0" ==
                x.s$"1"[i]]))/length(y[y == 0])
            for (i in 1:length(x.s$"0")) V01[i, j] <- (length(x.s$"1"[x.s$"0"[i] <
                x.s$"1"]) + 0.5 * length(x.s$"1"[x.s$"1" == x.s$"0"[i]]))/length(y[y ==
                1])
        }
        S10 <- (t(V10) %*% V10 - length(y[y == 1]) * th %*% t(th))/(length(y[y ==
            1]) - 1)
        S01 <- (t(V01) %*% V01 - length(y[y == 0]) * th %*% t(th))/(length(y[y ==
            0]) - 1)
        S <- S10/length(y[y == 1]) + S01/length(y[y == 0])
        contr <- L %*% th
        se <- sqrt((L %*% S %*% t(L)))
        test <- t(th) %*% t(L) %*% solve(L %*% (1/length(y[y ==
            1]) * S10 + 1/length(y[y == 0]) * S01) %*% t(L),
            t(t(th) %*% t(L)))
        p.value <- pchisq(test, df = qr(L %*% S %*% t(L))$rank,
            lower.tail = FALSE)
    }
    else {
        S <- NULL
        p.value <- NULL
        contr <- NULL
    }
    names(th) <- names(x)
    out <- list(th, p.value)
    return(out)

}
#______________________________________________________________________________#
#------------------------------------------------------------------------------#
################################################################################



if(test){
Blist<-list(NULL)
for(i in 1:length(pred)){
Blist[[i]]<-eval(parse(text=paste("data$",pred[i], sep='')))
}

testres <- a.roc.test.ade(eval(parse(text=paste("data$",event, sep=''))), x = Blist, L = t(contr.sum(length(Blist), contrasts = TRUE)))
pvaltest<-format_p.ade(testres[[2]], pdigs)
}


# Colors
if(is.null(tcol)  & wall==0)   tcol<-1
if(is.null(tcol)  & wall!=0)   tcol<-rgb(0.1,0.1,0.25)
if(is.null(bgcol) & wall==0)   bgcol<-1
if(is.null(bgcol) & wall!=0)   bgcol<-'#DBE0E8'
if(is.null(col) & a.N==1) col <- tcol
if(is.null(col) & a.N>1)  col <- a.getcol.ade(a.N)

################################################################################
################################################################################
################################################################################


#______________________________________________________________________________#
#------------------------------------------------------------------------------#
################################################################################
a.roc.calcs<-function(logistic.model)
{

ppfit<-logistic.model$fitted.values
if(length(logistic.model$coefficients)==2){
if(logistic.model$coefficients[2]<0)  ppfit<- (-1*ppfit)
}
    if (length(grep("cbind", names(model.frame(logistic.model)))) >  0) {
        firsttable1 <- cbind(ppfit, model.frame(logistic.model)[,1][, 2:1])
        firsttable1 <- firsttable1[order(firsttable1[, 1]), ]
    }
    else {
        if (length(grep("(weights)", names(model.frame(logistic.model)))) >
            0) {
            firsttable <- xtabs(as.vector(model.frame(logistic.model)[,
                ncol(model.frame(logistic.model))]) ~ ppfit +
                logistic.model$y)
        }
        else {
            firsttable <- table(ppfit,
                logistic.model$y)
        }
        colnames(firsttable) <- c("Non-diseased", "Diseased")
        firsttable1 <- cbind(as.numeric(rownames(firsttable)),
            firsttable)
    }

    colnames(firsttable1)[1] <- "predicted.prob"
    firsttable <- firsttable1[, 2:3]
    secondtable <- firsttable

sum1.ade<-sum(firsttable[, 1])
sum2.ade<-sum(firsttable[, 2])
run.sum1.ade<-0
run.sum2.ade<-0
    for (i in 1:length(secondtable[, 1])) {
run.sum1.ade<-run.sum1.ade + firsttable[i, 1]
run.sum2.ade<-run.sum2.ade + firsttable[i, 2]
        secondtable[i, 1] <- (sum1.ade - run.sum1.ade)/sum1.ade
        secondtable[i, 2] <- (sum2.ade - run.sum2.ade)/sum2.ade
    }

    secondtable <- rbind((c(1, 1)), secondtable)
    colnames(secondtable) <- c("1-Specificity", "Sensitivity")
    model.des <- deparse(logistic.model$formula)
    auc <- 0
    for (i in 1:(nrow(secondtable) - 1)) {
        auc <- auc + (secondtable[i, 1] - secondtable[(i + 1),
            1]) * 0.5 * (secondtable[i, 2] + secondtable[(i +
            1), 2])
    }

    list(model.description = model.des, auc = auc, predicted.table = firsttable1, diagnostic.table = secondtable)
}
#______________________________________________________________________________#
#------------------------------------------------------------------------------#
################################################################################





#####  Style 0 #################################################################
if(wall==0){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0) )
a1<-axis(1, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), col.ticks=bgcol, lwd.ticks=1)
if(diag) abline(0, 1, lty='dashed', col=bgcol)
box(col=bgcol)
}


#  Legend  #
legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc)  letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
a.adj <- c(0, 0.5)
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
legend(x=1, y=0, legend=letext , col=c(col, 0), adj=a.adj, box.lwd=1, lty = lty, lwd=3, merge =FALSE, yjust=0, xjust=1, box.col=bgcol, text.col=tcol, text.width=max(strwidth(letext,font = 2)))
}
}
################################################################################



#####  Style 1 #################################################################
if(wall==1){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0) )
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=tcol, labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
if(diag) abline(0, 1, lty='dashed', col=rgb(1,1,1))
box(col=rgb(1,1,1))
}


#  Legend  #
legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc) letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
a.adj <- c(0, 0.5)
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
legend(x=1, y=0, legend=letext , col=c(col, 0), adj=a.adj, lty = lty, lwd=3, merge =FALSE, yjust=0, xjust=1, border=rgb(1,1,1), box.lwd=2,  box.col=rgb(1,1,1), bg=bgcol, text.col=tcol, text.width=max(strwidth(letext,font = 2)))
}

}
################################################################################




#####  Style 2 #################################################################
if(wall==2){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0) )
a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -75), lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=a.coladd.ade(bgcol, -75), labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=1)
abline(v=a1, h=a2, lty=1,   col=bgcol, lwd=1)
if(diag) abline(0, 1, lty='dashed', col=bgcol)
box(col=a.coladd.ade(bgcol, -75))
}


#  Legend  #
legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc) letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
a.adj <- c(0, 0.5)
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
legend(x=1, y=0, legend=letext , col=c(col, 0), adj=a.adj, lty = lty, bg=rgb(1,1,1), lwd=3, merge =FALSE, yjust=0, xjust=1, border=rgb(1,1,1), box.lwd=1,  box.col=a.coladd.ade(bgcol, -75), text.col=tcol, text.width=max(strwidth(letext,font = 2)))
}

}
################################################################################




#####  Style 3 #################################################################
if(wall==3){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0) )
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -50), lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=a.coladd.ade(bgcol, -50), labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=1)
abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -50), lwd=1)
if(diag) abline(0, 1, lty='dashed', col=a.coladd.ade(bgcol, -50))
box(col=a.coladd.ade(bgcol, -50))
}


#  Legend  #
legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc) letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
a.adj <- c(0, 0.5)
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
legend(x=1, y=0, legend=letext , col=c(col, 0), adj=a.adj, lty = lty, bg=rgb(1,1,1), lwd=3, merge =FALSE, yjust=0, xjust=1, border=a.coladd.ade(bgcol, -50), box.lwd=1,  box.col=a.coladd.ade(bgcol, -50), text.col=tcol, text.width=max(strwidth(letext,font = 2)))
}

}
################################################################################




#####  Style 4 #################################################################
if(wall==4){
par(col.axis=tcol)
par(col.lab=rgb(1,1,1))
par(col.main=rgb(1,1,1))
par(font=2)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim , col=rgb(1,1,1,0))
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=tcol, labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=1)

abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
if(diag) abline(0, 1, lty='dashed', col=rgb(1,1,1))


# Outer
par(xpd=TRUE)
if(!spec) polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0, 2.75,  2.75, 0)), col=tcol, border=rgb(1,1,1))
if(spec)  polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(1.3, 3,  3, 1.3)), col=tcol, border=rgb(1,1,1))
if(!spec) text(a.glc(side=0), a.glc(side=3, line=1),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))
if(spec)  text(a.glc(side=0), a.glc(side=3, line=1.8),    labels=main, cex = 1.25, font=2, col=rgb(1,1,1), adj=c(0.5,0))


if(ylab!='' & ylab!=' ') polygon( a.glc(side=2, line=c(3.5, 3.5, 2, 2)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=rgb(1,1,1))
if(xlab!='' & xlab!=' ') polygon( a.glc(side=c(2, 2, 4, 4), line=0),     a.glc(side=1, line=c(4, 2.5, 2.5, 4)), col=bgcol, border=rgb(1,1,1))
text(a.glc(side=0), a.glc(side=1, line=3.5), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,  font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(lwd=1, col=rgb(1,1,1))
}


#  Legend  #
legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc) letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
if(test)  blcol<-c(rep(rgb(1,1,1), length(col)-1), rgb(1,1,1,0))
if(!test) blcol<-c(rep(rgb(1,1,1), length(col)))
legend('bottomright', legend=letext, lwd=5  ,  col = blcol,  box.lwd=1, lty = 1,    merge =FALSE, yjust=0, xjust=1, box.col=rgb(1,1,1), text.col=rgb(1,1,1), bg=tcol, text.width=max(strwidth(letext,font = 2)))
legend('bottomright', legend=letext, lwd=3  ,  col = col,  box.lwd=1, lty = lty,  merge =FALSE, yjust=0, xjust=1, box.col=rgb(1,1,1), text.col=rgb(1,1,1, 0), bg=a.alpha.ade(1,0), text.width=max(strwidth(letext,font = 2)))
}

}
################################################################################



#####  Style 5 #################################################################
if(wall==5){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)
par(font=2)
newmai<-rep(0, 4)
oldmai<-par('mai')
if(oldmai[2]<1) newmai[2]<- 1 - oldmai[2]
if(oldmai[3]>0.75 & oldmai[3]<=0.82) newmai[3]<- 0.75-oldmai[3]
if(oldmai[4]>0.25 & oldmai[4]<=0.42) newmai[4]<- 0.25-oldmai[4]
par(mai=(oldmai+newmai))


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', xlab='', ylab='', xlim=xlim, ylim=ylim , col=rgb(1,1,1,0))
a1<-axis(1, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=tcol, lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=tcol, labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=1)

if(diag) abline(0, 1, lty='dashed', col=bgcol)




# Outer
par(xpd=TRUE)
if(!spec) polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
if(!spec) polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=rgb(1,1,1,0), border=tcol)
if(!spec) polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(0.6, 3, 3, 0.6)), col=bgcol,        border=tcol)
if(!spec) polygon(a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 0.6, 0.6, 2.6)), col=bgcol,  border=tcol)

if(spec)  polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=3, line=c(1.6, 3.25, 3.25, 1.6)), col=bgcol,        border=tcol)
if(spec)  polygon(a.glc(side=c(2,2,4,4), line=c(0,0,0,0)), a.glc(side=3, line=c(1.6, 3.25, 3.25, 1.6)), col=rgb(1,1,1,0), border=tcol)
if(spec)  polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)),   a.glc(side=3, line=c(1.6, 3.25, 3.25, 1.6)), col=bgcol,        border=tcol)
if(spec)  polygon(a.glc(side=2, line=c(4.25, 4.25 ,3.65, 3.65)),  a.glc(side=c(1,3,3,1), line=c(2.6, 1.6, 1.6, 2.6)), col=bgcol,  border=tcol)


polygon(a.glc(side=4, line=c(0, 0 ,0.6, 0.6)), a.glc(side=c(1, 3, 3, 1), line=0), col=bgcol, border=tcol)
polygon(a.glc(side=2, line=c(4.25, 4.25, 0, 0)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)
polygon(a.glc(side=c(2, 2, 4, 4), line=0), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=rgb(1,1,1,0), border=tcol)
polygon(a.glc(side=4, line=c(0, 0, 0.6, 0.6)), a.glc(side=1, line=c(2.6, 4.5, 4.5, 2.6)), col=bgcol, border=tcol)

if(!spec) text(a.glc(side=0), a.glc(side=3, line=1.5),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))
if(spec)  text(a.glc(side=0), a.glc(side=3, line=2),  labels=main, cex = 1.25, font=2, col=tcol, adj=c(0.5,0))

text(a.glc(side=0), a.glc(side=1, line=3.75), labels=xlab, cex = 1.1,  font=2, col=tcol, adj=c(0.5,0))
text(a.glc(side=2, line=2.5), a.glc(side=5),  labels=ylab, cex = 1.1,   font=2,  col=tcol, adj=c(0.5,0), srt=90)
par(xpd=FALSE)
box(lwd=1, col=tcol)
}

legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc) letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
legend('bottomright', legend=letext, lwd=3  ,  col = col, box.lwd=1, lty = lty,  merge =FALSE, yjust=0, xjust=1, border=tcol, box.col=tcol, text.col=tcol, bg=rgb(1,1,1, 0), text.width=max(strwidth(letext,font = 2)))
}

}
################################################################################


#####  Style 6 #################################################################
if(wall==6){
par(col.axis=tcol)
par(col.lab=tcol)
par(col.main=tcol)


#  Plot  #
plot.box.ade<-function(xlab, ylab, main, xlim, ylim, lwd){
plot(0, 0, type='s', main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=rgb(1,1,1,0) )
polygon( c(par('usr')[c(1,1,2,2)]), par('usr')[c(3,4,4,3)], col=bgcol, border=FALSE)
a1<-axis(1, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a1<-axis(1, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1)
a2<-axis(2, col=rgb(1,1,1), col.ticks=a.coladd.ade(bgcol, -35), lwd.ticks=3)
a2<-axis(2, col=rgb(1,1,1), col.ticks=rgb(1,1,1), lwd.ticks=1)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=a.coladd.ade(bgcol, -35), labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=3)
if(spec) a3<-axis(3, padj =0.9, tck=-0.01, at=pretty(c(1,0)), col.ticks=rgb(1,1,1)              , labels=pretty(c(1,0))[length(pretty(c(1,0))):1], col=rgb(1,1,1), lwd.ticks=1)


abline(v=a1, h=a2, lty=1, col=a.coladd.ade(bgcol, -35), lwd=3)
abline(v=a1, h=a2, lty=1, col=rgb(1,1,1), lwd=1)
if(diag) abline(0, 1, lty='dashed', col=rgb(1,1,1))
box(lwd=3, col=rgb(1,1,1))
box(lwd=1, col=a.coladd.ade(bgcol, -35))
}


#  Legend  #
legens.ade<-function(vnames, groups, auc, p, col){
letext<-NULL
if(drawauc) letext=paste(vnames,' ',groups," (AUC: ", auc, ")", sep='')
if(!drawauc) letext=paste(vnames,' ',groups, sep='')
if(test) letext<-c(letext, paste('p-Value:', p))
a.adj <- c(0, 0.5)
if(length(groups)==1 & length(auc)==1)  a.adj <- c(0.1, 0.5)
legend(x=1, y=0, legend=letext , col=c(col, 0), adj=a.adj, lty = lty, lwd=3, merge =FALSE, yjust=0, xjust=1, border=rgb(1,1,1), box.lwd=3,  box.col=rgb(1,1,1),  bg=bgcol, text.col=tcol, text.width=max(strwidth(letext,font = 2)))
legend(x=1, y=0, legend=letext , col=c(col, 0), adj=a.adj, lty = lty, lwd=3, merge =FALSE, yjust=0, xjust=1, border=rgb(1,1,1), box.lwd=1,  box.col=a.coladd.ade(bgcol, -35), bg=rgb(0,0,0,0), text.col=tcol, text.width=max(strwidth(letext,font = 2)))

}

}
################################################################################



################################################################################
################################################################################
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

################################################################################
################################################################################
################################################################################

################################################################################
# 1 Gruppe
if(is.null(group)){


plot.box.ade( xlab, ylab, main, xlim=c(0 ,1), ylim=c(0, 1) , lwd=lwd)

if(length(pred)==1){  
if(any(is.na(col))) col <- tcol

roc.model<-a.roc.calcs(glm(eval(parse(text=paste("data$",event)))~eval(parse(text=paste("data$",pred))), family='binomial'))
x<- roc.model$diagnostic.table[ , 1]
y<- roc.model$diagnostic.table[ , 2]
auc <- round(roc.model$auc, digits=digits)


legens.ade(vnames, '', auc, pvaltest, rgb(0,0,0,0))
points(x, y, type='l', xlim=c(0 ,1), ylim=c(0, 1) , lty=lty, lwd=lwd, col=col)
}


if(length(pred)>1){
g<-NULL  
ltexts<-NULL
if(length(lty)<length(pred))  lty<-rep(lty, length(pred))
for(i in 1:length(pred)){
g[i]<-vnames[i]

# Calcs
roc.model<-a.roc.calcs(glm(eval(parse(text=paste("data$",event)))~eval(parse(text=paste("data$",pred[i]))), family='binomial'))
x<- roc.model$diagnostic.table[ , 1]
y<- roc.model$diagnostic.table[ , 2]
auc <- round_n.ade(roc.model$auc, digits=digits)

ltexts<-c(ltexts, auc)
points(x, y, type='l', xlim=c(0 ,1), ylim=c(0, 1) , lty=lty[i], lwd=lwd, col=col[i])
}
if(test) col<-c(col, rgb(1,1,1,0))
legens.ade(vnames, '', ltexts, pvaltest, col)
}

}
################################################################################         


################################################################################
# Mehrere Gruppen
if(!is.null(group)){
if(nlevels(g)>1){

plot.box.ade( xlab, ylab, main, xlim=c(0 ,1), ylim=c(0, 1), lwd=lwd)
ltexts<-NULL
if(length(lty)<nlevels(g))  lty<-rep(lty, nlevels(g))
for(i in 1:nlevels(g)){
subdata <-  subset(data, eval(parse(text=group))==levels(g)[i])
roc.model<-a.roc.calcs(glm(eval(parse(text=paste("subdata$",event)))~eval(parse(text=paste("subdata$",pred))), family='binomial'))
x<- roc.model$diagnostic.table[ , 1]
y<- roc.model$diagnostic.table[ , 2]
auc <- round_n.ade(roc.model$auc, digits=digits)
ltexts<-c(ltexts, auc)
points(x, y, type='l', xlim=c(0 ,1), ylim=c(0, 1) , lty=lty[i], lwd=lwd, col=col[i])
}
legens.ade(vnames, levels(g), ltexts, '', col)

}


################################################################################ 
################################################################################ 

}

}
