mypostscript=function (file="temp", mfrow=c(1,1), mfcol=NULL, width=NULL, height=NULL, ext="eps", oma=NULL, mar=NULL,main.outer=FALSE, save2file=TRUE, ...) {    
    
    print(paste(getwd(),"/",file,sep=""))
    
    if (!is.null(mfcol)) {
        nrow=mfcol[1]; ncol=mfcol[2]        
    } else {
        nrow=mfrow[1]; ncol=mfrow[2]
    }
    
    #if (nrow>4) warning ("nrow > 4 will not fit a page without making the figures hard to see")
    
    # sca controls how much to scale down for use in a paper
    if(is.null(width) | is.null(height))  {
        if (nrow==1 & ncol==1) {width=6.7; height=6.7
        } else if (nrow==1 & ncol==2) {width=9.7; height=5.2
        } else if (nrow==1 & ncol==3) {width=9.7; height=3.4
        } else if (nrow==1 & ncol==4) {width=14; height=3.4
        } else if (nrow==2 & ncol==3) {width=9.7; height=6.7
        } else if (nrow==4 & ncol==6) {width=15; height=10
        } else if (nrow==2 & ncol==4) {width=13; height=6.7
        } else if (nrow==3 & ncol==6) {width=17.5; height=9
        } else if (nrow==4 & ncol==8) {width=17.5; height=9
        } else if (nrow==4 & ncol==9) {width=20; height=9
        } else if (nrow==3 & ncol==5) {width=15; height=9.6
        } else if (nrow==3 & ncol==4) {width=12; height=9.6
        } else if (nrow==4 & ncol==5) {width=15; height=12.5
        } else if (nrow==5 & ncol==6) {width=9; height=8.3
        } else if (nrow==2 & ncol==2) {width=8; height=8.5
        } else if (nrow==3 & ncol==3) {width=9.7; height=10.3
        } else if (nrow==4 & ncol==4) {width=9.7; height=10.3
        } else if (nrow==6 & ncol==5) {width=18; height=17
        } else if (nrow==5 & ncol==5) {width=15; height=15
        } else if (nrow==5 & ncol==3) {width=9; height=15
        } else if (nrow==4 & ncol==2) {width=6; height=13
        } else if (nrow==6 & ncol==3) {width=9; height=19
        } else if (nrow==7 & ncol==3) {width=9; height=22
        } else if (nrow==8 & ncol==5) {width=10; height=16
        } else if (nrow==6 & ncol==4) {width=12; height=19
        } else if (nrow==7 & ncol==5) {width=18; height=19
        } else if (nrow==5 & ncol==4) {width=12; height=15
        } else if (nrow==2 & ncol==1) {width=6.7; height=9.7
        } else if (nrow==3 & ncol==1) {width=10; height=9.7
        } else if (nrow==5 & ncol==1) {width=5; height=13
        } else if (nrow==3 & ncol==2) {width=6.7; height=10.3
        } else if (nrow==4 & ncol==3) {width=9; height=12
        } else stop ("nrow x ncol not supported: "%+%nrow%+%" x "%+%ncol)
    }
    
    if(save2file){      
        if (ext=="pdf") pdf (paper="special", file=file%+%"."%+%ext, width=width, height=height, ...)
        else postscript (paper="special", horizontal=FALSE, file=file%+%"."%+%ext, width=width, height=height, ...)
    } else {
        print("not saving to file")
    }
    
    if (!is.null(mfcol)) par(mfcol=mfcol)
    else par(mfrow=mfrow)    
    
    if (!is.null(oma)) par(oma=oma)
    if (!is.null(mar)) {
        par(mar=mar)
    }
    
    if (main.outer) {
        tmp=par()$oma
        tmp[3]=tmp[3]+1
        par(oma=tmp)
    }
    
}
mypdf=function (...) {mypostscript(ext="pdf",...)} # cannot print both to screen and pdf
##test
#mypdf(mfrow=c(1,1),file="test1x1");plot(1:10,main="LUMX",xlab="t",ylab="y");dev.off()
#mypdf(mfrow=c(1,2),file="test1x2");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);dev.off()
#mypdf(mfrow=c(2,2),file="test2x2");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(1,3),file="test1x3");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(2,3),file="test2x3");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(4,4),file="test4x4");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10,main="Luminex");dev.off()



# if lty is specified, a line will be drawn
mylegend=function(legend, x, lty=NULL,bty="n", ...) {
    x=switch(x, "topleft", "top", "topright", "left", "center" , "right", "bottomleft", "bottom", "bottomright")
    legend(bty=bty,x=x, legend=legend, lty=lty, ...)
}


# copied from pairs help page
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor,  ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method="pearson", use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.nothing=function(x, ...) {}
mypairs=function(dat, ...){
    pairs(dat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist, ...)
}


getMfrow=function (len) {
    ret=NULL
    if (len==1) { 
        ret=c(1,1)
    } else if (len==2) { 
        ret=c(1,2)
    } else if (len==3) { 
        ret=c(1,3)
    } else if (len<=4) { 
        ret=c(2,2)
    } else if (len<=6) { 
        ret=c(2,3)
    } else if (len<=9) { 
        ret=c(3,3)
    } else if (len<=12) { 
        ret=c(3,4)
    } else if (len<=16) { 
        ret=c(4,4)
    } else if (len<=20) { 
        ret=c(4,5)
    } else if (len<=25) { 
        ret=c(5,5)
    }
    ret
}


empty.plot=function () {
    plot(1,1,type="n",xlab="",ylab="",xaxt="n", yaxt="n", bty="n")
}

myforestplot=function(dat, xlim=NULL, xlab="", main="", col.1="red", col.2="blue",plot.labels=TRUE,order=FALSE,decreasing=FALSE) {
    if (order) dat=dat[order(dat[,1],decreasing=decreasing),] 
    p=nrow(dat)    
    # makes no plot, but helps set the x axis later
    plot(c(dat[,2], dat[,3]),rep(1:p,2), xlim=xlim, yaxt="n", xaxt="s", xlab=xlab, ylab="", main="", type="n", cex.main=1.4, axes=F)
    mtext(side=3, line=3, adj=0, text=main, cex=1.4, font=2, xpd=NA)
    if (range(dat[,2:3])[1]>0) null.val=1 else null.val=0 # if all values are greater than 0, 1 is probably the null, otherwise, 0 is probably the null
    abline(v=null.val, col="gray")
    cols=ifelse(dat[,4]<0.05, col.1, col.2)
    points(dat[,1], nrow(dat):1, pch=19, col=cols)
    segments(dat[,2], nrow(dat):1, dat[,3], nrow(dat):1, lwd=2, col=cols)
    axis(1, cex.axis=1.4)
    # add labels
    if (plot.labels) axis(4, at=p:1, rownames(dat), tick=F, las=2, col=1, cex.axis=1, xpd=NA, line=-.5)
}


# can use after myboxplot
# both dat must have two columns, each row is dat from one subject
# x.ori=0; xaxislabels=rep("",2); cex.axis=1; add=FALSE; xlab=""; ylab=""; pcol=NULL; lcol=NULL
my.interaction.plot=function(dat, x.ori=0, xaxislabels=rep("",2), cex.axis=1, add=FALSE, xlab="", ylab="", pcol=NULL, lcol=NULL, ...){
    if (!add) plot(0,0,type="n",xlim=c(1,2),ylim=range(dat), ylab=ylab, xlab=xlab, xaxt="n", ...)
    cex=.25; pch=19
    if (is.null(lcol)) lcol=ifelse(dat[,1]>dat[,2],"red","black")
    for (i in 1:nrow(dat)) {
        points (1+x.ori, dat[i,1], cex=cex, pch=pch, col=ifelse(is.null(pcol), 1, pcol[i,1]))
        points (2+x.ori, dat[i,2], cex=cex, pch=pch, col=ifelse(is.null(pcol), 1, pcol[i,2]))
        lines (1:2+x.ori, dat[i,], lwd=.25, col=lcol[i])
    }
    axis(side=1, at=1:2+x.ori, labels=xaxislabels, cex.axis=cex.axis)
}

myboxplot <- function(object, ...) UseMethod("myboxplot") 

# myboxplot.formula and myboxplot.list make a boxplot with data points and do inferences for two group comparions. 
# cex=.5; ylab=""; xlab=""; main=""; box=FALSE; highlight.list=NULL; at=NULL;pch=1;col=1;
myboxplot.formula=function(formula, data, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL,
    pch=1, col=1, test="", reshape.formula=NULL, jitter=TRUE, add.interaction=FALSE, ...){
    
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {        
        set.seed(1)
        save.seed <- get(".Random.seed", .GlobalEnv)
    }                        
    set.seed(1)
    
    
    if (box) {
        res=boxplot(formula, data, range=0, ylab=ylab, xlab=xlab, at=at, main=main, col=NULL, cex=cex, ...)
    } else {
        res=boxplot(formula, data, range=0, ylab=ylab, xlab=xlab, at=at, main=main, col=NULL, cex=cex, border="white", ...)
    }
        
    dat.tmp=model.frame(formula, data)
    xx=interaction(dat.tmp[,-1])
    if(is.null(at)){        
        xx=as.numeric(xx)
    } else{
        xx=at[xx]
    }      
    if (add.interaction) jitter=FALSE
    if (jitter) xx=jitter(xx)  
    points(xx, dat.tmp[[1]], cex=cex,pch=pch,col=col)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)     
    
    # inference
    x.unique=unique(dat.tmp[[2]])
    if (length(test)>0) {
        sub=""
        if ("t" %in% test) sub=sub%+%" Student's t "%+%signif(t.test(formula, data)$p.value,2) 
        if ("w" %in% test) sub=sub%+%" Wilcoxon "%+%signif(wilcox.test(formula, data)$p.value,2)
        if ("k" %in% test) sub=sub%+%" Kruskal "%+%signif(kruskal.test(formula, data)$p.value,2)
        if ("f" %in% test) {
            dat.wide=myreshapewide (reshape.formula, data, idvar = NULL)
            #str(dat.wide)# show this so that we know we are using the right data to do the test
            ftest = friedman.test (as.matrix(dat.wide[,-(1)]))
            if (add.interaction) my.interaction.plot(as.matrix(dat.wide[,-1]), add=T)
            sub=sub%+%" Friedman "%+%signif(ftest$p.value,2)
        }
        title(sub=sub)
    }
    
    res
    
}

myboxplot.data.frame=function(object, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL, pch=1, col=1, test="", ...){
    myboxplot.list(as.list(object), cex=cex, ylab=ylab, xlab=xlab, main=main, box=box, at=at, pch=pch, col=col, test=test, ...)
}
myboxplot.matrix=function(object, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL, pch=1, col=1, test="", ...){
    myboxplot.list(as.list(as.data.frame(object)), cex=cex, ylab=ylab, xlab=xlab, main=main, box=box, at=at, pch=pch, col=col, test=test, ...)
}

myboxplot.list=function(object, ...){

    # make a dataframe out of list object
    dat=NULL
    if(is.null(names(object))) names(object)=1:length(object)
    for (i in 1:length(object)) {
        dat=rbind(dat,data.frame(y=object[[i]], x=names(object)[i]))
    }
    myboxplot(y~x, dat, ...)
    
}


# called butterfly.plot, because it is meant to plot two treatment arms at two time points, the two arms are plotted in a mirror fashion, see "by analyte.pdf" for an example
# if dat2 is null: dat is matrix with four columns. each row is one subject, the columns will be plotted side by side, with lines connecting values from one ptid
# if dat2 is not null, dat has two columns, which are plotted side by side with lines connecting them, same for dat2
# if add is true, no plot function will be called
butterfly.plot=function (dat, dat2=NULL, add=FALSE, xaxislabels=rep("",4), x.ori=0, xlab="", ylab="", cex.axis=1, ...){
    if (!add) plot(0,0,type="n",xlim=c(1,4),ylim=range(dat), xaxt="n", xlab=xlab, ylab=ylab, ...)
    for (i in 1:nrow(dat)) {
        lines (1:2+x.ori, dat[i,1:2], lwd=.25, col=ifelse(dat[i,1]<=dat[i,2],"red","black"))
        if (is.null(dat2)) {
            lines (2:3+x.ori, dat[i,2:3], lwd=.25, col="lightgreen")
            lines (3:4+x.ori, dat[i,3:4], lwd=.25, col=ifelse(dat[i,3]<=dat[i,4],"black","red"))
        }
    }
    if (!is.null(dat2)) {
        for (i in 1:nrow(dat2)) {
            lines (3:4+x.ori, dat2[i,1:2], lwd=.25, col=ifelse(dat2[i,1]<=dat2[i,2],"black","red"))
        }
    }
    axis(side=1, at=1:4+x.ori, labels=xaxislabels, cex.axis=cex.axis)
}


corplot <- function(object, ...) UseMethod("corplot") 

corplot.default=function(object,y,...){
    dat=data.frame(object,y)
    names(dat)=c("x1", "x2")
    corplot(x2~x1, dat, ...)
}

# col can be used to highlight some points
corplot.formula=function(formula,data,main="",method=c("pearson","spearman"),col=1,cex=.5,add.diagonal.line=TRUE,add.lm.fit=FALSE,col.lm=2,add.deming.fit=FALSE,col.deming=4,add=FALSE,
    log="",same.xylim=FALSE,xlim=NULL,ylim=NULL, ...){
    vars=dimnames(attr(terms(formula),"factors"))[[1]]
    cor.=NULL
    if (length(method)>0) {
        cor.=sapply (method, function (method) {
            cor(data[,vars[1]],data[,vars[2]],method=method,use="p")
        })
        main=main%+%ifelse(main=="","",", ")
        main=main%+%"cor: "%+%concatList(round(cor.,2),"|")
    }

    if (!add) {
        if (same.xylim) {
            xlim=range(model.frame(formula, data))
            ylim=range(model.frame(formula, data))
        }
        plot(formula,data,main=main,col=col,cex=cex,log=log,xlim=xlim,ylim=ylim,...)
    } else {
        points(formula,data,main=main,col=col,cex=cex,log=log,...)
    }
    
    if(add.diagonal.line) abline(0,1)
    if(add.lm.fit) {
        fit=lm(formula, data)
        abline(fit,untf=log=="xy", col=col.lm)
    }
    if(add.deming.fit) {
        # this implementation is faster than the one by Therneau, Terry M.
        fit=Deming(model.frame(formula, data)[[2]], model.frame(formula, data)[[1]]) # this function is in Deming.R copied from MethComp package by Bendix Carstensen
        abline(fit["Intercept"], fit["Slope"], untf=log=="xy", col=col.deming)   
        # Therneau, Terry M.'s implementation in a loose R file that is in 3software folder, slower than Deming, but may be more generalized?
        #fit <- deming(model.frame(formula, data)[[2]], model.frame(formula, data)[[1]], xstd=c(1,0), ystd=c(1,0))
        #abline(fit, untf=log=="xy", col=col.deming)        
    }
    
    cor.
}

abline.pts=function(pt1, pt2=NULL){
    if (is.null(pt2)) {
        if (nrow(pt1)>=2) {
            pt2=pt1[2,]
            pt1=pt1[1,]
        } else {
            stop("wrong input")
        }
    }
    slope=(pt2-pt1)[2]/(pt2-pt1)[1]
    intercept = pt1[2]-slope*pt1[1]
    abline(intercept, slope)
}
#abline.pts(c(1,1), c(2,2))

abline.pt.slope=function(pt1, slope,...){
    intercept = pt1[2]-slope*pt1[1]
    abline(intercept, slope,...)
}
#abline.pt.slope(c(1,1), 1)
mymatplot=function(x, y, type="b", lty=1:5, pch=NULL, col=1:6, xlab=NULL, ylab="", 
    draw.x.axis=TRUE, bg=NA, lwd=1, at=NULL, make.legend=TRUE, legend=NULL, 
    legend.x=9, legend.title=NULL, legend.cex=1, legend.inset=0, ...) {
    
    missing.y=FALSE
    if (missing(y)) {
        missing.y=TRUE
        y=x
        x=1:nrow(y)
    } 
    
    if (is.null(xlab)) xlab=names(dimnames(y))[1]
    if (is.null(legend.title)) legend.title=names(dimnames(y))[2]
    matplot(x, y, lty=lty, pch=pch, col=col, xlab=xlab, xaxt="n", ylab=ylab, bg=bg, lwd=lwd, type=type, ...)
    if(missing.y & draw.x.axis) axis(side=1, at=x, labels=rownames(y)) else if (draw.x.axis) axis(side=1, at=x, labels=x)
    if (make.legend) {
        if (is.null(legend)) legend=colnames(y)
        if (length(unique(pch))>1) {
            mylegend(legend, x=legend.x, lty=lty, title=legend.title, col=col, pt.bg=bg, cex=legend.cex, lwd=lwd, inset=legend.inset, pch=pch)
        } else {
            mylegend(legend, x=legend.x, lty=lty, title=legend.title, col=col, pt.bg=bg, cex=legend.cex, lwd=lwd, inset=legend.inset)
        }
    }
}
