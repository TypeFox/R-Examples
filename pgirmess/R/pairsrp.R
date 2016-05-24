"pairsrp" <-
function(dataframe,meth="spearman",pansmo=FALSE,abv=FALSE,lwt.cex=NULL,...){

# Giraudoux - 1.11.04
# dataframe = data.frame
# draw a scatterplot + regression line in the upper triangle of a matrix
# of every pairwise variables of a dataframe
# and gives the r and p values in the lower triangle
# dataframe = a data.frame of numerics
# meth = "pearson", "spearman" or "kendal"
# abv = should the diag labels be abbreviated?
# lwt.cex = character size expansion in the lower triangle

    if (abv) nm<-abbreviate(names(dataframe)) else nm<-names(dataframe)
    tit<-cor.test(dataframe[,1],dataframe[,2],method=meth)$method
# if (pansmo)  lines(lowess(x,y,...)
# if (pansmo)   lines(sort(x),predict(loess(y~x),data=sort(x),type=response),col="green")
    panel1<-function(x,y,...){points(x,y,...);abline(lm(y~x),col="blue");if (pansmo) lines(lowess(x,y,...),col="green")}
    panel2<-function(x,y,...) {
        res<-cor.test(x,y,method=meth)
        corre<-round(res$estimate,2);cort<-round(res$p.value,3)
        met<-res$method
        resf<-paste("\nr = ",format(abs(corre)*-1,nsmall=2),"\np=",format(cort,nsmall=3),sep="")
        res<-paste("\nr = ",corre,"\np=",cort,sep="")
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        if(cort<=0.05) cols="red" else if (cort<0.1) cols="black" else cols=grey(0.7)
        if (is.null(lwt.cex)) lwt.cex<-0.8/strwidth(resf)
        text(0.5, 0.6, res,col=cols,cex=lwt.cex)}

    pairs(dataframe, labels= nm, panel = points, cex.labels=min(0.9/strwidth(nm)), main=tit,gap=0,
            upper.panel= panel1,
            lower.panel = panel2,...)
            
}
