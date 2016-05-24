make.panel.svysmooth<-function(design,bandwidth=NULL){
  function(x,y,span=NULL,col.smooth="red",col=par("col"),bg=NA,pch=par("pch"),cex=1,...){
    if(!is.null(span))
      bandwidth<-diff(range(x))*span/3
    s<-svysmooth(y~x,design=design,bandwidth=bandwidth)
    points(x,y,pch=pch,bg=bg,col=col)
    lines(s[[1]],col=col.smooth,...)
  }
}


svyplot<-function(formula, design,...) UseMethod("svyplot",design)
svyplot.default<-function(formula,
                  design,
                  style=c("bubble","hex","grayhex","subsample","transparent"),
                  sample.size=500, subset=NULL,legend=1,inches=0.05,
                  amount=NULL,basecol="black",alpha=c(0,0.8), xbins=30,...){
  
  style<-match.arg(style)
  if (style %in% c("hex","grayhex") && !require(hexbin)){
    stop(style," plots require the hexbin package")
  }

  subset<-substitute(subset)
  subset<-with(design$variables, subset)
  if(length(subset)>0)
    design<-design[subset,]

  W<-weights(design, "sampling")

  mf<-model.frame(formula, design$variables,na.action=na.pass)  
  Y<-model.response(mf)
  X<-mf[,attr(attr(mf,"terms"),"term.labels")]
  
  switch(style, 
         bubble={
           if(is.function(basecol)) basecol<-basecol(model.frame(design))
           symbols(X,Y,circles=sqrt(W),inches=inches,fg=basecol,...)
         },
         hex={
           ## CRAN will be happier if we stop supporting the old version of hexbin
             ## new version
             rval<-hexbin(X,Y,IDs=TRUE,xbins=xbins)
             cell<-rval@cID
             rval@count<-as.vector(tapply(W,cell,sum))
             rval@xcm<-as.vector(tapply(1:length(X), cell,
                              function(ii) weighted.mean(X[ii],W[ii])))
             rval@ycm<-as.vector(tapply(1:length(Y), cell,
                              function(ii) weighted.mean(Y[ii],W[ii])))
             gplot.hexbin(rval, legend=legend, style="centroids",...)
           
           
         },
         grayhex={
             ## new version
             rval<-hexbin(X,Y,IDs=TRUE,xbins=xbins)
             cell<-rval@cID
             rval@count<-as.vector(tapply(W,cell,sum))
             gplot.hexbin(rval, legend=legend,...)
        
         },
         subsample={
           index<-sample(length(X),sample.size,replace=TRUE, prob=W)
           if (is.numeric(X))
             xs<-jitter(X[index],factor=3,amount=amount$x)
           else
             xs<-X[index]
           if (is.numeric(Y))
             ys<-jitter(Y[index],factor=3,amount=amount$y)
           else
             ys<-Y[index]
           plot(xs,ys,...)
         },
         transparent={
           transcol<-function(base,opacity){
             rgbs<-col2rgb(base)/255
             rgb(rgbs[1,],rgbs[2,], rgbs[3,], alpha=opacity)
           }
           if(is.function(basecol)) basecol<-basecol(model.frame(design))
           w<-weights(design)
           maxw<-max(w)
           minw<-0
           alphas<- (alpha[1]*(maxw-w)+alpha[2]*(w-minw))/(maxw-minw)
           plot(X,Y,col=transcol(basecol,alphas),...)
         })

}

svyboxplot<-function(formula, design,all.outliers=FALSE,...) UseMethod("svyboxplot",design)
svyboxplot.default<-function(formula, design,  all.outliers=FALSE,col=NULL,names,...){
    
    formula<-as.formula(formula)
    if(length(formula)!=3) stop("need a two-sided formula")
    ##if(length(formula[[3]])>2) stop("only one rhs variable allowed")
    
    outcome<-eval(bquote(~.(formula[[2]])))
    outcome.values<-model.frame(outcome, model.frame(design),na.action=na.pass)
    
    if (length(attr(terms(formula),"term.labels"))){
        groups<-eval(bquote(~.(formula[[3]])))
        qs <- svyby(outcome,groups,design,svyquantile,ci=FALSE,
                    keep.var=FALSE,
                    quantiles=c(0,0.25,0.5,0.75,1),na.rm=TRUE)
        group.values<-model.frame(groups, model.frame(design),na.action=na.pass)[[1]]
        n<-NCOL(qs)
        iqr<- qs[,n-1]-qs[,n-3]
        low<-pmax(qs[,n-4],qs[,n-2]-1.5*iqr)
        hi<-pmin(qs[,n],qs[,n-1]+1.5*iqr)
        stats<-t(as.matrix(cbind(low,qs[,n-(3:1)],hi)))
        z<-list(stats=stats,n=coef(svytotal(groups,design,na.rm=TRUE)))
        for(i in 1:ncol(stats)){
            out<-c(if(qs[i,n]!=hi[i]) qs[i,n],
                   if(qs[i,n-4]!=low[i])qs[i,n-4])
            if (all.outliers){
              outlo<-sort(outcome.values[!is.na(outcome.values) & (as.numeric(group.values) %in% i) & outcome.values<low[i] ])
              outhi<-sort(outcome.values[!is.na(outcome.values) & (as.numeric(group.values) %in% i) & outcome.values>hi[i] ])
              out<-na.omit(unique(c(outlo,outhi)))
            }
            z$out<-c(z$out,out)
            z$group<-c(z$group,rep(i,length(out)))
            z$names<-as.character(qs[,1])
        }
    } else {
        qs<-svyquantile(outcome,design,ci=FALSE,
                        quantiles=c(0,0.25,0.5,0.75,1),na.rm=TRUE)
        iqr<-qs[4]-qs[2]
        z<-list(stats=matrix(c(max(qs[1],qs[2]-1.5*iqr),
                qs[2:4],min(qs[5],qs[4]+1.5*iqr))),
                n=sum(weights(design,"sampling")))
        z$out<-c(if(qs[5]!=z$stats[5]) qs[5],
                 if(qs[1]!=z$stats[1]) qs[1])
        if (all.outliers){
          outlo<-sort(outcome.values[!is.na(outcome.values) &  outcome.values<qs[2]-1.5*iqr ])
          outhi<-sort(outcome.values[!is.na(outcome.values) & outcome.values>qs[4]+1.5*iqr])
          z$out<-na.omit(unique(c(outlo,outhi)))
        }
        z$group<-rep(1,length(z$out))
    }
    if (is.null(col)) col<-par("bg")
    if (!missing(names)) z$names<-names
    bxp(z,boxfill=col,...)
}



svycoplot<-function(formula, design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),...) UseMethod("svycoplot",design)
svycoplot.default<-function(formula, design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),xbins=15,...){
  require(lattice)
  style<-match.arg(style)
  wt<-weights(design,"sampling")
  
  switch(style,
         hexbin={
           require(hexbin) || stop("hexbin package is required (from Bioconductor)")
           hexscale<-match.arg(hexscale)
           xyplot(formula, data=model.frame(design), xbins=xbins,
                  panel=function(x,y,style="centroids",xbins,subscripts,...) {
                    if (!length(x)) return(panel.xyplot(x,y,...))
                    vp<-current.viewport()
                    wd<-convertWidth(vp$width,unitTo="cm",valueOnly=TRUE)
                    ht<-convertHeight(vp$height,unitTo="cm",valueOnly=TRUE)
                    W<-wt[subscripts]
                    rval<-hexbin(x,y,IDs=TRUE,xbins=xbins,shape=ht/wd,xbnds=vp$xscale,ybnds=vp$yscale)
                    cell<-rval@cID
                    rval@count<-as.vector(tapply(W,cell,sum))
                    rval@xcm<-as.vector(tapply(1:length(x), cell,
                                               function(ii) weighted.mean(x[ii],W[ii])))
                    rval@ycm<-as.vector(tapply(1:length(y), cell,
                                               function(ii) weighted.mean(x[ii],W[ii])))
                    grid.hexagons(rval,style=style, maxarea=switch(hexscale, relative=0.8,
                                                      absolute=0.8*sum(W)/sum(wt)))
                  },...)
         }, transparent={
           if(is.function(basecol)) basecol<-basecol(model.frame(design))
           transcol<-function(base,opacity){
             rgbs<-col2rgb(base)/255
             rgb(rgbs[1,],rgbs[2,], rgbs[3,], alpha=opacity)
           }
           maxw<-max(wt)
           minw<-0
           alphas<- (alpha[1]*(maxw-wt)+alpha[2]*(wt-minw))/(maxw-minw)
           cols<-transcol(basecol,alphas)
           xyplot(formula, data=model.frame(design), 
                  panel=function(x,y,basecol="black",subscripts,...) {
                    a<-alphas[subscripts]
                    panel.xyplot(x,y,col=cols[subscripts],pch=19,...)
                  },...)
         }
         )
}


barplot.svystat<-function(height,...) barplot(coef(height),...)
barplot.svrepstat<-function(height,...) barplot(coef(height),...)

plot.svystat<-function(x,...) barplot(coef(x),...)
plot.svrepstat<-function(x,...) barplot(coef(x),...)

barplot.svyby<-function(height,beside=TRUE,...){    
  aa <- attr(height, "svyby")
  rval <- height[, max(aa$margins) + (1:aa$nstats)]
  if (is.null(dim(rval))) {
    if (length(aa$margins)<2){
      names(rval) <- row.names(height)
    } else {
      rval<-matrix(rval, nrow=length(unique(height[,aa$margins[1]])))
      rownames(rval) <- unique(height[,aa$margins[1]])
      colnames(rval)<-levels(do.call(interaction, height[,aa$margins[-1],drop=FALSE]))
    }
  } else {
    rval <- as.matrix(rval)
    colnames(rval)<-sub("statistics\\.","",colnames(rval))
    rval<-t(rval)
  }
  barplot(rval,beside=beside,...)
}

plot.svyby<-function(x,...) barplot.svyby(x,...)


dotchart.default<-graphics::dotchart
dotchart<-function(x,...,pch=19) UseMethod("dotchart")
dotchart.svystat<-function(x,...,pch=19) dotchart(coef(x),...,pch=pch)
dotchart.svrepstat<-function(x,...,pch=19) dotchart(coef(x),...,pch=pch)

dotchart.svyby<-function(x,...,pch=19){
  height<-x
  aa <- attr(height, "svyby")
  rval <- height[, max(aa$margins) + (1:aa$nstats)]
  if (is.null(dim(rval))) {
    if (length(aa$margins)<2){
      names(rval) <- row.names(height)
    } else {
      rval<-matrix(rval, nrow=length(unique(height[,aa$margins[1]])))
      rownames(rval) <- unique(height[,aa$margins[1]])
      colnames(rval)<-levels(do.call(interaction, height[,aa$margins[-1],drop=FALSE]))
    }
  } else {
    rval <- as.matrix(rval)
    colnames(rval)<-sub("statistics\\.","",colnames(rval))
    rval<-t(rval)
  }
  dotchart(rval,...,pch=pch)
}
