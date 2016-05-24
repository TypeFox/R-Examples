getRpartModel <- function(model, dataset) {
    m <- getRpart(model);
    ee <-list();
    class(ee)<-"rpart";
    if (model$model == "regTree") 
       ee$method <- "anova" 
    else ee$method<-"class";
    ff<-m[[1]];
    ff<-matrix(as.numeric(ff),m[[4]],m[[5]],TRUE);
    dim(m[[3]])<-c(m[[4]],3);
    attrDesc<-matrix(m[[12]],m[[8]],m[[13]], TRUE)
    var<-m[[3]][,2];
	suppressWarnings(idv<-as.integer(var))
    var[!is.na(idv)]<-attrDesc[,2];
    dim(var)<-c(m[[4]],1);
    dfr<-data.frame(var,ff);
    
    desc<-m[[3]][,3]!= ""
    yvalue<-dfr[,9];
    yvalue[desc]<-m[[3]][desc,3]
    classLev<-length(model$class.lev);
    #special case if class has levels - 
    #it is not a regression tree.
    if(classLev > 0){
        stat<-matrix(m[[14]], m[[4]],classLev,TRUE);
        stat<-unlist(stat);
        dim(stat)<-c(m[[4]],classLev);
        stat1<-stat/dfr[,2];
        dim(stat)<-c(1,m[[4]]*classLev);
        dim(stat1)<-c(1,m[[4]]*classLev);
		y2 <- c()
		for (i in 1:length(yvalue)) {
			if (yvalue[i] %in% model$class.lev)
				y2[i] <- as.integer(factor(yvalue[i],levels=model$class.lev))
		    else y2[i] <- as.integer(yvalue[i])
		}
		dfryval2<-c(y2,stat,stat1);
		# dfryval2<-c(yvalue,stat,stat1);
        dim(dfryval2)<-c(m[[4]],2*classLev+1);
		method <- "class"
    }
    else{
        av<-mat.or.vec(m[[4]],1);
        ps<-strsplit(yvalue, " ", fixed=TRUE);
        isOldType<-length(unlist(ps)) == length(yvalue);
        if(isOldType){
            yvalue<-as.numeric(yvalue);
        }
        dfryval2<-c(yvalue);
        dim(dfryval2)<-c(m[[4]],1);
    }
    dfr[,9]<-dfryval2;
    ee$frame<-dfr;
    attr(ee$frame, "row.names")<-as.integer(c(m[[3]][,1]));
    attr(ee$frame, "names")<-m[[2]];
    splits <-matrix(as.numeric(m[[6]]),m[[8]],m[[9]],TRUE);
    splits<-splits[,2:dim(splits)[2]];
    dim(splits)<-c(m[[8]],as.integer(m[[9]])-1);
    splCol<-attrDesc[,2]
    attr(splits, "dimnames")<-list(splCol, m[[7]]);
    #replace count column with values from frame$n
    attrIndex<-ee$frame$var!="<leaf>"
    leftchild<-attr(ee$frame[attrIndex,],"row.names")*2
    i<-1
    while(i<=length(leftchild)){
        selected<-attr(ee$frame,"row.names")==leftchild[i]
        splits[i,1]<-ee$frame$n[selected]
        i<-i+1;
    }
    ee$splits<-splits;
    search <- all.vars(model$formula)[1]
    b<-lapply(attr(dataset, 'names'),function(x){x!=search});
    #we can use levels onto dataset, 
    #because rpart.labels matches attributes by name.
    cs<-lapply(dataset,levels);
    #cs is a list of attributes each of them has a list of values.
    cs<-cs[unlist(b)]
    d<-max(unlist(lapply(cs, length)));
    d<-min(d, as.numeric(m[[11]]));
    if(d > 0){
        csplit<-matrix(as.numeric(m[[10]]),m[[8]],m[[11]],TRUE);
        #csplit==3 is the default, 
        #if min is 1 there is at least one discrete attribute
        if(min(csplit) == 1){
            csplit<-csplit[,1:d];
            dim(csplit)<-c(m[[8]],d);
            ee$csplit<-csplit;
        }
    }
    attr(ee, "xlevels")<-cs;
    attr(ee, "ylevels")<-model$class.lev;
    if(classLev > 0){
        #function used for formatting discrete class
        textFunction <- function (yval, dev, wt, ylevel, 
                digits, n, use.n){
            nclass <- (ncol(yval) - 1L)/2
            group <- yval[, 1L]
            counts <- yval[, 1L + (1L:nclass)]
            dimCounts<-dim(counts);
            counts<-matrix(as.numeric(counts),dimCounts[1],dimCounts[2])
            temp1 <- rpart.formatg(counts, digits)
            if (nclass > 1) {
                temp1 <- apply(matrix(temp1, ncol = nclass), 
                        1, paste, collapse = "/")
            }
            if (use.n) {
                out <- paste(format(group, justify = "left"), 
                        "\n", temp1, sep = "")
            }
            else {
                out <- format(group, justify = "left")
            }
            return(out);
        }
    }
    else{
        #function used for formatting continuous class
        textFunction <- function (yval, dev, wt, ylevel, 
                digits, n, use.n){
            if(!is.numeric(yval)){
                ps<-strsplit(yval, "+", fixed=TRUE);
                for(ij in 1:length(ps)){
                    yval[ij]<-paste(ps[[ij]], collapse = "+\n");
                }
                if (use.n) {
                    yval<-paste(yval, "\nn=", n, sep = "")
                }
                return(yval);
            }
            else{
                if (use.n) {
                    paste(rpart.formatg(yval, digits), 
                            "\nn=", n, sep = "")
                }
                else {
                    paste(rpart.formatg(yval, digits))
                }
            }
        }
    }
    #rpartNamespace <- asNamespace("rpart");
    #environment(textFunction) <- rpartNamespace;
    ee$functions$text <- textFunction;
    ee;
}
getRpart <- function(model) {
    #regression tree
    if(model$model == "regTree"){
        .Call("exportModelRT", as.integer(model$modelID), PACKAGE="CORElearn")
    }
    #classification tree
    else if(model$model == "tree"){
        .Call("exportModelT", as.integer(model$modelID), PACKAGE="CORElearn")
    }
    else{
        stop("The model must be a regresion or a decision tree.");
    }
}
rfProximity <- function(model, outProximity=TRUE){
    if (model$model == "rf"){
        .Call("exportProximity", as.integer(model$modelID),
              as.integer(outProximity==FALSE), PACKAGE="CORElearn")
    }
    else{
        stop("The model must be a random forest.");
    }
}

getVarImportanceCluster <-function(model, cluster){
    modelID <- model$modelID;
    tmp <- .C("exportVarImportanceCluster",
            as.integer(modelID),
            clusterData = as.integer(cluster),
            var = double(model$noNumeric + model$noDiscrete-1),
            NAOK=TRUE,
            PACKAGE="CORElearn")
}
spaceScale <- function(pr, component){
    space<-cmdscale(pr, component, add=TRUE);
}
rfClustering <- function(model, noClusters=4){
    covMatrix<-rfProximity(model, outProximity=F);
    cluster<-cluster::pam(covMatrix, noClusters, diss=TRUE)
}

rfOutliers <- function(model, dataset){
    pr <- rfProximity(model, outProximity=TRUE);
    search <- all.vars(model$formula)[1];
    b<-lapply(attr(dataset, 'names'),function(x){x==search});
    attrClass<-c(1:length(b))[b==TRUE]
    caseIndex<-matrix(c(1:dim(pr)[2]), dim(pr)[1], dim(pr)[2], TRUE)
    el<-caseIndex;
    i<-1;
    prSum<-matrix(0);
    while(i<=dim(pr)[2]){
        rclass<-dataset[el[i,], attrClass]==dataset[i,attrClass]
        prSumMedian<-median(pr[i,rclass])
        prSumSd<-sd(pr[i,rclass]);
        prSum[i]<-sum(pr[i,rclass]^2)^(-1)-prSumMedian
        if(is.na(prSumSd)){
            prSum[i] = 10;
            warning(paste("Element ",i," is the only one in his class. Setting output to 10."))
        }
        else if(is.numeric(prSumSd) && prSumSd>0){
            prSum[i]<-prSum[i]/prSumSd
        }
        i<-i+1
    }
    outliers<-prSum
}
classPrototypes<-function(model, dataset, noPrototypes=10){
    search <- all.vars(model$formula)[1];
    b<-unlist(lapply(attr(dataset, 'names'),function(x){x==search}));
    lev<-levels(dataset[, b]);
    nclass<-length(lev);
    n<-length(dataset[, 1]);
    pre<-predict(model, dataset);
    cluster<-NULL;
    index<-matrix(1:n, n, nclass);
    out<-matrix(0, nclass, noPrototypes);
    for(j in 1:nclass){
        p<-pre$probabilities;
        maxelem<-NULL;
        i<-1;
        while(i<=n && length(p[p!=0])>0 && length(maxelem) < noPrototypes){
            tmp<-index[p==max(p)];
            ntmp<-length(tmp);
            outPos<-list();
            for(k in 1:ntmp){
                outPos[length(outPos)+1]<-c(1:nclass)[lev==dataset[tmp[k],b]]
            }
            outPos<-unlist(outPos);
            #       if(length(outPos==j) > 0)
            #       {
            maxelem<-c(maxelem, tmp[outPos==j]);
            #       }
            p[tmp,]<-0;
            i<-i+1;
        }
        nmaxelem<-length(maxelem);
        if(nmaxelem>=noPrototypes) {
            out[j, ]<-maxelem[1:noPrototypes];
        }
        else{
            out[j, c(1:nmaxelem)]<-maxelem[1:nmaxelem];
        }
        cluster<-c(cluster, (c(1:noPrototypes)*0+j)[out[j,]>0]);
    }
    out<-t(out);
    dim(out)<-c(1, noPrototypes*nclass);
    o<-list();
    o$prototypes<-as.numeric(out[out>0]);
    o$clustering<-as.numeric(cluster);
    o$levels<-as.character(lev);
    impPredictedExample<-o;
}
# the data.frame set is converted to a form such that all the attributes have values between 0 and 1
# this is useful in visualization
varNormalization<-function(md, set){
    #d-discrete, a-attribure, n-names
    column<-length(set[1,]);
    n<-length(set[,1]);
    colPos<-matrix(FALSE, column, column);
    dan<-md$discAttrNames;
    nd<-length(dan);
    ian<-0;
    if(nd>0){
        int<-vector("numeric",nd);
        for(ian in 1:nd)
        {
            search<-dan[ian];
            colPos[ian,]<-unlist(lapply(attr(set, 'names'),function(x){x==search}));
            int[ian]<-1/(length(levels(set[, colPos[ian, ]]))-1);
        }
    }
    nan<-md$numAttrNames;
    nn<-length(nan);
    if(nn > 0){
        offset<-ian;
        mi<-vector("numeric", nn);
        sigma<-vector("numeric", nn);
        maxnorm<-vector("numeric", nn);
        moveup<-vector("numeric", nn);
        for(ian in 1:nn){
            search<-nan[ian];
            tmp<-unlist(lapply(attr(set, 'names'),function(x){x==search}));
            colPos[ian+offset,]<-tmp
            mi[ian]<-as.numeric(mean(set[,tmp]));
            sigma[ian]<-as.numeric(sd(set[,tmp]));
            allcurval<-(set[, tmp]-mi[ian])/sigma[ian];
            moveup[ian]<-min(c(0,allcurval));
            maxnorm[ian]<-max(allcurval-moveup[ian])
        }
    }
    out<-NULL;
    classV<- all.vars(md$formula)[1];
    classV<-unlist(lapply(attr(set, 'names'),function(x){x!=classV}));
    for(ex in 1:n){
        pos<-vector("numeric", column);
        if(nd>0){
            for(da in 1:nd) {
                lev<-levels(set[, colPos[da, ]]);
                val<-set[ex, colPos[da, ]];
                if(is.na(val)){
                    index<-1;
                }
                else{
                    index<-c(1:length(lev))[lev==val];
                }
                pos[colPos[da, ]]<-(index-1)*int[da];
            }
        }
        if(nn>0){
            for(na in 1:nn){
                normal<-((set[ex, colPos[na+offset, ]]-mi[na])/sigma[na]);
                val<-(normal-moveup[na])/maxnorm[na];
                pos[colPos[na+offset, ]]<-val;
            }
        }
        out<-c(out,pos[classV]);
    }
    out<-matrix(out, n, column-1, TRUE);
}

getQuartils<-function(examples){
    int<-0.25;
    i<-1
    minOld<-FALSE
    while(i<5){
        minIndex<-examples<i*int
        aMedian<-median(examples[!minOld & minIndex]);
        examples[!minOld & minIndex]<-aMedian;
        minOld<- minOld | minIndex
        i<-i+1
    }
    getQuartils<-examples;
}

rfAttrEvalClustering<-function(model, dataset, clustering=NULL){
    search <- all.vars(model$formula)[1];
    b<-lapply(attr(dataset, 'names'),function(x){x==search})
    b<-unlist(b);
    i<-1;
    imp<-NULL;
    lev<-NULL;
    if(is.null(clustering)){
        cl2<-as.numeric(dataset[,b]);
        levSet<-dataset[,b];
    }
    else{
        cl2<-clustering;
        levSet<-cl2;
    }
    cl<-cl2;
    stopme<-length(levels(levSet));
    while(i<=stopme){
        cl[cl!=i]<-0;
        cl[cl==i]<-1;
        lev<-c(lev, as.character((levSet[cl==1])[1]));
        d<-getVarImportanceCluster(model, cl);
        imp<-c(imp,d$var);
        cl<-cl2;
        i<-i+1;
    }
    temp<-list();
    ncolumn<-length(b[b==FALSE]);
    imptemp<-matrix(imp, stopme, ncolumn, TRUE);
  	colnames(imptemp) <- names(dataset)[b==FALSE]
	  rownames(imptemp) <- levels(levSet)
    temp$imp<-imptemp
    temp$levels<-lev;
    rfAttrEvalClustering<-temp;
	
}

plotRFStats <- function(point, cluster=FALSE, plotLine=FALSE, 
        lOffset=0, myCount=7, myAxes=FALSE)
{
    pointLen <- length(point);
    if(is.null(dim(point)) || length(dim(point)) == 1)
    {
        tmpPoint <- point;
        point<-matrix(0, pointLen, 2);
        point[,1] <- c(1:pointLen);
        point[,2] <- tmpPoint;
    }
    noVar <- pointLen;
    ylim<-c(min(point[,2]), max(point[,2])+lOffset)
    xlim<-c(min(point[,1]), max(point[,1]))
    if(is.logical(myAxes)){
        axesShow<-TRUE;
    }
    else{
        axesShow<-FALSE;
    }
    plot(1, 1, xlim=xlim, ylim=ylim, type="n", ann=FALSE, frame=TRUE, axes=axesShow);
    if(!is.logical(myAxes) && length(myAxes) > 0){
        axis(2);
        axis(1, at=1:noVar, labels=myAxes);
    }
    if(cluster!=FALSE && length(cluster) > 0)
    {
        tmpCluster <- cluster;
        clusterLevelNames <- list();
        clusterLevels <- 0;
        i <- 1;
        while(length(tmpCluster) > 0 && i < 13)
        {
            clusterLevels[i] <- i;
            clusterLevelNames[i]<-tmpCluster[1]
            cluster[cluster==tmpCluster[1]]<-i;
            tmpCluster <- tmpCluster[tmpCluster!=tmpCluster[1]];
            i <- i+1;
        }
        clusterLevels <- clusterLevels[clusterLevels>0];
        myPch <-0;
        myColor <-0;
        for(value in clusterLevels)
        {
            #mod
            myPch[value]<-floor(value/myCount);
            myColor[value]<- 1+value - myCount*floor(value/myCount);
            
            points(point[cluster==value,], col=myColor[value], pch=myPch[value]);
            if(plotLine == TRUE)
            {
                lines(point[cluster==value,], col=myColor[value], pch=myPch[value]);
            }
        }
        prefix<-""
        if(is.integer(clusterLevelNames[1]))
        {
            prefix<-"skupina "
        }
        clusterNames <- paste(prefix, clusterLevelNames, sep = "");
        legend(xlim[1], ylim[2], clusterNames, cex=0.8, col=myColor, pch=myPch);
    }
    else
    {
        points(point);
        if(plotLine == TRUE)
        {
            lines(point);
        }
    }
}

plotRFMulti<-function(point, legendNames=FALSE, lOffset=0, 
        myCount=7, myHoriz=FALSE, myAxes=FALSE)
{
    noVar<-dim(point)[2]
    noCluster<-dim(point)[1]
    ylim<-c(min(point), max(point)+lOffset)
    xlim<-c(1, noVar)
    if(is.logical(myAxes)){
        axesShow<-TRUE;
    }
    else{
        axesShow<-FALSE;
    }
    plot(1, 1, xlim=xlim, ylim=ylim, type="n", ann=FALSE, frame=TRUE, axes=axesShow);
    if(!is.logical(myAxes) && length(myAxes) > 0){
        axis(2);
        axis(1, at=1:noVar, labels=myAxes);
    }
    myPch <-0;
    myColor <-0;
    pPoint<- matrix(0.0, noVar, 2)
    pPoint[,1]<-c(1:noVar)
    for(i in c(1:noCluster))
    {
        myPch[i]<-floor(i/myCount);
        myColor[i]<- 1+i - myCount*floor(i/myCount);
        
        pPoint[,2]<-point[i,]
        points(pPoint, col=myColor[i], pch=myPch[i]);
        lines(pPoint, col=myColor[i], pch=myPch[i]);
    }
    prefix<-""
    if(is.logical(legendNames) && legendNames==FALSE)
    {
        legendNames<-c(1:noCluster)
        prefix <-"cluster"
    }
    clusterNames <- paste(prefix, legendNames, sep = "");
    color<-c(1:noCluster);
    legend(xlim[1], ylim[2], clusterNames, cex=0.8, col=myColor, pch=myPch, horiz=myHoriz);
}

plotRFNorm<-function(point, cluster, somnames, lOffset, 
        myHoriz=FALSE, myAxes=FALSE)
{
    noVar<-dim(point)[2];
    ylim<-c(min(point), max(point)+lOffset)
    xlim<-c(1, noVar)
    if(is.logical(myAxes)){
        axesShow<-TRUE;
    }
    else{
        axesShow<-FALSE;
    }
    plot(1, 1, xlim=xlim, ylim=ylim, type="n", ann=FALSE, frame=TRUE, axes=axesShow);
    if(!is.logical(myAxes) && length(myAxes) > 0){
        axis(2);
        axis(1, at=1:noVar, labels=myAxes);
    }
    tmpCluster <- cluster;
    clusterLevelNames <- list();
    clusterLevels <- 0;
    i <- 1;
    while(length(tmpCluster) > 0 && i < 13){
        clusterLevels[i] <- i;
        clusterLevelNames[i]<-tmpCluster[1]
        cluster[cluster==tmpCluster[1]]<-i;
        tmpCluster <- tmpCluster[tmpCluster!=tmpCluster[1]];
        i <- i+1;
    }
    clusterLevels <- clusterLevels[clusterLevels>0];
    myPch <-0;
    myColor <-0;
    myCount<-7
    nexamples<-length(point[,1]);
    for(value in 1:nexamples){
        #mod
        myPch[cluster[value]]<-floor(cluster[value]/myCount);
        myColor[cluster[value]]<- 1+cluster[value] - myCount*floor(cluster[value]/myCount);
        points(point[value,], col=myColor[cluster[value]], pch=myPch[cluster[value]]);
        lines(point[value,], col=myColor[cluster[value]], pch=myPch[cluster[value]]);}
    legend(xlim[1], ylim[2], somnames, cex=0.8, col=myColor, pch=myPch, horiz=myHoriz);
}

# taken from rpart package, because it does not export it anymore
## format a set of numbers using C's "g" format
rpart.formatg <- function(x, digits = getOption("digits"),
		format = paste0("%.", digits, "g"))
{
	if (!is.numeric(x)) stop("'x' must be a numeric vector")
	
	temp <- sprintf(format, x)
	if (is.matrix(x)) matrix(temp, nrow = nrow(x)) else temp
}

