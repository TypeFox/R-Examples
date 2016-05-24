sigora <-
    function(GPSrepo,level,markers=FALSE,queryList=NULL,saveFile=NULL,weighting.method="invhm"){
        ##` GPSrepo:Output of makeGPS
        ##` fn:not used
        ##` queryList: query list
        ##`
        idmap<-get(data(idmap,envir=as.environment(parent.frame())))
        invhm <- function(a, b) {
            0.5 * ((1/a) + (1/b))
        }
        reciprod <- function(a, b) {
            1/(a * b)
        }
        jac <- function(a, b) {
            1/(a + b - 1)
        }
        topov <- function(a, b) {
            1/(apply(as.matrix(cbind(a, b)), FUN = min, MARGIN = 1))
        }
        cosine <- function(a, b) {
            1/(sqrt(a * b))
        }
        noweights <- function(a, b) {
            rep(1, length(a))
        }
        justPUGs <- function(a, b) {
            rep(0, length(a))
        }
        hh<-NULL
        jj<-grep("^L",names(GPSrepo),value=T)
        weights<-0
        ## mapping
        if(length(intersect(queryList,GPSrepo$origRepo[[2]]))==0){
            t1<-which.max(sapply(idmap,function(x)length(intersect(x,GPSrepo$origRepo[[2]]))))
            t2<-which.max(sapply(idmap,function(x)length(intersect(x,queryList))))
            queryList<-idmap[which(idmap[,t2]%in%queryList),t1]
            print(paste("Mapped identifiers from" ,colnames(idmap)[t2]," to ",colnames(idmap)[t1],"..."))
        }
        for(ind in 1:level){
            v1<-GPSrepo[[eval(jj[ind])]]
            com1 <- paste("weights<-", weighting.method, "(v1$degs[v1$gs[v1$GPS[,1]]],
v1$degs[v1$gs[v1$GPS[,2]]])")
            eval(parse(text = com1))
            hh<-rbind(hh,cbind(v1$gs[v1$GPS[,1]],
                               v1$gs[v1$GPS[,2]],
                               v1$ps[v1$GPS[,3]],
                               round(weights*100)/100))
        }
        colnames(hh)<-c("gene1","gene2","pathway","weight")      
        if(markers==TRUE){
            for(ind in 1:level){
                v1<-GPSrepo[[eval(jj[ind])]]
                hh<-rbind(hh,cbind(v1$gs[v1$PU[,1]],
                                   v1$gs[v1$PUG[,1]],
                                   v1$ps[v1$PUG[,2]],
                                   1))
            }
        }
        if(is.null(queryList)){
            ##TODO
        }
        hhd<-(hh[hh[,1]%in%queryList & hh[,2]%in%queryList,,drop=FALSE])
        k1<-(aggregate(as.numeric(hhd[,4]),by=list(hhd[,3]),FUN=sum))
        kN<-(aggregate(as.numeric(hh[,4]),by=list(hh[,3]),FUN=sum))
        sum(kN[,2])
        sum(k1[,2])
        ps<-phyper(k1[,2]-1,kN[match(k1[,1],kN[,1]),2],
                   sum(kN[,2])-kN[match(k1[,1],kN[,1]),2],sum(k1[,2]),lower.tail=F)
        ps<-signif(ps,digits = 4)
        Bonfer<-signif(p.adjust(ps,n=length(unique(hh[,3])),method='bonfer'),digits = 4)
        sp1<-GPSrepo$pathwaydescriptions
        summary_results<-data.frame(k1[,1],sp1[match(k1[,1],sp1[,1]),2],ps,Bonfer,k1[,2],
                                    kN[match(k1[,1],kN[,1]),2],sum(kN[,2]),sum(k1[,2]))
        detailes_results<-hhd
        colnames(summary_results)<-c("pathwy.id","description","pvalues","Bonferroni",
                                     "successes","PathwaySize","N","sample.size")
        summary_results<- summary_results[with(summary_results, order(pvalues)),]
        print(summary_results[which(summary_results$Bonfer<0.01),])
        res<-list()
        res[["summary_results"]]<-summary_results
        res[["detailes_results"]]<-detailes_results
        if(!is.null(saveFile)){
            Genes<-vector(mode='character',length=nrow(summary_results))
            for(i in 1:nrow(summary_results)){
                v11<-getGenes(res,i)
                if("Symbol" %in% colnames(v11)){
                    Genes[i]<-paste(v11$Symbol,sep="" ,collapse=";")}else{
                        Genes[i]<-paste(v11[,1],sep="",collapse=";")
                    }                            
            }
            write.table(cbind(summary_results,Genes),file=saveFile,quote=F,row.names=F,sep='\t')
                    }
            invisible(res)
            ##take a look at weget.cmbi.umcn.nl
        }
