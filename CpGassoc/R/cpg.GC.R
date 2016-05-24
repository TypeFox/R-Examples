cpg.GC <-
function(x) {
      if(class(x)%in% c("cpg","cpg.perm")) {
            holm.values<-p.adjust(x$results$gc.p.value,"holm")
            num.holm<-sum(holm.values<.05,na.rm=TRUE)
            holm.ind<-ifelse(holm.values<.05,TRUE,FALSE)
            if(sum(is.na(x$results[,3]))>0) {
              holm.ind[which(is.na(x$results[,3]))]<-FALSE
                  }
            fdr.method<-x$info$FDR.method
            if (fdr.method=="qvalue") {
                 if (!requireNamespace("qvalue", quietly = TRUE)) {
            stop("qvalue needed for this to work. Please install it.",
                      call. = FALSE)
                      }                   
                fdr.adj<-tryCatch(qvalue::qvalue(x$results$gc.p.value), error = function(e) NULL)
             if(is.null(fdr.adj)) {
                fdr.adj <- tryCatch(qvalue::qvalue(x$results$gc.p.value, pi0.method = "bootstrap"), 
                                      error = function(e) NULL)
                if(is.null(fdr.adj)) {
                  fdr.method="BH"
                }}}
            if(fdr.method!="qvalue") {
              fdr.adj<-p.adjust(x$results$gc.p.value,fdr.method)
              }
            num.fdr<-sum(fdr.adj<.05,na.rm=TRUE)
            beta.col<-nrow(x$results)
            levin<-is.factor(x$indep)
            n1<-coef(x)[1,1]
            gcvalue<-median(n1*ifelse(rep(levin,beta.col),x$results[,2],x$results[,2]**2),na.rm=TRUE)/qchisq(.5,n1)
            gcvalue<-ifelse(gcvalue<1,1,gcvalue)
            if(levin) {
              adj.test.stat<-x$results[,2]/gcvalue
                            }
            else{
                adj.test.stat<-sqrt(x$results[,2]**2/gcvalue)*sign(x$results[,2])
                }
            gc.results<-data.frame(x$results[,1],adj.test.stat,x$results$gc.p.value,holm.ind,fdr.adj,stringsAsFactors=FALSE)
            names(gc.results)<-c("CPG.Labels","GC.Adjusted","Adjust.P.value","Adj.Holm","Adj.FDR")
            gc.info<-data.frame(num.holm=num.holm,FDR.method=fdr.method,num.fdr=num.fdr,gcvalue=gcvalue,stringsAsFactors=FALSE)
            gc.ev<-list(gc.results=gc.results,gc.info=gc.info,coefficients=coef(x))

            if(class(x) %in% "cpg.perm") {
              perm.p.gc<-sum(x$gc.permutation.matrix[,1] <= min(gc.results[,3],na.rm=TRUE))/x$perm.p.values$nperm
              perm.holm.gc<-sum(x$gc.permutation.matrix[,2] >= gc.info[,1])/x$perm.p.values$nperm
              perm.fdr.gc<-sum(x$gc.permutation.matrix[,3] >= gc.info[,3])/x$perm.p.values$nperm
              gc.info<-data.frame(gc.info,perm.p.gc,perm.holm.gc,perm.fdr.gc,stringsAsFactors=FALSE)
              gc.ev<-list(gc.results=gc.results,gc.info=gc.info,coefficients=coef(x))
                             }
              class(gc.ev)<-ifelse(class(x)=="cpg","cpg.gc","cpg.perm.gc")
            gc.ev

            }

          }
