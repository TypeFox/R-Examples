`diagnostic_series` <-
function(
measured,
modelled,
window.size,
step.size=1,
integral_correction=FALSE, use_qualV=FALSE
){
t.pos<- seq(1, (NROW(modelled)-window.size), by=step.size)
dim(t.pos)<-c(length(t.pos),1)

t.diff<-measured-modelled
diff.ecdf <- ecdf(t.diff)

res <- lapply(t.pos, FUN=diagnostic_window, window.size=window.size, measured=measured, modelled=modelled,use_qualV=use_qualV, diff.ecdf=diff.ecdf)
r.diag<-do.call("rbind", res)

if(length(res)!=NROW(r.diag)){
    print("unexpected length(res)!=NROW(r.diag)")
    the.names <- colnames(res[[1]])
    all.elements <- sapply(res, FUN=function(x){all(colnames(x)==the.names)})
    browser()
    ele.comp <- sapply(1:length(res), FUN=function(x){bla <- try(all(res[[x]]==r.diag[x,],na.rm=TRUE));if(inherits(bla,"try-error")) print(paste("Not possible for", x));return(bla)})
}



#multiplicative correction of integrated flow
#time series need the same number of NA
if(integral_correction){
    correctA<-measured
    correctA[is.na(modelled)]<-NA
    correctB<-modelled
    correctB[is.na(correctA)]<-NA

    correct.integral<-sum(correctB,na.rm=TRUE)/sum(correctA,na.rm=TRUE)
    r.diag$I<-r.diag$I*correct.integral
}

toAppend<-data.frame(matrix(nrow=floor(window.size/step.size),ncol=NCOL(r.diag)) )
names(toAppend) <- names(r.diag)
r.diag<-rbind(toAppend,r.diag)

stopifnot(NROW(r.diag)==NROW(modelled))

return(r.diag)

}

