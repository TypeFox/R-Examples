join.results<-function(...,type = NULL, genenames = NULL)
{
args<-list(...)
N<-length(args)
if (!(is.null(type)) & N!=length(type)) stop ("Vector type has not correct length")
results<-list()

if (is.null(type)) {
for (i in 1:N) {
if ("metaMA.res" %in% class(args[[i]]) ) {results[[i]]<-metalist.to.matrix(args[[i]],genenames); class(results[[i]])<-"metaMA"}
if ("ES.GeneMeta.res" %in% class(args[[i]]) ) {results[[i]]<- args[[i]]; class(results[[i]])<-"ES.GeneMeta"}
if ("RankProduct.res" %in% class(args[[i]]) ) {results[[i]]<-rbind(args[[i]]$Table1, args[[i]]$Table2); class(results[[i]])<-"RankProduct"}
if ("SOGLresult" %in% class(args[[i]]) ) {
 genenames<-args[[i]]$all.genes
 dum<-as.data.frame(genenames %in% args[[i]]$genes)
 rownames(dum)<-genenames
 results[[i]]<-dum
 class(results[[i]])<-c("data.frame","SOGL")}
if ("posterior.mean"  %in% class(args[[i]]) ) {results[[i]]<-args[[i]]; class(results[[i]])<-c("data.frame","post.mean")}
if ("MAP.Matches.res" %in% class(args[[i]]) ) {results[[i]]<-probs.to.matrix(args[[i]]$genes, args[[i]]$all.genes);class(results[[i]])<-"MAP.Matches"}
if ("METRADISC.res"   %in% class(args[[i]]) ) {results[[i]]<-args[[i]];class(results[[i]])<-"METRADISC"}
}
} else {
if (is.null(genenames)) stop("The 'genenames' argument can not be missing.") 
for (i in 1:N)
{
if (type[i]==1) {results[[i]]<-metalist.to.matrix(args[[i]],genenames)}
if (type[i]==2) {dum<-as.data.frame(genenames %in% args[[i]]$genes)
 rownames(dum)<-genenames
 results[[i]]<-dum
names(results[[i]])<-"SOGL"}
if (type[i]==3) {results[[i]]<-rbind(args[[i]]$Table1, args[[i]]$Table2)}
if (type[i]==4) {results[[i]]<-probs.to.matrix(args[[i]], genenames)}
if (type[i]==5) {results[[i]]<-args[[i]]}
}
}
return(results)
}

