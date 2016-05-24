mod_probs <-
function(object,n.burnin=0,scale=0.1,best=NULL,thin=1){

if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(scale<0 | scale>1){
stop("scale should be between 0 and 1")}
if(!is.null(best)){
if(best<=0){
stop("best should be positive")}}

if(n.burnin>0){
innerBETA<-object$BETA[-(1:n.burnin),]
innerMODEL<-object$MODEL[-(1:n.burnin)]} else{
innerBETA<-object$BETA
innerMODEL<-object$MODEL}

n.sample<-length(innerMODEL)
every<-seq(from=thin,to=n.sample,by=thin)
innerBETA<-innerBETA[every,]
innerMODEL<-innerMODEL[every]

tab<-sort(table(innerMODEL)/length(innerMODEL),decreasing=TRUE)

tab.names<-dimnames(tab)[[1]]

if(is.null(best)){
ref.tab<-tab[tab>tab[1]*scale]} else{
if(best<=length(tab)){
ref.tab<-tab[1:best]} else{
ref.tab<-tab}}

if(length(ref.tab)>1){
ref.tab.names<-dimnames(ref.tab)[[1]]
dimnames(ref.tab)<-NULL} else{
ref.tab.names<-names(ref.tab)
names(ref.tab)<-NULL}

if(length(ref.tab)>1){
forms<-c()
for(j in 1:length(ref.tab)){
forms[j]<-paste0("~",as.character(index2formula(index=model2index(ref.tab.names,dig=dim(innerBETA)[2])[j,],maximal.mod=object$maximal.mod))[3])}
} else{
forms<-paste0("~",as.character(index2formula(index=model2index(ref.tab.names,dig=dim(innerBETA)[2]),maximal.mod=object$maximal.mod))[3])}

est<-list(table=data.frame(prob=ref.tab,model_formula=forms),totmodsvisit=length(tab),thin=thin)

class(est)<-"modprobs"

est}
