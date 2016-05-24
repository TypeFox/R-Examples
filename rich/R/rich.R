rich <-
function(matrix, verbose=FALSE, nrandom=NULL)
 {
    SRobs <- function(D, d){E=D[d,]
      return(specpool(E)$Species)}
    ASRobs <- function(D, d){E=D[d,]
	vec.rich<-apply(X=E, MARGIN=1,FUN=sum)
	return(mean(vec.rich))}
    
    bspfm <- function (a){
	a<-as.matrix(a)
	matrice.AS.boot<-boot(a, ASRobs, R=nrandom)
	matrice.AS.boot.ci<-boot.ci(matrice.AS.boot, conf=0.95, type="norm")
	s.moy.boot<-mean(matrice.AS.boot$t)
	s.boot.corr<-2*matrice.AS.boot$t0 - s.moy.boot
	bias<-s.moy.boot-matrice.AS.boot$t0
	stddev.valboot<-sd(as.vector(matrice.AS.boot$t)) 
	vec.boot.res<-vector(length=7)
	vec.boot.res[1]<-matrice.AS.boot$t0
	vec.boot.res[2]<-s.moy.boot
	vec.boot.res[3]<-s.boot.corr
	vec.boot.res[4]<-bias
	vec.boot.res[5]<-stddev.valboot
	if(is.numeric(matrice.AS.boot.ci$normal[2])==TRUE) 
	    vec.boot.res[6]<-matrice.AS.boot.ci$normal[2]
	else vec.boot.res[6]<-NA
	if(is.numeric(matrice.AS.boot.ci$normal[3])==TRUE)
	    vec.boot.res[7]<-matrice.AS.boot.ci$normal[3]
	else vec.boot.res[7]<-NA
	vec.boot.res.sm<-as.data.frame(t(vec.boot.res));names(vec.boot.res)<-c("mr.obs","mr.boot","mr.bcorr","mr.bias","mr.se",
	"mr.lbn","mr.ubn")
    return(vec.boot.res)}

    bspf <- function (a){
	a<-as.matrix(a)
	matrice.S.boot<-boot(a, SRobs, R=nrandom)
	matrice.S.boot.ci<-boot.ci(matrice.S.boot, conf=0.95, type="norm")
	s.moy.boot<-mean(matrice.S.boot$t)
	s.boot.corr<-2*matrice.S.boot$t0 - s.moy.boot
	bias<-s.moy.boot-matrice.S.boot$t0
	stddev.valboot<-sd(as.vector(matrice.S.boot$t)) 
	vec.boot.res<-vector(length=7)
	vec.boot.res[1]<-matrice.S.boot$t0
	vec.boot.res[2]<-s.moy.boot
	vec.boot.res[3]<-s.boot.corr
	vec.boot.res[4]<-bias
	vec.boot.res[5]<-stddev.valboot
	if(is.numeric(matrice.S.boot.ci$normal[2])==TRUE) 
	    vec.boot.res[6]<-matrice.S.boot.ci$normal[2]
	else vec.boot.res[6]<-NA
	if(is.numeric(matrice.S.boot.ci$normal[3])==TRUE)
	    vec.boot.res[7]<-matrice.S.boot.ci$normal[3]
	else vec.boot.res[7]<-NA
	vec.boot.res<-as.data.frame(t(vec.boot.res),row.names = "")
	names(vec.boot.res)<-c("cr.obs","cr.boot","cr.bcorr","cr.bias",
				"cr.se","cr.lbn","cr.ubn")
    return(vec.boot.res)}
	cat("computing...", "\n")
	matrix<-as.matrix(matrix)
	richcum<-SRobs(matrix)   
	b<-as.vector(matrix) ; b[which(b>=1)]<-1
	a<-matrix(b, nrow=dim(matrix)[1], ncol=dim(matrix)[2])
	colnames(a)<-colnames(matrix)
	vec.rich<-apply(X=a, MARGIN=1,FUN=sum)
	vec.abond<-apply(X=matrix, MARGIN=2,FUN=sum)
	singletons<-length(vec.abond[vec.abond==1])    
	doubletons<-length(vec.abond[vec.abond==2])     
    unique <- function(v){
	ind<-length(v)-length(v[v==0])
	if(ind==1) res<-1 
	else
	res<-0
    return(res)}

    uniques<-sum(apply(X=matrix, MARGIN=2,FUN=unique))

    duplicate <- function(v){
	ind<-length(v)-length(v[v==0])
	if(ind==2) res<-1 
	else
	res<-0
    return(res)} 

    duplicates<-sum(apply(X=matrix, MARGIN=2,FUN=duplicate))
	richmoy<-mean(vec.rich)
	richsd<-sd(as.vector(vec.rich))
	if(is.null(nrandom)==FALSE) {
	    if(nrandom<10) {nrandom<-99 ; resboot<-bspf(matrix)
	      resboot2<-data.frame(t(bspfm(matrix)),row.names = "")}
	    else {resboot<-bspf(matrix); resboot2<-data.frame(t(bspfm(matrix)),
		  row.names = "")}
	  }
	else {resboot<-NULL ; resboot2<-NULL}
    
    if(verbose==TRUE) 
	out<-list(cr=richcum, mr=richmoy, mrsd=richsd,
	singletons=singletons, doubletons=doubletons, uniques=uniques, 
	duplicates=duplicates, bootCR=resboot, bootMR=resboot2, nrandom=nrandom, 
	richvec=vec.rich, matrix=matrix, matrixbin=a,
	sumrow=apply(X=a, MARGIN=1,FUN=sum), sumcol=apply(X=a, MARGIN=2,FUN=sum),
	zeroes=length(which(as.vector(a)==0)))
    else
	out<-list(cr=richcum, mr=richmoy, mrsd=richsd,
	bootCR=resboot,bootMR=resboot2,nrandom=nrandom)

return(out)
}
