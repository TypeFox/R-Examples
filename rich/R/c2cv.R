c2cv <-
function(com1,com2,nrandom=99,pr1=0.025,pr2=0.975,verbose=TRUE)
 {
    SRobs <- function(D, d){E=D[d,]
    return(specpool(E)$Species)}

    n.groupe1 <- dim(com1)[1];n.groupe2 <- dim(com2)[1]
    obs.S1<-0 ; obs.S2<-0  
    obs.S1<-SRobs(com1) ; obs.S2<-SRobs(com2) 
    obs_dif1_2 <- obs.S1 - obs.S2       
    observation.vec<-rbind(com1,com2)
    GrandN<-n.groupe1 + n.groupe2 
    vecteur_numero_ordre<-c(1:GrandN)
    mean.scores.simulated1 <-vector(length=nrandom)
    mean.scores.simulated2 <-vector(length=nrandom)
    mean.scores.simulated_dif1_2 <-vector(length=nrandom)
    for (random in 1:(nrandom)) { 
	vecteur_ordre <- sample(vecteur_numero_ordre,GrandN,replace=FALSE)
	vecteur_ordrep1<-vecteur_ordre[1:n.groupe1]
	vecteur_ordrep2<-vecteur_ordre[(n.groupe1+1):(n.groupe1+n.groupe2)]
	vecteur.random.1<-observation.vec[vecteur_ordrep1,]
	vecteur.random.2<-observation.vec[vecteur_ordrep2,]
	rand.S1<-SRobs(vecteur.random.1) ; rand.S2<-SRobs(vecteur.random.2)
	mean.scores.simulated1[random] <- rand.S1
	mean.scores.simulated2[random] <- rand.S2}
    mean.scores.simulated_dif1_2 <- mean.scores.simulated1 - mean.scores.simulated2 
    mean.scores<-c(obs_dif1_2,mean.scores.simulated_dif1_2[1:nrandom])
    inf<-mean.scores[mean.scores<obs_dif1_2]
    sup<-mean.scores[mean.scores>obs_dif1_2]
    egal<-mean.scores[mean.scores==obs_dif1_2] 
    outtable<-as.data.frame(matrix(ncol=1,nrow=8));names(outtable)<-" " 
    row.names(outtable)<- c("cv1","cv2", "cv1-cv2", "p", paste("quantile",pr1,sep=" ")
    ,paste("quantile",pr2,sep=" "), "randomized cv1-cv2", "nrandom")   
    outtable[1,1]<- obs.S1
    outtable[2,1]<- obs.S2
    outtable[3,1]<-obs_dif1_2   
    if(obs_dif1_2>0) 
	outtable[4,1]<-(length(sup)+length(egal))/(nrandom+1)
    if(obs_dif1_2<0)
	outtable[4,1]<-(length(inf)+length(egal))/(nrandom+1)
    if(obs_dif1_2==0)
	outtable[4,1]<-"NC"
    outtable[5,1]<-quantile(mean.scores,probs = c(pr1, pr2),names=FALSE)[1]
    outtable[6,1]<-quantile(mean.scores,probs = c(pr1, pr2),names=FALSE)[2]
    outtable[7,1]<-mean(mean.scores)
    outtable[8,1]<-nrandom
    
    if(verbose==TRUE){
	out<-list(res=outtable, rand=mean.scores)}
    if(verbose==FALSE){
	out<-list(res=outtable)}    
return(out)
}
