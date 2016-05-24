sim2.bd.mrca <- function(mrca,numbsim,lambda,mu,K){
	phy2 <- list()
	if (K==0){
	for (j in 1:numbsim) {
		
		stop = 0
		while (stop==0){
			t1 <- sim2.bd.origin(0,mrca,lambda,mu,K)
			if (class(t1) == "phylo"){
				stop =1
			} else if (t1 == 1) {
				stop = 1	
			}
		}
		stop = 0
		while (stop==0){
			t2 <- sim2.bd.origin(0,mrca,lambda,mu,K)
			if (class(t2) == "phylo"){
				stop =1
			} else if (t2 == 1) {
				stop = 1	
			}
		}		
		
		if (class(t1) == "phylo" && class(t2) == "phylo"){
			#altered on July 8, 2014
			t1$root.edge<-t1$edge.length[1]
			t2$root.edge<-t2$edge.length[1]
			t1<-collapse.singles(t1)			
			t2<-collapse.singles(t2)	
			t<-t1+t2
			#t<- bind.tree(t1,t2,where="root")
			t$tip.label <- paste("t", sample(t$Nnode+1), sep = "")
		} else if (class(t1) == "phylo" || class(t2) == "phylo")  {
			if (class(t2) == "phylo") {
				t1<-t2
				}
			t<-t1
			t$edge <- t1$edge + 1
			root<- t$edge[1,1]
			t$edge <- rbind(c(root,1),t$edge)
			t$edge.length <- c(mrca,t$edge.length)
			t$tip.label <- paste("t", sample(t$Nnode+1), sep = "")
		} else {
			t<- sim2.bd.fast(2,1,1,0,1)[[1]][[1]]
			t$edge.length <- t$edge.length*0+mrca
		} 
		phy2 <- c(phy2,list(t))
	}} else {
		#K>0
		n<-0
		age<-mrca
	lambda0<-lambda
edge <- c(-1,-2)		#matrix of edges
edge<-rbind(edge,c(-1,-3))
leaves <- c(-2,-3)			#list of extant leaves
timecreation <-c(0,0,0)		#time when species -2, -3, ... -n was created after origin
extinct <- vector()		#list of extinct leaves
time <-0				#time after origin
maxspecies <- -3		#smallest species
edge.length <- c(0,0)		#edge length. if 0: leaf which didn't speciate /extinct yet
extincttree = 0
stop = 0

while (stop == 0 ){
if (length(leaves) == 0){
		#phy2 = 0
		if (age >0) {
			phy2=0
			}
		extincttree=1
		stop = 1
		} else {	

if (K>0) {lambda<-lambda0*(1-length(leaves)/K)
	if (lambda<0) {lambda<-0}
	}			
if (mu==0 && lambda==0) {
	#print("no events can occur any more as K reached. no death.")
	stop<-1
	} else {
					
timestep <- rexp(1,(length(leaves)*(lambda+mu)))    #time since last event
if ((age > 0 && (time+timestep) < age) || age == 0) {

time = time+timestep			#time after origin
species <- sample(leaves,1)  	#the leaf undergoing the next event
del <- which(leaves == species)
specevent <- runif(1,0,1)		#event speciation or extinction
edgespecevent <- which(edge == species) - length(edge.length)
if ((lambda/(lambda+mu)) > specevent) {
	edge.length[edgespecevent] <- time-timecreation[- species]
	edge <- rbind(edge,c(species,maxspecies-1))
	edge <- rbind(edge,c(species,maxspecies-2))
	edge.length <- c(edge.length,0,0)
	leaves <- c(leaves,maxspecies-1,maxspecies-2)
	maxspecies <- maxspecies-2
	leaves <- leaves[- del]
	timecreation <- c(timecreation,time,time)
	if (length(leaves) == n){
		stop = 1
		}
	} else {
	extinct <- c(extinct,leaves[del])
	leaves <- leaves[- del]
	edge.length[edgespecevent] <- time-timecreation[- species]
	}
	} else {
	stop = 1	
		}
}
}}

if (extincttree==0){	#length pendant edge
if (age==0){
	timestep <- rexp(1,(length(leaves)*(lambda+mu)))    #time since last event
	time = time+timestep			#time after origin
} else {
	time = age
	}
}


if (extincttree==0 || age ==0){
#assign pendant edge length
for (j in (1:length(leaves))){
	k = which( edge == leaves[j]  ) - length(edge.length)
	edge.length[ k ] <- time - timecreation[- leaves[j]]
	}

nodes <- (length(leaves)+length(extinct))*2
leaf=1
interior=length(leaves)+length(extinct)+1
edgetemp<-edge

if (nodes == 2) {
	#edge = c(2,1)
	#leaves = c(1,2)
	phy2 <-1
	} else {
		
			
for (j in (1:nodes)){
	if (sum(match(leaves,- j,0)) + sum(match(extinct,- j,0)) == 0) {
		# replace all -j values in edge by interior
		posvalue <- interior
		interior <- interior +1
		} else {
		posvalue <- leaf
		leaf <- leaf +1
		}
	replacel <- which(edge == - j)
	if (length(replacel)>0)  {
	for (k in 1:length(replacel)) {
		if ((replacel[k]-1) < length(edge.length)) {
			edge[replacel[k],1] <- posvalue
		} else {
			edge[(replacel[k]-length(edge.length)),2] <- posvalue
			}
	}
	}
}


phy <- list(edge = edge)
phy$tip.label <- paste("t", sample(length(leaves)+length(extinct)), sep = "")
phy$edge.length <- edge.length
phy$Nnode <- length(leaves)+length(extinct)-1
class(phy) <- "phylo"
phy2 <- c(phy2,list(phy))
}
}		
}  #end K>0
phy2	
}	
