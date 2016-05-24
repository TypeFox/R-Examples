sim.bdtypes.stt.taxa <- function(n, lambdavector,deathvector,sampprobvector,init=-1,EI=FALSE,eliminate=0){

	muvector<-deathvector*(1-sampprobvector)
	psivector<-deathvector*sampprobvector
	
	extincttree = 1

	if ((init == -1) && (length(deathvector)==2)){
		init<-2
		lamb<-lambdavector[1,1]-lambdavector[2,2]-deathvector[1]+deathvector[2]
		c<- sqrt(lamb^2 + 4*lambdavector[1,2]*lambdavector[2,1])
		f1 <- (c + lamb)/(c+lamb+2*lambdavector[1,2])
		r<-runif(1,0,1)
		if (r<f1){init <-1}
	}
	if ((init==-1) && (length(deathvector)!=2))  {init<-sample(1:length(deathvector),1)}

	if (EI==TRUE){init<-1}

	while(extincttree==1){
	edge <- c(-1,-2)		#matrix of edges
	leaves <- c(-2)			#list of extant leaves
	types <- c(init)		#list of types of extant leaves
	sampled <- vector()
	typessampled<-vector()
	timecreation <-c(0,0)	#time when species -1 (origin), -2, -3, ... -n was created after origin
	extinct <- vector()		#list of extinct leaves
	time <-0				#time after origin
	maxspecies <- -2		#smallest species
	edge.length <- c(0)		#edge length. if 0: leaf which didn't speciate /extinct yet
	extincttree = 0
	stop = 0

	while (stop == 0 ){
		if (length(leaves) == 0){
			phy2=0
			extincttree=1
			print("extinct tree")
			stop = 1
		} else {	
			sumrates<-vector()
			for (i in 1:length(lambdavector[,1])){
				sumrates<-c(sumrates,(length(which(types==i))*(sum(lambdavector[i,])+muvector[i]+psivector[i])))
			}
			timestep <- rexp(1,sum(sumrates))    #time since last event
			time = time+timestep			#time after origin
			r<-runif(1,0,sum(sumrates))
			chosentype<-min(which(cumsum(sumrates)>r))
			species <- sample(leaves[which(types==chosentype)],1)  	#the leaf undergoing the next event
			lambda<-sum(lambdavector[chosentype,])
			#7.9.14 EI model
			gamma <- 0
			if (EI==TRUE){
			if (chosentype == 1) {gamma<-lambda
				lambda <-0
				}}
			mu<-muvector[chosentype]
			psi<-psivector[chosentype]
			del <- which(leaves == species) #leaves[del]=species
			specevent <- runif(1,0,1)		#event speciation, sampling or extinction
			edgespecevent <- which(edge == species) - length(edge.length)
			if ((lambda/(lambda+gamma+mu+psi)) > specevent) {
				edge.length[edgespecevent] <- time-timecreation[- species]
				edge <- rbind(edge,c(species,maxspecies-1))
				edge <- rbind(edge,c(species,maxspecies-2))
				edge.length <- c(edge.length,0,0)
				r<-runif(1,0,lambda)
				newtype<-min(which(cumsum(lambdavector[chosentype,])>r))
				leaves <- c(leaves,maxspecies-1,maxspecies-2)
				types<- c(types, chosentype, newtype)
				maxspecies <- maxspecies-2
				leaves <- leaves[- del]
				types <- types[- del]
				timecreation <- c(timecreation,time,time)}
			else if (((lambda+gamma)/(lambda+gamma+mu+psi)) > specevent) {	#EI
				types[del] <- 2
			}
			else if (((lambda+gamma+psi)/(lambda+gamma+mu+psi)) > specevent) {	
				sampled<-c(sampled,leaves[del])
				if (EI == TRUE && length(typessampled)<eliminate){
					typessampled<-c(typessampled,1)
				} else {typessampled<-c(typessampled,chosentype)}
				leaves <- leaves[- del]
				types <- types[- del]
				edge.length[edgespecevent] <- time-timecreation[- species]
				if (length(sampled) == n){
					stop = 1
				}
			} else {
				extinct <- c(extinct,leaves[del])
				leaves <- leaves[- del]
				types <- types[- del]
				edge.length[edgespecevent] <- time-timecreation[- species]
			}
		}
	}}
	
	while (length(leaves)>0) {
		del<-1
		extinct <- c(extinct,leaves[del])
		k = which( edge == leaves[del]  ) - length(edge.length)
		edge.length[k] <- time-timecreation[- leaves[del]]
		leaves <- leaves[- del]
	}

	for (j in 1:length(extinct)){
		del<-which (edge==extinct[j]) - length(edge.length)
		surpress<-edge[del,1]
		edge.length<-edge.length[- del]
		edge<-edge[- del,]
		del2<-which (edge[,1]==surpress)
		modify <-which (edge[,2]==surpress)
		edge[modify,2]<-edge[del2,2]
		edge.length[modify]<-edge.length[modify]+edge.length[del2]
		edge.length<-edge.length[- del2]
		edge<-edge[- del2,]
	}

	#nodes <- (length(sampled))*2
	leaf=1							#leaf labels 1...length(sampled)
	interior=length(sampled)+1		#int node labels length(sampled)+1 ...
	edgetemp<-edge
		
	typessampledNew<-vector()	
	temp<-unique(c(edge))				
	temp<-sort(temp,decreasing=TRUE)    #all nodes
	for (j in temp){
		#if (sum(match(sampled,j,0))==0 || j==-1) {
		if (sum(match(sampled,j,0))==0 || j==-1) {  #Tanja 19.4.12: why j==-1 special case?
			# replace all -j values in edge by interior
			posvalue <- interior
			interior <- interior +1
		} else {
			typessampledNew<-c(typessampledNew,typessampled[which(sampled==j)])
			posvalue <- leaf
			leaf <- leaf +1
		}
		replacel <- which(edge == j)
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
	phy$tip.label <- paste("t", sample(length(sampled)), sep = "")
	phy$edge.length <- edge.length
	phy$Nnode <- length(sampled)
	phy$states<-typessampledNew
	class(phy) <- "phylo"
	br<-sort(getx(phy,sersampling=1),decreasing=TRUE)
        phy$root.edge<-br[1]-br[2]
        phy<-collapse.singles(phy)
        phy2 <- phy

phy2
}
