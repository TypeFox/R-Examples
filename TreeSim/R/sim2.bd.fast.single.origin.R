sim2.bd.fast.single.origin <-
function(n,lambda,mu,rho,origin){
edge <- c(-1,-2)		#matrix of edges
leaves <- c(-2)			#list of extant leaves
timecreation <-c(0,0)		#time when species -2, -3, ... -n was created after origin
time <-0				#time after origin
maxspecies <- -2		#smallest species
edge.length <- c(0)		#edge length. if 0: leaf which didn't speciate /extinct yet
index = 2

specevents = vector()
specevents <-c(specevents,origin)

for (j in 1:(n-1)){
	r <- runif(1,0,1)
	tor=origin
	lamb1= rho * lambda
    mu1 = mu - lambda *(1-rho)
    if (lambda>mu) {
    temp <- 1/(lamb1-mu1)*log((lamb1-mu1* exp((-lamb1+mu1)*tor) -mu1*(1-exp((-lamb1+mu1)*tor)) *r )/(lamb1-mu1* exp((-lamb1+mu1)*tor) -lamb1*(1-exp((-lamb1+mu1)*tor)) *r )   )  
	} else {
	temp<- -((tor* r)/(-1 - lambda* rho* tor + lambda* rho* tor*r ))
		}
	specevents <- c(specevents,temp)
	}
specevents <- sort(specevents,decreasing=TRUE)
#specevents <- c(specevents,0)
while (index<=n) {
	
timestep <- specevents[(index-1)]-specevents[index]    #time since last event. specevent[j] time of (j-1)th spec event

time = time+timestep			#time after origin
species <- sample(leaves,1)  	#the leaf undergoing the next event
del <- which(leaves == species)
edgespecevent <- which(edge == species) - length(edge.length)
	edge.length[edgespecevent] <- time-timecreation[- species]
	edge <- rbind(edge,c(species,maxspecies-1))
	edge <- rbind(edge,c(species,maxspecies-2))
	edge.length <- c(edge.length,0,0)
	leaves <- c(leaves,maxspecies-1,maxspecies-2)
	maxspecies <- maxspecies-2
	leaves <- leaves[- del]
	timecreation <- c(timecreation,time,time)
	index <- index+1
}

time <- specevents[1]
# assign pendant edge length
for (j in (1:length(leaves))){
	k = which( edge == leaves[j]  ) - length(edge.length)
	edge.length[ k ] <- time - timecreation[- leaves[j]]
	}

# relabel
nodes <- (length(leaves))*2
leaf=1
interior=length(leaves)+1
edgetemp<-edge

if (nodes == 2) {
	# edge = c(2,1)	
	phy2 <-1	
	} else {
for (j in (1:nodes)){
	if (sum(match(leaves,- j,0))  == 0) {
		# replace all -j values in edge by interior
		posvalue <- interior
		interior <- interior +1
		} else {
		posvalue <- leaf
		leaf <- leaf +1
		}
	replacel <- which(edge == - j)
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
phy$tip.label <- paste("t", sample(length(leaves)), sep = "")
phy$edge.length <- edge.length
phy$Nnode <- length(leaves)
class(phy) <- "phylo"
storage.mode(phy$edge) <- "integer"

phy
}

