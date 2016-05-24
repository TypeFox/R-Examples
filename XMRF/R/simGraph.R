simGraph <-
function(p=50, type="scale-free"){
	if(type=="scale-free"){
		g <- igraph::barabasi.game(n=p, power=0.01, zero.appeal=p, directed=F)
		theta <- as.matrix(igraph::get.adjacency(g, type="both"))
	}
	#
	if(type=="hub"){
		if(p > 40) {g = ceiling(p/20)}
		if(p <= 40) {g = 2} 
		
		g.list = c(rep(floor(p/g), (g - p%%g)), rep(floor(p/g)+1, p%%g))
		g.ind = rep(c(1:g), g.list)
		u = 0.1
		v = 0.3
		theta = matrix(0, p, p)
		for (i in 1:g) {
            tmp = which(g.ind == i)
            theta[tmp[1], tmp] = 1
            theta[tmp, tmp[1]] = 1
        }
	}
	#
	if(type=="lattice"){
		
		latticeStructure <- function(p){
			nnodes <- sqrt(p)
			if((nnodes %% 1) == 0){
				rnodes = nnodes
				cnodes = nnodes
			}
			if((nnodes %% 1) != 0){
				cnodes = 3
				rnodes = p / cnodes
				while((rnodes %% 1) != 0){
					cnodes = cnodes + 1
					rnodes = p / cnodes
				}
			}
			return(c(rnodes, cnodes))
		}
		
		temp = latticeStructure(p)
		rnodes = temp[1]
		cnodes = temp[2]
		
		theta=diag(1,nrow=rnodes*cnodes,ncol=cnodes*rnodes)
		alterSign = FALSE
		for(i in 1:(rnodes*cnodes)){
			#print(i);
			if( i%% rnodes !=0 ){
				theta[i,i+1]=1
			}
			if((i-1)%%rnodes!=0){
				theta[i,i-1]=1
			}
			if((i-rnodes)>0){
				theta[i,i-rnodes]=1
				if(alterSign) 	theta[i,i-rnodes] = -1

			}			
			if((i+rnodes)<=rnodes*cnodes){
				theta[i,i+rnodes]=1
				if(alterSign) 	theta[i,i+rnodes] = -1

			}			
		}
	}
	diag(theta) <- 0
	return(theta)
}
