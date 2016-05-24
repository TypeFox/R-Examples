#' eco.space scape simulation with a macro-ecological focus
#' 
#' \code{eco.scape} is a modified version of the Helmus et al. method
#' implemented in \code{\link{scape}}. It produces phylogenetically
#' structured communities. It allows phylogenetic signals in niche
#' optima, but unlike \code{\link{scape}}, does not include the
#' ability to specify niche optima signal type (attraction/repulsion)
#' or phylogenetic signal in range size. Instead, the focus is on
#' having more control over the macroecological characteristics of the
#' resulting landscapes. In particular, eco.scape produces landscapes
#' with fixed mean range sizes, reasonable range size and abundance
#' distributions, and control over whether species present on a tree
#' must be present in the landscape.
#'
#' Simulates a landscape with species (i.e., tree tips) distributions
#' dependent on a supplied phylogenetic tree. The amount and type of
#' structure is determened by the signal parameter
#' \code{g.center}. Parameters are based on an Ornstein-Uhlenbeck
#' model of evolution with stabilizing selection. Values of g=1
#' indicate no stabilizing selection and correspond to the Brownian
#' motion model of evolution; 01 corresponds to disruptive selection
#' where phylogenetic signal for the supplied tree is amplified. See
#' corBlomberg. Communities are simulated along two gradients where
#' the positions along those gradients, \code{g.center}, can exhibit
#' phylogenetic signal.
#'
#' The function returns a landscape where the average range size is
#' equivalent to the wd.all parameter - in the scape function, this
#' parameter is not necessarily returned in the resulting
#' landscape. To do this, the probability of presence (th) that
#' returns the wd.all parameter is solved for. If there is no solution
#' that can produce the wd.all given, the error "Error in uniroot(f,
#' lower = 0, upper = max(X.), tol = 10^-200): f() values at end
#' points not of opposite sign" will occur. This seems to mostly arise
#' for extreme or unlikely parameter values (small species pools, low
#' carrying capacities). Try adjusting parameter values first.
#'
#' The \code{extinction} parameter specifies whether all of the
#' species on the tree should be present in the final landscape. Some
#' species will have probabilities of presence less than those
#' required for presence. If extinctions is \code{TRUE}, these species
#' will not be present. If \code{FALSE}, these species will be present
#' in 1 site, that in which they have the highest probability of
#' presence.
#'
#' @param tree \code{\link[ape:phylo]{phylo}} object; must have branch
#' lengths and be ultrametric
#' @param scape.size edge dimension of square landscape
#' @param g.center strength of phylogenetic signal in species range
#' centers. See \code{\link[ape:corBlomberg]{corBlomberg}},
#' 1=brownian,<1=rates of evol accelerate, >1=rates decelerate.
#' @param wd.all niche width, larger values simulate broader range
#' sizes
#' @param signal.center simulate with phylosignal in range centers 
#' @param center.scale adjust strength of phylogenetic attraction in
#' range centers independent of g.center
#' @param site.stoch.scale adjust strength of random variation in
#' species richness across sites
#' @param sd.center sd in \code{\link[stats:rnorm]{rnorm}}rnorm for
#' the range centers, increase to get more variation in center values
#' across species
#' @param sd.range sd in rnorm for the range sizes, increase to get
#' more variation in range sizes across gradients
#' @param K carrying capacity of a site in terms of maximum
#' individuals that can be present. Currently a constant value. Used
#' to scale the presence-absence matrix to include abundances.
#' @param extinction TRUE/FALSE can species on the tree go extinct on
#' the landscape? If the number of species present on the landscape
#' should equal the number of tips on the tree, choose FALSE. See
#' Details.
#' @param rho Grafen branch adjustment of phylogenetic tree see
#' \code{\link[ape:corGrafen]{corGrafen}}
#' @return
#' \item{cc}{\code{\link{comparative.comm}} object with
#' presence/absence results of simulations. The site names are the
#' row.columns of the cells in the original grid cells that made up
#' the data, and these co-ordinates are also given in the \code{env}
#' slot of the object along with the environmental gradient
#' information.}
#' \item{Y}{presence/absence matrix}
#' \item{Yab}{abundance matrix}
#' \item{index}{spatial coordinates for X and Y (stacked columns)}
#' \item{X.joint}{full probabilities of species at sites, used to
#' construct Y}
#' \item{X1}{probabilities of species along gradient 1}
#' \item{X2}{probabilities of species along gradient 2}
#' \item{gradient1, gradient2}{environmental gradient values}
#' \item{nichewd}{average niche width of the assemblage}
#' \item{K}{carrying capacity of each cell}
#' \item{environ}{matrix depicting environmental values over the 2D
#' landscape}
#' \item{sppXs}{full probabilities of each species as an array
#' arranged in a scape.size X scape.size matr ix}
#' \item{V.phylo}{initial phylogenetic covariance matrix from tree,
#' output of vcv.phylo(tree, corr=T)}
#' \item{V.phylo.rho}{phylogenetic covariance matrix from tree scaled
#' by Grafen if rho is provided, other wise just an output of
#' vcv.phylo(tree, corr=T)}
#' \item{V.center}{scaled (by g.center) phylo covariance matrix used
#' in the simulations}
#' \item{bspp1}{species optima for gradient 1}
#' \item{bspp2}{pecies optima for gradient 2}
#' @author Matt Helmus, Caroline Tucker, cosmetic edits by Will Pearse
#' @examples
#' # Simulations
#' tree <- rcoal(64)
#'
#' scape1 <- eco.scape(tree, scape.size=25, g.center=1,
#'     signal.center=FALSE, K=100, extinction=TRUE)
#' scape2 <- eco.scape(tree, scape.size=16, g.center=0.2,
#'     signal.center=TRUE, K=100, extinction=FALSE)
#' scape3 <- eco.scape(tree, scape.size=16, g.center=20,
#'     signal.center=TRUE, K=100, extinction=TRUE)
#' 
#' # Plotting distributions and landscape patterns
#' original_landscape <- scape1
#' abundmax <- original_landscape$K
#' PA_mat <- as.matrix(original_landscape$Y)
#' abund_mat <- original_landscape$Yab
#' site.size <- nrow(PA_mat)
#' species <- ncol(PA_mat)
#' mx <- original_landscape$gradient
#' env <- original_landscape$environ$env.gradient
#' par(mfrow=c(2,2), oma=c(0,0,2,0))
#' heatcol <- (colorRampPalette(c("yellow","red")))
#' 
#' image(matrix(env,sqrt(site.size),sqrt(site.size),byrow=TRUE),
#'     col=heatcol(max(env)),xaxt="n",yaxt="n",main="Env gradient")
#' 
#' image(matrix(rowSums(PA_mat),sqrt(site.size),sqrt(site.size),byrow=TRUE),
#'     col=heatcol(16),xaxt="n",yaxt="n",main="Species Richness")
#' 
#' hist(colSums(PA_mat),ylab="Number of species",xlab="Number of sites",
#'     main="Species Area Relationship",col="lightgrey")
#' 
#' hist(colSums(abund_mat),ylab="Number of species",xlab="Number of individuals",
#'     main="Species Abundance Relationship",col="lightgrey")
#' mtext("Env random, phy.signal=0.2, 32 species", outer=TRUE, side=3, cex=1.25)
#' @importFrom ape vcv.phylo is.ultrametric
#' @importFrom vegan decostand
#' @importFrom stats uniroot rnorm
#' @seealso \code{\link{scape}} \code{\link{sim.phy}} \code{\link{sim.meta}}
#' @export
eco.scape <- function(tree, scape.size=10, g.center=1, wd.all=0.2*(scape.size+1)^2, signal.center=TRUE, center.scale = 1, site.stoch.scale = 0, sd.center=1, sd.range=1, K=100, extinction = FALSE, rho=NULL){
    #Argument checking and assertion
    if(!inherits(tree, "phylo"))
        stop("'tree' must be of class 'phylo'")
    if(is.null(tree$edge.length))
        stop("'tree' must have branch lengths")
    if(!is.ultrametric(tree))
        stop("'tree' must be ultrametric")
    V <- vcv.phylo(tree, corr = TRUE)
    Vinit<-V
    
    #initialize
    nspp <- dim(V)[1]

    Xscale <- 2     #scale the strength of the probability matrix X
    Mscale <- site.stoch.scale    #scale stochasticity in niche distributions 
    Vscale1 <- Vscale2 <- center.scale  #scale the strength of the optimal values on axis one
    
    # Grafen's rho adjust strength of phylogenetic signal overall.
    if(!is.null(rho)){
        V <- 1-(1-Vinit)^rho
        V <- V/max(V)
    }

    #################################
    #SIMULATION######################
    #################################
    nsites<-scape.size  #number of sites for the square landscape
    #Establish environmental gradient
    mx2<- mx <- t(as.matrix((-(nsites)/2):(nsites/2)))  #env gradient
    m <- length(mx) #new number of sites (equal to nsites + 1)    
    ############
    ####Establish range centers/niche optima 
    
    if(signal.center){
        g<-abs(g.center)
        V.a<-vcv(corBlomberg(g, tree), corr=T)   #adjust phylogenetic signal for niche optima
        iD <- t(chol(V.a))
    } else {
        V.a <- V
        iD <- diag(nspp)
        dimnames(iD) <- dimnames(V.a)
    }
    iD.w <- diag(nspp)

    ##Species response to environmental gradient 1   
    e <- iD %*% rnorm(nspp, sd=sd.center)                                                         #assign optimal values as related to branch lengths, includes variation about mean value. Absolute values meaningless
    e <- Vscale1 * (e - mean(e))/apply(e, 2, sd)  #z-scores and scaling of the optimal values based on environmental signal in phylogeny
    bspp1 <- e

    spmx <- t((array(1, c(nspp, 1))) %*% mx) #env array as matrix
    mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e)) #max response to env for all sp
    wd <- iD.w %*% rnorm(nspp, sd=sd.range) #add small variation
    wd <- wd + (abs(min(wd))) #remove negative values
    
    X <- exp(-((spmx - mxsp)^2)/t(matrix(wd, nspp, m))) #Niche distributions  
    X1 <- Xscale * X  #Scales this initial species x site probability matrix 

    dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
    rat<-mean(dif/sort(wd)[-1])
    wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat    
    #Assign the zero with the mean ratio of nearest neighbor distances over the larger item
    wd <- wd.all*wd
    
    ##Species response to environmental gradient 2   
    e <- iD %*% rnorm(nspp, sd=sd.center)
    e <- Vscale2 * (e - mean(e))/apply(e, 2, sd)
    bspp2 <- e

    spmx2 <- t((array(1, c(nspp, 1))) %*% mx2)
    mxsp2 <- max(mx2)*((array(1, c(length(mx2), 1))) %*% t(e))
    wd <- iD.w %*% rnorm(nspp, sd=sd.range)
    wd <- wd + (abs(min(wd)))
    
    X <- exp(-((spmx2 - mxsp2)^2)/t(matrix(wd, nspp, m))) #Niche distributions     
    X2 <- Xscale * X

    dif <- sort(wd)[-1]-sort(wd)[-length(wd)] #Assign the zero to the nonzero minimum
    rat <- mean(dif/sort(wd)[-1])
    wd[wd==0] <- sort(wd)[2]-sort(wd)[2]*rat 
    wd <- wd.all*wd       #Assign the zero with the mean rato of nearist neighbor distances over the larger item


    ###CALCULATE JOINT PROBABILITY MATRIX, by combining X1 & X2
    X. <- NULL
    spp.Xs <- array(NA, dim=c(m, m, nspp))
    for(i in 1:nspp){
        sppX <- matrix((X1[,i]) %*% t(X2[,i]))
        spp.Xs[,,i] <- sppX
        X. <- cbind(X., matrix(sppX))
    }
    colnames(X.) <- colnames(X2)


    ## SOLVE for probability cutoff that gives mean range = wd.all  
    th <- uniroot(
        function(y) mean(apply(X., 2, function(x) sum(x>y))) - wd.all,
                     lower=0, upper=max(X.), tol=10^-200)$root
    
    
    ##SUMMARIZE PA matrix
    m. <- dim(X.)[1]
    Y <- matrix(0, ncol = nspp, nrow = m.)
    Y[th < X.] <- 1
    if(extinction==FALSE)
        for(i in  which(colSums(Y)==0))
            Y[sample((which(X.[,i]==max(X.[,i]))),1),i] <- 1
    
    colnames(Y) <- colnames(X.)
    index <- cbind(matrix(sapply(1:m, rep, times=m)), matrix(rep(1:m, times=m)))
    colnames(index) <- c("X1", "X2")
    
    ##WEIGHT by K for abundance matrix- new range of values for each site bounded by 0 (absent), each site sums to K
    Yab <- X.
    Yab[th > X.] <- 0
    
    Yab <- t(apply(Yab, 1, function(x) if(sum(x)>0) ceiling(x*K/sum(x)) else x)) #scale by carrying capacity K			
    Yab[Yab>0] <- floor(sapply(Yab[Yab>0], function(x){runif(1,min=(x-1),max=(x+5))}))
    for(i in which(colSums(Yab)==0))
        Yab[(which(Y[,i]==1)),i] <- 1
    colnames(Yab) <- tree$tip.label

    ##CREATE full environmental matrix (add mx1+mx2)
    env <- matrix(mx, nrow=length(mx), ncol=length(mx), byrow=TRUE)
    for(i in 1:nrow(env))
	env[i,] <- sapply(env[i,],function(x){x+mx2[i]})

    ########### OUTPUT
    rownames(Yab) <- paste(index[,1], index[,2], sep=".")
    env <- data.frame(env.gradient=env[index], row=index[,1], column=index[,2])
    rownames(env) <- rownames(Yab)
    cc <- comparative.comm(tree, Yab, env=env)
    return(list(cc=cc, Y=Y, Yab=Yab, gradient1 = mx, gradient2 = mx2, X.joint = X., X1 = X1, X2 = X2, nichewd = wd.all, K = K, environ = env, sppXs = spp.Xs, V.phylo = Vinit, V.phylo.rho = V, V.center = V.a, bsp1 = bspp1, bspp2 = bspp2))
}
