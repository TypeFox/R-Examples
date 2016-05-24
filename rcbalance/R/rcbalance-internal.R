.onAttach <- function(libname, pkgname){
	packageStartupMessage('optmatch (>= 0.9-1) needed to run the rcbalance command.  Please load optmatch and agree to its academic license before calling rcbalance.')	
}

dist2net.matrix <- function(dist.struct, k, exclude.treated = FALSE){
	ntreat <- nrow(dist.struct)
	ncontrol <- ncol(dist.struct)
		
	z <- c(rep(1, ntreat), rep(0, ncontrol))
	
    nobs <- length(z)			


	fb.structure <- data.frame('root.layer' = rep(1, nobs))	

    nums <- 1:nobs
    treated <- nums[z == 1]
    startn <- NULL
    endn <- NULL
    cost <- NULL
    ucap <- NULL

	if (inherits(dist.struct, 'InfinitySparseMatrix')){
		lookup.obj <- data.frame('treated' = attributes(dist.struct)$rows, 'control' = attributes(dist.struct)$cols, 'distance' = as.vector(dist.struct))	
	}

    #build treatment-control arcs		
    for (i in c(1:ntreat)) {
    	
	    	if (inherits(dist.struct, 'InfinitySparseMatrix')){
	    	  index.i <- which(lookup.obj$treated == i & 
	        is.finite(lookup.obj$distance))	
	   	  controls <- lookup.obj$control[index.i]
	   	  match.costs <- lookup.obj$distance[index.i]			
	    	} else {
	    	  controls <- which(is.finite(dist.struct[i,]))
	    	  match.costs <- dist.struct[i,controls] 	
	    	}
    
	    	if(!exclude.treated && length(controls) == 0){
	    		stop('Match is infeasible: some treated units have no potential control partners')
	    	}else if(length(controls) == 0){
	    		next	
	    	} 
	    	controls <- controls + ntreat #offset control node numbers so they don't overlap with treated node numbers
        startn <- c(startn, rep(treated[i], length(controls)))
        endn <- c(endn, controls)
		cost <- c(cost, match.costs)
        ucap <- c(ucap, rep(1, length(controls)))
    }
    b <- z*k #give each treated unit a supply of k
    tcarcs <- length(startn)
	
	#set up penalty vector
	
	p <- max(cost)
	theta <- 2
	S.i <- NULL
	layers <- list()	
	b <- c(b, -sum(k*z))
	layers[[1]] <- data.frame('input' = length(b), 'out' = NA, 'bypass' = NA)
	rownames(layers[[1]]) <- '1'	
	
	#now connect lowest layer to controls

	#find control nodes 
	ctrl.idx <- which(z == 0) #since T and C are first created, Cs should be these nodes

	#connect controls to lowest layer of fine balance treated
	low.layer <- ncol(fb.structure)
	parent.categories <- fb.structure[ctrl.idx,low.layer]
	#look up indices of parent nodes
	parent.nodes <- layers[[low.layer]]$input[match(parent.categories, rownames(layers[[low.layer]]))]
	
	#add edges between controls and fine balance nodes
	startn <- c(startn, ctrl.idx)
	endn <- c(endn, parent.nodes)
	ucap <- c(ucap, rep(1, length(ctrl.idx)))
	cost <- c(cost, rep(0, length(ctrl.idx)))
	
	#if we are excluding treated units, add bypass edges from treated units to lowest level of fine balance	
	if(exclude.treated){
		treat.idx <- which(z==1)
		#look up indices of lowest-level nodes for treated
		ll.categories <- fb.structure[treat.idx,low.layer]
		#look up indices of lowest-layer nodes
		ll.nodes <- layers[[low.layer]]$input[match(ll.categories, rownames(layers[[low.layer]]))]
		
		#add edges
		startn <- c(startn, treat.idx)
		endn <- c(endn, ll.nodes)
		ucap <- c(ucap, rep(1, length(treat.idx)))

		#set bypass cost to sum of ntreat largest t-c distances
		ntreat.largest <- sum(sort(cost, decreasing = TRUE)[1:ntreat])
		cost <- c(cost, rep(ntreat.largest, length(treat.idx)))
	}
	
	
	
	
	net.layers <- list(startn = startn, endn = endn, ucap = ucap, b = b, 
        cost = cost, tcarcs = tcarcs, layers = layers, z = z, fb.structure = fb.structure, penalties = S.i, theta = theta, p = p)
	net.layers
}

	


dist2net <-
function(dist.struct, k, exclude.treated = FALSE, ncontrol = NULL){
	ntreat <- length(dist.struct)
	if(is.null(ncontrol)) ncontrol <- max(laply(dist.struct, function(x) max(c(as.numeric(names(x)),0))))
		
	z <- c(rep(1, ntreat), rep(0, ncontrol))
	
    nobs <- length(z)			


	fb.structure <- data.frame('root.layer' = rep(1, nobs))	

    nums <- 1:nobs
    treated <- nums[z == 1]
    startn <- NULL
    endn <- NULL
    cost <- NULL
    ucap <- NULL
    

    #build treatment-control arcs		
    for (i in c(1:ntreat)) {
    	controls <- as.numeric(names(dist.struct[[i]]))
    	if(!exclude.treated && length(controls) == 0){
    		stop('Match is infeasible: some treated units have no potential control partners')
    	} else if(length(controls) == 0){
    		next	
    	} 
    	controls <- controls + ntreat #offset control node numbers so they don't overlap with treated node numbers
        startn <- c(startn, rep(treated[i], length(controls)))
        endn <- c(endn, controls)
        cost <- c(cost, dist.struct[[i]])
        ucap <- c(ucap, rep(1, length(controls)))
    }
    b <- z*k #give each treated unit a supply of k
    tcarcs <- length(startn)

	#set up penalty vector
	
	p <- max(cost)
	theta <- 2
	S.i <- NULL

	layers <- list()	
	b <- c(b, -sum(k*z))
	layers[[1]] <- data.frame('input' = length(b), 'out' = NA, 'bypass' = NA)
	rownames(layers[[1]]) <- '1'	
	
	#now connect lowest layer to controls

	#find control nodes 
	ctrl.idx <- which(z == 0) #since T and C are first created, Cs should be these nodes

	#connect controls to lowest layer of fine balance treated
	low.layer <- ncol(fb.structure)
	parent.categories <- fb.structure[ctrl.idx,low.layer]
	#look up indices of parent nodes
	parent.nodes <- layers[[low.layer]]$input[match(parent.categories, rownames(layers[[low.layer]]))]
	
	#add edges between controls and fine balance nodes
	startn <- c(startn, ctrl.idx)
	endn <- c(endn, parent.nodes)
	ucap <- c(ucap, rep(1, length(ctrl.idx)))
	cost <- c(cost, rep(0, length(ctrl.idx)))

	#if we are excluding treated units, add bypass edges from treated units to lowest level of fine balance	
	if(exclude.treated){
		treat.idx <- which(z==1)
		#look up indices of lowest-level nodes for treated
		ll.categories <- fb.structure[treat.idx,low.layer]
		#look up indices of lowest-layer nodes
		ll.nodes <- layers[[low.layer]]$input[match(ll.categories, rownames(layers[[low.layer]]))]
		
		#add edges
		startn <- c(startn, treat.idx)
		endn <- c(endn, ll.nodes)
		ucap <- c(ucap, rep(1, length(treat.idx)))
		#set bypass cost to sum of ntreat largest t-c distances
		ntreat.largest <- sum(sort(cost, decreasing = TRUE)[1:ntreat])
		cost <- c(cost, rep(ntreat.largest, length(treat.idx)))
	}
	
	net.layers <- list(startn = startn, endn = endn, ucap = ucap, b = b, 
        cost = cost, tcarcs = tcarcs, layers = layers, z = z, fb.structure = fb.structure, penalties = S.i, theta = theta, p = p)
	net.layers
}

add.layer <-
function(net.layers, new.layer){
	#net.layers is a layered network object
	startn <- net.layers$startn
	endn <- net.layers$endn
	ucap <- net.layers$ucap
	b <- net.layers$b
	cost <- net.layers$cost
	tcarcs <- net.layers$tcarcs
	layers <- net.layers$layers
	z <- net.layers$z
	fb.structure <- net.layers$fb.structure
	S.i <- net.layers$penalties
	theta <- net.layers$theta
	my.p <- net.layers$p
	#figure out from previous network structure what control ratio is
	k <- b[min(which(z == 1))]

	#find where new layer fits in fine balance structure
	n.levels <- ncol(fb.structure)
	parent.layer <- NA
	for(i in c(n.levels:1)){
		nest.tab <- table(fb.structure[,i], new.layer)
		ncoarse.by.fine <- apply(nest.tab, 2, function(x)sum(x > 0))
		#check if new.layer nests inside this layer
		if(all(ncoarse.by.fine <= 1)){
			if(i == n.levels){
				parent.layer <- i
				fb.structure <- cbind(fb.structure, new.layer)
				colnames(fb.structure)[i+1] <- paste('f',i+1,sep ='.') 
				layers[[i+1]] <- list()
				n.levels <- n.levels + 1
				break
			}
			#for i < n.levels, check nesting in lower layer
			nest.tab2 <- table(new.layer, fb.structure[,i+1])
			ncoarse.by.fine2 <- apply(nest.tab2, 2, function(x)sum(x > 0))			
			if(all(ncoarse.by.fine2 <= 1)){
				parent.layer <- i
				fb.structure <- cbind(fb.structure[,c(1:i)], new.layer,fb.structure[,c((i+1):n.levels)])
				temp.layers <- list(length = n.levels + 1)
				temp.layers[c(1:i,(i+2):(n.levels + 1))] <- layers[c(1:i,(i+1):n.levels)]
				temp.layers[[i+1]] <- list()
				layers <- temp.layers
				n.levels <- n.levels + 1
				colnames(fb.structure)<- paste('f',c(1:n.levels), sep = '.')
				break
			}
		}
	}
	stopifnot(!is.na(parent.layer)) #stop if new layer doesn't nest correctly
	
	#change index of interest to newly added layer
	i <- parent.layer + 1
	
	#update penalty vector S.i
	#also update bypass penalties if we have them
	if(length(S.i) == 0){
		S.i <- my.p*theta
	}else{
		S.i <- c(theta*S.i[1], S.i)
		for(j in c(2:length(S.i))){
			if(j >= i) break #only need to update penalties in layers above new one
			cost[cost == S.i[j]] <- S.i[j-1] #step up penalty to higher level
		}		
	}
	#check if there are bypass edges from treated layer; if so update their penalties too
	if(any(startn[-c(1:tcarcs)] %in% which(z ==1)) && length(S.i) > 0){
		bypass.edges <- startn %in% which(z == 1)
		bypass.edges[1:tcarcs] <- FALSE
		new.byp.pen <- theta*max(S.i)
		cost[which(bypass.edges)] <- new.byp.pen
	}
	
	#find parent nodes for nodes in current layer
	ztab <- table(fb.structure[,i], z) 
	zerobins <- which(apply(ztab, 1,function(x)all(x == 0))) 
	if(length(zerobins) > 0){	
		ztab <- ztab[-zerobins,]
	}
	nnodes <- nrow(ztab)
	node.nums <- length(b) + c(1:nnodes)

	nest.tab <- table(fb.structure[,i-1], fb.structure[,i])
	parent.categories <- apply(nest.tab, 2, function(x) rownames(nest.tab)[which (x > 0)] )
	if(length(zerobins) > 0){			
		parent.categories <- parent.categories[-zerobins]
	}

	#put fine balance structures in the network
	#start with output nodes
	out.nums <- length(b) + c(1:nnodes)
	b <- c(b, rep(0, nnodes))
	node.lookup <- match(parent.categories, rownames(layers[[i-1]]))
	parent.nodes <- layers[[i-1]]$input[node.lookup]
	
	#detach the parent node layer from its former child 
	drop.edges <- which(endn %in% parent.nodes)
	#check if we are using bypass edges from treated layer so we can add them back in later
	exclude.treated <- any(which(z==1) %in% startn[drop.edges])
	
	if(length(drop.edges) > 0){
		startn <- startn[-drop.edges]
		endn <- endn[-drop.edges]
		ucap <- ucap[-drop.edges]
		cost <- cost[-drop.edges]
	}
	
	startn <- c(startn, out.nums)
	endn <- c(endn, parent.nodes)
	ucap <- c(ucap, rep(sum(k*z), nnodes))
	cost <- c(cost, rep(0,nnodes))
			
	#now do input nodes
	in.nums <- length(b) + c(1:nnodes)
	b <- c(b, rep(0, nnodes))

	startn <- c(startn, in.nums)
	endn <- c(endn, out.nums)		
	ucap <- c(ucap, k*ztab[,2]) #counts in treated populations times k are capacities here
	cost <- c(cost, rep(0,nnodes))
			
	#finally do bypass nodes
	bypass.nums <- length(b) + c(1:nnodes)
	b <- c(b, rep(0, nnodes))
			
	startn <- c(startn, in.nums)
	endn <- c(endn, bypass.nums)
	startn <- c(startn, bypass.nums)
	endn <- c(endn, out.nums)			
	ucap <- c(ucap, rep(sum(k*z), 2*nnodes)) 
	cost <- c(cost, rep(S.i[i-1],2*nnodes)) #give these edges high costs
			
	layers[[i]] <- data.frame('input' = in.nums, 'out' = out.nums, 'bypass' = bypass.nums)
	rownames(layers[[i]]) <- rownames(ztab)

	#need to attach new layer to controls or child.layer
	if(i == n.levels){
		#connect new lowest layer to controls
		#find control nodes 
		ctrl.idx <- which(z == 0) #since T and C are first created, Cs should be these nodes

		#connect controls to lowest layer of fine balance treated
		low.layer <- ncol(fb.structure)
		parent.categories <- fb.structure[ctrl.idx,low.layer]
		#look up indices of parent nodes
		parent.nodes <- layers[[low.layer]]$input[match(parent.categories, rownames(layers[[low.layer]]))]
		
		#add edges between controls and fine balance nodes
		startn <- c(startn, ctrl.idx)
		endn <- c(endn, parent.nodes)
		ucap <- c(ucap, rep(1, length(ctrl.idx)))
		cost <- c(cost, rep(0, length(ctrl.idx)))
		
		
		#Add bypass edges from treated to lowest fine balance layer
		if(exclude.treated){
				treat.idx <- which(z==1)
				#look up indices of lowest-level nodes for treated
				ll.categories <- fb.structure[treat.idx,low.layer]
				#look up indices of lowest-layer nodes
				ll.nodes <- layers[[low.layer]]$input[match(ll.categories, rownames(layers[[low.layer]]))]
		
				#add edges
				startn <- c(startn, treat.idx)
				endn <- c(endn, ll.nodes)
				ucap <- c(ucap, rep(1, length(treat.idx)))
				#EDIT
				#set bypass cost to max t-c distance plus one
				if(length(S.i) > 0){
					new.byp.pen <- theta*max(S.i)
				}else{
					new.byp.pen <- sum(sort(cost, decreasing = TRUE)[1:sum(z)])
				}
				cost <- c(cost, rep(new.byp.pen, length(treat.idx)))
		}
						
	}else{
		i <- i + 1 #now i indexes child layer
		
		ztab <- table(fb.structure[,i], z) 
		zerobins <- which(apply(ztab, 1,function(x)all(x == 0))) 
		if(length(zerobins) > 0){	
			ztab <- ztab[-zerobins,]
		}
		nnodes <- nrow(ztab)
		node.nums <- layers[[i]]$out
		
		nest.tab <- table(fb.structure[,i-1], fb.structure[,i])
		parent.categories <- apply(nest.tab, 2, function(x) rownames(nest.tab)[which (x > 0)] )
		pc.found <- sapply(parent.categories, length)
		zerobins <- which(pc.found == 0)
		if(length(zerobins) > 0){			
			parent.categories <- parent.categories[-zerobins]
		}
		node.lookup <- match(parent.categories, rownames(layers[[i-1]]))
		parent.nodes <- layers[[i-1]]$input[node.lookup]
	
		#now add in new edges connecting new layer to parent
		startn <- c(startn, node.nums)
		endn <- c(endn, parent.nodes)
		ucap <- c(ucap,rep(sum(k*z),nnodes))
		cost <- c(cost, rep(0,nnodes))
	}
	net.layers <- list(startn = startn, endn = endn, ucap = ucap, b = b, 
        cost = cost, tcarcs = tcarcs, layers = layers, z = z, fb.structure = fb.structure, penalties = S.i, theta = theta, p = my.p)
	net.layers
}


callrelax <- function (net) {
	if (!requireNamespace("optmatch", quietly = TRUE)) {
	  stop('Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.')
	}
    	startn <- net$startn
    	endn <- net$endn
    	ucap <- net$ucap
    	b <- net$b
    	cost <- net$cost
    	stopifnot(length(startn) == length(endn))
    	stopifnot(length(startn) == length(ucap))
    	stopifnot(length(startn) == length(cost))
    	stopifnot(min(c(startn, endn)) >= 1)
    	stopifnot(max(c(startn, endn)) <= length(b))
    	stopifnot(all(startn != endn))
    	nnodes <- length(b)
    my.expr <- parse(text = '.Fortran("relaxalg", nnodes, as.integer(length(startn)), 
    	    as.integer(startn), as.integer(endn), as.integer(cost), 
    	    as.integer(ucap), as.integer(b), x1 = integer(length(startn)), 
    	    crash1 = as.integer(0), large1 = as.integer(.Machine$integer.max/4), 
    	    feasible1 = integer(1), NAOK = FALSE, DUP = TRUE, PACKAGE = "optmatch")')
	fop <- eval(my.expr)	
   	x <- fop$x1
    	feasible <- fop$feasible1
    	crash <- fop$crash1
    	list(crash = crash, feasible = feasible, x = x)
}

penalty.update <-
function(net.layers, newtheta, newp = NA){
	oldpen <- net.layers$penalties
	if(length(oldpen)== 0) return(net.layers)
	if(is.na(newp)) newp <- net.layers$p #rev(oldpen)[1]/net.layers$theta #if p is not supplied set it to old value of p
	newpen <- newtheta^c(length(oldpen):1)*newp
	oldcost <- net.layers$cost
	newcost <- net.layers$cost
	for(i in c(1:length(oldpen))){
		newcost[which(oldcost == round(oldpen[i]))] <- newpen[i]
	}
	#check if there are bypass edges from treated layer; if so update their penalties too
	if(any(net.layers$startn[-c(1:net.layers$tcarcs)] %in% which(net.layers$z ==1))){
		bypass.edges <- net.layers$startn %in% which(net.layers$z == 1)
		bypass.edges[1:net.layers$tcarcs] <- FALSE
		new.byp.pen <- newtheta*max(newpen)
		newcost[which(bypass.edges)] <- new.byp.pen
	}	
	net.layers$cost <- newcost
	net.layers$penalties <- newpen
	net.layers$theta <- newtheta
	return(net.layers)
}

penalize.near.exact <- function(net.layers, near.exact){
	oldcost <- net.layers$cost
	startn <- net.layers$startn
	endn <- net.layers$endn
	tcarcs <- net.layers$tcarcs
	theta <- net.layers$theta
	z <- net.layers$z
	
	newcost <- oldcost
	if(any(startn[-c(1:tcarcs)] %in% which(z ==1))){	
		bypass.edges <- startn %in% which(z == 1)
		bypass.edges[1:tcarcs] <- FALSE
		old.byp.pen <- oldcost[which(bypass.edges)[1]]
		newcost[which(bypass.edges)] <- old.byp.pen*theta
		near.exact.pen <- old.byp.pen
	}else{
		near.exact.pen <- theta*max(net.layers$penalties,0)
		#correct penalty if net.layers$penalties was empty
		if(near.exact.pen == 0) near.exact.pen <- theta*net.layers$p
	}
	newcost[which(near.exact[startn] != near.exact[endn])] <- newcost[which(near.exact[startn] != near.exact[endn])] + near.exact.pen
	net.layers$cost <- newcost
	return(net.layers)
}

.Random.seed <-
c(403L, 11L, -1638544231L, 1706315048L, -1076233258L, -1306066863L, 
1371861763L, 465437330L, 1377005604L, 1240702695L, 822719293L, 
455789492L, 619198490L, 613632869L, 1477367295L, 112285830L, 
-961543152L, 923367091L, -1610177199L, 2057944320L, -1349498914L, 
481084457L, -1843797029L, -1198393110L, -932587892L, -1882977841L, 
-1661665915L, -458496452L, 1904077586L, -211482611L, -1828160697L, 
-561303250L, 1246165128L, -584606325L, 1743676777L, 1888880440L, 
-995915290L, 373370241L, -1298235117L, 162373314L, -89876108L, 
155304311L, -1160949139L, 1722178820L, -1533495254L, 590943637L, 
1535982127L, -1695447242L, 1867315040L, -179593821L, -151185727L, 
-301502800L, -236795442L, 3644729L, -1345236661L, -185224454L, 
-78994436L, -2079144897L, -1349305067L, -1289214420L, 1288407810L, 
-815884323L, 1707492055L, 372006334L, -1977526280L, -350989349L, 
-90575303L, 1399453576L, -1684666570L, 1327347057L, 1617827L, 
-1118726414L, 1774804548L, 594411783L, 85604573L, -1668633132L, 
-813096710L, -276481979L, 118538399L, 1772907046L, 913104560L, 
-92961965L, 669609585L, 1911837600L, -1856277698L, 1683937865L, 
-905646341L, 1926524618L, -352391124L, -557742033L, -1875560219L, 
-1416834660L, 1486519282L, 2007598701L, 2013860775L, -660271666L, 
-749450456L, -757415637L, 1367014217L, -177143720L, 552252678L, 
165781729L, 952211571L, 1110865122L, 1298515028L, -1457493673L, 
1353615693L, -645419804L, 1282969802L, -772473547L, 264868751L, 
551291030L, -724264896L, -643661565L, -1699058399L, -751256048L, 
-1567322642L, 1484685721L, 1519668779L, 1132222234L, 844279260L, 
1910412255L, -427093067L, 259157324L, 1412369314L, -1731468547L, 
1166408695L, -2087754466L, 837863384L, -156542661L, 819293913L, 
675267176L, 1096098710L, -683276527L, -238034749L, -1421295406L, 
-1552013980L, -1257547609L, 723381757L, 391581684L, 1193310426L, 
-1676543451L, -1153018049L, 1902531782L, -1887844016L, -348341901L, 
-1797569647L, 1229936064L, 1058918302L, -100437655L, -128319077L, 
1198306474L, 1028780748L, 829881231L, -703807547L, -1858837636L, 
-175121454L, 624188237L, 1879091719L, 1022720494L, 246289608L, 
-1394623541L, -1150734807L, 1793975416L, 554148262L, -1271774271L, 
537784403L, -1188940030L, -103185868L, 610105271L, -7837779L, 
1745261764L, 139794666L, 339394389L, -1311911697L, -903670026L, 
929219104L, -479934493L, -1203458943L, 1777593072L, -1478827634L, 
-724168967L, 561822859L, -2060540742L, 1437384636L, -1933722369L, 
1044779093L, -2010070164L, 1534671810L, 316042269L, 958188695L, 
192244222L, 1887037880L, 1564108187L, 133909497L, 67523400L, 
-822919562L, 2075916209L, 2128284131L, -651977038L, 900387588L, 
-67104697L, 1986901789L, -197309036L, 923251002L, 1884294533L, 
-61766049L, -2071039770L, 273691504L, 319886483L, 1102844209L, 
-1232406560L, -106431362L, -102174967L, 1812246843L, -1288039158L, 
-660255508L, 1087813999L, 422985381L, 1916768348L, -237025230L, 
-568361171L, -1408168473L, -415330034L, 1087960808L, 369105131L, 
-1070147191L, 876816664L, 1991294022L, -1241329726L, -2138524888L, 
-753010260L, 1690279228L, 1740495730L, -332223184L, 1190618916L, 
-165001928L, 1378584282L, -685072496L, 1414765308L, -839459948L, 
708670530L, -1409944320L, -1652914116L, -1199719088L, -1356077918L, 
-1775555768L, -440508660L, -749011236L, -216354990L, 109672064L, 
-1852884908L, 808092648L, -828331270L, 819582288L, 1979907884L, 
1667621748L, -1862851694L, 1146253920L, 1716066188L, 1120319328L, 
1744290786L, -518726744L, -838087476L, -702211588L, -1847953486L, 
-440039952L, -1447647804L, -1089013416L, -391639686L, -4801104L, 
-77470788L, -438611884L, 691695554L, 31514656L, 960262332L, 1520294480L, 
-860313630L, 1409252968L, -1657090484L, -1310805764L, -512821262L, 
1549133056L, -999742604L, 6725640L, -415370822L, 187522512L, 
2058646508L, 190466964L, 885742546L, -1511807136L, -528973492L, 
579924256L, -710733950L, -1153890584L, -615196436L, -145729988L, 
2102162290L, 749033072L, 502359076L, -1673842888L, 1576016026L, 
1927573200L, -884912964L, 1489619604L, 1090708098L, -621211456L, 
1291158844L, 1904040080L, 1024185890L, -1413467640L, -1365951220L, 
-1660337508L, -482150894L, -1634592384L, -1657595372L, 1221090856L, 
-1392123590L, 1737075664L, 1098010220L, -1847016652L, 490660498L, 
1475175776L, -435626676L, 251261920L, 184307618L, 2055241640L, 
72705932L, -314602820L, 1039483762L, -296900368L, 1537725508L, 
342848856L, 1297281082L, 239096688L, 2116830268L, 2034846484L, 
1234802114L, -580432672L, -1461346436L, 2004281040L, -1821754718L, 
234533928L, -1591189236L, -2109400132L, 246869298L, -546764864L, 
251692596L, -1707486264L, -861049670L, -500408432L, 466403820L, 
356492372L, 1477319314L, 1046301472L, 234154444L, 1549796320L, 
1612467906L, -1224972760L, 589184684L, 764604092L, -255913486L, 
1721880880L, -1402627932L, 1402624696L, 1173718874L, -2039390832L, 
384494972L, 459278612L, -1472506558L, -1630096640L, -1973795908L, 
1145320272L, -788411102L, -2018395832L, 810691596L, 1444045916L, 
-1061468974L, -1672072192L, 38647636L, -692438552L, -1177091078L, 
2147427024L, -1190040276L, -581073676L, 1830044562L, -1496492960L, 
1838319884L, 1552302816L, -1146679966L, -204650840L, 1511464524L, 
1260938108L, 513339186L, 1867090544L, -563560124L, -241535016L, 
-2123262854L, -496039376L, 704689852L, -1946224812L, -1570191806L, 
657282976L, 1111368892L, 1675289040L, -880152222L, -136253464L, 
-387329076L, 2047878908L, 1718698994L, 125938816L, 926189812L, 
810541192L, 718963002L, 386558800L, 1094565996L, -1905410924L, 
902022482L, -1180497184L, 1971850444L, -1258457952L, 77264258L, 
-144562968L, -934312724L, -1752795076L, 838233202L, 456694896L, 
356807204L, 997000760L, 286302746L, 595053776L, 1311092284L, 
-206268268L, -830094206L, 599589824L, 848898364L, -1582150000L, 
1606997026L, 957644808L, -4905844L, 1273912348L, 1479565458L, 
615046272L, -1684583404L, -2083432664L, 1568849978L, -1622419248L, 
-388116372L, 104400564L, -1209578478L, 1505190880L, 1976103884L, 
-600679968L, -696702942L, -205607640L, 1145811212L, 1326061919L, 
-884117271L, 907724798L, 360729580L, -94175427L, -879795449L, 
-1689109776L, -143547762L, -500040517L, -197388131L, -911612854L, 
-1801257368L, -948309423L, -580308669L, 833726532L, -186160398L, 
2058400455L, -1379233647L, 445685270L, -1592828812L, 1617942277L, 
1425569103L, 1292306248L, 1328090854L, 1722265619L, 264193909L, 
-474372270L, 1424433056L, 510700169L, 1496764091L, -1688404724L, 
1102119066L, -1378912753L, -1275141543L, 1399971566L, -1825178244L, 
1181990381L, -1060267369L, 1077835104L, -1019612866L, 752638475L, 
-1547608435L, 1738621882L, -1087284744L, -1985895007L, 1391573843L, 
-506801548L, 199947266L, 766269079L, 207601761L, 2136333478L, 
-961434460L, 1554279765L, 1783186943L, 930696344L, 404060598L, 
782713411L, 957785989L, -475134942L, 1097631504L, -1358196039L, 
-946894997L, 1462830684L, 860769162L, 1428923263L, -232796471L, 
369395550L, -1986159604L, -1401185315L, 1175281383L, -43229424L, 
-587420114L, -962669349L, 803223997L, 505313066L, 1203703752L, 
-1485314831L, -362309021L, -1822458268L, 674137490L, -2090831257L, 
-851580367L, 1122469750L, -1830605676L, 1846723237L, -1869492305L, 
1159351400L, 1337401414L, 217875955L, 302764373L, -1404886222L, 
-1308641536L, -1814764567L, 48201307L, 1622974636L, 52815546L, 
-933735313L, 1669796345L, -2071747442L, -1680011940L, -1731142195L, 
-2019286729L, -1218049088L, 675161246L, -732837525L, 1720168685L, 
-1793472166L, 766448152L, -1734881407L, -277406413L, 216077140L, 
1235879266L, 944166391L, 1384103489L, 827947718L, -1384739836L, 
-721253451L, -2128894817L, -1048604424L, -1637111594L, 1959706979L, 
-227914459L, -352232766L, 1999591344L, -1763750951L, 1977607499L, 
-602913860L, 561761258L, -1721414625L, -974894423L, -1425424834L, 
1445589804L, 1532059005L, 481915847L, 496313008L, 576856270L, 
-1068156421L, -1883252131L, -1862830710L, 1145359784L, 1645802257L, 
-18611709L, 1606751492L, 922382642L, 102064007L, -613159471L, 
-1814552746L, 1755643188L, 249091141L, -840745841L, -1129113080L, 
-1952254426L, 2145808083L, 1886249013L, 1330212754L, -1303195808L, 
-1246299447L, -1478193669L, -592373300L, 667914330L, 774147407L, 
-30915559L, -2027386322L, 1525409852L, -176644691L, 935910999L, 
1134374816L, -881685250L, 604057547L, 877485219L)
