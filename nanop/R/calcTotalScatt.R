###################################################################################################################################
#  X-RAY SCATTERING FACTORS
#
#     x      a1        b1         a2        b2        a3        b3          a4         b4         a5         b5        c
#    'C' : [2.657506, 14.780758, 1.078079, 0.776775, 1.490909, 42.086843,  -4.241070, -0.000294, 0.713791, 0.239535, 4.297983],
#    'S' : [6.372157, 1.514347, 5.154568, 22.092528, 1.473732, 0.061373,   1.635073,  55.445176, 1.209372, 0.646925, 0.154722],
#    'Pd': [6.121511, 0.062549,  4.784063, 0.784031, 16.631683, 8.751391,  4.318258, 34.489983, 13.246773, 0.784031, 0.883099],
#    'Ag': [6.073874, 0.055333, 17.155437, 7.896512, 4.173344, 28.443739,  0.852238, 110.376108, 17.988685, 0.716809, 0.756603],
#    'Au': [16.777389, 0.122737, 19.317156, 8.621570, 32.979682, 1.256902, 5.595453, 38.008821, 10.576854, 0.000601, -6.279078]
#  Table (1) of 
#      D. WAASMAIER AND A. KIRFEL, Acta Cryst. (1995). A51, 416-431 


calcTotalScatt <- function(nanop, dQ=.01, minQ=0.771, maxQ=35, type="neutron", 
                           scatterFactor=NA, scatterLength=NA, sigma=NA, 
						   n=0, delta=0, kind="fastHist", 
						   dr = 0.001,  del = 0.01, eps=1e-3) {
						   
#######################################################################################
# Q = 4*pi*sin(theta)/lambda
# kind="fast_av" ## for polydisperse particles
# kind="fast" 
# kind="exact"
# kind="fastHist"
# rmax, dr = histogram params
# rmax should be greater than particle diameter!
# del, eps = fast Cervillino approach params
# 
#
#
  diam_coef <- 1
  if(attributes(nanop)$dimer)
    diam_coef <- diam_coef*2

  if( kind!="fast" & kind!="fast_av" & kind!="exact" & kind!="fastHist"){
    warning("switch error; calculating exact scattering function \n", immediate. = TRUE)
    kind <- "exact"  
  }
  
  if(kind != "fastHist")
    dr <- 0
  
  useN <- !(n==0 && delta == 0)
  if (delta != 0){
    useN <- TRUE
	if(is.na(n) || n<=0 )
	  n <- 2 
  }

  nAtomTypes <- attributes(nanop)$nAtomTypes
  if (is.na(scatterLength[1]))
    scatterLength <- attributes(nanop)$scatterLength
	
  if (length(scatterLength) < nAtomTypes){
    warning("vector <scatterLength> length is less than number of atom types nAtomTypes\n", immediate. = TRUE)
	scatterLength <- rep(scatterLength[1], nAtomTypes)
  }
  if (is.na(sigma[1]))
    sigma <- attributes(nanop)$sigma

  if (is.na(sigma[1]))
    sigma <- 0
	
  if (length(sigma) < nAtomTypes){
    warning("vector <sigma> length is less than number of atom types nAtomTypes\n", immediate. = TRUE)
	sigma <- rep(sigma[1], nAtomTypes)
  }
  
  if (is.na(scatterFactor[[1]][1]))
    scatterFactor <- attributes(nanop)$scatterFactor
	
  if (length(scatterFactor[[1]]) < nAtomTypes){
    warning("list <scatterFactor> length is less than number of atom types nAtomTypes\n", immediate. = TRUE)
	for(i in 1:11)
	  scatterFactor[[i]] <- rep(scatterFactor[[i]][1], nAtomTypes)
  }
  a1=scatterFactor$a1
  b1=scatterFactor$b1
  a2=scatterFactor$a2
  b2=scatterFactor$b2
  a3=scatterFactor$a3
  b3=scatterFactor$b3
  a4=scatterFactor$a4
  b4=scatterFactor$b4
  a5=scatterFactor$a5
  b5=scatterFactor$b5
  c = scatterFactor$c
  
  if(type == "neutron"){
    type <- 0
  }else{ 
    type <- 1 
  }

   
  atomType <- attributes(nanop)$atomType
  layers_num <- length(attributes(nanop)$layer_end)
  tt <- 1:max(atomType)
  if(min(atomType) < 0)
    tt <- c(tt, -1:min(atomType))
#  tt <- sort(tt)
  layer_start <- attributes(nanop)$layer_start	
  layer_end <- attributes(nanop)$layer_end  
  
 
  iATL <- list()  #  positions of atoms in whole array 'nanop' for each atom type (first index) for each layer (second index) 
  for(i in 1:nAtomTypes){
    iATL[[i]] <- list()
	for(j in 1:layers_num){
      rs <- which(atomType == tt[i]) # only specific atom type
	  if(tt[i]>0){
	    rs <- rs[which(rs >= layer_start[j])]
	    rs <- rs[which(rs <= layer_end[j])]
	  }
	  else{
 	    layerS_start <- attributes(nanop)$layerS_start	
        layerS_end <- attributes(nanop)$layerS_end
	    rs <- rs[which(rs >= layerS_start[j] + attributes(nanop)$rowcore)]
	    rs <- rs[which(rs <= layerS_end[j] + attributes(nanop)$rowcore)]	  
	  }
	  iATL[[i]][[j]] <- rs	
    }
  }
  
  nAt <- 0 # number of atoms of each type
  for(i in 1:nAtomTypes){
    nAt[i] <- length( which(atomType == tt[i]) )  
  }
  iAtomType <- list() # positions of atoms in whole array 'nanop' for each atom type
  for(i in 1:nAtomTypes){
    iAtomType[[i]] <-  which(atomType == tt[i])
  }
  
 
  fAv <- 0  # average sc. length for each layer
  NN <- 0
  sx <- 0   # scale  factor 1/N<f>^2 for a given layer
  Q <- seq(minQ, maxQ, by=dQ)
  
  for(k in 1:length(Q)){ 
    for(i in 1:layers_num){
      fAv[i] <- 0
	  NN[i] <- 0	  
	  if(type==0){
	    ff <- scatterLength
	  }else{
	    q4 = (Q[k]/(4*pi))^2
        ff = a1*exp(-b1*q4) + a2*exp(-b2*q4) + a3*exp(-b3*q4) + 
		        a4*exp(-b4*q4) + a5*exp(-b5*q4)+c		
	  }
	  for(j in 1:nAtomTypes){
	    fAv[i] <- fAv[i] + length(iATL[[j]][[i]])*ff[j]
	    NN[i] <- NN[i] + length(iATL[[j]][[i]])
  	  }
      fAv[i] <- fAv[i]/NN[i]
	  sx[(i-1)*length(Q)+k] <- 1 / (NN[i]*fAv[i]^2)  # 1/N<f>^2 for given layer
    }
  }
#########################################	
#     uniform model  
  if(nAtomTypes==1){
    x2 <- nrow(nanop)+1
	x1 <- 1
    if(attributes(nanop)$shape=="ellipse")
	  diam <- 2*max(attributes(nanop)$r)
	else if(attributes(nanop)$shape=="box")
      diam <- 2*sqrt(sum(attributes(nanop)$r^2))
    else
	  diam <- 2*attributes(nanop)$r
 	if( (x2-x1) <= 1) 
	  res <- list(Q=Q, gQ = rep(0,length(Q)))
	else{
	  if(kind=="fast_av")
        res <- list(Q=Q, gQ=.C("fastCalcTotalScattAv",
			layer_start = as.integer(layer_start),			
			layer_end = as.integer(layer_end),
			layers_num = as.integer(layers_num),
			res = as.double(Q),
			Q = as.double(Q), 
			len = as.integer(length(Q)),
			np = as.double(as.vector(t(nanop[x1:(x2-1), ]))),
			nrow = as.integer(x2-x1),
			scale = as.double(sx),
			sigma = as.double(sigma[1]),					
			type=as.integer(type),					
			a1 = as.double(as.vector(c(a1[1], a1[1]))),
			b1 = as.double(as.vector(c(b1[1], b1[1]))),
			a2 = as.double(as.vector(c(a2[1], a2[1]))),
			b2 = as.double(as.vector(c(b2[1], b2[1]))),
			a3 = as.double(as.vector(c(a3[1], a3[1]))),
			b3 = as.double(as.vector(c(b3[1], b3[1]))),
			a4 = as.double(as.vector(c(a4[1], a4[1]))),
			b4 = as.double(as.vector(c(b4[1], b4[1]))),
			a5 = as.double(as.vector(c(a5[1], a5[1]))),
			b5 = as.double(as.vector(c(b5[1], b5[1]))),
			c = as.double(as.vector(c(c[1], c[1]))),	
			scatterLengths=as.double(as.vector(c(scatterLength[1],scatterLength[1]))),
			useN = useN,
			n = as.double(n),
			delta = as.double(delta),
			dr = as.double(dr),
			diam = as.double(as.vector(diam_coef*diam)),
			del = as.double(del),
			eps = as.double(eps),
			PACKAGE="nanop")$res)
      if(kind=="fast" || kind=="fastHist"){
        res <- list(Q=Q, gQ=.C("fastCalcTotalScatt",
			res = as.double(Q),
			Q = as.double(Q), 
			len = as.integer(length(Q)),
			np = as.double(as.vector(t(nanop))),
			nrow = as.integer(nrow(nanop)),
			scale = as.double(sx),
			sigma = as.double(sigma[1]),					
			type=as.integer(type),					
			a1 = as.double(as.vector(c(a1[1], a1[1]))),
			b1 = as.double(as.vector(c(b1[1], b1[1]))),
			a2 = as.double(as.vector(c(a2[1], a2[1]))),
			b2 = as.double(as.vector(c(b2[1], b2[1]))),
			a3 = as.double(as.vector(c(a3[1], a3[1]))),
			b3 = as.double(as.vector(c(b3[1], b3[1]))),
			a4 = as.double(as.vector(c(a4[1], a4[1]))),
			b4 = as.double(as.vector(c(b4[1], b4[1]))),
			a5 = as.double(as.vector(c(a5[1], a5[1]))),
			b5 = as.double(as.vector(c(b5[1], b5[1]))),
			c = as.double(as.vector(c(c[1], c[1]))),
			scatterLengths=as.double(as.vector(c(scatterLength[1],scatterLength[1]))),
			useN = useN,
			n = as.double(n),
			delta = as.double(delta),
			dr = as.double(dr),
			diam = as.double(diam_coef*diam),
			del = as.double(del),
			eps = as.double(eps),
			PACKAGE="nanop")$res)
			}
						
      if(kind == "exact")	
	   res <- list(Q=Q, gQ=.C("calcTotalScatt",
			res = as.double(Q),
			Q = as.double(Q), 
			len = as.integer(length(Q)),
			np = as.double(as.vector(t(nanop))),
			nrow = as.integer(nrow(nanop)),
			scale = as.double(sx),
			sigma = as.double(sigma[1]),					
			type=as.integer(type),					
			a1 = as.double(as.vector(c(a1[1], a1[1]))),
			b1 = as.double(as.vector(c(b1[1], b1[1]))),
			a2 = as.double(as.vector(c(a2[1], a2[1]))),
			b2 = as.double(as.vector(c(b2[1], b2[1]))),
			a3 = as.double(as.vector(c(a3[1], a3[1]))),
			b3 = as.double(as.vector(c(b3[1], b3[1]))),
			a4 = as.double(as.vector(c(a4[1], a4[1]))),
			b4 = as.double(as.vector(c(b4[1], b4[1]))),
			a5 = as.double(as.vector(c(a5[1], a5[1]))),
			b5 = as.double(as.vector(c(b5[1], b5[1]))),
			c = as.double(as.vector(c(c[1], c[1]))),
			scatterLengths=as.double(as.vector(c(scatterLength[1],scatterLength[1]))),	
			useN = useN,
			n = as.double(n),
			delta = as.double(delta),
			PACKAGE="nanop")$res)		
	}

  }
  else { 
    res_gQ_CCSS <- 0
    for(i in 1:nAtomTypes){
	  if( nAt[i] <= 1) 
	    res_gQ_CCSS <- res_gQ_CCSS + rep(0,length(Q))
	  else{
	    if(attributes(nanop)$shape=="ellipse"){
          diam <- 2*max(attributes(nanop)$r)
	      if(tt[i]>0 && min(tt) < 0 )
	        diam <- 2*max(attributes(nanop)$rcore)
        }else if(attributes(nanop)$shape=="box"){
          diam <- 2*sqrt(sum(attributes(nanop)$r^2))
	      if(tt[i]>0 && min(tt) < 0 )
	        diam <- 2*sqrt(sum(attributes(nanop)$rcore^2))         
        }else{
          diam <- 2*attributes(nanop)$r
	      if(tt[i]>0 && min(tt) < 0 )
	        diam <- 2*attributes(nanop)$rcore
        }
		  
	    ff1 <- scatterLength[i]
		ff2 <- scatterLength[i]	  
		layer_start <- 0
		layer_end <- 0
	    for(k in 1:layers_num){
		  layer_start[k] = which(iAtomType[[i]] == min(iATL[[i]][[k]]))
		  layer_end[k] = which(iAtomType[[i]] == max(iATL[[i]][[k]]))
        }
        if(kind=="fast_av")		
          res_gQ_CCSS <- res_gQ_CCSS + .C("fastCalcTotalScattAv",
					layer_start = as.integer(layer_start),			
					layer_end = as.integer(layer_end),
					layers_num = as.integer(layers_num),
					res = as.double(Q),
					Q = as.double(Q), 
					len = as.integer(length(Q)),
					np = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),
					nrow = as.integer(length(iAtomType[[i]])),
					scale = as.double(sx),   # change for X-ray scattering
					sigma = as.double(sigma[i]),
					type=as.integer(type),
					a1 = as.double(as.vector(c(a1[i], a1[i]))),
					b1 = as.double(as.vector(c(b1[i], b1[i]))),
					a2 = as.double(as.vector(c(a2[i], a2[i]))),
					b2 = as.double(as.vector(c(b2[i], b2[i]))),
					a3 = as.double(as.vector(c(a3[i], a3[i]))),
					b3 = as.double(as.vector(c(b3[i], b3[i]))),
					a4 = as.double(as.vector(c(a4[i], a4[i]))),
					b4 = as.double(as.vector(c(b4[i], b4[i]))),
					a5 = as.double(as.vector(c(a5[i], a5[i]))),
					b5 = as.double(as.vector(c(b5[i], b5[i]))),
					c = as.double(as.vector(c(c[i], c[i]))),
					scatterLengths=as.double(as.vector(c(ff1,ff2))),
					useN = useN,
					n = as.double(n),
					delta = as.double(delta),
					dr = as.double(dr),
					diam = as.double(as.vector(diam_coef*diam)),
					del = as.double(del),
					eps = as.double(eps),
					PACKAGE="nanop")$res
					
		if(kind=="fast" || kind=="fastHist")
          res_gQ_CCSS <- res_gQ_CCSS + .C("fastCalcTotalScatt",	
					res = as.double(Q),
					Q = as.double(Q), 
					len = as.integer(length(Q)),
					np = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),
					nrow = as.integer(length(iAtomType[[i]])),
					scale = as.double(sx),   # change for X-ray scattering
					sigma = as.double(sigma[i]),
					type=as.integer(type),
					a1 = as.double(as.vector(c(a1[i], a1[i]))),
					b1 = as.double(as.vector(c(b1[i], b1[i]))),
					a2 = as.double(as.vector(c(a2[i], a2[i]))),
					b2 = as.double(as.vector(c(b2[i], b2[i]))),
					a3 = as.double(as.vector(c(a3[i], a3[i]))),
					b3 = as.double(as.vector(c(b3[i], b3[i]))),
					a4 = as.double(as.vector(c(a4[i], a4[i]))),
					b4 = as.double(as.vector(c(b4[i], b4[i]))),
					a5 = as.double(as.vector(c(a5[i], a5[i]))),
					b5 = as.double(as.vector(c(b5[i], b5[i]))),
					c = as.double(as.vector(c(c[i], c[i]))),
					scatterLengths=as.double(as.vector(c(ff1,ff2))),
					useN = useN,
					n = as.double(n),
					delta = as.double(delta),
					dr = as.double(dr),
					diam = as.double(diam_coef*diam),
					del = as.double(del),
					eps = as.double(eps),
					PACKAGE="nanop")$res			
		
		if(kind=="exact")
          res_gQ_CCSS <- res_gQ_CCSS + .C("calcTotalScatt",	
					res = as.double(Q),
					Q = as.double(Q), 
					len = as.integer(length(Q)),
					np = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),
					nrow = as.integer(length(iAtomType[[i]])),
					scale = as.double(sx),   # change for X-ray scattering
					sigma = as.double(sigma[i]),
					type=as.integer(type),
					a1 = as.double(as.vector(c(a1[i], a1[i]))),
					b1 = as.double(as.vector(c(b1[i], b1[i]))),
					a2 = as.double(as.vector(c(a2[i], a2[i]))),
					b2 = as.double(as.vector(c(b2[i], b2[i]))),
					a3 = as.double(as.vector(c(a3[i], a3[i]))),
					b3 = as.double(as.vector(c(b3[i], b3[i]))),
					a4 = as.double(as.vector(c(a4[i], a4[i]))),
					b4 = as.double(as.vector(c(b4[i], b4[i]))),
					a5 = as.double(as.vector(c(a5[i], a5[i]))),
					b5 = as.double(as.vector(c(b5[i], b5[i]))),
					c = as.double(as.vector(c(c[i], c[i]))),
					scatterLengths=as.double(as.vector(c(ff1,ff2))),
					useN = useN,
					n = as.double(n),
					delta = as.double(delta),
					PACKAGE="nanop")$res		
					
#		cat("CCSS total time", a1, "i",i, "\n")
	  }  
	}
	res_gQ_CS <- 0
    for(i in 2:nAtomTypes){
	  for(j in 1:(i-1)){		
	    if( nAt[i] <= 1 || nAt[j] <= 1 ) 
	      res_gQ_CS <- res_gQ_CS + rep(0,length(Q))
	    else{		
		  if(attributes(nanop)$shape=="ellipse"){
  		    diam <- 2*max(attributes(nanop)$r)
	        if(min(tt) < 0){
	          if(tt[i]*tt[j]<0)
			    diam <- max(attributes(nanop)$rcore) + max(attributes(nanop)$r)
	          if(tt[i]>0 && tt[j]>0)
			    diam <- 2*max(attributes(nanop)$rcore)
		    }
          }else if(attributes(nanop)$shape=="box"){
  		    diam <- 2*sqrt(sum(attributes(nanop)$r^2))
	        if(min(tt) < 0){
	          if(tt[i]*tt[j]<0)
			    diam <- sqrt(sum(attributes(nanop)$rcore^2)) + sqrt(sum(attributes(nanop)$r^2))
	          if(tt[i]>0 && tt[j]>0)
			    diam <- 2*sqrt(sum(attributes(nanop)$rcore^2))
		    }          
          }else{
  		    diam <- 2*attributes(nanop)$r
	        if(min(tt) < 0){
	          if(tt[i]*tt[j]<0)
			    diam <- attributes(nanop)$rcore + attributes(nanop)$r
	          if(tt[i]>0 && tt[j]>0)
			    diam <- 2*attributes(nanop)$rcore
            }
          }			
		  
		  layer_start <- layer_end <- layerS_start <- layerS_end <- 0
	      for(k in 1:layers_num){
		    layerS_start[k] = which(iAtomType[[i]] == min(iATL[[i]][[k]]))  # i==shell
		    layerS_end[k] = which(iAtomType[[i]] == max(iATL[[i]][[k]])) 
    	    layer_start[k] = which(iAtomType[[j]] == min(iATL[[j]][[k]]))   # j==core
		    layer_end[k] = which(iAtomType[[j]] == max(iATL[[j]][[k]]))
          }  

		  if(kind=="fast_av")
		    res_gQ_CS <- res_gQ_CS + .C("fastCalcTotalScattAv_CS",
					layer_start = as.integer(as.vector(layer_start)),			
					layer_end = as.integer(as.vector(layer_end)),
					layerS_start = as.integer(as.vector(layerS_start)),			
					layerS_end = as.integer(as.vector(layerS_end)),
					layers_num = as.integer(layers_num),		
					res = as.double(Q),
					Q = as.double(Q),                
					len = as.integer(length(Q)),
					scale = as.double(sx),   # change for X-ray scattering
					np_mu = as.double(as.vector(t(nanop[iAtomType[[j]], ]))),  #core
					np_nu = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),  #shell
					nrow_mu = as.integer(length(iAtomType[[j]])),
					nrow_nu = as.integer(length(iAtomType[[i]])),
					type=as.integer(type),
					sigma_mu = as.double(sigma[j]),
					sigma_nu = as.double(sigma[i]),
					a1 = as.double(as.vector(c(a1[j],a1[i]))),
					b1 = as.double(as.vector(c(b1[j],b1[i]))),
					a2 = as.double(as.vector(c(a2[j],a2[i]))),
					b2 = as.double(as.vector(c(b2[j],b2[i]))),
					a3 = as.double(as.vector(c(a3[j],a3[i]))),
					b3 = as.double(as.vector(c(b3[j],b3[i]))),
					a4 = as.double(as.vector(c(a4[j],a4[i]))),
					b4 = as.double(as.vector(c(b4[j],b4[i]))),
					a5 = as.double(as.vector(c(a5[j],a5[i]))),
					b5 = as.double(as.vector(c(b5[j],b5[i]))),
					c = as.double(as.vector(c(c[j],c[i]))),	
					scatterLengths=as.double(as.vector(c(scatterLength[j], scatterLength[i]))), 
					useN = useN,
					n = as.double(n),
					delta = as.double(delta),	
					del = as.double(del),					
					PACKAGE="nanop")$res

    	  if(kind=="fast" || kind=="fastHist")
		    res_gQ_CS <- res_gQ_CS + .C("fastCalcTotalScatt_CS",
					res = as.double(Q),
					Q = as.double(Q), 
					len = as.integer(length(Q)),
					scale = as.double(sx),   # change for X-ray scattering
					np_mu = as.double(as.vector(t(nanop[iAtomType[[j]], ]))),  #core
					np_nu = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),  #shell
					nrow_mu = as.integer(length(iAtomType[[j]])),
					nrow_nu = as.integer(length(iAtomType[[i]])),
    				type=as.integer(type),
					sigma_mu = as.double(sigma[j]),
					sigma_nu = as.double(sigma[i]),
					a1 = as.double(as.vector(c(a1[j],a1[i]))),
					b1 = as.double(as.vector(c(b1[j],b1[i]))),
					a2 = as.double(as.vector(c(a2[j],a2[i]))),
					b2 = as.double(as.vector(c(b2[j],b2[i]))),
					a3 = as.double(as.vector(c(a3[j],a3[i]))),
					b3 = as.double(as.vector(c(b3[j],b3[i]))),
					a4 = as.double(as.vector(c(a4[j],a4[i]))),
					b4 = as.double(as.vector(c(b4[j],b4[i]))),
					a5 = as.double(as.vector(c(a5[j],a5[i]))),
					b5 = as.double(as.vector(c(b5[j],b5[i]))),
					c = as.double(as.vector(c(c[j],c[i]))),	
					scatterLengths=as.double(as.vector(c(scatterLength[j], scatterLength[i]))), 
					useN = useN,
					n = as.double(n),
					delta = as.double(delta),	
					dr = as.double(dr),
					diam = as.double(diam_coef*diam),
					del = as.double(del),
					eps = as.double(eps),			
					PACKAGE="nanop")$res		
		
          if(kind=="exact") 
	        res_gQ_CS <- res_gQ_CS + .C("calcTotalScatt_CS",
					res = as.double(Q),
					Q = as.double(Q), 
					len = as.integer(length(Q)),
					scale = as.double(sx),   # change for X-ray scattering
					np_mu = as.double(as.vector(t(nanop[iAtomType[[j]], ]))),  #core
					np_nu = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),  #shell
					nrow_mu = as.integer(length(iAtomType[[j]])),
					nrow_nu = as.integer(length(iAtomType[[i]])),
					type=as.integer(type),
					sigma_mu = as.double(sigma[j]),
					sigma_nu = as.double(sigma[i]),
					a1 = as.double(as.vector(c(a1[j],a1[i]))),
					b1 = as.double(as.vector(c(b1[j],b1[i]))),
					a2 = as.double(as.vector(c(a2[j],a2[i]))),
					b2 = as.double(as.vector(c(b2[j],b2[i]))),
					a3 = as.double(as.vector(c(a3[j],a3[i]))),
					b3 = as.double(as.vector(c(b3[j],b3[i]))),
					a4 = as.double(as.vector(c(a4[j],a4[i]))),
					b4 = as.double(as.vector(c(b4[j],b4[i]))),
					a5 = as.double(as.vector(c(a5[j],a5[i]))),
					b5 = as.double(as.vector(c(b5[j],b5[i]))),
					c = as.double(as.vector(c(c[j],c[i]))),	
					scatterLengths=as.double(as.vector(c(scatterLength[j], scatterLength[i]))), 
					useN = useN,
					n = as.double(n),
					delta = as.double(delta),				
					PACKAGE="nanop")$res		  
		  
#			cat("C/S total time", a1, "\n")								
        }	  
	  }
	}
    res_gQ <- res_gQ_CCSS + res_gQ_CS
    res <- list(Q=Q, gQ=res_gQ)
  }
  res
}
