calcPDF <- function(nanop, dr=0.01, minR=0.01, maxR=20, 
                    scatterLength=NA, scatterFactor=NA,
                    type="neutron", Qmin=1e-16) {
					
					
		
  if(type == "neutron") 
    type <- 0
  else 
    type <- 1 
	
  nAtomTypes <- attributes(nanop)$nAtomTypes
  if (is.na(scatterLength[1]))
    scatterLength <- attributes(nanop)$scatterLength
  if (length(scatterLength) < nAtomTypes){
    warning("vector <scatterLength> length is less than number of atom types nAtomTypes\n", immediate. = TRUE)
	scatterLength <- rep(scatterLength[1], nAtomTypes)
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
  
  
  atomType <- attributes(nanop)$atomType
  layers_num <- length(attributes(nanop)$layer_end)
  tt <- 1:max(atomType)
  if(min(atomType) < 0)
    tt <- c(tt, -1:min(atomType))
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
  for(i in 1:layers_num){
    fAv[i] <- 0
	NN[i] <- 0
    if(type==0){
	  ff <- scatterLength
	}else{
	  q4 = (Qmin/(4*pi))^2
      ff = a1*exp(-b1*q4) + a2*exp(-b2*q4) + a3*exp(-b3*q4) + 
		        a4*exp(-b4*q4) + a5*exp(-b5*q4)+c		
	}
	
	for(j in 1:nAtomTypes){
	  fAv[i] <- fAv[i] + length(iATL[[j]][[i]])*ff[j]
	  NN[i] <- NN[i] + length(iATL[[j]][[i]])
	}
    fAv[i] <- fAv[i]/NN[i]
	sx[i] <- 1 / (NN[i]*fAv[i]^2)  # 1/N<f>^2 for given layer
  }
   
  
  r <- seq(minR, maxR, by = dr)
  if(length(nanop) == 0)
    return(list(r=r,gr=rep(0,length(r))))
	
  if(nAtomTypes==1) {
  
	layer_start <- attributes(nanop)$layer_start	
    layer_end <- attributes(nanop)$layer_end
	layers_num <- length(attributes(nanop)$layer_end)
	nrow_tot <- layer_end-layer_start+1
		
    ret <- list(r=r, gr=.C("calcPDF",
					 res = as.double(rep(0,length(r))),
					 r = as.double(r), 
					 len = as.integer(length(r)),
					 layer_start = as.integer(layer_start),			
                     layer_end = as.integer(layer_end),
			         layers_num = as.integer(layers_num),
					 np = as.double(as.vector(t(nanop))),
					 nrow = as.integer(nrow(nanop)),
					 scale = as.double(as.vector(sx)),
					 a1 = as.double(as.vector(c(scatterFactor$a1[1], scatterFactor$a1[1]))),
					 b1 = as.double(as.vector(c(scatterFactor$b1[1], scatterFactor$b1[1]))),
					 a2 = as.double(as.vector(c(scatterFactor$a2[1], scatterFactor$a2[1]))),
					 b2 = as.double(as.vector(c(scatterFactor$b2[1], scatterFactor$b2[1]))),
					 a3 = as.double(as.vector(c(scatterFactor$a3[1], scatterFactor$a3[1]))),
					 b3 = as.double(as.vector(c(scatterFactor$b3[1], scatterFactor$b3[1]))),
					 a4 = as.double(as.vector(c(scatterFactor$a4[1], scatterFactor$a4[1]))),
					 b4 = as.double(as.vector(c(scatterFactor$b4[1], scatterFactor$b4[1]))),
					 a5 = as.double(as.vector(c(scatterFactor$a5[1], scatterFactor$a5[1]))),
					 b5 = as.double(as.vector(c(scatterFactor$b5[1], scatterFactor$b5[1]))),
					 c = as.double(as.vector(c(scatterFactor$c[1], scatterFactor$c[1]))),	
					 scatterLengths=as.double(as.vector(c(scatterLength[1],scatterLength[1]))),
					 type=as.integer(type),
					 Qmin=Qmin, 
					 dr = as.double(dr),
					 minR = as.double(minR),
					 PACKAGE="nanop")$res,
				facs=0)
					 	   
  }
   else {
    ffav <- 0
    res_gr_CCSS <- list()  
    for(i in 1:nAtomTypes){
	  if( nAt[i] <= 1 ) 
	    res_gr_CCSS[[i]] <- rep(0,length(r))
	  else{
	  	ff1 <- scatterLength[i]
		ff2 <- scatterLength[i]	  
		layer_start <- 0
		layer_end <- 0
	    for(k in 1:layers_num){
		  layer_start[k] = which(iAtomType[[i]] == min(iATL[[i]][[k]]))
		  layer_end[k] = which(iAtomType[[i]] == max(iATL[[i]][[k]]))
        }
	    res_gr_CCSS[[i]] <- .C("calcPDF",
								 res = as.double(rep(0,length(r))),
								 r = as.double(r), 
								 len = as.integer(length(r)),
								 layer_start = as.integer(layer_start),			
								 layer_end = as.integer(layer_end),
								 layers_num = as.integer(layers_num), 
								 np = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),
								 nrow = as.integer(length(iAtomType[[i]])),
								 scale = as.double(as.vector(sx)),   # change for X-ray scattering
								 a1 = as.double(as.vector(c(scatterFactor$a1[i], scatterFactor$a1[i]))),
								 b1 = as.double(as.vector(c(scatterFactor$b1[i], scatterFactor$b1[i]))),
								 a2 = as.double(as.vector(c(scatterFactor$a2[i], scatterFactor$a2[i]))),
								 b2 = as.double(as.vector(c(scatterFactor$b2[i], scatterFactor$b2[i]))),
								 a3 = as.double(as.vector(c(scatterFactor$a3[i], scatterFactor$a3[i]))),
								 b3 = as.double(as.vector(c(scatterFactor$b3[i], scatterFactor$b3[i]))),
								 a4 = as.double(as.vector(c(scatterFactor$a4[i], scatterFactor$a4[i]))),
								 b4 = as.double(as.vector(c(scatterFactor$b4[i], scatterFactor$b4[i]))),
								 a5 = as.double(as.vector(c(scatterFactor$a5[i], scatterFactor$a5[i]))),
								 b5 = as.double(as.vector(c(scatterFactor$b5[i], scatterFactor$b5[i]))),
								 c = as.double(as.vector(c(scatterFactor$c[i], scatterFactor$c[i]))),
								 scatterLengths=as.double(as.vector(c(ff1,ff2))),
								 type=as.integer(type),
								 Qmin=Qmin, 								 
								 dr = as.double(dr),
								 minR = as.double(minR),
								 PACKAGE="nanop")$res

	  }			   
	}  
    res_gr_CS <- list()
	
    for(i in 2:nAtomTypes){
	  for(j in 1:(i-1)){
		k <- ceiling(i*(i-3)/2+j+1)
		if ( nAt[i] <= 1 || nAt[j] <= 1 )
		  res_gr_CS[[k]] <-  rep(0,length(r))
	    else{
		  layer_start <- layer_end <- layerS_start <- layerS_end <- 0
	      for(m in 1:layers_num){
		    layerS_start[m] = which(iAtomType[[i]] == min(iATL[[i]][[m]]))  # i==shell
		    layerS_end[m] = which(iAtomType[[i]] == max(iATL[[i]][[m]])) 
    	    layer_start[m] = which(iAtomType[[j]] == min(iATL[[j]][[m]]))   # j==core
		    layer_end[m] = which(iAtomType[[j]] == max(iATL[[j]][[m]]))
          }   
	      res_gr_CS[[k]] <- .C("calcPDF_CS",
				  layer_start = as.integer(as.vector(layer_start)),			
				  layer_end = as.integer(as.vector(layer_end)),
				  layerS_start = as.integer(as.vector(layerS_start)),			
				  layerS_end = as.integer(as.vector(layerS_end)),
				  layers_num = as.integer(layers_num),		  
				  res = as.double(rep(0,length(r))),
				  r = as.double(r), 
				  len = as.integer(length(r)),
				  scale = as.double(as.vector(sx)),   # change for X-ray scattering
				  np_mu = as.double(as.vector(t(nanop[iAtomType[[j]], ]))),  #core
				  nrow_mu = as.integer(length(iAtomType[[j]])),
				  np_nu = as.double(as.vector(t(nanop[iAtomType[[i]], ]))),  #shell
				  nrow_nu = as.integer(length(iAtomType[[i]])),		  
				  a1 = as.double(as.vector(c(scatterFactor$a1[j],scatterFactor$a1[i]))),
				  b1 = as.double(as.vector(c(scatterFactor$b1[j],scatterFactor$b1[i]))),
				  a2 = as.double(as.vector(c(scatterFactor$a2[j],scatterFactor$a2[i]))),
				  b2 = as.double(as.vector(c(scatterFactor$b2[j],scatterFactor$b2[i]))),
				  a3 = as.double(as.vector(c(scatterFactor$a3[j],scatterFactor$a3[i]))),
				  b3 = as.double(as.vector(c(scatterFactor$b3[j],scatterFactor$b3[i]))),
				  a4 = as.double(as.vector(c(scatterFactor$a4[j],scatterFactor$a4[i]))),
				  b4 = as.double(as.vector(c(scatterFactor$b4[j],scatterFactor$b4[i]))),
				  a5 = as.double(as.vector(c(scatterFactor$a5[j],scatterFactor$a5[i]))),
				  b5 = as.double(as.vector(c(scatterFactor$b5[j],scatterFactor$b5[i]))),
				  c = as.double(as.vector(c(scatterFactor$c[j],scatterFactor$c[i]))),	
				  scatterLengths=as.double(as.vector(c(scatterLength[j], scatterLength[i]))),
				  type=as.integer(type),
				  Qmin=Qmin, 
				  dr = as.double(dr),
				  minR = as.double(minR),
				  PACKAGE="nanop")$res		

        }				  
	  }
	}	
	res_gr <-  0
	for(i in 1:nAtomTypes){
	  res_gr <- res_gr + res_gr_CCSS[[i]]
	}
	for(i in 2:nAtomTypes){
      for(j in 1:(i-1)){
	    k <- ceiling(i*(i-3)/2+j+1)
		res_gr <- res_gr + res_gr_CS[[k]]
	  }
	}
	
    ret <- list(r=r,
	            gr=res_gr,
                gr_CCSS = res_gr_CCSS,
				gr_CS = res_gr_CS, 
				facs = 1)
    
  }
  
  ret
}
