###------DST fmt functions--------------#####

### from matlab codes by Philippe Smets. FMT = Fast Mobius Transform
### depend: package library wle ####

##' @export
mtobetp  <-  function(InputVec) {
    # computing BetP on the signal points from the m vector (InputVec) out = BetP
    # vector beware: not optimize, so can be slow for >10 atoms
    
    # the length of the power set, f
    mf = length(InputVec)
    # the number of the signal point clusters
    natoms = round(log2(mf))
    if (2^natoms == mf) {
        if (InputVec[1] == 1) {
            # bba of the empty set is 1
            stop("warning: all bba is given to the empty set, check the frame\n")
            out = matrix(1, 1, natoms)/natoms
        } else {
            betp = matrix(0, 1, natoms)
            for (i in 2:mf) {
                # x , the focal sets InputVec the dec2bin form
                x = dec2bin(i - 1, natoms)
                # m_i is assigned to all the signal points equally
                betp = betp + InputVec[i]/sum(x) * x
            }
            out = betp/(1 - InputVec[1])
        }
        return(out)
    } else {
        stop("Error: the length of the InputVec vector should be power set of 2\n")
    }
}

##' @export
mtob  <-  function(InputVec) {
    
    # comput InputVec from m to b function.  belief function + m(emptset) InputVec = m
    # vector out = b vector
    
    mf = length(InputVec)
    natoms = round(log2(mf))
    if (2^natoms == mf) {
        # InputVec=1:4/sum(1:4) browser()
        for (myi in 1:natoms) {
            ## the first time, sum the items -2 (each focal sets that includes w1, sum
            ## the item just not include w1) the second time, sum the items -4(each
            ## focal sets that includes w2, sum the item just not include w2) at this
            ## time, the item just not include w2 has the item not include (w1,w2)
            ## continue until for that just not include wK
            i124 = 2^(myi - 1)
            i842 = 2^(natoms + 1 - myi)
            i421 = 2^(natoms - myi)
            InputVec = matrix(InputVec, i124, i842)
            InputVec[, (1:i421) * 2] = InputVec[, (1:i421) * 2] + InputVec[, (1:i421) * 2 - 
                1]
            # print(InputVec)
        }
        out = matrix(InputVec, 1, mf)
        return(out)
    } else {
        stop("Error: the length of the InputVec vector should be power set of 2\n")
    }
}


##' @export
mtopl  <-  function(InputVec) {
    # computing from m to pl in = m vector out = pl vector
    
    InputVec = mtob(InputVec)
    out = btopl(InputVec)
    return(out)
}


##' @export
mtoq <- function(InputVec) {
    # computing FMT from m to q.
	# input:
	#    InputVec = m vector 
	# output:
	#    out = q vector
    mf = length(InputVec)
    natoms = round(log2(mf))
    if (2^natoms == mf) {
        for (i in 1:natoms) {
            i124 = 2^(i - 1)
            i842 = 2^(natoms + 1 - i)
            i421 = 2^(natoms - i)
            InputVec = matrix(InputVec, i124, i842)
            InputVec[, (1:i421) * 2 - 1] = InputVec[, (1:i421) * 2 - 1] + InputVec[, (1:i421) * 2]
        }
        out = matrix(InputVec, 1, mf)
		return(out)
    } else {
        stop("ACCIDENT in mtoq: length of input vector not OK: should be a power of 2\n")
    }
} 

##' @export
mtonm <- function(InputVec) {
    # transform bbm m into normalized bbm
	# make the InputVec in the closed world
    if (InputVec[1] < 1) {
        out = InputVec/(1 - InputVec[1])
    }
    out[1] = 0
    return(out)
} 

##' @export
mtobel <- function(InputVec) {
    out = mtob(mtonm(InputVec))
	return(out)
}

##' @export
mtov <- function(InputVec) {
    #  not checked 
    # computing FMT from m to v. 
	# InputVec = m vector out = v vector
    out = btov(mtob(InputVec))
	return(out)
}

##' @export
mtow <- function(InputVec) {
    # not checked
    # computing FMT from m to w. 
	## InputVec = m vector
	## out = w vector
    out = qtow(mtoq(InputVec))
	return(out)
}

##' @export
btopl  <-  function(InputVec) {
    # compute pl from b InputVec ,
	## InputVec vector f*1 
	## out = pl
    
    mf = length(InputVec)
    natoms = round(log2(mf))
    if (2^natoms == mf) {
        # browser()
        InputVec = InputVec[mf] - InputVec[seq(length(InputVec), 1, by = -1)]
        # InputVec[1] = 0;
        out = InputVec
        return(out)
    } else {
        stop("Error: the length of the InputVec vector should be power set of 2\n")
    }
}

##' @export
btom <- function(InputVec) {
    # computing FMT from b to m.  
	# input: InputVec = b vector
	## out = m vector
    
    mf = length(InputVec)
    natoms = round(log2(mf))
    if (2^natoms == mf) {
        for (i in 1:natoms) {
            i124 = 2^(i - 1)
            i842 = 2^(natoms + 1 - i)
            i421 = 2^(natoms - i)
            InputVec = matrix(InputVec, i124, i842)
            InputVec[, (1:i421) * 2] = InputVec[, (1:i421) * 2] - InputVec[, (1:i421) * 2 - 1]
        }
        out = matrix(InputVec, 1, mf)
        return(out)
    } else {
        stop("ACCIDENT InputVec btom: length of InputVecput vector not OK: should be a power of 2\n")
    }
}

##' @export
btoq <- function(InputVec) {
    # computing FMT from b to q. Compute thru pl InputVec = b vector out = q vector
    out = pltoq(btopl(InputVec))
    return(out)
} 

##' @export
btobel <- function(InputVec) {
    # computing bel from b 
	# input: bufn = b
	# out = bel
    
    if (InputVec[1] == 1) 
        stop("ACCIDENT in btobel: you try to normalize with a 1 on m\n") else {
        
        mf = length(InputVec)
        natoms = round(log2(mf))
        
        if (2^natoms == mf) {
            InputVec = (InputVec - InputVec[1])
            InputVec[1] = 0
            out = InputVec
            return(out)
        } else {
            stop("ACCIDENT InputVec btobel: length of InputVecput vector not OK: should be a power of 2\n")
        }
    }
}

##' @export
pltoq <- function(InputVec) {
    # computing FMT from pl to q. 
	# input: InputVec = pl vector 
	# out = q vector
    
    out = abs(btom(InputVec))
    out[1] = 1
    return(out)
}

##' @export
pltob <- function(InputVec) {
    # compute b from pl 
	# input:InputVec = pl 
	# out = b
    
    mf = length(InputVec)
    natoms = round(log2(mf))
    if (2^natoms == mf) {
        InputVec = 1 - InputVec[length(InputVec):1]
        out = InputVec
		return(out)
    } else {
        stop("ACCIDENT InputVec btopl: length of InputVecput vector not OK: should be a power of 2\n")
    }

	return(out)
}

##' @export
pltobel <- function(InputVec) {
    
    out = pltob(InputVec)
    if (out[1] < 1) {
        out = out/(1 - out[1])
    }
    out[1] = 0
    return(out)
}

##' @export
pltom <- function(InputVec) {
    # compute m from pl 
	# input InputVec = pl out = m
    
    out = btom(pltob(InputVec))
    return(out)
}


##' @export
beltob <- function(InputVec){

   out = InputVec;
   return(out)
}


##' @export
beltom <- function(InputVec){
  out = btom(InputVec);
  return(out)
}

##' @export
beltopl <- function(InputVec){
  out = btopl(InputVec);
  return(out)
}

##' @export
beltoq <- function(InputVec){
  out = btoq(InputVec);
  return(out)
}

##' @export
qtom <- function(InputVec){
	# computing FMT from q to m.
	# in = q vector
	# out = m vector

	lm = length(InputVec);
	natoms = round(log2(lm)); 		
	if(2^natoms == lm){ 
	for (step in 1:natoms){ 
		i124 = 2^(step-1); 			
		i842 = 2^(natoms+1-step); 	
		i421 = 2^(natoms - step); 	
		InputVec = matrix(InputVec,i124,i842);
		InputVec[,(1:i421)*2-1] = InputVec[,(1:i421)*2-1] - InputVec[,(1:i421)*2];
	}	
	    out = matrix(InputVec,1,lm);
	}
	else{
		stop('ACCIDENT in qtom: length of input vector not OK: should be a power of 2\n');
	}
	return(out)
}

##' @export
qtow <- function(InputVec){

	# computing FMT from q to w. Use algorithm qtom on log q
	# in = q vector
	# out = w vector

	lm = length(InputVec);
	natoms = round(log2(lm)); 		
	if(2^natoms == lm){ 
		if(InputVec[lm] > 0){
			out = exp(-(qtom(log(InputVec))));
			out[lm] = 1;
		}else{
			cat('Accident in qtow: algorithm works only if q(lm) > 0\n')
			cat('add an epsilon to m(lm)\n')
			cat('No garantee it is really OK\n')
			mini = 1;
			for (i in 1:lm){
				if (InputVec[i] > 0){
					mini = min(mini,InputVec[i]);
				}
			}
			mini = mini / 10000000000;
			for(i in 1:lm){
			  InputVec[i] = max(InputVec[i],mini);
			}
			out = exp(-(qtom(log(InputVec))));
			out[lm] = 1;
		}
	}else{
		stop('ACCIDENT in qtow: length of input vector not OK: should be a power of 2\n');
	}
    return(out)
}
	
##' @export
wtom <- function(InputVec){
	# computing FMT from w to m. 
	# in = w vector
	# out = m vector
	out = qtom(wtoq(InputVec));
	return(out)
}

##' @export
wtoq <- function(InputVec){
	# computing FMT from w to q. use algortihm mtoq
	# in = w vector
	# out = q vector

	lm = length(InputVec);
	natoms = round(log2(lm)); 		
	if(2^natoms == lm){ 
		if(min(InputVec)>0){
			out = prod(InputVec)/exp(mtoq(log(InputVec)));
		}else{
    		stop('accident in wtoq: one of the weigths are non positive\n')
		}
	}else{
		stop('ACCIDENT in wtoq: length of input vector not OK: should be a power of 2\n')
	}

	return(out)
}


##' @export
btov <- function(InputVec){
	# computing FMT from b to v. Use algorithm btom on log b
	# in = b vector
	# out = v vector

	lm = length(InputVec);
	natoms = round(log2(lm)); 		
	if(2^natoms == lm){
		if(InputVec[lm] > 0){
			out = exp(-(btom(log(InputVec))));
			out[1] = 1;
		}else{
			stop('Accident in btov: algorithm works only if b(lm) > 0\n')
			stop('I nevertheless try to compute it, I add an epsilon to m(lm)\n')
			stop('No garantee it is really OK\n')
			mini=1;
			for(i in 1:lm){
				if(InputVec[i]>0){
					mini=min(mini,InputVec[i]);
				}
			}
			mini=mini/10000000000;
			for(i in 1:lm){
				InputVec[i]=max(InputVec[i],mini);
			}
			out = exp(-(btom(log(InputVec))));
			out[1] = 1;
		}
	}else{
		stop('ACCIDENT in btov: length of input vector not OK: should be a power of \n')
	 }
	return(out)
}


##' @export
vtom <- function(InputVec){
	# computing FMT from v to m. 
	# in = v vector
	# out = m vector
	 out = btom(vtob(InputVec));
	 return(out)
}

##' @export
vtob <- function(InputVec){
	# computing FMT from v to b. use algortihm mtoq
	# in = v vector
	# out = b vector

	lm = length(InputVec);
	natoms = round(log2(lm)); 		
	if(2^natoms == lm){ 
		if(min(InputVec)>0){
			out = prod(InputVec)/exp( mtob(log(InputVec)));
		}else{
		 stop('accident in vtob: one of the weigths are non positive\n')
		}
	}else{
		stop('ACCIDENT in vtob: length of input vector not OK: should be a power of 2\n')
	}

	return(out)
}

