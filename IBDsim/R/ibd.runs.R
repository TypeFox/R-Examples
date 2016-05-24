sap.segments <- #the new ibd.runs!!!
function(twostrands_list, sap) { #genedrops for several chromosomes. sap = single allele pattern	
	do.call(rbind, lapply(twostrands_list, sap.finder, sap=sap))
}

sap.finder <-
function(x, sap) {#x a list of haplo-objects (i.e. each element is a list of two matrices); sap = single allele pattern
	zero = sap[['0']]; one = sap[['1']]; two = sap[['2']]; atm1 = sap[['atmost1']]; atl1 = sap[['atleast1']]; not1=sap[['not1']]
	stopifnot(!is.null(id1 <- c(two, atl1, one)[1]))
	sap.indivs = c(zero,one,two,atm1,atl1,not1)
	breaks = unlist(lapply(unlist(x[sap.indivs], recursive=FALSE), function(m) m[-1,1]))
	breaks = c(0, .sortDouble(breaks[!duplicated(breaks)]))
	running_alleles = numeric(0); res=matrix(nrow=0,ncol=3); starts=numeric(); 
	ALLalleles = list(); ALLalleles[sap.indivs] <- lapply(x[sap.indivs], .getAlleles, posvec=breaks)
	for(i in seq_along(breaks)) { 
		a = ALLalleles[[id1]][1, i]; b = ALLalleles[[id1]][2, i]
		for (A in running_alleles[!running_alleles %in% c(a,b)]) {#reached end of segment, write to result matrix
			res = rbind(res, c(starts[A], breaks[i], A))
			running_alleles = running_alleles[running_alleles != A]
		}
		for(al in {if(a==b) a else c(a,b)}) {
			ibd_nr = numeric(length(x))
			ibd_nr[sap.indivs] = unlist(lapply(ALLalleles[sap.indivs], function(als) sum(als[,i]==al)))
			test = all(ibd_nr[two]==2) && all(ibd_nr[atl1]>0) && all(ibd_nr[one]==1) && all(ibd_nr[atm1]<2) && all(ibd_nr[zero]==0) && all(ibd_nr[not1]!=1)  #longer but faster
		
			ind = match(al, running_alleles, nomatch=0) #
			if(ind>0 && !test)  {#cat(al,"pos",breaks[i],"stop(fail)\n")
				res = rbind(res, c(starts[al], breaks[i], al))
				running_alleles = running_alleles[-ind] #cat(running_alleles, "running\n")
			}
			else if(ind==0 && test) {#cat(al,"pos",breaks[i],"start\n")
				starts[al]=breaks[i]; 
				running_alleles=c(running_alleles, al)
			}
		}
	}
	for (A in running_alleles) res = rbind(res, c(starts[A], attr(x,'length_Mb'), A)) #last entries
	
	if(!is.null(disloc <- attr(x, 'dis.locus'))) 
		disreg = as.numeric(res[,1] <= attr(x, 'dis.locus') & res[,2] > attr(x, 'dis.locus') & res[,3] == attr(x, 'dis.allel'))
	else disreg = rep.int(0,nrow(res))
	
	chrom = rep.int(attr(x, 'chromosome'), nrow(res))
	res = cbind(chrom, res, disreg)
	colnames(res) = c('chrom', 'start', 'end', 'allele', 'disreg') 
	res
}


.getAlleles = function(chromdata, posvec) {
	str1 = chromdata[[1]]; str2 = chromdata[[2]]
	posvec[posvec<0]=0
    index1 = findInterval(posvec, str1[,1])
    index2 = findInterval(posvec, str2[,1])
    rbind(str1[index1, 2], str2[index2, 2])
}

