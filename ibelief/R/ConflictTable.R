##' Computing the conflict table
##'
##' Computing the table of conflict for \eqn{nbexperts} masses  and \eqn{natoms = round(\log2(lm))} classes. 
##' This function gives the conflict focal set combinations for the \eqn{nbexperts} masses. The focal sets are labeled in natural order, e.g, 
##' number 2 denotes \eqn{\omega_1}, and number 4 donoets \eqn{\{\omega_1,\omega_2\}} if the discernment frame is \eqn{\{\omega_1,\omega_2,\ldots,\omega_n\}}. Note that only one case of 
##' conflict is given. For example,  if expert 1 says 3, and expert 2 says 2 the function returns \code{matrix(c(2,3),,1)}
##'  and if expert 1 says 2, and expert 2 says 3 the function also returns \code{matrix(c(2,3),,1)}. 
##' @export
##'
##' @param lm The length of the power set of the discernment frame, i.e., \eqn{2^{natoms}}
##' @param nbexperts The number of experts (masses)
##' @return  Matrix with \eqn{nbexperts} rows and number of conflict focal set combinations columns. 
##'
##' @seealso \code{\link{PCR6}}, \code{\link{decisionDST}} 
##' @examples
##' ## The conflict table for two experts in a discernment frame with three elements
##'  ConflictTable(2^3,2) 
##' ##The conflict table for three experts in a discernment frame with four elements
##'  ConflictTable(2^4,3) 
##'
ConflictTable <- function(lm,nbexperts){
    # depending programe
	AddConflict <- function(natoms,Int,ind,in1,indAjout,experts,nbexperts){
		# recursive function for add the conflict for 3 and more experts

		if(experts<(nbexperts+1)){
			for(k in 2:(2^natoms-1)){
				Int2=intersect(Int,ind[[k]]);
				if(!length(Int2)){
					m=rbind(indAjout, k, matrix(0,nbexperts-experts,1));
					in1=cbind(in1,m[order(m[,1]),]);
				}else{
					expertsRec=experts+1;
					indAjoutRec=rbind(indAjout,k);
					in1=AddConflict(natoms,Int,ind,in1,indAjoutRec,expertsRec,nbexperts);
				}
			}
		}
		out=in1;
	}

	natoms = round(log2(lm)); 		
	ind=list()
	if(2^natoms == lm){ 
		ind[[1]]=c(1);;
		ind[[2]]=c(2);
		mystep =3;
		while(mystep<2^natoms){
			ind[[mystep]]=c(mystep);
			mystep=mystep+1;
			indatom=mystep;
			for(mystep2 in 2:(indatom-2)){
				ind[[mystep]]=sort(unique(c(ind[[indatom-1]],ind[[mystep2]])));
				mystep=mystep+1;
		    }
		}
	
		out=c();
		temp=1;
		for(i in 2:(2^natoms-2))
			for (j in (i+1):(2^natoms-1)){
				Int=intersect(ind[[i]],ind[[j]]);
				if(!length(Int)){
					m=rbind(i,j,matrix(0,nbexperts-2,1));
					out=cbind(out,m[order(m[,1]),]);
					#browser()
				
				}else{
					experts=3;
					indAjout=rbind(i,j);
					out=AddConflict(natoms,Int,ind,out,indAjout,experts,nbexperts);
				    #browser()
				}
			}
	}else{
		stop('ACCIDENT in ConflictTable: length of input vector not OK: should be a power of 2\n');
    }
      # browser()
      mm=unique(t(out));
	  colnum=ncol(mm);

	  ## order the matrix according to column 1,2,...
	  temp="mm[,1]";
	  for(j in 2:(colnum)){
	    temp1=paste("mm[,",eval(j),"]",sep="")
	    temp=paste(temp,temp1,sep=",")
	  }
	  temp2=paste("order(",temp,")",sep="")
	  myorder=eval(parse(text=temp2))
	  mm=mm[myorder,];
	  out=t(mm);
	  return(matrix(out,nbexperts))

}
