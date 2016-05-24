#berend krebs master thesis
#import/export function



import_file <- function(modus, filename, locnum, popnum){

	#import data, put it together as a list and return it
	p<-array(0,c(locnum,1))
	
	population<-list();
	if(modus == 0){ 
		for(i in 1:popnum){
			tmp<-matrix(as.integer(scan(filename,what="",skip=5+(i-1)*(locnum+2),nlines=locnum)),ncol=3,byrow=TRUE);
			#colnames(tmp)<-(c("Index","Max","Hap"));
			tmp<-matrix(append(tmp,p),ncol=4)
			population[[i]]<-tmp;
		}
	} else {
		hapnum<-as.integer(scan(filename,what="",skip=5,nlines=1)[3]);
		for(i in 1:popnum){
			tmp<-matrix(as.integer(scan(filename,what="",skip=5+(i-1)*(locnum+2),nlines=locnum)),ncol=(hapnum+3),byrow=TRUE);
			population[[i]]<-tmp;
		}
	}

	print("Data successfully imported!");



	return(population)
}
