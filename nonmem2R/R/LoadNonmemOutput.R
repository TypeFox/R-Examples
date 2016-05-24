####################################################################################################
#' Load final parameter values from a nonmem .ext output file.
#'
#' @description
#' Load parameter values from a nonmem .ext output file.
#' For output files from FO, FOCE and IMP only the final parameter values are loaded.
#' For output tables of MCMC and SAEM all parameter values from ITERATIONS>0 will be loaded
#' unless positive.iterations=FALSE in case all the output from all iterations are loaded.
#'
#' For ext-files with multiple table results either only the last table result is loaded (last.table.only=TRUE, default)
#' or all table result are loaded (last.table.only=FALSE).
#' All but the last table results are then returned as sub-lists to the last table result.
#'
#' @param model
#' name of the ext file with or without the .ext extension.
#' model may include full or relative path to the ext file. See examples.
#' @param positive.iterations.only
#' Include only rows with ITERATIONS>0 from MCMC and SEAM table results (default=TRUE)
#' @param last.table.only
#' Include only the last table result for ext files with multiple table results
#' @return
#' Named list including theta, omega, sigma, and ofv. For MCMC output file each object are matrixes.
#' @export
#' @importFrom utils read.table
#' @examples
#' ##### Load the .ext file "run001.ext"
#' # 1) Get path to the example file included in nonmem2R package
#' file1 <- system.file("extdata", "run001.ext", package = "nonmem2R")
#' # 2) Load the file using the extload function
#' extload(file1)
extload <- function(model,positive.iterations.only=T,last.table.only=T) {
  ## Make sure model is without .ext, remove extension if needed
  if(substr(model,nchar(model)-3,nchar(model))==".ext") {
    model=substr(model,1,nchar(model)-4)
  }

	## First check number of TABLE in the ext file
	tmp<-read.table(paste(model,".ext",sep=""),sep="?",header=F,stringsAsFactors = FALSE)
	skip.rows<-grep("TABLE",tmp[,1])
	n.rows<-c(skip.rows[-1],nrow(tmp))-skip.rows-2

	###First do last table
	i<-	length(skip.rows)
	ret<-extload.sub.table(model,skip.rows[i],n.rows[i],positive.iterations.only)

	### Then add other tables as sublists
	if(!last.table.only & length(skip.rows)>1){
		for(i in 1:(length(skip.rows)-1)){
			reti<-extload.sub.table(model,skip.rows[i],n.rows[i],positive.iterations.only)
			eval(parse(text=paste("ret$table",i,"<-reti",sep="")))
		}
	}
	ret
}


####################################################################################################
#' Internal package function
#'
#' @param model
#' model
#' @param skip
#' number of rows to skip
#' @param nrows
#' number of rows to read
#' @param positive.iterations.only
#' positive iterations only
#' @return
#' Named list including theta, omega, sigma, and ofv. For MCMC output file each object are matrixes.
#' @importFrom utils read.table
extload.sub.table <- function(model,skip,nrows,positive.iterations.only) {
	ext<-read.table(file=paste(model,".ext",sep=""),skip=skip,nrows=nrows,header=T)
	### Check for type based on name of last column
	last.col<-colnames(ext)[ncol(ext)]
	type<-match(last.col,c("OBJ","SAEMOBJ","MCMCOBJ"))
	if(is.na(type)){
			cat("Unknown nonmem table type:",last.col,"\n")
			cat("Returning empty list\n")
			ret<-list(theta=NULL,omega=NULL,sigma=NULL)
	}
	else{
		if(type==1){  ## OBJ type, return first row with ITERATION<0
			jtheta<-grep("THETA",colnames(ext))
			jomega<-grep("OMEGA",colnames(ext))
			jsigma<-grep("SIGMA",colnames(ext))
			ii<-which(ext$ITERATION<0)[1]
			theta<-t(ext[ii,jtheta])[,1]

			omegaV<-t(ext[ii,jomega])[,1]
			k<-length(omegaV)
			n<-sqrt(2*k+1/4)-1/2
			omega<-matrix(0,n,n)
			omega[upper.tri(omega,diag=TRUE)]<-omegaV
			omega[lower.tri(omega,diag=TRUE)]<-t(omega)[lower.tri(omega,diag=TRUE)]
			colnames(omega)<-rownames(omega)<-paste("OMEGA",1:n,sep="")

			sigmaV<-t(ext[ii,jsigma])[,1]
			k<-length(sigmaV)
			n<-sqrt(2*k+1/4)-1/2
			sigma<-matrix(0,n,n)
			sigma[upper.tri(sigma,diag=TRUE)]<-sigmaV
			sigma[lower.tri(sigma,diag=TRUE)]<-t(sigma)[lower.tri(sigma,diag=TRUE)]
			colnames(sigma)<-rownames(sigma)<-paste("sigma",1:n,sep="")
			ret<-list(theta=theta,omega=omega,sigma=sigma)

		}
		else{
			jtheta<-c(1,grep("THETA",colnames(ext)))
			jomega<-c(1,grep("OMEGA",colnames(ext)))
			jsigma<-c(1,grep("SIGMA",colnames(ext)))
			ii<-1:nrow(ext)
			if(positive.iterations.only){
				ii<-which(ext$ITERATION>0)
			}
			theta<-ext[ii,jtheta]

			omega<-ext[ii,jomega]

			sigma<-ext[ii,jsigma]

			ret<-list(theta=theta,omega=omega,sigma=sigma)
		}
	}
	ret
}



####################################################################################################
#' Load the covariance matrix from a nonmem .cov output file.
#'
#' @description
#' Load the covariance matrix from a nonmem .cov output file.
#' Either the covariance matrix of all THETA parameters (default) or the covariance matrix of all parameters, THETA, OMEGA and SIGMA.
#'
#'
#' @param model
#' name of the cov file without the .cov extension. model may include full or relative path to the cov file. See examples.
#' @param theta.only
#' return covariance matrix of theta's only (default)
#' @return
#' the covariance matrix
#' @export
#' @importFrom utils read.table
#'
#' @examples
#' ##### Load the .cov file "run001.cov"
#' # 1) Get path to the example file included in nonmem2R package
#' file1 <- system.file("extdata", "run001.cov", package = "nonmem2R")
#' # 2) Load the file using the covload function
#' covload(file1)
covload <- function(model,theta.only=T){
  ## Make sure model is without .cov, remove extension if needed
  if(substr(model,nchar(model)-3,nchar(model))==".cov") {
    model=substr(model,1,nchar(model)-4)
  }
  cov<-read.table(file=paste(model,".cov",sep=""),skip=1,header=T)
	ii<-grep("THETA",cov$NAME)
	jj<-grep("THETA",colnames(cov))

	theta.cov<-t(cov[ii,jj])
	theta.cov
}
