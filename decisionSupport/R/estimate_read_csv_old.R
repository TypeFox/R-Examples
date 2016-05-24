#
# file: estimate_read_csv_old.R
#
# This file is part of the R-package decisionSupport
#
# Authors:
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF)
#	http://www.worldagroforestry.org
#
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' @include estimate.R
NULL
###############################################################################################
# estimate_read_csv_old(fileName, strip.white=TRUE, ...)
##############################################################################################
#' Read an Estimate from CSV - File (deprecated standard).
#'
#' \code{estimate_read_csv_old} reads an estimate from CSV file(s) according to the deprecated 
#' standard. This function is for backward compatibility only.
#' @rdname estimate_read_csv
#' @details 
#' \subsection{Deprecated input standard (\code{estimate_read_csv_old})}{
#'    File name structure of the correlation file: \code{<marginal-filename>.csv_correlations.csv}\cr
#' }
#' @seealso \code{\link{estimate_read_csv}}, \code{\link[utils]{read.csv}}, \code{\link{estimate}}
#' @export
estimate_read_csv_old<-function(fileName, strip.white=TRUE, ...){
  .Deprecated("estimate_read_csv")
	marginal<-NULL
	#	correlation_matrix<-NULL
	marginalFilename=fileName
	# Read marginal data:
	marginal<-read.csv(marginalFilename, strip.white=strip.white, stringsAsFactors=FALSE, ...)
	marginal<-subset(marginal,variable!="")
	marginal<-data.frame(marginal,row.names="variable")
	# Transform number format: remove ',' as thousand separator:
	marginal[,"lower"]<-as.numeric(as.numeric(gsub(",","", as.character(marginal[,"lower"]))))
	marginal[,"upper"]<-as.numeric(as.numeric(gsub(",","", as.character(marginal[,"upper"]))))


	# Read correlation file:
	if( length(marginal$distribution[marginal$distribution=="correlated"]) > 0 ){
		correlated_variables<-estimate_read_csv_old_correlation(paste(fileName,"_correlations.csv",sep=""))
	}

	# Merge marginal information and correlation information (ToDo: check):
	if( !setequal(intersect(row.names(marginal), row.names(correlated_variables$marginal)),
					 row.names(subset(x=marginal,subset=(distribution=="correlated")))) )
		stop("Variables in marginal file indicated as correlated are not the same as in correlation file.")


	marginal<-subset(x=marginal,subset=(distribution!="correlated"))
	marginalCor<-correlated_variables$marginal
#	mergedmarginal<-data.frame(merge(marginal,marginalCor,by=c("row.names",intersect(names(marginal),names(marginalCor)))),row.names="Row.names")
	mergedmarginal<-data.frame(merge(marginal,marginalCor,by=c("row.names",intersect(names(marginal),names(marginalCor))),all.y=TRUE),row.names="Row.names")
	marginal<-rbind(subset(x=marginal,subset=(distribution!="correlated")),
							mergedmarginal)

	# Transform old to new syntax
	marginal$distribution[marginal$distribution=="constant"]<-"const"
	marginal$distribution[marginal$distribution=="normal"]<-"norm"
	marginal$distribution[marginal$distribution=="pos_normal"]<-"posnorm"
	marginal$distribution[marginal$distribution=="normal_0_1"]<-"tnorm_0_1"

	# Return:
	as.estimate(marginal,correlation_matrix=correlated_variables$correlation_matrix)
}

###############################################################################################
# estimate_read_csv_old_correlation(inFileName)
# ToDo: review documentation (if pre and postconditions are correct)
##############################################################################################
# Read correlated variables:
estimate_read_csv_old_correlation<-function(inFileName){
	# Initialization:
	correlationMatrix<-NULL
	marginal<-NULL

	tabl<-read.csv(inFileName, header=FALSE,stringsAsFactors=FALSE)
	block_starts<-which(tabl[,1]=="variable")
	if(length(block_starts)==1) block_ends<-nrow(tabl) else
		block_ends<-c(block_starts[2:length(block_starts)]-1,nrow(tabl))

	for(i in 1:length(block_starts)){
		col_names<-as.character(unlist(tabl[block_starts[i],]))
		col_names<-sapply(tabl[block_starts[i],],as.character)
		block<-tabl[(block_starts[i]+1):block_ends[i],]
		colnames(block)<-as.character(col_names)
		block[,"lower"]<-as.numeric(as.numeric(gsub(",","", as.character(block[,"lower"]))))
		block[,"upper"]<-as.numeric(as.numeric(gsub(",","", as.character(block[,"upper"]))))
		block<-block[which(!is.na(block[,"upper"])),]

		block<-data.frame(block,row.names="variable")

		marginalBlock<-block[!(names(block) %in% row.names(block))]
		correlationMatrixBlock<-block[(names(block) %in% row.names(block))]
		# Replace '#' by the transposed numeric value
		correlationMatrixBlock[which(correlationMatrixBlock=="#", arr.ind=TRUE)]<-
			t(correlationMatrixBlock)[which(correlationMatrixBlock=="#", arr.ind=TRUE)]

		marginal<-rbind(marginal,marginalBlock)
		if( is.null(correlationMatrix) ){
			correlationMatrix<-correlationMatrixBlock
		}
		else{
			correlationMatrix	<- merge(correlationMatrix, correlationMatrixBlock, by = "row.names", all = TRUE)
			correlationMatrix[is.na(correlationMatrix)] <- 0
			correlationMatrix	<- data.frame(correlationMatrix, row.names="Row.names")
		}
	}
	correlationMatrix<-data.matrix(correlationMatrix)[row.names(marginal), row.names(marginal)]
	# Return:
	as.estimate(marginal,correlation_matrix=correlationMatrix)
}
