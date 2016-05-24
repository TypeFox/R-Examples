readData <- function( Y.original, ratio=NULL, range=NULL, balance=NULL, eps.bal=0.6 ){
	
	if (is.null(Y.original)){
		stop("Y.original is NULL")
	}
	D_EditIn_checked = .Check_D_EditIn(Y.original)
		
	var_names = dimnames(D_EditIn_checked)[[2]] ; n_var = dim(D_EditIn_checked)[[2]]
	
	RatioEdit_S1 = NULL
	RangeRest_S1 = NULL
	BalanceEdit_Final = NULL
	n_balance = 0
	
	if (!is.null(ratio)){
		
		if (class(ratio)=="data.frame"){
			RatioEdit_S1 = .read_Ratio_matrix(ratio,var_names,n_var)
		} else {
			if (class(ratio)=="editmatrix"){
				RatioEdit_S1 = .read_Ratio_editmatrix(ratio,var_names,n_var)
			} else {
			  stop("class(ratio) needs to be 'data.frame' or 'editmatrix'")		
			} # if (class(ratio)=="editmatrix") ... else ...
		} # if (class(ratio)=="matrix") ... else ...
		
	} # if (!is.null(ratio)) 
	
	if (!is.null(range)){

		if (class(range)=="data.frame"){
			RangeRest_S1 = .read_Range_matrix(range,var_names,n_var)
		} else {
			if (class(range)=="editmatrix"){
				RangeRest_S1 = .read_Range_editmatrix(range,var_names,n_var)
			} else {
			  stop("class(range) needs to be 'data.frame' or 'editmatrix'")		
			} # if (class(RangeRestriction)=="editmatrix") ... else ...
		} # if (class(RangeRestriction)=="matrix") ... else ...

	} # if (is.null(RangeRestriction)) 
	
	if (!is.null(balance)){

		BalanceEdit_Temp = NULL
		
		if (class(balance)=="data.frame"){
			if ( dim(balance)[[2]] != n_var ){
				stop("The number of columns of 'balance' does not match with that of 'Y.original'.")
			} else {
				BalanceEdit_Temp = .read_Balance_matrix(balance,n_var)
			}
			
		} else {
			if (class(balance)=="editmatrix"){
				BalanceEdit_Temp = .read_Balance_editmatrix(balance,var_names,n_var)
			} else {
			  stop("class(range) needs to be 'data.frame' or 'editmatrix'")		
			} # if (class(balance)=="editmatrix") ... else ...
		} # if (class(balance)=="matrix") ... else ...
		
		n_balance = dim(BalanceEdit_Temp)[[1]]
		
		if (!is.null(BalanceEdit_Temp)){
				BalanceEdit_Final = .make_BE_final(BalanceEdit_Temp,n_var,eps.bal)
		}
		
	} # if (is.null(balance)) 
	
	if ( is.null(ratio)*is.null(range)*is.null(balance)==1 ){
		Edit_editmatrix = NULL
	} else {
		Edit_editmatrix = .make_editmatrix(RatioEdit_S1,RangeRest_S1,BalanceEdit_Final,n_var)
	}
	
	###########
	
	Bound_L_U = .make_Bound_LU(D_EditIn_checked,RangeRest_S1,var_names,n_var)
	
	if ( is.null(ratio)*is.null(range)*is.null(balance)==1 ){
		Edit_matrix = .make_matrix(NULL,Bound_L_U,NULL,var_names,n_var)
	} else {
		Edit_matrix = .make_matrix(RatioEdit_S1,Bound_L_U,BalanceEdit_Final,var_names,n_var)
	}
	
	###########
	
	FaultyRecordID = NULL
	for (i_sample in 1:dim(D_EditIn_checked)[[1]]){
		if (max( Edit_matrix[,1:n_var]%*%t(D_EditIn_checked[i_sample,])-Edit_matrix[,"CONSTANT"] )>0) FaultyRecordID = c(FaultyRecordID,i_sample);
	}
	
	if (is.null(FaultyRecordID)) stop("There is no record violating edit rules. Check if the edit rules are correctly specified")
	
	EditIn_format = list(Y.input=D_EditIn_checked,Edit.editmatrix=Edit_editmatrix,Edit.matrix=Edit_matrix,Bound.LU=Bound_L_U,ratio=RatioEdit_S1,n.balance=n_balance,FaultyRecordID=FaultyRecordID)
	class(EditIn_format) <- "EditIn.data"
	return(EditIn_format)	
	
} # readData <- function


###############################################################################################
###############################################################################################

.read_Ratio_matrix <- function(RatioEdit.import,var.names,n.var){
	
	dimnames(RatioEdit.import)[[2]] = c("A","B","UpperLimit(A/B)")
	n.RatioEdit = dim(RatioEdit.import)[[1]] 

	RatioEdit.Mat1_return = array(0,c(n.RatioEdit,(n.var+1)))
	dimnames(RatioEdit.Mat1_return)[[2]] = c(var.names,"CONSTANT")
	for (i_row in 1:n.RatioEdit){
		RatioEdit.Mat1_return[i_row,which(RatioEdit.import[i_row,1]==var.names)] = 1
		RatioEdit.Mat1_return[i_row,which(RatioEdit.import[i_row,2]==var.names)] = -1 * RatioEdit.import[i_row,3]
	}
	
	return( RatioEdit.Mat1_return )
	
} # read.Ratio.matrix

.read_Ratio_editmatrix <- function(RatioEdit.expression,var.names,n.var){
	
	Temp_RE.exp.mat = as.matrix(RatioEdit.expression)
	names_Temp_RE = dimnames(Temp_RE.exp.mat)[[2]]
	n.RatioEdit = dim(Temp_RE.exp.mat)[[1]]
	
	RatioEdit.Mat2_return = array(0,c(n.RatioEdit,(n.var+1)))
	dimnames(RatioEdit.Mat2_return)[[2]] = c(var.names,"CONSTANT")
	for (i_temp in 1:length(names_Temp_RE)){
		RatioEdit.Mat2_return[,names_Temp_RE[i_temp]] = Temp_RE.exp.mat[,i_temp]
	}
	
	return( RatioEdit.Mat2_return )
	
} # read.Ratio.editmatrix


.read_Range_matrix <- function(RangeRest.import,var.names,n.var){
	
	dimnames(RangeRest.import)[[2]] = c("VariableName","LowerLimit","UpperLimit")
	n.RangeRest_temp = 2*dim(RangeRest.import)[[1]] 

	RangeRest.Mat1_return = array(0,c(n.RangeRest_temp,(n.var+1)))
	dimnames(RangeRest.Mat1_return)[[2]] = c(var.names,"CONSTANT")
	for (i_row in 1:(n.RangeRest_temp/2)){
		RangeRest.Mat1_return[(2*i_row-1),which(RangeRest.import[i_row,1]==var.names)] = -1 
		RangeRest.Mat1_return[(2*i_row-1),"CONSTANT"] = - RangeRest.import[i_row,2]
		RangeRest.Mat1_return[(2*i_row),which(RangeRest.import[i_row,1]==var.names)] = 1 
		RangeRest.Mat1_return[(2*i_row),"CONSTANT"] = RangeRest.import[i_row,3]
	}

	isInf = which(abs(RangeRest.Mat1_return[,"CONSTANT"])==Inf)
	if ( length(isInf)>0 ) RangeRest.Mat1_return = RangeRest.Mat1_return[-isInf,] ;
	
	return( RangeRest.Mat1_return ) ;  
	
} # read.Range.matrix <- function

.read_Range_editmatrix <- function(RangeRest.expression,var.names,n.var){
	
	Temp_RE.exp.mat = as.matrix(RangeRest.expression)

	names_Temp_RE = dimnames(Temp_RE.exp.mat)[[2]]
	n.RangeRest = length(names_Temp_RE)
	RangeRest.Mat2_return = array(0,c(n.RangeRest,(n.var+1)))
	dimnames(RangeRest.Mat2_return)[[2]] = c(var.names,"CONSTANT")
	for (i_temp in 1:length(names_Temp_RE)){
		RangeRest.Mat2_return[,names_Temp_RE[i_temp]] = Temp_RE.exp.mat[,i_temp]
	}

	isInf = which(abs(RangeRest.Mat2_return[,"CONSTANT"])==Inf)
	if ( length(isInf)>0 ) RangeRest.Mat2_return = RangeRest.Mat2_return[-isInf,] ;
	
	return( RangeRest.Mat2_return ) ;
	
} # read.Range.editmatrix <- function

.read_Balance_matrix <- function(Balance.Edit,n.var){

	CONSTANT = rep(0,dim(Balance.Edit)[[1]])
	BalanceEdit.Mat1_return = cbind(Balance.Edit,CONSTANT)
	
	return( as.matrix(BalanceEdit.Mat1_return) ) ;
	
} # read.Balance.matrix <- function


.read_Balance_editmatrix <- function(Balance.Edit,var.names,n.var){
	
	Temp_BE.exp.mat = as.matrix(Balance.Edit)

	names_Temp_BE = dimnames(Temp_BE.exp.mat)[[2]]
	
	BalanceEdit.Mat2_return = array(0,c(dim(Temp_BE.exp.mat)[[1]],(n.var+1)))
	dimnames(BalanceEdit.Mat2_return)[[2]] = c(var.names,"CONSTANT")
	for (i_temp in 1:length(names_Temp_BE)){
		BalanceEdit.Mat2_return[,names_Temp_BE[i_temp]] = Temp_BE.exp.mat[,i_temp]
	}
	
	return( BalanceEdit.Mat2_return ) ;
	
} # read.Balance.editmatrix <- function

.make_BE_final <- function(Balance.Edit.Temp,n.var,Eps){
	
	Balance.Edit.Final_return = rbind(Balance.Edit.Temp,Balance.Edit.Temp)
	for (i_row in 1:dim(Balance.Edit.Temp)[[1]]){
		Balance.Edit.Final_return[(2*i_row-1),] = Balance.Edit.Temp[i_row,]
		Balance.Edit.Final_return[(2*i_row),] = -Balance.Edit.Temp[i_row,]
	}
	Balance.Edit.Final_return[,(n.var+1)] = Eps
		
	return( Balance.Edit.Final_return )	
	
} # make.BE.final <- function

.make_editmatrix <- function(RatioEdit.Mat1,RangeRest.Mat1,Balance.Edit.Final,n.var){
	
	row_name_temp = NULL
	
	if (!is.null(RatioEdit.Mat1)){
		for (i in 1:dim(RatioEdit.Mat1)[[1]]) row_name_temp = c(row_name_temp,paste("Rt",i,sep="")) ;
	}
	if (!is.null(RangeRest.Mat1)){
		for (i in 1:dim(RangeRest.Mat1)[[1]]) row_name_temp = c(row_name_temp,paste("Rn",i,sep="")) ;
	}
	if (!is.null(Balance.Edit.Final)){
		for (i in 1:(dim(Balance.Edit.Final)[[1]]/2)){
			 row_name_temp = c(row_name_temp,paste("B",i,"a",sep=""),paste("B",i,"b",sep=""))
		};
	}
	
	EditCombined_editrules_ = rbind(RatioEdit.Mat1,RangeRest.Mat1,Balance.Edit.Final)
	dimnames(EditCombined_editrules_)[[1]] = row_name_temp

	return( as.editmatrix(as.matrix(EditCombined_editrules_[,c(1:n.var)]), b=EditCombined_editrules_[,(n.var+1)], ops = c(rep("<=",dim(EditCombined_editrules_)[[1]])) ) )
	
} # make.editmatrix <- function

.make_Bound_LU <- function(Y.in,RangeRest.Mat1,var.names,n.var){
	
	n_RangeE_temp = dim(RangeRest.Mat1)[[1]] 

	Bound_L_U.return = array(NA,c(n.var,2))
	dimnames(Bound_L_U.return)[[1]] = var.names ; dimnames(Bound_L_U.return)[[2]] = c("Lower","Upper")

	if (!is.null(RangeRest.Mat1)){
		
		for (i_row in 1:n_RangeE_temp){	
			i_var = which(RangeRest.Mat1[i_row,(1:n.var)] != 0)
				if (RangeRest.Mat1[i_row,i_var]==1){
					Bound_L_U.return[var.names[i_var],"Upper"] = RangeRest.Mat1[i_row,"CONSTANT"]
				} 
				if (RangeRest.Mat1[i_row,i_var]==-1){
					Bound_L_U.return[var.names[i_var],"Lower"] = -RangeRest.Mat1[i_row,"CONSTANT"]
				}	
		} ## for (i_row in 1:n_RangeE_temp)
		
	} # if (!is.null(RangeRest.Mat1))

	for (i_var in 1:n.var){	
		if ( is.na(Bound_L_U.return[var.names[i_var],"Upper"]) ){
			Bound_L_U.return[var.names[i_var],"Upper"] = 10 * max(Y.in[,i_var])
		}
		if ( is.na(Bound_L_U.return[var.names[i_var],"Lower"]) ){
			Bound_L_U.return[var.names[i_var],"Lower"] = max( ( 0.1 * min(Y.in[,i_var]) ), 1e-5 ) 
		}
	} ## for (i_var in 1:n_var)
	
	return(Bound_L_U.return)
	
} # make.Bound.LU <- function

.make_matrix <- function(RatioEdit.Mat1,Bound.L.U,Balance.Edit.Final,var.names,n.var){
	
	n.RangeFinal = 2*n.var
	RangeFinal = array(0,c(n.RangeFinal,(n.var+1)))
	dimnames(RangeFinal)[[2]] = c(var.names,"CONSTANT")
	for (i_row in 1:(n.RangeFinal/2)){
		RangeFinal[(2*i_row-1),var.names[i_row]] = -1 
		RangeFinal[(2*i_row-1),"CONSTANT"] = - Bound.L.U[i_row,1]
		RangeFinal[(2*i_row),var.names[i_row]] = 1 
		RangeFinal[(2*i_row),"CONSTANT"] = Bound.L.U[i_row,2]
	}

	return(as.matrix(rbind(RatioEdit.Mat1,RangeFinal,Balance.Edit.Final)))

} # make.matrix <- function

