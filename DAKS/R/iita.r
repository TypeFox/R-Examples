###########################################
# Inductive item tree analysis algorithms #
###########################################

#############################################
#                                           #
# This function can be used to perform one  #
# of the three inductive item tree analysis #
# algorithms (original, corrected, and      # 
# minimized corrected) selectively.         #
#                                           # 
############################################# 

iita<-function(dataset, v){
if ((!is.data.frame(dataset) & !is.matrix(dataset)) || ncol(dataset) == 1){
stop("data must be either a numeric matrix or a data.frame, with at least two columns.\n")
}

if(sum(!(dataset == 0 | dataset == 1) | is.na(dataset)) != 0){
stop("data must contain only 0 and 1")
}

if(v != 1 && v != 2 && v !=3){
stop("IITA version must be specified")
}

# call to the chosen algorithm
if(v == 3){
i<-ind_gen(ob_counter(dataset))
ii<-orig_iita(dataset, i)
}

if(v == 2){
i<-ind_gen(ob_counter(dataset))
ii<-corr_iita(dataset, i)
}

if(v == 1){
i<-ind_gen(ob_counter(dataset))
ii<-mini_iita(dataset, i)
}

obj<-list(diff = ii[[1]], implications = i[which.min(ii[[1]])][[1]], error.rate = ii[[2]][which.min(ii[[1]])], selection.set.index = which.min(ii[[1]])[[1]], v=v)
class(obj)<-"iita"
return(obj)
}
