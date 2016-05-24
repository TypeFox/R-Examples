surfaceForward<-
function(otree, odata, starting_list=NULL, starting_shifts=NULL, exclude=0, aic_threshold=0, max_steps=NULL, save_steps=FALSE, filename="temp_out_list.R", verbose=FALSE, plotaic=FALSE, error_skip=FALSE, sample_shifts=FALSE, sample_threshold=2){

if(any(duplicated(as(otree,"data.frame")$labels)))
	stop("Each node in 'tree' must have a unique label for back-compatibility between formats. Use 'nameNodes(tree)' prior to converting to ouch format")

n<-otree@nterm;nt<-dim(odata)[2]
if(is.null(max_steps))max_steps<-otree@nnodes

if(is.null(starting_list)){
	out_list<-startingModel(otree,odata,starting_shifts)
}else{
	out_list<-starting_list
	}
for(j in (length(out_list)+1):max_steps){
	out_list[[j]]<-addRegime(otree,odata,oldshifts=out_list[[j-1]]$savedshifts,oldfit=out_list[[j-1]]$fit,oldaic=out_list[[j-1]]$aic, alloldaic=out_list[[j-1]]$all_aic,exclude=exclude, aic_threshold=aic_threshold, verbose=verbose,plotaic=plotaic,error_skip=error_skip,sample_shifts=sample_shifts, sample_threshold=sample_threshold)
	if(save_steps)save(out_list,file=filename)
	if((out_list[[j]]$aic-out_list[[j-1]]$aic)>=(aic_threshold))break
	}
if((out_list[[j]]$aic-out_list[[j-1]]$aic)>=(aic_threshold)){
	out_list<-out_list[1:(j-1)]
	}
out_list	
}
