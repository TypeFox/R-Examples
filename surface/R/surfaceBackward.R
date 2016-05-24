surfaceBackward<-
function(otree,odata,starting_model,aic_threshold = 0, max_steps =  NULL,save_steps = FALSE,filename = "temp_back_list.R",verbose = FALSE, only_best = FALSE, plotaic = FALSE, error_skip = FALSE, sample_shifts = FALSE, sample_threshold = 2){

if(any(duplicated(as(otree,"data.frame")$labels)))
	stop("Each node in 'tree' must have a unique label for back-compatibility between formats. Use 'nameNodes(tree)' prior to converting to ouch format")

n<-otree@nterm;nt<-dim(odata)[2]
if(is.null(max_steps))max_steps<-starting_model$n_regimes[1]

back_list<-list()
back_list[[1]]<-list(fit=starting_model$fit,delta_aic=NULL, aic=starting_model$aic,savedshifts=starting_model$savedshifts,n_regimes=starting_model$n_regimes)
for(j in 2:max_steps){
	back_list[[j]]<-collapseRegimes(otree,odata,oldshifts=back_list[[j-1]]$savedshifts,oldaic=back_list[[j-1]]$aic,oldfit=back_list[[j-1]]$fit, aic_threshold=aic_threshold,verbose=verbose,only_best=only_best,plotaic=plotaic, error_skip=error_skip,sample_shifts=sample_shifts, sample_threshold=sample_threshold)
	if(save_steps)save(back_list,file=filename)
	if(back_list[[j]]$aic==back_list[[j-1]]$aic)break
	}
if(back_list[[j]]$aic==back_list[[j-1]]$aic){
	back_list<-back_list[1:(j-1)]
	}
back_list
}
