listToJSONStr <- function(list, direct = FALSE){

	# CONVERTS PARAMETERS INTO STRING FOR JAVASCRIPT TO READ

	json_str <- '{'

	if(direct){

		#out_str <- '{"update_status":"8 landmarks saved. 2 curves saved (control and curve points)."}'

		for(param_name in names(list)){

			if(!is.list(list[[param_name]])){
				json_str <- paste0(json_str, paste0('"', param_name, '":["', paste0(list[[param_name]], collapse='", "'), '"]', sep=', '))
			}else{
				if(length(list[[param_name]]) == 0) next
				json_str <- paste0(json_str, paste0('"', param_name, '":['))
				for(i in 1:length(list[[param_name]])){
					json_str <- paste0(json_str, paste0('["', paste0(list[[param_name]][[i]], collapse='", "'), '"]', sep=', '))
				}			
				json_str <- paste0(json_str, '],')
			}
		}
	}else{

		for(param_name in names(list)){

			if(!is.list(list[[param_name]])){
				json_str <- paste0(json_str, paste0('\"', param_name, '\":[\"', paste0(list[[param_name]], collapse='\", \"'), '\"]', sep=', '))
			}else{
				if(length(list[[param_name]]) == 0) next
				json_str <- paste0(json_str, paste0('\"', param_name, '\":['))
				for(i in 1:length(list[[param_name]])){
					json_str <- paste0(json_str, paste0('[\"', paste0(list[[param_name]][[i]], collapse='\", \"'), '\"]', sep=', '))
				}			
				json_str <- paste0(json_str, '],')
			}
		}
	}

	json_str <- paste0(json_str, "}")
	json_str <- gsub("(,[ ]*)(]|})", "\\2", json_str)

	json_str
}