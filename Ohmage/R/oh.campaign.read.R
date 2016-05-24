#' List campaigns available on the server
#' @param output_format short or long
#' @param to.data.frame should response be converted to a dataframe
#' @param ... other arguments passed to ohmage
#' @export
oh.campaign.read <- function(output_format="short", to.data.frame=TRUE, ...){
	xhr <- oh.call("/campaign/read", output_format=output_format, ...);
	if(output_format=="xml"){
		xhr$data$xmlTree <- XML::xmlTreeParse(xhr$data$configuration, useInternalNodes=T);
	}

	if(!to.data.frame){
		return(xhr);
	}

	if(length(xhr$data) == 0){
		message("No active campaigns founds.")
	} else {

		if(output_format=="short"){
			#unlist some stuff
			for(i in 1:length(xhr$data)){
				xhr$data[[i]] <- unlist(lapply(xhr$data[[i]][c("description","name","privacy_state","creation_timestamp","running_state")], unlist)); #note that the [-1] excludes the user role list.
			}

			mydataframe <- as.data.frame(do.call(rbind, xhr$data));
			mydataframe$UUID <- row.names(mydataframe);
			#row.names(mydataframe) <- NULL;
			return(mydataframe);
		}

		if(output_format=="long"){
			#unlist some stuff
			for(i in 1:length(xhr$data)){
				xhr$data[[i]]$classes <- unlist(xhr$data[[i]]$classes);
				xhr$data[[i]]$user_role_campaign <- lapply(xhr$data[[i]]$user_role_campaign, unlist)
				xhr$data[[i]]$user_roles <- unlist(xhr$data[[i]]$user_roles);
			}

			xhr$metadata$items <- unlist(xhr$metadata$items);
		}

		return(xhr)
	}
}

