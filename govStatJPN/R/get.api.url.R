get.api.url <-
function(app.id, stat.type = 'stat-list', param, stat.id='', n=''){
	file.uri  <-  'http://proory.com/govStatJPN.php'
	if(stat.type == 'stat-list'){
		api.param <- rep("", (length(param) + 2))
		api.param[1] <- paste('t', stat.type, sep="=")
		api.param[2] <- paste('appid', app.id, sep='=')
		for(i in 1:length(param))
		{
			api.param[i+2] <- paste(names(param)[i], param[i], sep="=")
		}
	
	}else if(stat.type == 'meta-info'){
		api.param <- c(
			paste('t', stat.type, sep="="),
			paste('appid', app.id, sep='='),
			paste('statsDataId', stat.id, sep='=')
		)
	}else{
		api.param <- rep("", (length(param) + 3))
		api.param[1] <- paste('t', stat.type, sep="=")
		api.param[2] <- paste('appid', app.id, sep='=')
		api.param[3] <- paste('statsDataId', stat.id, sep='=')
		for(i in 1:length(param))
		{
			api.param[i+3] <- paste(names(param)[i], param[i], sep="=")
		}
	}
	return(paste(file.uri, paste(api.param, collapse="&"), sep="?"))
}

