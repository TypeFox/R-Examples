zendesk <- function(username, password, url){
    if(!is.null(username) & !is.null(password) & !is.null(url)){
        .ZendeskEnv$data$username <- username
        .ZendeskEnv$data$password <- password
        .ZendeskEnv$data$url <- gsub("\\/$","", url)
    }
    else{
        warning("Username, Password and URL must be provided in order to access your organization's data")
    }
}

