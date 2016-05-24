#' Get Login credentials for NBN access
#' 
#' Opens a dialog box for the user to enter credentials if
#' username and password are not given, also manages cookies.
#' 
#' @export
#' @import httr
#' @param username optional character giving username
#' @param password optional character giving password
#' @param verbose logical, if \code{TRUE} successful logins are reported via console
#' @details This function is used within the getOccurrences function and should not
#' normally need to be used by users. It automatically handles cookies and only
#' prompts users for username and password if no cookies exist. The facility to provide
#' a username and password has been added for people who may have issues with cookies
#' or perhaps share an account with a number of others, however hardcoding usernames and
#' passwords has obvious security issues.
#' @return NULL
#' @author Tom August, CEH \email{tom.august@@ceh.ac.uk}


nbnLogin <- function(username = NULL, password = NULL, verbose = FALSE){
    
    # If we are not providing a username and password see if we have cookies
    if(is.null(username) | is.null(password)){    
        
        # See if we are known
        whoamI <- "https://data.nbn.org.uk/api/user"
        
        resp_who <- GET(whoamI)
       
        # check ststus
        if(!grepl("^2[[:digit:]]{2}$", resp_who$status_code)) stop(paste('Error accessing https://data.nbn.org.uk/api/user -', http_status(resp_who)$message))
                   
        #Login in with dialog box
        if(content(resp_who)$id == 1){
            
            # Get username and password
            dusername <- readline("Enter Username:")
            dpassword <- readline("Enter Password:")
            
            # Create login URL
            urlLogin <- paste("https://data.nbn.org.uk/api/user/login?username=", gsub(' ', '%20', dusername),
                              "&password=", gsub(' ', '%20', dpassword), sep='')
            
            # Check that login was a success (if not stop)
            resp_login <- GET(urlLogin)
            
            if(!grepl("^2[[:digit:]]{2}$", resp_login$status_code)) stop(paste('Error accessing https://data.nbn.org.uk/api/user/login? -', http_status(resp_login)$message))
                        
            if(!content(resp_login)$success){
                stop('Username and password invalid, have another go')
            } else {
                if(verbose) print('Login successful')
            }
        } else {
            if(verbose) print(paste('Logged in as', content(resp_who)$username, 'using cookies'))
        }
    } else { # If we have specified a username and password
        
        # Create login URL
        urlLogin <- paste("https://data.nbn.org.uk/api/user/login?username=", gsub(' ', '%20', username),
                          "&password=", gsub(' ', '%20', password), sep='')
        
        # Check that login was a success (if not stop)
        resp_login <- GET(urlLogin)
        
        # If 'Unauthorised'
        if(resp_login$status_code == 401) stop('Username and password invalid, have another go')
        
        # If some other error
        if(!grepl("^2[[:digit:]]{2}$", resp_login$status_code)) stop(paste('Error accessing https://data.nbn.org.uk/api/user/login? -', http_status(resp_login)$message))
        
        if(!content(resp_login)$success){
            stop('Username and password invalid, have another go')
        } else {
            if(verbose) print('Login successful')
        }        
    }   
}