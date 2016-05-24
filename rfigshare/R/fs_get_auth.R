# Helper function to get authentication
#
# The authentication environment is created by new.env function in the zzz.R file.  The authentication token 'oauth' is created by the figshare_auth function.  This helper function lets all other functions load the authentication.  
# @keywords internal
 fs_get_auth <- function(){
   if(!exists("oauth", envir=FigshareAuthCache))
     tryCatch(fs_auth(), error= function(e) 
       stop("Requires authentication. 
       Are your credentials stored in options? 
       See fs_auth function for details."))
   get("oauth", envir=FigshareAuthCache)
 }
# 
# 
# 
