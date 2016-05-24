#' List all models deployed.
#'
#' @param auth_details authorization parameters (must be stored in list format: list(user_name='exampleusername',api_key='exampleapikey')).
#'
#' @return A data.frame object with columns model_name, model_size, last_update, deploy_type.
#' @seealso \code{\link{duke_deploy}} to deploy a model or \code{\link{duke_predict}} for making predictions with deployed models.
#' @export
#' @examples
#' 
#' duke_auth <- list(
#' user_name = "try_it",
#' api_key   = "db1542b66f16aba5768d8a19c27dec4facf9168a",
#' endpoint = "/api/v1.0"
#' )
#' duke_list(auth_details=duke_auth)

duke_list <- function(auth_details){
  if (!("user_name" %in% names(auth_details))) return("No user_name supplied as part of the auth_details parameter.")
  if (!("api_key" %in% names(auth_details))) return("No api_key supplied as part of the auth_details parameter.")
  if (!("endpoint" %in% names(auth_details))) return("No endpoint supplied as part of the auth_details parameter.")
  r <- httr::GET(url= paste("http://deploy.dukeanalytics.com",auth_details$endpoint,"/list/",auth_details$user_name,"/",auth_details$api_key,sep=""))
  result <- httr::content(r,"parsed","application/json")
  result_df <- data.frame(matrix(unlist(result), ncol=4, byrow=T));result_df<- result_df[,c(2,3,1,4)];names(result_df) <- c("model_name","model_size","last_update","deploy_type")
  return(result_df)
}
