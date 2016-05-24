#' Predict using deployed model using REST API.
#'
#' @param auth_details authorization parameters (must be stored in list format: list(user_name='exampleusername',api_key='exampleapikey'))
#' @param model_name the name of the deployed model for predictions.
#' @param new_data new feature values for predictions.
#' @param host_name optional parameter.
#'
#' @return On successful prediction, a character vector with model predictions. An error message is returned on failure.
#' @seealso \code{\link{duke_deploy}} for deploying a model or \code{\link{duke_list}} to display a list of all deployed models.
#' @export
#' @examples
#' 
#' duke_auth <- list(
#' user_name = "try_it",
#' api_key   = "db1542b66f16aba5768d8a19c27dec4facf9168a",
#' endpoint = "/api/v1.0"
#' )
#' example_data <- as.data.frame(cbind(gl(3,50),rnorm(150)));names(example_data) <- c("x","y")
#' duke_predict(auth_details=duke_auth,model_name="example_fit",example_data)

duke_predict <- function(auth_details,model_name,new_data,host_name=NULL){
  if (!("user_name" %in% names(auth_details))) return("No user_name supplied as part of the auth_details parameter.")
  if (!("api_key" %in% names(auth_details))) return("No api_key supplied as part of the auth_details parameter.")
  if (!("endpoint" %in% names(auth_details))) return("No endpoint supplied as part of the auth_details parameter.")
  if (!is.null(host_name)){
    if(grepl("*.dukeanalytics.com",host_name)>0) return("Wrong host address. Model not deployed. Please try with the correct host_name parameter in the auth_details list.")
  }else{host_name <- "deploy.dukeanalytics.com"}
  if(is.null(model_name)) return("Please specify an appropriate model_name parameter.")
  if(grepl("[^[:alnum:]_]", model_name)) return("Model names may only contain alphanumeric and underscore characters.")
  
  body <- list("Model"=model_name,
               "New_Data"=jsonlite::toJSON(new_data))  
  r <- httr::POST(paste("http://",host_name,auth_details$endpoint,"/predict/",auth_details$user_name,
                        "/",auth_details$api_key,sep=""), 
            body = body,
            encode="form")
  result <- httr::content(r, "parsed", "application/json")
  #return(jsonlite::fromJSON(result))
  return(as.character(unlist(result)))
}
