#' Deploy model using REST API
#'
#' @param auth_details authorization parameters (must be stored in list format: list(user_name='exampleusername',api_key='exampleapikey')).
#' @param model_object the model to deploy e.g. lm() or glm() object.
#' @param predict_udf user defined function as alternative to the built-in model object predict function.
#' @param host_name optional parameter.
#'
#' @return Deployment status, success or failure with message.
#' @seealso \code{\link{duke_predict}} for making predictions with deployed models or \code{\link{duke_list}} to display a list of all deployed models.
#' @export
#' @examples
#' 
#' duke_auth <- list(
#' user_name = "try_it",
#' api_key   = "db1542b66f16aba5768d8a19c27dec4facf9168a",
#' endpoint = "/api/v1.0"
#' )
#' example_data <- as.data.frame(cbind(gl(3,50),rnorm(150)));names(example_data) <- c("x","y")
#' example_fit  <- lm(y~x,data=example_data)
#' duke_deploy(auth_details=duke_auth,model_object=example_fit)

duke_deploy <- function(auth_details,model_object,predict_udf=NULL,host_name=NULL){
  if (!("user_name" %in% names(auth_details))) return("No user_name supplied as part of the auth_details parameter.")
  if (!("api_key" %in% names(auth_details))) return("No api_key supplied as part of the auth_details parameter.")
  if (!("endpoint" %in% names(auth_details))) return("No endpoint supplied as part of the auth_details parameter.")
  if (!is.null(host_name)){
    if(grepl("*.dukeanalytics.com",host_name)>0) return("Wrong host address. Model not deployed. Please try with the correct host_name parameter in the auth_details list.")
  }else{host_name <- "deploy.dukeanalytics.com"}
  model_name <- deparse(substitute(model_object))
  if(is.null(model_name)) return("Please specify an appropriate model_name parameter.")
  if(grepl("[^[:alnum:]_]", model_name)) return("Model names may only contain alphanumeric and underscore characters.")
  if(!is.null(predict_udf)) {
      model_object$duke_predict_udf <- predict_udf
    }#else {return("Invalid predict function. Please retry.")}
#function(model_object,new_data) {out <- predict(model_object,newdata=new_data); return(out);}  
  temp_file <- base::tempfile(pattern=paste(model_name,"-",sep=""),fileext=".rds") 
  saveRDS(model_object,file=temp_file)
  model_object_size <- format(utils::object.size(model_object),units="Mb")
  model_upload_size <- round(base::file.info(temp_file)$size/1000000,2)
  if(model_upload_size>20)
  {
    return(paste("Model is ",model_upload_size," Mb in size, the limit is 20Mb.",sep=""))
  }

  #print(paste("Model upload size is :", model_upload_size, "Mb, Model size is :",model_object_size,sep=""))
  body <- list("files"=httr::upload_file(temp_file))  
  
  r        <- httr::POST(paste("http://",host_name,auth_details$endpoint,"/deploy","/",auth_details$user_name,
                               "/",auth_details$api_key,sep=""), 
                         body = body,
                         encode="multipart")
  
  result <- httr::content(r, "parsed", "application/json")
  base::unlink(temp_file)
  return(result)
}
