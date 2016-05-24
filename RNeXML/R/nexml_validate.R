ONLINE_VALIDATOR <- "http://162.13.187.155/nexml/phylows/validator"
CANONICAL_SCHEMA <- "http://162.13.187.155/nexml/xsd/nexml.xsd"
#ONLINE_VALIDATOR <- "http://www.nexml.org/nexml/phylows/validator"
#CANONICAL_SCHEMA <- "http://www.nexml.org/2009/nexml.xsd"

#' validate nexml using the online validator tool
#' @param file path to the nexml file to validate
#' @param schema URL of schema (for fallback method only, set by default).  
#' @details Requires an internet connection.  see http://www.nexml.org/nexml/phylows/validator for more information in debugging invalid files
#' @return TRUE if the file is valid, FALSE or error message otherwise
#' @export
#' @import httr XML
#' @examples \dontrun{
#' data(bird.orders)
#' birds <- nexml_write(bird.orders, "birds_orders.xml")
#' nexml_validate("bird_orders.xml")
#' unlink("bird_orders.xml") # delete file to clean up
#' }
nexml_validate <- function(file, schema=CANONICAL_SCHEMA){
  a = POST(ONLINE_VALIDATOR, body=list(file = upload_file(file)))
  if(a$status_code %in% c(200,201)){
    TRUE
  } else if(a$status_code == 504){
    warning("Online validator timed out, trying schema-only validation.")
    nexml_schema_validate(file, schema=schema)

  } else if(a$status_code == 400){
    warning(paste("Validation failed, error messages:",
         xpathSApply(htmlParse(content(a, "text")), 
                     "//li[contains(@class, 'error') or contains(@class, 'fatal')]", xmlValue)
         ))
    FALSE
  } else {
    warning(paste("Unable to reach validator. status code:", a$status_code, ".  Message:\n\n", content(a, "text")))
    NULL
  }
}




nexml_schema_validate <- function(file, schema=CANONICAL_SCHEMA){
  a = GET(schema)
  if(a$status_code == 200){
    if(is.null(xmlSchemaParse(schema))){
        warning(paste("Schema not accessible at", schema))
        NULL
    } else {
      result <- xmlSchemaValidate(schema, file) 
      if(length(result$errors) == 0){
        TRUE
      } else {
        warning(paste(result$errors))
        FALSE
      }
    }
  } else {
    warning("Unable to obtain schema, couldn't validate")
    NULL
  }
    
}
#xmlSchemaValidate(xmlSchemaParse(content(a, "text"), asText=TRUE), file)   # fails to get other remote resources



