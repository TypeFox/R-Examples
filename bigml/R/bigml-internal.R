#' @import RJSONIO plyr RCurl

.BIGML_URL <-
"https://bigml.io/andromeda/"
.DATASET_URL <-
"https://bigml.io/andromeda/dataset"
.MODEL_URL <-
"https://bigml.io/andromeda/model"
.PREDICTION_URL <-
"https://bigml.io/andromeda/prediction"
.SOURCE_URL <-
"https://bigml.io/andromeda/source"

.fixid <- 
function(name)
{
	stripthis = paste("^\\w+/",sep='')
	gsub(stripthis, '',name)
}

.basic_api <-
function (resource) 
{
    list(get = function(id, ...) {
		id = .fixid(id)
        resource = paste(resource, "/", id, sep = "")
        result = fromJSON(getURL(.build_url(resource, ...)))
        .check_for_code(result)
        result
    }, list = function(...) {
        result = fromJSON(getURL(.build_url(resource, ...)))
        .check_for_code(result)
        result
    }, postJson = function(postfields, ...) {
		if ('debug' %in% names(list(...))){
			writeLines(toJSON(postfields))
		}
        result = fromJSON(getURL(.build_url(resource, ...), customrequest = "POST", 
            httpheader = c(`Content-Type` = "application/json"), 
            postfields = toJSON(postfields)))
        .check_for_code(result)
        result
    }, putJson = function(id, postfields, ...) {
		if ('debug' %in% names(list(...))){
			writeLines(toJSON(postfields))
		}
        resource = paste(resource, "/", id, sep = "")
        result = fromJSON(getURL(.build_url(resource, ...), customrequest = "POST", 
            httpheader = c(`Content-Type` = "application/json"), 
            postfields = toJSON(postfields)))
        .check_for_code(result)
        result
    }, upload = function(file, name = NULL, source_parser = NULL, ...) {
        opts = list(file = file)
        if (!is.null(source_parser)) opts$source_parser = source_parser
        if (!is.null(name)) opts$name = name
        result = postForm(.build_url(resource, ...), .params = opts, 
            .checkParams = FALSE)
        result = fromJSON(result)
        .check_for_code(result)
        result
    }, delete = function(resource_id, ...) {
		id = .fixid(resource_id)
		resource = paste(resource, "/", id, sep = "")
        getURL(.build_url(resource, ...), customrequest = "DELETE")
    }, update = function(resource_id, ...) {
		id = .fixid(resource_id)
		resource = paste(resource, "/", id, sep = "")
	    getURL(.build_url(resource, ...), customrequest = "UPDATE")
	}

)

}
.build_url <-
function (request, ...) 
{
	val = paste(request, formEncodeURL(...), sep = "")
	if ('debug' %in% names(list(...))){
		message(val)
	}
    val
}
.check_for_code <-
function (result) 
{
    if ("code" %in% names(result) && result$code %in% .error_codes) {
        code = names(.error_codes[which(.error_codes == result$code)])
    }
    else if ("status" %in% names(result) && "code" %in% result$status && 
        result$status$code %in% .error_codes) {
        code = names(.error_codes[which(.error_codes == result$status$code)])
    }
    else {
        return(TRUE)
    }
    msg = paste("Error: BigML returned code", code, toJSON(result))
    stop(msg)
}
.error_codes <-
structure(list(HTTP_NO_CONTENT = 204, HTTP_BAD_REQUEST = 400, 
    HTTP_UNAUTHORIZED = 401, HTTP_PAYMENT_REQUIRED = 402, HTTP_FORBIDDEN = 403, 
    HTTP_NOT_FOUND = 404, HTTP_METHOD_NOT_ALLOWED = 405, HTTP_LENGTH_REQUIRED = 411, 
    HTTP_INTERNAL_SERVER_ERROR = 500, MISSING_PARAMETER = -1200, 
    INVALID_ID = -1201, FIELD_ERROR = -1203, BAD_REQUEST = -1204, 
    VALUE_ERROR = -1205, VALIDATION_ERROR = -1206, UNSUPPORTED_FORMAT = -1207, 
    INVALID_SORT_ERROR = -1208), .Names = c("HTTP_NO_CONTENT", 
"HTTP_BAD_REQUEST", "HTTP_UNAUTHORIZED", "HTTP_PAYMENT_REQUIRED", 
"HTTP_FORBIDDEN", "HTTP_NOT_FOUND", "HTTP_METHOD_NOT_ALLOWED", 
"HTTP_LENGTH_REQUIRED", "HTTP_INTERNAL_SERVER_ERROR", "MISSING_PARAMETER", 
"INVALID_ID", "FIELD_ERROR", "BAD_REQUEST", "VALUE_ERROR", "VALIDATION_ERROR", 
"UNSUPPORTED_FORMAT", "INVALID_SORT_ERROR"))
.resolve_field_id <-
function (idx, fields) 
{
    id = NULL
    if (class(idx) == "numeric") {
        field = Filter(function(y) y$column_number + 1 == idx, 
            fields)
        id = names(field)[1]
    }
    else if (class(idx) == "character" || class(idx) == "factor") {
        field = Filter(function(y) y$name == idx, fields)
        id = names(field)[1]
    }
    return(id)
}
.resolve_resource_id <-
function (x, name) 
{
    id_str = ""
    if (class(x) == "list") {
        if ("resource" %in% names(x)) {
            return(x$resource)
        }
        else {
            stop("argument doesn't apear to be a valid BigML response")
        }
    }
    else if (class(x) == "character") {
		x = .fixid(x)
        return(paste(name, "/", x, sep = ""))
    }
    else {
        stop("argument is not a string or BigML response")
    }
}
.success_codes <-
structure(list(HTTP_OK = 200, HTTP_CREATED = 201, HTTP_ACCEPTED = 202), .Names = c("HTTP_OK", 
"HTTP_CREATED", "HTTP_ACCEPTED"))

