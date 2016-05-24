#' @method qread yaml
#' @export
qread.yaml <- function(file, type, ...) {
	.check_package("yaml");

	yaml::yaml.load_file(file, ...)
}

#' @method qwrite yaml
#' @export
qwrite.yaml <- function(x, file, type, append=FALSE, ...) {
	.check_package("yaml");

	cat(yaml::as.yaml(x, ...), sep="", file=file, append=append)
}
