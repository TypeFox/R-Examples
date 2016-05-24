.RZabbixEnv <- new.env()

.onAttach <- function(...) {
   packageStartupMessage("RZabbix ", utils::packageVersion("RZabbix"))
}

.onLoad <- function(...) {
	assign('id', value = 1, envir = .RZabbixEnv)
}