#' Setup and initialize some scripts.
#'
#' @param bin path to bin folder
#'
#' @details Will add more to this to identify cluster and aid in other things
#' @export
setup <- function(bin = "~/bin"){
	pkg = "funr"
	if(!file.exists(bin))
		dir.create(bin) ## create bin, if it does not exist
	script = file.path(bin, pkg)
	if(!file.exists(script)){
		#message("Adding flowr executable to ~/bin")
		file.symlink(system.file(package = pkg, "scripts/funr"), bin)
	}
	tmp <- c("A small script ", pkg, " has been copied to ", bin,
					 ".\nConsider adding ", bin, " to your PATH variable, by running this:",
					 "\necho 'export PATH=$PATH:", bin,"' >> ~/.bashrc",
					 "\nAfter opening a new shell session, and then running ", pkg, " from the shell should work.")
	message(tmp)
	message("\n\nShowing a example ", pkg, " output: ", file.path(bin, pkg))
	cat(system(file.path(bin, pkg), intern = TRUE))
}
