# Create symlink to file with a simplified filename
# @param file  a \code{filename} object
.symlink_simplified <- function(file) {
	if (is.character(file)) {
		file <- as.filename(file);
	}
	if (is.filename(file)) {
		x <- as.character(file, simplify=FALSE);
		y <- as.character(file, simplify=TRUE);
		if (x != y) {
			info <- Sys.readlink(y);
			# strip away all but the last directory (date directory)
			x.rel <- as.character(
				set_fpath(file, file$path[length(file$path)]), simplify=FALSE);
			if (is.na(info)) {
				# target does not exist: no preparation needed
				# create symlink to the relative path
				file.symlink(x.rel, y);
			} else if (info == "") {
				# target exists and it is not a symbolic link
				warning("Cannot symlink ", x.rel, " to ", y,
					", because a file exists with the same name.");
			} else {
				# target exists as a symlink: remove it
				file.remove(y);
				file.symlink(x.rel, y);
			}
		}
	}
}
