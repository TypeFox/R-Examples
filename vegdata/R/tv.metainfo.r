tv.metadata <- function (db, refl, tv_home, filename = 'metadata.txt', ...)
{
    if (missing(tv_home)) tv_home <- tv.home()
    if (db[1] == "eco") {
	if(missing(refl)) refl <- tv.refl(db = db[1], tv_home = tv_home, ...)
	shell.exec(file.path(tv_home, "Species", refl, "metadata-eco.txt"))
      } else {
	for(i in 1:length(db)) {
	  meta <- file.path(tv_home, "Data", db[i], filename)
	  if (file.access(meta)) stop('No metainfo file "',filename, '" available in directory "', db[1], '".') else 
	  if(.Platform$OS.type == "windows") shell.exec(meta)  else file.show(meta)
  }}
}
