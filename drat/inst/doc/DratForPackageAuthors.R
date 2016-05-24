## ---- eval=FALSE---------------------------------------------------------
#  ## insert given package into default drat repo on local file system
#  drat::insertPackage("myPkg_0.5.tar.gz")

## ---- eval=FALSE---------------------------------------------------------
#  ## insert given package into given repo on local file system
#  drat::insertPackage("myPkg_0.5.tar.gz", "/srv/projects/git/drat")

## ---- eval=FALSE---------------------------------------------------------
#  ## insert given package into given repo on a network-local file system
#  drat::insertPackage("myPkg_0.5.tar.gz", "file:/nfs/groups/groupABC/R/drat")

