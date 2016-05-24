# ' Upgrades a package to the lastest version on R-forge
# ' @param package name of the packages to upgrade
# ' @param update boolean to say if it is an update (detach the package)
# ' @export
rforge<-function(package="CorReg",update=FALSE){
   if(package=="CorReg" |update){detach(unload=TRUE)}
   install.packages(package, repos="http://R-Forge.R-project.org")
}

