installPckg <- function(requiredPackages){
# Set your mirror of CRAN
old.repos <- getOption("repos") 
on.exit(options(repos = old.repos)) #this resets the repos option on exit
new.repos <- old.repos
new.repos["CRAN"] <- "http://cran.r-project.org" #set your favorite CRAN Mirror here
options(repos = new.repos) # overwrite the repos options

cat("\n----Loading required packages----\n")
cat("\n")

loadedPackages <- try(for (i in 1:length(requiredPackages)){

if (requiredPackages[i] %in% row.names(installed.packages())  == FALSE){

install.packages(requiredPackages[i],dependencies=TRUE)

} else {

	require(requiredPackages[i], quietly = FALSE,
        character.only = TRUE)
        }
     })
}