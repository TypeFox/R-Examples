# install.packages("devtools")
#devtools::install_github("hadley/staticdocs")

library(staticdocs)
build_site(pkg = "~/Dropbox/public/github_params", "../github_paramspages")

## stuff for MAC ONLY
if(Sys.info()['sysname'] == "Darwin"){
	setwd("~/Dropbox/public/github_paramspages/")
	system("git commit -a -m 'update website'")
	system("git push")
}
