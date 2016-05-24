message("TESTING: Package unloading...")

library("R.oo")

pkg <- Package("datasets")
load(pkg)
print(isLoaded(pkg))
unload(pkg)
print(isLoaded(pkg))

message("TESTING: Package unloading...DONE")
