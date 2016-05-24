library("R.cache")

root <- getCacheRootPath()
print(root)

path <- getCachePath(dirs=c("abc", "def"))
print(path)

path <- getCachePath(dirs=c("abc", "def"), path="subdir")
print(path)
