getPathname <- R.devices:::getPathname
getPath <- R.devices:::getPath
view <- R.devices:::view
getData <- R.devices:::getData

message("*** DevEvalFileProduct ...")

# Return the DevEvalProduct object by default
R.devices::devOptions("*", field="*")

message("*** DevEvalFileProduct - image file ...")

p <- R.devices::toEPS("foo", tags=c("a", "b"), aspectRatio=0.7, {
  plot(1:10)
})
print(p)
str(p)

library("R.devices")

fields <- c("name", "fullname", "filename", "pathname", "path", "mime", "dataURI")
for (ff in fields) {
  cat(sprintf("%s: %s\n", ff, substring(p[[ff]], 1, 50)))
}

pathnameA <- getPathname(p, relative=FALSE)
cat(sprintf("Absolute pathname: %s\n", pathnameA))
path <- getPath(p, relative=FALSE)
cat(sprintf("Path: %s\n", path))
pathnameR <- getPathname(p, relative=TRUE)
cat(sprintf("Relative pathname: %s\n", pathnameR))

data <- getData(p, mode="character")
str(data)

data <- getData(p, mode="raw")
str(data)


## Call view() but use void browser
view(p, browser="false")

message("*** DevEvalFileProduct - image file ... DONE")

message("*** DevEvalFileProduct - missing file ...")

# An empty file product
na <- DevEvalFileProduct()
print(na)
cat(sprintf("Pathname: %s\n", getPathname(na)))

message("*** DevEvalFileProduct - missing file ... DONE")

message("*** DevEvalFileProduct ... DONE")
