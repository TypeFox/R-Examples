library("R.rsp")

path <- system.file("exData", package="R.rsp")
pathname <- rfile("random.txt.rsp", path=path)
print(pathname)

lines <- readLines(pathname, warn=FALSE)
cat(lines, collapse="\n")

file.remove(pathname)


# Passing arguments
path <- system.file("exData", package="R.rsp")
pathname <- rfile("random-args.txt.rsp", path=path, args=list(K=50))
print(pathname)

lines <- readLines(pathname, warn=FALSE)
cat(lines, collapse="\n")

file.remove(pathname)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Various host content types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- system.file("exData", package="R.rsp")
filenames <- list.files(pattern="LoremIpsum.*.rsp$", path=path)
for (filename in filenames) {
  print(filename)
  pathname <- rfile(filename, path=path)
  print(pathname)
}
