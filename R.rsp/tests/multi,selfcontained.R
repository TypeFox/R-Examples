library("R.rsp")

path <- system.file("rsp_tests", package="R.rsp")
res <- rfile("multi,selfcontained.md.rsp", path=path, workdir="multi,selfcontained/")
print(res)
