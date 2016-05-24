library("R.rsp")

pathname <- system.file("exData", "slowcounting.txt.rsp", package="R.rsp")
rstr <- rclean(pathname, verbose=TRUE)
print(rstr)
