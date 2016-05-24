message("*** devAll() ...")

library("R.devices")
devAll <- R.devices:::devAll
print(devAll())
print(devAll(force=TRUE))

message("*** devAll() ... DONE")
