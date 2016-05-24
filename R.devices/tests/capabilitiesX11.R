message("*** capabilitiesX11() ...")

res <- R.devices::capabilitiesX11()
print(res)

res2 <- R.devices::capabilitiesX11()
stopifnot(identical(res2, res))

res3 <- R.devices::capabilitiesX11(reset=TRUE)
print(res3)

message("*** capabilitiesX11() ... DONE")
