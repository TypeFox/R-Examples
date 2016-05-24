library("R.methodsS3")

message("TESTING: isGenericS3/S4()...")

FUNs <- list(
  isGenericS3=isGenericS3,
  isGenericS4=isGenericS4
)

for (name in names(FUNs)) {
  cat(sprintf("%s():\n", name))
  FUN <- FUNs[[name]]
  print(FUN("print"))
  print(FUN("show"))
  print(FUN("unknown"))
  print(FUN(print))
  print(FUN(sum))
  print(FUN(function() NULL))
}

message("TESTING: isGenericS3/S4()...DONE")
