sysinfo<-function()
  {


A = .Machine

B = .Platform


print(paste(sep=" ", "Platform=", B$OS.type))
print(paste(sep=" ", "Endian=", B$endian))
print(paste(sep=" ", "Size of LONG=", A$sizeof.long))

invisible(list(A=A, B=B))
  }


