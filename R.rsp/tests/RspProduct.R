library("R.rsp")

p <- RspProduct()
print(p)

p <- RspStringProduct()
print(p)

p <- RspFileProduct("dummy", mustExist=FALSE)
print(p)

p <- RspSourceCode()
print(p)

p <- RspRSourceCode()
print(p)

p <- RspShSourceCode()
print(p)

