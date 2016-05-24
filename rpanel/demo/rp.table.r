panel <- rp.control()
rp.table(
  panel=panel,
  theexampletable, 
  matrix(
    c(1,2,3, 11,12,13, 4,5,6, 7,8,9), nrow=4, ncol=3, byrow=TRUE, 
    dimnames=list(c("row1", "row2", "row3", "row4"), c("C.1", "C.2", "C.3"))),
  action_click=function(panel) { print("click"); print(panel$theexampletable); panel },
  action_press=function(panel) { print("press"); print(panel$theexampletable); panel },
  action_return=function(panel) { print("return"); print(panel$theexampletable); panel },
  name="exampletable"
)
print(rp.table.value.get(panel, "exampletable", 1, 2))
rp.table.value.set(panel, "exampletable", 1,2,10)
print(rp.table.value.get(panel, "exampletable", 1, 2))
print(rp.table.get(panel, "exampletable"))