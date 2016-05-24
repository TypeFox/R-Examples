message("TESTING: Static methods and classes...")

clazz <- R.oo::Class
print(clazz)

clazz <- R.oo::Class$forName("Object")
print(clazz)

message("TESTING: Static methods and classes...DONE")
