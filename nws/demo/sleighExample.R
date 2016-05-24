s = sleigh()

# create workers on different machines
# s1 = sleigh(c('n5', 'n6', 'n7', 'n7'))

paste('add one to a list: 1 + (1:10)')
# add one to a list of 10 numbers
eachElem(s, function(x) {x+1}, list(1:10))


# add two numbers
cat('add two lists: (1:10) + (11:20)\n')
add = function(x, y) {x+y}
elemList = list(1:10, 11:20)
eachElem(s, add, elementArgs=elemList)


# pass fixed argument to the function
cat('use fixed argument to subtract a fixed number: (1:10)-(11:20)-1111\n')
sub = function(x, y, z) {x-y-z}
elemList = list(1:10, 11:20)
fixedArgs = list(1111)
eachElem(s, sub, elementArgs=elemList, fixedArgs=fixedArgs)

# permute the order of function parameters
cat('permute the order of arguments\n')
eo = list(argPermute=c(3, 2, 1))
print_func = function(x, y, z) {paste(x, y, z)}
cat('original order: ')
eachElem(s, print_func, elementArgs=elemList, fixedArgs=fixedArgs, eo=eo)
paste('after permute: ')
eachElem(s, sub, elementArgs=elemList, fixedArgs=fixedArgs, eo=eo)

# stop sleigh
cat('stop sleigh\n')
stopSleigh(s)


