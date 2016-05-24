
 library(rcdd)

 foo <- z2q(c(1:5), c(6:10))

 bar <- z2q(c(0:4), c(11:15))

 print(foo)

 print(bar)

 ##### addition #####

 .Call("qoq", foo, bar, as.integer(1))

 qpq(foo, bar)

 ##### subtraction #####

 .Call("qoq", foo, bar, as.integer(2))

 qmq(foo, bar)

 ##### multiplication #####

 .Call("qoq", foo, bar, as.integer(3))

 qxq(foo, bar)

 ##### division #####

 .Call("qoq", bar, foo, as.integer(4))

 qdq(bar, foo)

 options(error=dump.frames)

 .Call("qoq", foo, bar, as.integer(4))

 qdq(foo, bar)

