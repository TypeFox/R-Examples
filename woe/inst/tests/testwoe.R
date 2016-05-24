
test1<-woe(mtcars,"mpg",TRUE,"am",10,Bad=0,Good=1)

test2<-woe(mtcars,"cyl",FALSE,"am",10,Bad=0,Good=1)

test3<-woe(mtcars,"mpg_wrong",TRUE,"am",10,Bad=0,Good=1)

test4<-woe(mtcars,"cyl_wrong",FALSE,"am",10,Bad=0,Good=1)

expect_equal(class(test1),"data.frame")
expect_equal(class(test2),"data.frame")


expect_equal(class(test3),"numeric")
expect_equal(class(test4),"numeric")


