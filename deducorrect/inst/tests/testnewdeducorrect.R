
context("Deducorrect object")



test_that('newdeducorrect works if no corrections have been applied',{
    dd <- data.frame(x=1:3,y=rep(0,3))
    st <- data.frame(status=status(3,ini='invalid'))
    cr <- data.frame(row=numeric(0),variable=character(0),old=numeric(0),new=numeric(0))
    newdeducorrect(corrected=dd,status=st,corrections=cr)
})







