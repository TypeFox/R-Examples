## library(appell)
## The following code compares results with those published in
## Colavecchia et al. (2001), tables 2-5.
## It also illustrates the very minor differences between results
## obtained with the hyp2f1 algorithm by R. Forrey and or that by
## N. Michel and M. Stoitsov, with the exception of no convergence
## problems with the latter. Therefore it is used as default algorithm. 

## --------------------
## read the original table 2
table2orig <- read.table(file=
                         system.file(file.path("extdata", "table2.dat"),
                                     package="appell"),
                         col.names=c("x", "y", "absf1", "exact", "relError"))
table2orig

## compute the values here
table2orig <- cbind(table2orig,
                    absf1.forrey=0,
                    absf1.michel.stoitsov=0)

for(case in seq_len(nrow(table2orig)))
{
    table2orig[case, "absf1.forrey"] <-
        abs(appellf1(a=1,
                     b1=2+1i,
                     b2=1.5-0.5i,
                     c=1,
                     x=table2orig[case, "x"],
                     y=table2orig[case, "y"],
                     debug=TRUE,        # test debugging info as well
                     userflag=1L,
                     hyp2f1="forrey")$val)
    table2orig[case, "absf1.michel.stoitsov"] <-
        abs(appellf1(a=1,
                     b1=2+1i,
                     b2=1.5-0.5i,
                     c=1,
                     x=table2orig[case, "x"],
                     y=table2orig[case, "y"],
                     userflag=1L,
                     hyp2f1="michel.stoitsov")$val)
}
table2orig

## look at the (small) differences:
table2orig$absf1 - table2orig$absf1.forrey
table2orig$absf1 - table2orig$absf1.michel.stoitsov

## here no difference between the hyp2f1 choices:
identical(table2orig$absf1.forrey,
          table2orig$absf1.michel.stoitsov)

## --------------------
## read the original table 3
table3orig <- read.table(file=
                         system.file(file.path("extdata", "table3.dat"),
                                     package="appell"),
                         col.names=
                         c("x", "y", "f1ser", "f1int", "f1exact", "relErrorser",
                           "relErrorint"))

## compute the values here
table3orig <- cbind(table3orig,
                    f1ser.forrey=0,
                    f1int.forrey=0,
                    f1ser.michel.stoitsov=0,
                    f1int.michel.stoitsov=0)

for(case in seq_len(nrow(table3orig)))
{
    ## first everything with Forrey's algorithm for hyp2f1
    f1ser.forrey <- try(appellf1(a=1,
                                 b1=3+1i,
                                 b2=2-0.5i,
                                 c=5+0.5i,
                                 x=table3orig[case, "x"],
                                 y=table3orig[case, "y"],
                                 userflag=2L,
                                 debug=TRUE,
                                 hyp2f1="forrey"))
    
    table3orig[case, "f1ser.forrey"] <-
        if(inherits(f1ser.forrey, "try-error"))
            NA
        else
            abs(f1ser.forrey$val)

    table3orig[case, "f1int.forrey"] <- abs(appellf1(a=1,
                                                     b1=3+1i,
                                                     b2=2-0.5i,
                                                     c=5+0.5i,
                                                     x=table3orig[case, "x"],
                                                     y=table3orig[case, "y"],
                                                     userflag=1L,
                                                     debug=TRUE,
                                                     hyp2f1="forrey")$val)

    ## then everything with the algorithm by Michel and Stoitsov for hyp2f1
    f1ser.michel.stoitsov <- try(appellf1(a=1,
                                 b1=3+1i,
                                 b2=2-0.5i,
                                 c=5+0.5i,
                                 x=table3orig[case, "x"],
                                 y=table3orig[case, "y"],
                                 userflag=2L,
                                 debug=TRUE,
                                 hyp2f1="michel.stoitsov"))
    
    table3orig[case, "f1ser.michel.stoitsov"] <-
        if(inherits(f1ser.michel.stoitsov, "try-error"))
            NA
        else
            abs(f1ser.michel.stoitsov$val)

    table3orig[case, "f1int.michel.stoitsov"] <- abs(appellf1(a=1,
                                                     b1=3+1i,
                                                     b2=2-0.5i,
                                                     c=5+0.5i,
                                                     x=table3orig[case, "x"],
                                                     y=table3orig[case, "y"],
                                                     userflag=1L,
                                                     debug=TRUE,
                                                     hyp2f1="michel.stoitsov")$val)
}
table3orig

## look at the (small) differences:
table3orig$f1ser - table3orig$f1ser.michel.stoitsov
table3orig$f1int - table3orig$f1int.michel.stoitsov
## so we have no missing values for Michel & Stoitsov

## besides that, only very small differences between the two methods:
table3orig$f1ser.michel.stoitsov - table3orig$f1ser.forrey
table3orig$f1int.michel.stoitsov - table3orig$f1int.forrey

## --------------------
## read the original table 4
table4orig <- read.table(file=
                         system.file(file.path("extdata", "table4.dat"),
                                     package="appell"),
                         col.names=
                         c("x", "y", "f1", "hypergeo", "exact", "relErrorf1",
                           "relErrorhypergeo"))

## compute the values here
table4orig <- cbind(table4orig,
                    f1.forrey=0,
                    f1.michel.stoitsov=0,
                    flag=0)

for(case in seq_len(nrow(table4orig)))
{
    ## get Forrey result and flag
    thisRes <- appellf1(a=-0.5,
                        b1=2,
                        b2=1,
                        c=3,
                        x=table4orig[case, "x"],
                        y=table4orig[case, "y"],
                        hyp2f1="forrey")
    
    table4orig[case, "f1.forrey"] <- abs(thisRes$val)
    table4orig[case, "flag"] <- thisRes$algoflag

    ## get Michel & Stoitsov result
    table4orig[case, "f1.michel.stoitsov"] <-
        abs(appellf1(a=-0.5,
                     b1=2,
                     b2=1,
                     c=3,
                     x=table4orig[case, "x"],
                     y=table4orig[case, "y"],
                     hyp2f1="michel.stoitsov")$val)
}
table4orig

## look at the (small) differences:
table4orig$f1 - table4orig$f1.forrey
## very small errors all over the place!

## and extremely small differences between Forrey and Michel & Stoitsov:
table4orig$f1.michel.stoitsov - table4orig$f1.forrey


## look at the flags
subset(table4orig,
       select=c(x, y, flag))
