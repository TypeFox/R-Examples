# test.printcall.R
#
# TODO we don't test use of printcall in a namespace

library(plotmo)

options(warn=1) # print warnings as they occur

printf <- function(format, ...) cat(sprintf(format, ...), sep="") # like c printf

cat0 <- function(...) cat(..., sep="")

# test that we got an error as expected from a try() call
expect.err <- function(object, expected.msg="")
{
    if(class(object)[1] == "try-error") {
        msg <- attr(object, "condition")$message[1]
        if(length(grep(expected.msg, msg, fixed=TRUE)))
            cat0("Got error as expected from ",
                deparse(substitute(object)), "\n")
        else
            stop(sprintf("Expected: %s\n  Got:      %s",
                         expected.msg, substr(msg, 1, 1000)))
    } else
        stop("did not get expected error ", expected.msg)
}
for(all in c(FALSE, TRUE)) {
    for(EVAL in c(FALSE, TRUE)) {
        printf("=== Test printcall with all=%s EVAL=%s ===\n", all, EVAL)

        foo30 <- function() { plotmo:::printcall(all=all) }
        foo30()

        foo32 <- function(...) { plotmo:::printcall(all=all); plotmo:::printdots(..., EVAL=EVAL) }
        foo32()
        foo32(a=31)


        foo34 <- function(aa=1, ...) { plotmo:::printcall(all=all); plotmo:::printdots(..., EVAL=EVAL) }
        foo34()
        foo34(a=31) # argname a will be expanded to aa
        foo34(a=31, x=1:10, y=NULL)
        foo34(a=31, y=NULL)
        foo34(x=stopifnot(TRUE), y=NULL)

        foo36 <- function(aa=NULL, ...) { plotmo:::printcall(all=all); plotmo:::printdots(..., EVAL=EVAL) }
        foo36()
        foo36(a=NULL)
        foo36(a=1)
        foo36(a=1:3)
        foo36(a=1:3, x=NULL)

        # check formatting of various argument types
        # note that we correctly don't call stopifnot(FALSE) (which would call stop)

        foo38 <- function(aa=1:3, bb=4:6, cc=print.default,
                        dd=stopifnot(FALSE),
                        ee=function(m=1) cat(m), ff=7, ...)
            { plotmo:::printcall(all=all); plotmo:::printdots(..., EVAL=EVAL) }
        foo38(x=matrix(ncol=1, nrow=3))

        list1 <- list(aa=1:3, bb=4:6, cc=print.default,
                      dd=stopifnot(TRUE),
                      ee=function(m=1) cat(m), ff=7)

        cat("list1 ", plotmo:::list.as.char(list1), "\n", sep="")

        list2 <- list(lmmod=lm(Volume~Girth, data=trees),
              boolean=c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE), env=parent.frame(),
              chars=c("a", "b", "c", "a", "b", "c"),
              trees=trees, l=list(x=1, y="2", z=foo38))

        cat("list2 ", plotmo:::list.as.char(list2), "\n", sep="")

        # test unnamed arguments

        foo40 <- function(aa, ...) { plotmo:::printcall(all=all); plotmo:::printdots(..., EVAL=EVAL) }
        foo40()
        foo40(aa=b, c)
        foo40(b, c)

        # test printcall when called in an S3 method

        foo.s3 <- function(a=NULL, ...) { UseMethod("foo.s3") }
        foo.s3.list <- function(a=NULL, ...) {
            cat("in foo.s3.list: "); plotmo:::printcall(all=all)
            plotmo:::printdots(..., EVAL=EVAL)
        }
        foo.s3.default <- function(a=NULL, ...) {
            cat("in foo.s3.default: "); plotmo:::printcall(all=all)
            plotmo:::printdots(..., EVAL=EVAL)
        }
        foo.s3(a=list(m=1, n=2))
        foo.s3(a=NULL)
        foo.s3(a=list(m=1, n=2, o=3, p=4, q=5, r=6, s=7, t=8, u=9), b=30)

        # test formatting with long argument list

        foo46 <- function(mmmmmmmmmmm=1000, nnnnnnnnnnn=2000, ooooooooooo=3000, ppppppppppp=4000,
                        qqqqqqqqqqq=5000, rrrrrrrrrrr=6000, sssssssssss=7000, ttttttttttt=8000,
                        uuuuuuuuuuu=9000, vvvvvvvvvvv=1000, wwwwwwwwwww=2000, xxxxxxxxxxx=3000,
                        ...) { plotmo:::printcall(all=all); plotmo:::printdots(..., EVAL=EVAL) }
        foo46(a=30)

        # test call.as.char

        foo47 <- function(aa=1, ...) { s <- plotmo:::call.as.char(all=all); cat(s, "\n", sep="") }
        foo47(b=30)

        # create a variable named foo48 in foo48
        foo48 <- function(aa=1, ...) { foo48 <- 99; s <- plotmo:::call.as.char(all=all); cat(s, "\n", sep="") }
        foo48(b=30)

        # Note that the following doesn't do what you might expect.
        # The calling function is print(), not foo50() as you may expecty.

        foo50 <- function(...) { print(plotmo:::call.as.char(all=all)) }
        foo50(a=1)
    }
}

if(!interactive())
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
