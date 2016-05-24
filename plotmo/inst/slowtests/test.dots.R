# test.dots.R

library(plotmo)

options(warn=1) # print warnings as they occur

if(!interactive())
    postscript(paper="letter")

strip.space <- function(s) gsub("[ \t\n]", "", s)

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
cat0("=== test callers.name\n")

test.callers.name <- function(x) {
    caller0  <- plotmo:::callers.name(0)  # test.callers.name
    caller1  <- plotmo:::callers.name(1)  # caller of test.callers.name
    caller99 <- plotmo:::callers.name(99) # sys.call(-n) : not that many frames on the stack
    s <- sprintf("0 %s 1 %s 99 %s", caller0, caller1, caller99)
    cat(s, "\n", sep="")
    s
}
print(plotmo:::callers.name()) # "eval"

stopifnot(test.callers.name() == "0 test.callers.name 1 stopifnot 99 unknown")

myfunc <- function(func) func()
stopifnot(myfunc(function(x) test.callers.name(99)) == "0 test.callers.name 1 func 99 unknown")

cat0("=== test dotindex\n")

test.dotindex <- function(expected, ARGNAME, ..., EX=FALSE)
{
    dotindex <- plotmo:::dotindex(ARGNAME=ARGNAME, EX=EX, ...)
    stopifnot(all.equal(dotindex, expected))
}
test.dotindex(NA, "x") # empty dots
test.dotindex(NA, "x",  a=10, b=20)
test.dotindex(1,  "a",  a=10, b=20)
test.dotindex(2,  "b",  a=10, b=20)
test.dotindex(1,  "a1", a=10, b=20)
test.dotindex(NA, "a",  a1=10, a2=20)
expect.err(try(test.dotindex(1, nonesuch, a=10, a=20)), "object 'nonesuch' not found")
expect.err(try(test.dotindex(1, "a1", a=10, a=20)), "argument 'a' for test.dotindex() is duplicated")
expect.err(try(test.dotindex(1, "aa1", a=10, aa=20)), "arguments 'a' and 'aa' both match 'aa1' in test.dotindex")
stopifnot(is.na(plotmo:::dotindex("a", EX=1, a1=10, a2=20)))
stopifnot(plotmo:::dotindex("a2", EX=1, a1=10, a2=20) == 2)

# multiple argnames
test.dotindex(NA, c("a", "b")) # empty dots
test.dotindex(1,  c("a", "b"), a=2, c=3)
test.dotindex(1,  c("a", "b"), a=5, b=6)
test.dotindex(2,  c("a", "b"), x=1, a=5, b=6)
test.dotindex(3,  c("b,a"), x=1, a=5, b=6)
test.dotindex(1,  c("a b"), b=3, c=4)
test.dotindex(2,  c(" a b "), c=3, b=4)
test.dotindex(NA, c("a", "b"), c=3)
stopifnot(plotmo:::dotindex(c("x", "a1"), EX=1, a1=10, a2=20) == 1)

test.dot <- function(expected, ARGNAME, ..., DEF=NA, EX=FALSE)
{
    if(is.na(DEF))
        dot <- plotmo:::dot(ARGNAME, EX=EX, ...)
    else
        dot <- plotmo:::dot(ARGNAME, EX=EX, DEF=DEF, ...)
    stopifnot(all.equal(dot, expected))
}
cat0("=== test dot\n")
test.dot(NA, "x") # empty dots
test.dot(NA, "x",  a=10, b=20)
test.dot(10, "a",  a=10, b=20)
test.dot(20, "b",  a=10, b=20)
test.dot(99, DEF=99, "nonesuch", a=10, b=20)
test.dot(NA, "a", a1=10, a2=20)
expect.err(try(test.dot(1, "a1", a=10, a=20)), "argument 'a' for test.dot() is duplicated")
expect.err(try(test.dot(1, 99, a=10, a=20)), "is.character(argname) is not TRUE")
expect.err(try(test.dot(1, test.dot, a=10, a=20)), "is.character(argname) is not TRUE")
expect.err(try(test.dot(1, "", a=10, a=20)), "empty string in ARGNAME")
expect.err(try(test.dot(1, "x^x", a=10, a=20)), "illegal character \"^\" in ARGNAME")

test.dot(10, "abc", EX=T, abc=10)
test.dot(NA, "a",   EX=T, a1=10, a2=20)
expect.err(try(test.dot(1, "a1", a1=10, a1=20)), "argument 'a1' for test.dot() is duplicated")

stopifnot(is.na(plotmo:::dot("a", EX=1, a1=1, a2=2)))
stopifnot(plotmo:::dot("a2", EX=1, a1=10, a2=20, a3=30) == 20)

foo <- function(func, x) func(x)
foo(mean, 33)
foo(function(...) plotmo:::dot("x", ...), 33)
foo(function(...) plotmo:::dot("x99", ...), 33)
foo(function(...) { plotmo:::dot("nonesuch", ...) }, 33)

test.dot(1,  "a", EX=T, a=1)
test.dot(2,  "b", EX=T, a=1, b=2, c=3)
test.dot(NA, "x", EX=T, a=1, b=2, c=3)
test.dot(2,  "a", EX=T, ab=1, a=2)
test.dot(2,  "a", EX=T, aa=1, a=2)
test.dot(NA, "a", EX=T, aa=1, ab=2)
expect.err(try(test.dot(2, "a", EX=T, aa=1, a=2, a=3)), "argument 'a' for test.dot() is duplicated")

expect.err(try(test.dot(2, "a", EX=T, a=none.such)), "cannot evaluate 'a'")

# multiple argnames
test.dot(2,  c("a", "b"), a=2, c=3)
test.dot(5,  c("a", "b"), a=5, b=6)
test.dot(5,  c("a", "b"), x=1, a=5, b=6)
test.dot(3,  c("a", "b"), b=3, c=4)
test.dot(4,  c("a", "b"), c=3, b=4)
test.dot(NA, c("a", "b"), c=3)
expect.err(try(test.dot(1, c("b", "aa1"), a=10, aa=20)), "arguments 'a' and 'aa' both match 'aa1' in test.dot")
expect.err(try(test.dot(1, c("x", ""), a=10, b=20)), "empty string in ARGNAME")
stopifnot(plotmo:::dot(c("x", "a2", "y"), EX=1, a1=10, a2=20, a3=30) == 20)

test.dot(NA, c("a", "b"), aa=2, cc=3, EX=T)
test.dot(2,  c("aa", "b"), aa=2, cc=3, EX=T)
test.dot(3,  c("bb", "b"), bb=3, cc=4, EX=T)
test.dot(NA, c("a", "b"), c=3, EX=T)

foo.x <- function(...) { plotmo:::dot("x", ..., DEF="default", EX=FALSE) }
stopifnot(foo.x(x=3) == 3)
stopifnot(foo.x(y=3) == "default")

foo2 <- function(funcarg, ...) funcarg(...)
stopifnot(is.na(foo2(function(...) plotmo:::dot("x", ...), 3))) # 3 is unnamed
stopifnot(foo2(function(...) plotmo:::dot("x", EX=0, ...), x=3) == 3)
stopifnot(foo2(function(...) plotmo:::dot("x99", EX=0, ...), x=3) == 3)
stopifnot(foo2(function(...) { plotmo:::dot("x", DEF="default", EX=FALSE, ...) }, x=3) == 3)
stopifnot(foo2(function(...) { plotmo:::dot("y", DEF="default", EX=FALSE, ...) }, x=3) == "default")
# expect.err(try(foo2(function(...) { plotmo:::dot("y", DEF="default", EX=FALSE, ...) }, 3)), "unnamed arguments in ... are not allowed for funcarg()")

stopifnot(foo2(foo.x, x=3) == 3)
stopifnot(foo2(foo.x, y=3) == "default")

test.is.dot <- function(expected, ARGNAME, ...)
{
    present <- plotmo:::is.dot(ARGNAME, ...)
    stopifnot(all.equal(present, expected))
}
cat0("=== test is.dot\n")
test.is.dot(FALSE, "x") # empty dots
test.is.dot(FALSE, "x",  EX=0, a=10, b=20)
test.is.dot(TRUE,  "a",  EX=0, a=10, b=20)
test.is.dot(TRUE,  "b",  EX=0, a=10, b=20)
test.is.dot(TRUE,  "a1", EX=0, a=10, b=20)
test.is.dot(FALSE, "a",  EX=0, a1=10, a2=20)
expect.err(try(test.is.dot(TRUE, "a1", EX=0, a=10, a=20)), "argument 'a' for test.is.dot() is duplicated")
expect.err(try(test.is.dot(TRUE, "a", EX=0, a=10, a=20)), "argument 'a' for test.is.dot() is duplicated")
stopifnot(plotmo:::is.dot("a",  EX=1, a1=10, a2=20, a3=30) == FALSE)
stopifnot(plotmo:::is.dot("x",  EX=1, a1=10, a2=20, a3=30) == FALSE)
stopifnot(plotmo:::is.dot("a3", EX=1, a1=10, a2=20, a3=30) == TRUE)

# multiple argnames
test.is.dot(TRUE,  EX=0, c("a1", "b1"), a=2, c=3)
test.is.dot(TRUE,  EX=0, c("a1", "b1"), b=3, c=4)
test.is.dot(TRUE,  EX=0, c("a1", "b1"), c=3, b=4)
test.is.dot(FALSE, EX=0, c("a1", "b1"), c=3)
expect.err(try(test.is.dot(FALSE, c("aa1", "b"), EX=0, a=10, aa=20)), "arguments 'a' and 'aa' both match 'aa1' in test.is.dot")
stopifnot(plotmo:::is.dot(c("x", "a", "y"), EX=1, a1=10, a2=20, a3=30) == FALSE)
stopifnot(plotmo:::is.dot(c("x", "a2", "y"), EX=1, a1=10, a2=20, a3=30) == TRUE)

cat0("=== test expand.drop\n")

# nchar is used an example func, it has formals "x", "type", "allowNA"

stopifnot(is.null(plotmo:::expand.drop(NULL, prefix="prefix.", func=nchar)))

stopifnot(plotmo:::expand.drop("a", prefix="prefix.", func=nchar) == ">PREFIX|>EXPLICIT|^a")

stopifnot(plotmo:::expand.drop("a", prefix="prefix.", func=nchar, include.standard.prefixes=TRUE) == ">STANDARDPREFIXES|^force\\.|^def\\.|^drop\\.|>PREFIX|^prefix\\.|>EXPLICIT|^a")

stopifnot(plotmo:::expand.drop("FORMALS", prefix="prefix.", func=base::nchar) == ">FORMALS|^x|^type|^allowNA|^keepNA|>PREFIX|>EXPLICIT")

stopifnot(plotmo:::expand.drop("FORMALS", prefix="prefix.", func=base::nchar, include.standard.prefixes=TRUE) == ">FORMALS|^x|^type|^allowNA|^keepNA|>STANDARDPREFIXES|^force\\.|^def\\.|^drop\\.|>PREFIX|^prefix\\.|>EXPLICIT")

expect.err(try(plotmo:::expand.drop("FORMALS", prefix="prefix.", func=NULL)), "\"FORMALS\" specified in DROP, but FUNC is NULL")

expect.err(try(plotmo:::expand.drop("FORMALS", prefix="prefix.", func=base::c)), "\"FORMALS\" specified but formals(FUNC) returned no formal arguments")

foo99 <- function(...) NULL
expect.err(try(plotmo:::expand.drop("FORMALS", prefix="prefix.", func=foo99)), "\"FORMALS\" specified but formals(FUNC) returned only \"...\"")

stopifnot(plotmo:::expand.drop("a,FORMALS", prefix="prefix.", func=base::nchar) == ">FORMALS|^x|^type|^allowNA|^keepNA|>PREFIX|>EXPLICIT|^a")

stopifnot(plotmo:::expand.drop("a,FORMALS", prefix="prefix.", func=base::nchar, include.standard.prefixes=TRUE) == ">FORMALS|^x|^type|^allowNA|^keepNA|>STANDARDPREFIXES|^force\\.|^def\\.|^drop\\.|>PREFIX|^prefix\\.|>EXPLICIT|^a")

expect.err(try(plotmo:::expand.drop("", prefix="prefix.", func=base::nchar)), "DROP is an empty string")

stopifnot(plotmo:::expand.drop("a", prefix="lines.", func=base::nchar) == ">PREFIX|>EXPLICIT|^a")

stopifnot(plotmo:::expand.drop("a", "lines.a", prefix="lines.", func=base::nchar, include.standard.prefixes=TRUE) == ">STANDARDPREFIXES|^force\\.|^def\\.|^drop\\.|>PREFIX|^lines\\.|>EXPLICIT|^a")

stopifnot(plotmo:::expand.drop("a*", prefix="lines.", func=base::nchar) == ">PREFIX|>EXPLICIT|^a.*")

stopifnot(plotmo:::expand.drop("a.*", prefix="lines.", func=base::nchar) == ">PREFIX|>EXPLICIT|^a\\..*")

stopifnot(plotmo:::expand.drop("a$", prefix="lines.", func=base::nchar) == ">PREFIX|>EXPLICIT|^a$")

stopifnot(plotmo:::expand.drop("a$,b*,c*$", prefix="lines.", func=base::nchar) == ">PREFIX|>EXPLICIT|^a$|^b.*|^c.*$")

stopifnot(plotmo:::expand.drop(c("a", "b,c", " d e$ f ", "g h$, i"), prefix="lines.", func=base::nchar) ==
">PREFIX|>EXPLICIT|^a|^b|^c|^d|^e$|^f|^g|^h$|^i")

stopifnot(plotmo:::expand.drop("PLOT.ARGS", prefix="lines.", func=base::nchar) ==
">PREFIX|>EXPLICIT|>PLOT_ARGS|^add$|^adj$|^bty$|^cex$|^cex\\.axis$|^cex\\.lab$|^cex\\.main$|^cex\\.sub$|^col$|^col\\.axis$|^col\\.lab$|^col\\.main$|^col\\.sub$|^crt$|^family$|^font$|^font$|^font\\.axis$|^font\\.lab$|^font\\.main$|^font\\.sub$|^lend$|^ljoin$|^lmitre$|^lty$|^lwd$|^main$|^pch$|^srt$|^xaxp$|^xaxs$|^xaxt$|^xlab$|^xlim$|^xlog$|^xpd$|^yaxp$|^yaxs$|^yaxt$|^ylab$|^ylim$|^ylog$")

stopifnot(plotmo:::expand.drop("abc,PLOT.ARGS", prefix="lines.", func=base::nchar) ==
">PREFIX|>EXPLICIT|^abc|>PLOT_ARGS|^add$|^adj$|^bty$|^cex$|^cex\\.axis$|^cex\\.lab$|^cex\\.main$|^cex\\.sub$|^col$|^col\\.axis$|^col\\.lab$|^col\\.main$|^col\\.sub$|^crt$|^family$|^font$|^font$|^font\\.axis$|^font\\.lab$|^font\\.main$|^font\\.sub$|^lend$|^ljoin$|^lmitre$|^lty$|^lwd$|^main$|^pch$|^srt$|^xaxp$|^xaxs$|^xaxt$|^xlab$|^xlim$|^xlog$|^xpd$|^yaxp$|^yaxs$|^yaxt$|^ylab$|^ylim$|^ylog$")

stopifnot(plotmo:::expand.drop("abc,FORMALS,PLOT.ARGS", prefix="lines.", func=base::nchar) ==
">FORMALS|^x|^type|^allowNA|^keepNA|>PREFIX|>EXPLICIT|^abc|>PLOT_ARGS|^add$|^adj$|^bty$|^cex$|^cex\\.axis$|^cex\\.lab$|^cex\\.main$|^cex\\.sub$|^col$|^col\\.axis$|^col\\.lab$|^col\\.main$|^col\\.sub$|^crt$|^family$|^font$|^font$|^font\\.axis$|^font\\.lab$|^font\\.main$|^font\\.sub$|^lend$|^ljoin$|^lmitre$|^lty$|^lwd$|^main$|^pch$|^srt$|^xaxp$|^xaxs$|^xaxt$|^xlab$|^xlim$|^xlog$|^xpd$|^yaxp$|^yaxs$|^yaxt$|^ylab$|^ylim$|^ylog$")

stopifnot(plotmo:::expand.drop("abc,FORMALS,PAR.ARGS", prefix="lines.", func=base::nchar) ==
">FORMALS|^x|^type|^allowNA|^keepNA|>PREFIX|>EXPLICIT|^abc|>PAR_ARGS|^adj$|^ann$|^ask$|^bg$|^bty$|^cex$|^cex\\.axis$|^cex\\.lab$|^cex\\.main$|^cex\\.sub$|^col\\.axis$|^col\\.lab$|^col\\.main$|^col\\.sub$|^crt$|^err$|^family$|^fg$|^fig$|^fin$|^font$|^font\\.axis$|^font\\.lab$|^font\\.main$|^font\\.sub$|^lab$|^las$|^lend$|^lheight$|^ljoin$|^lmitre$|^lty$|^mai$|^mar$|^mex$|^mfcol$|^mfg$|^mfrow$|^mgp$|^mkh$|^new$|^oma$|^omd$|^omi$|^pch$|^pin$|^plt$|^ps$|^pty$|^srt$|^tck$|^tcl$|^usr$|^xaxp$|^xaxs$|^xaxt$|^xlog$|^xpd$|^yaxp$|^yaxs$|^yaxt$|^ylbias$|^ylog$")

stopifnot(plotmo:::expand.drop("abc,FORMALS,PLOTMO.ARGS", prefix="lines.", func=base::nchar) ==
">FORMALS|^x|^type|^allowNA|^keepNA|>PREFIX|>EXPLICIT|^abc|>PLOTMO_ARGS|^caption\\.|^cex\\.|^col\\.|^contour\\.|^cum\\.|^degree1\\.|^degree2\\.|^density\\.|^filled\\.contour\\.|^font\\.|^func\\.|^grid\\.|^heatmap\\.|^image\\.|^jitter\\.|^legend\\.|^label\\.|^level\\.|^line\\.|^lines\\.|^lty\\.|^lty\\.|^lwd\\.|^main\\.|^mtext\\.|^nresiduals|^par\\.|^pch\\.|^persp\\.|^plot\\.|^plotmath\\.|^qq\\.|^qqline\\.|^pt\\.|^response\\.|^rug\\.|^smooth\\.|^text\\.|^title\\.|^vfont\\.")

test.deprefix <- function(expected, ..., FNAME="test.deprefix", KEEP=NULL)
{
    args <- plotmo:::deprefix(..., FNAME=FNAME, KEEP=KEEP, CALLARGS="")
    # can't use all.equal because it complains about names
    # cat("args:\n")
    # print(args)
    # cat("expected:\n")
    # print(expected)
    stopifnot(length(args) == length(expected))
    for(i in seq_len(length(expected))) {
        stopifnot(names(args)[i] == names(expected)[i])
        stopifnot(args[[i]] == expected[[i]])
    }
}
cat0("=== test deprefix\n")

test.deprefix(
    expected=list(a=1, b=2), DROP="*",
    PREFIX="predict.", def.a=1, predict.b=2, c=3)

test.deprefix(TRACE=2,
    expected=list(b="predict.b", d="def.d", c="predict.c", e="predict.e"),
    PREFIX="predict.", DROP="*",
    a="a", b="b", c="c", w1.xlab="xlab",
    def.b="def.b", def.d="def.d",
    predict.b="predict.b", predict.c="predict.c", predict.e="predict.e")

test.deprefix(TRACE=2,
    expected=list(b="predict.b", d="def.d", a="a", c="predict.c", e="predict.e"),
    KEEP=NULL, PREFIX="predict.", DROP="w1.",
    a="a", b="b", c="c", w1.xlab="xlab",
    def.b="def.b", def.d="def.d",
    predict.b="predict.b", predict.c="predict.c", predict.e="predict.e")

test.deprefix(
    expected=list(a="predict.a"),
    KEEP=NULL, PREFIX="predict.", DROP="w1.",
    a="plain.a", predict.a="predict.a")

test.deprefix(expected=list(a="aa1"),
    KEEP=NULL, PREFIX="predict.", a="aa1")

test.deprefix(expected=list(a="aa2"),
    KEEP=NULL, PREFIX="predict.", def.a="aa2")

test.deprefix(expected=list(a="aa3", b="bb3"),
    KEEP=NULL, PREFIX="predict.", def.a="aa3", b="bb3")

test.deprefix(expected=list(10, 20), TRACE=2,
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon2=20, force.anon1=10)

test.deprefix(expected=list(10, 20, a=3), TRACE=2,
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon2=20, force.anon1=10,
                a=3)

expect.err(try(test.deprefix(expected=list(10, 20, a=4),
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon=10, force.anon=20,
                a=3, predict.a=4)),
                "argument 'force.anon' for test.deprefix() is duplicated")

expect.err(try(test.deprefix(expected=list(10, 20, a=4),
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.", FNAME="foobar",
                force.anon=10, force.anon=20,
                a=3, predict.a=4)),
                "argument 'force.anon' for foobar() is duplicated")

test.deprefix(expected=list(10, 20, a=4),
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon1=10, force.anon2=20,
                a=3, predict.a=4)

test.deprefix(expected=list(10, 20, b=3, a=4),
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon1=10, force.anon2=20, def.b=3,
                a=3, predict.a=4)

test.deprefix(expected=list(10, 20, b=5, a=3),
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon1=10, force.anon2=20, def.b=3,
                a=3, predict.b=5)

test.deprefix(expected=list(10, 20, b=6, a=3),
                KEEP=NULL, DROP="w1.,persp.,xlab.",
                PREFIX="predict.",
                force.anon1=10, force.anon2=20, def.b=3,
                a=3, b=6)

expect.err(try(test.deprefix(expected=NULL, KEEP=NULL, PREFIX="predict.", DROP="w1\\.")), "illegal character \"\\\" in DROP = \"w1\\.\"")

test.deprefix(expected=list(b="predict.b", d="def.d", a="a", c="predict.c", w1.xl="xlab2", e="predict.e"),
    PREFIX="predict.", DROP="w1.xlab$",
    a="a", b="b", c="c",
    w1.xlab="xlab1", # will be dropped (exact match)
    w1.xl="xlab2",   # will be kept (not an exact match)
    def.b="def.b", def.d="def.d",
    predict.b="predict.b", predict.c="predict.c", predict.e="predict.e")

# expect.err(try(plotmo:::deprefix(FNAME="test.deprefix", PREFIX="predict.", UPPER.CASE123=99,
#   def.a=1, predict.b=2, c=3)),
#   "uppercase argument names like \"UPPER.CASE123\" are not allowed for test.deprefix()")

test.expand.dotnames <- function(expected, PREFIX, FUNC=NULL,
                                 FNAME="test.expand.dotnames", FORMALS=NULL, ...)
{
    dots <- as.list(match.call(expand.dots=FALSE)$...)
    args <- plotmo:::expand.dotnames(dots, PREFIX, FUNC, FNAME, FORMALS)
    # can't use all.equal because it complains about named list versus unnamed list
    stopifnot(length(args) == length(expected))
    for(i in seq_len(length(expected))) {
        stopifnot(names(args)[i] == names(expected)[i])
        stopifnot(eval(args[[i]]) == expected[[i]])
    }
}
cat0("=== test expand.dotnames\n")

test.expand.dotnames(expected=list(x=9, persp.shade=3),
    "persp.", graphics:::persp.default, "persp.default", x=9, persp.sh=3)

test.expand.dotnames(expected=list(x=9, persp.shade=3, persp.nonesuch=4),
    "persp.", graphics:::persp.default, "persp.default", x=9, persp.sh=3, persp.nonesuch=4)

test.expand.dotnames(expected=list(x=9, persp.col=3),
    "persp.", graphics:::persp.default, "persp.default", x=9, persp.c=3)

# TODO not sure why this works as it does
test.expand.dotnames(expected=list(x=9, persp.x=3),
    "persp.", graphics:::persp.default, "persp.default", x=9, persp.x=3)

expect.err(try(test.expand.dotnames(expected=NULL,
    "persp.", graphics:::persp.default, "persp.default", x=9, persp.l=3)),
    "'l' matches both the 'ltheta' and 'lphi' arguments of persp.default()")

test.expand.dotnames(expected=list(x=9, plot.foo=3, plot.xlim=c(1,2)),
    "plot.", graphics:::plot.default, "plot.default", x=9, plot.foo=3, plot.xlim=c(1,2))

test.expand.dotnames(expected=list(x=9, plot.foo=3, plot.xlim=c(1,2)),
    "plot.", graphics:::plot.default, "plot.default", x=9, plot.foo=3, plot.xli=c(1,2))

expect.err(try(test.expand.dotnames(expected=NULL,
    "plot.", graphics:::plot.default, "plot.default", x=9, plot.foo=3, plot.xl=c(1,2))),
    "'xl' matches both the 'xlim' and 'xlab' arguments of plot.default()")

foo3 <- function(aaa=1, aa=2, bb=3, bba=4, cca=5, ccb=6, def=7)
    cat0("foo3: aaa=", aaa, " aa=", aa, ", bb=", bb, " bba=", bba,
         " cca=", cca, " ccb=", ccb, " def=", def, "\n")

# --- above tests again but using formals ---

# formal args for graphics:::persp.default (R version 3.2.0)
formals <- c( "x", "y", "z", "xlim", "zlim", "xlab", "ylab", "zlab",
    "main", "sub", "theta", "phi", "r", "d", "scale", "expand", "col",
    "border", "ltheta", "lphi", "shade", "box", "axes", "nticks",
    "ticktype")

test.expand.dotnames(expected=list(x=9, persp.shade=3),
    "persp.", graphics:::persp, "persp", FORMALS=formals, x=9, persp.sh=3)

test.expand.dotnames(expected=list(x=9, persp.shade=3, persp.nonesuch=4),
    "persp.", graphics:::persp, "persp", FORMALS=formals, x=9, persp.sh=3, persp.nonesuch=4)

test.expand.dotnames(expected=list(x=9, persp.col=3),
    "persp.", graphics:::persp, "persp", FORMALS=formals, x=9, persp.c=3)

# TODO not sure why this works as it does
test.expand.dotnames(expected=list(x=9, persp.x=3),
    "persp.", graphics:::persp, "persp", FORMALS=formals, x=9, persp.x=3)

expect.err(try(test.expand.dotnames(expected=NULL,
    "persp.", graphics:::persp, "persp", FORMALS=formals, x=9, persp.l=3)),
    "'l' matches both the 'ltheta' and 'lphi' arguments of persp()")

# done formals tests

test.expand.dotnames(expected=list(x=9, plot.foo=3, plot.xlim=c(1,2)),
    "plot.", graphics:::plot.default, "plot.default", x=9, plot.foo=3, plot.xlim=c(1,2))

test.expand.dotnames(expected=list(x=9, plot.foo=3, plot.xlim=c(1,2)),
    "plot.", graphics:::plot.default, "plot.default", x=9, plot.foo=3, plot.xli=c(1,2))

expect.err(try(test.expand.dotnames(expected=NULL,
    "plot.", graphics:::plot.default, "plot.default", x=9, plot.foo=3, plot.xl=c(1,2))),
    "'xl' matches both the 'xlim' and 'xlab' arguments of plot.default()")

test.expand.dotnames(expected=list(foo3.aa=99),
    "foo3.", foo3, "foo3", foo3.aa=99)
expect.err(try(plotmo:::call.plot(foo3, "foo3.", foo3.aa=99)), "Unnamed arguments are not allowed here\n       The argument's value is \"foo3.\"")
expect.err(try(plotmo:::call.plot(foo3, foo, foo3.aa=99)),
"Unnamed arguments are not allowed here\n       The argument's value is function.object")
expect.err(try(plotmo:::call.plot(foo3, NULL, foo3.aa=99)), "Unnamed arguments are not allowed here\n       The argument's value is NULL")
expect.err(try(plotmo:::call.plot(foo3, stop("stop was called"), foo3.aa=99)), "Unnamed arguments are not allowed here (argument ..1 is unnamed)")
expect.err(try(plotmo:::call.plot(foo3, cat("side effect\n"), foo3.aa=99)), "Unnamed arguments are not allowed here\n       The argument's value is NULL")
expect.err(try(plotmo:::call.plot(foo3, nonesuch1=1, nonesuch2, foo3.aa=99)), "Unnamed arguments are not allowed here (argument ..2 is unnamed)")
plotmo:::call.plot(foo3, PREFIX="foo3.", foo3.aa=99)

test.expand.dotnames(expected=list(foo3.aaa=99),
    "foo3.", foo3, "foo3", foo3.aaa=99)
plotmo:::call.plot(foo3, foo3.aaa=99)

expect.err(try(test.expand.dotnames(expected=list(foo3.aaa=99),
    "foo3.", foo3, "foo3", foo3.aa=88, foo3.aa=99)),
    "'foo3.aa' for foo3() is duplicated")

expect.err(try(test.expand.dotnames(expected=list(foo3.aaa=99),
    "foo3.", foo3, "foo3", foo3.a=88, foo3.aa=99)),
    "'a' matches both the 'aaa' and 'aa' arguments of foo3()")

expect.err(try(test.expand.dotnames(expected=list(foo3.aaa=99),
    "foo3.", foo3, "foo3", foo3.aaa=88, foo3.aaa=99)),
    "'foo3.aaa' for foo3() is duplicated")

test.expand.dotnames(expected=list(foo3.bbb=88, foo3.bba=99),
     "foo3.", foo3, "foo3", foo3.bbb=88, foo3.bba=99)
expect.err(try(plotmo:::call.plot(foo3, foo3.bbb=88, foo3.bba=99)),
     "unused argument (bbb = 88)")

# same as above but with TRACE (so don't use try in call.dots)
expect.err(try(plotmo:::call.plot(foo3, foo3.bbb=88, foo3.bba=99, TRACE=T)),
     "unused argument (bbb = 88)")

test.expand.dotnames(expected=list(foo3.bb=88),
     "foo3.", foo3, "foo3", foo3.bb=88)
plotmo:::call.plot(foo3, foo3.bb=88)

# test with FUNC=NULL

test.expand.dotnames(expected=list(foo3.aa=99),
    "foo3.", NULL, "foo3", foo3.aa=99)
plotmo:::call.plot(foo3, foo3.aa=99)

test.expand.dotnames(expected=list(foo3.aaa=99),
    "foo3.", NULL, "foo3", foo3.aaa=99)
plotmo:::call.plot(foo3, foo3.aaa=99)

expect.err(try(test.expand.dotnames(expected=list(foo3.aaa=99),
    "foo3.", NULL, "foo3", foo3.aa=88, foo3.aa=99)),
    "argument 'foo3.aa' for foo3() is duplicated")

test.expand.dotnames(expected=list(foo3.a=88, foo3.aa=99),
     "foo3.", NULL, "foo3", foo3.a=88, foo3.aa=99)
expect.err(try(plotmo:::call.plot(foo3, foo3.a=88, foo3.aa=99)),
     "'a' matches both the 'aaa' and 'aa' arguments of foo3()")

expect.err(try(test.expand.dotnames(expected=list(foo3.aaa=99),
     "foo3.", NULL, "foo3", foo3.aaa=88, foo3.aaa=99)),
     "argument 'foo3.aaa' for foo3() is duplicated")

test.expand.dotnames(expected=list(foo3.bbb=88, foo3.bba=99),
      "foo3.", NULL, "foo3", foo3.bbb=88, foo3.bba=99)
expect.err(try(plotmo:::call.plot(foo3, PREFIX="foo3.", foo3.bbb=88, foo3.bba=99)),
      "unused argument (bbb = 88)")

test.expand.dotnames(expected=list(foo3.bb=88),
      "foo3.", NULL, "foo3", foo3.bb=88)
plotmo:::call.plot(foo3, foo3.bb=88)

test.expand.dotnames(expected=list(foo3.bbx=88),
      "foo3.", NULL, "foo3", foo3.bbx=88)
expect.err(try(plotmo:::call.plot(foo3, foo3.bbx=88)),
      "unused argument (bbx = 88)")

test.expand.dotnames(expected=list(foo3.cc=77),
      "foo3.", NULL, "foo3", foo3.cc=77)
expect.err(try(plotmo:::call.plot(foo3, foo3.cc=77)),
      "'cc' matches both the 'cca' and 'ccb' arguments of foo3()")

# following two directly compare FUNC=NULL to FUNC=foo3
test.expand.dotnames(expected=list(foo3.cc=77),
           "foo3.", FUNC=NULL, "foo3", foo3.cc=77)
expect.err(try(test.expand.dotnames(expected=NULL,
           "foo3.", FUNC=foo3, "foo3", foo3.cc=77)),
           "'cc' matches both the 'cca' and 'ccb' arguments of foo3()")

test.expand.dotnames(expected=list(), "foo3.", foo3, "foo3", d=88, de=99)

expect.err(try(plotmo:::call.plot(graphics::plot, x=1:3, y=1:3, 99)),
"Unnamed arguments are not allowed here\n       The argument's value is 99\n       plotmo:::call.plot via try called call.dots(FUNC=plot, PREFIX=PREFIX, ...")

# test TRACE
print(plotmo:::deprefix(FUNC=nchar, PREFIX="foo3.", TRACE=TRUE, FNAME="nchar", allowN=1, b=2, foo3.c=3))
print(plotmo:::deprefix(FUNC=nchar, PREFIX="foo3.", TRACE=2,    allowN=1, b=2, foo3.c=3))
print(plotmo:::deprefix(FUNC=nchar, PREFIX="foo3.", TRACE=3,    allowN=1, b=2, foo3.c=3))

expect.err(try(plotmo:::call.plot(foo3, foo3.d=88, foo3.de=99)),
        "'foo3.d' and 'foo3.de' both match the 'def' argument of foo3()")

cat0("=== test stop.if.dots\n")

foo3 <- function(x=1, ...)  plotmo:::stop.if.dots(...)
foo3(1) # ok
expect.err(try(foo3(10, y=2)), "foo3: unrecognized argument 'y'")
expect.err(try(foo3(10, 99)), "foo3: unrecognized unnamed argument\n       The call was foo3(x=10, 99)")
expect.err(try(foo3(10, y=plot)), "foo3: unrecognized argument 'y'")
expect.err(try(foo3(10, plot)),
"foo3: unrecognized unnamed argument\n       The call was foo3(x=10, plot)")

expect.err(try(foo3(20, c(1,2,3), plot)),
"foo3: unrecognized unnamed argument\n       The call was foo3(x=20, c(1,2,3), plot)")

expect.err(try(foo3(20, c(1,2,3,4,5,6,7,8,9,10,11,12), plot)),
"foo3: unrecognized unnamed argument\n       The call was foo3(x=20, c(1,2,3,4,5,6,7,8,9,10,11,12), plot)")

# test that we don't crash because we eval the argument
expect.err(try(foo3(20, y=stop("stop was called"))), "foo3: unrecognized argument 'y'")
expect.err(try(foo3(20, stop("stop was called"))), "foo3: unrecognized unnamed argument")
expect.err(try(foo3(20, cat("side effect\n"))),
"foo3: unrecognized unnamed argument\n       The call was foo3(x=20, cat(")
foo2 <- function(...)  plotmo:::stop.if.dots(...)
foo2() # ok
expect.err(try(foo2(y=2)), "foo2: unrecognized argument 'y'")
expect.err(try(foo2(2)), "foo2: unrecognized unnamed argument\n       The call was foo2(2)")
expect.err(try(foo2(y=plot)), "foo2: unrecognized argument 'y'")
expect.err(try(foo2(plot)),
"foo2: unrecognized unnamed argument\n       The call was foo2(plot)")

foo2a <- function(funcarg, ...) funcarg(...)
expect.err(try(foo2a(function(x=1, ...) plotmo:::stop.if.dots(...), x=1, y=2)), "funcarg: unrecognized argument 'y'")

cat0("=== test warn.if.dots\n")

old.warn <- options("warn")
options(warn=2) # treat warnings as errors

foo3 <- function(x=1, ...)  plotmo:::warn.if.dots(...)
foo3(1) # ok
expect.err(try(foo3(1, y=2)), "foo3 ignored argument 'y'")
expect.err(try(foo3(1, 2)), "foo3 ignored unnamed argument\n       The call was foo3(x=1, 2)")
expect.err(try(foo3(1, y=plot)), "foo3 ignored argument 'y'")
# TODO would like to improve this error messsage
expect.err(try(foo3(1, plot)),
"(converted from warning) foo3 ignored unnamed argument\n       The call was foo3(x=1, plot)")
foo4 <- function(...)  plotmo:::warn.if.dots(...)
foo4() # ok
expect.err(try(foo4(y=2)), "foo4 ignored argument 'y'")
expect.err(try(foo4(2)), "foo4 ignored unnamed argument\n       The call was foo4(2)")
expect.err(try(foo4(y=plot)), "foo4 ignored argument 'y'")
expect.err(try(foo4(plot)),
"(converted from warning) foo4 ignored unnamed argument\n       The call was foo4(plot)")

options(warn=old.warn$warn)
foo3(1, nonesuch=12, nonesuch2=12, 999) # expect three warnings

cat0("=== test using sample functions that invoke call.dots\n")

x <- 1:10
y <- x * x
lmfit <- lm(y~x)
par(mfrow=c(3, 2))
par(oma=c(0, 0, 3, 0))

# plot1: simple example
# we choose to use predict() here rather than fitted() because nearly all
# models have a fitted() method, but many don't have a fitted() method.

plot1 <- function(object, ...)
{
    residuals <- residuals(object, ...)

    fitted <- predict(object, ...)

    plot(fitted, residuals, ...)
}
plot1(lmfit)
mtext("example plot functions using prefixed dots", outer=TRUE, font=2, line=1, cex=1)

# Following causes error in predict.lm().  The type argument meant for
# residuals() is also sent to predict.lm(), where it is rejected.

expect.err(try(plot1(lmfit, type="pearson")), "'arg' should be one of \"response\", \"terms\"")

# plot2: use prefixed args

plot2 <- function(object, ..., TRACE=2)
{
    resids <- plotmo:::call.dots(residuals, object=object, ..., TRACE=TRACE)

    fitted <- plotmo:::call.dots(predict, object=object, ..., TRACE=TRACE)

    plotmo:::call.plot(plot, x=fitted, y=resids, ..., TRACE=TRACE)
}
# we can now direct args using the prefixes "residuals.", "predict.", or "plot.")

plot2(lmfit, residuals.type="pearson")

# We can also use the usual plot arguments like ylab: call.dots drops
# them; call.plot recognizes them and passes them to lines().

plot2(lmfit, residuals.type="pearson", ylab="pearson residuals", main="plot2")

# plot3: further refinements
#   o namespace added to FUNC arg
#   o full name for plot.default
#   o force. and def. prefixes
#   o explicit xlab and ylab for call.plot
#   o unprefixed args are passed to residuals()

plot3 <- function(object, ..., TRACE=2)
{
    resids <- plotmo:::call.dots(stats::residuals,
                                 DROP="plotmo:::PLOTARGS,predict.,plot.", 
                                 TRACE=TRACE, force.object=object, ...)

    fitted <- plotmo:::call.dots(stats::predict,
                                 force.object=object, TRACE=TRACE, ...)

    plotmo:::call.plot(graphics::plot.default, force.x=fitted, force.y=resids,
                       def.xlab="fitted", def.ylab="residuals",
                       TRACE=TRACE, ...)
}
plot3(lmfit, type="pearson", main="plot3a") # type goes only to pearson, no prefix needed
plot3(lmfit, type="pearson", predict.type="response", main="plot3b")

if(!interactive()) {
    dev.off()        # finish postscript plot
    q(runLast=FALSE) # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
