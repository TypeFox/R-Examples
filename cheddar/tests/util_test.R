TestPredationMatrixToLinks <- function()
{
    # Square matrices
    n <- paste('S',1:10)   # Names for testing

    m <- matrix(0, ncol=10, nrow=10, dimnames=list(n,n))
    AssertEqual(0, nrow(PredationMatrixToLinks(m)))

    m <- matrix(0, ncol=10, nrow=10, dimnames=list(n,n))
    m[1,1] <- 1
    AssertEqual(data.frame(resource='S 1', consumer='S 1', stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    m <- matrix(0, ncol=10, nrow=10, dimnames=list(n,n))
    m[1,1] <- m[10,10] <- 1
    AssertEqual(data.frame(resource=c('S 1','S 10'),
                           consumer=c('S 1','S 10'), stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    m <- matrix(0, ncol=10, nrow=10, dimnames=list(n,n))
    m[1,1] <- m[10,10] <- m[1,10] <- 1
    AssertEqual(data.frame(resource=c('S 1','S 1','S 10'), 
                           consumer=c('S 1','S 10','S 10'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    t1 <- PredationMatrixToLinks(PredationMatrix(TL84))
    t2 <- TLPS(TL84)[,c('resource', 'consumer')]
    AssertEqual(t1, t2)

    # Logical values
    m <- matrix(0, ncol=10, nrow=10, dimnames=list(n,n))
    m[1,1] <- m[10,10] <- TRUE
    AssertEqual(data.frame(resource=c('S 1','S 10'),
                           consumer=c('S 1','S 10'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    # NA
    m <- matrix(NA, ncol=10, nrow=10, dimnames=list(n,n))
    m[1,1] <- m[10,10] <- TRUE
    AssertEqual(data.frame(resource=c('S 1','S 10'),
                           consumer=c('S 1','S 10'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    # Values other than 1
    m <- matrix(NA, ncol=10, nrow=10, dimnames=list(n,n))
    m[1,1] <- m[10,10] <- -10
    AssertEqual(data.frame(resource=c('S 1','S 10'),
                           consumer=c('S 1','S 10'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    # A non-square matrix
    m <- matrix(0, ncol=3, nrow=2)
    colnames(m) <- letters[1:3]
    rownames(m) <- letters[4:5]
    m[1,1] <- 0.3
    m[2,1] <- 0.7
    m[1,2] <- 1
    m[1,3] <- 0.1
    m[2,3] <- 0.9
    AssertEqual(data.frame(resource=c('d','e','d','d','e'), 
                           consumer=c('a','a','b','c','c'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    # The same non-square matrix with the link property extracted
    AssertEqual(data.frame(resource=c('d','e','d','d','e'), 
                           consumer=c('a','a','b','c','c'), 
                           diet.fraction=c(0.3,0.7,1,0.1,0.9), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m, link.property='diet.fraction'))

    # A data.frame as input
    m <- data.frame(a=c(1,1), b=c(1,0), c=c(1,1), row.names=c('d', 'e'))
    AssertEqual(data.frame(resource=c('d','e','d','d','e'), 
                           consumer=c('a','a','b','c','c'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    # A data.frame with a link property
    m <- data.frame(a=c(0.3,0.7), b=c(1,0), c=c(0.1,0.9), row.names=c('d', 'e'))
    AssertEqual(data.frame(resource=c('d','e','d','d','e'), 
                           consumer=c('a','a','b','c','c'), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m))

    # The same data.frame with the link property extracted
    AssertEqual(data.frame(resource=c('d','e','d','d','e'), 
                           consumer=c('a','a','b','c','c'), 
                           diet.fraction=c(0.3,0.7,1,0.1,0.9), 
                           stringsAsFactors=FALSE), 
                PredationMatrixToLinks(m, link.property='diet.fraction'))

    # Not a matrix
    AssertRaises(PredationMatrixToLinks(NA))
    AssertRaises(PredationMatrixToLinks(1:10))
    AssertRaises(PredationMatrixToLinks(NA))

    # No names
    AssertRaises(PredationMatrixToLinks(matrix(0, ncol=10, nrow=10)))
}

TestStripWhiteSpace <- function()
{
    AssertEqual('', cheddar:::.StripWhitespace(''))
    AssertEqual('', cheddar:::.StripWhitespace(' '))
    AssertEqual('', cheddar:::.StripWhitespace('   '))
    AssertEqual('a', cheddar:::.StripWhitespace('a'))
    AssertEqual('a', cheddar:::.StripWhitespace('a '))
    AssertEqual('a', cheddar:::.StripWhitespace('a  '))
    AssertEqual('a', cheddar:::.StripWhitespace(' a'))
    AssertEqual('a', cheddar:::.StripWhitespace('  a'))
    AssertEqual('a', cheddar:::.StripWhitespace(' a '))
    AssertEqual('a', cheddar:::.StripWhitespace('  a  '))
    AssertEqual('a b c', cheddar:::.StripWhitespace('a b c'))
    AssertEqual('a b c', cheddar:::.StripWhitespace(' a b c'))
    AssertEqual('a b c', cheddar:::.StripWhitespace(' a b c '))
    AssertEqual('\\.[]a b c.$^-+.;/"',  
               cheddar:::.StripWhitespace(' \\.[]a b c.$^-+.;/" '))
}

TestFormatLM <- function()
{
    # Values to 5 dp, no r squared
    models <- NvMLinearRegressions(TL84)
    res <- sapply(models, FormatLM, dp=5, r.squared=FALSE)
    expected <- expression(all = "y" == "-2.68628" ~ structure("-", .Names = "x") ~ "0.82711" * "x" * ""*"" * "", 
        producer = "y" == "2.55834" ~ structure("-", .Names = "x") ~ "0.40715" * "x" * "" * "" * "", 
        invertebrate = "y" == "1.46561" ~ structure("-", .Names = "x") ~ "0.32432" * "x" * "" * "" * "", 
        vert.ecto = "y" == "-34.66097" ~ structure("-", .Names = "x") ~ "11.62787" * "x" * "" * "" * "")
    AssertEqual(expected, res)

    # Values to 2 dp, lots of info
    models <- NvMLinearRegressions(TL84)
    res <- sapply(models, FormatLM, r=TRUE, slope.95.ci=TRUE, 
                  ci.plus.minus.style=TRUE)
    expected <- expression(all = "y" == "-2.69" ~ structure("-", .Names = "x") ~ "0.83" * "x" * ("" %+-% 0.1 ~ "(95% CI, n=56)") * ("," ~ r == "-0.92") * ("," ~ r^2 == "0.84"), 
        producer = "y" == "2.56" ~ structure("-", .Names = "x") ~ "0.41" * "x" * ("" %+-% 0.23 ~ "(95% CI, n=31)") * ("," ~ r == "-0.56") * ("," ~ r^2 == "0.32"), 
        invertebrate = "y" == "1.47" ~ structure("-", .Names = "x") ~ "0.32" * "x" * ("" %+-% 0.24 ~ "(95% CI, n=22)") * ("," ~ r == "-0.54") * ("," ~ r^2 == "0.29"), 
        vert.ecto = "y" == "-34.66" ~ structure("-", .Names = "x") ~ "11.63" * "x" * ("" %+-% 63.36 ~ "(95% CI, n=3)") * ("," ~ r == "-0.92") * ("," ~ r^2 == "0.84"))
    AssertEqual(expected, res)

    # names other than x and y
    m <- lm(Log10N(TL84) ~ Log10M(TL84))
    res <- FormatLM(m, r=TRUE, slope.95.ci=TRUE, ci.plus.minus.style=TRUE)
    expected <- expression("Log10N(TL84)" == "-2.69" ~ structure("-", .Names = "Log10M(TL84)")~ "0.83" * "Log10M(TL84)" * ("" %+-% 0.1 ~ "(95% CI, n=56)") * ("," ~ r == "-0.92") * ("," ~ r^2 == "0.84"))
    AssertEqual(expected, res)
}

