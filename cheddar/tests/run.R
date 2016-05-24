# Runs all of cheddar's tests. You must cd cheddar/tests before running 
# Rscript run_all.R
options(warn=2)
library(cheddar)

# Data for test plans
data(Benguela, BroadstoneStream, ChesapeakeBay, SkipwithPond, TL84, TL86, 
     YthanEstuary, pHWebs, Millstream)

# One species. No properties. No trophic links.
c1 <- Community(nodes=data.frame(node='S'), properties=list(title='c1'))

# One cannibalistic species
c2 <- Community(nodes=data.frame(node='S'), 
                trophic.links=cbind(resource='S', consumer='S'), 
                properties=list(title='c2'))

# Resource-consummer
c3 <- Community(nodes=data.frame(node=c('R','C')), 
                trophic.links=cbind(resource='R', consumer='C'), 
                properties=list(title='c3'))

# Three-species chain
c4 <- Community(nodes=data.frame(node=c('R','C','P')), 
                trophic.links=cbind(resource=c('R','C'), 
                                    consumer=c('C','P')), 
                properties=list(title='c4'))

# IGP
c5 <- Community(nodes=data.frame(node=c('R','C','O')), 
                trophic.links=cbind(resource=c('R','R','C'), 
                                    consumer=c('C','O','O')), 
                properties=list(title='c5'))

# Three species chain with M, N, taxonomy and trophic link properties
c6 <- Community(nodes=data.frame(node=c('R','C','P'), 
                                 M=c(1.5, 5, 100), 
                                 N=c(1000, 10, 5.5), 
                                 order=c('Order 1', 'Order 2', 'Order 2'), 
                                 family=c('Family 1', 'Family 2', 'Family 3'), 
                                 stringsAsFactors=FALSE), 
                trophic.links=data.frame(resource=c('R','C'), 
                                         consumer=c('C','P'),  
                                         link.evidence=c('Inferred', 'Known'), 
                                         link.strength=c(0.5, 0.2), 
                                         stringsAsFactors=FALSE), 
                properties=list(title='c6', M.units='g', 
                                N.units='m^-3'))

# A community with three nodes:
#   A is a cannibal with no resources other than itself
#   B consumes A
#   C is a cannibal that consumes B
#   D is a cannibal and has no other resources or consumers
#   E has no resources or consumers
# This community is biologically silly but mathematically interesting for 
# testing node status: basal, intermediate, top-level etc
c7 <- Community(properties=list(title='Test'),
                nodes=data.frame(node=c('A', 'B', 'C', 'D', 'E')), 
                trophic.links=data.frame(resource=c('A', 'A', 'B', 'C', 'D'), 
                                         consumer=c('A', 'B', 'C', 'C', 'D')))

# 10 species with different categories
c8 <- Community(nodes=data.frame(node=paste('Species', 1:8), 
                                 M=c(10,  9,  5, 0.1, 0.2, 0.3, 1, 2),
                                 N=c(10, 1,  21, 900, 500, 100, 1, 5),
                                 category=c(rep('vert.endo', 3), 
                                            rep('', 3),
                                            rep('vert.ecto', 2))),
                properties=list(title='c7', M.units='g', N.units='m^-3'))

AssertEqual <- function(a, b, ...)
{
    res <- all.equal(a, b, ...)
    if(!isTRUE(res))
    {
        stop(paste(res, collapse="\n"))
    }
}

AssertTrue <- function(v)
{
    AssertEqual(TRUE, v)
}

AssertFalse <- function(v)
{
    AssertEqual(FALSE, v)
}

AssertNull <- function(v)
{
    AssertEqual(NULL, v)
}

AssertRaises <- function(ex)
{
    # A function that expects an exception to be raise when ex is evalutated
    res <- tryCatch(eval(ex), error=function(e) e)
    if(!"error" %in% class(res))
    {
        stop('Did not raise error\n')
    }
}

RunTests <- function(tests)
{
    # tests should be a vector of function names
    options(error=function()
            {
                traceback(3)
                stop()
            }, 
            stringsAsFactors=FALSE)

    if(0==length(tests))
    {
        stop('No tests to run! Is the working directory cheddar/tests ?')
    }
    else
    {
        for(test in tests)
        {
            cat(paste('Running [', test, ']\n', sep=''))
            do.call(test, args=list())
       }
    }
}

# Source all files in this dir except this one
files <- list.files(getwd(), pattern='*R$')
files <- setdiff(files, 'run.R')
junk <- sapply(file.path(getwd(), files), source)
tests <- commandArgs(trailingOnly=TRUE)
if(0==length(tests))
{
    tests <- ls(pattern=glob2rx('^Test*'))
}

RunTests(tests)
