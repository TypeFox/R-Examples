### R code from vignette source 'latex_fdt.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: tab
###################################################
library(fdth)

## Example 1: The simplest possible
t1 <- fdt(rnorm(n=1e3,
                mean=10,
                sd=2))

t1x <- latex.fdt(t1)

t1x


## Example 2
(t1x <- latex.fdt(t1,
                  replace.breaks=FALSE,
                  columns=c(1:2, 4, 6)))


## Example 3
t2 <- fdt(rnorm(n=1e3,
                mean=10,
                sd=2),
          right=TRUE)

t2x <- latex.fdt(t2,
                 algtable='\\centering',
                 caption='Frequency distribution table 1',
                 label='tbl-2',
                 pattern='%.1f')

t2x


## Example 4
t3 <- fdt(rnorm(n=1e3,
                mean=10,
                sd=2))

t3x <- latex.fdt(t3,
                 algtable='\\flushright',
                 caption='Frequency distribution table 2',
                 label='tbl-3',
                 pattern='%.1e')

t3x

## Example 5: ftd.multiple
t4 <- fdt(iris,
          by='Species')

t4x <- lapply(t4,
              function(x) latex.fdt(x,
                                    caption='Frequency distribution table'))

cat(unlist(t4x),
    sep='\n')

t5 <- fdt_cat(sample(LETTERS[1:3], 
                     replace=TRUE,
                     size=30))

t5x <- latex.fdt_cat(t5,
                     caption='Frequency distribution table 15',
                     label='tbl-4')

cat(unlist(t5x),
    sep='\n')

t6 <- fdt_cat(data.frame(c1=sample(LETTERS[1:3],
                                   replace=TRUE,
                                   size=10),
                         c2=sample(letters[4:5],
                                   replace=TRUE,
                                   size=10)))

captions <- c('Frequency distribution table 16',
             'Frequency distribution table 17')

t6x <- latex.fdt_cat.multiple(t6,
                              caption=captions,
                              algtable='\\flushleft',
                              label=c('tbl-5','tbl-6'))
cat(unlist(t6x),
    sep='\n')


