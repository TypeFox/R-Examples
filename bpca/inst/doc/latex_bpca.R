### R code from vignette source 'latex_bpca.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: tab
###################################################
library(bpca)

##Example 1: The simplest possible
bp <- bpca(iris[-5],
           d=1:2)
lbp <- latex.bpca(bp)
lbp
summary(lbp)

##Example 2: With caption 
bp2 <- bpca(gabriel1971) 

latex.bpca(bp2,
           caption='Biplot gabriel1971 (alignment = centering)',
           label='example_2',
           algtable='\\centering')

## Example 3: Almost worked
latex.bpca(bp,
           caption='Biplot da base de dados iris (aligment = flushright)',
           label='example_3',
           v.partial='Parcial',
           v.accumulated='Acumulada',
           eigenvalues='Autovalores',
           eigenvectors='Autovetores',
           v.retained='Var. retida',
           algtable='\\flushright')

## Example 4: Changing the alignment of the column first 
bp3 <- bpca(gge2003,
            d=1:3)

latex.bpca(bp3,
           caption='Biplot gge2003 (changing the alignment of the column first)',    
           alg1='r')

## Example 5: Changing the alignment of the second column 
latex.bpca(bp3,
           caption='Biplot gge2003 (changing the alignment of the second column)',    
           #alg2='>{\\\\raggedright}p{0.1cm}')
           alg2='r')

## Example 6: Changing the column alignment with numbers
latex.bpca(bp3,
           caption='Biplot gge2003 (changing the column alignment with numbers)',    
           #algnumbers='>{\\\\raggedleft}p{2.2cm}')
           alg1='r',
           alg2='r',
           algnumbers='l')

## Example 7: Changing the header alignment 
latex.bpca(bp3,
           caption='Biplot gge2003 (changing the header alignment)',    
           algheader='r')

## Example 8: Two decimals
latex.bpca(bp3,
           round=2,
           caption='Biplot gge2003 (two decimals)', 
           algnumbers='>{\\\\centering}p{2.2cm}',
           pc.label='Principal Component-')

## Example 9: With bold in the header, subheader and variables
latex.bpca(bp3,
           round=2,
           caption='Biplot gge2003 (bold in the header, subheader and variables)',
           eigenvalues='\\\\textbf{Eigenvalues}',
           eigenvectors='\\\\textbf{Eigenvectors}',
           v.retained='\\\\textbf{Variance retained (%)}',
           v.partial='\\textbf{Partial}',
           v.accumulated='\\textbf{Accumulated}',
           ft.variable='bold',
           ft.components='bold')  

## Example 10: Changing the font
latex.bpca(bp3,
           round=2,
           caption='Biplot gge2003 (changing the font)', 
           pc.label='\\\\textbf{P. Component-}',           
           size='\\scriptsize')                     

## Example 11: Italic in the variables names
latex.bpca(bp2,
           round=2,
           caption='Biplot gabriel1971 (italic in the variables names)',     
           pc.label='\\\\textbf{Principal Component-}',
           algnumbers='>{\\\\centering}p{2.5cm}',
           ft.variable='italic') 

## Example 12: With footnote
latex.bpca(bp2,   
           round=3,
           caption='Biplot gabriel1971 (with footnote)',     
           footnote='$^1$\\\\scriptsize Example with footnote')

## Example 13: With others principal components
bp4 <- bpca(gabriel1971,
            d=2:4)

latex.bpca(bp4,
           round=5,
           caption='Biplot gabriel1971 (with principal components 2, 3 and 4)')

## Example 14: More than one bpca object
data(marina)

y_2007 <- subset(marina,
                 year==2007)

y_2008 <- subset(marina,
                 year==2008)

y_2009 <- subset(marina,
                 year==2009)   

bp_2007 <- bpca(y_2007[,-c(1:2)],
                d=1:3)

bp_2008 <- bpca(y_2008[,-c(1:2)],
                d=1:2) 

bp_2009 <- bpca(y_2009[,-c(1:2)],
                d=1:2)      

latex.bpca(list(bp_2007,
                bp_2008),
           round=2,
           caption='Biplot Marina (more than one bpca)',
           bpca.label=paste('Year',
                            2007:2008,
                            sep='-'),
           size='\\scriptsize')

## Example 15: With two lines in the table
latex.bpca(list(bp_2007,
                bp_2008,
                bp_2009),
           round=3,
           caption='Biplot Marina (with two lines)',
           label='example_15',
           bpca.label=c('2007','2008','2009'),
           hline1='\\hline \\hline',
           hline2='\\hline \\hline',
           algnumbers='>{\\\\raggedleft}p{1.3cm}',
           size='\\scriptsize',
           footnote='Note: F - Movie; D - Documentary; DH - Documentary directed by men; DF - Documentary directed by women.')  


