SensSpec.demo <- function(sens, spec, prev, n=100000, step=11) {
    mat <- matrix(NA, ncol=4, nrow=4)
    dimnames(mat) <- list( Test=c('Positive','Negative','','Total'),
                           Disease=c('   Yes','    No','  ','  Total') )
    pplines <- c(' ',
                 'PPV =',
                 'NPV =')
    mat[4,4] <- n
    if(step>1){
        mat[4,1] <- round(n*prev)
    }
    if(step>2){
        mat[4,2] <- n-mat[4,1]
    }
    if(step>3){
        mat[1,1] <- round( sens*mat[4,1] )
    }
    if(step>4){
        mat[2,1] <- mat[4,1] - mat[1,1]
    }
    if(step>5){
        mat[2,2] <- round( spec*mat[4,2] )
    }
    if(step>6){
        mat[1,2] <- mat[4,2]-mat[2,2]
    }
    if(step>7){
        mat[1,4] <- mat[1,1]+mat[1,2]
    }
    if(step>8){
        mat[2,4] <- mat[2,1]+mat[2,2]
    }
    if(step>9){
        pplines[2] <- paste( 'PPV = ', mat[1,1], '/', mat[1,4],
                             ' = ', round(mat[1,1]/mat[1,4], 4), sep='')
    }
    if(step>10){
        pplines[3] <- paste( 'NPV = ',mat[2,2], '/', mat[2,4],
                             ' = ', round(mat[2,2]/mat[2,4], 4), sep='')
    }

    print(mat, na.print='')
    cat(paste(pplines, collapse='\n'),"\n\n")

    invisible(mat[-3,-3])
}

