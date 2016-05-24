latex.bpca <- function(x,
                       round=4,
                       where='!tbp',
                       caption=NULL,
                       label=NULL,
                       size='', 
                       v.retained='Variance retained (%)',
                       v.partial='Partial',
                       v.accumulated='Accumulated',
                       eigenvalues='Eigenvalues',
                       eigenvectors='Eigenvectors',
                       algtable='\\flushleft',
                       alg1='l',
                       alg2='l',
                       algnumbers='r',
                       algheader='c',
                       algsubheader='c',
                       pc.label='PC',
                       bpca.label=NULL,
                       ft.variable='',
                       ft.components='',
                       hline1='\\hline',
                       hline2='\\hline',
                       footnote='')
{
  percent <- length(grep('%',
                         v.retained))

  if(percent)
    v.retained <- sub('%',
                      '\\\\\\\\%',
                      v.retained)     

  if(class(x)[1]=='list'){

    eigvec  <- lapply(lapply(x,
                             function(x) x$eigenvectors[,
                                                        c(x$number)]),
                      function(x)round(x,
                                       round))

    eigval  <- lapply(lapply(x,
                             function(x) x$eigenvalues[c(x$number)]),
                      function(x)round(x,
                                       round))

    xnumbers <- lapply(x,
                       function(x) x$number)

    lambdas <- list()
    for(i in 1:length(x)){
      lambdas[[i]] <- paste('$(', 
                            paste(paste('\\lambda_',
                                        xnumbers[[i]],
                                        sep=''),
                                  eigval[[i]],
                                  sep='='),
                            ')$', 
                            sep='')
    }

    lambdas <- unlist(lambdas)

    vret <- lapply(lapply(x,
                          function(x) x$eigenvalues[x$number[1]:x$number[length(x$number)]]^2 / sum(x$eigenvalues^2)),
                   function(x) round(x,
                                     round))

    if(percent)
      vret <-  lapply(lapply(lapply(x,
                                    function(x) x$eigenvalues[x$number[1]:x$number[length(x$number)]]^2 / sum(x$eigenvalues^2)),
                             function(x) x * 100),
                      function(x) round(x, 
                                        round))

    vacum <- lapply(vret,
                    cumsum)

    varis <- lapply(eigvec,
                    rownames)[[1]]

    if(ft.variable=='bold')
      varis <- paste('\\textbf{',
                     paste(varis,'}',
                           sep=''),
                     sep='')

    if(ft.variable=='italic')
      varis <- paste('\\textit{',
                     paste(varis,
                           '}',
                           sep=''),
                     sep='')  

    y <-cbind(c(varis,
                v.partial,
                v.accumulated),
              rbind(do.call(cbind,
                            eigvec),
                    unlist(vret),
                    unlist(vacum)))

    comp <- sub('PC',
                pc.label,
                colnames(y)[-1])   

    if(is.null(bpca.label))
      bpca.label <- paste('bp',
                          1:length(x),
                          sep='-')

    meansheader <- paste('&',
                         paste('&',
                               sub('algsubheader',
                                   algsubheader,
                                   paste(paste(paste(paste(rep('\\multicolumn{',
                                                               length(x)),
                                                           unlist(lapply(xnumbers,
                                                                         length)),
                                                           sep=''),
                                                     '}{algsubheader}{',
                                                     sep='')   ,
                                               bpca.label,
                                               sep=''),
                                         '}',
                                         sep='')),
                               sep='',
                               collapse=''),
                         '\\tabularnewline',
                         sep='')
  } else {   

    eigvec <- round(x$eigenvectors[, c(x$number)],
                    round)

    eigval <- round(x$eigenvalues[c(x$number)],
                    round)

    lambdas <- paste('$(',
                     '\\lambda_',
                     x$number,
                     '=',
                     eigval,
                     '$)',
                     sep='')

    comp <- sub('PC',
                pc.label,
                colnames(eigvec))

    vret <- round(x$eigenvalues[x$number[1]:x$number[length(x$number)]]^2 / sum(x$eigenvalues^2),
                  round)

    if(percent) 
      vret <- round(x$eigenvalues[x$number[1]:x$number[length(x$number)]]^2 / sum(x$eigenvalues^2) * 100,
                    round)

    vacum <- cumsum(vret)     

    varis <- rownames(eigvec)

    if(ft.variable=='bold')
      varis <- paste('\\textbf{',
                     paste(varis,
                           '}',
                           sep=''),
                     sep='')

    if(ft.variable=='italic')
      varis <- paste('\\textit{',
                     paste(varis,
                           '}',
                           sep=''),
                     sep='') 

    y <- cbind(c(varis,
                 v.partial,
                 v.accumulated),
               rbind(eigvec,
                     vret,
                     vacum))

  }

  y <- t(apply(y,
               1,
               function(x) paste('&',
                                 x,
                                 sep='')))

  dy <- dim(y)[2]+1

  cline <- sub('dy',
               dy,
               '\\cline{3-dy}')  

  rowname <- sub('numb1',
                 length(varis),
                 sub('v.retained',
                     v.retained,
                     sub('eigenvectors',
                         eigenvectors,
                         c('\\multirow{numb1}{*}{eigenvectors}',
                           rep('',
                               length(varis) - 1),
                           '\\hline\\multirow{2}{*}{v.retained}',
                           ''))))

  y <- cbind(rowname,y)

  rownames(y) <- NULL 
  colnames(y) <- NULL    

  algnumb <- rep(algnumbers, 
                 dim(y)[2]-2)   

  y <- paste(apply(y,
                   1,
                   paste,
                   collapse=''),
             '\\tabularnewline',
             sep='')                               

  smallheader <- paste(comp,
                       lambdas,
                       sep='\\quad')  

  if(ft.components=='bold')
    smallheader <- paste('\\textbf{',
                         smallheader,
                         '}',
                         sep='')

  if(ft.components=='italic')
    smallheader <- paste('\\textbf{',
                         smallheader,
                         '}',
                         sep='')      

  res <- list(start=NULL,
              algtable=algtable,
              size=size,
              caption=NULL,
              label=NULL,
              begintabular=NULL,
              hline1=hline1,
              greaterheader=NULL,
              cline=cline,
              meansheader= if(class(x)[1]=='list') meansheader else '',
              cline=if(class(x)[1]=='list') cline else '',
              smallheader=paste(c('&',
                                  paste('&',
                                        smallheader,
                                        sep=''),
                                  '\\tabularnewline'),
                                collapse=''),
              hline2=hline2,
              bpcaobject=y,
              endtabular=NULL,
              footnote=footnote,
              end=NULL)                    

  res$start <- sub('where',
                   where,
                   '\\begin{table}[where]')

#  if(!is.null(algtable))
#    res$algtable <- algtable    

  if(!is.null(caption))
    res$caption <- paste('\\caption{', 
                         caption,
                         '}',
                         sep='')

  if(!is.null(label))
    res$label <- paste('\\label{', 
                       label,
                       '}',
                       sep='')

  algnumber <- paste(algnumb,
                     collapse='')      

  begintabular <-sub('alg1',
                     alg1,
                     sub('alg2',
                         alg2,
                         sub('algnumb',
                             algnumber,
                             '\\begin{threeparttable} \n \\begin{tabular}{alg1alg2algnumb}'))) 

  res$begintabular <- begintabular

  if(class(x)[1]=='list')
    xnumbers <- sum(unlist(lapply(xnumbers,
                                  length)))
  else
    xnumbers <- length(x$number)

  res$greaterheader <- sub('eigenvalues',
                           eigenvalues,
                           sub('xnumbers',
                               xnumbers,
                               sub('algheader',
                                   algheader,
                                   '&&\\multicolumn{xnumbers}{algheader}{eigenvalues}\\tabularnewline')))

  res$endtabular <- '\\hline \n \\end{tabular}'

  res$footnote <- sub('footnote',
                      footnote,
                      '\\begin{tablenotes}[para,flushleft] \n \\item footnote \n \\end{tablenotes} \n \\end{threeparttable}')

  res$end <- paste('\n',
                   '\\end{table}',
                   sep=' ')

  class(res) <- c('latex.bpca', 'list')

  return(res)
}
