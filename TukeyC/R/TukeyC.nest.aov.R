##
## S3 method to 'aov' object
##

TukeyC.nest.aov <- function (x,
                             which, 
                             fl1, 
                             fl2=0,
                             sig.level=.05,
                             round=2,
                             dispersion=c('mm', 's', 'se'), ...)
{
  mt <- model.tables(x,
                     "means")           # summary tables for model fits

  if(is.null(mt$n))
    stop("No factors in the fitted model!")

  nfa <- names(mt$tables)[-1] # nomes dos fatores

  nf1 <- unlist(strsplit(which,
                         split=':'))[1] # nome do primeiro fator do which

  nf2 <- unlist(strsplit(which,
                         split=':'))[2] # nome do segundo fator do which

  nf3 <- unlist(strsplit(which,
                         split=':'))[3] # nome do terceiro fator do which

  MSE <- deviance(x)/df.residual(x)

  group <- NULL

  group2 <- NULL

  if(fl2 == 0){
    # MODELO SOMENTE COM DOIS FATORES
    if(length(nfa[grep('[[:punct:]]',
                       nfa)]) == 1 && # condição necessária para certificar-se que no modelo há somente uma interação!
       which != nfa[grep('[[:punct:]]',
                         nfa)]) {
      whichn <- paste(nf2,
                      nf1,
                      sep=':') # necessário pois no modelo original estamos em uma ordem inversa!

      r <- mt$n[names(mt$tables)][[whichn]] # groups and its number of replicates

      bal <- ifelse(length(r) == 1,
                    TRUE,
                    FALSE) # is (or not) balanced          

      m <- as.vector(mt$tables[whichn][[whichn]][,fl1]) # pegando as médias de interesse

      which1 <- names(dimnames(mt$tables[whichn][[whichn]]))[2] # corresponde ao primeiro fator do seu 'which'

      which2 <- names(dimnames(mt$tables[whichn][[whichn]]))[1] # corresponde ao segundo fator do seu 'which'

      m.inf <- m.inf.2a(x,
                        which1,
                        which2,
                        dispersion)
      
      f1 <- levels(x$model[,which2]) # correspondem aos fatores que se quer comparar!

      f2 <- levels(x$model[,which1])[fl1] # corresponde ao fator onde se está fazendo o desdobramento!

      m.inf <- subset(m.inf, 
                      group == f2)[,2]

      rownames(m.inf) <- paste(f2,
                               f1,
                               sep='/') 

      ord <- order(m,
                   decreasing=TRUE)

      m.inf <- cbind(m.inf[,1][ord],
                     m.inf[,2][ord],
                     m.inf[,3][ord])

      tab <- mt$tables[whichn][[whichn]]  
    } else if(length(nfa[grep('[[:punct:]]',
                              nfa)]) == 1 &&
              which == nfa[grep('[[:punct:]]',
                                nfa)]) { 

      r <- mt$n[names(mt$tables)][[which]] # groups and its number of replicates

      bal <- ifelse(length(r) == 1,
                    TRUE,
                    FALSE) # is (or not) balanced          

      m <- as.vector(mt$tables[which][[which]][fl1,]) # pegando as médias de interesse

      which1 <- names(dimnames(mt$tables[which][[which]]))[1] # corresponde ao primeiro fator do seu 'which'

      which2 <- names(dimnames(mt$tables[which][[which]]))[2] # corresponde ao segundo fator do seu 'which'

      m.inf <- m.inf.2a(x,
                        which1,
                        which2,
                        dispersion)

      f1 <- levels(x$model[,which2]) # correspondem aos fatores que se quer comparar!

      f2 <- levels(x$model[,which1])[fl1] # corresponde ao fator onde se está fazendo o desdobramento!

      m.inf <- subset(m.inf, 
                      group == f2)[,2]

      rownames(m.inf) <- paste(f2,
                               f1,
                               sep='/') 

      ord <- order(m,
                   decreasing=TRUE)

      m.inf <- cbind(m.inf[,1][ord],
                     m.inf[,2][ord],
                     m.inf[,3][ord])

      tab <- mt$tables[which][[which]] 
    } else
      # MODELO COM TRES FATORES
      if(length(nfa[grep('[[:punct:]]',
                         nfa)]) != 1 &&
         which == nfa[grep('[[:punct:]]',
                           nfa)][1] |
         which == nfa[grep('[[:punct:]]',
                           nfa)][2] |  
         which == nfa[grep('[[:punct:]]',
                           nfa)][3] ) {

        r <- mt$n[names(mt$tables)][[which]] # groups and its number of replicates

        bal <- ifelse(length(r) == 1,
                      TRUE,
                      FALSE) # is (or not) balanced          

        m <- as.vector(mt$tables[which][[which]][fl1,])

        which1 <- names(dimnames(mt$tables[which][[which]]))[1] 

        which2 <- names(dimnames(mt$tables[which][[which]]))[2]

        m.inf <- m.inf.2a(x,
                          which1,
                          which2,
                          dispersion)

        f1 <- levels(x$model[,which2])

        f2 <- levels(x$model[,which1])[fl1] 

        m.inf <- subset(m.inf, 
                        group == f2)[,2]

        rownames(m.inf) <- paste(f2,
                                 f1,
                                 sep='/') 

        ord <- order(m, decreasing=TRUE)

        m.inf <- cbind(m.inf[,1][ord],
                       m.inf[,2][ord],
                       m.inf[,3][ord]) 

        tab <- mt$tables[which][[which]]
      } else if(length(nfa[grep('[[:punct:]]',
                                nfa)]) != 1 &&
                which != nfa[grep('[[:punct:]]',
                                  nfa)][1]) {

        t1 <- unlist(strsplit(which,
                              split=':'))[1] 

        t2 <- unlist(strsplit(which,
                              split=':'))[2]

        whichn <- paste(t2,
                        t1,
                        sep=':')

        r <- mt$n[names(mt$tables)][[whichn]] # groups and its number of replicates

        bal <- ifelse(length(r) == 1,
                      TRUE,
                      FALSE) # is (or not) balanced          

        m <- as.vector(mt$tables[whichn][[whichn]][,fl1])

        which1 <- names(dimnames(mt$tables[whichn][[whichn]]))[2] 

        which2 <- names(dimnames(mt$tables[whichn][[whichn]]))[1]

        m.inf <- m.inf.2a(x,
                          which1,
                          which2,
                          dispersion)

        f1 <- levels(x$model[,which2])

        f2 <- levels(x$model[,which1])[fl1] 

        m.inf <- subset(m.inf, 
                        group == f2)[,2]

        rownames(m.inf) <- paste(f2,
                                 f1,
                                 sep='/')

        ord <- order(m,
                     decreasing=TRUE)

        m.inf <- cbind(m.inf[,1][ord],
                       m.inf[,2][ord],
                       m.inf[,3][ord])

        tab <- mt$tables[whichn][[whichn]] 
      }
  } else if(fl2 != 0) {
    natri <- nfa[grep('[[:punct:]].{,100}[[:punct:]]',
                      nfa)]

      nt1 <- unlist(strsplit(natri,
                             split=':'))[1]

      nt2 <- unlist(strsplit(natri,
                             split=':'))[2]

      nt3 <- unlist(strsplit(natri,
                             split=':'))[3]
      if(which == natri){
        r <- mt$n[names(mt$tables)][[which]] # groups and its number of replicates

        bal <- ifelse(length(r) == 1,
                      TRUE,
                      FALSE) # is (or not) balanced          

        which1 <- names(dimnames(mt$tables[which][[which]]))[1]

        which2 <- names(dimnames(mt$tables[which][[which]]))[2]

        which3 <- names(dimnames(mt$tables[which][[which]]))[3]

        m.inf <- m.inf.3a(x,
                          which1=which3,
                          which2=which2,
                          which3=which1,
                          dispersion)

        f1 <- levels(x$model[,which3])

        f2 <- levels(x$model[,which2])[fl2] 

        f3 <- levels(x$model[,which1])[fl1]

        m.inf <- subset(m.inf, 
                        group == f2 & 
                        group2 == f3)[,3]

        rownames(m.inf) <- paste(f3,
                                 f2,
                                 f1,
                                 sep='/')   

        ord <- order(as.vector(m.inf[,1]),
                     decreasing=TRUE)

        m.inf <- cbind(m.inf[,1][ord],
                       m.inf[,2][ord],
                       m.inf[,3][ord])

        tab <- mt$tables[which][[which]]
      } else if(which != natri){
        if(nf1 == nt2 && nf2 == nt3) {

          r <- mt$n[names(mt$tables)][[natri]] # groups and its number of replicates

          bal <- ifelse(length(r) == 1,
                        TRUE,
                        FALSE) # is (or not) balanced          

          m <- as.vector(mt$tables[natri][[natri]][,fl1,fl2])

          which1 <- names(dimnames(mt$tables[natri][[natri]]))[1] 

          which2 <- names(dimnames(mt$tables[natri][[natri]]))[2]

          which3 <- names(dimnames(mt$tables[natri][[natri]]))[3]

          m.inf <- m.inf.3a(x,
                            which1=which1,
                            which2=which3,
                            which3=which2,
                            dispersion)

          f1 <- levels(x$model[,which1])

          f2 <- levels(x$model[,which3])[fl2] 

          f3 <- levels(x$model[,which2])[fl1] 

          m.inf <- subset(m.inf, 
                          group == f2 &
                          group2 == f3)[,3]

          rownames(m.inf) <- paste(f3,
                                   f2,
                                   f1,
                                   sep='/')

          ord <- order(m,
                       decreasing=TRUE)

          m.inf <- cbind(m.inf[,1][ord],
                         m.inf[,2][ord],
                         m.inf[,3][ord]) 

          tab <- mt$tables[natri][[natri]]   
        } else if(nf1 == nt2 && nf2 == nt1) {
          r <- mt$n[names(mt$tables)][[natri]] # groups and its number of replicates

          bal <- ifelse(length(r) == 1,
                        TRUE, 
                        FALSE) # is (or not) balanced          

          which1 <- names(dimnames(mt$tables[natri][[natri]]))[1] 

          which2 <- names(dimnames(mt$tables[natri][[natri]]))[2]

          which3 <- names(dimnames(mt$tables[natri][[natri]]))[3]

          m.inf <- m.inf.3a(x,
                            which1=which3,
                            which2=which1,
                            which3=which2,
                            dispersion)

          f1 <- levels(x$model[,which3])

          f2 <- levels(x$model[,which1])[fl2] 

          f3 <- levels(x$model[,which2])[fl1] 

          m.inf <- subset(m.inf, 
                          group == f2 &
                          group2 == f3)[,3]

          rownames(m.inf) <- paste(f3,
                                   f2,
                                   f1,
                                   sep='/')

          ord <- order(as.vector(m.inf[,1]), 
                       decreasing=TRUE)

          m.inf <- cbind(m.inf[,1][ord],
                         m.inf[,2][ord],
                         m.inf[,3][ord])

          tab <- mt$tables[natri][[natri]]      
        } else if(nf1 == nt1 && nf2 == nt3) {
          r <- mt$n[names(mt$tables)][[natri]] # groups and its number of replicates

          bal <- ifelse(length(r) == 1,
                        TRUE,
                        FALSE) # is (or not) balanced          

          m <- as.vector(mt$tables[natri][[natri]][fl1,,fl2])

          which1 <- names(dimnames(mt$tables[natri][[natri]]))[1] 

          which2 <- names(dimnames(mt$tables[natri][[natri]]))[2]

          which3 <- names(dimnames(mt$tables[natri][[natri]]))[3]

          m.inf <- m.inf.3a(x,
                            which1=which2,
                            which2=which3,
                            which3=which1,
                            dispersion)

          f1 <- levels(x$model[,which2])

          f2 <- levels(x$model[,which3])[fl2] 

          f3 <- levels(x$model[,which1])[fl1] 

          m.inf <- subset(m.inf, 
                          group == f2 &
                          group2 == f3)[,3]

          rownames(m.inf) <- paste(f3,
                                   f2,
                                   f1,
                                   sep='/')

          ord <- order(m,
                       decreasing=TRUE)

          m.inf <- cbind(m.inf[,1][ord],
                         m.inf[,2][ord],
                         m.inf[,3][ord]) 

          tab <- mt$tables[natri][[natri]] 
        } else if(nf1 == nt3 && nf2 == nt2) {
          r <- mt$n[names(mt$tables)][[natri]] # groups and its number of replicates

          bal <- ifelse(length(r) == 1,
                        TRUE,
                        FALSE) # is (or not) balanced          

          m <- as.vector(mt$tables[natri][[natri]][,fl2,fl1])	

          which1 <- names(dimnames(mt$tables[natri][[natri]]))[1] 

          which2 <- names(dimnames(mt$tables[natri][[natri]]))[2]

          which3 <- names(dimnames(mt$tables[natri][[natri]]))[3]

          m.inf <- m.inf.3a(x,
                            which1=which1,
                            which2=which2,
                            which3=which3,
                            dispersion)

          f1 <- levels(x$model[,which1])

          f2 <- levels(x$model[,which2])[fl2] 

          f3 <- levels(x$model[,which3])[fl1] 

          m.inf <- subset(m.inf, 
                          group == f2 &
                          group2 == f3)[,3]

          rownames(m.inf) <- paste(f3,
                                   f2,
                                   f1,
                                   sep='/')

          ord <- order(m,
                       decreasing=TRUE)

          m.inf <- cbind(m.inf[,1][ord],
                         m.inf[,2][ord],
                         m.inf[,3][ord])

          tab <- mt$tables[natri][[natri]]   
        } else {
          r <- mt$n[names(mt$tables)][[natri]] # groups and its number of replicates

          bal <- ifelse(length(r) == 1,
                        TRUE,
                        FALSE) # is (or not) balanced          

          m <- as.vector(mt$tables[natri][[natri]][fl2,,fl1])

          which1 <- names(dimnames(mt$tables[natri][[natri]]))[1] 

          which2 <- names(dimnames(mt$tables[natri][[natri]]))[2]

          which3 <- names(dimnames(mt$tables[natri][[natri]]))[3]

          m.inf <- m.inf.3a(x,
                            which1=which2,
                            which2=which1,
                            which3=which3,
                            dispersion)

          f1 <- levels(x$model[,which2])

          f2 <- levels(x$model[,which1])[fl2] 

          f3 <- levels(x$model[,which3])[fl1] 

          m.inf <- subset(m.inf, 
                          group == f2 &
                          group2 == f3)[,3]

          rownames(m.inf) <- paste(f3,
                                   f2,
                                   f1,
                                   sep='/')

          ord <- order(m,
                       decreasing=TRUE)

          m.inf <- cbind(m.inf[,1][ord],
                         m.inf[,2][ord],
                         m.inf[,3][ord])

          tab <- mt$tables[natri][[natri]] 
        }
      }
  }
  colnames(m.inf) <- c('mean',
                       'min',
                       'max')

  mMSE <- MSE/r           

  dfr <- df.residual(x)  # residual degrees of freedom 

  out <- make.TukeyC.test(r=r,
                          MSE=MSE,
                          m.inf=m.inf,
                          ord=ord,
                          sig.level=sig.level,
                          dfr=dfr,
                          bal=bal,
                          mt=mt,
                          round)  

  class(out) <- c('TukeyC.nest',
                  'TukeyC',
                  'list')

  return(out)                             
}    
