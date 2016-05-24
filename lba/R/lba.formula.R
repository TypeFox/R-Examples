lba.formula <- function(formula,
                        data,
                        A           = NULL, # mixing parameters
                        B           = NULL, # latent components
                        K           =  1L,  # integer  
                        cA          = NULL, # position of the constraints mixing parameters
                        cB          = NULL, # position of the constraints latent components 
                        logitA      = NULL, # design matrix for row covariates IxS 
                        logitB      = NULL, # design matrix for row covariates JxT 
                        omsk        = NULL, # matrix of logit parameters SxK
                        psitk       = NULL, # matrix of logit parameters TxK
                        S           = NULL, # integer
                        T           = NULL, # integer
                        row.weights = NULL, # row weights matrix
                        col.weights = NULL, # column weights matrix
                        tolG        = 1e-10, 
                        tolA        = 1e-05, 
                        tolB        = 1e-05,
                        itmax.unide = 1e3, # With and without constraint
                        itmax.ide   = 1e3,
                        trace.lba   = TRUE,# only when K > 3
                        toltype     = 'all',
                        method      = c("ls", 
                                        "mle"),
                        what        = c('inner',
                                        'outer'),
                        ...)
{

 data <- as.data.frame(apply(data,
                             2,
                             factor))

 aux.form <- strsplit(as.character(formula),
                      '~')

 var.row <- aux.form[[2]]

 var.col <- aux.form[[3]]

 var.row1 <- unlist(strsplit(var.row,
                             ' \\+ '))

 var.col1 <- unlist(strsplit(var.col,
                             ' \\+ '))

 aux.tables <- rep(list(list()),
                   length(var.row1))

 aux.tables1 <- list()

 for(i in 1:length(var.row1)){

  for(j in 1:length(var.col1)){

   aux.tables[[i]][[j]] <- table(data[[var.row1[i]]],
                                 data[[var.col1[j]]])

  }

  aux.tables1[[i]]  <- do.call('cbind',
                               aux.tables[[i]]) 

 }

 tabs <- do.call('rbind',
                 aux.tables1)

 aux.namesR <- ifelse(length(length(var.row1))==1,
                      sapply(var.row1,
                             function(x)levels(data[[x]]),
                             simplify=F),
                      sapply(var.row1,
                             function(x)levels(data[[x]])))

 aux.namesR1 <- sapply(aux.namesR,
                       length)

 aux.namesR2 <- rep(var.row1,
                    aux.namesR1)

 rownames(tabs) <- paste(aux.namesR2,
                         dimnames(tabs)[[1]],
                         sep='')

 aux.namesC <- sapply(var.col1,
                      function(x)levels(data[[x]]),
                      simplify=F)

 aux.namesC1 <- sapply(aux.namesC,
                       length)

 aux.namesC2 <- rep(var.col1,
                    aux.namesC1)

 colnames(tabs) <- paste(aux.namesC2,
                         dimnames(tabs)[[2]],
                         sep='')

 switch(match.arg(what),
        inner = what <- 'inner',
        outer = what <- 'outer')

 switch(match.arg(method),
        ls = method <- 'ls',
        mle = method <- 'mle')

 if(is.null(cA) & is.null(cB) & is.null(logitA) & is.null(logitB)){

  class(tabs)  <- method

  result <- lba(tabs,
                A           =  A,           
                B           =  B,
                K           =  K,
                row.weights =  row.weights, 
                col.weights =  col.weights, 
                tolG        =  tolG,        
                tolA        =  tolA,        
                tolB        =  tolB,        
                itmax.unide =  itmax.unide,
                itmax.ide   =  itmax.ide,
                trace.lba   =  trace.lba,  
                toltype     =  toltype,
                what        =  what,
                ...) 

 } else 

  if((!is.null(cA) | !is.null(cB)) & is.null(logitA) & is.null(logitB)){

   class(tabs) <- paste(method,
                        'fe',
                        sep='.')
   result <- lba(tabs,
                 A           =  A,           
                 B           =  B,
                 K           =  K,
                 cA          =  cA,          
                 cB          =  cB,           
                 row.weights =  row.weights, 
                 col.weights =  col.weights, 
                 tolG        =  tolG,        
                 tolA        =  tolA,        
                 tolB        =  tolB,        
                 itmax.unide =  itmax.unide,
                 itmax.ide   =  itmax.ide,
                 trace.lba   =  trace.lba,  
                 toltype     =  toltype,
                 what        =  what,
                 ...)

  } else {

   class(tabs) <- paste(method,
                        'logit',
                        sep='.')

   result <- lba(tabs,
                 A           =  A,           
                 B           =  B,
                 K           =  K,
                 cA          =  cA,
                 cB          =  cB,
                 logitA      =  logitA,      
                 logitB      =  logitB,      
                 omsk        =  omsk,        
                 psitk       =  psitk,       
                 S           =  S,           
                 T           =  T,           
                 row.weights =  row.weights, 
                 col.weights =  col.weights, 
                 tolG        =  tolG,        
                 tolA        =  tolA,        
                 tolB        =  tolB,        
                 itmax.unide =  itmax.unide,
                 itmax.ide   =  itmax.ide,
                 trace.lba   =  trace.lba,  
                 toltype     =  toltype,
                 what        =  what,
                 ...)

  }

 cl <- match.call()

 result$call <- cl
 result$what <- what
 result$tab <- tabs

 class(result) <- c(class(result),
                    'lba.formula',
                    'lba')

 invisible(result)
}
