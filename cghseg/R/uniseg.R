setMethod(f = "uniseg",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax=NULL){

            uniKmax    = getuniKmax(.Object,CGHo,uniKmax)
            CGHr       = new("CGHresults",CGHd=.Object,CGHo=CGHo)
            from(CGHr) = "uniseg"

            if (CGHo["wavenorm"]!="none"){
              cat("[uniseg] CGHo[\"wavenorm\"] is not considered in uniseg \n")
            }
            if (CGHo["GCnorm"]!="none"){
              cat("[uniseg] CGHo[\"GCnorm\"] is not considered in uniseg \n")
            }
            
            if (CGHo["calling"]){
              instr = parse(text = "out = unisegclust(Ym,CGHo,Kmax)")
            } else {
              instr = parse(text = "out = unisegmean(Ym,CGHo,Kmax)")
            }
            instr2 = parse(text = "invisible(list(mu=out$mu,loglik=out$loglik))" )
            
            Res = lapply(names(.Object@Y), FUN = function(m){      
              Ym   = .Object@Y[[m]]
              Kmax = uniKmax[[m]]
              eval(instr)
              eval(instr2)
            })
            names(Res)    = names(.Object@Y)
            out.mu        = lapply(Res,FUN = function(x){x$mu})
            
            if (!is.null(.Object["genomic.position"])){
              x      = .Object["genomic.position"]
              out.mu = lapply(out.mu,FUN = function(y){
                       y$begin = x[y$begin]
                       y$end   = x[y$end]
                       return(y)
                     })
            }
            
            mu(CGHr)      = out.mu
            out.ll        = lapply(Res,FUN = function(x){x$loglik})
            loglik(CGHr)  = out.ll
            return(CGHr)
          } # end of uniseg
          )











