setMethod(f = "golden.search",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){
            
            P          = CGHo["nblevels"]
            M          = length(names(.Object@Y))
            n          = lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) 
            n          = sum(unlist(n))    
            select.tmp = CGHo["select"]  
            Kseq       = (M:multiKmax)+1
            instr      = fun2run(CGHo)

            Kmax         = multiKmax-M+1
            per          = 2  #percent
            Kseq1        = seq(0,Kmax)+M  #K=Kmin:Kmax
            Kseq1[1]     = M+1 
            select(CGHo) = select.tmp
            
            a         = 1
            b         = which(Kseq==Kseq1[length(Kseq1)])
            r         = (sqrt(5)-1)/2
            Kseqab    = c(a:b)
            s         = which.min(abs(a+(1-r)*(b-a)-Kseqab))
            x         = which.min(abs(a+r*(b-a)-Kseqab))
            
            multiKmax = Kseq[s]-1
            res       = eval(instr)
            Js        = -getmBIC(multiKmax,res$loglik,res$mu,CGHo)

            multiKmax = Kseq[x]-1
            res       = eval(instr)
            Jx        = -getmBIC(multiKmax,res$loglik,res$mu,CGHo) 

            xseq = c(x,s)
            Jseq = c(Jx,Js)

            
            while ((x-s)>0){    
              if (Js<Jx){ 
                b         = x
                x         = s
                Jx        = Js
                s         = which.min(abs(a+(1-r)*(b-a)-Kseqab))
                multiKmax = Kseq[s]-1
                res       = eval(instr)
                Js        = -getmBIC(multiKmax,res$loglik,res$mu,CGHo)

              } else {
                if (Js>Jx){
                  a         = s
                  s         = x
                  Js        = Jx
                  x         = which.min(abs(a+r*(b-a)-Kseqab))
                  
                  multiKmax = Kseq[x]-1
                  res       = eval(instr)
                  Jx        = -getmBIC(multiKmax,res$loglik,res$mu,CGHo)

                } else {
                  a         = s
                  b         = x
                  s         = which.min(abs(a+(1-r)*(b-a)-Kseqab))
                  x         = which.min(abs(a+r*(b-a)-Kseqab))                  
                  multiKmax = Kseq[s]-1
                  res       = eval(instr)
                  Js        = -getmBIC(multiKmax,res$loglik,res$mu,CGHo)

                  multiKmax = Kseq[x]-1
                  res       = eval(instr)
                  Jx        = -getmBIC(multiKmax,res$loglik,res$mu,CGHo)

                }    
              }
              xseq=c(xseq,c(x,s))
              Jseq=c(Jseq,c(Jx,Js))
            }
            gc()
            rg  = 1:length(Kseq)
            Kh  = Kseq[which((rg>=s) & (rg<=x))]
            return(Kh)
          }
          )
