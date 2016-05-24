# Verified 1.3.18
disaggregate.ts <-
function(x, ...) {
        if ( class(x) != "ts" ) { stop("stream must be of class <ts>.") }
        if ( frequency(x) != 12 ) { stop("stream must be of frequency 12.") }
        
        ast = -9999
        fun = median
        l = list(...)
        if ( length(l) > 0 ) {
                if ( ( "asterisk" %in% names(l) ) == T ) {
                        ast = l$asterisk
                } 
                if ( ( "fun" %in% names(...) ) == T ) {
                        fun = l$fun
                }
        }
        data.ts = x
        p.eps = 1e-2
        if ( !any(data.ts==ast, na.rm = T) ) { 
                warning("Nothing to do, data complete in given range.")
                return(x)
        }
        
        data.ts[is.na(data.ts)] = ast
        table.match = match(data.ts, ast)
        table.match[is.na(table.match)] = 0L
        END = length(data.ts)
        ENDTM = length(table.match)
        
        nidx = (1:END)[table.match == 0]
        if (nidx[length(nidx)] < END ) {
                data.ts = ts(data.ts[-((nidx[length(nidx)]+1):END)], start=start(data.ts),frequency=frequency(data.ts))
                table.match = match(data.ts, ast)
                table.match[is.na(table.match)] = 0L
                END = length(data.ts)
                ENDTM = length(table.match)
        }
        
        m = start(data.ts)[2]
        m = 0:11 + m # January? not always
        m [ m > 12 ] = m[ m > 12 ] - 12
        p = rep(NA, 12)
        
        idx = (1:ENDTM)[table.match == 1]
        d.ts = data.ts
        d.ts[ idx ] = NA
        if ( length(idx) > 0 ) {
                if ( idx[length(idx)] == ENDTM ) { idx = idx[-length(idx)] }
                d.ts[ idx + 1 ] = NA
        }
        
        if ( length(data.ts) < 12 ) {
                p = rep(1, 12)
        } else {
                for ( i in  m ) {
                        y = d.ts[seq(i, length(d.ts), 12)]
                        p[i] = apply(matrix(y, nrow=1), 1, fun, na.rm=T) # fun=median
                }
        }
        # Sanity check
        if ( any(is.na(p)) ) {
                warning("Some weights are NA due to too short ts. Equal weights assigned.\n  *** Discarding of this series is recomenden.")
                p = rep(1, 12)
        }
        if ( any(p < p.eps) ) {
                warning(paste0("Some weights are too small (<",p.eps,") !"))
                p[p<p.eps] = median(p)
        }
        if ( any(p < p.eps) ) {
                warning(paste0("Some weights are still too small (<",p.eps,") !\nEqual weights assigned."))
                p = rep(1, 12)
        }
        
        idx_s = 1
        while ( idx_s < END ) {
                idx_s = idx_s - 1 + match(1L, table.match[idx_s:END])
                if ( is.na(idx_s) ) { break }
                idx_e = (idx_s) + match(0L, table.match[(idx_s+1):END])
                sanity = data.ts[idx_e]
                mfix = m[m12(idx_s)] : m[m12(idx_e)]
                if ( m[m12(idx_s)] > m[m12(idx_e)] ) {
                        mfix = c(m[m12(idx_s)] : 12, 1 : m[m12(idx_e)]) 
                }
                mfix2 = m12(mfix[length(mfix)]+1) : m12(mfix[1]-1) 
                if ( idx_e - idx_s + 1 > 12 ) {
                        if ( m12(mfix[length(mfix)]+1) > m12(mfix[1]-1) ) {
                                mfix2 = c(m12(mfix[length(mfix)]+1) : 12, 1 : m12(mfix[1]-1)) 
                        }
                        mfix = c(mfix, mfix2)
                }        
                pp = rep(p[mfix], length.out=length(idx_s:idx_e))
                data.ts[idx_s:idx_e] = rep(pp/sum(pp)  * data.ts[idx_e], length.out=length(idx_s:idx_e))
                if ( abs(sum(data.ts[idx_s:idx_e]) - sanity ) > 1e-9) { stop("Sanity checksum error!") }
                idx_s = idx_e + 1
        }
        return(data.ts)
}
