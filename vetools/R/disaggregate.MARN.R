# Verified 1.3.18
disaggregate.MARN <-
function(stream=NULL, reference=NULL, na.action="error", asterisk=-9999, date.eps=0.004, float.eps=0.0001, return.incomplete=TRUE) {
        if ( class(stream) != "ts" ) { stop("stream must be of class <ts>.") }
        if ( class(reference) != "ts" ) { stop("reference stream must be of class <ts>.") }
        data.ts = stream
        
        if ( !any(data.ts==asterisk, na.rm = T) ) { 
                warning("Nothing to do, data complete in given range.")
                return(stream)
        }
        
        ref.ts = reference
        ast = asterisk
        f.eps = float.eps
        
        time.i = time(data.ts)
        table.match=match(data.ts, ast)
        table.match[is.na(table.match)] = 0L
        END = length(table.match)
        
        idx_s = 1
        while ( idx_s < END ) {
                idx_s = idx_s - 1 + match(1L, table.match[idx_s:END])
                if ( is.na(idx_s) ) { break }
                idx_e = (idx_s) + match(0L, table.match[(idx_s+1):END]) # Includes the TOTAL ACC. VALUE
                
                if ( time.i[idx_s] < start(ref.ts) ) {
                        pond = rep(NA, idx_e - idx_s + 1)
                } else {
                        pond = window(ref.ts, start=time.i[idx_s], end=time.i[idx_e]+date.eps, extend=T)
                }
                n = sum(pond[!is.na(pond)])
                if ( (all(is.na(pond))) | (abs(n) < f.eps) | any(is.nan(n)) ) {
                        if ( na.action %in% c("mean", "average", "warning", "continue")) {
                                warning(paste0("Reference series DOES NOT HAVE data at a needed time window: ", time.i[idx_s], ", [", idx_s, "] -> ", time.i[idx_e], ", [", idx_e, "]. Average will be used."))
                        } else {
                                if ( return.incomplete ) {
                                        cat(paste0("Error: Reference series DOES NOT HAVE data at a needed time window: ", time.i[idx_s], ", [", idx_s, "] -> ", time.i[idx_e], ", [", idx_e, "]. Returning partially fixed series."))
                                        return(data.ts)
                                } else {
                                        stop(paste0("Reference series DOES NOT HAVE data at a needed time window: ", time.i[idx_s], ", [", idx_s, "] -> ", time.i[idx_e], ", [", idx_e, "]"))
                                }
                                
                        }
                        n = length(pond)
                        pond = rep(1/n, n)
                } else {
                        pond = pond / n
                }
                if (length(window(data.ts, start=time.i[idx_s], end=time.i[idx_e])) != length(pond)){
                        stop(paste0("LENGTH DIFFER: stream [",length(window(data.ts, start=time.i[idx_s], end=time.i[idx_e])),"] vs. reference [",length(pond),"]. Modify date.eps!"))
                }
                data.ts[idx_s:idx_e] <- data.ts[idx_e] * pond
                idx_s = idx_e + 1
        }
        # Mass loss verification
        tmp = stream
        tmp[is.na(tmp)] = 0
        tmp[tmp == ast] = 0
        if (abs(sum(tmp) - sum(data.ts[!is.na(data.ts)])) > 1e-8){
                stop("Something wrong happened: Desagregated series sum mismatch.")
        }
        return(data.ts)
}
