enzCut <-
function(DNAsq, enznames, enzdata = enzdata) {
    
	if(length(DNAsq) >= 2){
	if(!inherits(DNAsq, "fasta"))
	{stop("The input DNAsq must be a \"fasta\" object.")}
	}
	
	REGEXPR <- 
    function(Echar) {
     loc = gsub("V", "[^T]", gsub("H", "[^G]{1}", gsub("D", "[^C]{1}", 
         gsub("B", "[^A]{1}", gsub("W", "[AT]{1}", gsub("S", "[GC]{1}", 
            gsub("K", "[GT]{1}", gsub("M", "[AC]{1}", gsub("Y", 
                "[CT]{1}", gsub("R", "[GA]{1}", gsub("N", ".", 
                  gsub("[',_]", "", Echar))))))))))))
       return(loc)
    }
	
    DNAsq = as.character(toupper(DNAsq))
	enz = selEnz(names = enznames, enzdata = enzdata)
    for (i in 1:length(enz$site)) {
        enzsq = toupper(enz$site[i])
        for (e in 1:2) {
            if (e == 1) {
                Enzchar = gsub("_", "", enzsq)
            }
            if (e == 2) {
                if (regexpr("_", enzsq) < 0) 
                  enzsq = gsub("'", "_", enzsq)
                Enzchar = gsub("-", "_", as.character(revComp(gsub("_", "-", gsub("'", "", enzsq)))))
            }
            up = min(unlist(gregexpr("[^N_']", Enzchar)))
            lo = max(unlist(gregexpr("[^N_']", Enzchar)))
            loc = REGEXPR(substr(Enzchar, up, lo))
            dis = regexpr("['_]", Enzchar) - up
            adFrag = ifelse(dis < 0, dis + 1, dis)
            if (regexpr(loc, DNAsq) > 0) 
                newsq = unlist(gregexpr(loc, DNAsq)) + adFrag
            if (regexpr(loc, DNAsq) < 0) 
                newsq = 0
            if (e == 1 & i == 1) 
                RFLP = newsq
            if (i != 1 | e != 1) 
                RFLP = sort(union(RFLP, newsq))
        }
    }
    if (max(RFLP) > 0) {
        RFLP.site = RFLP[RFLP > 0 & RFLP < nchar(DNAsq)]
        si = c(RFLP.site, (nchar(DNAsq) + 1))
        RFLP.frag = c()
        for (k in 1:length(si)) {
            if (k == 1) 
                RFLP.frag[k] = si[k] - 1
            if (k > 1) 
                RFLP.frag[k] = si[k] - si[k - 1]
        }
        T5 = range(RFLP.site)[1] - 1
        T3 = nchar(DNAsq) - max(RFLP.site) + 1
    }
    if (max(RFLP) == 0) {
        RFLP.site = 0
        RFLP.frag = nchar(DNAsq)
        T5 = nchar(DNAsq)
        T3 = nchar(DNAsq)
    }
    return(list(enz = enz, RFLP.site = RFLP.site, RFLP.frag = RFLP.frag, TRFLP = c(T5 = T5, T3 = T3)))
}

