paraMix <-
function (x, R, id.U, id.V = NULL, alleles, afreq = NULL, Xchrom = FALSE, 
    known_genotypes = list(), loop_breakers = NULL, eliminate = 0, check = TRUE, plot = TRUE, title=NULL) 
{
    if (inherits(x, "linkdat")) 
        x = list(x)
    all_typed = sapply(known_genotypes, "[", 1)
    contrib_typed = setdiff(all_typed, id.V)
    contrib_untyped = id.U
    K = unique(unlist(lapply(known_genotypes, function(triple) if (triple[1] %in% contrib_typed) triple[2:3])))
    R_not_masked = setdiff(R, K)
    if (length(alleles) == 1) 
        alleles = 1:alleles
    if (check) 
        checkInput(x, R, id.U, id.V, alleles, all_typed, K, R_not_masked)
    partialmarkers = lapply(x, function(ped) {
        m = marker(ped, alleles = alleles, afreq = afreq, chrom = if (Xchrom) 23 else NA)
        for (tup in known_genotypes) if (tup[1] %in% ped$orig.ids) 
            m = modifyMarker(ped, m, ids = tup[1], genotype = tup[2:3])
        m
    })
    if(!is.logical(plot) && plot!="plot_only") stop("'plot' must be either TRUE, FALSE or the character 'plot_only'.")
    if (isTRUE(plot) || plot=="plot_only") {
        topmargin = ifelse(isTRUE(title==""), 0, 2) 
        op = par(oma = c(0, 0, topmargin, 0), xpd = NA)
        layout(rbind(1:length(x)), widths = ifelse(sapply(x, 
            inherits, what = "singleton"), 1, 2))
        for (i in 1:length(x)) {
            ped = x[[i]]
            mm = if (length(known_genotypes) > 0) partialmarkers[[i]] else NULL
            cols = rep(1, ped$nInd)
            cols[ped$orig.ids %in% id.U] = 2
            cols[ped$orig.ids %in% id.V] = 4
            plot(ped, marker = mm, col = cols, mar = c(1, 2, 1, 2), title = "")
        }
        if (is.null(title)) {
         mixtext = paste("Mixture =", paste(R, collapse = "/"))
         idtext = paste("id.U =", if (length(id.U) > 0) paste(id.U, collapse = ",") else "None")
         idVtext = paste("id.V =", if (length(id.V) > 0) paste(id.V, collapse = ",   ") else "None")
         title = paste(mixtext, idtext, idVtext, sep = ";  ")
        }
        mtext(title, outer = T)
        par(op)
        if (plot=="plot_only") return()
    }
    
    R_INTERNAL = match(R, alleles)
    allgenos = allGenotypes(length(alleles))
    allgenos_in_R = which(allgenos[, 1] %in% R_INTERNAL & allgenos[, 2] %in% R_INTERNAL)
    id_split = lapply(x, function(ped) intersect(id.U, ped$orig.ids))
    id_order = do.call(c, id_split)
    a = lapply(1:length(x), function(i) {
        allgenos_row_index = geno.grid.subset(x[[i]], 
            partialmarkers[[i]], id_split[[i]], make.grid = F)
        lapply(allgenos_row_index, intersect, y = allgenos_in_R)
    })
    mixgrid = fast.grid(unlist(a, recursive = FALSE))
    if (length(R_not_masked) > 0) {
        R_not_masked_INTERNAL = match(R_not_masked, alleles)
        mixgrid = mixgrid[apply(mixgrid, 1, function(r) all(R_not_masked_INTERNAL %in% allgenos[r, ])), , drop = F]
    }
    if (!is.null(loop_breakers)) {
        for (i in 1:length(x)) {
            ped = x[[i]]
            lpb = loop_breakers[loop_breakers %in% ped$orig.ids]
            if (length(lpb) > 0) {
                y = breakLoops(setMarkers(ped, partialmarkers[[i]]), lpb)
                partialmarkers[[i]] = y$markerdata[[1]]
                y$markerdata = NULL
                x[[i]] = y
            }
        }
    }
    if (length(id.U) == 0) {
        lik = likelihood(x, locus1 = partialmarkers, eliminate = eliminate)
        return(list(likelihood = lik, all.likelihoods = lik))
    }
    if (nrow(mixgrid) == 0) {
        cat("No compatible genotype combinations\n")
        return(list(likelihood = 0, all.likelihoods = numeric(0)))
    }
    dat = lapply(id_order, function(i) {
        pedno = which(sapply(x, function(ped) i %in% ped$orig.ids))
        c(pedno = pedno, intid = match(i, x[[pedno]]$orig.ids))
    })
    all.likelihoods = apply(mixgrid, 1, function(allg_rows) {
        for (i in 1:length(id_order)) {
            genotyp = allgenos[allg_rows[i], ]
            pedno = dat[[i]]["pedno"]
            intid = dat[[i]]["intid"]
            partialmarkers[[pedno]][intid, ] = genotyp
        }
        likelihood(x, locus1 = partialmarkers, eliminate = eliminate)
    })
    list(likelihood = sum(all.likelihoods), all.likelihoods = all.likelihoods[all.likelihoods > 
        0])
}

checkInput= function (x, R, id.U, id.V, alleles, all_typed, K, R_not_masked) {
   if (any(duplicated(unlist(lapply(x, function(ped) ped$orig.ids)))))  # new
      stop("The input pedigrees must have disjoint sets of ID labels")
   if (length(intersect(id.U, id.V)) > 0) 
      stop("Contributors and noncontributors overlap")
   if (any(id.U %in% all_typed)) 
      stop("Input error. Unknown persons typed")
   if (!all(R %in% alleles))    # changed
      stop("Mixture contains non-existing alleles")
   if (!all(K %in% R)) 
      stop("Input error. Known contributors with alleles outside mixture")
   if (length(R_not_masked) > 2 * length(id.U)) 
      stop("Input error. Not enough unknown contributors")
   if (length(R_not_masked) > 8) 
      stop("No of unmasked alleles greater than 8 probably prohibitive. Set check=FALSE to try")
   if (length(length(id.U)) > 8) 
      stop("No of unknown greater than 8 probably prohibitive. Set check=FALSE to try")
}


# .checkInput=function(x,R,id.U,id.V,partialmarker,all_typed,K,R_not_masked){
#    if (length(intersect(id.U, id.V)) > 0)
#       stop("Contributors and noncontributors overlap") # endret teksten her
#    if (any(id.U %in% all_typed))
#       stop("Input error. Unknown persons typed")
#    if (!all(R %in% as.integer(attr(partialmarker, "alleles"))))
#       stop("Mixture contains non-existing alleles")
#    if (!all(K %in% R))
#       stop("Input error. Known contributors with alleles outside mixture")
#    if (length(R_not_masked) > 2*length(id.U))
#       stop("Input error. Not enough unknown contributors")
#    if (length(R_not_masked) > 8)
#       stop("No of unmasked alleles greater than 8 probably prohibitive. Set check=FALSE to try")	
#    if (length(length(id.U)) > 8)
#       stop("No of unknown greater than 8 probably prohibitive. Set check=FALSE to try")	  
#    	  
#    0
#    }

# #seed=177
# .Example2=function(id.victim=9,id.suspect=10,conditional=FALSE,nMarkers=23,seed=NULL,check=TRUE){
# if(!is.null(seed)) set.seed(seed)
# data(db)
# x=cousinPed(1)
# x=swapSex(addOffspring(x,father=7,mother=8,noffs=2),ids=9)
# x=addOffspring(x,father=1,noffs=1) #Pedigree has been defined
# index=c(setdiff(seq(2,10,2),id.victim),11)
# #loglik=rep(0,12)
# loglik=rep(0,length(index))
# names(loglik)=index
# nn=as.character(unique(db[,1]))
# ret=x
# N=min(length(nn),nMarkers)
# RR=vector("list",N)
# 
# for (k in 1:N){
#   
#  cat(k,"\n")
#  afreq1=as.double(db[db$Marker==nn[k],3])
#  afreq1=round(afreq1,7)
#  afreq1=afreq1/sum(afreq1)
#  alleles=as.integer(1:length(afreq1))
#  m1=marker(x,alleles=alleles,afreq=afreq1)
#  if(conditional){
#    y=markerSim(x,N=1,available=id.victim,partialmarker=m1,verbose=FALSE,loop_breakers=7,seed=NULL)
#    ret=addMarker(ret,y$markerdata[[1]])
#    m1=y$markerdata[[1]]
#    y=markerSim(x,N=1,available=id.suspect,partialmarker=m1,verbose=FALSE,loop_breakers=7,seed=NULL)
#  }
#  if(!conditional){
#    y=markerSim(x,N=1,available=c(id.victim,id.suspect),partialmarker=m1,verbose=FALSE,loop_breakers=7,seed=NULL)
#    ret=addMarker(ret,y$markerdata[[1]])
#    }
#  R=simMixParamlink(y,alleles)[[1]] #Finds mixtures
#  RR[[k]]=R
#  known=list(c(id.victim,as.integer(y$markerdata[[1]][id.victim,])))
#  for(i in 1:length(index)){
#    foo=paraMix(x,R,id.U=index[i],id.V=NULL,alleles=alleles,afreq=afreq1,loop_breakers=7,known_genotypes=known,check=check,plot=FALSE)
#    loglik[i]=loglik[i]+log(foo$likelihood)
#  }
# }
# #plot(x,marker=1:2,title="")
# #loglik=loglik[index]
# post=exp(loglik)/sum(exp(loglik))
# names(post)=index
# list(post=post,ret=ret,mixture=RR)
# }


