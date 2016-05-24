uniformMap <-
function(Mb=NULL, cM=NULL, M=NULL, cm.per.mb=1, chromosome=1) { #genL numeric of length 1 or 2: genetic length male & female
    if(!is.null(cM) && !is.null(M)) stop("Either 'cM' or 'M' must be NULL.")
    stopifnot(!is.null(cM) || !is.null(M) || !is.null(Mb))
    
    if(is.null(cM)) 
        if(!is.null(M)) cM = M*100 else  cM = cm.per.mb * Mb
    if(is.null(Mb)) Mb = cM / cm.per.mb
    
    if(is.character(chromosome) && tolower(chromosome)=="x") 
        chromosome = 23

    map = switch(max(length(Mb),length(cM)), {
        m = cbind(Mb=c(0, Mb), cM=c(0,cM))
        list(male=m, female=m)
    },{
        if (length(cM)==1) cM=c(cM,cM)
        if (length(Mb)==1) Mb=c(Mb,Mb)
        list(male=cbind(Mb=c(0, Mb[1]), cM=c(0,cM[1])), female=cbind(Mb=c(0, Mb[2]), cM=c(0,cM[2])))
    })
    female_phys = as.numeric(map$female[2,1])
    if(chromosome==23)
        map$male = NA
    else
        if(female_phys != map$male[2,1]) stop("Male and female chromosomes must have equal physical length.")  
    structure(map, length_Mb = female_phys, chromosome = chromosome, class = "chromosomeMap")
}
 
loadMap <-
function(map, chrom=NULL) {
  
    if (is.character(map)) {
        CHROM.LENGTH = cbind(male_morgan=c(1.9, 1.752, 1.512, 1.35, 1.302, 1.162, 1.238, 1.089, 1.047, 1.147,0.992, 1.154, 0.919, 0.857, 0.825, 0.88, 0.863, 0.737, 0.708, 0.563, 0.426, 0.45, NA), 
        female_morgan=c(3.34, 3.129, 2.687, 2.6, 2.49, 2.333, 2.21, 2.042, 1.879, 2.062, 1.886, 1.974, 1.468, 1.24, 1.393, 1.494, 1.529, 1.379, 1.152, 1.109, 0.67, 0.662, 1.733),
        Mb = c(247.2, 242.7, 199.3, 191.1, 180.6, 170.8, 158.7, 146.2, 140.1, 135.3, 134.4, 132.3, 114.1, 105.3, 100.2, 88.7, 78.6, 76.1, 63.8, 62.4, 46.9, 49.5, 154.6))

        if(is.null(chrom)) chrom=1:22
        if(is.character(chrom)) if(chrom=="AUTOSOMAL") chrom=1:22 else if(chrom=="X") chrom=23
        if (any(!chrom %in% 1:23)) stop(paste("Invalid chromosome number(s):", paste(setdiff(chrom, 1:23), collapse=", ")))
        
        maps = switch(tolower(map), 
        decode = IBDsim::DecodeMap[chrom], 
        uniform.sex.spec = lapply(chrom, function(chr) {dat=as.numeric(CHROM.LENGTH[chr,]); uniformMap(M=dat[1:2], Mb=dat[3], chromosome=chr)}),
        uniform.sex.aver = lapply(chrom, function(chr) {dat=as.numeric(CHROM.LENGTH[chr,]); uniformMap(M=mean(dat[1:2]), Mb=dat[3], chromosome=chr)}),
        stop("Invalid map name"))
    }
    else {
        maps = map
        if (inherits(maps, "chromosomeMap")) maps = list(maps)
    }
    attr(maps, 'length_Mb') = sum(unlist(lapply(maps, attr, "length_Mb")))
    maps
}

cm2phys <-
function(cM_locus, mapmat) {    # mapmat matrise med kolonner 'Mb' og 'cM'
    last = mapmat[nrow(mapmat), ]
    nontriv = cM_locus >= 0 & cM_locus <= last[['cM']]
    res = numeric(length(cM_locus))
    res[!nontriv] <- NA
    cm <- cM_locus[nontriv]
    interv = findInterval(cm, mapmat[, 'cM'], all.inside=TRUE) #.findInterval_quick not allowed
    res[nontriv] = mapmat[interv, 'Mb'] + (mapmat[interv+1, 'Mb'] - mapmat[interv, 'Mb']) * (cm - mapmat[interv, 'cM'])/(mapmat[interv+1, 'cM'] - mapmat[interv, 'cM'])
    res
}

#25.3.2014 (not used in IBDsim)
phys2cm <-
function(Mb_locus, mapmat) {    # mapmat matrise med kolonner 'Mb' og 'cM'
    last = mapmat[nrow(mapmat), ]
    nontriv = Mb_locus >= 0 & Mb_locus <= last[['Mb']]
    res = numeric(length(Mb_locus))
    res[!nontriv] <- NA
    mb <- Mb_locus[nontriv]
    interv = findInterval(mb, mapmat[, 'Mb'], all.inside=TRUE) #.findInterval_quick not allowed
    res[nontriv] = mapmat[interv, 'cM'] + (mapmat[interv+1, 'cM'] - mapmat[interv, 'cM']) * (mb - mapmat[interv, 'Mb'])/(mapmat[interv+1, 'Mb'] - mapmat[interv, 'Mb'])
    res
}

# Not allowed:
# .findInterval_quick = function(x, vec, all.inside=FALSE) {
    # nx <- length(x)
    # .C("find_interv_vec", xt = as.double(vec), n = length(vec), 
        # x = as.double(x), nx = nx, FALSE, 
        # as.logical(all.inside), index = integer(nx), DUP = FALSE, NAOK = TRUE, 
        # PACKAGE = "base")$index
# }

.genmapC <-
function(cM, mapmat) {    # mapmat matrise med kolonner 'pos' og 'cM'   ###NB gammel!!
    map_phys = as.double(mapmat[,1])
    map_cm = as.double(mapmat[,2])
    n = length(cM)
    .C("geneticMap", map_phys, map_cm, cM, n, res=numeric(n), PACKAGE="IBDsim")$res
}