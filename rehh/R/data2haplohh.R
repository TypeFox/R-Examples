setClass(Class = "haplohh",
	representation(haplo = "matrix",position = "numeric",snp.name="character",chr.name = "numeric",nhap = "numeric",nsnp = "numeric")
)


is.haplohh <- function (x) 
{
    res <- (is(x,"haplohh") & validObject(x))
    return(res)
}


make.example.files <- function() {
     file.copy(system.file('bta12_hapguess_switch.out',package='rehh'),'bta12_hapguess_switch.out')
     file.copy(system.file('map.inp',package='rehh'),'map.inp')
     file.copy(system.file('bta12_cgu.hap',package='rehh'),'bta12_cgu.hap')
}

data2haplohh<-function(hap_file,map_file,min_maf=0,min_perc_geno.hap=100,min_perc_geno.snp=100,chr.name=NA,popsel=NA,recode.allele=FALSE){
res<-new("haplohh")

#####Fichier map (verification)
map<-read.table(map_file,row.names=1) #snp name, chromosome, position, allele ancestral, allele derive
 if(ncol(map)!=4){
   cat("Wrong format for map data file: ",map_file,"\n Should contain 5 columns (Snp name, Chromosome Number, Position, Ancestral Allele and Derived Allele) and no header\n")
  stop("Conversion stopped")
 }else{
   tmp_chr=unique(as.numeric(map[,1]))
   if(length(tmp_chr)!=1 & is.na(chr.name)){
    cat("More than one chromosome name in Map file:",map_file,"\n")
       repeat {
         cat('Which chromosome should be considered among:',tmp_chr,"\n")
         chr.name <- scan(file = '',n = 1,what = integer(0),quiet = TRUE)
         if (chr.name %in% tmp_chr) break
         }
    }
    
    if(is.na(chr.name)){chr.name=tmp_chr[1]}
     map=map[as.character(map[,1])==chr.name,]
     tmp_nom<-rownames(map)
     if(nrow(map)==0){
      cat("No SNPs mapping to chromosome ",chr.name," are found in the map file\n")
      tmp_chr=unique(as.character(map[,1]))
      cat("Here is the list of chromosome names in the map file:",tmp_chr,"\n")
      stop("Conversion stopped")
      }
   #premier checks OK
    res@chr.name<-chr.name
    res@snp.name<-tmp_nom
    tmp_pos<-as.numeric(map[,2])

    if(sum(diff(tmp_pos)<0)>0){
      stop("SNP should be ordered on the map, check also that both haplotypes (column of haplo) and map (row of map) are ordered in the same way")}
    if(sum(diff(tmp_pos)==0)>0){
      warning("Some SNPs map to the same position")}
    res@position<-tmp_pos
  }

res@nsnp=nrow(map)
cat("Map file seems OK:",res@nsnp," SNPs declared for chromosome",res@chr.name,"\n")

###Fichier haplo

out_fphase<-scan(hap_file,what="character",sep="\n",quiet=TRUE,nlines=15)
test_fphase_1=grep("fastPHASE",out_fphase)
test_fphase_2=grep("BEGIN COMMAND_LINE",out_fphase)

 if(length(test_fphase_1)>0 & length(test_fphase_2)>0){ #fichier fastphase
  cat("Looks like a FastPHASE haplotype file\n")
  out_fphase<-scan(hap_file,what="character",sep="\n",quiet=TRUE)
  deb_geno=grep("BEGIN GENOTYPES",out_fphase)[1] + 1
  fin_geno=grep("END GENOTYPES",out_fphase)[1] - 1
  out_fphase<-out_fphase[deb_geno:fin_geno]

 test_poplabel=grep("subpop. label:",out_fphase)
 if(length(test_poplabel)>1){
   nom_hap_cplet=out_fphase[test_poplabel]
   nhap_tot=length(nom_hap_cplet)
   pop_label=numeric(nhap_tot)
   tmp_poplab=(strsplit(nom_hap_cplet,split="subpop. label:"))
   for(i in 1:nhap_tot){
     pop_label[i]=as.numeric(unlist(strsplit(tmp_poplab[[i]][2],split="\\(internally"))[1])
    }
    liste_pop=unique(pop_label)
    cat("Haplotypes originate from ",length(liste_pop)," different populations in the fastPhase output file\n")
   
    if(!(popsel %in% pop_label)){
     cat("Chosen pop. is not in the list of pop. number:",liste_pop,"\n")
     popsel=NA
    }
    if(is.na(popsel)){ #on demande interactivement la pop a selectionner
       repeat {
         cat('Which population should be considered among:',liste_pop,"\n")
         popsel <- scan(file = '',n = 1,what = integer(0),quiet = TRUE)
         if (popsel %in% pop_label) break
         }
    }
   hapsel=(which(pop_label==popsel) -1)*3 + 1
   hapsel=sort(as.numeric(cbind(hapsel,hapsel+1,hapsel+2)))
   out_fphase=out_fphase[hapsel]
   }else{
    cat("Haplotypes seems to originate from only one population\n")
    if(!(is.na(popsel))){ #bizarre car pas de sous population observe
         cat('No popsel are thus considered:',liste_pop,"\n")
    }
}

 nlignes=length(out_fphase) ; nind=nlignes/3 ; nhap=2*nind
 tmp_haplos=matrix(,nhap,res@nsnp)
   for(i in 1:nind){
    id_line=3*(i-1) + 1 ; hap1_line=id_line+1 ; hap2_line=id_line+2
    hap1_index=2*(i-1) +1 ; hap2_index=hap1_index + 1
  
    hap1=unlist(strsplit(out_fphase[hap1_line],split=" "))
    if(length(hap1)!=res@nsnp){
      stop("Number of snp in haplotypes differ from the number declared in the map file")
      }else{tmp_haplos[hap1_index,]=hap1}
  
    hap2=unlist(strsplit(out_fphase[hap2_line],split=" "))
    if(length(hap2)!=res@nsnp){
      stop("Number of snp in haplotypes differ from the number declared in the map file")
      }else{tmp_haplos[hap2_index,]=hap2}
   }

}else{ #fichier au format standard
   cat("Standard rehh input file assumed\n")
   tmp_haplos=as.matrix(read.table(hap_file,row.names=1))
   if(ncol(tmp_haplos)!=res@nsnp){
    cat("The number of snp in the haplotypes",ncol(tmp_haplos)," is not equal to the number of snps declared in the map file",res@nsnp,"\n")
    stop("Conversion stopped")
  }
}

 res@nhap=nrow(tmp_haplos)

##recodage des alleles si necessaire
   if(recode.allele){
    cat("Please Wait: alleles are recoded according to map file as: 0 (missing data), 1 (ancestral allele) or 2 (derived allele)\n")
    res@haplo=matrix(0,res@nhap,res@nsnp)
    anc_allele=as.character(map[,3]) ; der_allele=as.character(map[,4])
   for(i in 1:res@nsnp){
      all_anc= (as.character(tmp_haplos[,i])==anc_allele[i])
      all_der= (as.character(tmp_haplos[,i])==der_allele[i])
      all_miss= ((all_anc + all_der)==0)
      res@haplo[all_anc,i]=1 ; res@haplo[all_der,i]=2 ; res@haplo[all_miss,i]=0
      tmp_nmiss=sum(all_miss)
      if(tmp_nmiss>0){cat(tmp_nmiss," missing allele for SNP: ",res@snp.name[i],"\n")}
      if(i%%10==0){cat("SNP #",i," out of",res@nsnp,"\n")}
    }
  }else{
   if(sum(!(tmp_haplos==0 | tmp_haplos==1 | tmp_haplos==2 ))>0){
        cat("Alleles are not coded in the appropriate format: 0 (missing data), 1 (ancestral allele) or 2 (derived allele)\n")
        cat("Check your data or use recode.allele=TRUE option to recode according to the map information\n")
        stop("Conversion stopped")
      }
   res@haplo=matrix(as.numeric(tmp_haplos),res@nhap,res@nsnp)
  }

#selection des haplos d'apres les donnees manquantes
 if(min_perc_geno.hap<100){
   cat("Discard Haplotype with less than ",min_perc_geno.hap,"% of genotyped SNPs\n")
   hap_sel=(100*rowSums(res@haplo!=0)/res@nsnp)>=min_perc_geno.hap
   if(sum(hap_sel)==res@nhap){
    cat("No haplotype discarded\n")
   }else{
    cat(res@nhap-sum(hap_sel)," Haplos discarded\n")
    res@haplo=res@haplo[hap_sel,]
    res@nhap=sum(hap_sel)
    cat(res@nhap," Haplos remaining\n")
   }
 }

#selection des snps d'apres les donnees manquantes
 if(min_perc_geno.snp<100){
   cat("Discard SNPs genotyped on less than ",min_perc_geno.snp,"% of haplotypes\n")
   snp_sel=(100*colSums(res@haplo!=0)/res@nhap)>=min_perc_geno.snp
   if(sum(snp_sel)==res@nsnp){
    cat("No SNPs discarded\n")
   }else{
    cat(res@nsnp-sum(snp_sel)," SNPs discarded\n")
    res@haplo=res@haplo[,snp_sel]
    res@nsnp=sum(snp_sel)
    res@position=res@position[snp_sel]
    res@snp.name=res@snp.name[snp_sel]
    cat(res@nsnp," SNPs remaining\n")
   }
 }

#selection des snps sur MAF
 if(min_maf>0){
   cat("Discard SNPs with MAF below ",min_maf,"\n")
   tmp_n1=colSums(res@haplo==1) ; tmp_maf=tmp_n1/(tmp_n1 + colSums(res@haplo==2))
   tmp_maf[tmp_maf>0.5]=1-tmp_maf[tmp_maf>0.5]
   snp_sel=tmp_maf>min_maf
   if(sum(snp_sel)==res@nsnp){
    cat("No SNPs discarded\n")
   }else{
    cat(res@nsnp-sum(snp_sel)," SNPs discarded\n")
    res@haplo=res@haplo[,snp_sel]
    res@nsnp=sum(snp_sel)
    res@position=res@position[snp_sel]
    res@snp.name=res@snp.name[snp_sel]
    cat(res@nsnp," SNPs remaining\n")
   }
 }

    cat("Data consists of",res@nhap,"haplotypes and",res@nsnp,"SNPs\n")
    return(res)
}





