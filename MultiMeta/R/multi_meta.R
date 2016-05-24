#' Meta-analysis of multivariate GWAS results
#' 
#' \code{multi_meta} returns the meta-analysis results for multivariate GWAS across different cohorts.
#' 
#' This function applies an inverse-variance based method to meta-analyse multivariate GWAS results.
#' In particular, given \emph{n} different cohorts, for which \emph{p} phenotypes have been tested for genome-wide 
#' association, the results for each cohort will have \emph{p} different effect size coefficients i.e. beta values (one per each 
#' phenotype)
#' and a variance/covariance \emph{pxp} matrix representing beta's variances and covariances. In particular, the function is built to consider the output
#' from the GEMMA software multivariate association testing. If your output is not produced with GEMMA, the function works on
#' any results file containing the following column names: 
#' \itemize{
#' \item{\strong{chr} Chromosome} 
#' \item{\strong{ps} Position}
#' \item{\strong{rs} SNP name} 
#' \item{\strong{allele1} Effect allele} 
#' \item{\strong{allele0} Non-effect allele} 
#' \item{\strong{af} Effect-allele frequency} 
#' \item{\strong{beta_1}, \strong{beta_2}, ..., \strong{beta_p} Effect sizes for each of the \emph{p} traits} 
#' \item{\strong{Vbeta_1_1}, \strong{Vbeta_1_2}, ..., \strong{Vbeta_1_p}, \strong{Vbeta_2_2}, ...,
#'  \strong{Vbeta_2_p}, ..., \strong{Vbeta_p_p} 
#' variance-covariance matrix entries (diagonal and upper triangle values only, since this matrix is symmetric) 
#' }
#' }
#'
#' The function divides input files into chunks based on position. Only one chunk at a time is read and analysed;
#' thus a limited amount of data is loaded in the workspace at one given time. Default chunk dimension is 5 Mb for which
#' low memory is required (<250 MB for 2 cohorts). If you have larger RAM availability, sparse markers or a limited
#' number of cohorts, change chunks' dimension from the command line. 
#' 
#'@param files A vector containing the names of the results files to meta-analyse. These can be outputs from GEMMA multivariate analysis
#'or similar (see \strong{Details}). Furthermore they can be single-chromosome or genome-wide results.
#'@param N A vector containing sample sizes for each of the above files. This parameter is optional and is only
#'required for computing the overall allele frequency. 
#'@param output.file The name of the output file.
#'@param size.chunks Size of each chunk to be read and processed. Default is 5,000,000 (5 Mb). This size
#'will require very low memory usage. Increase this parameter if more memory is allocated or
#'if the number of cohorts is limited. Read more about the chunks in \strong{Details}.
#'@param min.pop Minimum number of populations required per SNP to compute meta-analysis. Default is 2, it can be any number up to the total number
#'of cohorts analysed.
#'@param sep Separator for reading input files. 
#'  
#'@examples
#'file1=system.file("extdata", "Example_file_1.txt", package="MultiMeta")
#'file2=system.file("extdata", "Example_file_2.txt", package="MultiMeta")
#'multi_meta(files=c(file1,file2), N=c(1200,600), sep=" ", 
#'output.file="Output_from_running_example.txt")
#'
#'
#'@export
#'
#'



multi_meta=function(files=c(), N=c(), output.file="Meta_Results.txt",size.chunks=5000000,min.pop=2,sep="\t")
{
    # ausilliary function to compute combined effect size coefficients (betas) and relevant variance
  mvt_mod=function (y, cov) {
    p <- dim(y)[1]
    n <- dim(y)[2]
    cov_i <- array(rep(NA, p * p * n), c(p, p, n))
    tmpy <- matrix(rep(NA, p * n), p, n)
    for (i in 1:n) {
      cov_i[, , i] <- solve(cov[, , i])
      tmpy[, i] <- cov_i[, , i] %*% unlist(y[, i])
    }
    covbeta_i <- apply(cov_i, 1:2, sum)
    covbeta <- solve(covbeta_i)
    beta <- as.vector(covbeta %*% apply(tmpy, 1, sum))
    return(list(beta = beta, cov = covbeta, cov_1=covbeta_i))
  }
  
  NF=length(files)
  if(length(files)!=length(N))
    stop("Error: The arguments have different length")
  dir=getwd()
  ####checking input files 
  for(n in 1:NF){
    if(length(grep("/", files[n])) ==0) files[n]=paste(dir, "/", files[n], sep="")
    he=scan(files[n], nlines=1, "char", quiet=T)
    if(!("chr" %in% he)) stop("Chromosome column missing or with wrong name, please check and rename \"chr\"")
    if(!("ps" %in% he)) stop("Position column missing or with wrong name, please check and rename \"ps\"")
    if(!("rs" %in% he)) stop("SNP column missing or with wrong name, please check and rename \"rs\"")
    if(!("allele1" %in% he)) stop("Effect allele column missing or with wrong name, please check and rename \"allele1\"")
    if(!("allele0" %in% he)) stop("Non-effect allele column missing or with wrong name, please check and rename \"allele0\"")
    if(!("af" %in% he)) stop("Effect allele frequency column missing or with wrong name, please check and rename \"af\"")
    if(length(grep("beta", he))==0) stop("Beta columns missing or with wrong name, please check manual and rename")
  }
  chr.ind=which(he=="chr")
  ps.ind=which(he=="ps")
  con=pipe(paste("cut -d \"",sep,"\" -f ",chr.ind," ",files[n],sep=""))
  chromosomes=unique(scan(con, "char", quiet=T))
  close(con)
  chromosomes=chromosomes[-1]
  
  #compute the number of phenotypes from the header line
  header=scan(files[1], nlines=1, "char", quiet=T)
  betas=grep(pattern="beta", header)
  L=length(betas)
  nphen=(-3+sqrt(9+8*L))/2
  he.ord=header[betas][order(header[betas])]
  
  to.chunk=ceiling(250000000/size.chunks)
  
  for(chr in chromosomes){  
    
    #cycle through the chunks
    for(k in c(1:to.chunk)){
      if(k==1){
        head_res=c("chr", "SNP", "Position", "allele1", "allele0","tot_af","n_pops","pops",he.ord, "p_value" )
        write.table(t(head_res),output.file, row.names=F, col.names=F, quote=F)
      }
      #read and assign names to files for chunk k
      
      j_1=(k-1)*size.chunks
      j_2=(k)*size.chunks-1
      for(fl in 1:NF){
        con=pipe(paste("awk '{if($",chr.ind,"==",chr," && ",j_1,"<=$",ps.ind," && $",ps.ind,"<",j_2,") print $0}' " ,files[fl],sep=""))
        assign(paste("res",fl,sep=""), try(read.table(con, stringsAsFactors=F), silent=T))
      }
      
      #check if all the files are non-empty for chunk k
      prova=0
      for(cl in 1:NF){
        if(class(get(paste("res",cl,sep=""))) != "try-error" ){ 
          prova=prova+1
        }
      }
      
      if(prova==NF){  #put all the tables into a list for simplicity
        res_list=list()
        snp=c()
        for(h in 1:NF){
          res_list[[h]]=get(paste("res",h,sep=""))  
          names(res_list[[h]]) <- scan(files[h], nlines=1, "char", quiet=T)
          res_list[[h]]=res_list[[h]][!duplicated(res_list[[h]]$rs),]
          row.names(res_list[[h]])=res_list[[h]]$rs
          snp=c(snp, res_list[[h]]$rs)    #create the vector of all the snps to be analyzed
        }
      }else{next}
      
      snp=unique(snp)
      #cycle over snps
      for(j in snp){
          count=c()
          for(r in 1:NF){
            assign(paste("j",r,sep=""), which(res_list[[r]]$rs==j))
            if(length(get(paste("j",r,sep=""))) > 0){
              count=c(count, 1)
            }else{
              count=c(count, 0)
            }
          }
         
          #check if the snp j is present in at least min.pop cohorts and which 
          tot_c=sum(count)
          good=which(count==1)
          if(tot_c >= min.pop){
          all=c()
          for(q in good){
            all=c(all, res_list[[q]]$allele1[get(paste("j",q,sep=""))] )
          }
          #switching the alleles and the betas
          which(count==1)[1]->pio
          if(dim(table(all))>=2){
            eff_all=res_list[[pio]][j ,"allele1"]
            alt_all=res_list[[pio]][j ,"allele0"]
            for(i in good[-1]){
              if(res_list[[i]][j, "allele1"]!=eff_all & res_list[[i]][j, "allele0"]==eff_all & res_list[[i]][j, "allele1"]==alt_all){
                res_list[[i]][j,"allele0"]=alt_all
                res_list[[i]][j, "allele1"]=eff_all
                res_list[[i]][j, paste("beta_",c(1:nphen),sep="")]=(-1)*res_list[[i]][j, paste("beta_",c(1:nphen),sep="")]
              }
            }
            all=c()
            for(q in good){
              all=c(all, res_list[[q]]$allele1[get(paste("j",q,sep=""))] )
            }
            if(dim(table(all))>=2){
              tab_strand=c("A"="T", "C"="G", "G"="C", "T"="A")
              rev_eff_all=tab_strand[eff_all]
              rev_alt_all=tab_strand[alt_all]
              if(res_list[[i]][j, "allele0"]==rev_alt_all & res_list[[i]][j,"allele1"]==rev_eff_all){
                res_list[[i]][j,"allele0"]=alt_all
                res_list[[i]][j, "allele1"]=eff_all
              }
              if(res_list[[i]][j, "allele0"]==rev_eff_all & res_list[[i]][j,"allele1"]==rev_alt_all){
                res_list[[i]][j,"allele0"]=alt_all
                res_list[[i]][j, "allele1"]=eff_all
                res_list[[i]][j, paste("beta_",c(1:nphen),sep="")]=(-1)*res_list[[i]][j, paste("beta_",c(1:nphen),sep="")]
              }
              all=c()
              for(q in good){
                all=c(all, res_list[[q]]$allele1[get(paste("j",q,sep=""))] ) 
              }
            }
          }
          if(dim(table(all))>1) {
            warning(paste("Check the alleles for snp ",j,sep=""))
            next
          }
          alt_all=c()
          for(q in good){
            alt_all=c(alt_all, res_list[[q]]$allele0[get(paste("j",q,sep=""))] ) 
          }
          #building the matrix of the beta values across cohorts and the array of the var/cov matrices
          if(dim(table(alt_all))==1 & table(all)==tot_c){
            Y=matrix(1, nrow=nphen, ncol=tot_c)
            cov=array(1, c(nphen,nphen,tot_c))
            for(q in 1:NF){
              if(count[q]==1){
                ind=sum(count[1:q])
                Y[,ind]=unlist(res_list[[q]][get(paste("j",q,sep="")), paste("beta_",c(1:nphen),sep="")])   
                cov[,,ind][lower.tri(cov[,,ind], diag=T)]=unlist(res_list[[q]][j, he.ord[-c(1:nphen)]])
                for(i in 1:nphen)  cov[,,ind][i,]=cov[,,ind][,i]
              }
            }
            #computing the meta-anlysis beta and variance/covariance matrix
            fin=mvt_mod(Y, cov)
            norm=suppressMessages(sqrtm(fin$cov_1) %*% fin$beta)
            zed=t(norm) %*% norm
            p_val=pchisq(as.numeric(zed), df=nphen, lower.tail=FALSE)
            res_final=c()
            tot_ac=0
            part_N=0
            for(g in 1:NF){
              if(count[g]==1){
                tot_ac=tot_ac+(res_list[[g]][j, "af"]*N[g])
                part_N=part_N+N[g]
              }
            }
            #writing the results
            tot_af=tot_ac/part_N
            res_final=unlist(c(res_list[[pio]][j, c("chr","rs","ps")], res_list[[pio]][j,"allele1"], res_list[[pio]][j,"allele0"], tot_af, tot_c, paste(count,collapse="") , fin$beta[1:nphen], fin$cov[lower.tri(fin$cov, diag=T)], p_val[1]))
            write.table(t(res_final), output.file, row.names=F, col.names=F, quote=F, append=T)
          }
        }
      }
    }
  }
}  
  
  