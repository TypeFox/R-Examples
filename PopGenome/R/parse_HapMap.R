parse_HapMap <- function(filepath){

cat("\n")
# first block to get basic informations

block         <- scan(file=filepath,nlines=500,what=character())
sub           <- substr(block,1,2)

# ind         <- sub=="NA" # FIXME 
startind      <- grep("QC",block)[1]
endind        <- grep("rs",block[2:length(block)])[1] # 2 da res auch ganz am Anfang, [1] nur erstes rs !

pp            <- (startind + 1):(endind)

# first position
# pp              <- which(ind)
pp_x            <- pp[length(pp)] 
first_pos       <- block[pp_x+4]
first_ref       <- block[pp_x+1]
first_type      <- block[pp_x+2]
first_chr       <- block[pp_x+3]

# print(first_pos)
# print(first_ref)
# print(first_type)
# print(first_chr)

# individuals
# individuals   <- block[ind]
individuals   <- block[pp]
# print(individuals)
n.individuals <- length(individuals)

# first block
# num.nucs      <- nchar(block)
# nuc           <- match(2,num.nucs)

nuc <- grep("^[A-Z]{2,2}$",block)[1]

first.block   <- block[nuc:(nuc+(n.individuals-1))]
sub_sub       <- block[(nuc+n.individuals):length(sub)]

 # print(length(first.block))
 # print(length(individuals))
 # print(sub_sub)
 # stop("")
# intermediate stuff
# stuff         <- match(2,nchar(sub_sub))

stuff         <- grep("^[A-Z]{2,2}$",sub_sub)[1]
inter         <- stuff - 1


# ---------------------------------------------------

# BLOCK       <- scan(file=filepath,skip=(nuc+(n.individuals-1)),what=character())
BLOCK         <- scan(file=filepath,what=character())
#print(BLOCK[(nuc+(n.individuals))])
#stop("")
BLOCK         <- BLOCK[(nuc+(n.individuals)):length(BLOCK)]


scan          <- seq(1,length(BLOCK),by=inter+n.individuals)     
Matr          <- matrix(,n.individuals,length(scan))

yy            <- 1

# init
  init           <- vector(,length(scan))
  ref            <- init
  type           <- init
  chr            <- init
  pos            <- init # ff(0,length(scan),vmode="double")

  cat("Scan Data ... \n")

for(xx in scan){
   #print(BLOCK[xx])
   ref[yy]          <- BLOCK[xx]
   type[yy]         <- BLOCK[xx+1]
   chr[yy]          <- BLOCK[xx+2]

 # print(xx)  

   pos[yy]      <- BLOCK[xx+3]
   blockx       <- BLOCK[(xx+inter):(xx+n.individuals+inter-1)] 

   # print(blockx)
   Matr[,yy]    <- blockx
   yy           <- yy + 1
}

 cat("Create Matrix ... \n")

# add first Block
Matr            <- cbind(first.block,Matr)
pos             <- as.numeric(c(first_pos,pos))
ref             <- c(first_ref,ref)
type            <- c(first_type,type)
chr             <- c(first_chr,chr)


MAT               <- apply(Matr,2,function(x){
                r <- unlist(strsplit(x,split=""))
		return(r)
})

 cat("Coding ... \n")

MAT <- .Call("code_nucs",MAT)           



# modify individuals
ind           <- sapply(individuals,function(x){return(c(x,paste(x,".2",sep="")))})


cat("Calculation ... \n")

# oft sind positionen doppelt im File !!!???
dupids    <- which(duplicated(pos))
if(length(dupids)>0){
 MAT       <- MAT[,-dupids]
 pos       <- pos[-dupids]
 ref       <- ref[-dupids]
}
 
rownames(MAT) <- as.vector(ind)
return(list(matrix=MAT,positions=pos,reference=ref,type=type,chr=chr,dupids=dupids))
#return(list(matrix=ff(MAT,dim=dim(MAT)),positions=ff(pos),reference=ref,type=type,chr=chr))

}



