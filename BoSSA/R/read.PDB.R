`read.PDB` <-
function(X)
{
  pdb <- scan(X,sep="\n",what="complex")
  category <- gsub(" ","",unlist(lapply(pdb,substring,1,6)))
  
  header <- unlist(lapply(pdb[category=="HEADER"],substring,11))
  
  compound <- unlist(lapply(pdb[category=="COMPND"],substring,11))
  id <- as.numeric(gsub(";","",gsub("MOL_ID:","",compound[grep("MOL_ID: ",compound)])))
  ch <- gsub(" ","",(gsub(";","",gsub("CHAIN:","",compound[grep("CHAIN: ",compound)]))))
  mol <- gsub(",","",gsub("MOLECULE:","",compound[grep("MOLECULE: ",compound)]))
  compound <- data.frame(id=id,chain=ch,molecule=mol)
  
  atomt <- unlist(lapply(pdb[category=="ATOM"],substring,14,16))
  aminoa <- as.character(unlist(lapply(pdb[category=="ATOM"],substring,18,20)))
  chain <- unlist(lapply(pdb[category=="ATOM"],substring,22,22))
  namino <- as.numeric(unlist(lapply(pdb[category=="ATOM"],substring,23,26)))
  Xcoor <- as.numeric(unlist(lapply(pdb[category=="ATOM"],substring,32,38)))
  Ycoor <- as.numeric(unlist(lapply(pdb[category=="ATOM"],substring,39,46)))
  Zcoor <- as.numeric(unlist(lapply(pdb[category=="ATOM"],substring,47,54)))
  atom <- data.frame(atom=atomt,aa=aminoa,chain=chain,naa=namino,X=Xcoor,Y=Ycoor,Z=Zcoor)
  
  sequence <- list(NULL)
  uch <- unique(atom$chain)
  for(i in 1:length(uch))
  {
    j <- i*2
    subatom <- atom[atom$chain==uch[i],]
    uaa <- unique(subatom$naa)
    sequence[[j-1]] <- uaa   
    sequence[[j]] <- as.character(subatom$aa[unlist(lapply(uaa,FUN <- function(X){grep(X,subatom$naa)[1]}))])
    names(sequence)[j-1] <- paste("ref_",uch[i],sep="")
    names(sequence)[j] <- paste("chain_",uch[i],sep="")
  }
  
  out <- list(header)
  out[[2]] <- compound
  out[[3]] <- atom
  out[[4]] <- sequence
  names(out) <- c("header","compound","atom","sequence")
  out
}

