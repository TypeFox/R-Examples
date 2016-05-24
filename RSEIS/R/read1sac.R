read1sac<-function(fname, Iendian=1 , HEADONLY=FALSE, BIGLONG=FALSE  )
  {

    #####   given a file name with one SAC, read it in and return a list
    
 if(missing(Iendian)) { Iendian=1 }
  if(missing(HEADONLY)) {HEADONLY=FALSE }
  if(missing(BIGLONG)) { BIGLONG=FALSE}

 DATAendian =c("little", "big", "swap")
  Kendian =c(1,2,3)

 

if(is.character(Iendian))
  {
    endianVEC =Iendian
  }
  else
    {
      endianVEC = DATAendian[match(Iendian , Kendian )]
    }


theENDIAN =  endianVEC
 ##

  SAChead.names  = c("delta", "depmin", "depmax", "scale", "odelta", "b", 
     "e", "o", "a", "internal1", "t0", "t1", 
     "t2", "t3", "t4", "t5", "t6", "t7", 
     "t8", "t9", "f", "resp0", "resp1", "resp2", 
     "resp3", "resp4", "resp5", "resp6", "resp7", "resp8", 
     "resp9", "stla", "stlo", "stel", "stdp", "evla", 
     "evlo", "evel", "evdp", "unused1", "user0", "user1", 
     "user2", "user3", "user4", "user5", "user6", "user7", 
     "user8", "user9", "dist", "az", "baz", "gcarc", 
     "internal2", "internal3", "depmen", "cmpaz", "cmpinc", "unused2", 
     "unused3", "unused4", "unused5", "unused6", "unused7", "unused8", 
     "unused9", "unused10", "unused11", "unused12", "nzyear", "nzjday", 
     "nzhour", "nzmin", "nzsec", "nzmsec", "internal4", "internal5", 
     "internal6", "npts", "internal7", "internal8", "unused13", "unused14", 
     "unused15", "iftype", "idep", "iztype", "unused16", "iinst", 
     "istreg", "ievreg", "ievtyp", "iqual", "isynth", "unused17", 
     "unused18", "unused19", "unused20", "unused21", "unused22", "unused23", 
     "unused24", "unused25", "unused26", "leven", "lpspol", "lovrok", 
     "lcalda", "unused27", "kstnm", "kevnm", "khole", "ko", 
     "ka", "kt0", "kt1", "kt2", "kt3", "kt4", 
     "kt5", "kt6", "kt7", "kt8", "kt9", "kf", 
     "kuser0", "kuser1", "kuser2", "kcmpnm", "knetwk", "kdatrd", 
     "kinst")

formsac = c( rep('float' , times=70), rep('long' , times=40), rep('char' , times=23) )

 
charlen = c(rep(NA, times=110), 8, 16, rep(8, 21) )

###   data.frame(SAChead.names, formsac, charlen)

 
HEAD = vector(mode="list")
 
 if(BIGLONG)
    {
      
    ishort = 2
  iint  = 4
  ilong = 8
  ifloat = 4
  idouble = 8

  }
  else
    {

    ishort = 2
  iint  = 4
  ilong = 4
  ifloat = 4
  idouble = 8
  }

sizes = c(ilong, ishort , iint, ifloat)

m1 = match( formsac, c("long",  "short", "char",  "float") )
isize = rep(iint, length(formsac))
isize = sizes[m1]

##  cat(paste("In read1sac:", fname), sep="\n")
 
zz <- file(fname , "rb")


 ###   the sac header is 70 floats, 40 longs and the rest chars
isign = TRUE
for(i in 1:70 )
  {

           A1 =  readBin(zz, numeric() , n = 1, size = isize[i] , signed = isign,
        endian = theENDIAN)
      HEAD[[i]] = A1
               
  }
 
 ## close(zz)



for(i in 71:110 )
  {

           A1 =  readBin(zz, integer() , n = 1, size = isize[i] , signed = TRUE,
        endian = theENDIAN)
      HEAD[[i]] = A1
               
  }

for(i in 111:length(formsac) )
  {
    fchar=readChar(zz, charlen[i], useBytes = FALSE)
    HEAD[[i]] = fchar
  }

names(HEAD) = SAChead.names

 sacLIST=list(HEAD=HEAD, amp=NULL)
 
 if(HEADONLY==TRUE)
   {
     close(zz)
     invisible(sacLIST)
   }
 else
   {

 
N = HEAD[[ which(SAChead.names== "npts")  ]]
D1 = readBin(zz, numeric()  , n = N , size =ifloat ,  endian = theENDIAN, signed = TRUE)

 sacLIST$amp = D1

    close(zz)

invisible(sacLIST)
}
}


