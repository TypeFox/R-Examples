
# Rdsm worker/manager code

# author:  N. Matloff

# note:  uses the portion of the "parallel" package derived from the
# "snow" package, to be referred to here as "snow" for brevity; the
# terms "manager" and "worker" will refer to the master and cluster
# nodes, respectively

# builds on R's "parallel" and "bigmemory packages" to form a threads
# type of environment on shared-memory machines (or on clusters, using
# file storage for sharing); each worker runs one thread

# shared variables are set up in physical shared memory (or in file
# storage); due to "bigmemory" restriction, all shared variables are
# matrices; they are accessed via name, no quotes, e.g.

#    mgrmakevar(cls,"w",1000,1000)  # at manager, create variable
#    w[2,5] <- 8  # run at worker, then 8 visible at other threads and mgr.

# initializing Rdsm:

#   create a snow cluster, cls 
#   run mgrinit(cls) at the manager to initialize Rdsm

# running an Rdsm app:
#   make shared variables at the manager for the app, using mgrmakevar()
#   possibly export to cluster additional objects
#   run Rdsm app code, using clusterEvalQ() at the manager to launch
#      threads running that code

# CRAN issue; no global variables are formed at the manager, so the code
# is CRAN-compliant, but CRAN check doesn't realize this
if(getRversion() >= "2.15.1") 
   globalVariables(c("myinfo","brlock","barrlock","gbl",
      "realrdsmlock","realrdsmunlock"))

# options(bigmemory.typecast.warning=FALSE)

# *************************  INITIALIZATION  ***********************

# mgrinit() starts Rdsm

# arguments:

#    cls:  snow cluster
#    boost:  if TRUE, then use synchronicity locks, else use 
#            backing store locks; 
#    barrback:  if TRUE, store certain variables involved with the 
#               barrier in backing store rather than in memory

mgrinit <- function(cls,boost=F,barrback=F) {
   # set up so that each worker node will have a global variable myinfo
   # that contains the thread ID and number of threads
   setmyinfo <- function(i,n) {
      assign("myinfo",list(id = i,nwrkrs = n),pos=tmpenv)
   }
   ncls <- length(cls)
   parallel::clusterEvalQ(cls,tmpenv <- new.env())
   parallel::clusterApply(cls,1:ncls,setmyinfo,ncls)
   parallel::clusterEvalQ(cls,myinfo <- get("myinfo",tmpenv))
   # we create global variables only at the workers, thus OK for CRAN,
   # but CRAN check complains anyway, so here is a workaround
   parallel::clusterEvalQ(cls,gbl <- globalenv())
   # set up the requested locking type, via synchronicity or in backing store 
   if (boost) {
      rdsmlock <- boostlock
      rdsmunlock <- boostunlock
   } else {
      rdsmlock <- backlock
      rdsmunlock <- backunlock
   }
   # set up the workers to use the proper lock type
   parallel::clusterExport(cls,"rdsmlock",envir=environment())
   parallel::clusterExport(cls,"rdsmunlock",envir=environment())
   parallel::clusterEvalQ(cls,realrdsmlock <- rdsmlock)
   parallel::clusterEvalQ(cls,realrdsmunlock <- rdsmunlock)
   # send the threads needed Rdsm functions
   parallel::clusterExport(cls,"barr")
   parallel::clusterExport(cls,"getidxs")
   parallel::clusterExport(cls,"getmatrix")
   parallel::clusterExport(cls,"readsync")
   parallel::clusterExport(cls,"writesync")
   # make a single barrier, set it up on the worker nodes
   makebarr(cls,boost,barrback)
}

# ******************  CREATING SHARED VARIABLES  *********************

# mgrmakevar() is used to created shared variables; executed at manager

# arguments:

#    cls:  snow cluster
#    varname:  name of variable
#    nr, nc:  number of rows, columns in matrix
#    vartype:  "double", "integer", etc.
#    fs:  if TRUE, shared variable will be stored in the file system,
#         not shared memory
#    mgrcpy:  if TRUE, the shared variable will be visible from the mgr
#    savedesc:  if TRUE, the descriptor of the variable will be save to
#               a file, e.g. "x.desc" if varname is "x"

# the variable is created in the global spaces of the worker nodes; at
# the manager, the variable is created in the environment of the caller,
# or not at all, depending on mgrcpy

mgrmakevar <- function(cls,varname,nr,nc,vartype="double",
      fs=FALSE,mgrcpy=TRUE,savedesc=TRUE) {
  if (!fs) {
     tmp <- bigmemory::big.matrix(nrow=nr,ncol=nc,type=vartype) 
     if (savedesc)
        dput(bigmemory::describe(tmp),paste(varname,".desc",sep=""))
  } else
     tmp <- bigmemory::filebacked.big.matrix(nrow=nr,ncol=nc,
         type=vartype,backingfile=varname,
         descriptorfile=paste(varname,".desc",sep=""))
   # make accessible to manager
   if (mgrcpy) assign(varname,tmp,pos=parent.frame())  
   # get the descriptor for this big.matrix object, to send to the
   # worker nodes
   parallel::clusterExport(cls,"varname",envir=environment())
   if (!fs) {
      desc <- bigmemory::describe(tmp)
      parallel::clusterExport(cls,"desc",envir=environment())
      parallel::clusterEvalQ(cls,
         tmp <- bigmemory::attach.big.matrix(desc))
   } else {
      parallel::clusterEvalQ(cls,
         tmp <- bigmemory::attach.big.matrix(paste(varname,".desc",sep="")))
   }
   parallel::clusterEvalQ(cls,assign(varname,tmp))  
   invisible(0)
}

# **********************  LOCKS, ETC.  ***********************

# LOCKS

# these are just stubs, need for CRAN check to be replaced at runtime 
rdsmlock <- function(lck) 0
rdsmunlock <- function(lck) 0

mgrmakelock <- function(cls,lockname,boost=F) {
   if (boost) {
      require(synchronicity)
      tmp <- synchronicity::boost.mutex()
      desc <- synchronicity::describe(tmp)
      parallel::clusterEvalQ(cls,require(synchronicity))
      parallel::clusterExport(cls,"desc",envir=environment())
      parallel::clusterEvalQ(cls,
         tmp <- synchronicity::attach.mutex(desc)) 
      parallel::clusterExport(cls,"lockname",envir=environment())
      parallel::clusterEvalQ(cls,assign(lockname,tmp))  
   } else {
      # make sure lock dir not left over from a previous run
      unlink(lockname,recursive=T)
   }
}

# synchronicity lock/unlock
boostlock <- function(lck) {
   if (is.character(lck)) lck <- get(lck)
   require(synchronicity)
   synchronicity::lock(lck)  
}
boostunlock <- function(lck) {
   if (is.character(lck)) lck <- get(lck)
   require(synchronicity)
   synchronicity::unlock(lck)  
}

# general lock/unlock (writing to a shared filesystem)
# lockname must be quoted
backlock <- function(lockname) {
   repeat {
      if (dir.create(lockname)) return()
   }
}
backunlock <- function(lockname) {
   unlink(lockname,recursive=T)
}

# BARRIER  

# sense-reversing barrier implementation; see mgrinit() above regarding
# boost and barrback
makebarr <- function(cls,boost=F,barrback=F) {
   mgrmakevar(cls,"barrnumleft",1,1,"integer",fs=barrback,mgrcpy=F)
   mgrmakevar(cls,"barrsense",1,1,"integer",fs=barrback,mgrcpy=F)
   mgrmakelock(cls,"barrlock",boost)
   if (boost) {
      clusterEvalQ(cls,brlock <- barrlock)
   } else {
      clusterEvalQ(cls,brlock <- "barrlock")
   }
   # first component is the count, second is "parity";
   # barrnumleft is number left before barrier done
   clusterEvalQ(cls,barrnumleft[1] <- myinfo$nwrkrs)  
   clusterEvalQ(cls,barrsense[1] <- 0)  # sense (0 or 1)
}

# barrier op
barr <- function() {
   realrdsmlock(brlock)
   # rdsmlock(brlock)
   count <- barrnumleft[1]
   sense <- barrsense[1]
   if (count == 1) {  # all done
      # reset count
      barrnumleft[1] <- myinfo$nwrkrs 
      # reverse sense
      barrsense[1] <- 1 - barrsense[1]
      realrdsmunlock(brlock)
      # rdsmunlock(brlock)
      return()
   } else {
      barrnumleft[1] <- barrnumleft[1] - 1
      realrdsmunlock(brlock)
      # rdsmunlock(brlock)
      repeat {
         if (barrsense[1] != sense) break
      }
   }
}

# SYNC ROUTINES

# these are needed only if the shared variables are in backing store

# typical usages:

#    critical section: 

#        readers call readsync() upon entry to the
#        section, writer calls writesync() before exit

#    barrier:

#        before barrier, one or more nodes write to a variable,
#        then call writesync()

#        after barrier, one or more nodes call readsync()

#    manager:

#       typically the manager calls readsync() to get final result

# sync so can read the variable after others have written to it;
readsync <- function(varname) {
   rm(list=varname,envir=gbl)
   gc()
   tmp <- bigmemory::attach.big.matrix(paste(varname,".desc",sep=""))
   assign(varname,tmp,envir=gbl)
}

# sync so write to variable is visible by others;
# callable only at the worker nodest
writesync <- function(varname) {
   backlock("synclock")
   bigmemory::flush(getmatrix(varname))
   # rm(list=varname,envir=.GlobalEnv)
   rm(list=varname,envir=gbl)
   tmp <- bigmemory::attach.big.matrix(paste(varname,".desc",sep=""))
   # assign(varname,tmp,envir=.GlobalEnv)
   assign(varname,tmp,envir=gbl)
   backunlock("synclock")
}

# *************************  UTILITIES  ***********************

# find indices among 1:m to be handled by the given node
getidxs <- function(m) {
   parallel::splitIndices(m,myinfo$nwrkrs)[[myinfo$id]]
}

# obtain access to the matrix m
# the matrix m can be either any of the following:
#    an ordinary R matrix
#    a bigmemory matrix
#    a bigmemory matrix descriptor
#    a bigmemory matrix name, i.e. a character string
# this function determines which, and returns an R matrix (first case
#    above) or a bigmemory matrix (last 3 cases above)
# in the last case, it is assumed that the descriptor is in the file
#    m.dexc
getmatrix <- function(m) {
   # require(bigmemory)
   cl <- class(m)
   if (cl %in% c("matrix","big.matrix")) return(m)
   if (cl == "big.matrix.descriptor") {
      return(bigmemory::attach.big.matrix(m))
   }
   if (cl == "character") {
      dfilename <- paste(m,".desc",sep="")
      if (!(dfilename %in% list.files())) 
         stop(paste('file ',dfilename,"missing"))
      return(bigmemory::attach.big.matrix(dfilename))
   }
}

# to shut down cluster (important!):
stoprdsm <- function(cls) {
   stopCluster(cls)
   rm(cls)
   # clean up .desc files
   unlink("*.desc")
}

# loads the examples file exfile, from the subdir directory within the
# tree for the package pkg ; if exfile=NA, the names of the files are
# displayed, without loading
loadex <- function(pkg,exfile=NA,subdir="examples") {
   sp <- searchpaths()
   idx <- grep(pkg,sp)
   dir <-  paste(sp[idx],"/",subdir,sep="")
   if (is.na(exfile)) return(dir(dir))
   source(paste(dir,"/",exfile,sep=""))
}

