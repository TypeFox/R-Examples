#Package parallel RSM

#constant variables
mytag=list(ready=1,result=2,closed=3,doJob=4,done=5)

# Initialize MPI
#library("Rmpi")


create_slaves <- function(slave_count)
{

if(mpi.comm.size()>1)
{
  print(sprintf("! %d slaves are running already",mpi.comm.size()))
  print("If you want change number of slaves please execute mpi.close.Rslaves()")
  return();
}
# Notice we just say "give us all the slaves you've got."
mpi.spawn.Rslaves(nslaves=slave_count)

if (mpi.comm.size() < 2) {
    print("More slave processes are required.")
    mpi.quit()
    }
}

.Last <- function(){
    if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 0){
            print("Please use mpi.close.Rslaves() to close slaves.")
            mpi.close.Rslaves()
        }
        print("Please use mpi.quit() to quit R")
        #.Call("mpi_finalize",PACKAGE="Rmpi")
    }
}

send_data <- function(y,X,m,initial_weights)
{
#p <- ncol(X)
#n <- nrow(X)

# Now, send the data to the slaves
mpi.bcast.Robj2slave(y)
mpi.bcast.Robj2slave(X)
#mpi.bcast.Robj2slave(p)
#mpi.bcast.Robj2slave(m)
#mpi.bcast.Robj2slave(n)
#mpi.bcast.Robj2slave(initial_weights)
mpi.bcast.Robj2slave(mytag)

}

#------------------------------ many models Static -----------------------------

# Function the slaves will call to perform a validation on the
# fold equal to their slave number.
# Assumes: thedata,fold,foldNumber,p
fitSubmodelStatic <- function(n,p,B,m,initial_weights) {
    # Note the use of the tag for sent messages:
    #     ready=1,result=2,closed=3,doJob=4,done=5

    master_id <- 0
    junk <-0
    done <- 0
    taskNumber = -1

    # Signal being ready to receive a new task
    mpi.send.Robj(junk,master_id,mytag$ready)
    # Receive a task
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]

    if (tag == mytag$doJob) {

        taskNumber <- task$taskNumber
        print(sprintf('slave %d got job (p=%d,n=%d,B=%d) at %s',taskNumber,p,n,B, Sys.time()))
        generator <- task$generator
        .Random.seed <- generator

        taskCount <- task$taskCount

        # Create vector of final scores
        finalScores <- numeric(p)
        # Create vector containing numbers of selection of each variable
        ns <- numeric(p)

  		  for(i in 1:taskCount)
  		  {
         submodel = sample(1:p,size=m,replace=FALSE,prob=initial_weights)
  		   lm1 = lm(y~X[,submodel])
         weights = as.numeric((summary(lm1)$coef[-1,3])^2)
  			 finalScores[submodel] <-  finalScores[submodel] + weights
         ns[submodel] <- ns[submodel] + 1
  		  }

        # Send a results message back to the master
        print(sprintf('slave %d sending result at %s',taskNumber,Sys.time()))
        results <- list(taskNumber=taskNumber, partialScores=finalScores,partialNS=ns)
        mpi.send.Robj(results,master_id,mytag$result)
   }

   #print(sprintf('slave %d is exiting at %s',taskNumber, Sys.time()))
   #mpi.send.Robj(junk,master_id,mytag$closed)
   rm(list = ls())
}

make_experimentStatic<-function(n,p,B,m,initial_weights)
{

# Send the function to the slaves
mpi.bcast.Robj2slave(fitSubmodelStatic)
#mpi.bcast.Robj2slave(initial_seed)
#mpi.bcast.Robj2slave(B)

# Call the function in all the slaves to get them ready to
# undertake tasks
mpi.bcast.cmd(cmd=fitSubmodelStatic, n=n,p=p,B=B,m=m,initial_weights=initial_weights)

# Create task list
nslaves=mpi.comm.size()-1
taskCount=floor(B/nslaves)
tasks <- vector('list')

#Create parallel random generator
requireNamespace("parallel",quietly = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(runif(1))
s <- .Random.seed[1:6]


for (i in 1:(nslaves-1)) {
    #print(s)
    s=nextRNGStream(s)
    tasks[[i]] <- list(taskNumber=i,taskCount=taskCount, generator=s)
    }
s=nextRNGStream(s)
tasks[[nslaves]] <- list(taskNumber=nslaves,taskCount=B-taskCount*(nslaves-1),  generator=s)
#print(tasks)

# Create vector of final scores
finalScores <- numeric(p)
# Create vector containing numbers of selection of each variable
ns <- numeric(p)

junk <- 0
closed_slaves <- 0
n_slaves <- mpi.comm.size()-1

while (closed_slaves < n_slaves) {
    #print('waiting for connection')
    # Receive a message from a slave
    message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    message_info <- mpi.get.sourcetag()
    slave_id <- message_info[1]
    tag <- message_info[2]
    #print(tag)

    if (tag == mytag$ready) {
        # slave is ready for a task. Give it the next task, or tell it tasks
        # are done if there are none.
      	if (length(tasks) > 0) {
                  # Send a task, and then remove it from the task list
                  mpi.send.Robj(tasks[[1]], slave_id, mytag$doJob);
                  tasks[[1]] <- NULL
      	         #print(slave_id)
      	           print(sprintf('job sent to slave %d at %s',slave_id,Sys.time()))

          }
        }
    if (tag == mytag$result) {
        # The message contains results. Do something with the results.
        # Store them in the data structure
        print(sprintf('job recived from slave %d at %s',slave_id, Sys.time()))
        taskNumber <- message$taskNumber
     	  partialScores <- message$partialScores
		    partialNS <- message$partialNS
     	  finalScores <-  finalScores + partialScores
   	    ns <- ns + partialNS
   	    closed_slaves <- closed_slaves + 1
   	  }
}


# Take an average of final scores
ns=ifelse(ns!=0,ns,1)
finalScores = finalScores/ns
return(finalScores)
}


