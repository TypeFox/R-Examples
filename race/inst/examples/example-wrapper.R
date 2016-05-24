# ------------------------------------------------------------------------- #
# example-wrapper.R                              An example of wrapper file #
# ------------------------------------------------------------------------- #
                                                                             
# ========================================================================= #
# Racing methods for the selection of the best                              #
# ------------------------------------------------------------------------- #
# Copyright (C) 2003 Mauro Birattari                                        #
# ========================================================================= #
# This program is free software; you can redistribute it and/or modify it   #
# under the terms of the GNU General Public License as published by the     #
# Free Software Foundation; either version 2 of the License, or (at your    #
# option) any later version.                                                #
#                                                                           #
# This program is distributed in the hope that it will be useful, but       #
# WITHOUT ANY WARRANTY; without even the implied warranty of                #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         #
# General Public License for more details.                                  #
#                                                                           #
# You should have received a copy of the GNU General Public License along   #
# with this program; if not, write to the Free Software Foundation, Inc.,   #
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.                  #
# ========================================================================= #

# ========================================================================= #
# Mauro BIRATTARI                                                           #
# IRIDIA - ULB, CP 194/6                                                    #
# Av. F. D. Roosevelt 50                                    mbiro@ulb.ac.be #
# 1050 Brussels, Belgium                     http://iridia.ulb.ac.be/~mbiro #
# ========================================================================= #

race.init<-function(){
  # Load libraries and data
  library(nnet) 
  data(iris)    

  # k-fold Cross-validation
  k<-15 
  n<-nrow(iris)
  # Set the random seed so that, if running under PVM, all the 
  # slaves use the same data partition for the cross-validation
  if (exists(".Random.seed"))
    save.seed<-.Random.seed
  set.seed(0)
  smpl<-sample(n)  
  # Restore the random seed
  if(exists("save.seed",inherits=FALSE))
    assign(".Random.seed",save.seed,envir=globalenv())
  else
    rm(".Random.seed",envir=globalenv())
  # As if nothing appened...
  msk<-k*c(0:(ceiling(n/k)-1))

  # Prepare a data.frame defining the candidates
  candidates<-expand.grid(size=c(1:10),decay=c(1:5)*1e-4)

  # Return init list
  return(list(no.candidates=nrow(candidates),
              no.tasks=k,
              msk=msk,
              smpl=smpl,
              iris=iris,
              candidates=candidates))
}


race.info<-function(data)
  return(list(race.name="Model Selection: nnet on iris",
              no.candidates=(data$no.candidates),
              no.tasks=data$no.tasks,
              extra=paste("This is a very simple example that shows",
                "how to use the 'race' package.",
                "By using a",
                paste(data$no.tasks,"-fold",sep=""),
                "cross-validation,",
                "we optimize the structure of a neural network.",
                "In this example, we optimize the number of nodes",
                "and the decay.  See the help of function 'nnet'",
                "in the homonymous package.",
                "We limit here the number of iteration of the",
                "training process for didactical reasons.")))


race.wrapper<-function(candidate,task,data){
  # We should check the parameters are of expected type,
  # are not out of range, etc... but this is just a quick example :)

  ### Evaluate the given candidate on the given task:
  # Index of examples to hold out
  idx<-data$smpl[task+data$msk]

  # Train the network
  training.set<-data$iris[-idx,]
  # nnet is too vocal :)
  # We better redirect output to /dev/null for a while
  # ONLY UNDER UNIX...
  if (.Platform$OS.type=="unix"){
    dev.null<-file("/dev/null",open="w");sink(dev.null)}
  nn<-nnet(Species~.,data=training.set,maxit=10,
           size=data$candidates[candidate,"size"],
           decay=data$candidates[candidate,"decay"])
  # ONLY UNDER UNIX...
  if (.Platform$OS.type=="unix"){
    sink(NULL);close(dev.null)}
  
  # Test the network on hold-out examples
  unseen.input<-data$iris[idx,1:4]
  target.output<-as.character(data$iris$Species[idx])
  predicted.output<-predict(nn,unseen.input,type="class")

  # Return error rate
  return(mean(predicted.output!=target.output))
}


race.describe<-function(candidate,data)
  return(data$candidates[candidate,])

