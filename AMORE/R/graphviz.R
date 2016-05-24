graphviz.MLPnet <- function(net,filename,digits=8) {
   if (class(net)!="MLPnet") {
      stop("Your net parameter does not belong to the MLPnet class. Are you aware that the result from the train function is now a list instead of a net? Check parameters and try again");
   }
   cat(file=filename," digraph AMOREnet { \n",append=FALSE);
   cat(file=filename,"rankdir=LR;         \n",append=TRUE);
   cat(file=filename,"ordering=out;       \n",append=TRUE);
   cat(file=filename,"ranksep=2;          \n",append=TRUE);
   cat(file=filename,"nodesep=1;          \n",append=TRUE);
   for (i in 1:length(net$layers[[1]])) {
      cat(file=filename,"node [shape = hexagon, color=\"green\"] ", paste("\"Input ",i,"\"",sep=""),";\n",append=TRUE);
   }
   for (ind.neuron in 1:length(net$neurons)) {
      neuron <- net$neuron[[ind.neuron]] ;
      cat(file=filename,"node [shape = record, color=\"blue\"] ",append=TRUE); 
      cat(file=filename,neuron$id,"[label = \"{<id> Id=\\N  | { ",append=TRUE);

      for ( ind.weight in 1:length(neuron$weights) ) {
         if (neuron$input.links[ind.weight] < 0 ) {
           cat(file=filename,"wi",-neuron$input.links[ind.weight],": ",round(neuron$weights[ind.weight],digits),"|",sep="",append=TRUE);
         } else {          
           cat(file=filename,"w",neuron$input.link[ind.weight],": ",round(neuron$weights[ind.weight],digits),"|",sep="",append=TRUE);
         }
      }
      cat(file=filename,"Bias:",round(neuron$bias,digits),"}|",neuron$activation.function,"|","<v0> v0:", round(neuron$v0,digits),"} \" ];\n",append=TRUE)
   }
 for (i in 1:length(net$layers[[length(net$layers)]])) {
      cat(file=filename,"node [shape = hexagon, color=\"red\"] ", paste("\"Output ",i,"\"",sep=""),";\n",append=TRUE);
   }

   for (ind.neuron in 1:length(net$neurons)) {
      neuron <- net$neurons[[ind.neuron]];
      for ( ind.weight in 1:length(neuron$weights)) {
         if (neuron$input.links[ind.weight] < 0 ) {
            cat(file=filename,"\"Input ",-neuron$input.links[ind.weight],"\" -> ",neuron$id," ;\n", sep="",append=TRUE);
         } else {
            cat(file=filename,neuron$input.links[ind.weight]," -> ",neuron$id," ;\n", sep="",append=TRUE);
         }
      }
      if (neuron$type=="output") {
         cat(file=filename,neuron$id," -> \"Output ",neuron$output.aims,"\" ;\n", sep="",append=TRUE);
      }





   }
cat(file=filename,"}\n",append=TRUE);


}
