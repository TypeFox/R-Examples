SampleTrees<-function(trees, burnin, final.number, format, prefix){           #nexus trees, %eg .25, 500
  trees<-read.nexus(trees)
  original.number<-as.numeric(length(trees))                           #NUMBER OF TREES IN ORIGINAL FILE
  post.burnin.trees<-trees[(burnin*original.number):original.number]   #THIS CREATES THE POST BURNIN PORTION
  final.trees<-sample(post.burnin.trees, final.number)                 #THIS DOWNSAMPLES THE COLLECTION OF TREES
  if(format=="new"){write.tree(final.trees, file=paste(prefix,".nwk"))
                    print("Your trees were saved in the Newick format")}         #THIS SAVES THEM AS NEWICK FORMAT
  if(format=="nex"){write.nexus(final.trees, file=paste(prefix,".nex"))
                    print("Your trees were saved in the Nexus format")}        #THIS SAVES THEM AS NEXUS FORMAT

}