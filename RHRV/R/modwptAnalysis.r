 modwptAnalysis=function(x,wavelet="d4", ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05, LFmin = 0.05, LFmax = 0.15, HFmin = 0.15, HFmax = 0.4, sampling=4,bandtolerance,relative)
{ 

  # Auxiliary variables
  ULFpower=0;
  VLFpower=0;
  LFpower=0;
  HFpower=0;
  
  # Working with a zero mean signal
  x=x-mean(x)
  # Getting tree nodes
  ULFnodes=getNodes(ULFmin,ULFmax,sampling,bandtolerance,relative)
  VLFnodes=getNodes(VLFmin,VLFmax,sampling,bandtolerance,relative)
  LFnodes=getNodes(LFmin,LFmax,sampling,bandtolerance,relative)
  HFnodes=getNodes(HFmin,HFmax,sampling,bandtolerance,relative)
  #All nodes 
  targets=unlist(c(ULFnodes[2:length(ULFnodes)],VLFnodes[2:length(VLFnodes)],LFnodes[2:length(LFnodes)],HFnodes[2:length(HFnodes)]))


  # Get the maximum depth of the nodes
  depth=max(ULFnodes[[1]],VLFnodes[[1]],LFnodes[[1]],HFnodes[[1]])
  # Wavelet analysis
  wx=BoundModwpt(x,wf=wavelet,n.levels=depth,targets);

  # Power calculation             
  ULFpower=getPower(wx,ULFnodes,wavelet)
  VLFpower=getPower(wx,VLFnodes,wavelet)
  LFpower=getPower(wx,LFnodes,wavelet)
  HFpower=getPower(wx,HFnodes,wavelet)


  # Store the results of the calculation
  return(list(ULF=ULFpower,VLF=VLFpower,LF=LFpower,HF=HFpower,depth=depth))
  
}
