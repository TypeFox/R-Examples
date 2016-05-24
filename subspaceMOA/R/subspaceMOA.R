#'@name subspaceMOA
#'@title Subspace and Projected stream clustering
#'
#'@description The \emph{subspaceMOA} package is an extension of the stream
#'package, focusing on stream clustering for high-dimensional data streams.
#'
#'To this end, two new data streams are provided: 
#'\link{DSD_RandomRBFSubspaceGeneratorEvents}, a synthetic data stream for high 
#'dimensional data with clusters that only exist in certain subspaces of the data and 
#'\link{DSD_SubspaceARFFStream}, a utility that reads ARFF files with subspace 
#'information.
#'
#'New subspace stream clustering algorithms for high-dimensional data can be
#'constructed using \link{DSC_ThreeStage} in combination with one of 
#'\link{DSC_subspaceDenStream} and \link{DSC_subspaceCluStream} and one of 
#'\link{DSC_clique}, \link{DSC_p3c}, \link{DSC_proclus}, \link{DSC_subclu}.
#'
#'Additionally, implementations of the existing subspace stream clustering
#'algorithms \link{DSC_PreDeConStream} and \link{DSC_HDDStream} are present.
#'
#'Lastly, functions for interactive visualization with the \link{shiny} package
#'are also included: A specific clustering can be plotted using
#'\link{plot_stream_interactive} and the clustering can also be animated with
#'the function \link{animate_stream_interactive}.
#'
#'
NULL