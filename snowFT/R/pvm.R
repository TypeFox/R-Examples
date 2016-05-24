#
# PVM Implementation
#

# This is commented out because the rpvm package was removed from CRAN

# recvOneDataFT.PVMcluster <- function(cl,type='b',time=0) {

    # rtag <- findRecvOneTag(cl, -1)

    # if (type == 'n') { # non-blocking receive
      # recv <- .PVM.nrecv(-1, rtag)
      # if (recv <= 0) return (NULL)
    # } else {
      # if (type == 't')  {# timeout receive
        # recv <- .PVM.trecv(-1, rtag, sec=time)
        # if (recv <= 0) return (NULL)
      # } else  # blocking receive
      # recv <- .PVM.recv(-1, rtag)
    # }
    # binfo <- .PVM.bufinfo(recv)
    # for (i in seq(along = cl)) {
      # if (cl[[i]]$tid == binfo$tid) {
        # n <- i
        # break
      # }
    # }
    # data <- .PVM.unserialize()
    # list(node = n, value = data)
  # }


# addtoCluster.PVMcluster <- function(cl, spec, ...,
                                    # options = defaultClusterOptions) {
  # options <- addClusterOptions(options, list(...))
  # n <- length(cl)
  # newcl <- vector("list",n+spec)
  # for (i in seq(along=cl))
    # newcl[[i]] <- cl[[i]]
  # for (i in (n+1):(n+spec)) {
    # newcl[[i]] <- newPVMnode(options = options)
    # newcl[[i]]$replic <- 0
  # }
  # class(newcl) <- class(cl)
  # newcl
# }

# repairCluster.PVMcluster <- function(cl, nodes, ...,
                                     # options = defaultClusterOptions) {
  # newcl <- vector("list",length(cl))
  # options <- addClusterOptions(options, list(...))
  # for (i in seq(along=cl)) {
    # if (length(nodes[nodes == i])>0) {
      # stopNode(cl[[i]])
      # newcl[[i]] <- newPVMnode(options = options)
      # clusterEvalQpart(newcl,i,require(snowFT))
      # newcl[[i]]$replic <- 0
    # } else 
    # newcl[[i]] <- cl[[i]]
  # }
  # class(newcl) <- class(cl)
	
  # newcl  
# }

# is.manageable.PVMcluster <- function(cl) {
	# return (c(cluster.size=TRUE, monitor.procs=TRUE, repair=TRUE))
# }

# processStatus.PVMnode <- function(node) {
  # status <- .PVM.pstats(node$tid)
  # if (status != "OK")
    # return(FALSE)
  # return(TRUE)
# }

# getNodeID.PVMnode <- function(node) {
  # return(node$tid)
# }

# do.administration.PVMcluster <- function(cl, clall, d, p, it, n, manage, mngtfiles, 
									# x, frep, freenodes, initfun, 
									# gentype, seed, ft_verbose, ...) {
	# free.nodes <- FALSE
	# if (length(d) <= 0) { # no results arrived yet
    	# while (TRUE) {
            # # do the administration in the waiting time
            # # ***************************************
            # if(manage['repair']) {
            # mfn <- findFailedNodes(cl) # look for failed nodes
            # if (nrow(mfn) > 0) { # failed nodes found
            	# fn<-mfn[,1]
                # frep <- c(frep, mfn[(mfn[,2]>0),2]) # only nodes where
                                        # # computation is running
				# if (ft_verbose) {
					# cat("   Failed slaves detected: ", fn,"\n")
					# cat("   Repair cluster ...\n")
				# }
                # cl <- repairCluster(cl,fn)
		
                # if (!is.null(initfun)){
		  			# if (ft_verbose) 
						# cat("   calling initfun ...\n")
                  	# clusterCallpart(cl,fn,initfun)
				# }
                # if (gentype != "None"){
		  			# if (ft_verbose) 
						# cat("   initializing RNG ...\n")
                  	# resetRNG(cl,fn,length(x),gentype,seed)
				# }
                # #keep all nodes in case
                # #messages from failed nodes arrive later
                # clall<-combinecl(clall,cl[fn])               
                # freenodes<-c(freenodes,fn)
                # if (it <= n) break # exit the loop only if there are
                                   				# # comput. to be started 
			# }
			# }
			# # read p from a file 
			# updated.values <- manage.replications.and.cluster.size(cl, clall, p, n, manage, mngtfiles, 
									# freenodes, initfun, gentype, seed, ft_verbose=ft_verbose)
			# newp <- updated.values$newp
			# if (updated.values$cluster.increased) {
				# p <- updated.values$p
				# cl <- updated.values$cl
				# clall <- updated.values$clall
				# freenodes <- updated.values$freenodes
				# break
			# }
            # p <- newp              
            # d <- recvOneResultFT(clall,'t',time=5) # wait for a result for
                                                  # # 5 sec
            # if (length(d) > 0) break # some results arrived, if not
                            	     # # do administration again
		# }  # end of while loop ****************************
        # if ((length(freenodes) > 0) && (it <= n)) free.nodes <- TRUE
	# } 
	
	# return(list(cl=cl, clall=clall, frep=frep, freenodes=freenodes, p=p, d=d, is.free.node=free.nodes))
# }