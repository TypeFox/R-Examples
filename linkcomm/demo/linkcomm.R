# Demonstration of linkcomm.

#require(linkcomm)

cat(paste(c("\n\n",rep("*",33),"\n*\tDemo of linkcomm\t*\n",rep("*",33),"\n\n"),collapse=""))

contFunc <- function()
	{
	cat("Press \"return\" to continue, \"m\" for the main menu, or \"q\" to quit:")
	cont <- readline()
	if(cont == "q"){
		return(TRUE)
	}else if(cont == "m"){
		linkcomm.demo.main()
		return(TRUE)
	}else{	return(FALSE)
		}
	}


linkcomm.demo.main <- function()
	{
	
	cat("\nExample networks to use:\n\n1. A set of 56 yeast proteins involved in 651 interactions related to transcription.\n2. The co-appearance network for Les Miserables.\n3. A social network of friendships between 34 members of a karate club at a US university in the 1970s.\n\n")

	choices <- c("1","2","3","q")
	choose <- NULL
	
	while(length(match(choose,choices)) == 0){
		cat("Enter the number of the network you wish to work with and press \"return\" or type \"q\" to quit:\n")
		choose <- as.character(readline())

		if(choose == "q"){
			return(invisible())
		}else if(choose == "1"){
			x <- pp_rnapol
			name <- "pp_rnapol"
		}else if(choose == "2"){
			x <- lesmiserables
			name <- "lesmiserables"
		}else if(choose == "3"){
			x <- karate
			name <- "karate"
			}
		}

	linkcomm.network(x = x, name = name)

	}
		

linkcomm.network <- function(x, name)
	{
	# x is the user-chosen network.

	cat(paste(c("\nThe input network is arranged as an edge list:\n\n> head(",name,")\n\n"),collapse=""))

	print(head(x))

	cat(paste(c("\nGet link communities for this network...\n\n> lc <- getLinkCommunities(",name,", hcmethod=\"average\")\n\n"),collapse=""))
	
	if(contFunc()){return(invisible())}

	lc <- getLinkCommunities(x)

	cat("\nDisplay community membership for the top 20 nodes that belong to the most communities...\n\n> plot(lc, type = \"members\")\n\n")

	if(contFunc()){return(invisible())}

	plot(lc,type="members")

	cat("\nDisplay the network with edges coloured according to community membership...\n\n> plot(lc, type = \"graph\")\n\n")

	if(contFunc()){return(invisible())}

	plot(lc,type="graph")

	cat("\nFind and display a subnetwork where the nodes of one link community are entirely nested within another link community...\n\n> getNestedHierarchies(lc, clusid = 1, plot = TRUE)\n\n")

	if(contFunc()){return(invisible())}

	getNestedHierarchies(lc,clusid=1)

	cat("\nPlot the network with a Spencer circle layout...\n\n> plot(lc, type = \"graph\", layout=\"spencer.circle\")\n\n")

	if(contFunc()){return(invisible())}

	plot(lc,type="graph",layout="spencer.circle")

	cat(paste(c("\nPlot Spencer circle for the top-connected node...\n\n> plot(lc, type = \"graph\", nodes = \"",names(lc$numclusters[1]),"\", layout=\"spencer.circle\", vertex.label.cex=0.8, jitter = 0.2))\n\n"),collapse=""))

	if(contFunc()){return(invisible())}

	plot(lc, type = "graph", nodes = names(lc$numclusters[1]), layout="spencer.circle")

	cat("\nDisplay the top modular networks...\n\n> plot(lc, type = \"commsumm\", summary = \"mod\")\n\n")

	if(contFunc()){return(invisible())}

	dev.off()

	plot(lc,type="commsumm",summary="mod")
	
	# Return to start.
	linkcomm.demo.main()

	}

linkcomm.demo.main()


