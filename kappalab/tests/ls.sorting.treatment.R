        library(kappalab)
	
	## n : number of criteria, here 4
	## k : search for a k-additive solution
	## d : minimal distance between 2 classes
	## t : number of prototypes        	
	## n.var.alt.A : number of elements of A

	## generate a random problem with "n.var.alt" alternatives and 4 criteria
	## n.var.alt <- 30 ## alternatives
	k <- 4
	d <- 0.1
	n.var.alt <- 10
	n.var.alt.A <- 10
  	n <- 4 ## criteria

	print("Number of prototypes: ")
	print(n.var.alt)
	print("Number of criteria: ")
	print(n)
	print("Number of elements of A: ")
	print(n.var.alt.A)
        print("Epsilon: ")
        print(d)
        print("k: ")
        print(k)
        
	print("*** Generating the data for the prototypes")
  	P <- matrix(runif(n.var.alt*n,0,1),n.var.alt,n)
  	cl.proto<-numeric(n.var.alt)

  	## the corresponding global scores
  	glob.eval <- numeric(n.var.alt)
  	a <- capacity(c(0:(2^n-3),(2^n-3),(2^n-3))/(2^n-3))
  	for (i in 1:n.var.alt)
    		glob.eval[i] <- Choquet.integral(a,P[i,])

  	cl.proto[glob.eval <= 0.33] <- 1

        ## decomment here if there should be errors in the
        ## classification of the prototypees 
        # cl.proto[glob.eval > 0.33 & glob.eval<=0.44] <-2
  	# cl.proto[glob.eval > 0.44 & glob.eval<=0.55] <-1
  	# cl.proto[glob.eval > 0.55 & glob.eval<=0.66] <-2
  	
	cl.proto[glob.eval>0.33 & glob.eval<=0.66] <-2
 
	cl.proto[glob.eval > 0.66] <- 3

        ## a Shapley preorder constraint matrix
        ## Sh(1) > Sh(2)
        ## Sh(3) > Sh(4)
        delta.S <-0.01
	Asp <- rbind(c(1,2,delta.S), c(3,4,delta.S))
	# Asp <- NULL
        
        ## a Shapley interval constraint matrix
        ## 0.3 <= Sh(1) <= 0.9 
        # Asi <- rbind(c(1,0.1,0.2))
         Asi <- NULL
        
        ## an interaction preorder constraint matrix
        ## such that I(12) > I(34)
        delta.I <- 0.01
        Aip <- rbind(c(1,2,3,4,delta.I))
        # Aip <- NULL
        
        ## an interaction interval constraint matrix
        ## i.e. 0.2 <= I(12) <= 0.4 
        ## delta.I <- 0.01
        # Aii <- rbind(c(1,2,0.2,0.4))
         Aii <- NULL
        
        ## an inter-additive partition constraint
        ## criteria 1,2 and criteria 3,4 are indepedent 
        # Aiap <- c(1,1,2,2)
         Aiap <- NULL
        
	print("*** Starting the calculations")
	## search for a capacity which satisfies the constraints
	lsc <- ls.sorting.capa.ident(n ,k, P, cl.proto, d,
                                   A.Shapley.preorder = Asp,
                                   A.Shapley.interval = Asi,
                                   A.interaction.preorder = Aip,
                                   A.interaction.interval = Aii,
                                   A.inter.additive.partition = Aiap)

        print("")
        print("Output of the QP")
        print(lsc$how)

	## analyse the quality of the model (classify the prototypes by the model)
        print("*** Starting the analysis of the results")

	lst <- ls.sorting.treatment(P,cl.proto,lsc$solution,P,cl.proto)

        barplot(Shapley.value(lsc$solution), main = "Learnt a", sub="Shapley")

        print("Assignments of the prototypes")
        print(lst$class.A)
        print("Assignment types")
        print(lst$correct.A)
        print("Evaluation")
        print(lst$eval.correct)

	stopifnot(lst$eval.correct[1] == 1)
	
	## Generate a second set of random alternatives (A)
	## Their "correct" class is determined as beforehand with the
	## randomly generated capacity
	## The goal is to see if we can reproduce this classification
	## by the capacity learnt from the prototypes

	## a randomly generated criteria matrix of alternatives
        print("*** Generating new alternatives")
	A <- matrix(runif(n.var.alt.A*n,0,1),n.var.alt.A,n)
        cl.A <-numeric(n.var.alt.A)
	
	## the corresponding global scores
        gA <- numeric(n.var.alt.A)
        for (i in 1:n.var.alt.A)
                gA[i] <- Choquet.integral(a,A[i,])

	cl.A[gA <= 0.33] <- 1
	
	## decomment this if there should be some errors in the classes of A
        # cl.A[gA > 0.33 & gA<=0.44] <-2
        # cl.A[gA > 0.44 & gA<=0.55] <-1
        # cl.A[gA > 0.55 & gA<=0.66] <-2

        cl.A[gA>0.33 & gA<=0.66] <-2

        cl.A[gA > 0.66] <- 3

	## let us now classify the alternatives of A according to the model
	## built on P
	print("*** Applying the model on these new alternatives")
	lst <- ls.sorting.treatment(P,cl.proto,lsc$solution,A,cl.A)


        ## decomment this in case the original class is not given
                
        ## Results
        print("Assignment of the alternatives of A")
        print(lst$class.A)
        print("Type of assignments")
        print(lst$correct.A)
        print("Evaluation")
        print(lst$eval.correct)

	
	
#       ## Show the original capacity
# 	x11()
# 	barplot(Shapley.value(a), main="Original a", sub="Shapley")
#         
#         print("Summary of the learnt capacity")
#         print(summary(lsc$solution))
# 
#         ## visualisation of the classes
#         choq.proto <- numeric(length(cl.proto))
#         for (i in 1:length(cl.proto))
#           choq.proto[i] <- Choquet.integral(lsc$solution,P[i,])
#         
#         max.cl.proto <- max(cl.proto)
#         min.cl.proto <- min(cl.proto)
#         pos <- 0.9*(max.cl.proto+min.cl.proto)/2
#         
#         x11()
#         plot(choq.proto, cl.proto, col="red")
#         q<-xy.coords(lst$choq.A, rep(pos,length(lst$choq.A)))
#         points(q,col="blue")
#         
#         print("Min and max of the prototype classes")
#         print(lst$minmax)
