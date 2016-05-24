cpa <- function(){
    cpa.function.env <- new.env()
    cpa.results <- new.env()

    X=V=A=NULL
    rm(X, V, A)
    
#    Read and store data

OpenX <- function()  {
	input <- tclvalue(tkgetOpenFile())
	if(input == ""){
            stop()
            }else{
                X <- read.table(input, header=TRUE, sep=',', dec='.')
                assign('X', X, cpa.function.env)
            }
}

#    Read and store the variable list

OpenV <- function()  {
	input <- tclvalue(tkgetOpenFile())
	if(input == ""){
            stop()
        }else{
            V <- read.table(input, header=FALSE, sep=',')
            assign('V', V, cpa.function.env)
        }
}

#    Read and store the direct interaction list

OpenA <- function()  {
	input <- tclvalue(tkgetOpenFile())
	if(input == ""){
            stop()
        }else{
            A <- read.table(input, header=FALSE, sep=',')
            assign('A', A, cpa.function.env)
        }
}

#    Close the program

Close <- function(){
	tkdestroy(tt)
}

#    Draw the model

PlotDAG <- function(){
    X <- cpa.function.env$X
    V <- cpa.function.env$V
    A <- cpa.function.env$A
    if(is.null(V)==TRUE | is.null(A)==TRUE | ncol(A)!=2 | ncol(V)!=1){
        stop()
    }else{
	theta <- seq(0,360)*(pi/180)
	x <- cos(theta)
	y <- sin(theta)
	p <- cbind(x, y)
	V <- as.vector(V[,1])
	n <- floor(nrow(p)/length(V))
	mat <- matrix(ncol=2, nrow=0)
	for(i in 0:length(V)){
		ind <- n*i
		mat1 <- p[ind,]
		mat <- rbind(mat, mat1)
	}
	mat <- data.frame(mat, V, check.rows=FALSE, row.names=NULL)
	s <- matrix(ncol=2, nrow=0)
	for(k in 1:nrow(A)){
		Vstart <- as.character(A[k,1])
		s1 <- subset(mat, V==Vstart, select=-c(V))
		s <- rbind(s, s1)
	}
	e <- matrix(ncol=2, nrow=0)
	for(k in 1:nrow(A)){
		Vend <- as.character(A[k,2])
		e1 <- subset(mat, V==Vend, select=-c(V))
		e <- rbind(e, e1)
	}
	par(mar=c(0,0,0,0))
	plot(mat[,1],mat[,2], cex=0, xaxt='n', yaxt='n', xlim=c(min(x)-.2,max(x)+.2), ylim=c(min(y)-.2, max(y)+.2), xlab='', ylab='', frame=FALSE)
	arrows(s[,1], s[,2], e[,1], e[,2], length=.15, col=palette())
	xt <- 1.1*cos(theta)
	yt <- 1.1*sin(theta)
	pt <- cbind(xt, yt)
	mt <- matrix(ncol=2, nrow=0)
	for(i in 0:length(V)){
		ind <- n*i
		mt1 <- pt[ind,]
		mt <- rbind(mt, mt1)
	}
	text(mt[,1], mt[,2], V)
    }
}

#    Main analyses

OnOK <- function(){
    X <- cpa.function.env$X
    V <- cpa.function.env$V
    A <- cpa.function.env$A
    if(is.null(X)==TRUE | is.null(V)==TRUE | is.null(A)==TRUE | ncol(A)!=2 | ncol(V)!=1){
        stop()
    }else{
	V <- as.vector(t(V))
	A <- as.matrix(A)

        #    (1) Find all the variable pairs that do not share an interaction and store them in a dataframe
        
	VV <- V
	AA <- 0
	BB <- 0
	bset <- data.frame(nrow=0, ncol=0, check.rows=FALSE, row.names=NULL)
	for(i in 1:length(V)){
		for(j in 1:length(VV)){
			coord <- c(V[i], VV[j])
			names(coord) <- c('V1','V2')
			for(k in 1:nrow(A)){
				if(sum(coord == as.vector(A[k,]))==2){
					AA <- 1
					break
				}
			}
			if(nrow(bset)>=1){
				for(y in 1:nrow(bset)){
					if(sum(coord == c(bset[y,2], bset[y,1]))==2){
						BB <- 1
						break
					}
				}
				for(h in 1:nrow(A)){
					if(sum(coord == c(A[h,2], A[h,1]))==2){
						BB <- 1
						break
					}
				}
			}
			if(coord[1]!=coord[2] & AA != 1 & BB != 1){
				bset <- rbind(bset, coord)
			}
			AA <- 0
			BB <- 0
		}
	}
	bset <- bset[-1,]

        #    (2) For each of the variable pairs found in step (1), find the direct causes ('parents') of both effects ('children') and store them in a list

	parent <- c()
	parents <- list()
	CC <- 0
	for(i in 1:nrow(bset)){
		parent <- c(parent, i)
		for(j in 1:nrow(A)){
			if(bset[i,1] == A[j,2] | bset[i,2] == A[j,2]){
				for(k in 1:length(parent)){
					if(A[j,1] == parent[k]){
						CC <- 1
						break
					}
				}
				if(CC != 1){
					parent <- c(parent, A[j,1])        
				}
				CC <- 0
			}
		}
		parents[[length(parents)+1]] <- parent
		parent <- c()
	}

        #    (3) Build a set of conditioning variables (using the list from step (3)) for the conditional independence tests

	conditions <- c()
	for(i in 1:length(parents)){
		a <- ""
		if(length(parents[[i]]) > 1){
			for(j in 2:length(parents[[i]])){
				a <- paste(a, parents[[i]][j], sep=' + ')
			}
			conditions <- c(conditions, a)
		}else{
			conditions <- c(conditions, "")
		}
	}

        #    (4) Test (with multiple linear models) the conditional independence claims and store the models and their results in lists

	tests <- list()
	tests_summary <- list()
	for(i in 1:nrow(bset)){
		a <- paste(bset[i,2], '~', bset[i,1], conditions[i], sep=' ')
		a <- as.formula(noquote(a))
		test <- lm(a, data=X)
		testsummary <- summary(test)
		tests[[length(tests)+1]] <- test
		tests_summary[[length(tests_summary)+1]] <- testsummary
	}

	p <- c()
	for(i in 1:length(tests_summary)){
		a <- coef(tests_summary[[i]])[,4][bset[i,1]]
		p <- c(p, a)
	}

        Ti <- tests
        Ti_summary <- tests_summary

        #    (5) Calculate relevant statistics (Fisher's C and P value)

	Cf <- -2*sum(log(p))

	P <- pchisq(Cf, 2*length(p), lower.tail=FALSE)

        #    (6) Build a fancy table with the d-separation statements and the null probabilities associated with

	Bu <- c()
	for(i in 1:length(parents)){
		a <- ""
		b <- ""
		if(length(parents[[i]]) > 1){
			for(j in 2:length(parents[[i]])){
				a <- paste(a, parents[[i]][j], sep=',')
			}
			a <- substr(a, 2, nchar(a))
			b <- paste('(', bset[i,1], ',', bset[i,2], ')', '|', '{', a, '}', sep='')
			Bu <- c(Bu, b)
		}else{
			Bu <- c(Bu, "{}")
		}
	}

	Bu <- as.data.frame(cbind(Bu, signif(p, 4)), check.rows=FALSE, row.names=NULL)
	names(Bu) <- c("Independence claim", "P-value")
	row.names(Bu) <- NULL

        #    (7) Find the direct causes of each variable in the dataset and store them in a list

        dclist <- list()
        dc <- c()
        for(i in 1:length(V)){
            I <- V[i]
            IA <- subset(A, A[,2]==I)
            if(nrow(IA) >= 1){
                for(k in 1:nrow(IA)){
                    dc <- c(dc, IA[k,1])
                }
            }else{
                dc <- NULL
            }
            dclist[[i]] <- dc
            dc <- c()
        }
        names(dclist) <- V
        dclist <- dclist

        #    (8) Calculate the AIC (the number of free parameters for each regression model is calculated considering the intercept and the residual error plus one additional parameter for each independent variable)

        Kt <- c()
        for(i in 1:length(dclist)){
            Ki <- 2 + length(dclist[[i]])
            Kt <- c(Kt, Ki)
        }

        K <- sum(Kt)

        AIC <- Cf + 2*K*(nrow(X)/(nrow(X)-K-1))

        #    (9) Build a set of predictors to test the strenght of the direct interactions by (multiple) linear regression

	dcconditions <- c()
	for(i in 1:length(dclist)){
		a <- ""
		if(is.null(dclist[[i]]) | length(dclist[[i]]) <= 1){
                    dcconditions <- c(dcconditions, "")
                }else{
                    for(j in 2:length(dclist[[i]])){
                        a <- paste(a, dclist[[i]][j], sep=' + ')
                    }
                    dcconditions <- c(dcconditions, a)
                }
	}
        
        #    (10) Test (with multiple linear models) the direct interactions

	dctests <- list()
	dctests_summary <- list()
        dctnames <- c()
	for(i in 1:length(V)){
            if(is.null(dclist[[i]])){
            }else{
                a <- paste(V[i], '~', dclist[[i]][1], dcconditions[i], sep=' ')
		a <- as.formula(noquote(a))
		dctest <- lm(a, data=X)
		dctestsummary <- summary(dctest)
		dctests[[length(dctests)+1]] <- dctest
		dctests_summary[[length(dctests_summary)+1]] <- dctestsummary
                dctnames <- c(dctnames, V[i])
            }
	}

        names(dctests) <- dctnames
        names(dctests_summary) <- dctnames

        #    (11) Assign the relevant objects to a new environment

        cpa.env <- new.env()

        assign('C', Cf, cpa.env)
        assign('P', P, cpa.env)
        assign('Bu', Bu, cpa.env)
        assign('Ti', Ti, cpa.env)
        assign('Ti_summary', Ti_summary, cpa.env)
        assign('AIC', AIC, cpa.env)
        assign('dctests', dctests, cpa.env)
        assign('dctests_summary', dctests_summary, cpa.env)

        #    (12) Print messages and save the results
                
	if(P <= 0.05){
		msg <- paste("\nC = ", round(Cf, 6), "\n\nP = ", round(P, 6), "\n\nAIC = ", round(AIC, 6), '', "\n\n\nThe hypothesized causal model is incompatible with the covariance structure of the data\n\nWould you like to save the results?\n\n")
	}else{
		msg <- paste("\nC = ", round(Cf, 6), "\n\nP = ", round(P, 6), "\n\nAIC = ", round(AIC, 6), '', "\n\n\nThe hypothesized causal model is compatible with the covariance structure of the data\n\nWould you like to save the results?\n\n")
            }
        mgsval <- tkmessageBox(message=msg, type='yesno', default='yes')
        if(tclvalue(mgsval)=='yes'){
            ttsave<-tktoplevel()
            Name <- tclVar("")
            entry.Name <-tkentry(ttsave,width="20",textvariable=Name)
            tkgrid(tklabel(ttsave,text="Enter the name of the environment that will contain the results"))
            tkgrid(entry.Name)
            Save <- function(){
                NameVal <- tclvalue(Name)
                tkdestroy(ttsave)
                assign(NameVal, cpa.env, cpa.results)
            }
            Save.but <-tkbutton(ttsave,text="   Save   ",command=Save)
            tkbind(entry.Name, "<Return>",Save)
            tkgrid(Save.but)
            tkfocus(ttsave)
        }else{
            return()
        }
    }
}

# Build the GUI

tt<-tktoplevel()
tkwm.title(tt, "Data")
entry.X <- tkbutton(tt, text = "  Browse...  ", command=OpenX)
entry.V <- tkbutton(tt, text = "  Browse...  ", command=OpenV)
entry.A <- tkbutton(tt, text = "  Browse...  ", command=OpenA)
entry.C <- tkbutton(tt, text = "  Exit  ", command=Close)
entry.DAG <- tkbutton(tt, text = "  Plot  ", command=PlotDAG)
tkgrid(tklabel(tt, text = "    "))
tkgrid(tklabel(tt, text = "  Select the file with the data matrix:"), entry.X, sticky='w')
tkgrid(tklabel(tt, text = "  Select the file with the variables list:"), entry.V, sticky='w')
tkgrid(tklabel(tt, text = "  Select the file with the direct interactions list:"), entry.A, sticky='w')
tkgrid(tklabel(tt, text = "    "))
tkgrid(tklabel(tt, text = "  Generate a graph (note: both the variable and the direct interaction lists must be already selected)"), entry.DAG, sticky='w')
tkgrid(tklabel(tt, text = "    "))
tkgrid(tklabel(tt, text = "    "))
tkgrid(tklabel(tt, text = "							"), entry.C, sticky='w')
tkgrid(tklabel(tt, text = "    "))
OK.but <-tkbutton(tt, text="  Run  ", command=OnOK)
tkgrid(  OK.but)
tkfocus(tt)
cpa.results
}
