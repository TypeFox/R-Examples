`esttab` <-
function(t.value=FALSE,p.value=FALSE,round.dec=3,caption=NULL,label=NULL,texfontsize=NULL,sig.levels=c(0.1,0.05,0.01),sig.sym=c("*","**","***"),filename=NULL,csv=FALSE,dcolumn=NULL,table="table",table.pos="htbp",caption.top=FALSE,booktabs=FALSE,var.order=NULL,sub.sections=NULL,var.rename=NULL,resizebox=c(0,0),colnumber=FALSE,store="default"){

# reading list from eststo
prev.list <- paste(store,".ccl",sep="")
model.coeff.list <- eval(lapply(prev.list,as.name)[[1]],estout:::estoutstorage)

# catching use if empty ccl
if(is.list(model.coeff.list)){}else{return("No values stored. I think you need to store some models first.")}

# setting dcolumn = NULL
if(is.null(dcolumn)){dcolumn <- "c"}else{dcolumn <- dcolumn}

# setting tablepos if non-NULL
if(is.null(table.pos)){}else{table.pos <- paste("[",table.pos,"]",sep="")}

# setting caption for TeX
if(is.null(caption)){texcaption <- caption}else{texcaption <- paste("\\caption{",caption,"}\n",sep="")}

# setting label if non-NULL
if(is.null(label)){}else{label <- paste("\\label{tab:",label,"}\n",sep="")}

# setting texfontsize if non-NULL
if(is.null(texfontsize)){}else{texfont <- paste(texfontsize,"\n",sep="")}
# converting sub.sections vector to list of type list[i] = c(position,text)
sub.sections.list <- NULL
if(is.null(sub.sections) == FALSE){
	if((length(sub.sections)%%2) == 0){
		sub.sections.list  <- list() 
		for(i in 1:(length(sub.sections)/2)){
			sub.sections.list[[i]]  <- c(sub.sections[[2*i-1]], sub.sections[[2*i]])
		}
	}
	else{
		return(cat("Your 'sub.section' vector contains a mistake. Please check and run again."))
	}
}
#print(sub.sections.list) #---- control

# converting var.names vector to list of type list[i] = c(name.old,name.new)
var.rename.list <- NULL
if(is.null(var.rename) == FALSE){
	if((length(var.rename)%%2) == 0){
		var.rename.list  <- list() 
		for(i in 1:(length(var.rename)/2)){
			var.rename.list[[i]]  <- c(var.rename[[2*i-1]], var.rename[[2*i]])
		}
	}
	else{
		return(cat("Your 'var.rename' vector contains a mistake. Please check and run again."))
	}
}
#print(var.rename.list)   #control --------!!!!!!!!!!!
# setting TeX commands for booktabs
if(booktabs == TRUE){
	toprule <- "\\toprule\n"
	midrule <- "\\midrule\n"
	bottomrule <- "\\bottomrule\n"
}
else{
	toprule <- "\\hline\\hline\n"
	midrule <- "\\hline\n"
	bottomrule <- "\\hline\\hline\n"
}

# starting list of used variables
var_list <- c()#model.coeff.list[[1]][[1]][[1]]

# for CSV
if(csv == TRUE){
	save.file <- paste(filename,".csv",sep="")
	delimiter <- ","
	R2 <- "R^2"
	aR2 <- "adj.R^2"
	N <- "N"
	om.end <- "\n"
	caption <- paste(caption,om.end,sep="")
	threestar <- paste(sig.sym[[3]],sep="")
	twostar <- paste(sig.sym[[2]],sep="")
	onestar <- paste(sig.sym[[1]],sep="")
}
# for TeX
else{
	save.file <- paste(filename,".tex",sep="")
	delimiter <- "\t\t&"
	R2 <- "$R^2$"
	aR2 <- "$adj.R^2$"
	N <- "$N$"
	om.end <- "\\\\\n"
	caption <- caption
	threestar <- paste("\\sym{",sig.sym[[3]],"}",sep="")
	twostar <- paste("\\sym{",sig.sym[[2]],"}",sep="")
	onestar <- paste("\\sym{",sig.sym[[1]],"}",sep="")
}
######### creating a vector of variable names ##########################
if(is.null(var.order) == FALSE){
	var_list <- var.order
}
else{
	for (j in 1:length(model.coeff.list)){
	        col_length <- length(model.coeff.list[[j]]) # length of one model column list
	        for (i in 1:(col_length-1)){
	                var.count <- 0
			for(var.exist in grepl(model.coeff.list[[j]][[i]][[1]],var_list)){
				if(var.exist == FALSE){
					var.count <- var.count + 1
					
				}
				else if(var.exist == TRUE){
					break
				}
			}
			if(var.count == length(var_list)){
	                          var_list <- c(var_list,model.coeff.list[[j]][[i]][[1]])
			}
		}
	}
}
var.count <- length(var_list)
#cat(var.count)
#cat("\n") #control
#cat(var_list)
#cat("\n")  #control
# creating code for resizebox in TeX
resizebox <- resizebox
if(resizebox[1] == 0 && resizebox[2] == 0){
    resize=FALSE
}
else{
    resizetop <- paste("\n\\resizebox{",resizebox[1],"}{",resizebox[2],"}{",sep="")
    resize <- TRUE
}
#########################################################################

######################## making a matrix  to be filled###################
index2 <- length(model.coeff.list[[1]])  # index length of columns in model 1
index3 <- length(model.coeff.list[[1]][[index2]]) # index length of values in last columns

if(model.coeff.list[[1]][[index2]][[index3]] == "lm"){     # checking if model is of lm class
adds <- 4
}
if(model.coeff.list[[1]][[index2]][[index3]] == "plm"){     # checking if model is of plm class
adds <- 4
}
if(model.coeff.list[[1]][[index2]][[index3]] == "glm"){     # checking if model is of glm class (here Zelig package)
adds <- 4
}

om.ncol = length(model.coeff.list) + 2	# om.ncol = number of columns of output matrix
om.nrow = var.count*2 + adds		# om.nrow = number of rows of output matrix
output_matrix <- matrix(delimiter,om.nrow,om.ncol)		# make matrix[delimiter,om.nrow,om.ncol]
output_matrix[,om.ncol] <- om.end							# fill last column
for(i in 1:var.count){
        output_matrix[i*2,1] <- var_list[i]  #insert var-name in 1st column
        output_matrix[i*2+1,1] <- " "           #clean & from first column
}
#-
if(model.coeff.list[[1]][[index2]][[index3]] == "lm"){     # checking if model is of lm class
output_matrix[length(output_matrix[,1])-2,1] <- R2
output_matrix[length(output_matrix[,1])-1,1] <- aR2
output_matrix[length(output_matrix[,1]),1] <- N
end.sep.line <- om.nrow - 3
}

if(model.coeff.list[[1]][[index2]][[index3]] == "plm"){     # checking if model is of plm class
output_matrix[length(output_matrix[,1])-2,1] <- R2
output_matrix[length(output_matrix[,1])-1,1] <- aR2
output_matrix[length(output_matrix[,1]),1] <- N
end.sep.line <- om.nrow - 3
}
if(model.coeff.list[[1]][[index2]][[index3]] == "glm"){     # checking if model is of glm class
output_matrix[length(output_matrix[,1]),1] <- N
end.sep.line <- om.nrow - 3
}
#print(output_matrix)  #control
#########################################################################
########## making stars and putting in matrix #### stars depend on p-values ##########

for (j in 1:length(model.coeff.list)){ 
        col_length <- length(model.coeff.list[[j]])
                for (i in 1:(col_length-1)){
                        if ( model.coeff.list[[j]][[i]][[5]] < sig.levels[3] ) {
                                sigs <- paste(delimiter,round(model.coeff.list[[j]][[i]][[2]],round.dec),threestar,sep="") #coefficient
                                std_err <- paste(delimiter,"(",round(model.coeff.list[[j]][[i]][[3]],round.dec),")",sep="")        #std.err
                                t_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[4]],round.dec),"]",sep="")          #t-value
                                p_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[5]],round.dec),"]",sep="")          #p-value
                        }
                        else if( model.coeff.list[[j]][[i]][[5]] < sig.levels[2] ) {
                                sigs <- paste(delimiter,round(model.coeff.list[[j]][[i]][[2]],round.dec),twostar,sep="")
                                std_err <- paste(delimiter,"(",round(model.coeff.list[[j]][[i]][[3]],round.dec),")",sep="")
                                t_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[4]],round.dec),"]",sep="")
                                p_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[5]],round.dec),"]",sep="")          #p-value
                        }
                        else if( model.coeff.list[[j]][[i]][[5]] < sig.levels[1] ) {
                                sigs <- paste(delimiter,round(model.coeff.list[[j]][[i]][[2]],round.dec),onestar,sep="")
                                std_err <- paste(delimiter,"(",round(model.coeff.list[[j]][[i]][[3]],round.dec),")",sep="")
                                t_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[4]],round.dec),"]",sep="")
                                p_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[5]],round.dec),"]",sep="")          #p-value
                        }
                        else{
                                sigs <- paste(delimiter,round(model.coeff.list[[j]][[i]][[2]],round.dec),"",sep="")
                                std_err <- paste(delimiter,"(",round(model.coeff.list[[j]][[i]][[3]],round.dec),")",sep="")
                                t_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[4]],round.dec),"]",sep="")
                                p_val <- paste(delimiter,"[",round(model.coeff.list[[j]][[i]][[5]],round.dec),"]",sep="")          #p-value
                        }
                        k <- 1
                        while(var_list[k] != model.coeff.list[[j]][[i]][[1]]){
                                k <- k +1
                       }
#print(k) #control-----------!!!!!!!!!!!!!
                        output_matrix[k*2,j+1] <- sigs #entry of coefficients
                        if(t.value == TRUE){
                                output_matrix[k*2+1,j+1] <- t_val  #if set entry of t-values
                        }
                        else if(p.value == TRUE){
                                output_matrix[k*2+1,j+1] <- p_val  #if set entry of t-values
                        }
                        else{
                                output_matrix[k*2+1,j+1] <- std_err  #if set entry of std.err.
                        }
                } # end for (i)
#print(j) #control-----------------!!!!!!!!!!!
#print(model.coeff.list[[j]][[col_length]][[1]]) #control-------------!!!!

		# rename dep.var
		dep.var <- deparse(model.coeff.list[[j]][[col_length]][[1]])
		if(is.null(var.rename) == FALSE){
			for(i in 1:length(var.rename.list)){
				if(dep.var == var.rename.list[[i]][[1]]){
					dep.var <- var.rename.list[[i]][[2]]
				}
			}
		}
#print(is.null(var.rename)) #control ---------------!!!!!!!!!

		#--------------------------------------------

                output_matrix[1,1] <- " "
                if(csv == TRUE){
                	output_matrix[1,j+1] <- paste(delimiter,dep.var,sep="") #dep.var
                }
                else{
	            	output_matrix[1,j+1] <- paste("&\\multicolumn{1}{c}{",dep.var,"}",sep="") #dep.var
                }
                if(model.coeff.list[[j]][[col_length]][[length(model.coeff.list[[j]][[col_length]])]] == "lm"){     # checking if model is of lm class
#print("Model is of class 'lm'")
    	           	output_matrix[length(output_matrix[,j+1])-2,j+1] <- paste(delimiter,round(model.coeff.list[[j]][[col_length]][[2]],round.dec),sep="") # Rsquared
        	    	output_matrix[length(output_matrix[,j+1])-1,j+1] <- paste(delimiter,round(model.coeff.list[[j]][[col_length]][[3]],round.dec),sep="") # adj.Rsquared
                	output_matrix[length(output_matrix[,j+1]),j+1] <- paste(delimiter,if(! csv==TRUE){"\\multicolumn{1}{c}{"},model.coeff.list[[j]][[col_length]][[4]],if(! csv==TRUE){"}"},sep="") # N
                }
                if(model.coeff.list[[j]][[col_length]][[length(model.coeff.list[[j]][[col_length]])]] == "plm"){     # checking if model is of plm class
    	           	output_matrix[length(output_matrix[,j+1])-2,j+1] <- paste(delimiter,round(model.coeff.list[[j]][[col_length]][[2]],round.dec),sep="") # Rsquared
        	    	output_matrix[length(output_matrix[,j+1])-1,j+1] <- paste(delimiter,round(model.coeff.list[[j]][[col_length]][[3]],round.dec),sep="") # adj.Rsquared
                	output_matrix[length(output_matrix[,j+1]),j+1] <- paste(delimiter,if(! csv==TRUE){"\\multicolumn{1}{c}{"},model.coeff.list[[j]][[col_length]][[4]],if(! csv==TRUE){"}"},sep="") # N
#print("Model is of class 'plm'")
                }
                if(model.coeff.list[[j]][[col_length]][[length(model.coeff.list[[j]][[col_length]])]] == "glm"){     # checking if model is of glm class
#print("Model is of class 'glm'")
                	output_matrix[length(output_matrix[,j+1]),j+1] <- paste(delimiter,if(! csv==TRUE){"\\multicolumn{1}{c}{"},model.coeff.list[[j]][[col_length]][[2]],if(! csv==TRUE){"}"},sep="") # N
                }
		#cat(output_matrix)  #control ----------!!!!!!!!!!
} # end for model (j)
#----------------------------------------------

# replacing regression var names with real var names
if(is.null(var.rename) == FALSE){
	for(j in 2:length(output_matrix[,1])){
		if(is.null(var.rename) == FALSE){
			for(i in 1:length(var.rename.list)){
				if(output_matrix[j,1] == var.rename.list[[i]][[1]]){
					output_matrix[j,1] <- var.rename.list[[i]][[2]]
				}
			}
		}
	}
}
#----------------------------------------------

if(csv == TRUE ){
	# begin sink to file if filename!=NULL
	if(! is.null(filename)){
		sink(save.file)
	}
	cat(caption)
	for(i in 1:length(output_matrix[,1])){
		cat(output_matrix[i,])
	}
	if(t.value==TRUE){
	        cat("t-values in brackets")
	}
	else if((p.value==TRUE) && (t.value==FALSE)){
	        cat("p-values in brackets")
	}
	else{
	        cat("Standard errors in parentheses\n")
	}
	# end sink
	if(! is.null(filename)){
		sink()
	}
}
### writing tex
else{
# begin sink
	if(! is.null(filename)){
		sink(save.file)
	}
# collate TeX formatted table
	cat("\\def\\sym#1{\\ifmmode^{#1}\\else\\(^{#1}\\)\\fi}\n\\begin{",table,"}",table.pos, if(caption.top==TRUE){texcaption},"
\\centering\n",texfontsize,
if(resize==TRUE){resizetop},"
\\begin{tabular}{l*{",om.ncol-2,"}{",dcolumn,"}}
",toprule,sep="")
	for (j in 1:om.nrow){
	                if(j == 1){
                            if(colnumber==TRUE){
                                cat(paste(paste("\t","&\\multicolumn{1}{c}{(",1:(om.ncol-2),")} ",collapse=" "),"\\\\",collapse=""))
                            }
	                }
	                cat(output_matrix[j,])			# writing output <- matrix
	                if(j == 1){
	                    cat(midrule)
	                }
	                if(j == (end.sep.line)){
	                    cat(midrule)
	                }
			if(!is.null(sub.sections)){
				for(k in 1:length(sub.sections.list)){
					if(j == (as.integer(sub.sections.list[[k]][[1]])*2+1)){
						cat("\\multicolumn{2}{c}{",sub.sections.list[[k]][[2]],"}",rep("&",om.ncol-3),"\\\\\n",sep="")
					}
				}
			}
	}
	cat(bottomrule)
	if(t.value==TRUE){
	        cat(paste("\\multicolumn{",om.ncol-1,"}{l}{\\footnotesize t-values in brackets}\\\\\n",sep=""))
	}
	else if((p.value==TRUE) && (t.value==FALSE)){
	        cat(paste("\\multicolumn{",om.ncol-1,"}{l}{\\footnotesize p-values in brackets}\\\\\n",sep=""))
	}
	else{
	        cat(paste("\\multicolumn{",om.ncol-1,"}{l}{\\footnotesize Standard errors in parentheses}\\\\\n",sep=""))
	}
	cat(paste("\\multicolumn{",om.ncol-1,"}{l}{\\footnotesize $^{",sig.sym[1],"}$ (p $\\le$ ",sig.levels[1],"), $^{",sig.sym[2],"}$ (p $\\le$ ",sig.levels[2],"), $^{",sig.sym[3],"}$ (p $\\le$ ",sig.levels[3],")}\\\\
\\end{tabular}\n",if(resize==TRUE){"}\n"},
if(caption.top==FALSE){texcaption},label,"\\end{",table,"}\n",sep=""))
	# end sink
	if(! is.null(filename)){
		sink()
	}
} 
# end writing tex
# esttab-end
} 
