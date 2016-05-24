# Nom fichier :        utilitaires.R
# Projet :             stacomiR


#############################################
# functions copied from Hmisc
#############################################
#' function used to print the html tables of output (see xtable documentation)
#' @param data a data frame
#' @param caption the caption
#' @param top of top the caption is placed on top
#' @param outfile outfile is the path to the file
#' @param clipboard if clipboard TRUE, a copy to the clipboard is made
#' @param append is the file appended to the previous one ?
#' @param digits 
#' @param ... 
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
funhtml=function(data,caption=NULL,top=TRUE,outfile=NULL,clipboard=FALSE,append=TRUE,digits=NULL,...){
	data[is.na(data)]<-""
	xt=xtable(data, caption=caption,digits=digits)
	xt=print(xt,type="html",caption.placement="top",file=outfile)
	# pour changer le defaut "bottom" des caption
	if (clipboard) writeClipboard(xt) 
} 
###########################################
# special functions (exported as they are usefull
#############################################
#' This function replaces the variable names in a data.frame
#' @param objet a data frame
#' @param old_variable_name 
#' @param new_variable_name 
#' @returnType data.frame
#' @return objet
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
chnames=function(objet,
		old_variable_name,
		new_variable_name){
		if (length(old_variable_name)!=length(new_variable_name)) stop("les variables de remplacement doivent avoir le meme nombre que les variables de depart")
		if (!all(!is.na(match(old_variable_name,colnames(objet))))) {
		   stop(paste("les noms",paste(is.na(match(old_variable_name,colnames(objet))),collapse="/"),"ne correspondent pas aux variables du tableau"))
    }
	colnames(objet)[match(old_variable_name,colnames(objet))]<- new_variable_name
	return(objet)
}

# fonction qui retourne l'index des valeurs repetees d'un vecteur
# see duplicated

#' fonction qui renvoit l'index des valeurs apparaissant une seule fois
#' @param a 
#' @returnType vector
#' @return the index unique  values within a vector
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
induk=function(a){
	sol=match(unique(a),a)     #index des valeurs uniques
	return(sol)   
}


#' very usefull function used to "kill" these bloody factors that appears, noticeably after loading with odbc
#' @param df a data.frame
#' @returnType data.frame
#' @return df
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
killfactor=function(df){
	for (i in 1:ncol(df))
	{
		if(is.factor(df[,i])) df[,i]=as.character(df[,i])
	}
	return(df)
}

#' ex fonction to write to excel, not used within the program but can still be used
#' @param d 
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
ex<-function(d=NULL){
	if (is.null(d)){
		xl=select.list(choices=ls(envir=globalenv()), preselect = NULL, multiple = FALSE, title = "choisir l'objet")
		write.table(get(xl),"clipboard",sep="\t",col.names=NA)
	} else {
		write.table(d,"clipboard",sep="\t",col.names=NA)
	}
}



#' id.odd function modified from package sma (which did not verify that the entry was indeed an integer)
#' @param x 
#' @returnType logical
#' @return a logical
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
is.odd=function (x) 
{
    if (x==as.integer(x)) {
        if (x%%2 == 0) {
            return(FALSE)
        }
        else {
            return(TRUE)
        }
    }
    else {
        stop("is.odd should be used with an integer")
    }
}
#' is.even function modified from package sma (which did not verified that the entry was indeed an integer)
#' @param x 
#' @returnType logical
#' @return a logical
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
is.even=function (x) 
{
    if (x==as.integer(x)) {
        if (x%%2 != 0) {
            return(FALSE)
        }
        else {
            return(TRUE)
        }
    }
    else {
        stop("is.even should be used with an integer")
    }
}

#' Function to transform a ftable into dataframe but just keeping the counts works with ftable of dim 2
#' @param tab 
#' @author Cedric Briand \email{cedric.briand00@@gmail.com}
#' @export
tab2df=function(tab){
	if (length((attributes(tab)$dim))>2) stop("only works with tables of dim 2")
	df=as.data.frame(matrix(as.vector(tab),nrow(tab),ncol(tab)))
	rownames(df)<-attributes(tab)$row.vars[[1]]
	colnames(df)<-attributes(tab)$col.vars[[1]]	
	return(df)
}

