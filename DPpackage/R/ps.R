### ps.R                   
###
### function for use in PSgam formulae to specify smooth term, e.g. 
### ps(x0,x1,x2,k=40,degree=3,prod=1) specifies a penalized spline 
### regression of x0,x1 and x2, using B-splines of degree 3, with
### 40 equally spaced knots and a first order difference penalty 
### when it enters the model.
###
### Copyright: Alejandro Jara, 2007-2012.
###
### Last modification: 23-07-2007.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The author's contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
###      Fax  : +56-2-3547729  Email: atjara@uc.cl
###

ps <- function (..., k=50,degree=3,pord=1)
{
   bs <- "bs"
   vars <- as.list(substitute(list(...)))[-1]
   d <- length(vars)
   term <- deparse(vars[[1]],backtick=TRUE,width.cutoff=500) # first covariate
   if (term[1]==".") stop("ps(.) not yet supported.")
   if (d>1) # then deal with further covariates
   for (i in 2:d)
   { term[i]<-deparse(vars[[i]],backtick=TRUE,width.cutoff=500)
     if (term[i]==".") stop("ps(.) not yet supported.")
   }
   for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),"term.labels")
   
   # check for repeated variables in function argument list
   if (length(unique(term))!=d) stop("Repeated variables as arguments of a smooth are not permitted")
   # assemble version of call with all options expanded as text
   full.call<-paste("ps(",term[1],sep="")
   if (d>1) for (i in 2:d) full.call <- paste(full.call,",",term[i],sep="")
   label <- paste(full.call,")",sep="") # used for labelling parameters
   full.call<-paste(full.call,",k=",k,",degree=",degree,
                    ",pord=",pord,")",sep="")
   ret<-list(term=term,bs=bs,bs.dim=k,dim=d,s.degree=degree,pord=pord,
             full.call=full.call,label=label)
   class(ret)<-paste(bs,".smooth.spec",sep="")
   ret
}
