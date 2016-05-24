#' Create global variables in the PopED database
#' 
#' Function takes design variables from input files
#' and converts them to the global variables needed
#' in PopED.  Typically not used by the user.  Instead 
#' use the function \code{\link{create.poped.database}}.
#' 
#' @param poped.db A PopED database 
#' @return A PopED database
#' @family poped_input
#' @export
#' @keywords internal

## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

convert_variables <- function(poped.db){
    design = poped.db$design
    design_space = poped.db$design_space

    #poped.db$ga = zeros(poped.db$design$m,size(poped.db$design$a,2))
    #poped.db$gx = zeros(poped.db$design$m,size(poped.db$design$x,2))
    #poped.db$gxt = zeros(poped.db$design$m,max(poped.db$design_space$maxni))
    #poped.db$gminxt = zeros(poped.db$design$m,max(poped.db$design_space$maxni))
    #poped.db$gmaxxt = zeros(poped.db$design$m,max(poped.db$design_space$maxni))
    #poped.db$gmaxa = zeros(poped.db$design$m,size(poped.db$design$a,2))
    #poped.db$gmina = zeros(poped.db$design$m,size(poped.db$design$a,2))

    #     if((test_mat_size(matrix(c(poped.db$design$m, 1),nrow=1,byrow=T),design$ni,'ni')==1)){
    #         poped.db$gni = design$ni
    #     }    
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, max(poped.db$design_space$maxni)),nrow=1,byrow=T),design$model_switch,'model_switch')==1)){
    #         if(is.null(poped.db$global_model_switch)){
    #             poped.db$global_model_switch=design$model_switch
    #         } else {
    #             poped.db$global_model_switch[1:poped.db$design$m,1:max(poped.db$design_space$maxni)]=design$model_switch
    #         }
    #     }

    #     if((test_mat_size(matrix(c(poped.db$design$m, max(poped.db$design_space$maxni)),nrow=1,byrow=T),design$xt,'xt')==1)){
    #         poped.db$gxt[1:poped.db$design$m,1:max(poped.db$design_space$maxni)]=design$xt
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, max(poped.db$design_space$maxni)),nrow=1,byrow=T),design_space$maxxt,'maxxt')==1)){
    #         poped.db$gmaxxt[1:poped.db$design$m,1:max(poped.db$design_space$maxni)] = design_space$maxxt
    #     }

    #     if((test_mat_size(matrix(c(poped.db$design$m, max(poped.db$design_space$maxni)),nrow=1,byrow=T),design_space$minxt,'minxt')==1)){
    #         poped.db$gminxt[1:poped.db$design$m,1:max(poped.db$design_space$maxni)] = design_space$minxt
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, size(poped.db$design$a,2)),nrow=1,byrow=T),design$a,'a')==1)){
    #       poped.db$ga[1:poped.db$design$m,1:size(poped.db$design$a,2)]=design$a
    #     }
       
    #     if((test_mat_size(matrix(c(poped.db$design$m, size(poped.db$design$a,2)),nrow=1,byrow=T),design_space$maxa,'maxa')==1)){
    #         poped.db$gmaxa[1:poped.db$design$m,1:size(poped.db$design$a,2)]=design_space$maxa
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, size(poped.db$design$a,2)),nrow=1,byrow=T),design_space$mina,'mina')==1)){
    #         poped.db$gmina[1:poped.db$design$m,1:size(poped.db$design$a,2)]=design_space$mina
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, size(poped.db$design$x,2)),nrow=1,byrow=T),design$x,'x')==1)){
    #         poped.db$gx[1:poped.db$design$m,1:size(poped.db$design$x,2)]=design$x
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, 1),nrow=1,byrow=T),design$groupsize,'groupsize')==1)){
    #         poped.db$groupsize=design$groupsize
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, max(poped.db$design_space$maxni)),nrow=1,byrow=T),
    #                       design_space$G_xt,'G')==1)){
    #         poped.db$G=design_space$G_xt
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, size(poped.db$design$a,2)),nrow=1,byrow=T),design_space$G_a,'Ga')==1)){
    #         if(is.null(dim(design_space$G_a))){
    #             design_space$G_a=matrix(data=design_space$G_a,nrow=1,byrow=TRUE)
    #         }
    #         poped.db$Ga=design_space$G_a
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, size(poped.db$design$x,2)),nrow=1,byrow=T),design_space$G_x,'Gx')==1)){
    #         poped.db$Gx=design_space$G_x
    #     }
    
    #     if((test_mat_size(matrix(c(poped.db$design$m, 1),nrow=1,byrow=T),design_space$maxgroupsize,'maxgroupsize')==1)){
    #         poped.db$maxgroupsize=design_space$maxgroupsize
    #     }

    #     if((test_mat_size(matrix(c(poped.db$design$m, 1),nrow=1,byrow=T),design_space$mingroupsize,'mingroupsize')==1)){
    #         poped.db$mingroupsize=design_space$mingroupsize
    #     }
    
    return( poped.db) 
}
