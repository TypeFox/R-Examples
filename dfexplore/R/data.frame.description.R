# Classe "data.frame.description"

# Définition de la classe
setClass(Class="Data.frame.description",
         representation=representation( 
           df_source = "data.frame",
           df_pour_plot = "data.frame",
           title = "character"
         )
)


#### Methods ####

# Generic method
setGeneric(name="produce_df_pour_plot",
           def=function(source_data_frame){
             standardGeneric("produce_df_pour_plot")
           }
)

#methode produisant un objet description
setMethod(f="produce_df_pour_plot",
          signature="data.frame",
          def=function(source_data_frame){
            
            # Check if there is matrix inside and if true, transform to standard
            matrix_position <- search_matrix_position(source_data_frame)
            if(length(matrix_position) > 0){
              source_data_frame <- expand_dfmatrix( dataframe=source_data_frame, matrix_position)
            }
            # Function to name to an human understandable way classes of data
            class_detect <- function(vector_tested,tested_list){
              #find the internal class of the vector. Keep only the first class.
              vector_class <- class(vector_tested)[1]
              
              # TODO : if the class is a matrix, then extend it with nameofmatrix.nameofvar
              
              #select if in tested list
              if(vector_class %in% tested_list){
                return(vector_class)
              }
              else{
                return("other")
              }
            }
            
            # List the names of the variable in an ordered factor (good for plot purposes)
            variables_names <- factor(x=names(source_data_frame),levels=names(source_data_frame),ordered=T)
            
            #Count the number of observations
            nb_obs <- nrow(source_data_frame)
            
            modalities <- c("NA", "integer",   "numeric",   "Date",         "character", "factor", "ordered",  "other")
            # Create a data.frame in a long format with the presence of NA for each variable and each obs
            
            df_na <- data.frame(
              #Repeat each variable's name
              "variable"=rep(x=variables_names,each=nb_obs),
              #Find the class of each variable and repeat
              "class"=rep(x=as.factor(sapply(source_data_frame,class_detect,modalities)),each=nb_obs),
              #give a number for each observation
              "observation"= rep(x=1:nb_obs, times=length(variables_names)),
              #produce a logical vector indicating if the cell (variable+obs) is NA
              "isNA"=as.vector(is.na(source_data_frame))
            )
            
            # Calculate the modality of each cell
            cell_content <- apply(X=df_na[,c("isNA","class")],MARGIN=1,
                                  FUN=function(x){
                                    #x[1]==" TRUE" is a workaround to fix 
                                    # (because "apply" transform everything in a matrix character )
                                    if(x[1]==" TRUE"){ 
                                      # If is NA, return "NA"
                                      return(as.character("NA"))
                                    }
                                    else{
                                      # else return human redeable class
                                      return(as.character(x[2]))
                                    }
                                  }
            )
            
            #use ordered factor to have it in good order in plot
            cell_content <- factor(x=cell_content,levels=modalities,ordered=T)
            
            #produce an data.frame with only useful data
            df_for_plot<-cbind(df_na[,c("variable","observation")],cell_content)
            
            return(df_for_plot)
          }
)

# Initiateur
setMethod(f="initialize",
          signature="Data.frame.description",
          definition=function(.Object,df_source,name){
            .Object@df_source <- df_source
            .Object@df_pour_plot <- produce_df_pour_plot(df_source)
      
            
            .Object@title <- name
            return(.Object)
          }
)
#test
#new(Class="Data.frame.description",cars)

# Constructeur
data.frame.description <- function(data_frame_source,data_frame_name = as.character(as.list(match.call())$data_frame_source)){
  #Constructeur de la classe Data.frame.description 
  new(Class="Data.frame.description",data_frame_source,data_frame_name)
}

# Methode plot pour data.frame.description
