require(caper)
require(ape)
#' A simplified interface to the main pastis function. 
#' 
#' This function assimilates sequences, taxonomic information and tree constraints 
#' into a mrBayes file. This permits the construction of trees that are 
#' compatible with all of these sources of tdata and contain all known taxa.
#' 
#' This is a simplified version of pastis_main which assumes that (i) the data are input as a pastisData object (e.g. created with read_input), or (ii) the input files 
#' all have the same base (specified by base_name) with diffferent extensions for each file
#' type:
#'  
#' .sequences: the sequence file (FASTA format)
#' 
#' .tree: the constraint tree (Newick format)
#' 
#' .taxa: a list of taxa (see below)
#' 
#' .missingclades: locations of missing clades (see below)
#' 
#' .template: the template for the mrBayes input file (see default_output_template() )
#' 
#' The taxa file consists of a header 'taxon,clade' with each subsequent line containing
#' a taxon,clade pair (separated by a comma)
#' 
#' Each line of the missingclades file consists of the missing clade, the word include or exclude
#' and a list of the reference clades (all separated by commas). Lines containing include specify
#' that a taxon is contained below the MRCA of the reference clades. Lines containing exclude
#' specify that the missing clade cannot attach below the MRCA of the reference clades.
#' 
#' The mrBayes input file is written to base_name.nexus.
#' 
#' @seealso \link{pastis_main} provides a more flexible interface.
#' 
#' \preformatted{
#' }
#' 
#' \strong{PASTIS}: Phylogenetic Assembly with Soft Taxonomic InferenceS? 
#' \preformatted{
#'A bright motmot was acting quite rowdy-
#'weaving and squawking quite loudly
#'"Pastis is delise"
#'he burped with a sneese
#'"but why is everything suddenly cloudy?"} -- \emph{Arne Mooers}
#'
#' @examples
#' # Generate MrBayes input files with constraints
#' \dontrun{
#' data(accipitridaeFullPastis)
#' pastis_simple(accipitridaeFullPastis, base_name="Accipitridae")
#'
#' data(accipitridaeBasicPastis)
#' pastis_simple(accipitridaeBasicPastis, base_name="AccipitridaeBasic")
#' }
#'
#' data(pastis_data_1)
#' pastis_simple(pastis_data_1, base_name="pastis_data_1")
#' unlink("pastis_data_1.nexus")
#'
#' data(pastis_data_2)
#' pastis_simple(pastis_data_2, base_name="pastis_data_2")
#' unlink("pastis_data_2.nexus")
#'
#' data(pastis_data_3)
#' pastis_simple(pastis_data_3, base_name="pastis_data_3")
#' unlink("pastis_data_3.nexus")
#'
#' @inheritParams pastis_main
#' @param base_name The base name for all input files may include a leading directory, but should 
#' not include a trailing .
#' @export
pastis_simple <- function(pastisData=NULL,base_name,paraphyly_constrains=TRUE,monophyly_constrains=TRUE,omit_sequences=FALSE)
{
  
  file_check <- function(base,extension,silent=FALSE)
  {
    file<-paste(base,extension,sep='.')
    if(file.exists(file))
    {
      return(file)
    } else
    {
      if (!silent)
      {
        warning(paste('No',extension,'file was found, proceeding without (this is fine if it is intentional)'))
      }
      return(NA)
    }
  }
  
  return(pastis_main(pastisData=pastisData,
				constraint_tree=file_check(base_name,'tree'),
                taxa_list=file_check(base_name,'taxa'),
                missing_clades=file_check(base_name,'missingclades'),
                sequences=file_check(base_name,'sequences'),
                output_file=paste(base_name,'nexus',sep='.'),
                output_template=file_check(base_name,'template'),
                paraphyly_constrains=paraphyly_constrains,
                monophyly_constrains=monophyly_constrains,
                omit_sequences=omit_sequences))
}

#' Phylogenetic Assembly with Soft Taxonomic Inferences
#' 
#' This function assimilates sequences, taxonomic information and tree constraints 
#' into a mrBayes file. This permits the construction of trees that are 
#' compatible with all of these sources of tdata and contain all known taxa.
#' 
#' This is the main function in pastis which assimilates sequences, taxonomic 
#' information and tree constraints and creates a mrBayes input file. This input 
#' file contains the tree structure specified by constraint_tree with missing taxa 
#' in taxa_list and missing clades in missing_clades added and placed loosely in the 
#' tree using the constraint logic outlined in Thomas et al. MEE (in review) and Jetz et al. (2012 Nature, 491, 444-448). 
#' 
#' See read_input for a description of the required format of the input files. At a 
#' minimum the constraining input tree and taxa list must be provided.
#' 
#' In addition to the input checks conducted by read_input this function also
#' checks for compatibility between the missing genus constraints and constraint
#' tree.
#' 
#' \preformatted{
#'  }
#' 
#' \strong{PASTIS}: Phylogenetic Assembly with Soft Taxonomic InferenceS? 
#' \preformatted{
#'A bright motmot was acting quite rowdy-
#'weaving and squawking quite loudly
#'"Pastis is delise"
#'he burped with a sneese
#'"but why is everything suddenly cloudy?"} -- \emph{Arne Mooers}
#' 
#' @param pastisData Input data object of class pastisData
#' @inheritParams read_input
#' @param output_file The filename for the mrBayes output file (will be overwritten if it exists)
#' @param paraphyly_constrains If TRUE, missing clades are prevented from entering
#' paraphyletic clades.
#' @param monophyly_constrains If TRUE, missing clades are precented from entering monophyletic clades.
#' @param omit_sequences If set to TRUE the sequence file (if any) will be ignored. This is useful for testing the constraints
#' created by pastis as mrBayes runs much quicker without sequence data!
#' @return NULL
#' @seealso \link{pastis_simple} provides a simplified interface to pastis_main.
#' 
#' \link{read_input} describes the required file formats
#' 
#' \link{default_output_template} provides an example of the output template (also the default)
#'
#' @examples
#' \dontrun{
#' # Generate MrBayes input files with constraints
#' data(accipitridaeFullPastis)
#' pastis_main(accipitridaeFullPastis, output_file="Accipitridae.nexus")
#'
#' data(accipitridaeBasicPastis)
#' pastis_main(accipitridaeBasicPastis, output_file="AccipitridaeBasic.nexus")
#' }
#'
#' data(pastis_data_1)
#' pastis_main(pastis_data_1, output_file="pastis_data_1")
#' unlink("pastis_data_1.nexus")
#'
#' data(pastis_data_2)
#' pastis_main(pastis_data_2, output_file="pastis_data_2")
#' unlink("pastis_data_2.nexus")
#'
#' data(pastis_data_3)
#' pastis_main(pastis_data_3, output_file="pastis_data_3")
#' unlink("pastis_data_3.nexus")
#' 
#' @export
pastis_main <- function(pastisData=NULL, constraint_tree,taxa_list,missing_clades=NA,sequences=NA,output_template=NA,output_file='output.nex',paraphyly_constrains=TRUE,monophyly_constrains=TRUE,omit_sequences=FALSE)
{
  #LOGIC
  #The constraints are created by collapsing the tree one cherry at a time and creating 
  #constraints corresponding to the internal node of each collapsed cherry.
  #Missing taxa need be added using this tree traversal process to determine whether 
  #the constraint for an internal node should be hard or partial.
  #Missing genus constraints could be added in a separate process, however by setting them
  #with this process it is possible to simultaneously check that they are compatible with 
  #the constraint tree.
  
  constraints<-c()
  reset_constraint_counters()
  
  #Read and sanitise input
  # GT: omit sequence data and template if omit_sequences=TRUE

  if (is.null(pastisData)) {
	if (omit_sequences)
	{
	input<-read_input(constraint_tree=constraint_tree,taxa_list=taxa_list,missing_clades=missing_clades, NA, NA)
	} else {
	input<-read_input(constraint_tree,taxa_list,missing_clades,sequences,output_template)
	}
	} else {
		if (class(pastisData)!="pastisData") {
		stop("Object is not of class pastisData")
		}
	input <- pastisData 
	}	
	
  constraint_tree<-input$constraint_tree
  taxa_list<-input$taxa_list
  missing_clades<-input$missing_clades
  sequences<-input$sequences
  output_template<-input$output_template
  
  
 # if (omit_sequences)
 # {
 #   sequences=NA
 # }
  #Build clade lists
  #A list with all taxa belonging to each clade
  clade_list = list()
  #A list with all taxa from each clade that are present in the constraint tree
  clade_list_constraining = list()
  for (clade in levels(taxa_list$clade))
  {
    clade_list[[clade]]<-as.character(taxa_list$taxon[taxa_list$clade==clade])
    clade_list_constraining[[clade]]<-intersect(taxa_list$taxon[taxa_list$clade==clade],constraint_tree$tip.label)
    #force missing clades to remain monophyletic
    if (length(clade_list_constraining[[clade]])==0 & length(clade_list[[clade]])>1)
    {
      constraints<-c(constraints,create_constraint('hard',clade_list[[clade]],name_class='missing_clade'))
    }
  }
  
  #Get a list of clades which will be reduced as clades are completed
  remaining_clades<-names(clade_list)
  
  constraining_taxa<-constraint_tree$tip.label
  all_taxa<-as.character(taxa_list$taxon)
  #missing_taxa<-setdiff(all_taxa,constraining_taxa)
  
  #Create a copy of the constraint tree which we will collapse
  tree<-constraint_tree
  
  st<-FALSE
  pendant_edge_index<-0
  #First we have to check pendant edges (in case there are clades with a single taxon
  #in the constraint tree), then we collapse the constraint tree one node at a time
  while(length(tree$tip.label)>1)
  {
    if (pendant_edge_index<length(tree$tip.label))
    {
      pendant_edge_index<-pendant_edge_index+1
      index<-pendant_edge_index
    } else
    {
      #Find a cherry
      edge <- find_cherry(tree)
      n_leaves <- sum(tree$edge[,1]==edge)
      
      #We check for any evidence to support breaking the polytom
      #NB: using a brute force approach for generality. Future versions could also permit
      #missing genera constraints to break polytomies. (i.e. would need to check more than the
      #completeness of each combination)
      found <- FALSE
      if (n_leaves>2)
      {
      for (s in 2:(n_leaves-1))
      {          
        for (combination in combn(seq(n_leaves),s,simplify=FALSE))
        {
          collapsed_temp<-collapse_cherry(tree,edge,combination)
          cherry_taxa <- collapsed_temp[[1]]$tip.sets[[collapsed_temp[[2]]]]
          if (is_complete_clade(cherry_taxa,clade_list_constraining))
          {
            collapsed<-collapse_cherry(tree,edge,combination)
            tree<-collapsed[[1]]
            index<-collapsed[[2]]              
            found<-TRUE
            break
          }
        }
        if (found==TRUE)
        {
          break
        }
      }
      }
      if (found == FALSE)
      {
        collapsed<-collapse_cherry(tree,edge)
        tree<-collapsed[[1]]
        index<-collapsed[[2]]
      }
    }
    
    cherry_taxa <- tree$tip.sets[[index]]
    #Check if the taxa corresponding to this cherry form a complete set of clades
    complete<-is_complete_clade(cherry_taxa,clade_list_constraining)
    
    #If complete we check if there are new clades and can add missing taxa and genera
    #as appropriate
    if(complete)
    {      
      #The clades for the cherry taxa
      clades<-get_clades(cherry_taxa,clade_list_constraining)
      
      #Any new clades (which as complete==TRUE must be complete)
      new_clades<-intersect(clades,remaining_clades)
      
      #Update the list of remaining (uncompleted) clades
      remaining_clades<-setdiff(remaining_clades,new_clades)
      
      #Check whether there are missing taxa to enter in at this point
      cherry_all_taxa<-unique(c(unlist(clade_list[clades]),cherry_taxa))
      includes_missing_taxa <- (length(cherry_taxa)!=length(cherry_all_taxa))
      
      #If a monophyletic or paraphyletic clade is completed and we are preventing the missing clades from
      #entering this part of the tree then we set a hard constraint
      if ((length(clades)==1 & monophyly_constrains) | (length(new_clades)>0 & length(clades) > 1 & paraphyly_constrains))
      {
        constraints<-c(constraints,create_constraint('hard',cherry_all_taxa,name_class='ctree'))
      } else
      {
        #otherwise a soft constraint which permits missing clades and taxa past
        exclude<-setdiff(constraining_taxa,cherry_taxa)
        if (length(exclude)>0)
        {
          constraints<-c(constraints,create_constraint('partial',cherry_all_taxa,exclude=exclude,name_class='ctree'))  
        }
      }      
      
      #Determine if there are any new missing genus constraints
      for (type in c('include','exclude'))
      {
        for (mc in names(missing_clades[[type]]))
        {

          #If all the clades for this missing clade constraint are present then we add constraints
          if (prod(is.element(missing_clades[[type]][[mc]],clades))==1)
          {
           if (type == 'exclude')
            {
              #Soft constraint to prevent taxa in this missing clade descending below this node
              constraints<-c(constraints,create_constraint('partial',
                                                           cherry_all_taxa,
                                                           exclude=clade_list[[mc]],
                                                           name_class=paste('missing_clade_exclude_',mc,sep='')))
              if (!is.element(mc,names(missing_clades[['include']])))
              {
                stop(paste('Missing clade',mc,'reached exclude statement without a remaining include'))
              }
            } else
            {
              #Missing clade taxa will be added at this point and the appropriate soft constraint implemented
              #as the presence of a new clade will be detected further on.
              cherry_all_taxa<-c(cherry_all_taxa,clade_list[[mc]])

              #Soft constraint with cherries + missing genus taxa versus remaining constrained taxa here
              constraints<-c(constraints,create_constraint('partial',
                                                           cherry_all_taxa,
                                                           exclude=setdiff(constraining_taxa,cherry_all_taxa),
                                                           name_class=paste('missing_clade_include_',mc,sep='')))
              if (is.element(mc,names(missing_clades[['exclude']])))
              {
                stop(paste('Missing clade',mc,'reached include statement while exclude statements remain'))
              }
            }
            missing_clades[[type]][[mc]]<-NULL
          }
        }
      }
      names(cherry_all_taxa)<-NULL
      tree$tip.sets[[index]]<-cherry_all_taxa
      
    } else
    {
      #otherwise a soft constraint to enforce the tree constraint but permit missing clades and taxa past
      constraints<-c(constraints,create_constraint('partial',cherry_taxa,exclude=setdiff(constraining_taxa,cherry_taxa),name_class='ctree'))  
    }
  }
  
  #Make a string for the constraints
  constraint_string<-paste(constraints,collapse=';\n')
  constraint_splits<-strsplit(constraints,split=' ')
  constraint_names<-paste(sapply(constraint_splits,function(x) { return(x[2])}),collapse=', ')
  
  constraint_string<-paste(constraint_string,paste('prset topologypr=constraints(',constraint_names,')',sep=''),sep=';\n')
  constraint_string<-paste(constraint_string,';',sep='')
  
  #Make a string for the sequences
  sequence_string = ''


  # GT: set number of characters
  if (!omit_sequences) 
  {
	n_char<-dim(sequences)[2] 
  } else {
	n_char <- 1
  }


  missing_sequence<-paste(rep('?',n_char),collapse='')
  for (taxon in all_taxa){
    if (is.element(taxon,row.names(sequences)))
    {
      sequence_string<-paste(sequence_string,'\n',taxon,paste(sequences[taxon,],collapse=''),sep=' ')
    } else
    {
      sequence_string<-paste(sequence_string,'\n',taxon,' ',missing_sequence,sep='')
    }
    
  }
  
  #Create an output file based on the template
  if (grep('<ntax>',output_template))
  {
    output<-output_template
  } else
  {
    output<-readChar(output_template, file.info(output_template)$size)
  } 
  output<-sub('<ntax>',length(all_taxa),output)
  output<-sub('<nchar>',n_char,output)
  output<-sub('<outputfile>',output_file,output)
  output<-sub('<constraints>',constraint_string,output)
  output<-sub('<sequences>',sequence_string,output)
  
  writeChar(output,output_file)
  
}
#' Called by create_job to read input files 
#' 
#' Reads the specified input files and performs some basic consistency checking
#' between the inputs. 
#' 
#' @param constraint_tree A tree with constraints that are forced to be present in 
#'   all output trees. Either a filename to a nexus file readable by read.tree or a 
#'   ape phylo object.
#' @param taxa_list A list of all taxa and their clades. Either a data frame 
#'   with columns "taxa" and "clade" or a filename for a file readable by 
#'   read.csv with those columns
#' @param missing_clades A file containing missing clades. Each line of the missing clades 
#' file consists of the missing clade, the word "include" or "exclude" and a list of the 
#' reference clades (all separated by commas). Lines containing "include" specify
#' that a taxon is contained below the MRCA of the reference clades. Lines containing "exclude"
#' specify that the missing clade cannot attach below the MRCA of the reference clades.
#' #' @param output_template The filename for a template for the output nexus file. This file should 
#' look like a regular mrBayes input file with special tags replacing content that will be filled by
#' pastis. In particular:
#' 
#' <sequences> will be replaced by the sequences (and should go below the MATRIX line)
#' 
#' <ntax> the number of taxa (i.e. "ntax=<ntax>" must be somewhere in your template)
#' 
#' <nchar> the number of characters
#' 
#' <constraints> the constraints will go here
#' 
#' <outputfile> where the summaries will be written, (i.e. \verb{"sumt filename=<outputfile> burnin ...."} should be in your template) 
#' 
#' see default_output_template for an example (which is used by default)
#' 
#' @param sequences A file with all the available sequence information in fasta format
#'   for details on that format see read.dna in the ape package.
#'   
#' @return A list with the loaded input
#' @export
read_input <- function(constraint_tree,taxa_list,missing_clades=NA,sequences=NA,output_template=NA)
{

  #Load the output template
  if(is.na((output_template)))
  {
    output_template<-default_output_template()
  } else
  {
    if (!file.exists(output_template))
    {
      stop(paste('the output template file (',output_template,') does not exist'))
    }
    #Read template file
    output_template<-readChar(output_template, file.info(output_template)$size)
    #Tags that are required
    tags<-c("<ntax>", "<nchar>", "<sequences>", "<constraints>", "<outputfile>")
    for (tag in tags)
    {
      if (!grep(tag,output_template))
      {
        stop(paste('the output template file is missing the tag: ',tag))
      }
    }
  }
  
  #Load and check the constraint tree
  if(class(constraint_tree)=='character') { constraint_tree=read.tree(constraint_tree) }
  if(class(constraint_tree)!='phylo') { stop('invalid constraint tree specified')}
  constraint_tree <- reorder(constraint_tree,'pruningwise')
  constraint_tree <- add_tip_sets(constraint_tree)
  #Load and check the taxa list
  if(class(taxa_list)=='character') { taxa_list=read.csv(taxa_list) }
  if(class(taxa_list)!='data.frame') {stop('invalid taxa list specified')}
  if(length(intersect(names(taxa_list),c('taxon','clade')))!=2) { stop('invalid column names in taxa list')}
  #A character version of the taxa list
  taxa_list$taxa_char <- as.character(taxa_list$taxon)
  
  #Check compatibility of taxa in constraint tree and taxa list
  if(length(intersect(constraint_tree$tip.label,taxa_list$taxon))!=length(constraint_tree$tip.label)) 
  { 
    missing_taxa = setdiff(constraint_tree$tip.label,intersect(constraint_tree$tip.label,taxa_list$taxon))
    stop(paste("the taxa:",missing_taxa,'is/are present in the constraint tree but not in the taxa list'))
  }
  
  #We need to work out which clades are represented in the constraint tree to sanitise the input
  constraint_taxa<-constraint_tree$tip.label
  constraint_clades<-unique(taxa_list$clade[is.element(taxa_list$taxon,constraint_taxa)])
  #Load the missing clades
  if (!is.na(missing_clades))
  {
    missing_clades_raw<-scan(missing_clades,what='character',sep='\n')
    missing_clades<-list(include=list(),exclude=list())
    
    #Turn the missing clades into a list. Keys are the missing clade names. 
    for (row in 1:length(missing_clades_raw))
    {
      split<-strsplit(missing_clades_raw[row],split=',')[[1]]
      this_clade<-split[1]
      
      #Check that this missing clade is really missing
      if (is.element(this_clade,constraint_clades))
      {
        stop(paste("the clade",this_clade,'is present in both the constraint tree and the missing clade file'))
      }
      
      type<-tolower(split[2])
      cogeners<-split[3:length(split)]
      #Remove this clade from the cogeners if it is included
      cogeners<-setdiff(cogeners,this_clade)
      #Check the missing clade contains taxa
      if (!is.element(this_clade,taxa_list$clade))
      {
        stop(paste('missing clade',this_clade,'is not present in the taxa list'))
      }
      #Check each of the cogeners contains taxa
      if (length(intersect(cogeners,taxa_list$clade))!=length(cogeners))
      {
        discrepancy<-setdiff(cogeners,intersect(cogeners,taxa_list$clade))
        stop(paste('missing clade',this_clade,'has sister clade(s) that are not in the taxa list'))
      }
      
      #Check the constraint type is valid
      if (!is.element(type,c('include','exclude')))
      {
        stop(paste('missing clade',this_clade,'has an invalid constraint type (valid types are include and exclude)'))
      }      
      
      #Check that cogeners aren't missing
      missing_cogeners<-setdiff(cogeners,constraint_clades)
      if (length(missing_cogeners)>0)
      {
        warning(paste('the missing clade constraint for',this_clade,'refers to missing clade(s):',paste(missing_cogeners,collapse=','),', this reference will be ignored when positioning this clade'))
      }
      cogeners<-intersect(cogeners,constraint_clades)
      missing_clades[[type]][[this_clade]]<-cogeners
    }
  } else
  {
    missing_clades<-list(include=list(),exclude=list())
  }
  #Check if all clades (in the taxa list) are represented in the missing_clades OR the constraint_tree
  
  
  
  #Load and check sequences
  
  if (is.na(sequences))
  {
    sequences<-matrix(c(1,1))
  } else
  {
    e<-try(sequences <- read.dna(sequences, format="fasta", as.character=TRUE))
    
    if(class(e)=='try-error')
    {
      stop(paste('could not load sequence data from',sequences,'. See read.dna for a description of the fasta format.'))
    }
    
  }

	pastis_object <- list(constraint_tree=constraint_tree, taxa_list=taxa_list, missing_clades=missing_clades, sequences=sequences, output_template=output_template)
	attr(pastis_object, "class") <- "pastisData"
  
  return(pastis_object)
}

# Clade completeness check
# 
# Determines all clades represented by a set of taxa and checks whether all 
# species in those clades are in that set of taxa
# 
# @inheritParams get_clades
# @return logical, whethe the clades were completely represented
is_complete_clade<-function(set,clades)
{
  present_clades<-get_clades(set,clades)
  required_taxa<-unlist(clades[present_clades])
  return(length(intersect(required_taxa,set))==length(required_taxa))
}

# Returns the clades for a set of taxa
# 
# This function determines which clades are represented by a set of taxa
# 
#  @param set A vector of taxa
#  @param clades A list of clades, each element contains a vector of taxa
#  belonging to that clade
#  @return A vector of clades represented by the taxa
get_clades <-function(set,clades)
{
  out=c()
  for (clade in names(clades))
  {
    if (length(intersect(set,clades[[clade]]))>0)
    {
      out <- c(out,clade)
    }
  }
  return(out)
}

e<-new.env()
#constraint_id<-0
#constraint_class_id<-list()

# Reset the constraint counters 
# 
# Resets the constraint counters that are used to give unique names to
# the constraints
#
# @return Nothing 
reset_constraint_counters <- function()
{
  e$constraint_id<-0
  e$constraint_class_id<-list()
}

# Constraint creator
# 
# Creates mrBayes style constraints for the supplied list of species. 
# Constraints are named according to either the supplied name or by the
# supplied name_class. If name_class is used each new constraint is 
# incremented.
# 
# @param type constraint type either "hard","negative" or "partial"
# @param include the first list of taxa
# @param exclude the second list of taxa (not used for hard constraints)
# @param name the name for the constraint. If none specified they are automatically assigned
# @param name_class alternatively automatically incremented names
# @return A string containing the constraint
create_constraint <- function(type,include,exclude=c(),name=NA,name_class=NA) {
  if (length(include) == 1)
  {
    return()
  }
  type<-tolower(type)
  if (!is.element(type,c('hard','negative','partial')))
  {
    stop('invalid constraint type specified (valid types are hard, negative and partial)')
  }
  if (!is.na(name_class))
  {
    if(!is.element(name_class,names(e$constraint_class_id)))
    {
      e$constraint_class_id[name_class]<-0
    }
    e$constraint_class_id[[name_class]]<-e$constraint_class_id[[name_class]]+1
    name<-paste(name_class,e$constraint_class_id[[name_class]],sep='_')
  }
  if (is.na(name))
  {
    name <- as.character(e$constraint_id)
    e$constraint_id<-e$constraint_id+1
  }
  out <- paste('constraint',name)
  out <- paste(c(out,type,'=',sort(include)),collapse=' ')
  if (length(exclude)>0)
  {
    out <- paste(c(out,':',sort(exclude)),collapse=' ')
  }
  return(out)  
}


# Add sets of tips to a phylo tree
# 
# @param tree an ape tree
# @return The modified ape tree
add_tip_sets <- function(tree)
{
  tip.sets <- list()
  for (species in 1:length(tree$tip.label))
  {
    tip.sets[[species]]<- tree$tip.label[species]
  }
  tree$tip.sets<-tip.sets
  return(tree)
}

# Find a cherry (may be a polytomy)
#
#
# @param tree An ape tree
# @return An index to the edge above the cherry
find_cherry<- function(tree)
{
  n_species <- length(tree$tip.label)
  for (edge in n_species+1:max(tree$edge))
  {
    select <- tree$edge[,1]==edge
    if (sum(tree$edge[select,2]>n_species)==0)
    {
      break
    }
  }  
  return(edge)
}

# Collapse a cherry (may be a polytomy)
#
#
# @param tree An ape tree
# @param edge The edge to collapse (NA to choose the first cherry)
# @param combine For polytomies use this to choose which tips to combine. Specified as a list
#   containing elements between 1 and the number of tips in the polytomy. Tips are numbered according
#   to their order in tree$edge. Omit or set to NA to combine all tips.
# @return A list containing: 1) the tree with a collapsed cherry 
#   2) an index to the edge that was collapsed
collapse_cherry <- function(tree,edge=NA,combine=NA)
{
  
  #If no edge specified, find the first cherry to collapse
  if (is.na(edge))
  {
    n_species <- length(tree$tip.label)
    edge<-find_cherry(tree)
    if (!is.na(combine))
    {
      stop('You can only specify which branches to combine if you have specified an edge')
    }
  }
  
  #Check that the tree has tip name sets, if not add them
  if (!is.element('tip.sets',names(tree)))
  {
    tree <- add_tip_sets(tree)
  }
  
  #The rows corresponding to the cherry edges
  rows<-tree$edge[,1]==edge
  
  #Check that valid edges are set to be combined
  if (!is.na(combine)[1])
  {
    if (max(combine)>sum(rows) | min(combine)<1)
    {
      stop('Invalid edges specified in combine')
    }
    if (length(unique(combine))!=length(combine))
    {
      stop('Elements in combine must be unique')
    }
    if (length(combine)<2)
    {
      stop('You must specify multiple edges to combine')
    }
  }
  
  #If we are combining all edges
  if (length(combine)==sum(rows) | is.na(combine)[1])
  {
    #The cherry leaves
    leaves<-sort(tree$edge[rows,2])
    #The new 'leaf' will have the id of the first cherry leaf
    new_leaf<-leaves[1]
    #Setup a new edge matrix
    new_edge<-tree$edge
    #Relabel the new 'leaf'
    new_edge[new_edge[,2]==edge,2]=new_leaf
    #Remove the cherry edges
    new_edge<-new_edge[!rows,]
    
    #Remove the superfluous internal edge
    select <- new_edge[,2]>= edge
    new_edge[select,2]<-new_edge[select,2]-1
    select <- new_edge[,1]>= edge
    new_edge[select,1]<-new_edge[select,1]-1
    
    #Renumber the leaves to account for the removed leaves
    for (removed_leaf in seq(length(leaves),2,by=-1))
    {
      select <- new_edge[,2]>= leaves[removed_leaf]
      new_edge[select,2]<-new_edge[select,2]-1
      select <- new_edge[,1]>= leaves[removed_leaf]
      new_edge[select,1]<-new_edge[select,1]-1
      
    }
    
    #New tip labels
    new_tip.label <- tree$tip.label[-leaves[2:length(leaves)]]
    new_tip.label[leaves[1]]<-'merge'
    
    #New tip sets
    new_tip.sets <- tree$tip.sets
    new_tip.sets[[leaves[1]]]<-unique(unlist(new_tip.sets[leaves]))
    new_tip.sets <- new_tip.sets[-leaves[2:length(leaves)]]
    
    new_tree <- tree
    new_tree$tip.label<-new_tip.label
    new_tree$edge <- new_edge
    new_tree$Nnode<-new_tree$Nnode-1
    new_tree$tip.sets <- new_tip.sets
  } else
  {
    rows<-which(tree$edge[,1]==edge)[combine]
    #The cherry leaves
    leaves<-sort(tree$edge[rows,2])
    #Setup a new edge matrix
    new_edge<-tree$edge
    #Remove the cherry edges
    new_edge<-new_edge[-rows[-1],]
    
    #Renumber the leaves to account for the removed leaves
    for (removed_leaf in seq(length(leaves),2,by=-1))
    {
      select <- new_edge[,2]>= leaves[removed_leaf]
      new_edge[select,2]<-new_edge[select,2]-1
      select <- new_edge[,1]>= leaves[removed_leaf]
      new_edge[select,1]<-new_edge[select,1]-1
      
    }
    
    #New tip labels
    new_tip.label <- tree$tip.label[-leaves[2:length(leaves)]]
    new_tip.label[leaves[1]]<-'merge'
    
    #New tip sets
    new_tip.sets <- tree$tip.sets
    new_tip.sets[[leaves[1]]]<-unique(unlist(new_tip.sets[leaves]))
    new_tip.sets <- new_tip.sets[-leaves[2:length(leaves)]]
    
    new_tree <- tree
    new_tree$tip.label<-new_tip.label
    new_tree$edge <- new_edge
    new_tree$Nnode<-new_tree$Nnode
    new_tree$tip.sets <- new_tip.sets
    
  }
  
  return(list(new_tree,leaves[1]))
}

#' mrBayes output interrogator (COnstraint CHecker)
#' 
#' This function examines mrBayes output from a pastis run to examine where in
#' the original constraint tree missing taxa have been placed
#' 
#' An analysis of the placement of each taxon not contained in the constraint tree is 
#' conducted. For each such taxon the edge lengths in the constraint tree are adjusted
#' according to the parameter edge_scaling. If simple_edge_scaling is TRUE, the edges
#' will have length 0 if a taxon is never descendant from a tree and 1 if it is descendant 
#' in at least one tree. If FALSE, indicate the proportion of sampled trees in which the 
#' taxon is descendant from that edge. 
#' 
#' The output from this function is useful for checking that missing taxa are placed in
#' appropriate positions relative to the original constraint tree.
#' 
#' Note that this routine has not been optimised and slow (possibly unacceptable) 
#' performance is to be expected with large trees and/or large posterior samples.
#' 
#' @param constraint_tree the constraint tree used with pastis (either the filename or ape phylo tree)
#' @param mrbayes_output the mrBayes output .t file (either the filename or ape multiPhylo tree object)
#' @param simple_edge_scaling boolean, see function description
#' @param species_set if specified, this should be a list of species and output trees will only be generated for these species (the default is all missing species)
#' @return NULL. A file for each missing taxon is created in the current directory
#' 
#' @examples
#' data(pastis_data_3_trees)
#'
#' \dontrun{
#'
#' # Check constraints for all missing taxa (takes ~6 seconds to run: sped up by Anonymous Reviewer 2)
#' conch(pastis_data_3_trees[[1]], pastis_data_3_trees[[2]])
#'
#' }
#'
#' # Check constraints for missing taxon "a_4"
#' conch(pastis_data_3_trees[[1]], pastis_data_3_trees[[2]], species_set="a_4")
#' unlink("taxonposition_a_4.tree")
#' 
#' @export
conch <- function(constraint_tree,mrbayes_output,simple_edge_scaling=TRUE,species_set=NA)
{
  
  #Load and check the constraint tree
  if(class(constraint_tree)=='character') { constraint_tree=read.tree(constraint_tree) }
  if(class(constraint_tree)!='phylo') { stop('invalid constraint tree specified')}
  constraint_tree <- reorder(constraint_tree,'cladewise')
  
  
  #Get list of all constraining taxa
  constraining_taxa<-constraint_tree$tip.label
  
  #Load trees
  if(class(mrbayes_output)=="character") { mrbayes_trees<-read.nexus(mrbayes_output) } else { mrbayes_trees<-mrbayes_output }
  if(class(mrbayes_trees)!='multiPhylo') { stop('invalid tree(s) specified')}
  
  missing_taxa<-setdiff(mrbayes_trees[[1]]$tip.label,constraining_taxa)
  if (is.na(species_set))
  {
    species_set<-missing_taxa
  }

  #Transform multiPhylo object into list (thanks to Anonymous MEE reviewer 2 for the suggestion)
  mrbayes_trees = .compressTipLabel(mrbayes_trees)
  mrbayes_label = attr(mrbayes_trees, "TipLabel")
  mrbayes_trees = unclass(mrbayes_trees)  

  #For each missing taxon M
  count<-1
  for (taxon in species_set)
  {
    print(paste('Processing species ',count,' of ',length(species_set),': ',taxon,sep=''))
    count<-count+1
    #Set edges to 1
    constraint_tree$edge.length<-c(1:dim(constraint_tree$edge)[1])
    constraint_tree$edge.length[]<-0
    constraint_tree$root.edge<-0

    for (tree_index in 1:length(mrbayes_trees))
    {
	  this_tree<-mrbayes_trees[[tree_index]]
	  this_tree$tip.label = mrbayes_label  
      #Copy & trim tree to include only constraining taxa and M
      this_tree<-drop.tip(this_tree,setdiff(this_tree$tip.label,c(constraining_taxa,taxon)))
      #Determine descendants of M's pendant edge
      edge<-this_tree$edge[,2]==which(this_tree$tip.label==taxon)
      this_tree<-extract.clade(this_tree,this_tree$edge[this_tree$edge[,2]==which(this_tree$tip.label==taxon),1])

      mrca_descendants <- intersect(this_tree$tip.label,constraining_taxa)
      if (length(mrca_descendants)==1)
      {
        mrca_edge<-which.edge(constraint_tree,mrca_descendants)[1]
      } else
      {
        mrca_node<-constraint_tree$edge[which.edge(constraint_tree,mrca_descendants)[1],1]
        mrca_edge<-which(constraint_tree$edge[,2]==mrca_node)
      }
      if (simple_edge_scaling)
      {
        constraint_tree$edge.length[mrca_edge]<-1
      } else
      {
        constraint_tree$edge.length[mrca_edge]<-constraint_tree$edge.length[mrca_edge]+1
      }
    }
    write.tree(constraint_tree,file=paste('taxonposition_',taxon,'.tree',sep=''))
  }
}
#conch('1.tree','1.nexus.t')  
  
#' A sample output template for pastis_main and pastis_simple
#' 
#' This is the default output template filled in by pastis_main and 
#' pastis_simple to create the input file for mrBayes. 
#' 
#' If you want to change the parameters in the mrBayes nexus file 
#' created by pastis then it is easier to create your own template rather than 
#' editing the created nexus files manually (which you would have to do again
#' if you reran pastis).
#' 
#' This is the default template used by pastis and is a good starting point for
#' creating your own template. To have 
#' a look at this template try:
#' 
#' \verb{
#' cat(default_output_template())
#' }
#'
#' Once you understand the format, you can create your own template as a string
#' and pass it to pastis.
#'
#' @return The default output template
#' @export
default_output_template <- function()
{
return('
Begin DATA; 
Dimensions ntax=<ntax> nchar=<nchar>;
Format datatype=DNA gap=- missing=?; 
Matrix 

<sequences>   
  ; 


begin MRBAYES; 

unlink shape=(all) tratio=(all) statefreq=(all) revmat=(all) pinvar=(all); 
<constraints>
  
  
  set usebeagle=no Beaglesse=no; 

prset brlenspr=clock:birthdeath; 
prset Extinctionpr = Fixed(0); 
prset Speciationpr=exponential(1); 
prset clockvarpr=ibr; 
prset ibrvarpr=exponential(10); 
mcmcp nruns=1 nchains=1 ngen=50000000 samplefreq=1000; 
mcmc; 

sumt filename=<outputfile> burnin=5000000 contype=halfcompat;

end; ')
}



#a<-pastis_simple('Accipitridae')

#' @name accipitridaeBasicPastis
#' @title Accipitridae basic data
#' @description Constraint tree and taxon list for Accipitridae
#' @docType data
#' @usage data(accipitridaeBasicPastis)
#' @author Gavin Thomas, 2013-07-22
NULL

#' @name accipitridaeFullPastis
#' @title Accipitridae full data
#' @description Constraint tree, taxon list, missing clades, sequence data and template for Accipitridae
#' @docType data
#' @usage data(accipitridaeFullPastis)
#' @author Gavin Thomas, 2013-07-22
NULL

#' @name pastis_data_1
#' @title Example data 1
#' @description Constraint tree and taxon list examples
#' @docType data
#' @usage data(pastis_data_1)
#' @author Gavin Thomas, 2013-07-22
NULL

#' @name pastis_data_2
#' @title Example data 2
#' @description Constraint tree and taxon list examples
#' @docType data
#' @usage data(pastis_data_2)
#' @author Gavin Thomas, 2013-07-22
NULL

#' @name pastis_data_3
#' @title Example data 3
#' @description Constraint tree, taxon and missing clade examples
#' @docType data
#' @usage data(pastis_data_3)
#' @author Gavin Thomas, 2013-07-22
NULL

#' @name pastis_data_3_trees
#' @title Example output
#' @description [[1]] Constraint tree and [[2]] distributions of trees produced in MrBayes using pastis_data_3
#' @docType data
#' @usage data(pastis_data_3_trees)
#' @author Gavin Thomas, 2013-07-22
NULL