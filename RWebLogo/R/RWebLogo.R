# TODO list:
# color TODO Specify symbol colors, e.g. --color black AG 'Purine' --color red TC 'Pyrimidine'
# default.color TODO Symbol color if not otherwise specified. e.g. 'red', '#FF0000'.

# require(jpeg)
# require(raster)

# Alignment formats supported
.DATATYPES = c('clustal', 'fasta', 'plain', 'msf', 'genbank', 'nbrf', 'nexus', 'phylip',
               'stockholm', 'intelligenetics', 'table', 'array','transfac')

# File output formats supported
.FORMATS = c('eps', 'png', 'pdf', 'jpeg', 'svg', 'logodata', 'png_print')

# Different sequence types that can be used 
.SEQTYPE = c('protein', 'rna', 'dna')

# Supported Y-axis units 
.UNITS = c('bits', 'nats', 'digits', 'kT', 'kJ/mol', 'kcal/mol', 'probability')

# Background sequence probabilities
.COMP =  c('auto', 'equiprobable', 'none', 
           'H. sapiens', 'E. coli', 'S. cerevisiae', 'C. elegans', 'D. melanogaster', 'M. musculus', 'T. thermophilus')

# Size of the plot
.SIZE = c('small', 'medium', 'large', 'xlarge', 'xxlarge')

# Different supported color schemes 
.COLOR_SCHEME = c('auto', 'base pairing', 'charge', 'chemistry','chemistry2', 'chemistry3', 'classic', 'hydrophobicity', 'monochrome')

# Write lines using sprintf
.writeLines = function(v, x, ...){
  if(v) writeLines(sprintf(x, ...))
}

.python.valid <- function(){
  ver = tryCatch({
    system2('python', '-V', stderr=T, stdout=T)
  }, error = function(e) {
    ""
  })
  if(ver == "") return(FALSE)
  
  return(TRUE)
}

# Add to command vector
.cmd.add = function(cmd, opt, val=NULL, class='other', has.val=T){
  vn =  as.character( as.list( sys.call() )[[4]] )
  .check.class(val, class, val.name=vn)
  if(is.logical(val)){
    if(val == TRUE) val = 'YES'
    if(val == FALSE) val = 'NO'
  }
  
  if(has.val){
    #val = gsub(' ', replacement='\\\\ ', val)
    val = paste0("'", val, "'")
    t = paste(opt, val)
  }else{
    t = opt
  }
  
  return(c(cmd, t))
}

# Check if value is in the allowed list of values
.check <- function(val, allowed){
  val.name = as.character( as.list( sys.call() )[[2]] )
  if(! val %in% allowed) 
    stop(sprintf('Invalid "%s" value, please specifity one of the following: %s', val.name, paste0(allowed, collapse=', ')))
}

# Ensure class of value is correct
.check.class <- function(val, class, val.name){
  v = c("numeric", "character", "logical", "integer")
  is.valid = TRUE
  if(class == v[1]){
    if(!is.numeric(val)) is.valid = FALSE
  }else if(class == v[2]){
    if(!is.character(val)) is.valid = FALSE
  }else if(class == v[3]){
    if(!is.logical(val)) is.valid = FALSE
  }else if(class == v[4]){
    if(! val%%1==0) is.valid = FALSE
  }
  if(!is.valid){
    stop(sprintf("Invalid '%s' value. Must be %s.", val.name, class ))
  }
}

#' Plot a sequence logo
#'
#' This function will plot a sequence logo given aligned sequences. For more details on the parameters, see the WebLogo user manual at \url{http://weblogo.threeplusone.com/manual.html}
#'
#' @param seqs Aligned sequences as an R character vector. Sequences must all have the same length. Alternatively, you can provide a file containing your sequence alignment using \code{file.in}.
#' @param file.in A file containing your sequence alignment in one of the following formats: clustal, fasta, plain, msf, genbank, nbrf, nexus, phylip, stockholm, intelligenetics, table, array, transfac. This file is only to be provided if you are not inputting data with 'seqs'. To set your data format, see the \code{datatype} option.
#' @param open Open the generated logo file? (default: TRUE).
#' @param verbose Write status messages to the R console? (default: TRUE).
#' @param return.cmd Logical indicating if RWebLogo should return the WebLogo command generated (default: FALSE).
#' @param datatype Type of multiple sequence alignment or position weight matrix file: ('clustal', 'fasta','plain', 'msf', 'genbank', 'nbrf', 'nexus', 'phylip', 'stockholm', 'intelligenetics', 'table', 'array', 'transfac'). You usually don't need to specify this, as weblogo will try figure out the format of your file.
#' @param file.out Output file. For example, /path/to/dir/mylogo.pdf. By default this is your working directory + RWebLogo. + selected \code{format}.
#' @param format Format of output: 'pdf' (default) 'eps', 'png', 'jpeg', 'svg'.
#' @param sequence.type The type of sequence data: 'protein', 'rna' or 'dna'.
#' @param alphabet The set of symbols to count, e.g. 'AGTC'. All characters not in the alphabet are ignored. If neither the alphabet nor sequence-type are specified then weblogo will examine the input data and make an educated guess. 
#' @param units  A unit of entropy ('bits' (default), 'nats', 'digits'), or a unit of free energy ('kT', 'kJ/mol', 'kcal/mol'), or 'probability' for probabilities.
#' @param composition The expected composition of the sequences: 'auto' (default), 'equiprobable', 'none' (do not perform any compositional adjustment), a CG percentage, a species name ('H. sapiens', 'E. coli', 'S. cerevisiae', 'C. elegans', 'D. melanogaster', 'M. musculus', 'T. thermophilus'), or an explicit distribution as a named numerical vector (e.g. c(A=10, C=40, G=40, T=10)). The automatic option uses a typical distribution for proteins and equiprobable distribution for everything else.
#' @param weight The weight of prior data. Default depends on alphabet length.
#' @param first.index Index of first position in sequence data (default: 1).
#' @param lower Lower bound index of sequence to display.
#' @param upper Upper bound index of sequence to display.
#' @param ignore.lower.case Disregard lower case letters and only count upper case letters in sequences.
#' @param reverse reverse sequences.
#' @param complement complement DNA sequences.
#' @param size Specify a standard logo size: 'small', 'medium', 'large' (default).
#' @param stacks.per.line Maximum number of logo stacks per logo line (default: 40).
#' @param title Logo title text.
#' @param label A figure label, e.g. '2a'.
#' @param show.xaxis Display sequence numbers along x-axis? (default: TRUE).
#' @param xlabel X-axis label.
#' @param annotate  A comma separated list or vector of custom stack annotations, e.g. '1,3,4,5,6,7' or c(1,3,4,5,6,7).  Annotation list must be same length as sequences.
#' @param yaxis Height of yaxis in units (default: maximum value with uninformative prior).
#' @param show.yaxis Display entropy scale along y-axis? (default: TRUE).
#' @param ylabel Y-axis label (default: depends on plot type and units).
#' @param show.ends Label the ends of the sequence? (default: FALSE).
#' @param fineprint The fine print text at the bottom right corner (default: blank).
#' @param ticmarks Distance between ticmarks (default: 1.0).
#' @param errorbars Display error bars? (default: FALSE).
#' @param reverse.stacks Draw stacks with largest letters on top? (default: TRUE).
#' @param color.scheme Specify a standard color scheme ('auto', 'base pairing', 'charge', 'chemistry', 'classic', 'hydrophobicity', 'monochrome').
#' @param stack.width Width of a logo stack (default: 10.8).
#' @param aspect.ratio Ratio of stack height to width (default: 5).
#' @param box Draw boxes around symbols? (default: FALSE).
#' @param resolution Bitmap resolution in dots per inch (DPI). Low resolution bitmaps with DPI<300 are antialiased  (default: 300 DPI).
#' @param scale.width Scale the visible stack width by the fraction of symbols in the column?  i.e. columns with many gaps of unknowns are narrow (default: TRUE).
#' @param rotate.numbers Rotate values of x-axis? (default: FALSE).
#' @param hide.tics Hide tic marks? (default: FALSE).
#' 
#' @export
#' 
#' @import findpython
#' 
#' @examples
#' # Make a sequence logo using an external alignment file format 
#' # In this example we'll use the EMBOSS alignment format or msf
#' # However, you can use any format supported by WebLogo e.g. fasta
#' fpath = system.file("extdata", "example_data.msf", package="RWebLogo")
#' weblogo(file.in=fpath)
#' # Now for an example using an alignment as an R character vector
#' aln <- c('CCAACCCAA', 'CCAACCCTA', 'AAAGCCTGA', 'TGAACCGGA')
#' # Simple WebLogo
#' weblogo(seqs=aln)
#' # Lets get rid of those ugly error bars and add some text!
#' weblogo(seqs=aln, errorbars=FALSE, title='Yay, No error bars!', 
#'         fineprint='RWebLogo 1.0', label='1a')
#' # We can also change the format of the output like this
#' weblogo(seqs=aln, format='png', resolution=500)
#' # You can change the axis labels like this
#' weblogo(seqs=aln, xlabel='My x-axis', ylabel='Awesome bits')
#' # You get the idea! See ?weblogo for more awesome options! 
weblogo <- function(seqs, file.in,
                open=TRUE, verbose=TRUE, return.cmd=F,
                datatype='plain',
                file.out,
                format='pdf', # eps (default), png, png_print, pdf, jpeg, svg, logodata
                sequence.type='protein', #The type of sequence data: 'protein', 'rna' or 'dna'.
                alphabet, units='bits', composition='auto',
                weight,
                first.index, lower, upper, 
                
                ignore.lower.case = TRUE, reverse, complement,
                
                size = 'large', stacks.per.line = 40, title, label,
                show.xaxis = TRUE, xlabel, annotate, 
                yaxis, show.yaxis=TRUE, ylabel,
                show.ends = FALSE, fineprint = '', ticmarks = 1.0, errorbars = FALSE,
                reverse.stacks = TRUE, 

                color.scheme = 'auto',
                stack.width=10.8, aspect.ratio=5, box=FALSE, resolution=300, scale.width=TRUE,
                
                # New options
                rotate.numbers = FALSE, hide.tics = FALSE
                ){
  
  
  if(!can_find_python_cmd(minimum_version='2.6', required_modules='numpy' )){
    writeLines('Please install and add Python (>=2.6) to your path')
    return(1)
  }
  
  v = verbose
  
  if(missing(file.in) & missing(seqs))
    stop('No inputed data detected. Please provide sequences or an input file')
  
  if(!missing(file.in) & !missing(seqs))
    stop('Please provide either sequences or an input file, but not both!')
  
  if(missing(file.out)){
    file.out=file.path(getwd(), sprintf('RWebLogo.%s', format))
  }
  
  if(!missing(seqs)){
    if(class(seqs) != "character") 
      stop("Invalid sequence data: 'seqs' should be a character vector, please try again")
    
    if(length(seqs) == 0)
      stop("Invalid sequence data: 'seqs' cannot be empty, please try again")
    
    len = sapply(seqs, nchar) 
    if( any(len != len[1]))
      stop("Invalid sequence data: all sequences must have the same number of characters")
  }
  
  if(missing(file.in)){
    file.in = file.path(tempdir(), 'tmp.logo')
    writeLines(seqs, file.in)
  }else if(!file.exists(file.in)){
    stop(sprintf('The input file %s does not exist, please use a valid input file!', file.in))
  }
  
  # /Users/omarwagih/Development/RWebLogo/inst/extdata/weblogo-3.3/weblogo'
  exec = file.path( system.file("extdata", "weblogo-3.3", package="RWebLogo"), 'weblogo')
  z = c(exec,
        sprintf('< %s > %s', file.in, file.out))
  
  if(!missing(datatype)){
    .check(datatype, .DATATYPES)
    z = .cmd.add(z, '--datatype', datatype)
  }
  
  .check(format, .FORMATS)
  z = .cmd.add(z, '--format', format, 'character')
  
  ##################################################################### 
  # Logo Data Options:
  ##################################################################### 
  
  if(!missing(sequence.type)){
    .check(sequence.type, .SEQTYPE)
    z = .cmd.add(z, '--sequence-type', sequence.type, 'character')
  }
  
  if(!missing(alphabet)){
    # Get unique chars 
    alphabet = toupper( unique(strsplit(alphabet, '')[[1]]) )
    # Find number of letters  
    nl = sum(grepl(pattern='[A-Za-z]', alphabet))
    if(nl == 0) stop('Alphabet must have at least one letter from A-Z')
    if(length(alphabet) < 2) stop('Alphabet must have at least two letters')
    
    alphabet = paste0(alphabet, collapse='')
    z = .cmd.add(z, '--alphabet', alphabet, 'character')
  }
  
  if(!missing(units)){
    .check(units, .UNITS)
    z = .cmd.add(z, '--units', units, 'character')
  }
  
  if(!missing(composition)){
    if(is.character(composition)){
      .check(composition, .COMP)
    }else if(!is.numeric(composition) | any(is.null(names(composition)))){
      stop("Composition must be a named numerical vector e.g. c(A=10, C=40, G=40, T=10)")
    }else{
      composition = paste0(sprintf("%s:%s", names(composition), composition), collapse=',')
      composition = paste0("{", composition, "}")
    }
    
    
    z = .cmd.add(z, '--composition', composition, 'character')
  }
  
  if(!missing(weight))
    z = .cmd.add(z, '--weight', weight, 'numeric')
  
  if(!missing(first.index))
    z = .cmd.add(z, '--first-index', first.index, 'numeric')
  
  if(!missing(lower))
    z = .cmd.add(z, '--lower', lower, 'integer')
  
  if(!missing(upper))
    z = .cmd.add(z, '--upper', upper, 'integer')
  
  ##################################################################### 
  # Transformations: Optional transformations of the sequence data.
  ##################################################################### 

  if(!missing(ignore.lower.case))
    z = .cmd.add(z, '--ignore-lower-case', ignore.lower.case, 'logical', has.val=F)
  
  if(!missing(reverse))
    z = .cmd.add(z, '--reverse', reverse, 'logical', has.val=F)
  
  if(!missing(complement))
    z = .cmd.add(z, '--complement', complement, 'logical', has.val=F)
  
  ##################################################################### 
  # Logo Format Options: These options control the format and display of the logo.
  ##################################################################### 

  .check(size, .SIZE)
  z = .cmd.add(z, '--size', size, 'character')
  
  if(!missing(stacks.per.line))
    z = .cmd.add(z, '--stacks-per-line', stacks.per.line, 'numeric')
  
  if(!missing(title))
    z = .cmd.add(z, '--title', title, 'character')
  if(!missing(label))
    z = .cmd.add(z, '--label', label, 'character')
  
  if(!missing(show.xaxis))
    z = .cmd.add(z, '--show-xaxis', show.xaxis, 'logical')
  
  if(!missing(xlabel))
    z = .cmd.add(z, '--xlabel', xlabel, 'character')
  
  if(!missing(annotate)){
    annotate = paste0(annotate, collapse=',')
    z = .cmd.add(z, '--annotate', annotate, 'character')
  }
  
  if(!missing(yaxis))
    z = .cmd.add(z, '--yaxis', yaxis, 'numeric')
  
  if(!missing(show.yaxis))
    z = .cmd.add(z, '--show-yaxis', show.yaxis, 'logical')
  
  if(!missing(ylabel))
    z = .cmd.add(z, '--ylabel', ylabel, 'character')
  
  if(!missing(show.ends))
    z = .cmd.add(z, '--show-ends', show.ends, 'logical')
  
  z = .cmd.add(z, '--fineprint', fineprint, 'character')
  
  if(!missing(ticmarks))
    z = .cmd.add(z, '--ticmarks', ticmarks, 'numeric')
  
  if(!missing(errorbars))
    z = .cmd.add(z, '--errorbars', errorbars, 'logical')
  
  if(!missing(reverse.stacks))
    z = .cmd.add(z, '--reverse-stacks', reverse.stacks, 'logical')
  

  ##################################################################### 
  # Color Options: Colors can be specified using CSS2 syntax. e.g. 'red', '#FF0000', etc.
  ##################################################################### 

  if(!missing(color.scheme)){
    .check(color.scheme, .COLOR_SCHEME)
    z = .cmd.add(z, '--color-scheme', color.scheme, 'character')
  }
  
#   if(!missing(color))
#     z = .cmd.add(z, '--color', color, 'character')
#   
#   if(!missing(default.color))
#     z = .cmd.add(z, '--default-color', default.color, 'character')
  
  # Advanced Format Options: These options provide fine control over the display of the logo.
  if(!missing(stack.width))
    z = .cmd.add(z, '--stack-width', stack.width, 'numeric')
  
  if(!missing(aspect.ratio))
    z = .cmd.add(z, '--aspect-ratio', aspect.ratio, 'numeric')
  
  if(!missing(box))
    z = .cmd.add(z, '--box', box, 'logical')
  
  if(!missing(resolution))
    z = .cmd.add(z, '--resolution', resolution, 'integer')
  
  if(!missing(scale.width))
    z = .cmd.add(z, '--scale-width', scale.width, 'logical')

  if(!missing(rotate.numbers))
    z = .cmd.add(z, '--rotate-numbers', rotate.numbers, 'logical')

  if(!missing(hide.tics)){
    if(hide.tics) z = .cmd.add(z, '--tic-length', 0.01, 'numeric')
  }
  
  command = paste(z, collapse=' ')
  rt = system(command)
  
  status = attr(rt, 'status')
  if(length(status) > 0){
    if(status == 2) stop('RWebLogo failed to plot due to an error from WebLogo')
  }
  
  .writeLines(v, 'Sequence logo saved to %s', file.out)
  if(open){
    if(interactive()){
      browseURL(file.out)
    }else{
      warning('Non interactive system, cannot open plotted sequence logo.')
    }
  }
  
  if(return.cmd){
    return(command)
  }
}