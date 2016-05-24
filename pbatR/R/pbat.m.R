####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     06/??/2005                                          #
# MODIFIED:    01/18/2006                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  The pbat.m(...) interface to pbat.                              #
####################################################################


PBAT.M.DEBUG <- FALSE;


strsplitFix <- function( x, split ) {
  if( length(x) > 1 ) stop( "strSplitFix(...) only works on a single string." );
  if( length(x)==0 || x=="" ) return("");

  ##print( paste( "TO BE SPLIT: (", x, "), (", split, ")", sep="") );

  res=unlist( strsplit( x, split, fixed=TRUE ) ); # split, return as vector of strings
  ##print( "SPLIT" );

  return( res[res!=""] ); # eliminate any empty strings!
}

####################################################################
# pbat.m(...)    STUPID PRINT METHOD IN CLASSES...                 #
# DESCRIPTION: Pretty model building for pbat!                     #
# PARAM  formula  best described in the examples :)                #
#        ...      other options for pbat (we will enforce only     #
#                 possible options later, so, for instance,        #
#                 you cannot redo haplos=()                        #
# EXAMPLES:                                                        #
#  pbatm( phenos1 ~ preds1 ); # basic model, does all snps (if none specified)
#  pbatm( phenos1 ~ mi(preds1) ); # does all snps, the mi() tells it should be a marker interaction
#  pbatm( time & censor ~ 1 * preds1^3 + preds2   | snp1 | snp2 | snp3 / group, temp );  # single snp analysis (b/c there's only one snp), stratified by group (presence of censor auto-indicates log-rank analysis)
#  pbatm( time & censor ~ mi(preds1^2) + mi(preds2^3) | block1snp1 + block1snp2 | block2snp1 + block2snp2 / group, temp );  # haplotype analysis, stratified by group
####################################################################
pbat.m <- function(
       formula, phe, ped, fbat="",
       max.pheno=1, min.pheno=1,
       null="no linkage, no association", alpha=0.05,
       trans.pheno="none", trans.pred="none", trans.inter="none",
       scan.pred="all", scan.inter="all",
       scan.genetic="additive",
       offset="gee",
       screening="conditional power", distribution="default",
       logfile="",
       max.gee=1,
       max.ped=14, min.info=0,
       incl.ambhaplos=TRUE, infer.mis.snp=FALSE,
       sub.haplos=FALSE, length.haplos=2, adj.snps=TRUE,
       overall.haplo=FALSE, cutoff.haplo=FALSE,
       output="normal",
       max.mating.types=10000,
       commandfile="",
       future.expansion=NULL,
       LOAD.OUTPUT=TRUE,
       monte=0,
       mminsnps=NULL, mmaxsnps=NULL,
       mminphenos=NULL, mmaxphenos=NULL,
       env.cor.adjust=FALSE,
       gwa=FALSE,
       snppedfile=FALSE,
       extended.pedigree.snp.fix=FALSE,
       new.ped.algo=FALSE,
       cnv.intensity=2, cnv.intensity.num=3
                   )
{
  #cat("entered pbat.m") ## debug hell

  # make sure some of the variables are of a certain format
  max.pheno <- as.numeric( max.pheno );
  min.pheno <- as.numeric( min.pheno );
  alpha <- as.numeric( alpha );
  max.gee <- as.numeric( max.gee );
  min.info <- as.numeric( min.info );
  length.haplos <- as.numeric( length.haplos );
  max.mating.types <- as.numeric( max.mating.types );

  # Note that most debugging of the formula should be done by the coding
  #  that I've already written, so we really only care to verify that
  #  the parsing works correctly!

  if( !is.character(formula) ) {
    #######################################
    # first, get a string for the formula #
    #######################################

    call <- match.call( expand.dots=FALSE );
    mf <- match( c("formula"), names(call) );
    formula <- call[[mf]];

    # We actually need to reconstruct the formula :(
    # I don't want to deal with infix when the coding is so simple here!
    ## Simple solution! W
    formula <- as.character( as.expression( formula ) );
  }

  ##print( "*** FORMULA ***" ); ## debug only
  ##print( formula );
  ##print( "*** FORMULA ***" );

  ########################################
  # Eliminate all spaces in the formula. #
  ########################################
  tmpVec <- strsplitFix( formula, " " );
  formula="";
  for( i in 1:length(tmpVec) )
    formula <- paste( formula, tmpVec[i], sep="" );

  ########################
  # split up the formula #
  ########################

  tmp <- strsplitFix( as.character(formula), "~" );
  if( length(tmp) > 2 )
    stop( "The 'formula' option must be specified with the phenotypes on the lhs, the model on the rhs, sepereated by a single '~' character." );
  lhs <- tmp[1];
  rhs <- ""; if( length(tmp)==2 ) rhs <- tmp[2];

  ###############################################
  # some other things for pathological cases... #
  ###############################################
  if( lhs=="ALL" ) lhs="";
  if( rhs=="NONE" ) rhs="";

  ###########################################################################
  ## ALL doesn't function right. Removing. Doesn't make much sense anyway. ##
  ###########################################################################
  if( lhs=="" )
    stop( "You must specify some phenotype. These are from the names of the columns in the phenotype file/object, and the keyword 'AffectionStatus' which stands for the affection status of an individual. ALL is no longer supported." );


  ########################
  # take care of the lhs #
  ########################

  tmp <- strsplitFix( lhs, "&" );
  censor=""; # PASSING
  time="";   # PASSING
  phenos=""; # PASSING
  if( length(tmp) == 2 ) {
    time <- tmp[1];
    censor <- tmp[2];
    fbat="logrank";
  }else if( length(tmp) == 1 ) {
    phenos <- strsplitFix( lhs, "+" );
  }else{
    stop( "usage: pbat.m( phenos1 + phenos2 ~ ... ) or pbat.m( time & censor )" );
  }

  ############################
  # now take care of the rhs #
  ############################

  # seperate if there's any grouping
  groups.var=""; # PASSING
  groups="";     # PASSING
  tmp <- strsplitFix( rhs, "/" );
  if( length(tmp)>2 )
    stop( "Only one grouping variable can be specified" );
  if( length(tmp)==2 ) {
    rhs <- tmp[1];
    groups.var <- tmp[2];
    tmp <- strsplitFix( groups.var, " " );
    if( length(tmp)>1 ) {
      groups.var=tmp[1];
      #######groups.order <- tmp[-1];   ### 09/11/06 - there is no such thing
    }
  }

  ##########################################
  # seperate the model and block structure #
  ##########################################

  snps <- c(); # PASSING
  haplos <- list(); # PASSING
  blocks <- strsplitFix(rhs,"|");
  model <- "";
  if( length(blocks)==1 ) {
    model <- blocks;
    blocks <- c();
  }else if( length(blocks)>1 ){
    model <- blocks[1];
    blocks <- blocks[-1];

    # now check - if each of the blocks are just one word long, then we've got single snps.
    # Otherwise, we need to create the haplotype structure...
    singleSnps <- TRUE;
    for( i in 1:length(blocks) ) {
      curSnps <- strsplitFix( blocks[i], "+" ); # NEED fixed=TRUE
      if( length(curSnps) > 1 )
        singleSnps <- FALSE;
      snps <- c(snps, curSnps);
      haplos[[i]] <- curSnps;
      names(haplos)[i] <- paste("block",i,sep="");
    }

    if( singleSnps==FALSE ) {
      snps=""; # we've got a blocking structure
    }else{
      haplos=NULL;  # nope
    }
  }

  ####################################
  # now, seperate the 'model' string #
  ####################################

  # first seperate out the model by '+' signs for the preds
  preds <- strsplitFix( model, "+" ); # PASSING

  # now check the preds for 'mi(pred1)' to stand for marker interaction
  #  put it in inters, and strip the mi
  inters <- c(); # PASSING
  for( i in 1:length(preds) ) {
    if( substring(preds[i],1,3)=="mi(" ) {
      # marker interaction present!
      preds[i] <- substring( preds[i],4,strlen(preds[i])-1 );
      inters <- c(inters, preds[i]);
    }
  }

  # now check the order of each of each by spliting each pred on '^'
  preds.order <- "";
  ##if( length(preds) > 1 ) {
  if( length(preds)>=1 && preds[1]!="" ) { ## 01/18/2006 bugfix
    preds.order <- rep(1,length(preds)); # PASSING - order 1 unless changed
    for( i in 1:length(preds) ) {
      tmp <- strsplitFix( preds[i], "^" );
      if( length(tmp)==2 ) {
        preds[i] <- tmp[1];
        preds.order[i] <- tmp[2];
      }
    }
  }

  # we need to strip off the '^' in the inters as well... sigh...
  for( i in 1:length(inters) ){
    tmp <- strsplitFix(inters[i],"^");
    if( length(tmp)==2 ) inters[i] <- tmp[1];
  }

  #######################
  # fix up a few things #
  #######################
  if( is.null(snps) ) snps <- "";
  if( is.null(inters) ) inters <- "";
  if( length(preds)==1 && preds=="NONE" ) preds=""; # we actually want 'NONE' here...
  if( !is.null(haplos) && length(haplos)==0 ) haplos <- NULL;


  # Let's just include some of the debugging in this file...
  ##if( PBAT.M.DEBUG==TRUE ) {
  ##  # print out everything that has the # PASSING comment;
  ##  #  make sure we haven't forgotten anything!
  ##
  ##  print( "CENSOR" ); print( censor );
  ##  print( "TIME" ); print( time );
  ##  print( "PHENOS" ); print( phenos );
  ##  print( "GROUPS.VAR" ); print( groups.var );
  ##  print( "GROUPS" ); print( groups );
  ##  print( "SNPS" ); print( snps );
  ##  print( "HAPLOS" ); print( haplos );
  ##  print( "PREDS" ); print( preds );
  ##  print( "INTERS" ); print( inters );
  ##  print( "PREDS.ORDER" ); print( preds.order );
  ##
  ##  return(NULL);
  ##}

  if( is.null(mminsnps) ) mminsnps <- NULL; ## WHAT THE HELL??

  # And that should be it!  Just run it baby!
  res <- pbat.obj( phe, ped, paste("pbat",getTimeStamp(),"data",sep=""),
       snps=snps,
       phenos=phenos, time=time,
       preds=preds, preds.order=preds.order,
       inters=inters,
       groups.var=groups.var, groups=groups,
       fbat=fbat,
       censor=censor,
       max.pheno=max.pheno, min.pheno=min.pheno,
       null=null, alpha=alpha,
       trans.pheno=trans.pheno, trans.pred=trans.pred, trans.inter=trans.inter,
       scan.pred=scan.pred, scan.inter=scan.inter,
       scan.genetic=scan.genetic,
       offset=offset,
       screening=screening, distribution=distribution,
       logfile=logfile,
       max.gee=max.gee,
       max.ped=max.ped, min.info=min.info,
       haplos=haplos, incl.ambhaplos=incl.ambhaplos, infer.mis.snp=infer.mis.snp,
       sub.haplos=sub.haplos, length.haplos=length.haplos, adj.snps=adj.snps,
       overall.haplo=overall.haplo, cutoff.haplo=cutoff.haplo,
       output=output,
       max.mating.types=max.mating.types,
       commandfile=commandfile,
       future.expansion=future.expansion,
       LOAD.OUTPUT=LOAD.OUTPUT,
       monte=monte,
       mminsnps=mminsnps, mmaxsnps=mmaxsnps,
       mminphenos=mminphenos, mmaxphenos=mmaxphenos,
       env.cor.adjust=env.cor.adjust,
       gwa=gwa,
       snppedfile=snppedfile,
       extended.pedigree.snp.fix=extended.pedigree.snp.fix,
       new.ped.algo=new.ped.algo,
       cnv.intensity=cnv.intensity, cnv.intensity.num=cnv.intensity.num
                  );

  ##print( names(res) ); ## DEBUG
  res$call <- formula;
  ##print( names(res) ); ## DEBUG
  return( res );
}
