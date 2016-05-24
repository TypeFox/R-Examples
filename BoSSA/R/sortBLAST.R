`sortBLAST` <-
function(blastout)
{
  def <- blastout['BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_def/#']
  accession <- blastout['BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_accession/#']
  evalue <- as.numeric(blastout['BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp/Hsp_evalue/#'][blastout['BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp/Hsp_num/#']==1])
  score <- as.numeric(blastout['BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp/Hsp_score/#'][blastout['BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp/Hsp_num/#']==1])
  if(length(def)==0) def <- "no hit"
  if(length(accession)==0) accession <- "no hit"
  if(length(evalue)==0) evalue <- "no hit"
  if(length(score)==0) score <- "no hit"
  out <- data.frame(accession=accession,definition=def,score=score,evalue=evalue,stringsAsFactors = FALSE)
  out
}

