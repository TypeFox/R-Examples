write.progress <-
function# write progress to the screen
### internal function for sisus
(text.to.cat
### internal variable
, time.start
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  time.sofar = progress.time(time.start);
  #cat("Progress:", time.sofar, text.to.cat);
  #capture.output(expr = cat("...", time.sofar, text.to.cat), append=TRUE, file="process_info.txt");
  #capture.output(expr = cat("...", sprintf("%8.2f",time.sofar), text.to.cat), append=TRUE, file="process_info.txt");
  capture.output(expr = cat(sprintf(">%8.2f",time.sofar), text.to.cat), append=TRUE, file="process_info.txt");
  flush.console();
  ### internal variable
}
