write.out <-
function# write directly to the screen
### internal function for sisus
(text.to.cat
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  #cat(text.to.cat);
  capture.output(expr = cat(text.to.cat), append=TRUE, file="process_info.txt");
  flush.console();
  ### internal variable
}
