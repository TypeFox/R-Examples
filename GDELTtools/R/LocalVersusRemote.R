LocalVersusRemote <- function(filelist, local.folder) {
  # Given a list of files, determines which ones are local and which are remote
  return( list(local=filelist[filelist %in% dir(local.folder)],
               remote=filelist[!(filelist %in% dir(local.folder))]) )
}
