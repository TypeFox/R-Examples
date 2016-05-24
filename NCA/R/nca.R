nca <-
function (data, nx, ny, 
          ols=TRUE, cols=FALSE, qr=FALSE, lh=FALSE, 
          ce_vrs=FALSE, cr_vrs=FALSE, ce_fdh=TRUE, cr_fdh=TRUE, sfa=FALSE,
          fast.fdh=TRUE, fast.vrs=TRUE, invisible=FALSE,
          title="NCA Plot", use.title=TRUE, pdf=FALSE, prefix="out",
          results=FALSE, bottleneck.x='percentage.range',
          bottleneck.y='percentage.range', steps=10, cutoff=0) {

  bottleneck.x <- p_validate_bottleneck(bottleneck.x, "x")
  bottleneck.y <- p_validate_bottleneck(bottleneck.y, "y")

  all_loopdata    <- matrix(list(), nrow = nx, ncol = ny)
  all_values      <- matrix(list(), nrow = nx, ncol = ny)
  all_bottlenecks <- matrix(list(), nrow = 1, ncol = ny)

  # Start looping the Ys
  for (id.y in 1:ny) {
    bottleneck.xy <- p_mp_mpy(data, id.y, nx, steps, bottleneck.y)
    mpx <- bottleneck.xy[[1]]
    mpy <- bottleneck.xy[[2]]

    bottlenecks <- list()
    if (lh)     { bottlenecks[['lh']]     <- mpx }
    if (cols)   { bottlenecks[['cols']]   <- mpx }
    if (qr)     { bottlenecks[['qr']]     <- mpx }
    if (cr_vrs) { bottlenecks[['cr_vrs']] <- mpx }
    if (ce_vrs) { bottlenecks[['ce_vrs']] <- mpx }
    if (ce_fdh) { bottlenecks[['ce_fdh']] <- mpx }
    if (cr_fdh) { bottlenecks[['cr_fdh']] <- mpx }
    if (sfa)    { bottlenecks[['sfa']]    <- mpx }

    # Start looping the Xs
    for (id.x in 1:nx) {
      loop.data <- loop_data(data, id.x, id.y, nx)
      all_loopdata[[id.x, id.y]] <- loop.data

      p_warn_percentage_max(bottleneck.y, loop.data, id.y)

      results.nca <- list()
      if (ols) {
        analysis <- nca_ols(loop.data, mpy, cutoff, bottleneck.x)
        analysis$bottleneck <- NULL
        results.nca$ols <- analysis
      }
      if (lh) {
        analysis <- nca_lh(loop.data, mpy, cutoff, bottleneck.x)
        bottlenecks$lh <- cbind(bottlenecks$lh, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$lh <- analysis
      }
      if (cols) {
        analysis <- nca_cols(loop.data, mpy, cutoff, bottleneck.x)
        bottlenecks$cols <- cbind(bottlenecks$cols, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$cols <- analysis
      }
      if (qr) {
        analysis <- nca_qr(loop.data, mpy, cutoff, bottleneck.x)
        bottlenecks$qr <- cbind(bottlenecks$qr, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$qr <- analysis
      }
      if (ce_vrs) {
        analysis <- nca_ce_vrs(loop.data, mpy, cutoff, bottleneck.x, fast.vrs)
        bottlenecks$ce_vrs <- cbind(bottlenecks$ce_vrs, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$ce_vrs <- analysis
      }
      if (cr_vrs) {
        analysis <- nca_cr_vrs(loop.data, mpy, cutoff, bottleneck.x, fast.vrs)
        bottlenecks$cr_vrs <- cbind(bottlenecks$cr_vrs, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$cr_vrs <- analysis
      }
      if (ce_fdh) {
        analysis <- nca_ce_fdh(loop.data, mpy, cutoff, bottleneck.x, fast.fdh)
        bottlenecks$ce_fdh <- cbind(bottlenecks$ce_fdh, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$ce_fdh <- analysis
      }
      if (cr_fdh) {
        analysis <- nca_cr_fdh(loop.data, mpy, cutoff, bottleneck.x, fast.fdh)
        bottlenecks$cr_fdh <- cbind(bottlenecks$cr_fdh, analysis$bottleneck)
        analysis$bottleneck <- NULL
        results.nca$cr_fdh <- analysis
      }
      if (sfa) {
        analysis <- nca_sfa(loop.data, mpy, cutoff, bottleneck.x)
        if (!is.na(analysis$bottleneck)) {
          bottlenecks$sfa <- cbind(bottlenecks$sfa, analysis$bottleneck)
        }
        analysis$bottleneck <- NULL
        results.nca$sfa <- analysis
      }

      # Output graph and tables for this X+Y
      if (!invisible) {
        p_display_graphs(results.nca, loop.data, title, use.title, pdf, prefix)
        p_display_tables(results.nca, loop.data, use.title, pdf)
        if (results) {
          p_display_results(results.nca, loop.data, prefix, pdf)
        }
      }

      all_values[[id.x, id.y]] <- results.nca
    }

    # Display bottleneck for this Y
    if (!invisible) {
      p_display_bottleneck(bottlenecks, colnames(data), id.y, prefix, use.title, pdf, bottleneck.x, bottleneck.y, nx, steps)
    }

    all_bottlenecks[[1, id.y]] <- bottlenecks
  }

  if (invisible) {
    return ( list(
      loopdatas=all_loopdata, values=all_values, bottlenecks=all_bottlenecks ) )
  }
}