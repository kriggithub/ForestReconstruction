plotStandReconstruction <- function(finalData) {
  ## ---- Setup multi-panel layout ----
  # one density plot + one BA plot per species
  df <- finalData[[1]]
  densCols <- grep("\\.density$", names(df), value = TRUE)
  
  allDF <- do.call(rbind, lapply(names(finalData), function(p) {
    cbind(percentile = p, finalData[[p]])
  }))
  baCols <- grep("\\.ba$", names(allDF), value = TRUE)
  
  nPlots <- 1 + length(baCols)   # 1 density plot + BA plots
  nCols <- ceiling(sqrt(nPlots))
  nRows <- ceiling(nPlots / nCols)
  op <- par(mfrow = c(nRows, nCols), mar = c(4,4,3,1))
  on.exit(par(op))  # reset layout when function exits
  
  ## ---- DENSITY plot (all species, once) ----
  if (length(densCols) > 0) {
    plot(df$refYear, df[[densCols[1]]],
         type = "l", lwd = 2,
         ylim = range(df[densCols], na.rm = TRUE),
         xlim = rev(range(df$refYear)),
         xlab = "Reference Year", ylab = "Density (trees/ha)",
         main = "Density vs Reference Year")
    
    if (length(densCols) > 1) {
      for (j in 2:length(densCols)) {
        lines(df$refYear, df[[densCols[j]]], col = j, lwd = 2)
      }
    }
    
    legend("topright", legend = gsub("\\.density", "", densCols),
           col = seq_along(densCols), lwd = 2, bty = "n", cex = 0.8)
  }
  
  ## ---- BASAL AREA plots (all percentiles on one plot per species) ----
  for (spCol in baCols) {
    spName <- gsub("\\.ba", "", spCol)
    
    plot(allDF$refYear, allDF[[spCol]], type = "n",
         ylim = range(allDF[[spCol]], na.rm = TRUE),
         xlim = rev(range(allDF$refYear)),
         xlab = "Reference Year", ylab = "Basal Area (m²/ha)",
         main = paste("Basal Area -", spName))
    
    for (p in unique(allDF$percentile)) {
      percentiles <- unique(allDF$percentile)
      lowHigh <- c(min(percentiles), max(percentiles))
      ltyVal <- ifelse(p %in% lowHigh, 2, 1)
      lines(allDF$refYear[allDF$percentile == p],
            allDF[[spCol]][allDF$percentile == p],
            lty = ltyVal, lwd = 2)
    }
    
    legend("topright", legend = unique(allDF$percentile),
           lty = c(2,1,2), lwd = 2, bty = "n", title = "Percentile", cex = 0.8)
  }
}




plotStandReconstruction(Alpine)
plotStandReconstruction(Gumo)
