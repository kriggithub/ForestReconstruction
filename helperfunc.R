normCols <- function(x) {
  if (is.character(x)) {
    x
  } else {
    deparse(substitute(x))
  }
}

liveTreeAge <- function(data = data,
                        speciesCol = Species,
                        ageCol = Age,
                        dbhCol = DBH,
                        statusCol = Status,
                        measYear = measYear){
  
  # Pull out column names as strings
  Species <- normCols(speciesCol)
  Age <- normCols(ageCol)
  DBH <- normCols(dbhCol)
  Status <- normCols(statusCol)
  
  # Loop for each unique species in data
  for (sp in unique(data[[Species]])) {
    
    # Set LM to predict Age by DBH
    formula <- stats::as.formula(paste(Age, "~", DBH))
    
    # Create new data frame for live trees with know ages for LM
    speciesData <- data[data[[Species]] == sp &
                          !is.na(data[[Age]]) &
                          data[[Status]] == 1, ]
    
    # Fit LM for known trees (at least 2 needed per species)
    if (nrow(speciesData) > 1) {
      fit <- stats::lm(formula, data = speciesData)
      
      # Grab live trees that have missing age data
      missingRows <- which(is.na(data[[Age]]) &
                             !is.na(data[[Species]]) &
                             data[[Species]] == sp &
                             data[[Status]] == 1)
      
      # Fill in missing age data for live trees
      if (any(missingRows)) {
        predictedAge <- stats::predict(fit, newdata = data[missingRows, ])
        data[[Age]][missingRows] <- predictedAge
      }
    }
  }
  
  # Round ages up to next whole number
  data[[Age]] <- ceiling(as.numeric(data[[Age]]))
  
  # Create new column for established year
  data$estabYear <- measYear - data[[Age]]
  
  # Return data
  data
}

fuleDecomp <- function(rate, percentiles){
  (log(percentiles) - log(1))/(log(1 + rate))
}

createConclass <- function(statusVec, decayVec){
  ifelse(statusVec == 2 & decayVec %in%  1:2, 3,
         ifelse(statusVec == 2 & decayVec %in%  3:4, 4,
                ifelse(statusVec == 2 & decayVec %in%  5:6, 5,
                       ifelse(statusVec == 3, 6,
                              ifelse(statusVec == 4 & decayVec %in%  1:2, 3,
                                     ifelse(statusVec == 4 & decayVec %in%  3:4, 4,
                                            ifelse(statusVec == 4 & decayVec %in%  5:6, 5,
                                                   ifelse(statusVec == 4 & decayVec == 7, 8,
                                                          ifelse(statusVec == 6, 7,
                                                                 ifelse(statusVec == 5, 7, NA))))))))))
}


plotStandReconstruction <- function(standOutput) {
  ## ---- Extract summary data ----
  # Now the summary is stored under standOutput$summary
  finalData <- standOutput$finalOutput
  
  # First percentile table (just to get refYears & column names)
  df <- finalData[[1]]
  densCols <- grep("\\.density$", names(df), value = TRUE)
  
  # Bind all percentiles into one big DF
  allDF <- do.call(rbind, lapply(names(finalData), function(p) {
    cbind(percentile = p, finalData[[p]])
  }))
  baCols <- grep("\\.ba$", names(allDF), value = TRUE)
  
  ## ---- Setup multi-panel layout ----
  nPlots <- 1 + length(baCols)   # 1 density plot + one BA plot per species
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
