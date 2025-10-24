liveTreeAge <- function(data = data,
                        speciesCol = "Species",
                        ageCol = "Age",
                        dbhCol = "DBH",
                        statusCol = "Status",
                        measYear = measYear){
  
  # Pull out column names as strings
  Species <- speciesCol
  Age <- ageCol
  DBH <- dbhCol
  Status <- statusCol
  
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

