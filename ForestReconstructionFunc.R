# Function for Forest Stand Reconstruction
# Assumes Cleaned Dataset (fullTreeData)


AlpineTreeData <- read.csv("AlpineFullTreeData.csv")

# Function to run linear regression on aged trees

ageLM <- function(data, speciesCol = Species, ageCol = Age, dbhCol = DBH, measYear) {
  
  # column names as strings
  Species <- deparse(substitute(speciesCol))
  Age <- deparse(substitute(ageCol))
  DBH <- deparse(substitute(dbhCol))
  
  # new dataframe for predicted ages
  ageData <- data
  
  for (sp in unique(data[[Species]])) {
    formula <- as.formula(paste(Age, "~", DBH))
    
    # dataframe for lm (has the target species with Age not NA)
    speciesData <- ageData[ageData[[Species]] == sp & !is.na(ageData[[Age]]), ]
    
    if (nrow(speciesData) > 1) {
      fit <- lm(formula, data = speciesData)
      
      missingRows <- is.na(ageData[[Age]]) & ageData[[Species]] == sp
      
      if (any(missingRows)) {
        predictedAge <- predict(fit, newdata = ageData[missingRows, ])
        ageData[[Age]][missingRows] <- predictedAge
      }
    } 
  }
  
  ageData[[Age]] <- ceiling(as.numeric(ageData[[Age]]))
  ageData$estabYear <- measYear - ageData[[Age]]
  
  
  return(ageData)
}



a <- ageLM(data = AlpineTreeData, measYear = 2013)
# b <- ageLM(data = GumoTreeData, species, age, dbh, measYear = 2004)




# function for subtracting DBH back to reference year
# avgIncVec: a named numeric vector, names must match species codes in data
# Ex: c(ABBI = 0.822818, PIEN = 0.85)

refLiveDBH <- function(data, refYear, measYear, speciesCol, dbhCol, avgIncVec) {
 
  Species <- deparse(substitute(speciesCol))
  DBH <- deparse(substitute(dbhCol))
  
  
  # growth to subtract for species (cm per year)
  growthSub <- (avgIncVec*2)/10
  
  
  # filter trees older than refYear
  data <- data[data$estabYear <= refYear & complete.cases(data), ]

  
  # match increments to species
  incs <- growthSub[match(data[[Species]], names(growthSub))]
  
  
  # create new columns for RefAge and RefDBH
  data$RefAge <- refYear - data$estabYear
  data$RefDBH <- data[DBH] - ((measYear - refYear) * incs)
  
  
  # keep non-negative DBH
  data <- data[data$RefDBH >= 0, ]
  
  return(data)
   
}

a2 <- refLiveDBH(a, 1960, 2013, speciesCol = Species, dbhCol = DBH, avgIncVec = c(PIEN = 0.85, ABBI = 0.822818))




# function to create conclass column

conclass <- function(data, statusCol, decayCol){
  
  # grab column names
  statusName <- deparse(substitute(statusCol))
  decayName <- deparse(substitute(decayCol))
  
  # pull out vectors
  status <- data[[statusName]]
  decay <- data[[decayName]]
  
  
  # create vector with conclasses
  conclassVec <- ifelse(status == 2 & decay %in%  1:2, 3, 
  ifelse(status == 2 & decay %in%  3:4, 4, 
  ifelse(status == 2 & decay %in%  5:6, 5, 
  ifelse(status == 3, 6, 
  ifelse(status == 4 & decay %in%  1:2, 3,
  ifelse(status == 4 & decay %in%  3:4, 4, 
  ifelse(status == 4 & decay %in%  5:6, 5, 
  ifelse(status == 4 & decay == 7, 8,
  ifelse(status == 6, 7, 
  ifelse(status == 5, 7, NA))))))))))
  
  
  # make vector into data column
  data$Conclass <- conclassVec
  
  return(data)
  
}


a3 <- conclass(a2, Status, Decay)





# function for reviving dead trees


decompRate <- function(data, speciesCol = Species, dbhCol = DBH, conclassCol = Conclass, percentiles = c(0.25, 0.5, 0.75),  avgIncVec, measYear, refYear){
  
  # column names as strings
  Species <- deparse(substitute(speciesCol))
  DBH <- deparse(substitute(dbhCol))
  Conclass <- deparse(substitute(conclassCol))
  
  
  # average DBH by species alive at reference year (cm)
  avgDBHcm <- tapply(
    data[[DBH]],
    data[[Species]],
    mean,
    na.rm = T
  )
  
  # convert to inches
  avgDBHin <- avgDBHcm/2.54
  
  
  snagLife <- 2 * avgDBHin
  snagFall <- 1 / snagLife
  
  
  
  
  # create reference table
  speciesnames <- names(avgDBHcm)
  conclasses <- 3:7
  decompreference <-  expand.grid(Species = speciesnames, 
                                  Conclass = conclasses)
  decompreference$Rate <- NA_real_
  
  
  # Fill rate by conclass
  decompreference$Rate[decompreference$Conclass == 3] <- 0
  decompreference$Rate[decompreference$Conclass == 4] <- 0.2
  decompreference$Rate[decompreference$Conclass == 5] <- 0.15
  decompreference$Rate[decompreference$Conclass == 6] <- snagFall
  decompreference$Rate[decompreference$Conclass == 7] <- snagFall
  
  
  
  # Create empty percentile columns
  percentnames <- paste0("p", percentiles*100)
  for (n in percentnames){
    decompreference[[as.character(n)]] <- NA
  }
  
  # Fule decomposition rate
  decompRate <- function(rate, percentiles){
    (log(percentiles) - log(1))/(log(1 + rate))
  }
  
  
  # loop through species
  
  for (sp in unique(decompreference$Species)) {
    spRows <- decompreference$Species == sp
    decompSubDF <- decompreference[spRows, ]
    decompSubDF <- decompSubDF[order((decompSubDF$Conclass))]
    
    for (i in seq_along(percentiles)) {
      p <- percentiles[i]
      cname <- percentnames[i]
      
      vals <- rep(NA_real_, nrow(decompSubDF))
      
      # Conclass == 3 always 0
      vals[decompSubDF$Conclass == 3] <- 0
      
      # Conclass >= 4: cumulative sum of decompRate
      concl4up <- which(decompSubDF$Conclass >= 4)
      if (length(concl4up) > 0) {
        vals[concl4up] <- cumsum(decompRate(decompSubDF$Rate[concl4up], p))
      }
      
      # write back to main dataset
      decompreference[spRows, cname] <- vals
    }
  }
  
  
  
  stitchedData <- merge(data, decompreference, by.x = c(Species, Conclass), by.y = c("Species", "Conclass"), all.x = T)
  
  
  for (n in percentnames){
    deathCol <- paste0(n, "DeathYear")
    pRefDBH <- paste0(n, "refDBH")

    # create death year columns
    stitchedData[[deathCol]] <- ceiling(measYear + stitchedData[[n]])


    # create reference DBH columns
    growthSub <- (avgIncVec*2)/10
    incs <- growthSub[stitchedData[[Species]]]
    stitchedData[[pRefDBH]] <- stitchedData[[DBH]] - (stitchedData[[deathCol]] - refYear) * incs
  }
  
  
  
  # 
  
  
  return(stitchedData)

  

}


a4 <- decompRate(a3, avgIncVec = c(ABBI = 0.822818, PIEN = 0.85), measYear = 2013, refYear = 1960)




# correct bark for dead trees


addBark <- function(data, conclassCol = Conclass, speciesCol, percentiles = c(0.25, 0.5, 0.75)) {
  
  Species <- deparse(substitute(speciesCol))
  Conclass <- deparse(substitute(conclassCol))
  
  percentnames <- paste0("p", percentiles*100)
  
  finalData <- list()
  
  for (n in percentnames){
    
    pRefDBH <- paste0(n, "refDBH")
    
    posRefDBHdat <- data[data[[pRefDBH]] > 0 | is.na(data[[pRefDBH]]), ]
    
    # bark correction
    posRefDBHdat$RefDBH <- NA_real_
    posRefDBHdat$RefDBH[posRefDBHdat$Conclass %in% 3:4] <- posRefDBHdat[[pRefDBH]][posRefDBHdat$Conclass %in% 3:4]
    posRefDBHdat$RefDBH[posRefDBHdat$Conclass %in% 5:7] <- (posRefDBHdat[[pRefDBH]][posRefDBHdat$Conclass %in% 5:7]*1.0508) + 0.2824
    
    
    
    # basal area (m^2)
    posRefDBHdat$BA <- 0.00007854 * (posRefDBHdat$RefDBH^2)
    
    
    # summarized data by species
    # BAsum <- aggregate(x = posRefDBHdat$BA,
    #                    by = list(Species = posRefDBHdat[[Species]]),
    #                    FUN = sum,
    #                    na.rm = TRUE)
    # names(BAsum)[2] <- "BA"
    
    
    BAsum <- rowsum(posRefDBHdat$BA, posRefDBHdat[[Species]], na.rm = TRUE)
    BAsum <- data.frame(Species = rownames(BAsum), BA = BAsum[,1], row.names = NULL)
    spCounts <- table(posRefDBHdat[[Species]])
    
    
    summaryDat <- merge(
      BAsum,
      data.frame(Species = names(spCounts),
                 treeDensity = as.numeric(spCounts) * 25),
      by = "Species",
      all = T
    )
    
    
    
    finalData[[n]] <- list(
      data = posRefDBHdat,
      summary = summaryDat
    )
      
      

    
  }
  
  

  
  return(finalData)
  
  
}


a5 <- addBark(a4, Conclass, Species)

a5$p25$summary
a5$p50$summary
a5$p75$summary




