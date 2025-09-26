#### Full function
#### Assumes cleaned data:
# Species, Age, DBH, Status, Decay


forestStandReconstruction <- function(data = data, measYear, refYear, avgIncVec, speciesCol = Species, ageCol = Age, dbhCol = DBH, statusCol = Status, decayCol = Decay, percentiles = c(0.25, 0.5, 0.75)) {
  
  # column names as strings
  Species <- deparse(substitute(speciesCol))
  Age <- deparse(substitute(ageCol))
  DBH <- deparse(substitute(dbhCol))
  
  # create frame for final data
  finalData <- list()
  
  
  for (rY in refYear){
    
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
    
    data <- ageData
    
    # growth to subtract for species (cm per year)
    growthSub <- (avgIncVec*2)/10
    
    
    # filter trees older than refYear
    data <- data[data$estabYear <= rY & complete.cases(data), ]
    
    
    # match increments to species
    incs <- growthSub[match(data[[Species]], names(growthSub))]
    
    
    # create new columns for RefAge and RefDBH
    data$RefAge <- rY - data$estabYear
    data$RefDBH <- data[[DBH]] - ((measYear - rY) * incs)
    
    
    # keep non-negative DBH
    data <- data[data$RefDBH >= 0, ]
    
    
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
    
    Conclass <- "Conclass"
    
    
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
      stitchedData[[pRefDBH]] <- stitchedData[[DBH]] - (stitchedData[[deathCol]] - rY) * incs
    }
    
    
    data <- stitchedData
    
    percentnames <- paste0("p", percentiles*100)
    
    for (n in percentnames){
      
      pRefDBH <- paste0(n, "refDBH")
      
      posRefDBHdat <- data[data[[pRefDBH]] > 0 | is.na(data[[pRefDBH]]), ]
      
      # bark correction
      posRefDBHdat$RefDBH[posRefDBHdat$Conclass %in% 3:4] <- posRefDBHdat[[pRefDBH]][posRefDBHdat$Conclass %in% 3:4]
      posRefDBHdat$RefDBH[posRefDBHdat$Conclass %in% 5:7] <- (posRefDBHdat[[pRefDBH]][posRefDBHdat$Conclass %in% 5:7]*1.0508) + 0.2824
      
      
      
      # basal area (m^2)
      posRefDBHdat$BA <- 0.00007854 * (posRefDBHdat$RefDBH^2)
      
      
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
      
      # add refYear column
      summaryDat$refYear <- rY
      
      
      
      if (is.null(finalData[[n]])) {
        finalData[[n]] <- summaryDat
      } else {
        finalData[[n]] <- rbind(finalData[[n]], summaryDat)
      }
      
      
    }
    
  }
  return(finalData)
  
}


GumoTreeData <- read.csv("GumoFullTreeData.csv")


Alpine <- forestStandReconstruction(AlpineTreeData, 
                                    measYear = 2013,
                                    refYear = c(1960),
                                    avgIncVec = c(PIEN = 0.85, ABBI = 0.822818),
                                    speciesCol = Species,
                                    ageCol = Age,
                                    dbhCol = DBH,
                                    statusCol = Status,
                                    decayCol = Decay)  

Alpine




Alpine1975 <- forestStandReconstruction(AlpineTreeData, 
                                        measYear = 2013,
                                        refYear = 1975,
                                        avgIncVec = c(PIEN = 0.85, ABBI = 0.822818),
                                        speciesCol = Species,
                                        ageCol = Age,
                                        dbhCol = DBH,
                                        statusCol = Status,
                                        decayCol = Decay)  

Alpine1975$p25$summary
Alpine1975$p50$summary
Alpine1975$p75$summary

# area of plots, so if plot is 20 x 20, you have to multiply by 25 for hectare
# regression should only run for live trees
# presettlement trees (all alive)




Alpine1960 <- forestStandReconstruction(AlpineTreeData, 
                                        measYear = 2013,
                                        refYear = 1960,
                                        avgIncVec = c(PIEN = 0.85, ABBI = 0.822818),
                                        speciesCol = Species,
                                        ageCol = Age,
                                        dbhCol = DBH,
                                        statusCol = Status,
                                        decayCol = Decay)  

Alpine1960$p25$summary
Alpine1960$p50$summary
Alpine1960$p75$summary


data <- Alpine1960$p25$data


Gumo1922 <- forestStandReconstruction(GumoTreeData, 
                                      measYear = 2004,
                                      refYear = 1922,
                                      avgIncVec = c(PIED = 0.102439024390244/5, PIPO = 0.207317073170732/5, PIST = 0.136585366/5, PSME = 0.190243902439024/5, QUGA = 0.075609756097561/5),
                                      speciesCol = sp,
                                      ageCol = a,
                                      dbhCol = d,
                                      statusCol = st,
                                      decayCol = dc)  



Gumo1922$p25$summary
Gumo1922$p50$summary
Gumo1922$p75$summary



Gumo1975 <- forestStandReconstruction(GumoTreeData, 
                                      measYear = 2004,
                                      refYear = 1975,
                                      avgIncVec = c(PIED = 0.102439024390244/5, PIPO = 0.207317073170732/5, PIST = 0.136585366/5, PSME = 0.190243902439024/5, QUGA = 0.075609756097561/5),
                                      speciesCol = sp,
                                      ageCol = a,
                                      dbhCol = d,
                                      statusCol = st,
                                      decayCol = dc)  



Gumo1975$p25$summary
Gumo1975$p50$summary
Gumo1975$p75$summary






