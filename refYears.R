#### Full function
#### Assumes cleaned data:
# Species, Age, DBH, Status, Decay


forestStandReconstruction <- function(data = data, measYear, refYear, avgIncVec, plotSize, barkEqList = list(), speciesCol = Species, ageCol = Age, dbhCol = DBH, statusCol = Status, decayCol = Decay, percentiles = c(0.25, 0.5, 0.75)) {
  
  # pull out column names as strings
  Species <- deparse(substitute(speciesCol))
  Age <- deparse(substitute(ageCol))
  DBH <- deparse(substitute(dbhCol))
  Status <- deparse(substitute(statusCol))
  Decay <- deparse(substitute(decayCol))
  
  # growth to subtract for species (cm per year)
  growthSub <- (avgIncVec*2)/10
  
  # Fule decomposition rate function
  decompRate <- function(rate, percentiles){
    (log(percentiles) - log(1))/(log(1 + rate))
  }
  
  # create names for percents
  percentnames <- paste0("p", percentiles*100)
  
  # create final dataframe
  finalData <- list()
  
  
  # lookup table for default bark equations
  defaultBarkEq <- list(
    PIEN = function(x) x*1.0508 + 0.2824,
    ABBI = function(x) x*1.0508 + 0.2824,
    ABCO = function(x) x*1.1238 + 0.1952,
    ABMA = function(x) ((x*1.1799)^0.9803) + 0.4403,
    CADE = function(x) x*1.1975 - 0.0427,
    PILA = function(x) x*1.1252 + 0.2488,
    PIST = function(x) x*1.1252 + 0.2488,
    PIED = function(x) x*1.1252 + 0.2488,
    PIPO = function(x) x*1.1029 + 0.7162,
    PIJE = function(x) x*1.1029 + 0.7162,
    PSME = function(x) x*1.1759 - 0.2721
  )
  
  
  barkEqList <- modifyList(defaultBarkEq, barkEqList)
  
  
  for (sp in names(barkEqList)) {
    if (is.character(barkEqList[[sp]]) && barkEqList[[sp]] %in% names(barkEqList)) {
      barkEqList[[sp]] <- barkEqList[[ barkEqList[[sp]] ]]
    }
  }
  
  
  # Step 1: fit LMs per species
  for (sp in unique(data[[Species]])) {
    formula <- as.formula(paste(Age, "~", DBH))
    
    # dataframe for lm (has the target species with Age not NA)
    speciesData <- data[data[[Species]] == sp & 
                          !is.na(data[[Age]]) &
                          data[[Status]] == 1, ]
    
    if (nrow(speciesData) > 1) {
      fit <- lm(formula, data = speciesData)
      
      missingRows <- is.na(data[[Age]]) & data[[Species]] == sp
      
      if (any(missingRows)) {
        predictedAge <- predict(fit, newdata = data[missingRows, ])
        data[[Age]][missingRows] <- predictedAge
      }
    } 
  }
  
  data[[Age]] <- ceiling(as.numeric(data[[Age]]))
  data$estabYear <- measYear - data[[Age]]
  
  
  
  # filter trees older than refYear
  for (rY in refYear){
  
  dat <- data
  
  dat <- dat[dat$estabYear <= rY & complete.cases(dat), ]
  
  
  # match increments to species
  # incs <- growthSub[match(dat[[Species]], names(growthSub))]
  incs <- growthSub[dat[[Species]]]
  
  
  # create new columns for RefAge and RefDBH
  dat$RefAge <- rY - dat$estabYear
  dat$RefDBH <- dat[[DBH]] - ((measYear - rY) * incs)
  
  
  # keep non-negative DBH
  dat <- dat[dat$RefDBH >= 0, ]
  
  
  # pull out vectors
  statusVec <- dat[[Status]]
  decayVec <- dat[[Decay]]
  
  
  # create vector with conclasses
  conclassVec <- ifelse(statusVec == 2 & decayVec %in%  1:2, 3, 
                        ifelse(statusVec == 2 & decayVec %in%  3:4, 4, 
                               ifelse(statusVec == 2 & decayVec %in%  5:6, 5, 
                                      ifelse(statusVec == 3, 6, 
                                             ifelse(statusVec == 4 & decayVec %in%  1:2, 3,
                                                    ifelse(statusVec == 4 & decayVec %in%  3:4, 4, 
                                                           ifelse(statusVec == 4 & decayVec %in%  5:6, 5, 
                                                                  ifelse(statusVec == 4 & decayVec == 7, 8,
                                                                         ifelse(statusVec == 6, 7, 
                                                                                ifelse(statusVec == 5, 7, NA))))))))))
  
  
  # make vector into dat column
  dat$Conclass <- conclassVec
  
  
  # average DBH by species alive at reference year (cm)
  avgDBHcm <- tapply(
    dat[[DBH]][data[[Status]] == 1],
    dat[[Species]][data[[Status]] == 1],
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
  for (n in percentnames){
    decompreference[[as.character(n)]] <- NA
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
      
      # write back to main datset
      decompreference[spRows, cname] <- vals
    }
  }
  
  
  
  dat <- merge(dat, decompreference, by.x = c(Species, "Conclass"), by.y = c("Species", "Conclass"), all.x = T)
  
  
  for (n in percentnames){
    deathCol <- paste0(n, "DeathYear")
    pRefDBH <- paste0(n, "refDBH")
    
    # create death year columns
    dat[[deathCol]] <- ceiling(measYear + dat[[n]])
    
    
    # create reference DBH columns
    incs <- growthSub[dat[[Species]]]
    dat[[pRefDBH]] <- dat[[DBH]] - (dat[[deathCol]] - rY) * incs
    
    pRefDBH <- paste0(n, "refDBH")
    
    dat <- dat[dat[[pRefDBH]] > 0 | is.na(dat[[pRefDBH]]), ]
    
    # bark correction
    dat$RefDBH[dat$Conclass %in% 3:4] <- dat[[pRefDBH]][dat$Conclass %in% 3:4]
    
    mask <- dat$Conclass %in% 5:7
    for (sp in unique(dat[[Species]][mask])) {
      spMask <- mask & dat[[Species]] == sp
      if (sp %in% names(barkEqList)) {
        dat$RefDBH[spMask] <- barkEqList[[sp]](dat[[pRefDBH]][spMask])
      } else {
        # if really no match, just leave unchanged
        dat$RefDBH[spMask] <- dat[[pRefDBH]][spMask]
      }
    }
    
    
    # dat$RefDBH[dat$Conclass %in% 5:7] <- (dat[[pRefDBH]][dat$Conclass %in% 5:7]*1.0508) + 0.2824
    
    
    
    # basal area (m^2)
    dat$BA <- 0.00007854 * (dat$RefDBH^2)
    
    
    
    BAsum <- rowsum(dat$BA, dat[[Species]], na.rm = TRUE)
    # BAsum <- dat.frame(Species = rownames(BAsum), BA = BAsum[,1], row.names = NULL)
    BAw <- setNames(as.list(BAsum[,1]), paste0(rownames(BAsum), ".ba"))
    
    
    
    
    spCounts <- table(dat[[Species]])
    Densw <- setNames(as.list(as.numeric(spCounts) * (10000/plotSize)),
                      paste0(names(spCounts), ".density"))
    

    rowOut <- data.frame(
      refYear = rY,
      BAw,
      Densw,
      check.names = FALSE
    )
    
    # append row to that percentile’s table
    if (is.null(finalData[[n]])) {
      finalData[[n]] <- rowOut
    } else {
      finalData[[n]] <- rbind(finalData[[n]], rowOut)
    }
    
    
  }
  
  
}

  
  return(finalData)
  
}



#### Function usage
AlpineTreeData <- read.csv("AlpineFullTreeData.csv")
GumoTreeData <- read.csv("GumoFullTreeData.csv")



Alpine <- forestStandReconstruction(AlpineTreeData, 
                                        measYear = 2013,
                                        refYear = seq(from = 1900, to = 2013, by = 1),
                                        avgIncVec = c(PIEN = 0.85, ABBI = 0.822818),
                                        plotSize = 400,
                                        speciesCol = Species,
                                        ageCol = Age,
                                        dbhCol = DBH,
                                        statusCol = Status,
                                        decayCol = Decay 
                                        )  

Alpine$p50
# plot(Alpine$p50$refYear, Alpine$p50$ABBI.density)
# plot(Alpine$p50$refYear, Alpine$p50$PIEN.density)
# plot(Alpine$p50$refYear, Alpine$p50$ABBI.ba)
# plot(Alpine$p50$refYear, Alpine$p50$PIEN.ba)




Gumo <- forestStandReconstruction(GumoTreeData, 
                                      measYear = 2004,
                                      refYear = seq(from = 1890, to = 2004, by = 1),
                                      avgIncVec = c(PIED = 0.102439024390244/5, PIPO = 0.207317073170732/5, PIST = 0.136585366/5, PSME = 0.190243902439024/5, QUGA = 0.075609756097561/5),
                                      plotSize = 400,
                                      speciesCol = sp,
                                      ageCol = a,
                                      dbhCol = d,
                                      statusCol = st,
                                      decayCol = dc)  



Gumo$p50
# plot(Gumo$p50$refYear, Gumo$p50$PIPO.density)


