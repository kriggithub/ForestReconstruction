corwinaCppData <- read.csv("corwinaCppData.csv")
corwinaC2Data <- read.csv("corwinaC2Data.csv")
GumoTreeData <- read.csv("GumoFullTreeData.csv")
AlpineTreeData <- read.csv("AlpineFullTreeData.csv")



Alpine <- standRecon(AlpineTreeData, 
                     measYear = 2013,
                     refYear = 1900,
                     avgIncVec = c(PIEN = 0.85, ABBI = 0.822818),
                     plotSize = 400,
                     nPlots = 23,
                     speciesCol = "Species",
                     ageCol = "Age",
                     dbhCol = "DBH",
                     statusCol = "Status",
                     decayCol = "Decay")  




corwinaCpp <- standRecon(corwinaCppData,
                         measYear = 2023,
                         refYear = 1912,
                         avgIncVec = c(PSME = 1.0282, PIPO = 0.6067),
                         plotSize = 400,
                         nPlots = 5,
                         speciesCol = "Species",
                         ageCol = "Adjusted.Age",
                         dbhCol = "DBH",
                         statusCol = "Status",
                         decayCol = "Decay")




corwinaC2 <- standRecon(corwinaC2Data,
                        measYear = 2023,
                        refYear = 1912,
                        avgIncVec = c(PSME = 0.8976, PIPO = 0.512),
                        plotSize = 400,
                        nPlots = 1,
                        speciesCol = "Species",
                        ageCol = "Adjusted.Age",
                        dbhCol = "DBH",
                        statusCol = "Status",
                        decayCol = "Decay")

Gumo <- standRecon(GumoTreeData,
                   measYear = 2004,
                   refYear = 1922,
                   avgIncVec = c(PIED = 0.102439024390244*5, 
                                 PIPO = 0.207317073170732*5, 
                                 PIST = 0.136585366*5, 
                                 PSME = 0.190243902439024*5, 
                                 QUGA = 0.075609756097561*5),
                   plotSize = 400,
                   nPlots = 159,
                   speciesCol = "sp",
                   ageCol = "a",
                   dbhCol = "d",
                   statusCol = "st",
                   decayCol = "dc")




# Extract outputs
finalData5 <- Alpine$finalOutput


# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData5)) {
  names(finalData5$measured) <- gsub("\\.measYear$", "", names(finalData5$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData5), function(p) {
  df <- finalData5[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2013",
  p25 = "1900 (25%)",
  p50 = "1900 (50%)",
  p75 = "1900 (75%)"
)

# choose density columns
cols <- grep("\\.density$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Density (trees/ha)",
        main = paste("Alpine Tree Density"))

legend("topright",
       legend = gsub("\\.density$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)






# Extract outputs
finalData <- corwinaCpp$finalOutput


# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData)) {
  names(finalData$measured) <- gsub("\\.measYear$", "", names(finalData$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData), function(p) {
  df <- finalData[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2023",
  p25 = "1912 (25%)",
  p50 = "1912 (50%)",
  p75 = "1912 (75%)"
)

# choose density columns
cols <- grep("\\.density$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Density (trees/ha)",
        main = paste("Corwina CPP Tree Density"))

legend("topright",
       legend = gsub("\\.density$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)











# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData)) {
  names(finalData$measured) <- gsub("\\.measYear$", "", names(finalData$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData), function(p) {
  df <- finalData[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2023",
  p25 = "1912 (25%)",
  p50 = "1912 (50%)",
  p75 = "1912 (75%)"
)

# choose basal area columns
cols <- grep("\\.ba$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Basal Area (m^2/ha)",
        main = paste("Corwina CPP Basal Area"))

legend("topright",
       legend = gsub("\\.ba$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)















# Extract outputs
finalData2 <- corwinaC2$finalOutput


# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData2)) {
  names(finalData2$measured) <- gsub("\\.measYear$", "", names(finalData2$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData2), function(p) {
  df <- finalData2[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2023",
  p25 = "1912 (25%)",
  p50 = "1912 (50%)",
  p75 = "1912 (75%)"
)

# choose density columns
cols <- grep("\\.density$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Density (trees/ha)",
        main = paste("Corwina C2 Tree Density"))

legend("topright",
       legend = gsub("\\.density$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)











# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData2)) {
  names(finalData2$measured) <- gsub("\\.measYear$", "", names(finalData2$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData2), function(p) {
  df <- finalData2[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2023",
  p25 = "1912 (25%)",
  p50 = "1912 (50%)",
  p75 = "1912 (75%)"
)

# choose basal area columns
cols <- grep("\\.ba$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Basal Area (m^2/ha)",
        main = paste("Corwina C2 Basal Area"))

legend("topright",
       legend = gsub("\\.ba$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)

















# Extract outputs
finalData3 <- Gumo$finalOutput


# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData3)) {
  names(finalData3$measured) <- gsub("\\.measYear$", "", names(finalData3$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData3), function(p) {
  df <- finalData3[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2004",
  p25 = "1922 (25%)",
  p50 = "1922 (50%)",
  p75 = "1922 (75%)"
)

# choose density columns
cols <- grep("\\.density$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Density (trees/ha)",
        main = paste("GUMO Tree Density"))

legend("topright",
       legend = gsub("\\.density$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)











# --- make sure measured columns match percentile columns ---
if ("measured" %in% names(finalData3)) {
  names(finalData3$measured) <- gsub("\\.measYear$", "", names(finalData3$measured))
}

# Combine into one dataframe
allDF <- do.call(rbind, lapply(names(finalData3), function(p) {
  df <- finalData3[[p]]
  df$scenario <- p
  df
}))

# Explicit order of scenarios
allDF$scenario <- factor(allDF$scenario, 
                         levels = c("measured", "p25", "p50", "p75"))

# ---- Rename scenarios for display ----
scenario_labels <- c(
  measured = "2004",
  p25 = "1922 (25%)",
  p50 = "1922 (50%)",
  p75 = "1922 (75%)"
)

# choose basal area columns
cols <- grep("\\.ba$", names(allDF), value = TRUE)

# Build matrix and reorder rows
mat <- as.matrix(allDF[, cols])
rownames(mat) <- allDF$scenario
mat <- mat[c("measured", "p25", "p50", "p75"), , drop = FALSE]  # enforce order

# Make labels in that same order
plot_labels <- scenario_labels[rownames(mat)]

# Plot
barplot(t(mat),
        beside = FALSE,
        col = gray.colors(ncol(mat), start = 0.9, end = 0.3),
        border = NA,
        names.arg = plot_labels,
        ylab = "Basal Area (m^2/ha)",
        main = paste("GUMO Basal Area"))

legend("topright",
       legend = gsub("\\.ba$", "", cols),
       fill = gray.colors(ncol(mat), start = 0.9, end = 0.3), bty = "n", cex = 0.8)


