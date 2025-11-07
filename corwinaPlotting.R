library(tidyverse)
corwinaCppData <- read.csv("corwinaCppData.csv")
corwinaC2Data <- read.csv("corwinaC2Data.csv")

# Stand Reconstruction
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
                        plotSize = 1000,
                        nPlots = 1,
                        speciesCol = "Species",
                        ageCol = "Adjusted.Age",
                        dbhCol = "DBH",
                        statusCol = "Status",
                        decayCol = "Decay")

# Reformat Final Output
corwinaCpp <- bind_rows(
  corwinaCpp$finalOutput$p25 %>% mutate(Sensitivity = "1912 (25%)"),
  corwinaCpp$finalOutput$p50 %>% mutate(Sensitivity = "1912 (50%)"),
  corwinaCpp$finalOutput$p75 %>% mutate(Sensitivity = "1912 (75%)"),
  corwinaCpp$finalOutput$measured %>% mutate(Sensitivity = "2023")
)

corwinaC2 <- bind_rows(
  corwinaC2$finalOutput$p25 %>% mutate(Sensitivity = "1912 (25%)"),
  corwinaC2$finalOutput$p50 %>% mutate(Sensitivity = "1912 (50%)"),
  corwinaC2$finalOutput$p75 %>% mutate(Sensitivity = "1912 (75%)"),
  corwinaC2$finalOutput$measured %>% mutate(Sensitivity = "2023")
)



# Basal Area Plots
# Corwina CPP

corwinaCppBA <- corwinaCpp %>% 
  select(Sensitivity, refYear, PIPO.ba, PSME.ba) %>% 
  pivot_longer(cols = c(PIPO.ba, PSME.ba),
               names_to = "Species",
               values_to = "BasalArea")

corwinaCppBA$Species <- sub("\\.ba$", "", corwinaCppBA$Species)

corwinaCppBA$Sensitivity <- factor(corwinaCppBA$Sensitivity,
                                   levels = c("2023", "1912 (25%)", "1912 (50%)", "1912 (75%)"))


corwinaCPPbaPlot <- ggplot(corwinaCppBA, aes(x = Sensitivity, y = BasalArea, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits = c(0, 43),
    breaks = seq(0, 43, by = 10),
    expand = c(0, 0)
  ) +
  labs(title = "Corwina CPP Basal Area",
       x = NULL,
       y = "Basal Area (m²/ha)") +
  scale_fill_manual(
    values = c("PIPO" = "#d8b87c",   # warm orange
               "PSME" = "#9eae4c")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_blank(),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(face = "bold", color = "black")
  )

ggsave("corwinaCPP_BA.png", corwinaCPPbaPlot, width = 8, height = 7, dpi = 300)


# Corwina C2

corwinaC2BA <- corwinaC2 %>% 
  select(Sensitivity, refYear, PIPO.ba, PSME.ba) %>% 
  pivot_longer(cols = c(PIPO.ba, PSME.ba),
               names_to = "Species",
               values_to = "BasalArea")

corwinaC2BA$Species <- sub("\\.ba$", "", corwinaC2BA$Species)

corwinaC2BA$Sensitivity <- factor(corwinaC2BA$Sensitivity,
                                   levels = c("2023", "1912 (25%)", "1912 (50%)", "1912 (75%)"))


corwinaC2baPlot <- ggplot(corwinaC2BA, aes(x = Sensitivity, y = BasalArea, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits = c(0, 43),
    breaks = seq(0, 43, by = 10),
    expand = c(0, 0)
  ) +
  labs(title = "Corwina C2 Basal Area",
       x = NULL,
       y = "Basal Area (m²/ha)") +
  scale_fill_manual(
    values = c("PIPO" = "#d8b87c",   # warm orange
               "PSME" = "#9eae4c")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_blank(),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(face = "bold", color = "black")
  )

ggsave("corwinaC2_BA.png", corwinaC2baPlot, width = 8, height = 7, dpi = 300)



# Density Plots ########################
# Corwina CPP

corwinaCppDensity <- corwinaCpp %>% 
  select(Sensitivity, refYear, PIPO.density, PSME.density) %>% 
  pivot_longer(cols = c(PIPO.density, PSME.density),
               names_to = "Species",
               values_to = "Density")

corwinaCppDensity$Species <- sub("\\.density$", "", corwinaCppDensity$Species)

corwinaCppDensity$Sensitivity <- factor(corwinaCppDensity$Sensitivity,
                                   levels = c("2023", "1912 (25%)", "1912 (50%)", "1912 (75%)"))


corwinaCPPdensityPlot <- ggplot(corwinaCppDensity, aes(x = Sensitivity, y = Density, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits = c(0, 1225),
    breaks = seq(0, 1225, by = 200),
    expand = c(0, 0)
  ) +
  labs(title = "Corwina CPP Tree Density",
       x = NULL,
       y = "Tree Density (stems/ha)") +
  scale_fill_manual(
    values = c("PIPO" = "#d8b87c",   # warm orange
               "PSME" = "#9eae4c")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_blank(),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(face = "bold", color = "black")
  )

ggsave("corwinaCPP_Density.png", corwinaCPPdensityPlot, width = 8, height = 7, dpi = 300)


# Corwina C2

corwinaC2Density <- corwinaC2 %>% 
  select(Sensitivity, refYear, PIPO.density, PSME.density) %>% 
  pivot_longer(cols = c(PIPO.density, PSME.density),
               names_to = "Species",
               values_to = "Density")

corwinaC2Density$Species <- sub("\\.density$", "", corwinaC2Density$Species)

corwinaC2Density$Sensitivity <- factor(corwinaC2Density$Sensitivity,
                                        levels = c("2023", "1912 (25%)", "1912 (50%)", "1912 (75%)"))


corwinaC2densityPlot <- ggplot(corwinaC2Density, aes(x = Sensitivity, y = Density, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits = c(0, 1225),
    breaks = seq(0, 1225, by = 200),
    expand = c(0, 0)
  ) +
  labs(title = "Corwina C2 Tree Density",
       x = NULL,
       y = "Tree Density (stems/ha)") +
  scale_fill_manual(
    values = c("PIPO" = "#d8b87c",   # warm orange
               "PSME" = "#9eae4c")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_blank(),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(face = "bold", color = "black")
  )

ggsave("corwinaC2_Density.png", corwinaC2densityPlot, width = 8, height = 7, dpi = 300)



