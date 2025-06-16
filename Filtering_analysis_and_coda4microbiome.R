library(MicrobiotaProcess) 
library(phyloseq) 
library(ggplot2) 
library(tidyverse) 
library(vegan)
library(coin) 
library(reshape2) 
library(ggnewscale) 
library(ggplot2)
# Data load
data <- read.delim("/Users/manuelfernandezlemos/Desktop/Fecal Transplant First Analysis/newDataSet20240920.tsv", header = TRUE, sep = "\t")
load('R_objects/metadataWithOrderTable.RData')
load('R_objects/metadataWithFamilyTable.RData')
load('R_objects/metadataWithGenusTable.RData')
load('R_objects/metadataWithSpeciesTable.RData')




data <- data %>%
  mutate(
    bmi_category = case_when(
      baseline_bmi < 18.5 ~ "underweight",
      baseline_bmi >= 18.5 & baseline_bmi <= 25 ~ "healthy range",
      baseline_bmi > 25 & baseline_bmi <= 29.9 ~ "overweight",
      baseline_bmi >= 30 & baseline_bmi <= 39.9 ~ "obesity",
      baseline_bmi >= 40 ~ "severe obesity"
    )
  )


df<- data %>%
  filter(age>=60) %>%
  filter(w1_ab_cdi==FALSE) %>%
  filter(w1_ab_noncdi==FALSE) %>%
  filter(w8_ab_cdi!= 'NA')

donor_efficiency <- df %>%
  group_by(Sample) %>%
  summarize(
    efficiency = mean(outcome == TRUE),
    number_of_donations = n()
  )

#Check everything looks right
print(donor_efficiency)


newTry <- donor_efficiency
newTry <- donor_efficiency %>%
  filter(number_of_donations>=10 )



newTry<-newTry %>%
  rename(patient=Sample)

#Check everything looks right
View(newTry)


# Join microbiome data with patient metadata
correlationMatrix<-left_join(order, newTry, 'patient')

# Take out troublesome columns
correlationMatrix<-correlationMatrix %>%
  filter(efficiency!='NA')
correlationMatrix$Fusobacteriales<-NULL
correlationMatrix$fmt_GROUPS<-NULL
correlationMatrix$depth<-NULL
correlationMatrix$Ursocholic.acid<-NULL
correlationMatrix$number_of_donations<-NULL

correlationMatrix <- correlationMatrix %>%
  rename(
    `3-oxoCDCA` = X3.Oxochenodeoxycholic.acid
  )

correlationMatrix <- correlationMatrix %>%
  rename(
    `GCA` = Glycocholic.acid
  )

correlationMatrix <- correlationMatrix %>%
  rename(
    `Efficiency` = efficiency
  )
#Create correlation matrix
correlation_matrix <- cor(correlationMatrix[,54:ncol(correlationMatrix)], method='spearman')

library(reshape2)

# Use melt from reshape2 to create a long format dataframe
cor_data_melted <- melt(correlation_matrix)

#Check everything looks right
efficiency_cor <- cor_data_melted %>%
  filter(Var1 == 'Efficiency' | Var2 == 'Efficiency')

View(efficiency_cor)


# Filter the melted data to include only strong correlations
strong_correlations <- cor_data_melted[abs(cor_data_melted$value) >= 0.6 & cor_data_melted$Var1 != cor_data_melted$Var2, ]

#Selecting variables with cor>=0.6
selected_vars <- c("3-oxoCDCA", "Efficiency", "Clostridiales", "GCA")
# Subset the correlation matrix
selected_cor_matrix <- correlation_matrix[selected_vars, selected_vars]
cor_data_melted<-melt(selected_cor_matrix)

# Plot strong correlations
library(ggplot2)
heatmap_plot<-ggplot(cor_data_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_fixed()


ggsave(filename = "correlation_heatmap.svg", plot = heatmap_plot, device = "svg", width = 6, height = 5)


# Order analysis

#Load data if not loaded yet
load('R_objects/metadataWithOrderTables.RData')

#Load packages if not loaded yet
library(ggpubr)
library(dplyr)
library(rstatix)


#For Clostridiales

# 1) Prepare the data
dat.clean <- order %>%
  mutate(
    week = case_when(
      grepl("week_00", fmt_GROUPS) ~ "Week 0",
      grepl("week_01", fmt_GROUPS) ~ "Week 1",
      grepl("week_08", fmt_GROUPS) ~ "Week 8",
      fmt_GROUPS == "donor"        ~ "Donor"
    ),
    status = case_when(
      grepl("failure", fmt_GROUPS) ~ "Failure",
      grepl("effect",  fmt_GROUPS) ~ "Effect",
      fmt_GROUPS == "donor"        ~ "Donor"
    )
  ) %>%
  mutate(
    week   = factor(week,   levels = c("Week 0", "Week 1", "Week 8", "Donor")),
    status = factor(status, levels = c("Failure", "Effect", "Donor"))
  )

# 2) Run Pairwise Wilcoxon Tests (Failure vs Effect by week)
stat.test <- dat.clean %>%
  filter(status %in% c("Failure","Effect")) %>%
  group_by(week) %>%
  wilcox_test(Clostridiales ~ status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

# Key Trick: use 'dodge=0' so the bracket center is exactly between the two boxes
stat.test <- stat.test %>%
  add_xy_position(x = "week", dodge = 0)

# Make a custom label with "p-value="
stat.test <- stat.test %>%
  mutate(pval_label = paste0("p-value=", p))

# 3) Build the boxplot
custom_palette <- c("Failure" = "red", "Effect" = "darkgreen", "Donor" = "blue")
bxp9 <- ggboxplot(
  data = dat.clean,
  x    = "week",
  y    = "Clostridiales",
  color= "status",
  palette = custom_palette
) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "gray40") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8)
  )

# 4) Add centered p-values, remove bracket lines
clostridiales <- bxp9 +
  stat_pvalue_manual(
    stat.test,
    label         = "pval_label",
    remove.bracket= TRUE,      # Hide the bracket line
    tip.length    = 0.01,
    step.increase = 0.1
  )

print(clostridiales)
ggsave(filename = "clostridiales.svg", plot = clostridiales, device = "svg", width = 6, height = 5)


#Lactobacillales

# 1) Prepare the data
dat.clean <- order %>%
  mutate(
    week = case_when(
      grepl("week_00", fmt_GROUPS) ~ "Week 0",
      grepl("week_01", fmt_GROUPS) ~ "Week 1",
      grepl("week_08", fmt_GROUPS) ~ "Week 8",
      fmt_GROUPS == "donor"        ~ "Donor"
    ),
    status = case_when(
      grepl("failure", fmt_GROUPS) ~ "Failure",
      grepl("effect",  fmt_GROUPS) ~ "Effect",
      fmt_GROUPS == "donor"        ~ "Donor"
    )
  ) %>%
  mutate(
    week   = factor(week,   levels = c("Week 0", "Week 1", "Week 8", "Donor")),
    status = factor(status, levels = c("Failure", "Effect", "Donor"))
  )

# 2) Run Pairwise Wilcoxon Tests (Failure vs Effect by week)
stat.test <- dat.clean %>%
  filter(status %in% c("Failure","Effect")) %>%
  group_by(week) %>%
  wilcox_test(Lactobacillales ~ status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

# Key Trick: use 'dodge=0' so the bracket center is exactly between the two boxes
stat.test <- stat.test %>%
  add_xy_position(x = "week", dodge = 0)

# Make a custom label with "p-value="
stat.test <- stat.test %>%
  mutate(pval_label = paste0("p-value=", p))

# 3) Build the boxplot
custom_palette <- c("Failure" = "red", "Effect" = "darkgreen", "Donor" = "blue")
bxp9 <- ggboxplot(
  data = dat.clean,
  x    = "week",
  y    = "Lactobacillales",
  color= "status",
  palette = custom_palette
) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "gray40") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8)
  )

# 4) Add centered p-values, remove bracket lines
lactobacillales <- bxp9 +
  stat_pvalue_manual(
    stat.test,
    label         = "pval_label",
    remove.bracket= TRUE,      # Hide the bracket line
    tip.length    = 0.01,
    step.increase = 0.1
  )
print(lactobacillales)

ggsave(filename = "lactobacillales.svg", plot = lactobacillales, device = "svg", width = 6, height = 5)
#coda4microbiome

library(coda4microbiome)

# Running the algorithms leads to different results, 
# although the algorithm is always able to create a higher efficiency group. 
# For reproducibility we provide the values that were used in the paper.

xOriginal<-correlationMatrix[72:108]
codaPredictions<-coda_glmnet(xOriginal, correlationMatrix$Efficiency)

codaPredictions$`predictions plot`

# saveRDS(codaPredictions, file = "coda_predictions.rds")


codaPredictions <- readRDS("coda_predictions.rds")



predicciones <-log(as.matrix(correlationMatrix[codaPredictions$taxa.name])) %*% codaPredictions$`log-contrast coefficients`
predicciones



efficiency<-correlationMatrix$Efficiency
efficiency


# Combine these into a data frame for ggplot
table <- data.frame(
  predicciones = predicciones,
  efficiency = efficiency
)

corr_test <- cor.test(table$predicciones, table$efficiency)
r_value <- corr_test$estimate          
r2_value <- r_value^2  
p_value <- corr_test$p.value
p_label <- format(p_value, scientific = TRUE, digits = 2)

# Create the plot
efficiency_coda<-ggplot(table, aes(x = predicciones, y = efficiency)) +
  geom_point(color = "#1f77b4", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "#d62728", fill = "#ffa07a", alpha = 0.3) +
  # Add the R² and p-value as text annotation
  annotate("text", x = min(table$predicciones) + 0.02, 
           y = max(table$efficiency) - 0.05,
           label = paste0("R² = ", round(r2_value, 2), ", p = ", p_label),
           hjust = 0, size = 5) +
  labs(x = "Predictions", y = "Efficiency score (Sₑ)") +
  theme_minimal(base_size = 14)


#Apply algorithm to another filter
ggsave(filename = "efficiency_coda.svg", plot = efficiency_coda, device = "svg", width = 6, height = 5)

df2<- df %>%
  filter(age>=60) %>%
  filter(outcomeFull != 'Death') %>%
  filter(outcomeFull != 'Suspected recurrence') %>%
  filter(charlson_score <=8) %>%
  filter(bmi_category != "underweight") %>%
  filter(w1_ab_cdi==FALSE) %>%
  filter(w1_ab_noncdi==FALSE) %>%
  filter(w8_ab_cdi!= 'NA')

donor_efficiency <- df2 %>%
  group_by(Sample) %>%
  summarize(
    efficiency = mean(outcome == TRUE),
    number_of_donations = n()
  )

View(donor_efficiency)

# Check everything looks right

newTry <- donor_efficiency
newTry <- donor_efficiency %>%
  filter(number_of_donations>=5 )
newTry<-newTry %>%
  rename(patient=Sample)

View(newTry)


# Join microbiome data with patient metadata
correlationMatrix2<-left_join(order, newTry, 'patient')

# Take out troublesome columns
correlationMatrix2<-correlationMatrix2 %>%
  filter(efficiency!='NA')
correlationMatrix2$Fusobacteriales<-NULL
correlationMatrix2$fmt_GROUPS<-NULL
correlationMatrix2$depth<-NULL
correlationMatrix2$Ursocholic.acid<-NULL
correlationMatrix2$number_of_donations<-NULL

correlationMatrix2 <- correlationMatrix2 %>%
  rename(
    `3-oxoCDCA` = X3.Oxochenodeoxycholic.acid
  )

correlationMatrix2 <- correlationMatrix2 %>%
  rename(
    `GCA` = Glycocholic.acid
  )

correlationMatrix2 <- correlationMatrix2 %>%
  rename(
    `Efficiency` = efficiency
  )



prediccionesNueva <- log(as.matrix(correlationMatrix2[codaPredictions$taxa.name])) %*% codaPredictions$`log-contrast coefficients`

correlationMatrix2$EfficiencyLabel <- ifelse(prediccionesNueva > 0, "High Efficiency", "Low Efficiency")

dat.clean <- correlationMatrix2


library(ggsci)

PREDICTOR <- "EfficiencyLabel"
OUTCOME <- "Efficiency"
SUBJECT <- "patient"



dat.overall <- dat.clean %>% mutate(EfficiencyLabel = "Overall")

dat.for.plot <- bind_rows(dat.overall, dat.clean) %>%
  mutate(EfficiencyLabel = factor(
    EfficiencyLabel,
    levels = c("High Efficiency", "Low Efficiency", "Overall")
  ))


## Make the plot

efficiency_bxp <- dat.for.plot %>%
  ggboxplot(
    x      = PREDICTOR,
    y      = OUTCOME,
    color  = PREDICTOR,
    palette = "jco"
  ) +

  # zoom the y-axis
  coord_cartesian(ylim = c(0.55, 0.95)) +

scale_color_jco(
  name   = "Efficiency Label",          # ← legend title
  labels = c(                           # ← legend entries
    expression("High " * S[e]),
    expression("Low "  * S[e]),
    "Overall"
  )
) +
  scale_fill_jco(guide = "none") +        # keep fills in sync, hide its guide

  # remove everything from the x-axis
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    legend.text  = element_text(size = 8)
  )

print(efficiency_bxp)

ggsave(filename = "efficiency_bxp.svg", plot = efficiency_bxp, device = "svg", width = 6, height = 5)

library(ggplot2)
library(ggpubr)



dat.overall <- dat.clean %>% mutate(EfficiencyLabel = "Overall")

dat.for.plot <- bind_rows(dat.overall, dat.clean) %>%
  mutate(EfficiencyLabel = factor(
    EfficiencyLabel,
    levels = c("High Efficiency", "Low Efficiency", "Overall")
  ))

dat.for.plot <- bind_rows(dat.overall, dat.clean) %>%
  mutate(EfficiencyLabel = factor(
    EfficiencyLabel,
    levels = c("High Efficiency", "Low Efficiency", "Overall")
  ))


efficiency_bxp <- ggplot(dat.for.plot, aes(x = EfficiencyLabel, y = Efficiency, color = EfficiencyLabel)) +
  geom_boxplot() +
  
  # zoom the y-axis
  coord_cartesian(ylim = c(0.55, 0.95)) +
  
  # Set the y-axis title
  ylab(expression("Efficiency Score (S"["e"]*")")) +
  
  # Apply the jco color palette and set legend properties
  scale_color_jco(
    name   = "Efficiency Label",          # legend title
    labels = c(                           # legend entries
      expression("High " * S[e]),
      expression("Low "  * S[e]),
      "Overall"
    )
  ) +
  
  # Apply a clean, white theme
  theme_classic() +
  
  # remove everything from the x-axis and format legend
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    
    legend.text  = element_text(size = 8)
  )

print(bxp9)
  
ggsave(filename = "efficiency_bxp.svg", plot = efficiency_bxp, device = "svg", width = 6, height = 5)




# --- Calculate and Print Summary Statistics ---
efficiency_summary <- dat.for.plot %>%
  group_by(EfficiencyLabel) %>%
  summarise(
    mean_efficiency = mean(Efficiency, na.rm = TRUE),
    median_efficiency = median(Efficiency, na.rm = TRUE),
    sd_efficiency = sd(Efficiency, na.rm = TRUE),
    count = n()
  )

View(efficiency_summary)


high_vs_low_data <- dat.for.plot %>%
  filter(EfficiencyLabel %in% c("High Efficiency", "Low Efficiency"))


stat_test_result <- t.test(Efficiency ~ EfficiencyLabel, data = high_vs_low_data)


print("Statistical comparison between High and Low Efficiency groups:")
print(stat_test_result)



library(dplyr)
library(ggpubr)
library(rstatix)



#Distances Bray Curtis
params <- readRDS("R_objects/params_betadiv.RDS")

# Choose metric
METRIC <- "bray"

# Choose variable 
VAR <- "fmt_GROUPS"

# Load data
load(params$input)

if (METRIC == "unif") {
  load("R_objects/bdiv_unif.RData")
  dist.used <- unif.dist
  nmds.used <- unif.nmds
  pcoa.used <- unif.pcoa
  phy.used <- phy.rare
  rm(unif.dist, unif.nmds, unif.pcoa, phy.rare)
} else if (METRIC == "wunif") {
  load("R_objects/bdiv_wunif.RData")
  dist.used <- wuf.dist
  nmds.used <- wuf.nmds
  pcoa.used <- wuf.pcoa
  phy.used <- phy.clean
  rm(wuf.dist, wuf.nmds, wuf.pcoa,phy.clean)
} else if (METRIC == "bray"){
  load("R_objects/bdiv_bray.RData")
  dist.used <- bray.dist
  nmds.used <- bray.nmds
  pcoa.used <- bray.pcoa
  phy.used <- phy.ra
  rm(bray.dist, bray.nmds, bray.pcoa, phy.ra)
} else if (METRIC == "jac"){
  load("R_objects/bdiv_jac.RData")
  dist.used <- jac.dist
  nmds.used <- jac.nmds
  pcoa.used <- jac.pcoa
  phy.used <- phy.rare
  rm(jac.dist, jac.nmds, jac.pcoa, phy.rare)
} else if (METRIC == "ait"){
  load("R_objects/bdiv_ait.RData")
  dist.used <- ait.dist
  nmds.used <- ait.nmds
  pcoa.used <- ait.pcoa
  phy.used <- phy.clean
  rm(ait.dist, ait.nmds, ait.pcoa, phy.clean)
}

# Extract metadata from phyloseq
mdat <-sample_data(phy.used)
mdat<-data.frame(mdat)
mdat[,VAR] <- as.factor(mdat[,VAR])


library(metagMisc)
library(dplyr)
library(stringr)
library(rstatix)

# Convert distance object to long format
dist.list <- dist2list(dist.used, tri=FALSE)


 # Annotate distance list with group info
dist.list2 <- dist.list %>%
  rename(SampleID_col = col, SampleID_row = row) %>%
  left_join(mdat %>% select(SampleID, fmt_GROUPS),
            by = c("SampleID_col" = "SampleID")) %>%
  rename(Group_col = fmt_GROUPS) %>%
  left_join(mdat %>% select(SampleID, fmt_GROUPS),
            by = c("SampleID_row" = "SampleID")) %>%
  rename(Group_row = fmt_GROUPS)

# Filter donor <-> non-donor pairs, compute mean distance per patient sample
dist_vs_donor <- dist.list2 %>%
  filter(
    (Group_col == "donor" & Group_row != "donor") |
      (Group_row == "donor" & Group_col != "donor")
  )

mean_dist_to_donors <- dist_vs_donor %>%
  mutate(
    patientSample = if_else(Group_col == "donor", SampleID_row, SampleID_col),
    patientGroup  = if_else(Group_col == "donor", Group_row, Group_col)
  ) %>%
  group_by(patientSample, patientGroup) %>%
  summarise(mean_dist = mean(value), .groups = "drop")

# 3) Boxplot of mean distance, grouped by fmt_GROUPS
boxplot_data <- mean_dist_to_donors %>%
  mutate(patientGroup = factor(patientGroup,
                               levels = c("week_00_failure","week_00_effect",
                                          "week_01_failure","week_01_effect",
                                          "week_08_failure","week_08_effect",
                                          "donor")))
# Separate week & status from patientGroup
boxplot_data2 <- boxplot_data %>%
  separate(patientGroup, into=c("week","status"), sep="_(?=failure|effect)",
           remove=FALSE, extra="merge") %>%
  filter(status %in% c("failure","effect"))

library(grid)   # for unit()

p <- ggboxplot(
  data    = boxplot_data2,
  x       = "week",
  y       = "mean_dist",
  color   = "status",
  palette = c("failure" = "red", "effect" = "darkgreen")
) +
  # rotate & relabel
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mean distance") +
  scale_x_discrete(labels = c(
    "week_00" = "Week 0",
    "week_01" = "Week 1",
    "week_08" = "Week 8"
  )) +
  
  coord_cartesian(clip = "off") +
  
  theme(
    plot.margin = unit(c(1, 1, 1.5, 1), "cm")
  )


stat.test <- boxplot_data2 %>%
  group_by(week) %>%
  wilcox_test(mean_dist ~ status) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p") %>%
  add_xy_position(x="week", dodge=0.8) %>%
  mutate(pval_label = paste0("p-value=", p))

distance_plot <- p +
  stat_pvalue_manual(
    stat.test,
    label      = "pval_label",
    tip.length = 0.01
  ) +
  labs(x = NULL)

print(distance_plot)

ggsave(filename = "distance_plot.svg", plot = distance_plot, device = "svg", width = 6, height = 5)

#Observed

library(ggpubr)
library(dplyr)
library(rstatix)



order <- order %>%
  rename(
    `Observed richness` = Observed
  )


# 1) Prepare the data
dat.clean <- order %>%
  mutate(
    week = case_when(
      grepl("week_00", fmt_GROUPS) ~ "Week 0",
      grepl("week_01", fmt_GROUPS) ~ "Week 1",
      grepl("week_08", fmt_GROUPS) ~ "Week 8",
      fmt_GROUPS == "donor"        ~ "Donor"
    ),
    status = case_when(
      grepl("failure", fmt_GROUPS) ~ "Failure",
      grepl("effect",  fmt_GROUPS) ~ "Effect",
      fmt_GROUPS == "donor"        ~ "Donor"
    )
  ) %>%
  mutate(
    week   = factor(week,   levels = c("Week 0", "Week 1", "Week 8", "Donor")),
    status = factor(status, levels = c("Failure", "Effect", "Donor"))
  )

# 2) Run Pairwise Wilcoxon Tests (Failure vs Effect by week)
stat.test <- dat.clean %>%
  filter(status %in% c("Failure","Effect")) %>%
  group_by(week) %>%
  wilcox_test(`Observed richness` ~ status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

# Key Trick: use 'dodge=0' so the bracket center is exactly between the two boxes
stat.test <- stat.test %>%
  add_xy_position(x = "week", dodge = 0)

# Make a custom label with "p-value="
stat.test <- stat.test %>%
  mutate(pval_label = paste0("p-value=", p))

# 3) Build the boxplot
custom_palette <- c("Failure" = "red", "Effect" = "darkgreen", "Donor" = "blue")
bxp9 <- ggboxplot(
  data = dat.clean,
  x    = "week",
  y    = "Observed richness",
  color= "status",
  palette = custom_palette
) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", color = "gray40") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8)
  )

# 4) Add centered p-values, remove bracket lines
observed_richness <- bxp9 +
  stat_pvalue_manual(
    stat.test,
    label         = "pval_label",
    remove.bracket= TRUE,      # Hide the bracket line
    tip.length    = 0.01,
    step.increase = 0.1
  )

print(observed_richness)




ggsave(filename = "observed_richness.svg", plot = observed_richness, device = "svg", width = 6, height = 5)

