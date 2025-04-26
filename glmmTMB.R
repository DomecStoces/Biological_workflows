library(glmmTMB)
library(lme4)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggeffects)
library(car)
library(emmeans)
library(sjPlot)
library(SYNCSA)
library(cluster)
library(effects)

# Histogram Number for each functional group of ground-dwelling arthropods
ggplot(dataset2[dataset2$Functional.group == "Detritivore", ], aes(x = Abundance)) +
  geom_histogram(binwidth = 1, color = "black", fill = "blue") +
  labs(
    title = "Histogram of Abundance",
    x = "Abundance",
    y = "Frequency"
  ) +
  theme_minimal()

# Calculate mean and variance for "Number" by Functional.group
group_stats <- dataset6 %>%
  filter(Functional.group %in% c("Detritivore", "Herbivore", "Omnivore", "Predator", "Saproxylic")) %>%
  group_by(Functional.group)

summarise(
  Mean_Number = mean(Number),
  Variance_Number = var(Number),
  Overdispersion = Variance_Number / Mean_Number
)

# Print the results
print(group_stats)

#What family to choose?
# Fit the Poisson model
poisson_model <- glmmTMB(SpeciesRichness ~ Treatment + Movement.pattern * Season + (1 | Trap), 
                         data = species_richness_data, 
                         family = poisson(link = "log"))
poisson_model <- glmmTMB(Number~ Treatment*Movement.pattern + (1 | Trap)+(1|Month), 
                         data = dataset7, 
                         family = poisson(link = "log"))

# Calculate overdispersion statistic
overdispersion_stat <- sum(residuals(poisson_model, type = "pearson")^2) / df.residual(poisson_model)

# Print the overdispersion statistic
print(overdispersion_stat)

# Categorical variable with Trap ID (unique for each trap) same as Clearing (1 or 2) and Movement pattern
dataset6$Clearing <- as.factor(dataset6$Clearing)
dataset6$Trap <- as.factor(dataset6$Trap)
dataset6$Movement.pattern <- as.factor(dataset6$Movement.pattern)
dataset6$Treatment <- as.factor(dataset6$Treatment)

#Total Abundance of ground-dwelling arthropods

dataset6$Treatment <- gsub("\\.", " ", dataset6$Treatment)
dataset6$Treatment <- factor(dataset6$Treatment, 
                             levels = c("Forest interior", "Ecotone", "Retention clearcut"))
levels(dataset6$Treatment)

model9 <- glmmTMB(Number ~ Treatment*Movement.pattern+ (1 | Trap)+(1|Month), 
                  data = dataset6, 
                  family = nbinom2(link = "sqrt"))
Anova(model9,type="III"))
emm <- emmeans(model9, ~ Movement.pattern*Treatment, type = "response")
emm_df <- as.data.frame(emm)

d<-ggplot(emm_df, aes(x = Treatment, y = response, 
                      color = Movement.pattern, 
                      group = Movement.pattern)) +
  # Add points with dodge for separation
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  # Add error bars with dodge and custom linewidth
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                linewidth = 0.8, 
                position = position_dodge(width = 0.5)) +
  # Add lines connecting points
  geom_line(aes(linetype = Movement.pattern), 
            position = position_dodge(width = 0.5),
            linewidth = 0.8) +
  # Customize labels
  labs(
    x = "Treatment",
    y = "Total predicted N° of individuals"
  ) +
  # Customize Y-axis
  scale_y_continuous(
    limits = c(0, 6), # Adjust range if needed
    labels = function(x) ifelse(x == 0, "0", x) 
  ) +
  # Customize theme
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text = element_text(size = 15, family = "Arial"), 
    strip.background = element_rect(fill = "grey90"), 
    axis.text.x = element_text(size = 12, family = "Arial", angle = -45, hjust = 0, vjust = 1),
    axis.text.y = element_text(size = 12, family = "Arial"),
    axis.title.x = element_text(size = 15, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial", margin = margin(r = 10)),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_blank(), 
    legend.position = "right", 
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black") 
  ) +
  # Customize color and line types
  scale_color_manual(values = c("Across" = "grey60", "Along" = "black")) + 
  scale_linetype_manual(values = c("Across" = "dashed", "Along" = "solid"))

tiff('Total_abundance.tiff',units="in",width=6,height=5,bg="white",res=300)
d
dev.off()

,expand = expansion(mult = c(0))
emmeans_results <- emmeans(model9, ~ Movement.pattern|Treatment)

# Apply pairwise contrasts with Sidak adjustment: To identify specific differences between levels of categorical predictors.
contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
summary(contrast_results)

##########################################################################################
# Check for sandardised residuals
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model9)

# Plot diagnostics
plot(sim_res)

# Run outlier test with bootstrap method
testOutliers(sim_res, type = "bootstrap")

#Likelihood Ratio Test (LRT) for Predictor Significance: For testing the overall significance of predictors or interactions
# Full model
model_full <- glmmTMB(Number ~ Season* Movement.pattern+ (1 | Trap), 
                      data = dataset6, 
                      family = nbinom2(link = "sqrt"))

Anova(model_full,type="III")
emmeans_results <- emmeans(model_full, ~ Movement.pattern|Treatment)

# Apply pairwise contrasts with Sidak adjustment: To identify specific differences between levels of categorical predictors.
contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
summary(contrast_results)

##########################################################################################
#MODEL variant for Treatment*Movement.pattern+(1|Month) separately for each Functional group
# List to store models
models <- list()

# Loop through each functional group
groups <- unique(dataset6$Alternative)

for (group in groups) {
  # Filter data for the current functional group
  group_data <- dataset6 %>% filter(Alternative == group)
  
  # Fit Negative Binomial GLMM (use group_data, not dataset)
  models[[group]] <- glmmTMB(Number ~ Treatment*Movement.pattern+ (1 |Trap)+(1|Month), 
                             data = group_data, 
                             family = nbinom2(link = "sqrt"))
  
  # Print a message after fitting the model for this group
  cat("Model fitted for Functional Group:", group, "\n")
}

# Access a summary for a specific group (e.g., "Detritivore")
summary(models[["Detritivore"]])
Anova(models[["omnivore"]],type="III")

# Estimated marginal means from the model
emm <- emmeans(models[["omnivore"]], ~ Movement.pattern * Treatment, type = "response")
emm_df <- as.data.frame(emm)

d <- ggplot(emm_df, aes(x = Treatment, y = response, 
                        color = Movement.pattern, 
                        group = Movement.pattern)) +
  # Add points
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  
  # Add error bars
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  
  # Add lines
  geom_line(aes(linetype = Movement.pattern), 
            position = position_dodge(width = 0.5),
            linewidth = 0.8) +
  
  # Axis labels
  labs(
    y = "Predicted N° of individuals",
    x = NULL
  ) +
  
  # Y-axis settings
  scale_y_continuous(
    limits = c(0, 14.5),
    labels = function(x) ifelse(x == 0, "0", x)
  ) +
  
  # Color and line type customizations
  scale_color_manual(values = c("Across" = "grey60", "Along" = "black")) + 
  scale_linetype_manual(values = c("Across" = "dashed", "Along" = "solid")) +
  
  # Theme settings
  theme_minimal(base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    strip.text = element_text(size = 25),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(size = 20, angle = -45, hjust = 0, vjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 3)
  )

# Print plot
print(d)

#Overall Species Richness
species_richness_data <- dataset6 %>%
  group_by(Trap, Treatment, Month, Movement.pattern) %>%
  summarize(
    SpeciesRichness = n_distinct(Species),  # Count of unique species
    Abundance = sum(Number, na.rm = TRUE),  # Sum of abundances (Number)
    .groups = "drop"
  )

# Fit the model
model5 <- glmmTMB(SpeciesRichness ~ Treatment + Movement.pattern + (1 | Trap)+(1|Month), 
                  data = species_richness_data, 
                  family = nbinom2(link = "log"))
summary(model5)

Anova(model2)

species_richness_data1 <- dataset3 %>%
  group_by(Trap, Treatment, Season, Movement.pattern) %>%
  summarize(
    SpeciesRichness = n_distinct(Species),  # Count of unique species
    Abundance = sum(Number, na.rm = TRUE),  # Sum of abundances (Number)
    .groups = "drop")
model5 <- glmmTMB(SpeciesRichness ~ Treatment + Movement.pattern*Season + (1 | Trap), 
                  data = species_richness_data1, 
                  family = nbinom2(link = "log"))

##########################################################################################
# MODEL for SpeciesRichness ~ Treatment * Movement.pattern+(1|Month) for each Functional group
species_richness_data <- dataset6 %>%
  group_by(Trap, Treatment, Month, Movement.pattern, Functional.group) %>%
  summarize(
    SpeciesRichness = n_distinct(Species),  # Unique species count
    Abundance = sum(Number, na.rm = TRUE),  # Sum of individuals
    Number = sum(Number, na.rm = TRUE),     # Retain Number
    .groups = "drop"
  )

# Calculating Species richness for Variant Treatment*Movement.pattern+(1|Month)
# Initialize a list to store models
models2 <- list()

# Get unique Functional Groups
functional_groups2 <- unique(species_richness_data$Functional.group)

# Loop through each Functional Group
for (group2 in functional_groups2) {
  # Filter data for the current Functional Group
  group_data2 <- species_richness_data %>% filter(Functional.group == group2)
  
  # Check for missing or insufficient data
  if (nrow(group_data2) < 5) {
    cat("Skipping Functional Group:", group2, "due to insufficient data\n")
    next
  }
  
  # Fit the Negative Binomial GLMM
  tryCatch({
    model <- glmmTMB(
      SpeciesRichness ~ Treatment * Movement.pattern + (1 | Trap) + (1 | Month),
      data = group_data2,
      family = nbinom2(link = "sqrt")
    )
    
    # Store the model
    models2[[group2]] <- model
    
    # Print a success message
    cat("Model fitted for Functional Group:", group2, "\n")
  }, error = function(e) {
    cat("Error fitting model for Functional Group:", group2, ": ", e$message, "\n")
  })
}

summary(models2[["Predator"]])
Anova(models2[["Herbivore"]],type="III")

emmeans_results <- emmeans(models2[["Herbivore"]], ~ Movement.pattern | Treatment)
contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
summary(contrast_results)

plot(allEffects(models[["Detritivore"]]))

# Check the number of observations for "Detritivore"
nobs(models[["Detritivore"]])

# Compare with the number of rows in the filtered dataset
nrow(dataset6 %>% filter(Functional.group == "Detritivore"))

# Extract the Detritivore model
predator_model <- models[["Predator"]]

# Generate predictions with confidence intervals
predictions <- ggpredict(predator_model, terms = c("Treatment", "Movement.pattern", "Season"))

# Convert to a data frame for plotting
results_data <- as.data.frame(predictions)

# Reorder Season levels
results_data$facet <- factor(results_data$facet, levels = c("Spring", "Summer", "Autumn"))

# Reorder Movement Pattern levels
results_data$group <- factor(results_data$group, levels = c("Along", "Across"))

d<-ggplot(results_data, aes(x = facet, y = predicted, color = group, group = group)) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +  # Points for predictions
  geom_line(aes(group = interaction(group, x)), position = position_dodge(width = 0.4)) +  # Lines connecting seasons
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +  # Error bars for confidence intervals
  facet_wrap(~ x, scales = "free_x") +  # Separate panels for each Treatment
  labs(
    title = "Predicted abundance for Herbivore by treatment and season",
    x = "Season",
    y = "Number of specimens from negive binomial model",
    color = "Movement pattern"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

tiff('Herbivore_negative_binomial.tiff',units="in",width=8,height=7,bg="white",res=600)
d
dev.off()


# Sum the Number column for each Functional.group
observation_counts <- dataset3 %>%
  group_by(Functional.group) %>%
  summarise(Total_Specimens = sum(Number, na.rm = TRUE))

# View the result
print(observation_counts)

# Verify the total count matches 9355
total_count <- sum(observation_counts$Total_Specimens)
print(total_count)

###############################################################################################
#MODEL for interaction of Number~Movement*Treatment each separated into Season for each Functional group
dataset6$Treatment <- gsub("\\.", " ", dataset6$Treatment)
dataset6$Treatment <- factor(dataset6$Treatment, 
                             levels = c("Forest interior", "Ecotone", "Retention clearcut"))
levels(dataset6$Treatment)

# List to store models
models <- list()

# Get unique Functional Groups
functional_groups <- unique(dataset6$Functional.group)

# Loop through each Functional Group
for (group in functional_groups) {
  # Filter data for the current Functional Group
  group_data <- dataset6 %>% filter(Functional.group == group)
  
  # Get unique Seasons
  seasons <- unique(group_data$Season)
  
  # Initialize a nested list for the group
  models[[group]] <- list()
  
  # Loop through each Season
  for (season in seasons) {
    # Filter data for the current Season
    season_data <- group_data %>% filter(Season == season)
    
    # Fit the Negative Binomial GLMM
    model <- glmmTMB(
      Number ~ Treatment*Movement.pattern + (1 | Trap),
      data = season_data,
      family = nbinom2(link = "sqrt")
    )
    
    # Store the model
    models[[group]][[season]] <- model
    
    # Print a message after fitting the model
    cat("Model fitted for Functional Group:", group, ", and Season:", season, "\n")
  }
}

# Example: Access a summary for a specific Functional Group and Season
# Replace "Detritivore" and "Spring" with your specific values
summary(models[["Predator"]][["Spring"]])

emm <- emmeans(models[["Predator"]][["Spring"]], ~ Movement.pattern*Treatment, type = "response")
emm_df <- as.data.frame(emm)
  
d <- ggplot(emm_df, aes(x = Treatment, y = response, 
                          color = Movement.pattern, 
                          group = Movement.pattern)) +
    # Add points with dodge for separation
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    # Add error bars with dodge and custom linewidth
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                  width = 0.2, 
                  linewidth = 0.8, # Use linewidth instead of size
                  position = position_dodge(width = 0.5)) +
    # Add lines connecting points
    geom_line(aes(linetype = Movement.pattern), 
              position = position_dodge(width = 0.5),
              linewidth = 0.8) + # Use linewidth for line thickness
    # Customize labels (Remove X-axis title)
    labs(
      y = "Predicted N° of individuals"
    ) +
    # Add a facet for a header-like layout
    facet_grid(~"Spring", scales = "free", space = "free") +
    # Customize Y-axis with proper labels and limits
    scale_y_continuous(
      limits = c(0, 11.5), # Replace with your desired range
      labels = function(x) ifelse(x == 0, "0", x) # Ensure 0 is not displayed as 0.0
    ) +
    # Customize theme
    theme_minimal(base_family = "Arial") +
    theme(
      panel.border = element_rect(color = "black", fill = NA), # Add border
      strip.text = element_text(size = 25, family = "Arial"), # Customize facet header text
      strip.background = element_rect(fill = "grey90"), # Light grey background for header
      axis.text.x = element_text(size = 20, family = "Arial", angle = -45, hjust = 0, vjust = 1),
      axis.text.y = element_text(size = 20, family = "Arial"),
      axis.title.x = element_blank(), # REMOVE X-axis title
      axis.title.y = element_blank(), # REMOVE X-axis title,
      legend.position = "none", # Remove legend completely
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      plot.margin = margin(t = 20, r = 20, b = 10, l = 3) # Adjust plot margins
    ) +
    # Customize color and line types
    scale_color_manual(values = c("Across" = "grey60", "Along" = "black")) + 
    scale_linetype_manual(values = c("Across" = "dashed", "Along" = "solid"))
  
tiff('Predator_Spring.tiff', units="in", width=8, height=5, res=500)
d
dev.off()
  
# Apply pairwise contrasts with Sidak adjustment: To identify specific differences between levels of categorical predictors.
contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
summary(contrast_results)

  ###############################################################################################
  #MODEL for interaction of SpeciesRichness~Movement*Treatment divided into Season for each Functional group
  species_richness_data <- dataset6 %>%
    group_by(Trap, Treatment, Season, Movement.pattern, Functional.group) %>%
    summarize(
      SpeciesRichness = n_distinct(Species),  # Count of unique species
      Abundance = sum(Number, na.rm = TRUE),  # Sum of abundances (Number)
      Number = sum(Number, na.rm = TRUE),     # Retain Number
      .groups = "drop"
    )
  
  # Initialize a list to store models
  models2 <- list()
  
  # Get unique Functional Groups
  functional_groups2 <- unique(species_richness_data$Functional.group)
  
  # Loop through each Functional Group
  for (group2 in functional_groups2) {
    # Filter data for the current Functional Group
    group_data2 <- species_richness_data %>% filter(Functional.group == group2)
    
    # Check for missing or inconsistent data
    if (nrow(group_data2) == 0) {
      cat("Skipping Functional Group:", group2, "due to no data\n")
      next
    }
    
    # Get unique Seasons
    seasons <- unique(group_data2$Season)
    
    # Initialize a nested list for the group
    models2[[group2]] <- list()
    
    # Loop through each Season
    for (season in seasons) {
      # Filter data for the current Season
      season_data2 <- group_data2 %>% filter(Season == season)
      
      # Check if sufficient data exists for modeling
      if (nrow(season_data2) < 5) {
        cat("Skipping Functional Group:", group2, "Season:", season, "due to insufficient data\n")
        next
      }
      
      # Fit the Negative Binomial GLMM
      tryCatch({
        model <- glmmTMB(
          SpeciesRichness ~ Treatment * Movement.pattern + (1 | Trap),
          data = season_data2,
          family = nbinom2(link = "sqrt")
        )
        
        # Store the model
        models2[[group2]][[season]] <- model
        
        # Print a message after successfully fitting the model
        cat("Model fitted for Functional Group:", group2, ", and Season:", season, "\n")
      }, error = function(e) {
        cat("Error fitting model for Functional Group:", group2, ", and Season:", season, ": ", e$message, "\n")
      })
    }
  }
  
  # Example: Access a summary for a specific Functional Group and Season
  # Replace "Predator" and "Summer" with your specific values
  if (!is.null(models2[["Predator"]][["Spring"]])) {
    summary(models2[["Predator"]][["Spring"]])
  } else {
    cat("Model for Functional Group 'Predator' and Season 'Summer' is not available.\n")
  }
  summary(models2[["Predator"]][["Spring"]])
  Anova(models2[["Predator"]][["Autumn"]],type = "III")
  
  # Full model
  model_full1 <- glmmTMB(SpeciesRichness ~ Treatment * Movement.pattern + (1 | Trap)+(1|Season), 
                         data = species_richness_data , 
                         family = nbinom2(link = "sqrt"))
  Anova(model_full1,type = "III")
  emmeans_results <- emmeans(model_full1, ~ Movement.pattern|Treatment)
  
  # Apply pairwise contrasts with Sidak adjustment: To identify specific differences between levels of categorical predictors.
  contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
  summary(contrast_results)
  
  emm <- emmeans(models2[["Detritivore"]][["Autumn"]], ~ Treatment*Movement.pattern, type = "response")
  emm_df <- as.data.frame(emm)
  
  d <- ggplot(emm_df, aes(x = Treatment, y = response, 
                          color = Movement.pattern, 
                          group = Movement.pattern)) +
    # Add points with dodge for separation
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    # Add error bars with dodge and custom linewidth
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                  width = 0.2, 
                  linewidth = 0.8, # Use linewidth instead of size
                  position = position_dodge(width = 0.5)) +
    # Add lines connecting points
    geom_line(aes(linetype = Movement.pattern), 
              position = position_dodge(width = 0.5),
              linewidth = 0.8) + # Use linewidth for line thickness
    # Customize labels (Remove X-axis title)
    labs(
      y = "Predicted N° of species"
    ) +
    # Add a facet for a header-like layout
    facet_grid(~"Autumn", scales = "free", space = "free") +
    # Customize Y-axis with proper labels and limits
    scale_y_continuous(
      limits = c(0, 15), # Replace with your desired range
      labels = function(x) ifelse(x == 0, "0", x) # Ensure 0 is not displayed as 0.0
    ) +
    # Customize theme
    theme_minimal(base_family = "Arial") +
    theme(
      panel.border = element_rect(color = "black", fill = NA), # Add border
      strip.text = element_text(size = 25, family = "Arial"), # Customize facet header text
      strip.background = element_rect(fill = "grey90"), # Light grey background for header
      axis.text.x = element_text(size = 20, family = "Arial", angle = -45, hjust = 0, vjust = 1),
      axis.text.y = element_text(size = 20, family = "Arial"),
      axis.title.x = element_blank(), # REMOVE X-axis title
      axis.title.y = element_blank(), # REMOVE X-axis title,
      legend.text = element_text(size = 15, family = "Arial"),
      legend.title = element_blank(), # Remove the legend title
      legend.position = "none", # Move legend to the top
      legend.direction = "horizontal", # Horizontal legend layout
      legend.justification = c(0, 1), # Align legend to the top-left corner
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      plot.margin = margin(t = 20, r = 20, b = 10, l = 3) # Adjust plot margins
    ) +
    # Customize color and line types
    scale_color_manual(values = c("Across" = "grey60", "Along" = "black")) + 
    scale_linetype_manual(values = c("Across" = "dashed", "Along" = "solid"))
  
  ,expand = expansion(mult = c(0))
  element_text(size = 20, family = "Arial", margin = margin(r = 10))
  tiff('Detritivore_Autumn.tiff',units="in",width=8,height=5,bg="white",res=500)
  d
  dev.off()
  
  emmeans_results <- emmeans(models2[["Detritivore"]][["Summer"]], ~ Movement.pattern|Treatment)
  
  # Apply pairwise contrasts with Sidak adjustment: To identify specific differences between levels of categorical predictors.
  contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
  summary(contrast_results)
  
  ############################################################################################
  #To support results from CANOCO5, we can calculate FD Rao as rao.diversity
  
  dataset6 <- as.data.frame(dataset6)
  
  # Remove leading/trailing spaces from Functional.group names
  dataset6$Functional.group <- trimws(dataset6$Functional.group)
  
  # Summarize abundance per functional group per site
  abundance_matrix <- dataset6 %>%
    group_by(Trap, Treatment, Movement.pattern, Season, Functional.group) %>%
    summarise(TotalNumber = sum(Number), .groups = "drop") %>%
    pivot_wider(names_from = Functional.group, values_from = TotalNumber, values_fill = 0)
  
  # Convert to matrix format (remove metadata columns)
  abundance_data <- as.matrix(abundance_matrix[, -c(1:4)])
  rownames(abundance_data) <- paste(abundance_matrix$Trap, 
                                    abundance_matrix$Treatment, 
                                    abundance_matrix$Movement.pattern, 
                                    abundance_matrix$Season, sep = "_")
  
  # Create a trait matrix (Functional groups as categorical variables)
  functional_groups <- sort(unique(dataset6$Functional.group))  # Unique functional groups
  trait_matrix <- diag(length(functional_groups))  # Identity matrix
  rownames(trait_matrix) <- functional_groups
  colnames(trait_matrix) <- functional_groups
  
  # Ensure column names in abundance_data match trait_matrix row names
  matching_groups <- intersect(colnames(abundance_data), rownames(trait_matrix))
  abundance_data <- abundance_data[, matching_groups, drop = FALSE]
  trait_matrix <- trait_matrix[matching_groups, matching_groups, drop = FALSE]
  
  # Compute Rao's Quadratic Entropy (RaoQ) using rao.diversity()
  raoQ_results <- rao.diversity(abundance_data, traits = trait_matrix)
  
  # Extract RaoQ values and format results
  raoQ_values <- data.frame(
    Site = names(raoQ_results$FunRao),  # Extract site names
    RaoQ = raoQ_results$FunRao  # Extract computed RaoQ values
  )
  
  # Split 'Site' into original grouping variables
  raoQ_values <- raoQ_values %>%
    separate(Site, into = c("Trap", "Treatment", "Movement.pattern", "Season"), sep = "_", convert = TRUE)
  
  # Convert Trap to character in both datasets for merging
  abundance_matrix$Trap <- as.character(abundance_matrix$Trap)
  raoQ_values$Trap <- as.character(raoQ_values$Trap)
  
  # Merge RaoQ values with metadata from abundance_matrix
  final_results <- left_join(abundance_matrix[, 1:4], raoQ_values, by = c("Trap", "Treatment", "Movement.pattern", "Season"))
  
  # Visualization: Boxplot of RaoQ Across Treatments and Movement Patterns
  ggplot(final_results, aes(x = Treatment, y = RaoQ, fill = Movement.pattern)) +
    geom_boxplot() +
    facet_wrap(~ Season) +
    theme_minimal() +
    labs(title = "Rao's Quadratic Entropy across Treatments",
         y = "RaoQ", x = "Treatment")
  
  # Statistical Analysis: Mixed-Effects Model to Test RaoQ Differences
  fit_raoQ <- lm(RaoQ ~ Treatment*Movement.pattern, data = final_results)
  summary(fit_raoQ)
  
  # Post-Hoc Pairwise Comparisons (Tukey's Test)
  emmeans_raoQ <- emmeans(fit_raoQ, pairwise ~ Treatment * Movement.pattern, adjust = "tukey")
  print(emmeans_raoQ$contrasts)
  
  # Save results to CSV file (optional)
  write.csv(final_results, "RaoQ_results.csv", row.names = FALSE)
  
  emmeans_results <- emmeans(fit_raoQ, ~ Movement.pattern | Treatment)
  
  # Apply pairwise contrasts with Sidak adjustment
  contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
  summary(contrast_results)
  
  # Compute the dissimilarity matrix between functional groups
  dissimilarity_matrix <- daisy(trait_matrix, metric = "gower")
  
  # Convert to dataframe
  dissim_df <- as.data.frame(as.matrix(dissimilarity_matrix))
  dissim_df$FunctionalGroup <- rownames(dissim_df)
  
  # Identify which functional groups have the highest dissimilarity
  dissim_long <- dissim_df %>%
    pivot_longer(cols = -FunctionalGroup, names_to = "ComparisonGroup", values_to = "Dissimilarity") %>%
    arrange(desc(Dissimilarity))
  
  # Filter top 10 most dissimilar pairs
  top_dissimilarities <- dissim_long %>% head(10)
  top_dissimilarities
  
  # Merge RaoQ results with functional groups
  raoQ_with_groups <- dataset6 %>%
    group_by(Functional.group) %>%
    summarise(Mean_RaoQ = mean(raoQ_values$RaoQ, na.rm = TRUE),
              SD_RaoQ = sd(raoQ_values$RaoQ, na.rm = TRUE),
              N = n()) %>%
    arrange(desc(Mean_RaoQ))
  
  # Print the summary
  print(raoQ_with_groups)
  #######################################################################################################

  
  