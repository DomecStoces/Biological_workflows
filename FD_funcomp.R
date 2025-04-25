library(FD)
library(glmmTMB)
library(emmeans)
library(dplyr)
library(lme4)
library(tidyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(adiv)

dataset6 <- as.data.frame(dataset6)

# Ensure there are no leading/trailing spaces in Functional.group names
dataset6$Functional.group <- trimws(dataset6$Functional.group)

# Summarize abundance per functional group per site
abundance_matrix <- dataset6 %>%
  group_by(Trap, Treatment, Movement.pattern, Season, Functional.group) %>%
  summarise(TotalNumber = sum(Number), .groups = "drop") %>%
  pivot_wider(names_from = Functional.group, values_from = TotalNumber, values_fill = 0)

# Convert the abundance matrix to the correct format
abundance_data <- as.matrix(abundance_matrix[, -c(1:4)])  # Remove non-numeric columns
rownames(abundance_data) <- paste(abundance_matrix$Trap, 
                                  abundance_matrix$Treatment, 
                                  abundance_matrix$Movement.pattern, 
                                  abundance_matrix$Season, sep = "_")

# Create trait matrix (Identity matrix since Functional Groups are categorical)
functional_groups <- sort(unique(dataset6$Functional.group))  # Sort to match order
trait_matrix <- diag(length(functional_groups))
rownames(trait_matrix) <- functional_groups
colnames(trait_matrix) <- functional_groups

# Ensure column names in abundance_data match row names in trait_matrix
matching_groups <- intersect(colnames(abundance_data), rownames(trait_matrix))

# Filter both matrices to only include matching groups
abundance_data <- abundance_data[, matching_groups, drop = FALSE]
trait_matrix <- trait_matrix[matching_groups, matching_groups, drop = FALSE]

# Verify if dimensions match before running functcomp
print(dim(trait_matrix))   # Should be (n functional groups × n functional groups)
print(dim(abundance_data)) # Should be (n sites × n functional groups)

# Calculate functional composition
functional_composition <- functcomp(trait_matrix, abundance_data, CWM.type = "all")

# Convert results to a dataframe and merge with metadata
results <- bind_cols(abundance_matrix[, 1:4], as.data.frame(functional_composition))

# View results
print(results)

fit_predator <- lmer(Predator_1 ~ Treatment * Movement.pattern + (1|Trap), data = results)
fit_herbivore <- lmer(Herbivore_1 ~ Treatment * Movement.pattern + (1|Trap), data = results)
fit_omnivore <- lmer(Omnivore_1 ~ Treatment * Movement.pattern + (1|Trap), data = results)
fit_detritivore <- lmer(Detritivore_1 ~ Treatment * Movement.pattern + (1|Trap), data = results)
fit_saproxylic <- lmer(Saproxylic_1 ~ Treatment * Movement.pattern + (1|Trap), data = results)

#Post-hoc Pairwise comparisons with applied Sidak correction for multiple comparisons to control the family-wise error rate
# Compute estimated marginal means
emmeans_results <- emmeans(fit_saproxylic, ~ Movement.pattern | Treatment)

# Apply pairwise contrasts with Sidak adjustment
contrast_results <- contrast(emmeans_results, method = "pairwise", adjust = "sidak")
summary(contrast_results)