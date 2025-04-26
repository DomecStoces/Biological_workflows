#To support results from CANOCO5, we can calculate FD Rao as rao.diversity

library(SYNCSA)
library(dplyr)
library(tidyr)
library(emmeans)

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