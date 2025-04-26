#indicspecies - multipatt: Multi-level pattern analysis
# Aggregate data
dataset7 <- dataset7 %>%
  group_by(Trap, Species, Movement.pattern, Season, Treatment) %>%
  summarize(Number = sum(Number, na.rm = TRUE), .groups = "drop")

# Create species-by-site matrix
species_matrix <- dataset7 %>%
  pivot_wider(
    names_from = Species,
    values_from = Number,
    values_fill = list(Number = 0)
  ) %>%
  as.data.frame() %>%
  select(where(is.numeric)) %>%
  as.matrix()

# Create grouping factor
grouping_factor <- dataset7 %>%
  distinct(Trap, Movement.pattern, Season, Treatment) %>%
  mutate(Group = interaction(Movement.pattern, Treatment)) %>%
  arrange(Trap, Season) %>%
  pull(Group)

# Debugging
stopifnot(nrow(species_matrix) == length(grouping_factor))

# Multi-level pattern analysis
set.seed(123)
results <- multipatt(species_matrix, grouping_factor, func = "r.g", control = how(nperm = 999))

#func="r.g." - Pearson's r calculates the correlation coefficient between species presence/absence (or abundance) and the grouping factor: works well for detecting correlations between species and gradients.

# Extract raw p-values
raw_p_values <- results$sign$p.value

# Adjust p-values using the desired method (e.g., "fdr", "holm", "bonferroni")
adjusted_p_values <- p.adjust(raw_p_values, method = "fdr")

# Replace raw p-values with adjusted p-values in the results object
results$sign$p.value.adj <- adjusted_p_values

# Summary with adjusted p-values
summary_results <- summary(results)
