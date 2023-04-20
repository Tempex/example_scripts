

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Start the analysis for group 3.2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# The main research question of this group is whether participants are able
# to remember words better if they have been produced by saying them out loud
# compared to only hearing them and whether such an effect would interact with
# a word's semantic context. Eventually, the group hypothesized that a
# production is especially beneficial for remembering words that are
# not semantically associated with any other words that need to be learned.

# Delete all variables from the environment
rm(list = ls() )

# Add toolboxes
require(tidyverse)
require(ggplot2)
require(ez)
require(readxl)
require(psycho)

# Define a working directory (here we only work in file directories in which
# the R script was placed)
wd <- ""


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Read in data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Read in the data of the No-Production Condition
df_horen <- 
  read_excel("data/23-02-08_fopra-group32_logfile-hören.xlsx")

# Read in the data of the Production Condition
df_sagen <- 
  read_excel("data/23-02-08_fopra-group32_logfile-sagen.xlsx")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Clean data columns
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Extract relevant columns from the No-Production Condition
df_horen <- df_horen %>%
  select("jatosStudyResultId", "Gelernt", "Semantik", "correct_response",
         "response") %>%
  mutate(prod = 0)

# Extract relevant columns from the Production Condition
df_sagen <- df_sagen %>%
  select("jatosStudyResultId", "Gelernt", "Semantik", "correct_response",
         "response") %>%
  mutate(prod = 1)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Combine datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Combine both tables
df <- rbind(df_horen, df_sagen)

# Relabel the columns
colnames(df) <- c("vpn", "learn", "sem", "corr", "resp", "prod")

# Rearrange the columns
df <- df %>% select("vpn", "prod", "sem", "learn", "resp", "corr")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Clean data rows
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Exclude rows from the calculation task
df <- df %>% filter(resp == "n" | resp == "a")

# Exclude participants without a participant number
df <- df %>% filter(!is.na(vpn))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Calculate Sensitivity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Determine Hits, False Alarms, Correct Rejections and Misses
df <- df %>%
  mutate(hit = if_else(corr == "a" & resp == "a",1,0),
         mis = if_else(corr == "a" & resp == "n",1,0),
         fa  = if_else(corr == "n" & resp == "a",1,0),
         cr  = if_else(corr == "n" & resp == "n",1,0))

# Drop columns that have become non-relevant
df <- df %>% select(-learn, -corr, -resp)

# Calculate the number of trials in each report category
df <- df %>%
  group_by(vpn, prod, sem) %>%
  summarise(n_hit = sum(hit),
            n_mis = sum(mis),
            n_fa = sum(fa),
            n_cr = sum(cr))

# Correct the proportions according to Hautus (1995)
df$hr <- (df$n_hit + .5) / ((df$n_hit + df$n_mis) + 1)
df$fr <- (df$n_fa  + .5) / ((df$n_fa  + df$n_cr ) + 1)

# Calculate the sensitivity (with correction)
df$dc <- round(qnorm(df$hr) - qnorm(df$fr),3)

# Alternatively we can use the following formula
dprime <- dprime(
  
  n_hit = df$n_hit,
  n_fa = df$n_fa,
  n_miss = df$n_mis,
  n_cr = df$n_cr,
  n_targets = 18,
  n_distractors = 18,
  adjusted = TRUE
  
)

# Extract the d values from the function results
df$dcf <- dprime$dprime


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Exclude participants
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Exclude participants that are too bad in both of their conditions
df <- df %>% group_by(vpn) %>% mutate(excl = sum(dc < .5)) %>% filter(excl < 2)

# Remove columns that have become irrelevant from the dataset
df <- df %>% select(vpn, prod, sem, dc)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Perform the ANOVA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Define factors
df$vpn <- as.factor(df$vpn)
df$prod <- as.factor(df$prod)
df$sem <- as.factor(df$sem)

# Rename variable expressions
levels(df$prod) <- c("no", "yes")
levels(df$sem) <- c("no", "yes")

# Perform the ANOVA (there are 31 vs 34 participants in each production group)
ezANOVA(
  
  data = df,
  dv = dc,
  wid = vpn,
  within = .(sem),
  between = .(prod),
  type = 3,
  detailed = TRUE
  
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot the results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# Create a plot
df_plot <- ezPlot(
  
  data = df,
  dv = .(dc),
  wid =.(vpn),
  within = .(sem),
  between = .(prod),
  type = 3,
  x = prod,
  split = sem,
  do_lines = TRUE,
  y_lab = "sensitivity",
  x_lab = "production",
  split_lab = "semantic"
  
) +

theme(
  
  legend.position = "top",
  panel.border = element_rect(color = "black",  fill = NA)
  
)

# Save the plot
ggsave(file = "images/group32_plot_anova_prod-sem.jpg", plot = df_plot,
       device = "jpeg", width = 5, height = 6, dpi = 250)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# End the analysis for group 3.2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

