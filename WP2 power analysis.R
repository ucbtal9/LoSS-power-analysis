# WP2 Power Analysis

#   1. Estimate expected ground reaction forces for selecting an
#      appropriate force plate.
#   2. Establish null allometric relationships from published
#      locomotor data (Astley 2016).
#   3. Perform simulation-based power analyses to optimise the
#      number of individuals and species required for the study.

#-------------------------------------------------------------------------------
# 0. Setting up
#-------------------------------------------------------------------------------
rm(list = ls())


# Load  packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom)
library(performance)
library(DHARMa)
library(patchwork)
library(picante)

#Import study taxa tree
tree <- read.nexus("LOTSS_tree.nex")
data <- read.csv("Study taxa.csv", stringsAsFactors = FALSE)

# Keep only taxa present in tree
study_taxa <- data %>%
  filter(
    Species %in% tree$tip.label
  )

# Prune to only those who will have biomechanical data
study_taxa <- study_taxa %>%
  filter(!is.na(Site))

# Keep field taxa only for the main field sampling analysis
study_taxa <- study_taxa %>%
  filter(Site != "LITERATURE")

# Prune tree to field taxa
study_tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, study_taxa$Species)
)

# Check final field tree
plot(study_tree,cex = 0.7)


# Field + literature taxa for H3 evolutionary analysis
comparative_taxa <- data %>%
  filter(
    Species %in% tree$tip.label,
    !is.na(Site)
  )

# Prune tree to field + literature taxa
literature_tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, comparative_taxa$Species)
)

# Check field + literature tree
plot(literature_tree,cex = 0.7)

#Import Astley (2016) supplementary data
astley <- read_csv("Astley 2016 - TableS2.csv")

# Define key variables
astley$species <- as.factor(astley$species)
mass_var <- "mass_g"


### Number of observations for each variable ----

data_summary <- tibble(
  Variable = names(astley),
  N = sapply(astley, function(x) sum(!is.na(x))),
  Missing = sapply(astley, function(x) sum(is.na(x))),
  Percent_Missing = round(
    100 * sapply(astley, function(x) mean(is.na(x))), 1)
)

print(data_summary)


### Histogram of body masses ----

ggplot(astley, aes(x = mass_g)) +
  geom_histogram(bins = 20) +
  labs(x = "Body mass (g)",
       y = "Number of individuals") +
  theme_classic()


### Log-transformed body masses ----

ggplot(astley, aes(x = log10(mass_g))) +
  geom_histogram(bins = 20) +
  labs(x = expression(log[10]*" Body mass (g)"),
       y = "Number of individuals") +
  theme_classic()


### Number of individuals per species ----

astley %>%
  count(species) %>%
  ggplot(aes(x = reorder(species, n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(x = "Species",
       y = "Number of individuals") +
  theme_classic()

#-------------------------------------------------------------------------------
# 1. Estimate the expected peak ground reaction force of a miniaturised frog under 
# the null hypothesis that it follows the same scaling relationship as the frogs 
# in Astley (2016)
#-------------------------------------------------------------------------------

# Do miniature frogs follow the established scaling relationship for jump force?
# We first test whether relative force (body weights) is independent of body mass,
# then estimate absolute force (mN) expected at miniature body sizes.

#-------------------------------------------------------------------------------
# 1.1 Test whether relative peak jump force is independent of body mass
#-------------------------------------------------------------------------------
jump_data <- astley %>%
  filter(!is.na(Peak_Jump_Force_BW),
         !is.na(mass_g))

force_bw_model <- lm(log10(Peak_Jump_Force_BW) ~ log10(mass_g),
                     data = jump_data)
summary(force_bw_model)

force_bw_results <- broom::tidy(force_bw_model,
                                conf.int = TRUE)
write.csv(force_bw_results,
          "Relative_Force_Scaling_Coefficients.csv",
          row.names = FALSE)

p1 <- ggplot(jump_data,
             aes(x = mass_g,
                 y = Peak_Jump_Force_BW)) +
  geom_point(size = 2) +
  stat_smooth(method = "lm",
              colour = "black",
              fill = "grey80") +
  scale_x_log10() +
  labs(x = "Body mass (g)",
       y = "Peak jump force (body weights)") +
  theme_classic()

p1


ggsave("Figure1_PeakForce_BW.pdf",
       p1,
       width = 6,
       height = 5)

#-------------------------------------------------------------------------------
# 1.2 Convert peak force from body weights to Newtons and milliNewtons
#-------------------------------------------------------------------------------
g <- 9.81

astley <- astley %>%
  mutate(
    Peak_Jump_Force_N = Peak_Jump_Force_BW * (mass_g / 1000) * g,
    Peak_Jump_Force_mN = Peak_Jump_Force_N * 1000
  )

jump_data <- astley %>%
  filter(
    !is.na(Peak_Jump_Force_BW),
    !is.na(Peak_Jump_Force_mN),
    !is.na(mass_g)
  )

#-------------------------------------------------------------------------------
# 1.3 Fit absolute force allometry
#-------------------------------------------------------------------------------
# Model: log10(force) ~ log10(body mass)
force_model <- lm(log10(Peak_Jump_Force_mN) ~ log10(mass_g),
                  data = jump_data)
summary(force_model)

force_results <- broom::tidy(force_model,
                             conf.int = TRUE)

write.csv(force_results,
          "Absolute_Force_Scaling_Coefficients.csv",
          row.names = FALSE)

#-------------------------------------------------------------------------------
# 1.4 Predict expected force for miniature frogs
#-------------------------------------------------------------------------------
prediction_data <- tibble(mass_mg = seq(5,100,1)) %>%
  mutate(mass_g = mass_mg/1000)

pred <- predict(force_model,
                newdata = prediction_data,
                interval = "prediction")

prediction_data <- bind_cols(prediction_data,
                             as.data.frame(pred)) %>%
  mutate(Predicted_mN = 10^fit,
         Lower_PI_mN = 10^lwr,
         Upper_PI_mN = 10^upr)

mini_frog_predictions <- prediction_data %>%
  filter(mass_mg %in% c(5,10,20,50,100))
print(mini_frog_predictions)

write.csv(prediction_data,
          "Predicted_Jump_Forces.csv",
          row.names = FALSE)

write.csv(mini_frog_predictions,
          "MiniFrog_Expected_Forces.csv",
          row.names = FALSE)

#-------------------------------------------------------------------------------
# 1.5 Prepare predicted force scaling relationship for combined figure
#-------------------------------------------------------------------------------
# Prediction grid in grams
prediction_plot_data <- tibble(
  mass_g = seq(
    0.01,
    100,
    length.out = 500
  )
)

# Generate predictions
pred <- predict(
  force_model,
  newdata = prediction_plot_data,
  interval = "prediction"
)

prediction_plot_data <- bind_cols(
  prediction_plot_data,
  as.data.frame(pred)
) %>%
  mutate(
    Predicted_mN = 10^fit,
    Lower_PI_mN = 10^lwr,
    Upper_PI_mN = 10^upr
  )

# Observed Astley data
observed_force_plot <- jump_data %>%
  select(
    mass_g,
    Peak_Jump_Force_mN
  )

# Predicted 10 mg miniature frog
mini_frog_force_prediction <- prediction_plot_data %>%
  slice(
    which.min(abs(mass_g - 0.01))
  )



#-------------------------------------------------------------------------------
# Section 2 - Define expected scaling relationships for locomotor performance
#-------------------------------------------------------------------------------

# Aim:
# Quantify expected scaling relationships across locomotor modes using
# published comparative data (Astley 2016).
#
# Model:
# log10(Trait) ~ log10(Body mass)
#
# These relationships define the null expectation that miniaturised frogs
# follow established locomotor scaling rules.

#-------------------------------------------------------------------------------
# 2.1 Define primary locomotor traits
#-------------------------------------------------------------------------------
traits <- c(
  # Jumping
  "Peak_Jump_Force_mN",
  "Jump_Takeoff_Velocity_SVL.second",
  
  # Swimming
  "Mean_Swim_Speed_SVL.second",
  
  # Walking
  "Walking_Speed_SVL.second"
)

#-------------------------------------------------------------------------------
# 2.2 Function to fit allometric relationships
#-------------------------------------------------------------------------------

fit_allometry <- function(data, trait, mass = "mass_g") {
  
  
  temp <- data %>%
    select(
      all_of(c(mass, trait))
    ) %>%
    filter(
      !is.na(.data[[mass]]),
      !is.na(.data[[trait]]),
      .data[[mass]] > 0,
      .data[[trait]] > 0
    ) %>%
    mutate(
      log_mass = log10(.data[[mass]]),
      log_trait = log10(.data[[trait]])
    )
  
  
  model <- lm(
    log_trait ~ log_mass,
    data = temp
  )
  
  
  model_summary <- broom::tidy(
    model,
    conf.int = TRUE
  )
  
  
  slope <- model_summary %>%
    filter(
      term == "log_mass"
    )
  
  
  intercept <- model_summary %>%
    filter(
      term == "(Intercept)"
    )
  
  
  tibble(
    
    Trait = trait,
    
    N = nrow(temp),
    
    Slope = slope$estimate,
    
    Slope_lower95 = slope$conf.low,
    
    Slope_upper95 = slope$conf.high,
    
    Intercept = intercept$estimate,
    
    R2 = summary(model)$r.squared,
    
    Residual_SD = summary(model)$sigma
    
  )
  
}

#-------------------------------------------------------------------------------
# 2.3 Run models
#-------------------------------------------------------------------------------

allometry_results <- map_dfr(
  traits,
  ~fit_allometry(
    astley,
    .x
  )
)
print(allometry_results)

# Save results
write.csv(
  allometry_results,
  "Allometric_relationships_summary.csv",
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 2.4 Generate prediction plots for each locomotor trait
#-------------------------------------------------------------------------------

plot_prediction <- function(data, trait, mass = "mass_g") {
  
  observed_plot <- data %>%
    select(
      mass_g = all_of(mass),
      Trait_value = all_of(trait)
    ) %>%
    filter(
      !is.na(mass_g),
      !is.na(Trait_value),
      mass_g > 0,
      Trait_value > 0
    )
  
  model_data <- observed_plot %>%
    mutate(
      log_mass = log10(mass_g),
      log_trait = log10(Trait_value)
    )
  
  model <- lm(
    log_trait ~ log_mass,
    data = model_data
  )
  
  model_summary <- broom::tidy(model)
  
  slope <- model_summary %>%
    filter(term == "log_mass") %>%
    pull(estimate)
  
  p_value <- model_summary %>%
    filter(term == "log_mass") %>%
    pull(p.value)
  
  prediction_plot_data <- tibble(
    mass_g = seq(
      0.01,
      100,
      length.out = 500
    )
  ) %>%
    mutate(
      log_mass = log10(mass_g)
    )
  
  pred <- predict(
    model,
    newdata = prediction_plot_data,
    interval = "prediction"
  )
  
  prediction_plot_data <- bind_cols(
    prediction_plot_data,
    as.data.frame(pred)
  ) %>%
    mutate(
      Predicted = 10^fit,
      Lower_PI = 10^lwr,
      Upper_PI = 10^upr
    )
  
  mini_frog_prediction <- prediction_plot_data %>%
    slice(
      which.min(abs(mass_g - 0.01))
    )
  
  trait_label <- case_when(
    trait == "Peak_Jump_Force_mN" ~ "Peak jump force (mN)",
    trait == "Jump_Takeoff_Velocity_SVL.second" ~ "Take-off velocity (SVL s-1)",
    trait == "Mean_Swim_Speed_SVL.second" ~ "Mean swimming speed (SVL s-1)",
    trait == "Walking_Speed_SVL.second" ~ "Walking speed (SVL s-1)",
    TRUE ~ trait
  )
  
  stats_label <- paste0(
    "Slope = ",
    round(slope,2),
    "\nP = ",
    format.pval(p_value, digits = 2),
    "\nN = ",
    nrow(observed_plot)
  )
  
  # Position statistics differently for positive relationship
  
  if(trait == "Peak_Jump_Force_mN"){
    
    stats_x <- 10
    stats_y <- min(prediction_plot_data$Lower_PI) * 15
    
  } else {
    
    stats_x <- 10
    stats_y <- max(prediction_plot_data$Upper_PI)
    
  }
  
  p <- ggplot() +
    
    geom_ribbon(
      data = prediction_plot_data,
      aes(
        x = mass_g,
        ymin = Lower_PI,
        ymax = Upper_PI,
        fill = "95% prediction interval"
      ),
      alpha = 0.2
    ) +
    
    geom_point(
      data = observed_plot,
      aes(
        x = mass_g,
        y = Trait_value,
        colour = "Astley (2016) observations"
      ),
      size = 2
    ) +
    
    geom_line(
      data = prediction_plot_data,
      aes(
        x = mass_g,
        y = Predicted,
        colour = "Predicted scaling relationship"
      ),
      linewidth = 1.1
    ) +
    
    geom_point(
      data = mini_frog_prediction,
      aes(
        x = mass_g,
        y = Predicted
      ),
      shape = 21,
      size = 4,
      fill = "white",
      colour = "black"
    ) +
    
    annotate(
      "text",
      x = stats_x,
      y = stats_y,
      label = stats_label,
      hjust = 0,
      vjust = 1,
      size = 3.2
    ) +
    
    scale_x_log10(
      breaks = c(
        0.01,
        0.1,
        1,
        10,
        100
      ),
      labels = c(
        "0.01",
        "0.1",
        "1",
        "10",
        "100"
      )
    ) +
    
    scale_y_log10(
      labels = scales::label_number()
    ) +
    
    scale_colour_manual(
      values = c(
        "Astley (2016) observations" = "black",
        "Predicted scaling relationship" = "steelblue"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "95% prediction interval" = "steelblue"
      )
    ) +
    
    labs(
      x = "Body mass (g)",
      y = trait_label,
      colour = "",
      fill = ""
    ) +
    
    theme_classic() +
    
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  return(p)
  
}

#-------------------------------------------------------------------------------
# 2.5 Generate combined Figure 3
#-------------------------------------------------------------------------------
p_force <- plot_prediction(
  astley,
  "Peak_Jump_Force_mN"
)

p_velocity <- plot_prediction(
  astley,
  "Jump_Takeoff_Velocity_SVL.second"
)

p_swim <- plot_prediction(
  astley,
  "Mean_Swim_Speed_SVL.second"
)

p_walk <- plot_prediction(
  astley,
  "Walking_Speed_SVL.second"
)

figure3 <- (
  p_force +
    p_velocity +
    p_swim +
    p_walk
) +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
figure3

ggsave("Figure3_Locomotor_Allometry_Predictions.pdf",
  figure3,width = 10,height = 8)



#-------------------------------------------------------------------------------
# Section 3 - Simulation-based sampling sensitivity analysis
#-------------------------------------------------------------------------------

# Aim:
# Determine the minimum number of individuals per species and species required
# to recover the expected allometric relationships.
#
# IMPORTANT:
# All biological parameters are derived from real data:
#
# 1. Astley (2016) provides measured SVL and body mass.
# 2. The Astley SVL-mass relationship is used to estimate mass for study taxa.
# 3. Sections 1-2 provide the empirical locomotor allometries.
# 4. Residual variation around those relationships is taken directly from Astley.
# 5. Study taxa are the actual candidate species in Study taxa.csv.
#
# No arbitrary biological effect sizes are introduced.

#-------------------------------------------------------------------------------
# 3.1 Estimate body mass from SVL
#-------------------------------------------------------------------------------

# Astley contains measured body mass (g) and SVL (mm).
# We fit the empirical allometric relationship:
#
# log10(mass) ~ log10(SVL)
#
# This allows body mass to be estimated for study species for which mean SVL
# is available but body mass has not yet been measured.

mass_model <- lm(log10(mass_g) ~ log10(SVL_mm),data=astley)
summary(mass_model)

# Estimate body mass for every study species.
# SVL is the published/known mean SVL in Study taxa.csv.
study_taxa <- study_taxa %>%
  mutate(
    Estimated_mass_g=10^predict(
      mass_model,
      newdata=tibble(SVL_mm=SVL)
    )
  )

# Check estimated masses
study_taxa[order(study_taxa$Estimated_mass_g),c("Species","SVL","Estimated_mass_g")]

#-------------------------------------------------------------------------------
# 3.2 Fit null locomotor models from empirical data
#-------------------------------------------------------------------------------
# These are the same null allometric relationships established in Sections 1-2.
# Each model uses the real Astley observations.
#
# The residual SD from each model represents the observed individual-level
# variation around the expected allometric relationship.

trait_models <- list()
for(trait in traits){
  
  temp <- astley %>%
    select(mass_g,all_of(trait)) %>%
    filter(!is.na(mass_g),!is.na(.data[[trait]]),mass_g>0,.data[[trait]]>0) %>%
    mutate(log_mass=log10(mass_g),log_trait=log10(.data[[trait]]))
  
  trait_models[[trait]] <- lm(log_trait~log_mass,data=temp)
  
}

trait_parameters <- tibble(
  Trait=traits,
  Intercept=map_dbl(trait_models,~coef(.x)[1]),
  Slope=map_dbl(trait_models,~coef(.x)[2]),
  Residual_SD=map_dbl(trait_models,sigma)
)
print(trait_parameters)

write.csv(
  trait_parameters,
  "Simulation_parameters.csv",
  row.names=FALSE
)

#-------------------------------------------------------------------------------
# 3.3 Predict expected locomotor performance for study species
#-------------------------------------------------------------------------------

# Predictions are based on:
# 1. estimated body mass from Astley SVL-mass relationship
# 2. empirical locomotor allometries from Astley
#
# The models are on log10 scale, so predictions are back-transformed.

study_taxa <- study_taxa %>%
  mutate(log_mass=log10(Estimated_mass_g))

for(trait in names(trait_models)){
  
  study_taxa[[paste0("Expected_",trait)]] <- 10^predict(
    trait_models[[trait]],
    newdata=tibble(log_mass=study_taxa$log_mass)
  )
  
}

study_taxa %>%
  select(
    Species,
    SVL,
    Estimated_mass_g,
    starts_with("Expected_")
  )

#-------------------------------------------------------------------------------
# 3.4 Simulate individual sampling
#-------------------------------------------------------------------------------
# Individuals within each species are simulated around that species' expected
# locomotor value.
#
# The amount of individual variation is NOT arbitrary: it is the residual SD
# estimated from the corresponding Astley allometric model.
#
# We therefore ask:
# "How much does increasing the number of individuals improve the precision
#  with which a species mean is estimated?"

#-------------------------------------------------------------------------------
# Simulate individuals around Astley-predicted species mean
#-------------------------------------------------------------------------------

simulate_individuals <- function(model,mass_g,n_ind){
  
  expected_log <- predict(
    model,
    newdata=tibble(log_mass=log10(mass_g))
  )
  
  residual_sd <- sigma(model)
  
  10^rnorm(
    n=n_ind,
    mean=expected_log,
    sd=residual_sd
  )
  
}

individual_numbers <- c(3,5,8,10,15)
n_sim <- 1000
set.seed(123)

# Use actual study-species body masses rather than invented masses.
# Only species with an available field site are included here because these
# are the taxa that could realistically contribute to the proposed dataset.

field_species <- study_taxa %>%
  filter(!is.na(Site),Site!="NA")

individual_sensitivity <- map_dfr(
  individual_numbers,
  function(n_ind){
    
    map_dfr(
      1:n_sim,
      function(i){
        
        sampled_species <- field_species %>%
          slice_sample(n=1)
        
        mass <- sampled_species$Estimated_mass_g
        
        values <- simulate_individuals(
          trait_models$Peak_Jump_Force_mN,
          mass,
          n_ind
        )
        
        tibble(
          Individuals=n_ind,
          Species=sampled_species$Species,
          Estimated_mean=mean(log10(values)),
          Estimated_SD=sd(log10(values))
        )
        
      }
      
    )
    
  }
)

# Summarise precision of species-level estimates.
individual_summary <- individual_sensitivity %>%
  group_by(Individuals) %>%
  summarise(
    Mean_SD=mean(Estimated_SD),
    SD_precision=sd(Estimated_SD),
    .groups="drop"
  )
print(individual_summary)

individual_summary <- individual_summary %>%
  mutate(
    Percent_improvement=round(
      100*(first(SD_precision)-SD_precision)/first(SD_precision),
      1
    )
  )

print(individual_summary)


#-------------------------------------------------------------------------------
# Section 3.5 - Simulation-based power and sensitivity analysis
#-------------------------------------------------------------------------------
# Aim:
#
# The evolutionary analysis will be conducted at the species level using mean
# locomotor performance calculated from multiple individuals.
#
# The purpose of this analysis is therefore to determine how many species and
# individuals per species are required to detect departures of study-species
# locomotor performance from the empirical Astley (2016) relationship.
#
# H0: study-species performance follows the Astley-predicted expectation.
# H1: study-species performance departs from the Astley-predicted expectation.
#
# No biologically justified alternative effect size currently exists for
# miniaturised frog locomotor performance. Therefore, ±20%, ±50% and ±80%
# departures are used only as sensitivity scenarios. They are not presented as
# biologically established effect sizes.
#
# The analysis uses:
# - empirical Astley locomotor models;
# - empirical Astley residual variation;
# - body masses estimated for the actual study taxa in Section 3.1;
# - the actual study species available for field sampling;
# - species means calculated from simulated individuals.
#
# The statistical test is performed on species-level departures from the
# Astley prediction, matching the level of the subsequent evolutionary analysis.
#-------------------------------------------------------------------------------

# Define sampling scenarios
species_numbers <- 4:10
individual_numbers <- c(5, 8, 10)
n_sim <- 1000
alpha <- 0.05

# Sensitivity scenarios for departure from the Astley-predicted performance.

# These are sensitivity scenarios, not biologically established effect sizes.
deviation_levels <- c(0.00, 0.20, 0.50, 0.80)
directions <- c("Negative", "Positive")

set.seed(123)


# Identify species available for the proposed study
field_species <- study_taxa %>%
  filter(!is.na(Site), Site != "NA", !is.na(Estimated_mass_g), Estimated_mass_g > 0) %>%
  distinct(Species, .keep_all = TRUE)

print(field_species %>% select(Species, SVL, Estimated_mass_g))

if(max(species_numbers) > nrow(field_species)){
  stop(paste0("The requested maximum number of species (", max(species_numbers),
              ") exceeds the number of available study species (", nrow(field_species), ")."))
}


# Function to simulate species-level performance under a departure from Astley
simulate_deviation_dataset <- function(model, study_species, n_species, n_ind, deviation, direction){
  
  # Randomly select the actual study species.
  
  sampled_species <- study_species %>%
    slice_sample(n = n_species, replace = FALSE)
  
  # Empirical Astley residual variation.
  
  residual_sd <- sigma(model)
  
  # Calculate Astley-predicted performance for each sampled species.
  
  sampled_species <- sampled_species %>%
    mutate(
      log_mass = log10(Estimated_mass_g),
      null_log_trait = as.numeric(predict(model, newdata = tibble(log_mass = log_mass))),
      null_trait = 10^null_log_trait
    )
  
  # Apply the sensitivity departure directly to predicted species performance.
  
  sampled_species <- sampled_species %>%
    mutate(
      alternative_trait = case_when(
        deviation == 0 ~ null_trait,
        direction == "Positive" ~ null_trait * (1 + deviation),
        direction == "Negative" ~ null_trait * (1 - deviation),
        TRUE ~ NA_real_
      ),
      alternative_log_trait = log10(alternative_trait)
    )
  
  # Simulate individuals around the departed species-level expectation.
  
  simulated_individuals <- map_dfr(seq_len(nrow(sampled_species)), function(i){
    
    simulated_log_trait <- rnorm(
      n = n_ind,
      mean = sampled_species$alternative_log_trait[i],
      sd = residual_sd
    )
    
    tibble(
      Species = sampled_species$Species[i],
      mass_g = sampled_species$Estimated_mass_g[i],
      log_mass = sampled_species$log_mass[i],
      null_trait = sampled_species$null_trait[i],
      alternative_trait = sampled_species$alternative_trait[i],
      log_trait = simulated_log_trait
    )

    
  })
  
  # Calculate species-level mean performance, matching the evolutionary analysis.
  
  simulated_species_means <- simulated_individuals %>%
    group_by(Species, mass_g, log_mass, null_trait, alternative_trait) %>%
    summarise(
      log_trait = mean(log_trait),
      trait = 10^log_trait,
      .groups = "drop"
    ) %>%
    mutate(
      departure = (trait - null_trait) / null_trait,
      log_departure = log10(trait / null_trait)
    )
  
  simulated_species_means
}

# Primary model
primary_model <- trait_models$Peak_Jump_Force_mN

print(paste("Astley null slope:", round(coef(primary_model)[2], 4)))

# Simulation-based power analysis for the primary trait
power_results <- map_dfr(species_numbers, function(n_species){
  map_dfr(individual_numbers, function(n_ind){
    map_dfr(deviation_levels, function(deviation){
      
      directions_use <- if(deviation == 0) "Null" else directions
      
      map_dfr(directions_use, function(direction){
        
        map_dfr(seq_len(n_sim), function(i){
          
          simulated <- simulate_deviation_dataset(
            model = primary_model,
            study_species = field_species,
            n_species = n_species,
            n_ind = n_ind,
            deviation = deviation,
            direction = direction
          )
          
          # Test whether species-level departures from the Astley prediction
          # have a mean significantly different from zero.
          departure_test <- t.test(
            simulated$departure,
            mu = 0
          )
          
          mean_departure <- mean(simulated$departure)
          sd_departure <- sd(simulated$departure)
          p_value <- departure_test$p.value
          detectable <- p_value < alpha
          
          # Check whether the detected departure is in the simulated direction.
          correct_direction <- case_when(
            deviation == 0 ~ NA,
            direction == "Positive" & detectable & mean_departure > 0 ~ TRUE,
            direction == "Negative" & detectable & mean_departure < 0 ~ TRUE,
            TRUE ~ FALSE
          )
          
          tibble(
            Species = n_species,
            Individuals = n_ind,
            Deviation = deviation * 100,
            Direction = direction,
            Mean_departure = mean_departure,
            SD_departure = sd_departure,
            P_value = p_value,
            Detectable = detectable,
            Correct_direction = correct_direction
          )
        })
      })
    })
    
  })
})

# Summarise power
power_summary <- power_results %>%
  group_by(Species, Individuals, Deviation, Direction) %>%
  summarise(
    Power = mean(Detectable, na.rm = TRUE),
    Directional_power = mean(Correct_direction, na.rm = TRUE),
    Mean_departure = mean(Mean_departure, na.rm = TRUE),
    SD_departure = sd(Mean_departure, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Power = round(Power, 3),
    Directional_power = round(Directional_power, 3),
    Mean_departure = round(Mean_departure, 4),
    SD_departure = round(SD_departure, 4)
  )

print(power_summary)


# Under the null, the true mean departure is zero. The proportion of simulations
# declared significant should therefore be approximately alpha = 0.05.
null_calibration <- power_summary %>%
  filter(Deviation == 0) %>%
  select(Species, Individuals, Direction, Power) %>%
  rename(Empirical_Type_I_Error = Power)

print(null_calibration)

# Identify scenarios achieving >=80% directional power
power_80 <- power_summary %>%
  filter(Deviation > 0) %>%
  mutate(Adequate_power = Directional_power >= 0.80)

print(power_80)

# Minimum design achieving >=80% power
minimum_power_design <- power_80 %>%
  filter(Adequate_power) %>%
  arrange(Species, Individuals)

print(minimum_power_design)

# Minimum design separately for each sensitivity scenario
minimum_by_effect <- power_80 %>%
  filter(Adequate_power) %>%
  group_by(Deviation, Direction) %>%
  arrange(Species, Individuals) %>%
  slice(1) %>%
  ungroup()

print(minimum_by_effect)

# Compact power table
power_table <- power_summary %>%
  select(Species, Individuals, Deviation, Direction, Power, Directional_power,
         Mean_departure, SD_departure) %>%
  arrange(Deviation, Direction, Species, Individuals)

print(power_table)

# Save primary-trait results
write.csv(power_summary, "Deviation_Power_Sensitivity_Analysis.csv", row.names = FALSE)
write.csv(minimum_by_effect, "Minimum_Design_By_Effect_Size.csv", row.names = FALSE)
write.csv(null_calibration, "Null_Type_I_Error_Calibration.csv", row.names = FALSE)

# Decision table
decision_table <- minimum_by_effect %>%
  select(Deviation, Direction, Species, Individuals, Directional_power,
         Mean_departure, SD_departure) %>%
  arrange(Deviation, Direction)

print(decision_table)

# Repeat the sensitivity analysis for all locomotor traits
all_trait_power <- map_dfr(names(trait_models), function(trait_name){
  
  current_model <- trait_models[[trait_name]]
  
  map_dfr(species_numbers, function(n_species){
    map_dfr(individual_numbers, function(n_ind){
      map_dfr(deviation_levels, function(deviation){
        directions_use <- if(deviation == 0) "Null" else directions
        
        map_dfr(directions_use, function(direction){
          
          map_dfr(seq_len(n_sim), function(i){
            
            simulated <- simulate_deviation_dataset(
              model = current_model,
              study_species = field_species,
              n_species = n_species,
              n_ind = n_ind,
              deviation = deviation,
              direction = direction
            )
            
            departure_test <- t.test(simulated$departure, mu = 0)
            
            mean_departure <- mean(simulated$departure)
            detectable <- departure_test$p.value < alpha
            
            correct_direction <- case_when(
              deviation == 0 ~ NA,
              direction == "Positive" & detectable & mean_departure > 0 ~ TRUE,
              direction == "Negative" & detectable & mean_departure < 0 ~ TRUE,
              TRUE ~ FALSE
            )
            
            tibble(
              Trait = trait_name,
              Species = n_species,
              Individuals = n_ind,
              Deviation = deviation * 100,
              Direction = direction,
              Mean_departure = mean_departure,
              Detectable = detectable,
              Correct_direction = correct_direction
            )
          })
        })
      })
    })
  })
})

# Summarise power across locomotor traits
all_trait_power_summary <- all_trait_power %>%
  group_by(Trait, Species, Individuals, Deviation, Direction) %>%
  summarise(
    Power = mean(Detectable, na.rm = TRUE),
    Directional_power = mean(Correct_direction, na.rm = TRUE),
    Mean_departure = mean(Mean_departure, na.rm = TRUE),
    SD_departure = sd(Mean_departure, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Power = round(Power, 3),
    Directional_power = round(Directional_power, 3),
    Mean_departure = round(Mean_departure, 4),
    SD_departure = round(SD_departure, 4)
  )

print(all_trait_power_summary)

# Save all-trait results
write.csv(all_trait_power_summary, "All_Traits_Deviation_Power_Analysis.csv", row.names = FALSE)





#-------------------------------------------------------------------------------
# 3.6 Integrated sampling design sensitivity (site + phylogeny)
#-------------------------------------------------------------------------------

# Evaluates realistic field sampling scenarios after accounting for the MVD.
#
# Only designs meeting the minimum viable dataset from Section 3.5 are retained.
# Candidate designs are ranked by:
# 1) Faith's phylogenetic diversity
# 2) number of species
# 3) representation of miniaturised and non-miniaturised taxa
# 4) fewer field sites required
#
# This identifies the strongest achievable sampling design for testing:
# H2 scaling vs deviation from scaling
# H3 convergence vs unique evolutionary solutions

# Minimum viable dataset established from Section 3.5
minimum_species <- 5

# Expand taxa occurring at multiple sites
# Example: Plethodontohyla_notosticta contributes to both Yellow and Red
field_species_expanded <- field_species %>%
  separate_rows(Site,sep=";")


# Available field sites
available_sites <- unique(
  field_species_expanded$Site
)


# Generate all possible site combinations
site_combinations <- map_dfr(
  2:length(available_sites),
  function(n_sites){
    
    map_dfr(
      combn(
        available_sites,
        n_sites,
        simplify=FALSE
      ),
      function(site_set){
        
        # Species available from chosen sites
        
        selected_species <- field_species_expanded %>%
          filter(Site %in% site_set) %>%
          distinct(Species,.keep_all=TRUE)
        
        
        # Prune phylogeny to sampled taxa
        
        pruned_tree <- drop.tip(
          study_tree,
          setdiff(
            study_tree$tip.label,
            selected_species$Species
          )
        )
        
        
        # Calculate Faith's phylogenetic diversity
        
        faith_pd <- picante::pd(
          matrix(
            1,
            nrow=1,
            ncol=nrow(selected_species),
            dimnames=list(
              "community",
              selected_species$Species
            )
          ),
          pruned_tree,
          include.root=TRUE
        )$PD
        
        tibble(
          Sites=paste(site_set,collapse="+"),
          Number_sites=n_sites,
          Species=nrow(selected_species),
          Mini_species=sum(selected_species$SVL<=20),
          Nonmini_species=sum(selected_species$SVL>20),
          Faith_PD=faith_pd,
          
          Mean_Difficulty=mean(
            selected_species$Difficulty,
            na.rm=TRUE
          ),
          
          Total_Difficulty=sum(
            selected_species$Difficulty,
            na.rm=TRUE
          )
        )
        
      }
      
    )
    
  }
  
)


# Remove designs below the minimum viable dataset
sampling_designs <- site_combinations %>%
  filter(
    Species>=minimum_species
  ) %>%
  arrange(
    desc(Faith_PD),
    desc(Species),
    desc(Mini_species),
    desc(Nonmini_species),
    Mean_Difficulty,
    Total_Difficulty,
    Number_sites
  ) %>%
  mutate(
    Rank=row_number()
  ) %>%
  select(
    Rank,
    everything()
  )


# View ranked designs
print(sampling_designs,n=Inf)


# Export final design ranking
write.csv(
  sampling_designs,
  "Integrated_sampling_designs_ranked.csv",
  row.names=FALSE
)




#-------------------------------------------------------------------------------
# 3.7 Final integrated sampling design optimisation
#-------------------------------------------------------------------------------
# Evaluate realistic field designs:
#
# Mandatory:
# - Seychelles included
#
# Compare:
# - Seychelles + 1 additional site
# - Seychelles + 2 additional sites
# - Seychelles + 3 additional sites (realistic maximum)
# - Seychelles + 4 additional sites (diminishing return comparison)
#
# Designs are evaluated by:
# - Faith's PD
# - MPD
# - miniaturised species
# - miniaturised families
# - miniaturised genera
#
# Species number is retained as descriptive information only.

# Add genus information
field_species_expanded <- field_species_expanded %>%
  mutate(
    Genus=sub("_.*","",Species)
  )


# Select candidate site combinations
candidate_designs <- site_combinations %>%
  filter(
    grepl("Seychelles",Sites),
    Number_sites <= 5, Species >= minimum_species
  )

# Calculate evolutionary returns
design_summary <- map_dfr(
  1:nrow(candidate_designs),
  function(i){
    
    sites <- strsplit(
      candidate_designs$Sites[i],
      "\\+"
    )[[1]]
    
    
    selected_species <- field_species_expanded %>%
      filter(Site %in% sites) %>%
      distinct(Species,.keep_all=TRUE)
    
    
    # Require minimum viable dataset
    
    if(nrow(selected_species) < minimum_species){
      
      return(NULL)
      
    }
    
    
    # prune full phylogeny
    
    selected_tree <- drop.tip(
      tree,
      setdiff(
        tree$tip.label,
        selected_species$Species
      )
    )
    
    
    if(length(selected_tree$tip.label)<2){
      
      return(NULL)
      
    }
    
    
    # Faith PD
    
    species_matrix <- matrix(
      1,
      nrow=1,
      ncol=length(selected_tree$tip.label)
    )
    
    colnames(species_matrix) <- selected_tree$tip.label
    
    
    faith_pd <- picante::pd(
      species_matrix,
      selected_tree,
      include.root=TRUE
    )$PD
    
    
    # Mean pairwise distance
    
    distance_matrix <- cophenetic(selected_tree)
    
    mpd_value <- mean(
      distance_matrix[
        upper.tri(distance_matrix)
      ]
    )
    
    
    tibble(
      
      Sites=candidate_designs$Sites[i],
      
      Number_sites=length(sites),
      
      Additional_sites=length(sites)-1,
      
      Species_total=nrow(selected_species),
      
      Mini_species=sum(
        selected_species$SVL <=20
      ),
      
      Mini_families=n_distinct(
        selected_species$Family[
          selected_species$SVL <=20
        ]
      ),
      
      Mini_genera=n_distinct(
        selected_species$Genus[
          selected_species$SVL <=20
        ]
      ),
      
      Faith_PD=faith_pd,
      
      MPD=mpd_value,
      
      Mean_Difficulty=mean(
        selected_species$Difficulty,
        na.rm=TRUE
      ),
      
      Total_Difficulty=sum(
        selected_species$Difficulty,
        na.rm=TRUE
      )
      
    )
    
  }
)


# Summarise returns by sampling effort
design_summary <- design_summary %>%
  group_by(Additional_sites) %>%
  mutate(
    
    PD_gain=Faith_PD-min(Faith_PD),
    
    MPD_gain=MPD-min(MPD)
    
  ) %>%
  ungroup() %>%
  arrange(
    Additional_sites,
    desc(Faith_PD)
  )

print(design_summary, n=Inf)

write.csv(
  design_summary,
  "Final_sampling_design_optimisation.csv",
  row.names=FALSE
)


# Integrated field design contingency analysis
#
# Purpose: Evaluate robustness of the proposed field sampling design.
#
# Scenarios:
#
# 1. Preferred design:
#    Seychelles + Blue + Red + Purple
#
# 2. Replacement:
#    One site unavailable, replaced by the best alternative site combination.
#
# 3. Loss without replacement:
#    One site unavailable and no additional site added.
#
# Designs are evaluated by:
#
# - miniaturised species replication
# - miniaturised genera
# - miniaturised families
# - Faith's phylogenetic diversity
# - mean pairwise distance
#
# Species number is retained as descriptive information only.


library(dplyr)
library(purrr)

# Define preferred design
preferred_sites <- c(
  "Seychelles",
  "Andasibe_Blue",
  "Marojejy_White",
  "MontagneDAmbre_Purple"
)

# Function to evaluate a design
evaluate_design <- function(site_set){
  
  
  selected_species <- field_species_expanded %>%
    filter(Site %in% site_set) %>%
    distinct(Species,.keep_all=TRUE) %>%
    filter(Species %in% tree$tip.label)
  
  if(nrow(selected_species)<2){
    return(NULL)
  }
  
  pruned_tree <- drop.tip(
    tree,
    setdiff(
      tree$tip.label,
      selected_species$Species
    )
  )
  
  community <- matrix(
    1,
    nrow=1,
    ncol=nrow(selected_species)
  )
  
  colnames(community) <- selected_species$Species
  
  
  faith_pd <- picante::pd(
    community,
    pruned_tree,
    include.root=TRUE
  )$PD
  
  
  mpd_value <- picante::mpd(
    community,
    cophenetic(pruned_tree)
  )
  
  tibble(
    
    Sites=paste(site_set,collapse="+"),
    
    Number_sites=length(site_set),
    
    Species_total=nrow(selected_species),
    
    Mini_species=sum(selected_species$SVL<=20),
    
    Mini_families=n_distinct(
      selected_species$Family[selected_species$SVL<=20]
    ),
    
    Mini_genera=n_distinct(
      selected_species$Genus[selected_species$SVL<=20]
    ),
    
    Faith_PD=faith_pd,
    
    MPD=mpd_value,
    
    Mean_Difficulty=mean(
      selected_species$Difficulty,
      na.rm=TRUE
    ),
    
    Total_Difficulty=sum(
      selected_species$Difficulty,
      na.rm=TRUE
    )
  )
  
}

# Scenario 1: preferred design
preferred_result <- evaluate_design(
  preferred_sites
) %>%
  mutate(
    Scenario="Preferred design",
    Lost_site="None"
  )

# Scenario 2: replacement if one site unavailable
all_sites <- unique(
  field_species_expanded$Site
)

replacement_results <- map_dfr(
  
  preferred_sites,
  
  function(site_lost){
    
    
    available_sites <- setdiff(
      all_sites,
      site_lost
    )
    
    
    candidate_designs <- combn(
      available_sites,
      4,
      simplify=FALSE
    )
    
    
    results <- map_dfr(
      candidate_designs,
      evaluate_design
    )
    
    
    results %>%
      arrange(
        desc(Mini_species),
        desc(Mini_genera),
        desc(Mini_families),
        desc(Faith_PD),
        desc(MPD)
      ) %>%
      slice(1) %>%
      mutate(
        Scenario="Replacement available",
        Lost_site=site_lost
      )
    
  }
)


# Scenario 3: loss without replacement
loss_results <- map_dfr(
  
  preferred_sites,
  
  function(site_lost){
    
    
    remaining_sites <- setdiff(
      preferred_sites,
      site_lost
    )
    
    
    evaluate_design(
      remaining_sites
    ) %>%
      mutate(
        Scenario="No replacement",
        Lost_site=site_lost
      )
    
  }
)


# Combine results
integrated_design_summary <- bind_rows(
  
  preferred_result,
  
  replacement_results,
  
  loss_results
  
) %>%
  select(
    Scenario,
    Lost_site,
    Sites,
    Number_sites,
    Species_total,
    Mini_species,
    Mini_families,
    Mini_genera,
    Faith_PD,
    MPD,
    Mean_Difficulty,
    Total_Difficulty
  )

print(integrated_design_summary)

write.csv(
  integrated_design_summary,
  "Integrated_field_design_contingency_analysis.csv",
  row.names=FALSE
)


