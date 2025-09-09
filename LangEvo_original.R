# Packages
library(igraph)

# === PARAMETERS ===
N <- 30  # number of agents
t_max <- 200  # number of timesteps
network_type <- "small-world"  # options: "random", "small-world", "regular", "fully-connected", "data-based"
p_link <- 0.04  # probability of connection (for random/small-world)
include_outgroup <- TRUE  # include outgroup connections
agent_turnover <- FALSE  # implement agent turnover over time
base_interaction_rate <- 5  # average interactions per timestep when weight = 1.0

# === INITIALIZATION ===
set.seed(65328)  # for reproducibility

## Initialize agents with linguistic variables
agent <- data.frame(
  id = 1:N,
  generation = sample(1:3, N, replace = TRUE, prob = c(0.8, 0.15, 0.05)), # 1=1st gen, 2=2nd gen, 3=3rd gen
  init_exp_greek = pmax(0, pmin(1, rnorm(N, mean = 0.4, sd = 0.2))), # exposure to Greek, Gaussian
  init_exp_french = pmax(0, pmin(1, rnorm(N, mean = 0.6, sd = 0.2))), # exposure to French, slightly higher mean of Gaussian
  theta_0 = rep(NA, N),  # initial prior (subject dropping probability)
  theta = rep(NA, N),    # current theta
  n = rep(0, N),         # total interactions
  k = rep(0, N),         # subject drops observed
  v = rep(0, N),          # velocity of change (Δθ/Δt)
  # Dynamic effect modifiers
  generation_learning_modifier = rep(NA, N),    # how generation affects learning
  greek_learning_modifier = rep(NA, N),         # how Greek exposure affects learning
  french_learning_modifier = rep(NA, N),        # how French exposure affects learning
  total_learning_modifier = rep(NA, N)          # combined learning rate modifier
)

## Set initial theta based on generation and Greek exposure
# 1st gen: higher theta (more dropping), 2nd gen: lower, 3rd gen: almost none
gen_effect <- c(0.7, 0.3, 0.1)  # base theta by generation
#Set prior distribution cases
p <- seq(0, 1, by = 0.001)
alpha1 <- 100
beta1 <- 1
prior1 <- dbeta(p, alpha1, beta1)
# Set each agent's theta
for(i in 1:N) {
  # Everyone starts with the same prior distribution
  agent$theta_0[i] <- rbeta(1, alpha1, beta1)
}

## Apply generation effects as post-initialization adjustments
generation_effects <- function(agent_df) {
  # Define how each generation differs in their linguistic behavior
  gen_theta_adjustment <- c(0.0, -0.3, -0.4)    # 1st gen: baseline, newer generations are more affected by French as the language of their childhood environment
  gen_learning_modifier <- c(0.9, 1.1, 1.2)    # 1st gen: slower learning, 3rd gen: faster learning
  
  for(i in 1:nrow(agent_df)) {
    gen <- agent_df$generation[i]
    
    # Adjust initial theta based on generation
    agent_df$theta_0[i] <- pmax(0, pmin(1, agent_df$theta_0[i] + gen_theta_adjustment[gen]))
    
    # Set generation-based learning modifier
    agent_df$generation_learning_modifier[i] <- gen_learning_modifier[gen]
  }
  
  return(agent_df)
}

## Apply Greek exposure effects
greek_exposure_effects <- function(agent_df) {
  for(i in 1:nrow(agent_df)) {
    greek_exp <- agent_df$init_exp_greek[i]
    
    # Greek exposure increases initial theta (more subject dropping)
    greek_theta_effect <- (greek_exp - 0.5) * 0.3  # centered at 0.5, scaled by 0.3
    agent_df$theta_0[i] <- pmax(0, pmin(1, agent_df$theta_0[i] + greek_theta_effect))
    
    # Greek exposure affects learning rate (higher exposure = more sensitive to Greek-like patterns)
    # This could make agents more likely to adopt dropping behavior they observe
    agent_df$greek_learning_modifier[i] <- 1 + (greek_exp * 0.5)  # 1.0 to 1.5 range
  }
  
  return(agent_df)
}

## Apply French exposure effects
french_exposure_effects <- function(agent_df) {
  for(i in 1:nrow(agent_df)) {
    french_exp <- agent_df$init_exp_french[i]
    
    # French exposure decreases initial theta (less subject dropping)
    french_theta_effect <- (french_exp - 0.5) * -0.2  # negative effect, scaled by 0.2
    agent_df$theta_0[i] <- pmax(0, pmin(1, agent_df$theta_0[i] + french_theta_effect))
    
    # French exposure makes agents more resistant to adopting dropping behavior
    # (lower learning rate for pro-dropping evidence)
    agent_df$french_learning_modifier[i] <- 1 - (french_exp * 0.3)  # 0.7 to 1.0 range
  }
  
  return(agent_df)
}

## Apply all effects in sequence
agent <- generation_effects(agent)
agent <- greek_exposure_effects(agent)
agent <- french_exposure_effects(agent)

## Calculate total learning modifier (combination of all effects)
for(i in 1:N) {
  agent$total_learning_modifier[i] <- agent$generation_learning_modifier[i] * 
    agent$greek_learning_modifier[i] * 
    agent$french_learning_modifier[i]
}

## Set current theta to initial theta
agent$theta <- agent$theta_0

## Initialize n(0) and k(0) based on theta(0)
initial_interactions <- 10
for(i in 1:N) {
  agent$n[i] <- initial_interactions
  agent$k[i] <- rbinom(1, initial_interactions, agent$theta[i])
}

# === NETWORK GENERATION ===
generate_network <- function(N, type, p_link = 0.04, include_outgroup = TRUE) {
  if(type == "random") {
    g <- sample_gnp(N, p_link, directed = FALSE)
  } else if(type == "small-world") {
    nei <- max(2, round(N * p_link / 2))  # ensure at least 2 neighbors
    g <- sample_smallworld(1, N, nei, p_link)
  } else if(type == "regular") {
    nei <- max(2, round(N * p_link))
    if(nei %% 2 == 1) nei <- nei - 1  # make even for regular graph
    g <- make_lattice(dimvector = N, periodic = TRUE, directed = FALSE)
  } else if(type == "fully-connected") {
    g <- make_full_graph(N, directed = FALSE)
  } else {
    # Default to small-world
    nei <- max(2, round(N * p_link / 2))
    g <- sample_smallworld(1, N, nei, p_link)
  }
  
  # Add outgroup node if requested
  if(include_outgroup) {
    g <- add_vertices(g, 1)  # add outgroup node
    N_total <- N + 1
    # Connect agents to outgroup
    for(i in 1:N) {
      g <- add_edges(g, c(i, N_total))
      }
    }
  
  return(g)
}

# Generate network
g <- generate_network(N, network_type, p_link, include_outgroup)

# Add edge weights (strength of linguistic influence)
grmin <- 0.1
grmax <- 1
E(g)$weight <- runif(length(E(g)), grmin, grmax)

# If outgroup exists, give it lower weights
if(include_outgroup) {
  frmin <- 0.05
  frmax <- 0.3 #
  outgroup_edges <- incident(g, vcount(g), mode = "all")
  E(g)$weight[outgroup_edges] <- runif(length(outgroup_edges), frmin, frmax)
}

# Visualize initial network without outgroup node
igraph.options(print.id = F)
outgroup_index <- vcount(g)
incident_edges <- incident(g, outgroup_index, mode = "all")
n_colors <- 100  # Number of color gradations
color_palette <- heat.colors(n_colors)

# Map theta values to color indices (clamped between 0 and 1)
theta_clamped <- pmax(0, pmin(1, agent$theta))  # Ensure values are between 0 and 1
color_indices <- round(theta_clamped * (n_colors - 1)) + 1  # Map to color palette indices

V(g)$color <- ifelse(1:vcount(g) <= N, color_palette[color_indices], "black")
V(g)$size <- 3
V(g)[outgroup_index]$color <- "transparent"  # make outgroup node invisible
V(g)$size <- 10  # default size
V(g)[outgroup_index]$size <- 0
incident_edges <- incident(g, outgroup_index)
E(g)$color <- "black"  # default edge color
E(g)[incident_edges]$color <- "transparent"  # make outgroup links invisible

plot(g, 
     vertex.label = ifelse(1:vcount(g) <= N, 1:N, "OUT"),
     vertex.label = NA,
     edge.arrow.size = 0,
     main = paste("Initial Network:", network_type, "with outgroup (links not depicted here)"))

legend("topright", 
       legend = c("0", "0.5", "1"), 
       col = "black",  # border color
       pt.bg = color_palette[c(1, n_colors/2, n_colors)],  # fill color
       pch = 21,  # filled circle with border
       title = "θ values",
       cex = 0.8)

# === SIMULATION FUNCTIONS ===

# Function to perform Bayesian update
bayesian_update <- function(current_theta, k_obs, n_obs, learning_modifier = 1.0, base_learning_rate = 0.1) {
  if(n_obs == 0) return(current_theta)
  
  observation <- k_obs / n_obs
  # Apply the learning modifier to the base learning rate
  adjusted_learning_rate <- base_learning_rate * learning_modifier
  
  # Bayesian update with modified learning rate
  new_theta <- current_theta * (1 - adjusted_learning_rate) + observation * adjusted_learning_rate
  
  # Keep theta in bounds [0, 1]
  return(pmax(0, pmin(1, new_theta)))
}

# Function to simulate one timestep
simulate_timestep <- function(agent_df, graph, h) {
  N_agents <- nrow(agent_df)
  new_agent_df <- agent_df
  
  # Store previous theta for velocity calculation
  prev_theta <- agent_df$theta
  
  # For each agent
  for(i in 1:N_agents) {
    # Reset interaction counters for this timestep
    n_hi <- 0 # try also 0.2 * previous n
    k_hi <- 0
    
    # Get neighbors (excluding outgroup for now)
    neighbors <- neighbors(graph, i)
    neighbors <- neighbors[neighbors <= N_agents]  # exclude outgroup
    
    if(length(neighbors) > 0) {
      # Circle through neighbors and collect their data
      for(j in neighbors) {
        # Get edge weight
        edge_id <- get.edge.ids(graph, c(i, j))
        weight <- E(graph)$weight[edge_id]
        
        # Read neighbor's n and k (weighted by connection strength)
        n_hij <- agent_df$n[j] * weight
        k_hij <- agent_df$k[j] * weight
        
        # Sum into agent i's observations
        n_hi <- n_hi + n_hij
        k_hi <- k_hi + k_hij
      }
      
      # Update theta based on neighbor observations
      if(n_hi > 0) {
        new_agent_df$theta[i] <- bayesian_update(agent_df$theta[i],
            k_hi, n_hi, learning_modifier = agent_df$total_learning_modifier[i])
      }
    }
    
    # Handle outgroup influence if connected
    if(include_outgroup && vcount(graph) > N_agents) {
      outgroup_id <- vcount(graph)
      if(are.connected(graph, i, outgroup_id)) {
        edge_id <- get.edge.ids(graph, c(i, outgroup_id))
        outgroup_weight <- E(graph)$weight[edge_id]
        
        # Outgroup has fixed theta (e.g., diaspora = low dropping)
        outgroup_theta <- 0
        outgroup_n <- 5
        outgroup_k <- outgroup_n * outgroup_theta
        
        # Apply outgroup influence (weighted)
        new_agent_df$theta[i] <- bayesian_update(new_agent_df$theta[i], 
                                                 outgroup_k * outgroup_weight, 
                                                 outgroup_n * outgroup_weight,
                                                 learning_modifier = agent_df$total_learning_modifier[i])
      }
    }
  }
  
  # Generate new interactions for next timestep based on updated theta
  for(i in 1:N_agents) {
    neighbors <- neighbors(graph, i)
    neighbors <- neighbors[neighbors <= N_agents]
    
    if(length(neighbors) > 0) {
      for(j in neighbors) {
        edge_id <- get.edge.ids(graph, c(i, j))
        weight <- E(graph)$weight[edge_id]
        
        # Generate interactions based on current theta and connection weight
        n_interactions <- rpois(1, 5 * weight)  # Poisson-distributed interactions, vary 5
        k_drops <- rbinom(1, n_interactions, new_agent_df$theta[i])
        
        # Update agent j's observations from agent i
        new_agent_df$n[j] <- new_agent_df$n[j] + n_interactions
        new_agent_df$k[j] <- new_agent_df$k[j] + k_drops
      }
    }
  }
  
  # Calculate velocity (change in theta)
  new_agent_df$v <- new_agent_df$theta - prev_theta ### try moving avg?
  
  return(new_agent_df)
}

# === TRACKING VARIABLES ===
totals <- data.frame(
  timestep = 1:t_max,
  k_total = rep(NA, t_max),
  theta_avg = rep(NA, t_max),
  theta_sd = rep(NA, t_max),
  v_avg = rep(NA, t_max),
  clustering = rep(NA, t_max)
)

# Store agent history for analysis
agent_history <- array(NA, dim = c(N, t_max, 4))  # agents x time x variables (theta, k, n, v)

# === MAIN SIMULATION LOOP ===
cat("Running simulation for", t_max, "timesteps...\n")

for(h in 1:t_max) {
  # Simulate one timestep
  agent <- simulate_timestep(agent, g, h)
  
  # Calculate totals
  totals$k_total[h] <- sum(agent$k)
  totals$theta_avg[h] <- mean(agent$theta)
  totals$theta_sd[h] <- sd(agent$theta)
  totals$v_avg[h] <- mean(abs(agent$v))
  totals$clustering[h] <- transitivity(g, type = "global")
  
  # Store agent states
  agent_history[, h, 1] <- agent$theta
  agent_history[, h, 2] <- agent$k
  agent_history[, h, 3] <- agent$n
  agent_history[, h, 4] <- agent$v
  
  # Optional: Print progress
  if(h %% 20 == 0) cat("Completed timestep", h, "\n")
}

# === ANALYSIS ===
cat("\n=== SIMULATION RESULTS ===\n")

# Final vs initial theta
cat("Correlation θ(final) ~ θ(0):", cor(agent$theta, agent$theta_0), "\n")

# Generation effects
for(gen in 1:3) {
  gen_agents <- agent$generation == gen
  cat("Generation", gen, "- Final θ mean:", mean(agent$theta[gen_agents]), "\n")
}

# Function to analyze generation effects
analyze_generation_effects <- function(gen_probs, scenario_name) {
  test_generation <- sample(1:3, 1000, replace = TRUE, prob = gen_probs) # Simulate with given probabilities
  # Calculate expected initial theta
  gen_effects <- c(0.0, -0.3, -0.4)  # as defined in your model
  expected_theta_adjustment <- sum(gen_probs * gen_effects)
  # Calculate expected learning rate
  gen_learning <- c(0.9, 1.1, 1.2)
  expected_learning_rate <- sum(gen_probs * gen_learning)
  
  cat("=== SCENARIO:", scenario_name, "===\n")
  cat("Generation distribution:", paste0(gen_probs * 100, "%"), "\n")
  cat("Expected theta adjustment:", round(expected_theta_adjustment, 3), "\n")
  cat("Expected learning rate:", round(expected_learning_rate, 3), "\n")
  cat("Interpretation: ")
  
  if(expected_theta_adjustment > 0.05) {
    cat("Community will start with higher dropping rates (heritage language maintenance)\n")
  } else if(expected_theta_adjustment < -0.05) {
    cat("Community will start with lower dropping rates (language shift toward host language)\n")
  } else {
    cat("Community will start with moderate dropping rates (balanced)\n")
  }
  
  if(expected_learning_rate > 1.05) {
    cat("Community will change quickly (young, flexible speakers)\n")
  } else if(expected_learning_rate < 0.95) {
    cat("Community will change slowly (older, more fixed speakers)\n")
  } else {
    cat("Community will change at moderate pace\n")
  }
  cat("\n")
}

analyze_generation_effects(c(0.8, 0.15, 0.05), "Greek communities in France")

# Greek exposure effects
cat("Correlation θ(final) ~ initial Greek exposure:", cor(agent$theta, agent$init_exp_greek), "\n")

analyze_greek_exposure_effects <- function(greek_mean = 0.4, greek_sd = 0.2, scenario_name = "") {
  # Simulate Greek exposure distribution
  n_sim <- 10000
  greek_exp_sim <- pmax(0, pmin(1, rnorm(n_sim, mean = greek_mean, sd = greek_sd)))
  
  # Calculate effects based on your model parameters
  # Greek theta effect: (greek_exp - 0.5) * 0.3
  greek_theta_effects <- (greek_exp_sim - 0.5) * 0.3
  
  # Greek learning modifier: 1 + (greek_exp * 0.5) [range 1.0 to 1.5]
  greek_learning_modifiers <- 1 + (greek_exp_sim * 0.5)
  
  cat("=== GREEK EXPOSURE ANALYSIS", ifelse(scenario_name != "", paste(":", scenario_name), ""), "===\n")
  cat("Distribution parameters: mean =", greek_mean, ", sd =", greek_sd, "\n")
  cat("Actual exposure range:", round(min(greek_exp_sim), 3), "to", round(max(greek_exp_sim), 3), "\n")
  cat("Mean exposure:", round(mean(greek_exp_sim), 3), "\n\n")
  
  cat("THETA EFFECTS (initial subject dropping probability):\n")
  cat("  Theta adjustment range:", round(min(greek_theta_effects), 3), "to", round(max(greek_theta_effects), 3), "\n")
  cat("  Mean theta adjustment:", round(mean(greek_theta_effects), 3), "\n")
  cat("  Agents with positive theta effect (> baseline):", round(sum(greek_theta_effects > 0) / length(greek_theta_effects) * 100, 1), "%\n")
  
  cat("\nLEARNING RATE EFFECTS:\n")
  cat("  Learning modifier range:", round(min(greek_learning_modifiers), 3), "to", round(max(greek_learning_modifiers), 3), "\n")
  cat("  Mean learning modifier:", round(mean(greek_learning_modifiers), 3), "\n")
  cat("  Agents with enhanced learning (> 1.0):", round(sum(greek_learning_modifiers > 1.0) / length(greek_learning_modifiers) * 100, 1), "%\n")
  
  cat("\nINTERPRETATION:\n")
  if(mean(greek_theta_effects) > 0.05) {
    cat("  - Community will start with HIGHER subject dropping (Greek heritage maintenance)\n")
  } else if(mean(greek_theta_effects) < -0.05) {
    cat("  - Community will start with LOWER subject dropping (Greek attrition)\n")
  } else {
    cat("  - Community will start with MODERATE subject dropping (balanced Greek influence)\n")
  }
  
  if(mean(greek_learning_modifiers) > 1.1) {
    cat("  - Community will be HIGHLY responsive to pro-dropping evidence\n")
  } else if(mean(greek_learning_modifiers) > 1.05) {
    cat("  - Community will be MODERATELY responsive to pro-dropping evidence\n")
  } else {
    cat("  - Community will have STANDARD responsiveness to pro-dropping evidence\n")
  }
  cat("\n")
}

# Function to analyze French exposure effects
cat("Correlation θ(final) ~ initial French exposure:", cor(agent$theta, agent$init_exp_french), "\n")

analyze_french_exposure_effects <- function(french_mean = 0.6, french_sd = 0.2, scenario_name = "") {
  # Simulate French exposure distribution
  n_sim <- 10000
  french_exp_sim <- pmax(0, pmin(1, rnorm(n_sim, mean = french_mean, sd = french_sd)))
  
  # Calculate effects based on your model parameters
  # French theta effect: (french_exp - 0.5) * -0.2
  french_theta_effects <- (french_exp_sim - 0.5) * -0.2
  
  # French learning modifier: 1 - (french_exp * 0.3) [range 0.7 to 1.0]
  french_learning_modifiers <- 1 - (french_exp_sim * 0.3)
  
  cat("=== FRENCH EXPOSURE ANALYSIS", ifelse(scenario_name != "", paste(":", scenario_name), ""), "===\n")
  cat("Distribution parameters: mean =", french_mean, ", sd =", french_sd, "\n")
  cat("Actual exposure range:", round(min(french_exp_sim), 3), "to", round(max(french_exp_sim), 3), "\n")
  cat("Mean exposure:", round(mean(french_exp_sim), 3), "\n\n")
  
  cat("THETA EFFECTS (initial subject dropping probability):\n")
  cat("  Theta adjustment range:", round(min(french_theta_effects), 3), "to", round(max(french_theta_effects), 3), "\n")
  cat("  Mean theta adjustment:", round(mean(french_theta_effects), 3), "\n")
  cat("  Agents with negative theta effect (< baseline):", round(sum(french_theta_effects < 0) / length(french_theta_effects) * 100, 1), "%\n")
  
  cat("\nLEARNING RATE EFFECTS:\n")
  cat("  Learning modifier range:", round(min(french_learning_modifiers), 3), "to", round(max(french_learning_modifiers), 3), "\n")
  cat("  Mean learning modifier:", round(mean(french_learning_modifiers), 3), "\n")
  cat("  Agents with reduced learning (< 1.0):", round(sum(french_learning_modifiers < 1.0) / length(french_learning_modifiers) * 100, 1), "%\n")
  
  cat("\nINTERPRETATION:\n")
  if(mean(french_theta_effects) < -0.05) {
    cat("  - Community will start with LOWER subject dropping (French influence dominant)\n")
  } else if(mean(french_theta_effects) > -0.02) {
    cat("  - Community will start with MODERATE subject dropping (weak French influence)\n")
  } else {
    cat("  - Community will start with MODERATE-LOW subject dropping (moderate French influence)\n")
  }
  
  if(mean(french_learning_modifiers) < 0.85) {
    cat("  - Community will be HIGHLY resistant to pro-dropping evidence\n")
  } else if(mean(french_learning_modifiers) < 0.95) {
    cat("  - Community will be MODERATELY resistant to pro-dropping evidence\n")
  } else {
    cat("  - Community will have STANDARD resistance to pro-dropping evidence\n")
  }
  cat("\n")
}

# Function to analyze combined Greek and French effects
analyze_combined_exposure_effects <- function(greek_mean = 0.4, greek_sd = 0.2, 
                                              french_mean = 0.6, french_sd = 0.2,
                                              correlation = 0, scenario_name = "") {
  n_sim <- 10000
  
  # Generate correlated exposures if specified
  if(correlation != 0) {
    library(MASS)
    sigma <- matrix(c(greek_sd^2, correlation * greek_sd * french_sd,
                      correlation * greek_sd * french_sd, french_sd^2), 2, 2)
    exposures <- mvrnorm(n_sim, mu = c(greek_mean, french_mean), Sigma = sigma)
    greek_exp_sim <- pmax(0, pmin(1, exposures[,1]))
    french_exp_sim <- pmax(0, pmin(1, exposures[,2]))
  } else {
    greek_exp_sim <- pmax(0, pmin(1, rnorm(n_sim, mean = greek_mean, sd = greek_sd)))
    french_exp_sim <- pmax(0, pmin(1, rnorm(n_sim, mean = french_mean, sd = french_sd)))
  }
  
  # Calculate combined effects
  greek_theta_effects <- (greek_exp_sim - 0.5) * 0.3
  french_theta_effects <- (french_exp_sim - 0.5) * -0.2
  combined_theta_effects <- greek_theta_effects + french_theta_effects
  
  greek_learning_modifiers <- 1 + (greek_exp_sim * 0.5)
  french_learning_modifiers <- 1 - (french_exp_sim * 0.3)
  combined_learning_modifiers <- greek_learning_modifiers * french_learning_modifiers
  
  cat("=== COMBINED EXPOSURE ANALYSIS", ifelse(scenario_name != "", paste(":", scenario_name), ""), "===\n")
  cat("Greek: mean =", greek_mean, ", sd =", greek_sd, "\n")
  cat("French: mean =", french_mean, ", sd =", french_sd, "\n")
  if(correlation != 0) cat("Correlation between exposures:", correlation, "\n")
  cat("\nCOMBINED THETA EFFECTS:\n")
  cat("  Combined adjustment range:", round(min(combined_theta_effects), 3), "to", round(max(combined_theta_effects), 3), "\n")
  cat("  Mean combined adjustment:", round(mean(combined_theta_effects), 3), "\n")
  cat("  Net effect direction:", ifelse(mean(combined_theta_effects) > 0, "PRO-DROPPING (Greek wins)", "ANTI-DROPPING (French wins)"), "\n")
  
  cat("\nCOMBINED LEARNING EFFECTS:\n")
  cat("  Combined modifier range:", round(min(combined_learning_modifiers), 3), "to", round(max(combined_learning_modifiers), 3), "\n")
  cat("  Mean combined modifier:", round(mean(combined_learning_modifiers), 3), "\n")
  
  cat("\nCONFLICT ANALYSIS:\n")
  greek_dominates <- (greek_theta_effects > 0) & (abs(greek_theta_effects) > abs(french_theta_effects))
  french_dominates <- (french_theta_effects < 0) & (abs(french_theta_effects) > abs(greek_theta_effects))
  balanced_conflict <- (!greek_dominates) & (!french_dominates)
  
  cat("  Agents where Greek dominates:", round(sum(greek_dominates) / length(greek_dominates) * 100, 1), "%\n")
  cat("  Agents where French dominates:", round(sum(french_dominates) / length(french_dominates) * 100, 1), "%\n")
  cat("  Agents with balanced/conflicting influences:", round(sum(balanced_conflict) / length(balanced_conflict) * 100, 1), "%\n")
  cat("\n")
}

# === SCENARIO TESTING FUNCTION ===
test_exposure_scenarios <- function() {
  cat("==================================================================\n")
  cat("                    EXPOSURE SCENARIOS ANALYSIS                   \n")
  cat("==================================================================\n\n")
  
  # Current scenario (your parameters)
  analyze_greek_exposure_effects(0.4, 0.2, "Current Model")
  analyze_french_exposure_effects(0.6, 0.2, "Current Model")
  analyze_combined_exposure_effects(0.4, 0.2, 0.6, 0.2, 0, "Current Model")
  
  cat("------------------------------------------------------------------\n\n")
  
  # Low French proficiency community
  analyze_greek_exposure_effects(0.6, 0.15, "High Greek Heritage")
  analyze_french_exposure_effects(0.3, 0.15, "Low French Integration")
  analyze_combined_exposure_effects(0.6, 0.15, 0.3, 0.15, 0, "Heritage-Dominant Community")
  
  cat("------------------------------------------------------------------\n\n")
  
  # Highly integrated community
  analyze_greek_exposure_effects(0.2, 0.15, "Low Greek Maintenance")
  analyze_french_exposure_effects(0.8, 0.1, "High French Integration")
  analyze_combined_exposure_effects(0.2, 0.15, 0.8, 0.1, 0, "Highly Integrated Community")
  
  cat("------------------------------------------------------------------\n\n")
  
  # Balanced bilingual community
  analyze_combined_exposure_effects(0.5, 0.2, 0.5, 0.2, -0.3, "Balanced Bilingual (Negative Correlation)")
}

# === PARAMETER SENSITIVITY ANALYSIS ===
sensitivity_analysis <- function() {
  cat("=== PARAMETER SENSITIVITY ANALYSIS ===\n\n")
  
  # Test different Greek effect magnitudes
  cat("Testing Greek theta effect magnitude (currently 0.3):\n")
  for(magnitude in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
    greek_effects <- (0.4 - 0.5) * magnitude  # using mean Greek exposure
    cat("  Magnitude", magnitude, "-> Mean theta effect:", round(greek_effects, 3), "\n")
  }
  
  cat("\nTesting French theta effect magnitude (currently -0.2):\n")
  for(magnitude in c(-0.1, -0.15, -0.2, -0.25, -0.3)) {
    french_effects <- (0.6 - 0.5) * magnitude  # using mean French exposure
    cat("  Magnitude", magnitude, "-> Mean theta effect:", round(french_effects, 3), "\n")
  }
  cat("\n")
}

# Run all analyses
cat("RUNNING EXPOSURE EFFECTS ANALYSIS...\n\n")
test_exposure_scenarios()
sensitivity_analysis()


# Network effects
cat("Final average θ:", mean(agent$theta), "\n")
cat("Network clustering coefficient:", transitivity(g, type = "global"), "\n")
cat("Network density:", edge_density(g), "\n")

# === VISUALIZATIONS ===
# Plot theta evolution over time
plot(1:t_max, totals$theta_avg, type = "l", 
     xlab = "Timestep", ylab = "Average θ", 
     main = "Evolution of Subject Dropping Rate",
     ylim = c(0, 1))
lines(1:t_max, totals$theta_avg + totals$theta_sd, lty = 2, col = "red")
lines(1:t_max, totals$theta_avg - totals$theta_sd, lty = 2, col = "red")
legend("topright", c("Mean", "±1 SD"), lty = c(1, 2), col = c("black", "red"))

plot(1:t_max, totals$v_avg, type = "l", col = "blue",
     xlab = "Timestep", ylab = "Average θ change rate", 
     main = "Evolution of Subject Dropping Rate",
     ylim = c(0, 0.2))

# Final network state
igraph.options(print.id = F)
outgroup_index <- vcount(g)
incident_edges <- incident(g, outgroup_index, mode = "all")
theta_clamped <- pmax(0, pmin(1, agent$theta))  # Ensure values are between 0 and 1
color_indices <- round(theta_clamped * (n_colors - 1)) + 1  # Map to color palette indices

V(g)$color <- ifelse(1:vcount(g) <= N, color_palette[color_indices], "black")
V(g)$size <- 3
V(g)[outgroup_index]$color <- "transparent"  # make outgroup node invisible
V(g)$size <- 10  # default size
V(g)[outgroup_index]$size <- 0
incident_edges <- incident(g, outgroup_index)
E(g)$color <- "black"  # default edge color
E(g)[incident_edges]$color <- "transparent"  # make outgroup links invisible

plot(g, 
     #vertex.size = ifelse(1:vcount(g) <= N, 5 + agent$theta * 10, 12),
     vertex.label = ifelse(1:vcount(g) <= N, round(agent$theta, 2), "OUT"),
     edge.arrow.size = 0,
     main = "Final Network State (colored and sized by θ)")

legend("topright", 
       legend = c("0", "0.5", "1"), 
       col = "black",  # border color
       pt.bg = color_palette[c(1, n_colors/2, n_colors)],  # fill color
       pch = 21,  # filled circle with border
       title = "θ values",
       cex = 0.8)

# Individual trajectories
matplot(t(agent_history[1:min(10, N), , 1]), type = "l", 
        xlab = "Timestep", ylab = "θ", 
        main = "Individual θ Trajectories (first 10 agents)")

cat("\n=== Key Variables Saved ===")
cat("\n- agent: final agent states")
cat("\n- totals: time series data") 
cat("\n- agent_history: full trajectory data")
cat("\n- g: network object")
