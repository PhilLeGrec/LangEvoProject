# Packages
library(igraph)

# === PARAMETERS ===
N <- 30  # number of agents
t_max <- 200  # number of timesteps
network_type <- "small-world"  # options: "random", "small-world", "regular", "fully-connected", "data-based"
p_link <- 0.04  # probability of connection (for random/small-world)
include_outgroup <- TRUE  # include outgroup connections
agent_turnover <- FALSE  # implement agent turnover over time
base_interaction_rate <- 15  # increased from 5 to compensate for probabilistic interactions

# === INITIALIZATION ===
set.seed(65328)  # for reproducibility

## Initialize agents with linguistic variables
agent <- data.frame(
  id = 1:N,
  generation = sample(1:3, N, replace = TRUE, prob = c(0.8, 0.15, 0.05)), # 1=1st gen, 2=2nd gen, 3=3rd gen
  init_exp_greek = pmax(0, pmin(1, rnorm(N, mean = 0.4, sd = 0.2))), # exposure to Greek, Gaussian
  init_exp_french = pmax(0, pmin(1, rnorm(N, mean = 0.6, sd = 0.2))), # exposure to French, slightly higher mean of Gaussian
  # Dynamic exposures (will change over time)
  current_exp_greek = rep(NA, N),
  current_exp_french = rep(NA, N),
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

# Initialize current exposures with initial values
agent$current_exp_greek <- agent$init_exp_greek
agent$current_exp_french <- agent$init_exp_french

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

## Apply Greek exposure effects (now using current exposure)
greek_exposure_effects <- function(agent_df) {
  for(i in 1:nrow(agent_df)) {
    greek_exp <- agent_df$current_exp_greek[i]
    
    # Greek exposure increases theta (more subject dropping)
    greek_theta_effect <- (greek_exp - 0.5) * 0.3  # centered at 0.5, scaled by 0.3
    
    # Greek exposure affects learning rate (higher exposure = more sensitive to Greek-like patterns)
    agent_df$greek_learning_modifier[i] <- 1 + (greek_exp * 0.5)  # 1.0 to 1.5 range
  }
  
  return(agent_df)
}

## Apply French exposure effects (now using current exposure, reduced initial effect)
french_exposure_effects <- function(agent_df) {
  for(i in 1:nrow(agent_df)) {
    french_exp <- agent_df$current_exp_french[i]
    
    # French exposure decreases theta (less subject dropping) - reduced from -0.2 to -0.1
    french_theta_effect <- (french_exp - 0.5) * -0.1  # reduced initial effect
    
    # French exposure makes agents more resistant to adopting dropping behavior
    agent_df$french_learning_modifier[i] <- 1 - (french_exp * 0.3)  # 0.7 to 1.0 range
  }
  
  return(agent_df)
}

## Apply all effects in sequence
agent <- generation_effects(agent)
agent <- greek_exposure_effects(agent)
agent <- french_exposure_effects(agent)

## Calculate initial theta adjustments
for(i in 1:N) {
  greek_theta_effect <- (agent$current_exp_greek[i] - 0.5) * 0.3
  french_theta_effect <- (agent$current_exp_french[i] - 0.5) * -0.1  # reduced effect
  agent$theta_0[i] <- pmax(0, pmin(1, agent$theta_0[i] + greek_theta_effect + french_theta_effect))
}

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

# Add edge weights using lognormal distribution (bounded between 0 and 1)
# First, get all edges excluding outgroup connections
all_edges <- E(g)
if(include_outgroup) {
  outgroup_id <- vcount(g)
  outgroup_edges <- incident(g, outgroup_id, mode = "all")
  ingroup_edges <- setdiff(all_edges, outgroup_edges)
} else {
  ingroup_edges <- all_edges
  outgroup_edges <- c()
}

# Assign weights to ingroup connections
if(length(ingroup_edges) > 0) {
  weights_raw <- rlnorm(length(ingroup_edges), meanlog = -0.5, sdlog = 0.8)
  # Normalize to 0-1 range but keep some variation
  weights_normalized <- pmax(0.1, pmin(1.0, weights_raw / quantile(weights_raw, 0.9)))
  E(g)[ingroup_edges]$weight <- weights_normalized
}

# Assign weights to outgroup connections (if they exist)
if(include_outgroup && length(outgroup_edges) > 0) {
  outgroup_weights_raw <- rlnorm(length(outgroup_edges), meanlog = -1.2, sdlog = 0.6)
  outgroup_weights_normalized <- pmax(0.05, pmin(0.5, outgroup_weights_raw / quantile(outgroup_weights_raw, 0.9)))
  E(g)[outgroup_edges]$weight <- outgroup_weights_normalized
}

# Function to update exposures based on network structure
update_exposures <- function(agent_df, graph, timestep) {
  N_agents <- nrow(agent_df)
  
  for(i in 1:N_agents) {
    # Calculate weighted in-degree centrality (excluding outgroup)
    neighbors <- neighbors(graph, i)
    ingroup_neighbors <- neighbors[neighbors <= N_agents]
    
    if(length(ingroup_neighbors) > 0) {
      weighted_indegree <- 0
      for(j in ingroup_neighbors) {
        # Check if edge exists and get weight safely
        if(are.connected(graph, i, j)) {
          edge_id <- get.edge.ids(graph, c(i, j))
          if(!is.na(edge_id) && edge_id > 0) {
            weight <- E(graph)$weight[edge_id]
            weighted_indegree <- weighted_indegree + weight
          }
        }
      }
      
      # Normalize by maximum possible weighted in-degree
      max_possible_indegree <- length(ingroup_neighbors)  # if all weights were 1
      if(max_possible_indegree > 0) {
        normalized_indegree <- weighted_indegree / max_possible_indegree
        
        # Update Greek exposure: higher centrality -> maintain/increase Greek exposure
        greek_change_rate <- 0.002  # small change per timestep
        centrality_effect <- normalized_indegree - 0.5  # centered effect
        agent_df$current_exp_greek[i] <- pmax(0, pmin(1, 
                                                      agent_df$current_exp_greek[i] + greek_change_rate * centrality_effect))
      }
    }
    
    # Update French exposure based on outgroup connection
    if(include_outgroup && vcount(graph) > N_agents) {
      outgroup_id <- vcount(graph)
      if(are.connected(graph, i, outgroup_id)) {
        edge_id <- get.edge.ids(graph, c(i, outgroup_id))
        if(!is.na(edge_id) && edge_id > 0) {
          outgroup_weight <- E(graph)$weight[edge_id]
          
          # French exposure increases over time, accelerated by outgroup connection strength
          french_base_increase <- 0.001  # base increase per timestep
          outgroup_effect <- outgroup_weight * 0.002  # additional increase based on connection
          
          agent_df$current_exp_french[i] <- pmax(0, pmin(1,
                                                         agent_df$current_exp_french[i] + french_base_increase + outgroup_effect))
        } else {
          # Still some increase even with problematic outgroup connection
          agent_df$current_exp_french[i] <- pmax(0, pmin(1,
                                                         agent_df$current_exp_french[i] + 0.0005))
        }
      } else {
        # Still some increase even without outgroup connection (societal influence)
        agent_df$current_exp_french[i] <- pmax(0, pmin(1,
                                                       agent_df$current_exp_french[i] + 0.0005))
      }
    }
  }
  
  return(agent_df)
}

# Function to recalculate exposure effects
recalculate_exposure_effects <- function(agent_df) {
  # Recalculate Greek and French learning modifiers based on current exposures
  for(i in 1:nrow(agent_df)) {
    greek_exp <- agent_df$current_exp_greek[i]
    french_exp <- agent_df$current_exp_french[i]
    
    # Update modifiers
    agent_df$greek_learning_modifier[i] <- 1 + (greek_exp * 0.5)
    agent_df$french_learning_modifier[i] <- 1 - (french_exp * 0.3)
    
    # Recalculate total learning modifier
    agent_df$total_learning_modifier[i] <- agent_df$generation_learning_modifier[i] * 
      agent_df$greek_learning_modifier[i] * 
      agent_df$french_learning_modifier[i]
  }
  
  return(agent_df)
}

# Visualize initial network with proper edge styling
igraph.options(print.id = F)

# Prepare node colors based on theta values
n_colors <- 100
color_palette <- heat.colors(n_colors)
theta_clamped <- pmax(0, pmin(1, agent$theta))
color_indices <- round(theta_clamped * (n_colors - 1)) + 1

# Set node properties
V(g)$color <- ifelse(1:vcount(g) <= N, color_palette[color_indices], "gray")
V(g)$size <- 15  # Make nodes larger to see better
V(g)$label <- ifelse(1:vcount(g) <= N, 1:N, "OUT")
V(g)$label.cex <- 0.8
V(g)$label.color <- "black"
V(g)$shape <- "circle"  # Set default shape for all vertices

# If outgroup exists, style it differently
if(include_outgroup) {
  outgroup_id <- vcount(g)
  V(g)[outgroup_id]$color <- "cyan"
  V(g)[outgroup_id]$size <- 20
  V(g)[outgroup_id]$shape <- "square"
}

# Set edge properties based on weights
edge_weights <- E(g)$weight
edge_weights[is.na(edge_weights)] <- 0.1  # Handle any remaining NAs

# Edge colors: darker = stronger connection
edge_colors <- gray(1 - edge_weights)  # Higher weight = darker color
E(g)$color <- edge_colors

# Edge widths: thicker = stronger connection  
E(g)$width <- edge_weights * 4 + 0.5  # Scale to reasonable width range

# Set default line type for all edges
E(g)$lty <- 1  # Solid lines as default

# If outgroup exists, style outgroup edges differently
if(include_outgroup) {
  outgroup_edges <- incident(g, outgroup_id, mode = "all")
  E(g)[outgroup_edges]$color <- "lightskyblue1"
  E(g)[outgroup_edges]$lty <- 2  # Dashed lines for outgroup connections
}

# Create the plot
plot(g, 
     vertex.label = V(g)$label,
     edge.arrow.size = 0,
     layout = layout_with_fr(g),  # Force-directed layout
     main = paste("Initial Network:", network_type, "with weighted connections"),
     sub = "Node color = θ value, Edge darkness/width = connection strength")

# Add legends
legend("topright", 
       legend = c("θ = 0", "θ = 0.5", "θ = 1"), 
       col = "black",
       pt.bg = color_palette[c(1, n_colors/2, n_colors)],
       pch = 21,
       title = "Subject Dropping",
       cex = 0.7)

legend("topleft",
       legend = c("Strong (1.0)", "Medium (0.5)", "Weak (0.1)"),
       col = gray(c(0, 0.5, 0.9)),
       lwd = c(4.5, 2.5, 0.6),
       title = "Connection Strength",
       cex = 0.7)

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

# Function to simulate one timestep with probabilistic interactions
simulate_timestep <- function(agent_df, graph, h) {
  N_agents <- nrow(agent_df)
  new_agent_df <- agent_df
  
  # Update exposures based on network structure and time
  new_agent_df <- update_exposures(new_agent_df, graph, h)
  
  # Recalculate exposure effects
  new_agent_df <- recalculate_exposure_effects(new_agent_df)
  
  # Store previous theta for velocity calculation
  prev_theta <- new_agent_df$theta
  
  # For each agent
  for(i in 1:N_agents) {
    # Reset interaction counters for this timestep
    n_hi <- 0
    k_hi <- 0
    
    # Get neighbors (excluding outgroup for learning)
    neighbors <- neighbors(graph, i)
    neighbors <- neighbors[neighbors <= N_agents]  # exclude outgroup
    
    if(length(neighbors) > 0) {
      # Calculate total available interaction budget
      total_budget <- base_interaction_rate
      
      # Distribute interactions probabilistically among neighbors
      neighbor_weights <- numeric(length(neighbors))
      for(j_idx in 1:length(neighbors)) {
        j <- neighbors[j_idx]
        # Check if edge exists and get weight safely
        if(are.connected(graph, i, j)) {
          edge_id <- get.edge.ids(graph, c(i, j))
          if(!is.na(edge_id) && edge_id > 0) {
            neighbor_weights[j_idx] <- E(graph)$weight[edge_id]
          } else {
            neighbor_weights[j_idx] <- 0.1  # default small weight
          }
        } else {
          neighbor_weights[j_idx] <- 0.1  # default small weight
        }
        
        # Ensure weight is not NA or negative
        if(is.na(neighbor_weights[j_idx]) || neighbor_weights[j_idx] < 0) {
          neighbor_weights[j_idx] <- 0.1
        }
      }
      
      # Normalize weights to probabilities
      if(sum(neighbor_weights, na.rm = TRUE) > 0) {
        # Remove any NA values from neighbor_weights
        neighbor_weights[is.na(neighbor_weights)] <- 0.1
        neighbor_probs <- neighbor_weights / sum(neighbor_weights, na.rm = TRUE)
        
        # Distribute total interactions among neighbors
        for(j_idx in 1:length(neighbors)) {
          j <- neighbors[j_idx]
          weight <- neighbor_weights[j_idx]
          
          # Probability of interaction at each "slot" (ensure valid probability)
          interaction_prob <- pmin(1, pmax(0, neighbor_probs[j_idx]))
          
          if(!is.na(interaction_prob) && interaction_prob > 0) {
            # Sample number of interactions (binomial with adjusted probability)
            n_interactions_with_j <- rbinom(1, total_budget, interaction_prob)
            
            if(n_interactions_with_j > 0) {
              # Read neighbor's linguistic behavior
              n_hij <- new_agent_df$n[j] * (n_interactions_with_j / total_budget)
              k_hij <- new_agent_df$k[j] * (n_interactions_with_j / total_budget)
              
              # Sum into agent i's observations
              n_hi <- n_hi + n_hij
              k_hi <- k_hi + k_hij
            }
          }
        }
      }
      
      # Update theta based on neighbor observations
      if(n_hi > 0) {
        new_agent_df$theta[i] <- bayesian_update(new_agent_df$theta[i],
                                                 k_hi, n_hi, learning_modifier = new_agent_df$total_learning_modifier[i])
      }
    }
    
    # Handle outgroup influence if connected (with dynamic French effect)
    if(include_outgroup && vcount(graph) > N_agents) {
      outgroup_id <- vcount(graph)
      if(are.connected(graph, i, outgroup_id)) {
        edge_id <- get.edge.ids(graph, c(i, outgroup_id))
        if(!is.na(edge_id) && edge_id > 0) {
          outgroup_weight <- E(graph)$weight[edge_id]
          
          # Outgroup influence now depends on current French exposure
          outgroup_theta <- 0.1 * (1 - new_agent_df$current_exp_french[i])  # less dropping with more French
          outgroup_n <- 5
          outgroup_k <- outgroup_n * outgroup_theta
          
          # Apply outgroup influence (weighted, with probabilistic interaction)
          if(!is.na(outgroup_weight) && outgroup_weight > 0 && 
             rbinom(1, 1, pmin(1, outgroup_weight * 0.5)) == 1) {  # probability of outgroup interaction
            new_agent_df$theta[i] <- bayesian_update(new_agent_df$theta[i], 
                                                     outgroup_k * outgroup_weight, 
                                                     outgroup_n * outgroup_weight,
                                                     learning_modifier = new_agent_df$total_learning_modifier[i])
          }
        }
      }
    }
  }
  
  # Generate new interactions for next timestep based on updated theta
  for(i in 1:N_agents) {
    neighbors <- neighbors(graph, i)
    neighbors <- neighbors[neighbors <= N_agents]
    
    new_n_total <- 0
    new_k_total <- 0
    
    if(length(neighbors) > 0) {
      # Each neighbor contributes to i's linguistic data
      for(j in neighbors) {
        # Check if edge exists and get weight safely
        if(are.connected(graph, i, j)) {
          edge_id <- get.edge.ids(graph, c(i, j))
          if(!is.na(edge_id) && edge_id > 0) {
            weight <- E(graph)$weight[edge_id]
          } else {
            weight <- 0.1  # default small weight
          }
        } else {
          weight <- 0.1  # default small weight
        }
        
        # Ensure weight is valid
        if(is.na(weight) || weight < 0) {
          weight <- 0.1
        }
        
        # Probabilistic interaction (ensure weight is valid)
        if(!is.na(weight) && weight > 0 && rbinom(1, 1, pmin(1, weight)) == 1) {
          # Generate interactions based on current theta and connection weight
          n_interactions <- rpois(1, base_interaction_rate * weight * 0.5)
          k_drops <- rbinom(1, n_interactions, new_agent_df$theta[j])  # observe j's behavior
          
          new_n_total <- new_n_total + n_interactions
          new_k_total <- new_k_total + k_drops
        }
      }
    }
    
    # Update with exponential decay of old observations
    decay_factor <- 0.9
    new_agent_df$n[i] <- new_agent_df$n[i] * decay_factor + new_n_total
    new_agent_df$k[i] <- new_agent_df$k[i] * decay_factor + new_k_total
  }
  
  # Calculate velocity (change in theta)
  new_agent_df$v <- new_agent_df$theta - prev_theta
  
  return(new_agent_df)
}

# === TRACKING VARIABLES ===
totals <- data.frame(
  timestep = 1:t_max,
  k_total = rep(NA, t_max),
  theta_avg = rep(NA, t_max),
  theta_sd = rep(NA, t_max),
  v_avg = rep(NA, t_max),
  clustering = rep(NA, t_max),
  greek_exp_avg = rep(NA, t_max),
  french_exp_avg = rep(NA, t_max)
)

# Store agent history for analysis
agent_history <- array(NA, dim = c(N, t_max, 6))  # agents x time x variables (theta, k, n, v, greek_exp, french_exp)

# === MAIN SIMULATION LOOP ===
cat("Running simulation for", t_max, "timesteps...\n")

for(h in 1:t_max) {
  # Simulate one timestep
  agent <- simulate_timestep(agent, g, h)
  
  # Calculate totals
  totals$k_total[h] <- sum(agent$k, na.rm = TRUE)
  totals$theta_avg[h] <- mean(agent$theta, na.rm = TRUE)
  totals$theta_sd[h] <- sd(agent$theta, na.rm = TRUE)
  totals$v_avg[h] <- mean(abs(agent$v), na.rm = TRUE)
  totals$clustering[h] <- transitivity(g, type = "global")
  totals$greek_exp_avg[h] <- mean(agent$current_exp_greek, na.rm = TRUE)
  totals$french_exp_avg[h] <- mean(agent$current_exp_french, na.rm = TRUE)
  
  # Store agent states
  agent_history[, h, 1] <- agent$theta
  agent_history[, h, 2] <- agent$k
  agent_history[, h, 3] <- agent$n
  agent_history[, h, 4] <- agent$v
  agent_history[, h, 5] <- agent$current_exp_greek
  agent_history[, h, 6] <- agent$current_exp_french
  
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

# Exposure correlations
cat("Correlation θ(final) ~ final Greek exposure:", cor(agent$theta, agent$current_exp_greek), "\n")
cat("Correlation θ(final) ~ final French exposure:", cor(agent$theta, agent$current_exp_french), "\n")

# Exposure changes
cat("Greek exposure change: ", mean(agent$current_exp_greek - agent$init_exp_greek), "\n")
cat("French exposure change: ", mean(agent$current_exp_french - agent$init_exp_french), "\n")

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

# Plot velocity over time
v_avg_clean <- totals$v_avg
v_avg_clean[is.na(v_avg_clean) | is.infinite(v_avg_clean)] <- 0
plot(1:t_max, v_avg_clean, type = "l", col = "blue",
     xlab = "Timestep", ylab = "Average |θ change|", 
     main = "Rate of Change in Subject Dropping",
     ylim = c(0, max(v_avg_clean, na.rm = TRUE) * 1.1))

# Plot exposure evolution over time
plot(1:t_max, totals$greek_exp_avg, type = "l", col = "darkgreen",
     xlab = "Timestep", ylab = "Average Exposure", 
     main = "Evolution of Language Exposures",
     ylim = c(0, 1))
lines(1:t_max, totals$french_exp_avg, type = "l", col = "darkblue")
legend("right", c("Greek", "French"), lty = c(1, 1), col = c("darkgreen", "darkblue"))

# Final network state
theta_clamped <- pmax(0, pmin(1, agent$theta))
color_indices <- round(theta_clamped * (n_colors - 1)) + 1

# Update node colors for final state
V(g)$color <- ifelse(1:vcount(g) <= N, color_palette[color_indices], "cyan")

# Keep edge styling the same (weights don't change during simulation)
# Edge colors and widths are already set based on weights

plot(g, 
     vertex.label = V(g)$label,
     edge.arrow.size = 0,
     layout = layout_with_fr(g),
     main = "Final Network State",
     sub = "Node color = final θ value, Edge darkness/width = connection strength")

legend("topright", 
       legend = c("θ = 0", "θ = 0.5", "θ = 1"), 
       col = "black",
       pt.bg = color_palette[c(1, n_colors/2, n_colors)],
       pch = 21,
       title = "Subject Dropping",
       cex = 0.6)

legend("topleft",
       legend = c("Strong (1.0)", "Medium (0.5)", "Weak (0.1)"),
       col = gray(c(0, 0.5, 0.9)),
       lwd = c(4.5, 2.5, 0.6),
       title = "Connection Strength",
       cex = 0.6)

# Individual trajectories plot
theta_trajectories <- t(agent_history[1:min(15, N), , 1])
# Replace NA values with the last known value or 0
theta_trajectories[is.na(theta_trajectories)] <- 0
# Ensure we have valid y-limits
ylim_range <- range(theta_trajectories, na.rm = TRUE, finite = TRUE)
if(any(is.infinite(ylim_range))) ylim_range <- c(0, 1)

matplot(theta_trajectories, type = "l", 
        xlab = "Timestep", ylab = "θ", 
        main = "Individual θ Trajectories (first 15 agents)",
        ylim = ylim_range)

# Additional visualization: Exposure trajectories for sample agents

# Greek exposure trajectories
greek_trajectories <- t(agent_history[1:10, , 5])
greek_trajectories[is.na(greek_trajectories)] <- 0
greek_ylim <- range(greek_trajectories, na.rm = TRUE, finite = TRUE)
if(any(is.infinite(greek_ylim))) greek_ylim <- c(0, 1)

matplot(greek_trajectories, type = "l", 
        xlab = "Timestep", ylab = "Greek Exposure", 
        main = "Greek Exposure Trajectories (first 10 agents)",
        ylim = greek_ylim)

# French exposure trajectories  
french_trajectories <- t(agent_history[1:min(10, N), , 6])
french_trajectories[is.na(french_trajectories)] <- 0
french_ylim <- range(french_trajectories, na.rm = TRUE, finite = TRUE)
if(any(is.infinite(french_ylim))) french_ylim <- c(0, 1)

matplot(french_trajectories, type = "l", 
        xlab = "Timestep", ylab = "French Exposure", 
        main = "French Exposure Trajectories (first 10 agents)",
        ylim = french_ylim)

par(mfrow = c(1, 1))

cat("\n=== Key Variables Saved ===")
cat("\n- agent: final agent states with dynamic exposures")
cat("\n- totals: time series data including exposure averages") 
cat("\n- agent_history: full trajectory data including exposures")
cat("\n- g: network object with weighted edges")