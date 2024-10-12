install.packages("igraph")
library(igraph)

install.packages("ggplot2")

# Load necessary library
if (!requireNamespace("mclust", quietly = TRUE)) {
  install.packages("mclust")
}

if (!require("bnlearn")) {
  install.packages("bnlearn", dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("gRain", "RBGL", "gRbase", "Rgraphviz"))

install.packages("mice")

# Load required libraries

library(bnlearn)
library(tidyverse)
library(mice)
library(bnlearn)
library(Rgraphviz)
library(ggplot2)
library(gRain)
library(RBGL)
library(gRbase)
library(readxl)
library(dplyr)
library(tidyr)
library(igraph)
library(mclust)

weather_data <- read.csv("AIMSWaterTempData.csv", header = FALSE, col.names = c("Temperature"))


#1.1
# Plotting the histogram of underwater temperature weather_data
ggplot(weather_data, aes(x = Temperature)) +
  geom_histogram(binwidth = 0.1, fill = "red", color = "yellow", alpha = 0.7) +
  labs(title = "Histogram of Underwater Temperature weather_data",
       x = "Temperature (°C)",
       y = "Frequency") +
  theme_minimal()


# I've already loaded the weather_data into a variable named weather_data
avg_temp <- mean(weather_data$Temperature)
std_temp <- sd(weather_data$Temperature)

# Creating a sequence of temperatures for plotting the Gaussian model
temperature_range <- seq(from = avg_temp - 4 * std_temp, to = avg_temp + 4 * std_temp, length.out = 1000)
# Calculating the density of the normal distribution at these points
density_values <- dnorm(temperature_range, mean = avg_temp, sd = std_temp)

# Plotting the histogram with the density plot
hist(weather_data$Temperature, breaks = 50, probability = TRUE, col = 'red', main = "Temperature Distribution with Gaussian Fit", xlab = "Temperature (°C)", ylab = "Density")
lines(temperature_range, density_values, col = "black", lwd = 2)
#################################

temperature <- weather_data[[1]]  # Extract temperature weather_data

# Checking if temperature weather_data is valid
if (length(temperature) > 0 && all(!is.na(temperature))) {
  # Fit a mixture of Gaussians model with 3 components
  model <- Mclust(temperature, G = 3)
  
  # Summary of the model
  summary(model, parameters = TRUE)
  
  # Extracting parameters for visualization
  parameters <- model$parameters
  means <- parameters$mean
  variances <- parameters$variance$sigmasq
  weights <- parameters$pro
  
  # Printing the parameters
  cat("Mixing Coefficients (weights):", weights, "\n")
  cat("Means (µ):", means, "\n")
  cat("Standard Deviations (σ):", sqrt(variances), "\n")
  
  # Creating a sequence of temperature values for plotting the densities
  temp_seq <- seq(min(temperature), max(temperature), length.out = 1000)
  
  # Calculating the densities for each Gaussian component
  densities <- sapply(1:length(weights), function(i) {
    dnorm(temp_seq, mean = means[i], sd = sqrt(variances[i])) * weights[i]
  })
  
  # Calculating the total density as the sum of individual densities
  total_density <- rowSums(densities)
  
  # Creating a weather_data frame for plotting
  plot_weather_data <- weather_data.frame(
    Temperature = temp_seq,
    Combined = total_density,
    Component1 = densities[, 1],
    Component2 = densities[, 2],
    Component3 = densities[, 3]
  )
  
  # Plotting the histogram and Gaussian densities
  ggplot() +
    geom_histogram(aes(x = temperature, y = ..density..), binwidth = 0.1, fill = "grey", alpha = 0.4) +
    geom_line(weather_data = plot_weather_data, aes(x = Temperature, y = Combined, color = "Combined"), size = 1.2) +
    geom_line(weather_data = plot_weather_data, aes(x = Temperature, y = Component1, color = "Component 1"), size = 1) +
    geom_line(weather_data = plot_weather_data, aes(x = Temperature, y = Component2, color = "Component 2"), size = 1) +
    geom_line(weather_data = plot_weather_data, aes(x = Temperature, y = Component3, color = "Component 3"), size = 1) +
    scale_color_manual(values = c("Combined" = "black", "Component 1" = "yellow", "Component 2" = "blue", "Component 3" = "pink")) +
    labs(title = "Mixture of Gaussians Fit to Temperature weather_data and Density of Gaussian component",
         x = "Temperature (°C)", y = "Density") + theme_minimal()
  
} else {
  cat("Temperature weather_data is empty or contains NA values.")
}

#1.4
# Extracting log-likelihood values from the model
log_likelihood_values <- model$loglik

# Creating a weather_data frame for plotting
likelihood_weather_data <- weather_data.frame(Step = 1:length(log_likelihood_values), LogLikelihood = log_likelihood_values)

# Identifying the maximum log-likelihood value
max_likelihood <- likelihood_weather_data[which.max(likelihood_weather_data$LogLikelihood), ]

# Generating the plot
ggplot(likelihood_weather_data, aes(x = Step, y = LogLikelihood)) +
  geom_line(color = "yellow") +
  geom_point(color = "blue") +
  geom_text(aes(label = ifelse(Step == max_likelihood$Step, "Max", "")),
            vjust = -1) +
  labs(title = "Log-Likelihood Values",
       x = "Step",
       y = "Log Likelihood") +
  theme_minimal()


# Defining states for each node
node_options <- list(
  Temperature = c("low", "medium", "high"),
  Location = c("urban", "rural"),
  HasAccident = c("yes", "no"),
  Weather = c("sunny", "cloudy", "rainy", "foggy"),
  AirQuality = c("good", "poor"),
  AccidentType = c("collision", "skid", "other"),
  Visibility = c("high", "medium", "low"),
  Severity = c("minor", "moderate", "major")
)

# Defining the network structure
nodes_list <- c("Tonnage", "Location", "HasAccident", "Weather", "AirQuality", "AccidentType", "Visibility", "Severity")
arcs_matrix <- matrix(c("", "AirQuality",
                        "Location", "AirQuality",
                        "HasAccident", "AccidentType",
                        "Weather", "Visibility",
                        "AirQuality", "AccidentType",
                        "AirQuality", "Severity",
                        "AccidentType", "Severity",
                        "Visibility", "Severity"), byrow = TRUE, ncol = 2)
network <- graph_from_edgelist(arcs_matrix, directed = TRUE)

# Function to calculate parameters for each node
calculate_parameters <- function(node, network, node_options) {
  parents_nodes <- neighbors(network, node, mode = "in")
  num_states_node <- length(node_options[[node]])
  # Checking if there are parent nodes and calculate the product of their states
  if (length(parents_nodes) > 0) {
    product_parent_states <- prod(sapply(parents_nodes, function(p) length(node_options[[p]])))
  } else {
    product_parent_states <- 1  # If no parent nodes, product of states is 1
  }
  (num_states_node - 1) * product_parent_states
}

# Calculating total parameters required for the network
total_parameters <- sum(sapply(nodes_list, calculate_parameters, network = network, node_options = node_options))

# Printing total parameters
cat("Total parameters:", total_parameters, "\n")

# Setting seed for reproducibility
set.seed(123)

# Generatinging simulated binary data for each variable
num_sam <- 1000  # Number of samples
sim_data <- data.frame(
  Temp = sample(0:1, num_sam, replace = TRUE),
  Loc = sample(0:1, num_sam, replace = TRUE),
  Accident = sample(0:1, num_sam, replace = TRUE),
  Weather = sample(0:1, num_sam, replace = TRUE),
  AirQ = sample(0:1, num_sam, replace = TRUE),
  Type = sample(0:1, num_sam, replace = TRUE),
  Vis = sample(0:1, num_sam, replace = TRUE),
  Severity = sample(0:1, num_sam, replace = TRUE)
)

# Calculating joint probability for the event where all variables are 1
joint_prob <- nrow(sim_data[sim_data$Temp == 1 & sim_data$Loc == 1 & sim_data$Accident == 1 & sim_data$Weather == 1 & sim_data$AirQ == 1 & sim_data$Type == 1 & sim_data$Vis == 1 & sim_data$Severity == 1, ]) / num_sam

# Printing the joint probability
cat("Joint Probability when all variables are 1:", joint_prob, "\n")

# Defining states for each variable
num_states <- c(
  Tonnage = 3, # Tonnage
  Length = 2, # Length:
  HumanFactors = 2, # Human factors
  Weather = 4, # Weather
  KilowattPower = 2, # Kilowatt Power
  AccidentType = 3, # Accident type
  Visibility = 3, # Visibility
  Severity = 3  # Severity of accident
)

# Function to calculate the total number of parameters required
calculate_parameters <- function(states) {
  total_params <- prod(states) - 1
  total_params
}

# Calculating the total number of parameters required
total_params_required <- calculate_parameters(num_states)

# Printing the result
cat("Total number of parameters required without independencies:", total_params_required, "\n")

#2.4A 

# Define the nodes in the network
node_list <- c("Tonnage", "Length", "HumanFactors", "Weather", "KilowattPower", "AccidentType", "Visibility", "Severity")

# Create an empty graph for the Bayesian network
bayesian_network <- empty.graph(node_list)

# Defining the arcs between nodes
arcs_matrix <- matrix(c("Tonnage", "KilowattPower",
                        "HumanFactors", "AccidentType",
                        "Weather", "Visibility",
                        "KilowattPower", "AccidentType",
                        "KilowattPower", "Severity",
                        "AccidentType", "Severity",
                        "Visibility", "Severity"), byrow = TRUE, ncol = 2)

# Set arcs in the Bayesian network
colnames(arcs_matrix) <- c("from", "to")
for (i in 1:nrow(arcs_matrix)) {
  bayesian_network <- set.arc(bayesian_network, from = arcs_matrix[i, "from"], to = arcs_matrix[i, "to"])
}

# Ploting the network to verify it visually
plot(bayesian_network, main = "Bayesian Network")

# Defining the number of states for each variable
states <- list(Tonnage = 3, Location = 2, HasAccident = 2, Weather = 4, AccidentType = 2, AirQuality = 3, Visibility = 3, Severity = 3)

# Calculating parameters for the original network
original_params_AccidentType <- (states$AccidentType - 1) * (states$Tonnage * states$Location)
original_total_params <- original_params_AccidentType  

# Calculating parameters for the modified network
modified_params_AccidentType <- (states$AccidentType - 1) * states$Tonnage
modified_total_params <- modified_params_AccidentType  

# Compute the difference
params_difference <- original_total_params - modified_total_params

# Printing the result
print(paste("Change in the number of parameters:", params_difference))


# Defining the nodes for the Bayesian network
network_nodes <- c("Tonnage", "Length", "HumanFactors", "Weather", "KilowattPower", "AccidentType", "Visibility", "Severity")

# Creating an empty graph with these nodes
bayesian_network <- empty.graph(network_nodes)

# Defining the arcs as per the network diagram
network_arcs <- matrix(c("Tonnage", "KilowattPower",
                         "Length", "KilowattPower",
                         "HumanFactors", "AccidentType",
                         "Weather", "Visibility",
                         "KilowattPower", "AccidentType",
                         "KilowattPower", "Severity",
                         "AccidentType", "Severity",
                         "Visibility", "Severity"), byrow = TRUE, ncol = 2)

# Setting the arcs in the Bayesian network
colnames(network_arcs) <- c("from", "to")
for (i in 1:nrow(network_arcs)) {
  bayesian_network <- set.arc(bayesian_network, from = network_arcs[i, "from"], to = network_arcs[i, "to"])
}

# List of all nodes to check against Tonnage
nodes_to_check <- setdiff(network_nodes, c("Tonnage", "AccidentType", "Visibility"))

# Check conditional independence
independent_nodes <- sapply(nodes_to_check, function(node) {
  dsep(bayesian_network, x = "Tonnage", y = node, z = c("AccidentType", "Visibility"))
})

# Output the nodes that are conditionally independent of Tonnage given AccidentType and Visibility
independent_nodes <- names(independent_nodes[independent_nodes == TRUE])

# Printing the results
cat("Nodes conditionally independent of Tonnage given AccidentType and Visibility:\n")
print(independent_nodes)

# Defining the nodes for the Bayesian network
network_nodes <- c("Tonnage", "Length", "HumanFactors", "Weather", "KilowattPower", "AccidentType", "Visibility", "Severity")

# Creating an empty graph with these nodes
bayesian_network <- empty.graph(network_nodes)

# Defining the arcs as per the network diagram
network_arcs <- matrix(c("Tonnage", "KilowattPower",
                         "Length", "KilowattPower",
                         "HumanFactors", "AccidentType",
                         "Weather", "Visibility",
                         "KilowattPower", "AccidentType",
                         "KilowattPower", "Severity",
                         "AccidentType", "Severity",
                         "Visibility", "Severity"), byrow = TRUE, ncol = 2)

# Setting the arcs in the Bayesian network
colnames(network_arcs) <- c("from", "to")
for (i in 1:nrow(network_arcs)) {
  bayesian_network <- set.arc(bayesian_network, from = network_arcs[i, "from"], to = network_arcs[i, "to"])
}

# Plot the Bayesian network
graphviz.plot(bayesian_network, main = 'Bayesian network using d-seperation')

# d-Separation test for Weather ⊥ KilowattPower | Severity
test_a <- dsep(bayesian_network, x = "Weather", y = "KilowattPower", z = "Severity")
cat("d-separation test for Weather independent of KilowattPower given Severity:", test_a, "\n")

# d-Separation test for HumanFactors ⊥ Visibility | {KilowattPower, AccidentType}
test_b <- dsep(bayesian_network, x = "HumanFactors", y = "Visibility", z = c("KilowattPower", "AccidentType"))
cat("d-separation test for HumanFactors independent of Visibility given KilowattPower and AccidentType:", test_b, "\n")

# Define the nodes for the Bayesian network
network_nodes <- c("Tonnage", "Length", "HumanFactors", "Weather", "KilowattPower", "AccidentType", "Visibility", "Severity")

# Creating an empty graph with these nodes
bayesian_network <- empty.graph(network_nodes)

# Defining the arcs as per the network diagram
network_arcs <- matrix(c("Tonnage", "KilowattPower",
                         "Length", "KilowattPower",
                         "HumanFactors", "AccidentType",
                         "Weather", "Visibility",
                         "KilowattPower", "AccidentType",
                         "KilowattPower", "Severity",
                         "AccidentType", "Severity",
                         "Visibility", "Severity"), byrow = TRUE, ncol = 2)

# Seting the arcs in the Bayesian network
colnames(network_arcs) <- c("from", "to")
for (i in 1:nrow(network_arcs)) {
  bayesian_network <- set.arc(bayesian_network, from = network_arcs[i, "from"], to = network_arcs[i, "to"])
}

# Finding the Markov blanket of KilowattPower
mb_k <- mb(bayesian_network, "KilowattPower")

# Print the nodes in the Markov blanket of KilowattPower
cat("Markov blanket of KilowattPower includes:", paste(mb_k, collapse = ", "), "\n")

# Creating a graphNEL object from the Bayesian network object
graph <- as.graphNEL(bayesian_network)

# Defining node colors
node_colors <- rep("lightblue", length(network_nodes))
names(node_colors) <- network_nodes
node_colors[network_nodes %in% mb_k] <- "red"
node_colors["KilowattPower"] <- "green"

# Setting node attributes for plotting
nodeRenderInfo(graph) <- list(fill = node_colors)

# Plot the Bayesian network with the Markov blanket of KilowattPower
plot(graph, main = "Bayesian Network with Markov Blanket of KilowattPower")



dataset <- data.frame(
  A = c(1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1),
  B = c(0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1),
  C = c(0, 0, 2, 2, 0, 2, 2, 2, 1, 0, 2, 1, 2, 1, 0, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2),
  D = c(0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1),
  E = c(1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0)
)
print(dataset)

# Compute MLE for alpha
alpha_mle <- sum(dataset$A == 0) / nrow(dataset)
cat("Alpha (α):", alpha_mle, "\n")

# Compute MLE for beta
beta_mle <- sum(dataset$E == 0 & dataset$C == 1) / sum(dataset$C == 1)
cat("Beta (β):", beta_mle, "\n")

# Calculate P(E = 0 | A = 0, B = 1) using the formula: 0.36 + 0.5 * beta
P_E0_given_A0_B1 <- 0.36 + 0.5 * beta_mle
cat("P(E = 0 | A = 0, B = 1):", P_E0_given_A0_B1, "\n")



# Define the CPTs based on the provided values
alpha_val <- 0.3
beta_val <- 0.6
lambda_val <- 0.5
gamma_val <- 0.2
sigma_val <- 0.3

# CPT for Var1
cpt_A <- array(c(alpha_val, 1 - alpha_val), dim = 2, dimnames = list(A = c("0", "1")))

# CPT for Var2
cpt_B <- array(c(gamma_val, 1 - gamma_val), dim = 2, dimnames = list(B = c("0", "1")))

# CPT for C given A and B
values_C <- c(0.1, 0.4, 0.5,   # A=0, B=0
              0.4, 0.5, 0.1,   # A=0, B=1
              0.2, lambda_val, 0.8 - lambda_val,   # A=1, B=0
              0.1, 0.1, 0.8)   # Var1=1, Var2=1
cpt_C <- cptable(~C | A:B, values = values_C, levels = c("0", "1", "2"))

# CPT for D given C
values_D <- c(0.7, 0.3,    # C=0
              sigma_val, 1 - sigma_val,    # C=1
              0.4, 0.6)    # Var3=2
cpt_D <- cptable(~D | C, values = values_D, levels = c("0", "1"))

# CPT for E given C
values_E <- c(0.4, 0.6,    # C=0
              beta_val, 1 - beta_val,    # C=1
              0.8, 0.2)    # C=2
cpt_E <- cptable(~E | C, values = values_E, levels = c("0", "1"))

# Compile the CPTs into a belief network
prob_list <- compileCPT(list(cpt_A, cpt_B, cpt_C, cpt_D, cpt_E))

# Create the belief network
belief_network <- grain(prob_list)

# Reset the graphics device
dev.off()

# Plot the belief network
plot(belief_network)

# Display the probability tables
cptA <- querygrain(belief_network, nodes = "A")
cptB <- querygrain(belief_network, nodes = "B")
cptC <- querygrain(belief_network, nodes = "C")
cptD <- querygrain(belief_network, nodes = "D")
cptE <- querygrain(belief_network, nodes = "E")

# Print CPTs
cat("CPT for A:\n")
print(cptA)
cat("CPT for B:\n")
print(cptB)
cat("CPT for C:\n")
print(cptC)
cat("CPT for D:\n")
print(cptD)
cat("CPT for E:\n")
print(cptE)

# Define the probability tables with different values
cpt_A2 <- cptable(~A, values=c(0.67, 0.33), levels=c("0", "1"))
cpt_B2 <- cptable(~B, values=c(0.67, 0.33), levels=c("0", "1"))
cpt_C2 <- cptable(~C | A + B, values=c(0.1, 0.9, 0.2, 0.8, 0.3, 0.7, 0.4, 0.6), levels=c("0", "1"))
cpt_D2 <- cptable(~D | C, values=c(0.5, 0.5, 0.6, 0.4), levels=c("0", "1"))
cpt_E2 <- cptable(~E | D, values=c(0.7, 0.3, 0.4, 0.6), levels=c("0", "1"))

# Compile the new network
cpt_list2 <- compileCPT(list(cpt_A2, cpt_B2, cpt_C2, cpt_D2, cpt_E2))
belief_network2 <- grain(cpt_list2)

# i) Marginal distribution of C
marginal_C <- querygrain(belief_network, nodes = "C")
cat("Marginal distribution of C:\n")
print(marginal_C)

# ii) Joint distribution of A, B, and C
joint_ABC <- querygrain(belief_network, nodes = c("A", "B", "C"), type = "joint")
cat("Joint distribution of A, B, and C:\n")
print(joint_ABC)

# iii) Conditional probability P(C=2 | A=0, D=1, E=0)
# Set evidence for A, D, and E
belief_network <- setEvidence(belief_network, evidence = list(A = "0", D = "1", E = "0"))

# Compute the conditional probability P(C=2 | A=0, D=1, E=0)
conditional_C_given_ADE <- querygrain(belief_network, nodes = "C")[["C"]]["2"]
cat("Conditional probability P(C=2 | A=0, D=1, E=0):", conditional_C_given_ADE, "\n")
##########################



# Loading the alarm dataset
data("alarm")

# Function to learn the network and plot it
learn_and_plot <- function(data, score_type, title) {
  bn <- hc(data, score = score_type)
  bn_score <- score(bn, data, type = score_type)
  plot_title <- paste(title, "score:", round(bn_score, 2))
  graphviz.plot(bn, main = plot_title, shape = "ellipse")
  return(bn_score)
}

# Subset the data for different sample sizes
data_500 <- alarm[1:500, ]
data_5000 <- alarm[1:5000, ]
data_10000 <- alarm[1:10000, ]

# Learn and plot networks for sample size 500
score_500_BIC <- learn_and_plot(data_500, "bic", "500 samples, BIC")
score_500_BDe <- learn_and_plot(data_500, "bde", "500 samples, BDe")

# Learn and plot networks for sample size 5000
score_5000_BIC <- learn_and_plot(data_5000, "bic", "5000 samples, BIC")
score_5000_BDe <- learn_and_plot(data_5000, "bde", "5000 samples, BDe")

# Learn and plot networks for sample size 10000
score_10000_BIC <- learn_and_plot(data_10000, "bic", "10000 samples, BIC")
score_10000_BDe <- learn_and_plot(data_10000, "bde", "10000 samples, BDe")

#4.2 Print the scores obtained
scores <- data.frame(
  Sample_Size = c(500, 5000, 10000),
  BIC_Score = c(score_500_BIC, score_5000_BIC, score_10000_BIC),
  BDe_Score = c(score_500_BDe, score_5000_BDe, score_10000_BDe)
)

print(scores)


# Use the full dataset
full_data <- alarm

# Learn and plot networks for the full dataset using BIC score
score_full_BIC <- learn_and_plot(full_data, "bic", "Full dataset, BIC")

# Learn and plot networks for the full dataset using BDe score
score_full_BDe <- learn_and_plot(full_data, "bde", "Full dataset, BDe")

# Print the scores obtained
scores <- data.frame(
  Dataset = c("Full dataset"),
  BIC_Score = c(score_full_BIC),
  BDe_Score = c(score_full_BDe)
)

print(scores)



# Define the true network structure
modelstring = paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
                     "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                     "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                     "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT]",
                     "[PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB]",
                     "[DISC][MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB]",
                     "[VALV|INT:VLNG][ACO2|VALV][CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL]",
                     "[CO|HR:STKV][BP|CO:TPR]")
true_dag <- model2network(modelstring)

# Check the class of true_dag
print(class(true_dag))  # Should be "bn"

# Function to learn the network
learn_network <- function(data, score_type) {
  bn <- hc(data, score = score_type)
  return(bn)
}

# Learn networks for the full dataset using BIC and BDe scores
bn_full_BIC <- learn_network(full_data, "bic")
bn_full_BDe <- learn_network(full_data, "bde")

# Check the class of the learned networks
print(class(bn_full_BIC))  # Should be "bn"
print(class(bn_full_BDe))  # Should be "bn"

# Compare the learned networks with the true network structure
comparison_BIC <- compare(bn_full_BIC, true_dag)
comparison_BDe <- compare(bn_full_BDe, true_dag)

# Print the comparison results
print("Comparison with BIC scoring method:")
print(comparison_BIC)

print("Comparison with BDe scoring method:")
print(comparison_BDe)


# Visualize the comparisons
graphviz.compare(bn_full_BIC, true_dag, main = c("BIC Learned Network", "True Network"))
graphviz.compare(bn_full_BDe, true_dag, main = c("BDe Learned Network", "True Network"))

# Fit the data to the learned network
fitted_bn <- bn.fit(bn_full_BIC, full_data)

#Extract and display the CPD table for the variable "HR"
cpd_hr <- fitted_bn$HR
print(cpd_hr)


# Query the conditional probability P(HR = "HIGH" | BP = "LOW", PRSS = "NORMAL")
query_result <- cpquery(fitted_bn, event = (HR == "HIGH"), evidence = (BP == "LOW" & PRSS == "NORMAL"))
print(query_result)




# Load the data
df <- read.csv("Book1.csv")

# Remove rows with any missing values
df <- na.omit(df)

# Perform multiple imputation
imputed_data <- mice(df, m = 1, method = 'pmm', seed = 500)
df <- complete(imputed_data, 1)


# Display the structure of the dataset
str(df)

colnames(df)

# Summary statistics for numeric variables
summary(df[, c("Speed.Limit", "Vehicle.Year", "Wind.Speed..MPH.", "TMP10.DK", "Relative.Humidity....", "Rain..Inches.")])

# Bar plot for categorical variables
categorical_vars <- c("ACRS.Report.Type", "Route.Type", "Collision.Type", "Weather", "Surface.Condition", "Light", "Traffic.Control", "Driver.At.Fault", "Injury.Severity", "Driver.Distracted.By", "Vehicle.Body.Type", "Vehicle.Movement", "Vehicle.Continuing.Dir", "Vehicle.Going.Dir", "Vehicle.Make", "Municipality")
for(var in categorical_vars) {
  print(
    ggplot(df, aes_string(x = var)) +
      geom_bar() +
      labs(title = var, x = var, y = "Frequency") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}

# Histogram plots for numeric variables
numeric_vars <- c("Speed.Limit", "Wind.Speed..MPH.", "TMP10.DK", "Relative.Humidity....", "Rain..Inches.")
for(var in numeric_vars) {
  print(
    ggplot(df, aes_string(x = var)) +
      geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
      labs(title = var, x = var, y = "Frequency") +
      theme_minimal()
  )
}

# Summary statistics for Weather and Surface Condition
summary_weather <- table(df$Weather)
summary_surface <- table(df$Surface.Condition) 

# Bar plot for Weather
barplot(summary_weather, main = "Weather Conditions", xlab = "Weather", ylab = "Number of Accidents", col = "skyblue")

# Bar plot for Surface Condition
barplot(summary_surface, main = "Surface Conditions", xlab = "Surface Condition", ylab = "Number of Accidents", col = "lightgreen")

# Cross-tabulation of Weather and Surface Condition with Accident Count
cross_tab <- table(df$Weather, df$Surface.Condition)

# Heatmap for cross-tabulation
heatmap(cross_tab, 
        Rowv = NA, Colv = NA, 
        col = terrain.colors(length(unique(cross_tab))), 
        main = "Accidents by Weather and Surface Condition",
        xlab = "Surface Condition", ylab = "Weather")


cross_tab1 <- table(df$Traffic.Control, df$Driver.At.Fault)

my_colors <- colorRampPalette(c("blue", "white", "red"))(length(unique(cross_tab1)))
# Create the heatmap with the color key
heatmap(cross_tab1, 
        Rowv = NA, Colv = NA, 
        col = my_colors, 
        main = "Accidents by Traffic control and drivers fault",
        xlab = "Traffic Control", ylab = "Drivers Fault",
        key = TRUE)

# Convert all columns to factor variables using lapply
df[] <- lapply(df, as.factor)


# Structure Learning using Hill-Climbing algorithm
hc_result <- hc(df)

# Parameter Learning
bn_fit <- bn.fit(hc_result, data = df)

# Print the learned network
print(bn_fit)

install.packages("graphviz")

# Plot the Bayesian network structure
graphviz.plot(hc_result)

#Tabu BN Algorithm
bn_tabu <- tabu(df, score = "bde", max.id = 2)

bn_tabu_mle <- bn.fit(bn_tabu, data = df, method = "mle")

print(bn_tabu)

graphviz.plot(bn_tabu)

#MMHC BN Algorithm
mmhc_result <- mmhc(df)

bn_fit_mmhc <- bn.fit(mmhc_result, data = df)

print(bn_fit_mmhc)
graphviz.plot(bn_fit_mmhc)

# Convert the bnlearn object to a gRain object
bn_grain <- as.grain(bn_fit)


# Set the evidence
evidence <- list("Surface.Condition" = "WET")

# Set evidence in the Bayesian network
bn_grain_with_evidence <- setEvidence(bn_grain, nodes = names(evidence), states = unlist(evidence))

# Query the conditional probability
query_result <- querygrain(bn_grain_with_evidence, nodes = "Injury.Severity")

# Print the result
print(query_result$`Injury.Severity`)


# Test for conditional independence
independence_test <- ci.test("Injury.Severity", "Speed.Limit", data = df)

# Print the result
print(independence_test)

#  Finding the Markov blanket of "Driver.Distracted.By"
markov_blanket <- mb(bn_fit, "Driver.Distracted.By")
print(markov_blanket)

############################END##############################


