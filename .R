# Load necessary libraries
library(bnlearn)
library(gRain)

# Load and process the dataset
# Assuming 'df' is your dataset, convert all columns to factors for Bayesian network analysis
df[] <- lapply(df, as.factor)

# Structure Learning using the Hill-Climbing algorithm
# This learns the structure of a Bayesian Network based on data.
hc_result <- hc(df)

# Parameter Learning: Fit the learned structure with the data
# This step estimates the probabilities (parameters) for each node based on the data.
bn_fit <- bn.fit(hc_result, data = df)

# Print the learned Bayesian Network structure with its parameters
print(bn_fit)

# Plot the learned Bayesian Network structure
graphviz.plot(hc_result)

# Apply the Tabu Search algorithm for structure learning with "bde" scoring
# Tabu search is another method for learning the network structure, avoiding local minima.
bn_tabu <- tabu(df, score = "bde", max.id = 2)

# Fit the network with parameters using Maximum Likelihood Estimation (MLE)
bn_tabu_mle <- bn.fit(bn_tabu, data = df, method = "mle")

# Print the Tabu search results
print(bn_tabu)

# Plot the Bayesian Network learned by Tabu search
graphviz.plot(bn_tabu)

# Structure learning using the MMHC (Max-Min Hill Climbing) algorithm
mmhc_result <- mmhc(df)

# Fit the network structure learned by MMHC with the data
bn_fit_mmhc <- bn.fit(mmhc_result, data = df)

# Print the MMHC algorithm results
print(bn_fit_mmhc)

# Plot the MMHC Bayesian Network structure
graphviz.plot(bn_fit_mmhc)

# Convert the bnlearn object to a gRain object for probabilistic reasoning
bn_grain <- as.grain(bn_fit)

# Setting evidence in the Bayesian Network
# Example: Set "Surface Condition" as "WET" to query the impact on other nodes.
evidence <- list("Surface.Condition" = "WET")

# Apply the evidence to the Bayesian Network
bn_grain_with_evidence <- setEvidence(bn_grain, nodes = names(evidence), states = unlist(evidence))

# Query the conditional probability for "Injury Severity" given the evidence
query_result <- querygrain(bn_grain_with_evidence, nodes = "Injury.Severity")

# Print the result of the conditional probability for "Injury Severity"
print(query_result$`Injury.Severity`)

# Test for conditional independence between "Injury Severity" and "Speed Limit"
# This tests if these two variables are independent, given the data.
independence_test <- ci.test("Injury.Severity", "Speed.Limit", data = df)

# Print the result of the conditional independence test
print(independence_test)

# Find the Markov blanket of "Driver Distracted By"
# The Markov blanket represents the minimal set of nodes that shield this node from the rest of the network.
markov_blanket <- mb(bn_fit, "Driver.Distracted.By")

# Print the Markov blanket for "Driver Distracted By"
print(markov_blanket)

# End of code
