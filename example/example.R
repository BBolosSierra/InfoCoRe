##reset
rm(list=ls())

library(infoCoRe)
library(igraph)
library(ggplot2)
library(scales)
library(Matrix)

set.seed(3)
## load graph
gg = read.graph("/Users/s2595589/WORK/InfoCoRe/example/PPI_Presynaptic.gml",format="gml")

## 
adj=get.adjacency(gg)

dg <- diag(degree(gg))

D = dg - as.matrix(adj)

#decompose(igraph::graph_from_adjacency_matrix(adj, mode = "undirected"))


result <- infoCoRe::driver(Adj=as(adj,"generalMatrix"),
                 weighted=0,
                 directed=0,
                 val_only=1,
                 norm=0,
                 custom_tau_vec=c(seq(0, 50, by = 0.01)))

# for visualization
df <- data.frame(tau = result$tau_vec, entropy = 1 - result$ent_vector, heat=result$specific_heat)



## Plot the full graph
# plot the entropy vs the tau
ggplot(df, aes(x = tau, y = heat)) +
  geom_line() +
  labs(title = "Heat vs. Tau",
       x = "Tau",
       y = "Heat")
# Create a ggplot object for the entropy
p <- ggplot(df, aes(x = tau, y = entropy)) +
  geom_line(aes(y = entropy), linetype = "dashed", color = "blue") +
  labs(title = "Entropy and Heat vs. Tau",
       x = "Tau",
       y = "1 - Entropy")
# Add the heat to the plot with a secondary y-axis
p + geom_line(aes(y = heat), color = "red") +
  geom_vline(xintercept = result$tau_peak_values, linetype = "dashed", color = "black") + 
  #scale_x_continuous(trans = log10_trans(),
  #                   breaks= trans_breaks("log10", function(x) 10^x),
  #                   labels= trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    name = "Heat",
    sec.axis = sec_axis(~., name = "1 - Entropy", breaks = scales::breaks_pretty())
  ) + 
  annotate("rect", xmin = 10, xmax=50, ymin=0, ymax=1, alpha=.1, fill = "blue")

# Find components
components <- decompose(graph_from_adjacency_matrix(result$adj2, mode = "undirected"))
# Find clusters
clus <- clusters(igraph::graph_from_adjacency_matrix(result$adj2, mode = "undirected"))

## to check if results are consistent we can compute it in R 
## it will be quite slow 
# Calculate the operator ρ(τ)
tau <- 8.15
K <- expm(-tau * result$L)
tr <-  sum(diag(K))
rho <- K / tr

n = length(gg)
rho = result$L
adj2 = matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) { # this is the approach in the python code
    adj2[i, j] = ifelse((rho[i, j] >= rho[j, j]) | (rho[i, j] >= rho[i, i]), 1, 0)
  }
}

components <- decompose(graph_from_adjacency_matrix(result$adj2, mode = "undirected"))

# Find connected components
clus <- clusters(igraph::graph_from_adjacency_matrix(result$adj2, mode = "undirected"))

# Contract nodes within each connected component
for (supernode in components$csize) {
  nodes <- sort(components$membership[components$membership == supernode])
  for (node in nodes[-1]) {
    # Contract nodes
    G1 <- contract.vertices(G1, nodes[1], node, self.loops = FALSE)
  }
}

decompose(graph_from_adjacency_matrix(result$adj2))

# Relabel nodes
new_labels <- seq_along(V(G1))
set_vertex_attr(G1, "name", value = as.character(new_labels))

# Print degree of each node
degree_sequence <- degree(G1)
print(degree_sequence)

# calculating degree distributions
degree_distribution <- degree_distribution(gg)

# Assuming 'degree_distribution' is your numeric vector
degree_df <- data.frame(degree = seq_along(degree_distribution), count = degree_distribution)

# Plot the degree distribution
ggplot(degree_df, aes(x = degree, y = count)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Degree Distribution",
       x = "Degree (log scale)",
       y = "Count (log scale)")


## Barabasi–Albert (BA) network

ba <- sample_pa(n=1780,
          power = 2,
          m = 1,
          directed = FALSE,
          algorithm = c("psumtree"))

adj_ba=get.adjacency(ba)


dg_ba <- diag(degree(ba))

D = dg_ba - as.matrix(adj_ba)


result <- infoCoRe::driver(Adj=as(adj_ba,"generalMatrix"),
                           weighted=0,
                           directed=0,
                           val_only=1,
                           norm=0,
                           custom_tau_vec=c(seq(0, 100, by = 1)))

L = as.matrix(result$L)


eigval <- eigen(L)[1]

e <- exp(-eigval$values * result$tau_peak_values[1])
mu <- e/sum(e)

S = sum(mu*log(mu))*(-1/log(dim(L)[1]))

K <- expm(-(result$tau_peak_values[1]) * result$L)
tr <-  sum(diag(K))
rho <- K / tr

mu2 <- eigen(rho)

df <- data.frame(tau = result$tau_vec, entropy = 1 - result$ent_vector, heat=result$specific_heat)


components <- decompose(graph_from_adjacency_matrix(result$adj2))

# Find connected components
clus <- clusters(igraph::graph_from_adjacency_matrix(result$adj2, mode = "undirected"))



# Contract nodes within each connected component
for (supernode in components$csize) {
  nodes <- sort(components$membership[components$membership == supernode])
  for (node in nodes[-1]) {
    # Contract nodes
    G1 <- contract.vertices(G1, nodes[1], node, self.loops = TRUE)
  }
}

# plot the entropy vs the tau
ggplot(df, aes(x = tau, y = entropy)) +
  geom_line() +
  labs(title = "Entropy vs. Tau",
       x = "Tau",
       y = "1 - Entropy")

# plot the entropy vs the tau
ggplot(df, aes(x = tau, y = heat)) +
  geom_line() +
  labs(title = "Heat vs. Tau",
       x = "Tau",
       y = "Heat")

# Create a ggplot object for the entropy
p <- ggplot(df, aes(x = tau, y = entropy)) +
  geom_line(aes(y = entropy), linetype = "dashed", color = "blue") +
  labs(title = "Entropy and Heat vs. Tau",
       x = "Tau",
       y = "1 - Entropy")

# Add the heat to the plot with a secondary y-axis
p + geom_line(aes(y = heat), color = "red") +
  scale_y_continuous(
    name = "Heat",
    sec.axis = sec_axis(~., name = "1 - Entropy", breaks = scales::breaks_pretty())
  )


## Directed Graph Example
## http://www-personal.umich.edu/~mejn/netdata/
##gg = sample_pa(n=100, power=1, m=1,  directed=T)
gg = read.graph("polblogs.gml",format="gml")

adj=get.adjacency(gg)

infoCoRe::driver(Adj=as(adj,"generalMatrix"),
                 weighted=0,
                 directed=1,
                 val_only=1,
                 norm=1)

