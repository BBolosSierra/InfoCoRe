##reset
rm(list=ls())

library(infoCoRe)
library(igraph)
library(ggplot2)

## load graph
gg = read.graph("/Users/s2595589/WORK/InfoCoRe/example/PPI_Presynaptic.gml",format="gml")

## 
adj=get.adjacency(gg)


result <- infoCoRe::driver(Adj=as(adj,"generalMatrix"),
                 weighted=0,
                 directed=0,
                 val_only=1,
                 norm=1,
                 custom_tau_vec=c(seq(1, 100, by = 0.05), seq(110, 210, by = 50)))

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
  geom_vline(xintercept = x_value_max_heat, linetype = "dashed", color = "black") + 
  scale_x_continuous(trans = log10_trans(),
                     breaks= trans_breaks("log10", function(x) 10^x),
                     labels= trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    name = "Heat",
    sec.axis = sec_axis(~., name = "1 - Entropy", breaks = scales::breaks_pretty())
  ) + 
  annotate("rect", xmin = 10, xmax=50, ymin=0, ymax=1, alpha=.1, fill = "blue")

## similar to HM-cp from https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.4.033196

degree_distribution <- degree_distribution(gg)

dg <- degree(gg)

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


#dg <- degree(ba)
#degree_distribution <- degree_distribution(ba)

result <- infoCoRe::driver(Adj=as(adj_ba,"generalMatrix"),
                           weighted=0,
                           directed=0,
                           val_only=1,
                           norm=1,
                           custom_tau_vec=c(seq(1, 100, by = 10), seq(110, 210, by = 50)))

df <- data.frame(tau = result$tau_vec, entropy = 1 - result$ent_vector, heat=result$specific_heat)

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

