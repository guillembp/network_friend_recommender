library(igraph)
library(dplyr)
library(reshape2)
getwd()
setwd('~')
edge_list <- read.csv2('Desktop/networks_data/email-Eu-core.txt', sep = ' ')
communities <- read.csv2('Desktop/networks_data/email-Eu-core-department-labels.txt', sep = ' ')
colnames(communities) <- c('Node','Community')
communities <- arrange(communities,communities$Node) # Order by node number Node

# Remove 0's from the node list
# edge_list <- edge_list-1
# communities <- communities

#Rank nodes by number of connections to identify hubs
node_connections <- as.data.frame(table(unlist(edge_list)))
arrange(desc(node_connections$Freq)) 
colnames(node_connections) <- c('Node','EdgeQuant')
node_connect_comm <- merge.data.frame(x=node_connections, y=communities, by = 'Node')

# Create a manageable dataset
Nodes <- 300
MinEdgesPerNode <- 10
MaxEdgesPerNode <- 600

# Delete nodes with less than 10 nodes and more than MaxEdgesPerNode 
node_connections_useful <- node_connect_comm[node_connect_comm$EdgeQuant>=MinEdgesPerNode & node_connect_comm$EdgeQuant<=MaxEdgesPerNode,]

# Create a new edge list only containing nodes in previous list
edgelist_new <- edge_list[edge_list$X0 %in% node_connections_useful$Node & edge_list$X1 %in% node_connections_useful$Node,]

# Create a strongly connected graph using breadth first search
G <- graph_from_data_frame(d = edgelist_new, directed = FALSE)
f <- function(graph, data, extra) {
 data['rank'] == Nodes
}
bfsy  <- bfs(G, root=120, "all", order=TRUE, callback = f)
suby <- bfsy$order[1:Nodes]
newgraph = induced.subgraph(G, suby)

# Generate adjacency matrix from the edge list
Adj <- as_adjacency_matrix(graph = newgraph, type = 'both', sparse=FALSE, names = TRUE)
Adj <- Adj[1:Nodes,1:Nodes]
Adj <- Adj*upper.tri(Adj, diag = FALSE)

newgraph <- simplify(graph=newgraph, remove.multiple = TRUE, remove.loops = TRUE)
e <- get.edgelist(newgraph, names = FALSE) # Rename edges from 1 to #Nodes
colnames(e) = c('X0','X1')

# Plot graph
pg <- graph_from_data_frame(d = e, directed = FALSE)
plot(pg)

# Ravasz algorithm for community clustering, dendrogram plotting
X <- matrix(0, nrow = Nodes, ncol = Nodes)
for(i in 1:Nodes){
  print(paste0("Row: ",i))
  for(j in 1:Nodes){
    if(i!=j){
      print(j)
      di <- 0
      dj <- 0
      J <- 0
      for(k in 1:Nodes){
        print(k)
        J <- J+Adj[i,k]*Adj[j,k]
        di <- di+Adj[i,k] 
        dj <- dj+Adj[j,k]
      }
      if(J>0){
        X[i,j] <- J/(min(di,dj))
        print('positive')
      } else {
        X[i,j] <- J/(min(di,dj)+1)
        print('negative')
      }
    }
  }
}

# Plot dendrogram
# prepare hierarchical cluster
hc = hclust(dist(X))
plot(hc)

# Remove edges at random
number_of_edges = length(e)/2
dropped <- round(runif(n=integer(0.1*number_of_edges), min = 1, max = number_of_edges))
new_edge_truth <- data.frame(e[dropped,])
e_missing <- data.frame(e[-dropped,])
Ge_m <- graph_from_data_frame(d = e_missing, directed = FALSE)
Ae_m <- as_adjacency_matrix(graph = Ge_m, type = 'both', sparse=FALSE, names = FALSE)
# colnames(Ae_m) <- c(1:Nodes)

# Predict the next edges that will be created
# Around 5%, very low performance
# 1. Cosine Similarity.                             
CosSim <- function(ix){
    Q = Ae_m[ix[1],]
    W = Ae_m[ix[2],]
    if(Q==0||W==0){
      return(0)
    }
    return(sum(Q*W)/sqrt(sum(Q^2)*sum(W^2)))
}   
cmb <- expand.grid(i=1:Nodes, j=1:Nodes)
similarity_matrix <- matrix(data = apply(X = cmb, MARGIN = 1, FUN = CosSim), nrow = Nodes, ncol = Nodes)

# 3. Jaccard. Up to 30% performance!
similarity_matrix <- similarity(Ge_m, method = "jaccard")

# 4. Iverse log-weighted. 
similarity_matrix <- similarity(Ge_m, method = "invlogweighted")


# Fix matrix
diag(similarity_matrix) <- 0
similarity_matrix[lower.tri(similarity_matrix)] <- 0


# Recover a list of edges with predicted high cosine similarity
candidate_edges <- melt(similarity_matrix, varnames = c('row', 'col'))
candidate_edges <- arrange(candidate_edges, desc(candidate_edges$value))

# Rating the performance
# Take the top link suggestions over number in truth
links <- length(new_edge_truth$X0)
top <- candidate_edges[1:links,]
r <- data.frame(X0=integer(), X1=integer())
for(i in 1:links){
  for(j in 1:links){
    if(top$row[i]==new_edge_truth$X0[j]){
      if(top$col[i]==new_edge_truth$X1[j]){
        r[nrow(r) + 1,] = list(top$row[i],top$col[i])
        break
      }
    }
    if(top$col[i]==new_edge_truth$X0[j]){
      if(top$row[i]==new_edge_truth$X1[j]){
        r[nrow(r) + 1,] = list(top$row[i],top$col[i])
        break
      }
    }
  }
}
score <- nrow(r)/links



