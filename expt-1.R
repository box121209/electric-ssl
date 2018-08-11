# EXPERIMENT 1
# Explore the harmonic SSL solution on synthetic data,
# including NN graph construction with distance threshold
# just enough to form connected graph.

source("functions.R")

# generate some data
dat <- make_2_normal_clusters(100,100, l1=2, l2=2, sd1=1.5)
attach(dat)

plot(X1,X2, col=cls+2, pch=19, xlab='', ylab='', axes=FALSE)
points(X1[lbl==1], X2[lbl==1], lwd=3, cex=1.5)

# build NN graph

library(igraph)

dmat <- as.matrix(dist(dat[,1:2], method='euclidean'))
n <- nrow(dmat)


threshold <- 1.0
cc <- 2
while(cc > 1){
  adj <- (dmat <= threshold)
  for(i in 1:n) adj[i,i] <- 0
  g <- graph.adjacency(adj, mode="undirected")
  cc <- clusters(g)$no
  threshold <- 1.01 * threshold
}

lout <- as.matrix(dat[,1:2])

#V(g)$color <- sapply(lbl[V(g)] + cls[V(g)], 
#                     function(i) 
#                       if(i==0) 'white' else if(i==1) 'red' else if(i==2) 'green'
#                     )
#plot(g, 
#     vertex.size=3, 
#     vertex.label=NA, 
#     layout=lout,
#     vertex.frame.color="grey", 
#     edge.color="lightgrey"
#)

# the harmonic solution

laplace <- diag(rowSums(adj)) - adj

A.uu <- laplace[lbl==0, lbl==0]
A.ul <- laplace[lbl==0, lbl==1]

lab.l <- cls[lbl==1]
lab.u <- -solve(A.uu) %*% A.ul %*% lab.l

voltage <- rep(0,n)
voltage[lbl==1] <- lab.l
voltage[lbl==0] <- lab.u

library(colorRamps)

V(g)$color <- rev(green2red(101))[1+round(voltage,2)*100]
V(g)$frame.color <- sapply(lbl, function(i) if(i==0) 'grey' else 'black')
V(g)$size <- sapply(lbl, function(i) if(i==0) 3 else 5)

plot(g, 
     vertex.label=NA, 
     layout=lout,
     edge.color="lightgrey"
)
hist(voltage, breaks=20, col='lightgrey', main="")

V(g)$color <- sapply(cls, function(i) if(i==0) 'red' else 'green')
V(g)$frame.color <- sapply(lbl, function(i) if(i==0) 'grey' else 'black')
V(g)$size <- sapply(lbl, function(i) if(i==0) 3 else 5)

plot(g, 
     vertex.label=NA, 
     layout=lout,
     edge.color="white"
)

