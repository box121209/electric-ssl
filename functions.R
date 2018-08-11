library(igraph)

make_2_normal_clusters <- function(n1, n2, 
                                   l1=1, 
                                   l2=1, 
                                   mu1=c(2,1),
                                   mu2=c(-2,-1),
                                   sd1=1.0,
                                   sd2=1.0){
  
  x1 <- rnorm(n1, mean=mu1[1], sd=sd1)
  y1 <- rnorm(n1, mean=mu1[2], sd=sd1)
  x2 <- rnorm(n2, mean=mu2[1], sd=sd2)
  y2 <- rnorm(n2, mean=mu2[2], sd=sd2)
  
  dat1 <- cbind(x1, y1, rep(0,n1), c(rep(1,l1), rep(0,n1-l1)))
  dat2 <- cbind(x2, y2, rep(1,n2), c(rep(1,l2), rep(0,n2-l2)))
  dat <- data.frame(rbind(dat1, dat2))
  names(dat) <- c('X1', 'X2', 'cls', 'lbl')
  
  # return:
  dat
}


make_graph <- function(dat, init_threshold=1.0){
  # dat should be data frame of Euclidean vectors
  dmat <- as.matrix(dist(dat, method='euclidean'))
  n <- nrow(dmat)
  
  # set nearest neighbour threshold just large enough 
  # to connect the graph:
  threshold <- init_threshold
  cc <- 2
  while(cc > 1){
    adj <- (dmat <= threshold)
    for(i in 1:n) adj[i,i] <- 0
    g <- graph.adjacency(adj, mode="undirected")
    cc <- clusters(g)$no
    threshold <- 1.01 * threshold
  }
  # return
  list(g=g, threshold=threshold)
}

permute_data <- function(dat, random=TRUE, labelsfirst=TRUE){
  # takes data frame with 'cls' and 'lbl' columns
  # and outstreams 'instream' frame beginning with two 
  # labelled rows and permuting the remainder randomly
  idx <- 1:nrow(dat)
  i1 <- which(dat$lbl==1 & dat$cls==0)[1]
  i2 <- which(dat$lbl==1 & dat$cls==1)[1]
  if(random==TRUE){
    # start with label from each class:
    acc <- c(i1,i2)
    # then the rest of the labels:
    if(labelsfirst){
      acc <- c(acc, 
               which(dat$lbl[idx]==1 & !(idx %in% acc)))
    }
    # randomly permute the rest:
    acc <- c(acc,
             sample(which(!(idx %in% acc))))
    instream <- dat[acc,]
    
  } else {
      ord <- order(1-dat$lbl, dat$cls)
      ord <- c(i1, i2, ord[!(idx[ord] %in% c(i1,i2))])
      instream <- dat[ord,]
  }
  # return:
  instream
}

# euclidean similarity function:
simfn <- function(p1, p2){ 1/sqrt(sum((p1 - p2)^2)) }

tlp_1 <- function(instream, winsize=6){
  # takes input instream (data frame) with occasional 
  # class binary labels and applies 'temporal lable propagation'
  # i.e. electric model on a sliding window NN graph,
  # to outstream a instream with soft class estimates added to all points
  
  n <- ncol(instream) - 2 # vector space dimension
  h <- graph.empty(0, directed=FALSE)

  # iteratively build outstream data frame
  outstream <- data.frame(instream[1:2,])
  
  h <- h %>% add_vertices(1, attr=cbind(instream[1,],timestamp=1))
  h <- h %>% add_vertices(1, attr=cbind(instream[2,],timestamp=2))
  
  for(thistime in 3:nrow(instream)){

    # read from instream:
    item <- instream[thistime,]
    
    if(item$lbl==1){ # if receive labelled point...
      
      b <- 1 + item$cls
      
      # update weights to the non-poles:
      for(v in V(h)){
        
        if(V(h)$timestamp[v] %in% c(1,2)) next
        
        pt <- outstream[V(h)$timestamp[v], 1:n]
        wt <- simfn(pt, item[1:n])

        if(are_adjacent(h,b,v)){
          e <- get.edge.ids(h, c(b,v))
          E(h)$weight[e] <- E(h)$weight[e] + wt
        } else {
          h <- add.edges(h, c(b,v), attr = list(weight=wt))
        }
      }
      
      # write to outstream:
      outstream[thistime,] <- item
      
    } else { # if receive unlabelled point...
      
      h <- add_vertices(h, 1, attr=cbind(instream[thistime,],
                                         timestamp=thistime))
      thisv <- which(V(h)$timestamp==thistime)
      
      # update weights to the poles:
      vert.1 <- outstream[outstream$lbl==1 & outstream$cls==0, 1:n]
      wt.1 <- sum(
        sapply(1:nrow(vert.1), function(j) 
          simfn(vert.1[j,], item[1,1:n])
        )
      )

      h <- add.edges(h, c(1, thisv), 
                     attr = list(weight=wt.1))
      vert.2 <- outstream[outstream$lbl==1 & outstream$cls==1, 1:n]
      wt.2 <- sum(
        sapply(1:nrow(vert.2), function(j) 
          simfn(vert.2[j,], item[1:n]))
      )
      h <- add.edges(h, c(2, thisv), 
                     attr = list(weight=wt.2))
      
      # update weights to the non-poles:
      for(v in V(h)){
        if(V(h)$timestamp[v] %in% c(1,2,thistime)) next
        pt <- sapply(1:n, function(i) 
          eval(parse(text = sprintf("V(h)[%d]$X%d", v, i))))
        wt <- simfn(pt, item[1:n])
        h <- add.edges(h, c(v, thisv), 
                       attr = list(weight=wt))
      }
      
      # star-mesh transform:
      if(length(V(h)) > winsize+2){
        
        nonpole <- V(h)$timestamp[V(h)$timestamp>2]
        v.oldest <- which(V(h)$timestamp==min(nonpole))
        deg.v.oldest <- degree(h, v.oldest)
        if(deg.v.oldest==0){ 
          print("Dividing by zero!")
          stop()
        }
        for(p in V(h)){
          if(!are_adjacent(h,p,v.oldest)) next
          for(q in V(h)){
            if(p==q | p==v.oldest | q==v.oldest) next
            if(!are_adjacent(h,q,v.oldest)) next
            e.old.1 <- get.edge.ids(h, c(p, v.oldest))
            e.old.2 <- get.edge.ids(h, c(q, v.oldest))
            delta <- E(h)$weight[e.old.1]*E(h)$weight[e.old.2] / deg.v.oldest
            if(delta==Inf | delta==-Inf){
              cat("Adding infinity!\n")
              stop()
            }
            if(are_adjacent(h,p,q)){
              e <- get.edge.ids(h, c(p,q))
              E(h)$weight[e] <- E(h)$weight[e] + delta
            } else {
              h <- add.edges(h, c(p,q), attr = list(weight=delta))
            }
          }
        }
        h <- delete_vertices(h,v.oldest)
        thisv <- which(V(h)$timestamp==thistime)
      }
      
      # harmonic solution applied to h:
      A <- as.matrix(as_adjacency_matrix(h, attr="weight"))
      laplace <- diag(rowSums(A)) - A
      
      lap.uu <- laplace[V(h)$lbl==0, V(h)$lbl==0]
      lap.ul <- laplace[V(h)$lbl==0, V(h)$lbl==1]
      
      lab.l <- V(h)$cls[V(h)$lbl==1]
      lab.u <- -solve(lap.uu) %*% lap.ul %*% lab.l
      
      if(thisv > winsize+2){
        cat("Inconsistent vertex number!\n")
        stop()
      }
      soft.cls <- lab.u[thisv - 2]
      
      # write to outstream instream:
      outstream[thistime,] <- item
      outstream[thistime,n+1] <- soft.cls
    }
    
    # renormalise the weights:
    E(h)$weight <- E(h)$weight/sum(E(h)$weight)
    countdown <- nrow(instream)-thistime
    cat(countdown, '\r')
  }
  # return
  outstream
}

tlp_2 <- function(instream, winsize=6, viz=FALSE){
  # takes input instream (data frame) with occasional 
  # class binary labels and applies 'temporal lable propagation'
  # i.e. electric model on a sliding window NN graph,
  # to outstream a instream with soft class estimates added to all points
  
  n <- ncol(instream) - 2 # vector space dimension
  h <- graph.empty(0, directed=FALSE)
  dmat <- 1/as.matrix(dist(instream[,1:n], method='euclidean'))
  
  # iteratively build outstream data frame
  outstream <- data.frame(instream[1:2,])
  
  h <- h %>% add_vertices(1, attr=list(timestamp=1, 
                                       lbl=instream$lbl[1],
                                       cls=instream$cls[1]))
  h <- h %>% add_vertices(1, attr=list(timestamp=2, 
                                       lbl=instream$lbl[2],
                                       cls=instream$cls[2]))
  
  for(thistime in 3:nrow(instream)){
    
    # read from instream:
    item <- instream[thistime,]
    
    if(item$lbl==1){ # if receive labelled point...
      
      b <- 1 + item$cls
      
      # update weights to the non-poles:
      for(v in V(h)){
        
        if(V(h)$timestamp[v] %in% c(1,2)) next
        
        wt <- dmat[thistime, V(h)$timestamp[v]]
        if(are_adjacent(h,b,v)){
          e <- get.edge.ids(h, c(b,v))
          E(h)$weight[e] <- E(h)$weight[e] + wt
        } else {
          h <- add.edges(h, c(b,v), attr = list(weight=wt))
        }
      }
      
      # write to outstream instream:
      outstream[thistime,] <- item
      
    } else { # if receive unlabelled point...
      
      h <- add_vertices(h, 1, attr=list(timestamp=thistime, 
                                        lbl=instream$lbl[thistime],
                                        cls=instream$cls[thistime]))
      thisv <- which(V(h)$timestamp==thistime)
      
      # update weights to the poles:
      vert.1 <- which(outstream$lbl==1 & outstream$cls==0)
      wt.1 <- sum( sapply(vert.1, function(j) dmat[j, thistime]) )
      h <- add.edges(h, c(1, thisv), 
                     attr = list(weight=wt.1))
      vert.2 <- which(outstream$lbl==1 & outstream$cls==1)
      wt.2 <- sum( sapply(vert.2, function(j) dmat[j, thistime]) )
      h <- add.edges(h, c(2, thisv), 
                     attr = list(weight=wt.2))
      
      # update weights to the non-poles:
      for(v in V(h)){
        if(V(h)$timestamp[v] %in% c(1,2,thistime)) next
        pt <- outstream[V(h)$timestamp[v], 1:n]
        wt <- simfn(pt, item[1:n])
        h <- add.edges(h, c(v, thisv), 
                       attr = list(weight=wt))
      }
      
      # star-mesh transform:
      if(length(V(h)) > winsize+2){
        nonpole <- V(h)$timestamp[V(h)$timestamp>2]
        v.oldest <- which(V(h)$timestamp==min(nonpole))
        deg.v.oldest <- degree(h, v.oldest)
        if(deg.v.oldest==0){ 
          stop()
          print("Dividing by zero!")
        }
        for(p in V(h)){
          if(!are_adjacent(h,p,v.oldest)) next
          for(q in V(h)){
            if(p==q | p==v.oldest | q==v.oldest) next
            if(!are_adjacent(h,q,v.oldest)) next
            e.old.1 <- get.edge.ids(h, c(p, v.oldest))
            e.old.2 <- get.edge.ids(h, c(q, v.oldest))
            delta <- E(h)$weight[e.old.1]*E(h)$weight[e.old.2] / deg.v.oldest
            if(delta==Inf | delta==-Inf){
              cat("Adding infinity!\n")
              stop()
            }
            if(are_adjacent(h,p,q)){
              e <- get.edge.ids(h, c(p,q))
              E(h)$weight[e] <- E(h)$weight[e] + delta
            } else {
              h <- add.edges(h, c(p,q), attr = list(weight=delta))
            }
          }
        }
        h <- delete_vertices(h,v.oldest)
        thisv <- which(V(h)$timestamp==thistime)
      }
      
      # harmonic solution applied to h:
      A <- as.matrix(as_adjacency_matrix(h, attr="weight"))
      laplace <- diag(rowSums(A)) - A
      
      lap.uu <- laplace[V(h)$lbl==0, V(h)$lbl==0]
      lap.ul <- laplace[V(h)$lbl==0, V(h)$lbl==1]
      
      lab.l <- V(h)$cls[V(h)$lbl==1]
      lab.u <- -solve(lap.uu) %*% lap.ul %*% lab.l
      
      if(thisv > winsize+2){
        cat("Inconsistent vertex number!\n")
        stop()
      }
      soft.cls <- lab.u[thisv - 2]
      
      # write to outstream instream:
      outstream[thistime,] <- item
      outstream[thistime,n+1] <- soft.cls
    }
    
    # renormalise the weights:
    E(h)$weight <- E(h)$weight/sum(E(h)$weight)
    countdown <- nrow(instream)-thistime
    cat(countdown, '\r')
    if(viz){
      col <- rev(green2red(101))[1+round(outstream$cls,2)*100]
      out <- outstream[1:thistime,]
      
      setwd("../tmp")
      png(sprintf("img_%d.png", thistime))
      plot(instream$X1,instream$X2, col='grey', pch=1, cex=2,
           xlab='', ylab='', axes=FALSE)
      for(u in V(h)){
        tu <- V(h)$timestamp[u]
        for(v in V(h)){
          if(u==v) next
          tv <- V(h)$timestamp[v]
          segments(instream[tu,1],
                   instream[tu,2],
                   instream[tv,1],
                   instream[tv,2],
                   col='grey')
        }
      }
      points(out$X1,out$X2,
             col=col[1:thistime], pch=19, cex=2)
      points(out$X1[out$lbl==1], 
             out$X2[out$lbl==1], 
             lwd=3, cex=2.5)
      dev.off()
      setwd("../electric-ssl")
    }
  }
  # return
  outstream
}

tlp <- tlp_2
