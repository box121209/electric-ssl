# EXPERIMENT 3
# Explore dependence of stream classification performance 
# and of running time on the window graph size.
# Fix dimensional reduction of zip data set, treated as 
# binary 0/not-0 classification problem, and randomly sample 
# data subsets and model hyperparameters, to build data set 
# of performance.

source("functions.R")

library(ElemStatLearn)
library(igraph)
library(Rtsne)

lowdim <- 3

#tsne <- Rtsne(zip.train[,2:257], dims=lowdim, verbose=TRUE)
#write.table(cbind(tsne$Y, zip.train[,1]), 
#            file=sprintf("zip_tsne_dim_%d.txt", lowdim), 
#            col.names=FALSE, row.names=FALSE)

zip <- read.table(sprintf("zip_tsne_dim_%d.txt", lowdim))
names(zip) <- c(sapply(1:lowdim, function(i) sprintf("X%d", i)), "cls")

cls <- sapply(zip$cls, function(i) i==0) 
all.dat <- data.frame(zip[,1:lowdim], cls)


if(lowdim==2){
  plot(all.dat$X1, all.dat$X2, 
       col=all.dat$cls+2, cex=0.2, 
       xlab='', ylab='', frame.plot=0)
}

# EXPERIMENTS
for(expt in 1:1000){

  # the main parameter choices

  datprob <- 0.5*runif(1)    # proportion of unlabelled data
  labprob <- 0.2*runif(1)  # proportion of labeled points
  winsize <- sample(seq(4,32,2), 1)      # size of sliding graph

  # make unlabelled data set

  dat <- all.dat[sample(1:nrow(all.dat), 
                      ceiling(datprob*nrow(all.dat))), ]
  dat <- data.frame(dat, lbl=rep(0, nrow(dat)))
  idx <- 1:nrow(dat)
  dat$lbl[ sample(idx[dat$cls==0], ceiling(labprob*sum(dat$cls==0))) ] <- 1
  dat$lbl[ sample(idx[dat$cls==1], ceiling(labprob*sum(dat$cls==1))) ] <- 1

  # make labels
  instream <- permute_data(dat, random=TRUE, labelsfirst=TRUE)

  if(lowdim==2){
    plot(instream$X1,instream$X2, 
       col=instream$cls+2, 
      pch=19, xlab='', ylab='', 
      frame.plot=0)
    points(instream$X1[instream$lbl==1], instream$X2[instream$lbl==1], lwd=3, cex=1.5)
  }

  # make output data stream

  runtime <- system.time( 
    outstream <- tlp(instream, winsize=winsize)
  )


  # we can now compae the output soft classes with the 
  # (known but withheld) actual stream classes:

  #boxplot(outstream$cls ~ instream$cls, 
  #      xlab="Actual class",
  #      ylab="Estimated class",
  #      col='yellow', frame.plot=0)

  # or scatterplot:

  if(lowdim==2){
    col <- rev(green2red(101))[1+round(outstream$cls,2)*100]
    plot(outstream$X1, outstream$X2, col=col, pch=19, xlab='', ylab='', frame.plot=0)
    points(outstream$X1[outstream$lbl==1], outstream$X2[outstream$lbl==1], lwd=3, cex=1.5)
  }

  # average error:
  err.soft <- mean(outstream$cls[instream$cls==0]) + 1 - mean(outstream$cls[instream$cls==1])
  err.hard <- (length(outstream$cls[instream$cls==0 & outstream$cls >= 0.5])
             + length(outstream$cls[instream$cls==1 & outstream$cls <= 0.5]))/length(outstream$cls)

  #cat(sprintf("Soft error: %g\nHard error: %g\nRun time: %g\n", 
  #          err.soft, 
  #          err.hard,
  #          runtime[3]))
  output <- data.frame(
              n=nrow(dat), 
              cls=sum(dat$cls), 
              lbl=sum(dat$lbl), 
              win=winsize, 
              softerr=err.soft, 
              harderr=err.hard, 
              time=as.numeric(runtime[3]))
  write.table(output, 
              sprintf("zip_expt_%d.txt", lowdim), 
              append=TRUE, row.names=FALSE, col.names=FALSE)
}

# ANALYSIS OF THE RESULTS

lowdim <- 3
res <- read.table(sprintf("zip_expt_%d.txt", lowdim))
names(res) <- c("n", 
                "cls", 
                "lbl", 
                "win", 
                "softerr", 
                "harderr", 
                "time")
attach(res)

# runtime by window size:
normtime <- time/res$n
boxplot(normtime ~ win, 
        xlab="Window size",
        ylab="Time per item",
        col='lightgrey', frame.plot=0)

# prediction performance by window size:
boxplot(harderr ~ win, 
        xlab="Window size",
        ylab="Prediction error",
        col='lightgrey', frame.plot=0)

# prediction performance by proportion of labels:
lblprob <- lbl/res$n
plot(softerr ~ lblprob, 
     col='blue', 
     xlab="Proportion of labels",
     ylab="Error",
     ylim=c(0,0.9),
     frame.plot=0)
points(harderr ~ lblprob, 
     col='red')
legend("topright",
       legend=c("Soft error", "Hard error"),
       col=c("blue","red"),
       pch=1)     

# prediction performance by data size:
plot(softerr ~ res$n, 
     col='blue', 
     xlab="Data size",
     ylab="Error",
     ylim=c(0,1),
     frame.plot=0)
points(harderr ~ res$n, 
       col='red')
legend("topright",
       legend=c("Soft error", "Hard error"),
       col=c("blue","red"),
       pch=1, bty='n')     


