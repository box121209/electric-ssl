# EXPERIMENT 2  
# Explore TLP implementation on synthetic data,
# for performance and visualisations.

source("functions.R")


# the main parameter choices

labprob1 <- 0.1 # proportion of labeled points
labprob2 <- 0.05 # proportion of labeled points
ssize1 <- 100    # sample size per class
ssize2 <- 50    # sample size per class
mixing <- 10.0   # relative sd of the two classes
winsize <- 6    # size of sliding graph

dat <- make_2_normal_clusters(ssize1,ssize2, 
                              l1=ceiling(labprob1*ssize1), 
                              l2=ceiling(labprob2*ssize2), 
                              sd1=mixing)

# input/output data streams...
# either astream through one class and then the other, 
# or permute randomly:
instream <- permute_data(dat, random=TRUE, labelsfirst=TRUE)

plot(instream$X1,instream$X2, col=instream$cls+2, pch=19, cex=2,
     xlab='', ylab='', axes=FALSE)
points(instream$X1[instream$lbl==1], instream$X2[instream$lbl==1], lwd=3, cex=2.5)

runtime <- system.time( 
  outstream <- tlp(instream, winsize=winsize, viz=TRUE)
)

# we can now compare the output soft classes with the 
# (known but withheld) actual stream classes:

boxplot(outstream$cls ~ instream$cls, 
        xlab="Actual class",
        ylab="Voltage",
        col='grey', frame.plot=0)

# or scatterplot:

col <- rev(green2red(101))[1+round(outstream$cls,2)*100]
plot(outstream$X1, outstream$X2, col=col, pch=19, cex=2, xlab='', ylab='', 
     axes=FALSE)
points(outstream$X1[outstream$lbl==1], outstream$X2[outstream$lbl==1], lwd=3, cex=2.5)

# average error:
err.soft <- mean(outstream$cls[instream$cls==0]) + 1 - mean(outstream$cls[instream$cls==1])
err.hard <- (length(outstream$cls[instream$cls==0 & outstream$cls >= 0.5])
+ length(outstream$cls[instream$cls==1 & outstream$cls <= 0.5]))/length(outstream$cls)

cat(sprintf("Soft error: %g\nHard error: %g\nRun time: %g\n", 
            err.soft, 
            err.hard,
            runtime[3]))


