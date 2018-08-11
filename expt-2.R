source("functions.R")


# the main parameter choices

labprob <- 0.03 # proportion of labeled points
winsize <- 6    # size of sliding graph
ssize <- 100    # sample size per class
mixing <- 1.0   # relative sd of the two classes

dat <- make_2_normal_clusters(ssize,ssize, 
                              l1=ceiling(labprob*ssize), 
                              l2=ceiling(labprob*ssize), 
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


