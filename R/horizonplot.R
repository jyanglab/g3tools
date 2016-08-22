
#http://rgraphgallery.blogspot.com/2013_05_01_archive.html


require(latticeExtra)

#example data
set.seed(123)
mydat <- ts(matrix(cumsum(rnorm(150 * 10)), ncol = 10))
colnames(mydat) <- paste("TS", letters[1:10], sep = "-")

#simple line plot
xyplot(mydat, scales = list(y = "same"))


# panel with different origin and scale:
horizonplot(mydat, layout = c(1,12), colorkey = TRUE) +  
    layer(panel.scaleArrow(x = 0.99, digits = 1, col = "grey",
                           srt = 90, cex = 0.7)) +
    layer(lim <- current.panel.limits(),
          panel.text(lim$x[1], lim$y[1], round(lim$y[1],1), font = 2,
                     cex = 0.7, adj = c(-0.5,-0.5), col = "blue")) 
