# This is an example of how to implement the Newton-Raphson method on complex numbers
# using R's built-in 'complex' number definition, which is written so that you can use +,-,*,/,^ in the usual way

library("viridis")
library("grDevices")

# These are some color palettes you can try. there are many others available.
#palette(viridis(64))
#palette(inferno(64))
palette(hcl.colors(64, palette = "Earth"))
#palette(rainbow(64))


# 1. CREATING A CLUSTER
library(parallel)
cl <- makePSOCKcluster(8)
# First define a function to work with - it will return two values.
# The first will be the value of the function at z (i.e., f(z)), the second will be the derivate of f
# at z (i.e., f'(z))
# F7<-function(z){ # this has three roots 1,-1/2 + sqrt(3)/2i and -1/2 - sqrt(3)/2i
#   # a complex function z^3-1
#   return(c(z^3 - 1, 3 * z^2)) # note that arithmetic ops on complex numbers work just like they should
# }

F7 <- function(x) {
  result <- (x+2)^3 - 2*x
  deriv_value <- deriv(~(x+2)^3 - 2*x, "x") |> eval() |> attr("gradient")
  return(c(result, deriv_value))
}



# Next we define a global variable to control whether we record roots or the number of iterations taken to find them when drawing our picture
# This is an inelegant way of doing it. It should really be part of the argument to the relevent functions.
bRootOrIterations <- 0   # Set <-1 to record which root is found, or <- 0 to record number of iterations needed
clusterExport(cl, "bRootOrIterations")


# Here's the function that performs Newton-Raphson
TwoDNewtonRaphson <- function(func, StartingValue, Tolerance, MaxNumberOfIterations) {
  i <- 0  # something to count the iterations
  NewZ <- StartingValue  # start the algorithm at the complex number 'StartingValue'
  Deviation = abs(func(StartingValue)[1])   # Work out how far away from (0,0) func(NewZ) is.
  
  #Set up a while loop until we hit the required target accuracy or the max. number of steps
  while ((i < MaxNumberOfIterations) && (Deviation > Tolerance)) {
    # Find the next Z-value using Newton-Raphson's formula
    Z <- func(NewZ)   # Remember, this is is a vector of two elements. Z[1] is the is the value of the function; Z[2] is its derivative
    if ((Z[1] == "NaN") || (Z[2] == "NaN")) {
      cat("Function or derivative not defined error.")
      cat("\n", NewZ, Z)
      break
    }
    
    # So we need to calculate the next value of Z using this formula Z(n+1) <- Z(n)-f(Z(n))/f'(z(n))
    NewZ <- NewZ - Z[1] / Z[2]
    
    # calculate how far f(z) is from 0
    NewVal <- func(NewZ)
    Deviation <- abs(NewVal[1])
    i <- i + 1
    #cat(paste("\nIteration ",i,":   Z=",NewZ,"  Devn=",Deviation))
  }
  
  # output the result
  if (Deviation > Tolerance) {
    cat(paste("\nConvergence failure. Deviation:", Deviation, "after ", i, 	"iterations"))
  }
  
  # what the function returns depends upon whether you are counting how many iterations it takes
  # to converge or checking which root it converged to...
  if (bRootOrIterations == 1) {
    return(NewZ)
  } else {
    return(c(i, i))
  }
}

# A function to check whether two points are close together
CloseTo <- function(x, y) {
  # returns 1 if x is close to y
  if (abs(x - y) < 0.1) {
    return(1)
  } else {
    return(0)
  }
}

# And now here's the function that will draw a pretty picture
Root_calculator <- function(Funcn, xmin, xmax, xsteps, ymin, ymax, ysteps) {
  # First define a grid of x and y coordinates over which to run Newton-Raphson.
  # When we run ut for the point (x,y) it will start with the complex number x+iy
  
  x <- seq(xmin, xmax, length.out = xsteps)
  y <- seq(ymin, ymax, length.out = ysteps)
  
  out_dat <- expand.grid(x = x, y = y)
  
  ThisZ <- complex(1, out_dat$x, out_dat$y)
  
  Root <- parSapply(cl = cl,
                    X = ThisZ,
                    FUN = TwoDNewtonRaphson,
                    func = Funcn,
                    Tolerance = 1e-1,
                    MaxNumberOfIterations = 100)
  
  if(bRootOrIterations == 0) {
    out_dat$color <- 261 + 5 * Root[1, ]
    out_dat$root1 <- Root[1, ]
    out_dat$root2 <- Root[2, ]
  } else {
    out_dat$color <- 261 + 5 * Root
    out_dat$root1 <- Root
  }
  return(out_dat)
}

root_plotter <- function(x, xmin, xmax, ymin, ymax, PtSize = 0.02) {
  # And now we have everything, so let's draw the picture.
  # plot(x[, 1:2], col = x[, 3],
  plot(x[, c(2,1)], col = x[, 3],
       xlim = c(xmin, xmax), ylim = c(ymin, ymax),
       pch = 16, cex = PtSize * 3,
       yaxt = "n",
       xaxt = "n",
       box = FALSE)
  # The graphics parameter pch controls what shape is used to plot each point
  # The graphics parameter cex is used to control the size of each point
  # Change these to improve the quality of your pictures
}

A <- Root_calculator(F7, -0.4, 0.3, 1440, -0.5, 0.5, 2560)

stopCluster(cl)
root_plotter(A,-0.5,0.5,-0.4,0.3,1)

png(file= "C:\\Users\\jwill\\Desktop\\image_assignment.png", width = 2560, height = 1440)


































