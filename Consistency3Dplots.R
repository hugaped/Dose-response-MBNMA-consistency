library(plot3D)
library(rgl)

####### INPUTS ######

d.1 <- c(0, 5, 9, 10, 14) / 10
d.2 <- c(0, 0.2, 0.25, 0.4, 0.7) / 5


data <- data.frame("studyID"=c(1,1,1,1,1,
                               2,2,2,
                               3,3,
                               4,4,
                               5,5,5,
                               6,6,
                               7,7,7,7
),
"treat"=c(1,3,3,3,4,
          1,2,2,
          2,2,
          3,4,
          1,4,5,
          1,5,
          3,3,3,3
),
"dose"=c(0,2,1.7,1,1,
         0,2,1,
         0,1.7,
         0.5,0.5,
         0,1,1,
         0,1,
         0.5,1,1.7,2
))


#data <- do.call("rbind",
#                by(data, INDICES=list(data$studyID),
#                   FUN=function(DF) if (DF$treat %in% 1) {} else {DF}))




####### ARGUEMENTS ######

res <- 10 # set resolution (smoothness of dose-response curves) - 10=default
studies <- c(1,2) # studies to inclulde (also decides which agents to plot therefore)
size <- 3
contrasts.show <- unique(data$studyID)
points.show <- unique(data$studyID)

emax.fun <- function(d.1,d.2,dose) {
  0 + ((d.1 * dose) / (exp(d.2) + dose))
}




####### FUNCTION ######
agent <- unique(data$treat)
curves <- data.frame("y"=0, "x"=0, "z"=0, "agent"=unique(data$treat))
axislim <- (max(data$dose) / res) + max(data$dose)
for (i in 2:length(agent)) {
  for (k in 1:res) {
    row <- c( emax.fun(d.1[i], d.2[i], k * (axislim/res)),
              k * (axislim/res),
              k * (axislim/res) * agent[i],
              agent[i])
    curves <- rbind(curves, row)
  }
}

# Plot curves
cols <- c(NA, "red", "green", "blue", "black")

i <- 2
subset <- subset(curves, agent==agent[i])
plot3d(x=subset$x, y=subset$y, z=subset$z, col=cols[agent[i]], type="l", lwd=3, surf=NULL,
       box=FALSE, axes=TRUE, xlab="", ylab="Response", zlab="Dose",
       plot=TRUE)
for (i in 3:length(agent)) {
  subset <- subset(curves, agent==agent[i])
  plot3d(x=subset$x, y=subset$y, z=subset$z, col=cols[agent[i]], type="l", lwd=3, add=TRUE)
  #plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE)
}


# Add points

if (!is.null(points.show)) {
  for (study in seq_along(points.show)) {
    points <- unique(subset(data, studyID==points.show[study])[,-1])
    y <- vector()
    z <- vector()

    for (i in seq_along(points$treat)) {
      y <- append(y, emax.fun(d.1[points$treat[i]], d.2[points$treat[i]], points$dose[i]))
      z <- append(z, points$dose[i] * points$treat[i])
    }

    points <- cbind(points, y, z)

    plot3d(x=points$dose, points$y, points$z, type="s", size=size, surf=NULL,
           plot=TRUE, add=TRUE)
  }
}


# Add contrasts (could be for a specific study given a name or number, or for all (as below))
if (!is.null(contrasts.show)) {
  for (study in seq_along(contrasts.show)) {
    contrasts <- data.frame("x"=NULL, "y"=NULL, "z"=NULL)
    subset <- subset(data, studyID==contrasts.show[study])
    for (i in seq_along(subset$studyID)) {
      k <- i + 1
      while (k<=nrow(subset)) {
        rows <- rbind(
          c(subset$dose[i],
            emax.fun(d.1[subset$treat[i]], d.2[subset$treat[i]], subset$dose[i]),
            subset$dose[i] * subset$treat[i]
          ),
          c(subset$dose[k],
            emax.fun(d.1[subset$treat[k]], d.2[subset$treat[k]], subset$dose[k]),
            subset$dose[k] * subset$treat[k]
          )
        )
        contrasts <- rbind(contrasts, rows)
        k <- k+1
      }
    }
    names(contrasts) <- c("x", "y", "z")

    plot3d(contrasts$x, contrasts$y, contrasts$z, type="l", lwd=1, surf=NULL,
           plot=TRUE, add=TRUE)
  }
}







######## ABC trial ###########


#### Dose-response curves

A <- data.frame("x"=0, "y"=0, "z"=0, "agent"=1)
B <- data.frame("x"=0, "y"=0, "z"=0, "agent"=2)
C <- data.frame("x"=0, "y"=0, "z"=0, "agent"=3)
row <- vector()

for (i in 1:10) {
  row <- c(i, exp(0.5*i), 1*i, 1)
  A <- rbind(A, row)

  row <- c(i, exp(0.4*i), 2*i, 2)
  B <- rbind(B, row)

  row <- c(i, exp(0.45*i), 1.5*i, 3)
  C <- rbind(C, row)
}

data <- rbind(A,B,C)

# All lines plotted at once
#plot3d(data$x, data$y, data$z, col=as.numeric(data$agent), type="l", surf=NULL,
#               box=FALSE, axes=FALSE, xlab="Time", ylab="Response", zlab="",
#                theta=100, phi=2,
#               plot=TRUE)

plot3d(A$x, A$y, A$z, col="blue", type="l", surf=NULL,
       box=FALSE, axes=TRUE, xlab="", ylab="Response", zlab="Dose",
       plot=TRUE)

plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE)
plot3d(C$x, C$y, C$z, col="green", type="l", add=TRUE, plot=TRUE)


#### Points

# Black
arm1 <- data.frame("x"=8,
                   "y"=exp(0.5*8),
                   "z"=1*8)

arm5 <- data.frame("x"=6,
                   "y"=exp(0.5*6),
                   "z"=1*6)

# Red
arm2 <- data.frame("x"=8,
                   "y"=exp(0.4*8),
                   "z"=2*8)

# Green
arm3 <- data.frame("x"=9,
                   "y"=exp(0.45*9),
                   "z"=1.5*9)

# Placebo
arm4 <- data.frame("x"=0,
                   "y"=0,
                   "z"=0)

#points <- rbind(arm1, arm2, arm3, arm5)
#points <- rbind(arm1, arm2, arm3)
points <- rbind(arm1, arm2, arm3, arm4)

plot3d(points$x, points$y, points$z, type="s", size=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)


#### Texts

texts <- points
texts$y <- texts$y+15
#texts$value <- c(expression(paste("f(x"["2"], ",", beta["A"], ")")),
#                 expression(paste("f(x"["1"], ",", beta["B"], ")")),
#                 expression(paste("f(x"["1"], ",", beta["C"], ")"))
#)
#texts$value <- c(expression(paste("f(x"["2"], ",", beta["A"], ")")),
#                 expression(paste("f(x"["1"], ",", beta["B"], ")")),
#                 expression(paste("f(x"["1"], ",", beta["C"], ")")),
#                 expression(paste("f(x"["1"], ",", beta["A"], ")"))
#)
texts$value <- c(expression(paste("f(x"["2"], ",", beta["A"], ")")),
                 expression(paste("f(x"["1"], ",", beta["B"], ")")),
                 expression(paste("f(x"["1"], ",", beta["C"], ")")),
                 expression(paste("E"["0"]))
)

text3d(texts$x, texts$y, texts$z, texts=texts$value, cex=1.5)


#### Contrasts

#contrasts <- rbind(points, points[1,])
contrasts <- rbind(points, points[1,], points[3,], points[4,], points[2,])

plot3d(contrasts$x, contrasts$y, contrasts$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)


# Alternative example (use code from before but exclude contrast bit - use this instead)

contrasts <- rbind(points[1,], points[4,])
plot3d(contrasts$x, contrasts$y, contrasts$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

contrasts <- rbind(points[2,], points[3,])
plot3d(contrasts$x, contrasts$y, contrasts$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)



######## ABPlacebo trial ###########

#data <- rbind(A,B)

#plot3d(data$x, data$y, data$z, col=as.numeric(data$agent), type="l", surf=NULL,
#       box=TRUE, axes=FALSE,
#       plot=TRUE)

plot3d(A$x, A$y, A$z, col="blue", type="l", surf=NULL,
       box=FALSE, axes=TRUE, xlab="", ylab="Response", zlab="Dose",
       plot=TRUE)

plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE)


points <- rbind(arm1, arm2, arm4)

plot3d(points$x, points$y, points$z, type="s", size=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)


contrasts <- rbind(points, points[1,])

plot3d(contrasts$x, contrasts$y, contrasts$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

texts <- points
texts$y <- texts$y+15
texts$value <- c(expression(paste("f(x"["2"], ",", beta["A"], ")")),
                 expression(paste("f(x"["1"], ",", beta["B"], ")")),
                 expression(paste("E"["0"]))
)

text3d(texts$x, texts$y, texts$z, texts=texts$value, cex=1.5)





################################################
########## FIGURE 2 - DOSE-RESPONSE CONSISTENCY ############
################################################

E0 <- 15

A <- data.frame("x"=0, "y"=0, "z"=0, "agent"=1)
B <- data.frame("x"=0, "y"=0, "z"=0, "agent"=2)
C <- data.frame("x"=0, "y"=0, "z"=0, "agent"=3)
row <- vector()

for (i in 1:10) {
  row <- c(i, exp(0.5*i), 1*i, 1)
  A <- rbind(A, row)

  row <- c(i, exp(0.4*i), 2*i, 2)
  B <- rbind(B, row)

  row <- c(i, exp(0.45*i), 1.5*i, 3)
  C <- rbind(C, row)
}

A$y <- A$y+E0
B$y <- B$y+E0


#### Points

# Black
arm1 <- data.frame("x"=8,
                   "y"=E0+exp(0.5*8),
                   "z"=1*8)

arm5 <- data.frame("x"=6,
                   "y"=E0+exp(0.5*6),
                   "z"=1*6)

# Red
arm2 <- data.frame("x"=8,
                   "y"=E0+exp(0.4*8),
                   "z"=2*8)

# Green
arm3 <- data.frame("x"=9,
                   "y"=E0+exp(0.45*9),
                   "z"=1.5*9)

# Placebo
arm4 <- data.frame("x"=0,
                   "y"=E0+0,
                   "z"=0)



#### Plot ####

plot3d(A$x, A$y, A$z, col="blue", type="l", surf=NULL,
       box=FALSE, axes=FALSE, xlab="", ylab="", zlab="",
       plot=TRUE)

plot3d(0, 0, 0, type="p", size=0.1, surf=NULL, # To ensure axes plot properly
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE)

points <- rbind(arm1, arm2, arm4)

plot3d(points$x, points$y, points$z, type="p", size=15, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

contrasts <- rbind(points, points[1,])
contrasts1 <- contrasts[-1,]

plot3d(contrasts1$x, contrasts1$y, contrasts1$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

# contrasts2 <- contrasts[1:2,]
# plot3d(contrasts2$x, contrasts2$y, contrasts2$z, type="l", lwd=3, lty=2, surf=NULL,
#        box=FALSE, axes=FALSE,
#        plot=TRUE, add=TRUE)

axis3d('y', pos = c(0, NA, 0))
axis3d('y', pos = c(10, NA, 20), labels=FALSE, tick = FALSE)
axis3d('y', pos = c(10, NA, 0), labels=FALSE, tick = FALSE)

axis3d('z', pos = c(0, 0, NA))
axis3d('z', pos = c(10, 150, NA), labels=FALSE, tick = FALSE)
axis3d('z', pos = c(10, 0, NA), labels=FALSE, tick = FALSE)
#axis3d('z', pos = c(0, 100, NA), labels=FALSE, tick = FALSE)

axis3d('x', pos = c(NA, 0, 0), labels=FALSE, tick=FALSE)
axis3d('x', pos = c(NA, 0, 20), labels=FALSE, tick = FALSE)
axis3d('x', pos = c(NA, 150, 0), labels=FALSE, tick = FALSE)






################################################
########## FIGURE 3 - DOSE-RESPONSE CONSISTENCY ############
################################################

E0 <- 15

A <- data.frame("x"=0, "y"=0, "z"=0, "agent"=1)
B <- data.frame("x"=0, "y"=0, "z"=0, "agent"=2)
C <- data.frame("x"=0, "y"=0, "z"=0, "agent"=3)
row <- vector()

for (i in 1:10) {
  row <- c(i, exp(0.5*i), 1*i, 1)
  A <- rbind(A, row)

  row <- c(i, exp(0.4*i), 2*i, 2)
  B <- rbind(B, row)

  row <- c(i, exp(0.45*i), 1.5*i, 3)
  C <- rbind(C, row)
}

A$y <- A$y+E0
B$y <- B$y+E0


#### Points

# Black
arm1 <- data.frame("x"=8,
                   "y"=E0+exp(0.5*8),
                   "z"=1*8)

arm2 <- data.frame("x"=6,
                   "y"=E0+exp(0.5*6),
                   "z"=1*6)

arm3 <- data.frame("x"=4,
                   "y"=E0+exp(0.5*4),
                   "z"=1*4)

# Red
arm4 <- data.frame("x"=8,
                   "y"=E0+exp(0.4*8),
                   "z"=2*8)

arm5 <- data.frame("x"=7,
                   "y"=E0+exp(0.4*7),
                   "z"=2*7)

arm6 <- data.frame("x"=5,
                   "y"=E0+exp(0.4*5),
                   "z"=2*5)


# Placebo
arm0 <- data.frame("x"=0,
                   "y"=E0+0,
                   "z"=0)



#### Plot ####

plot3d(A$x, A$y, A$z, col="blue", type="l", surf=NULL, lwd=2,
       box=FALSE, axes=FALSE, xlab="", ylab="", zlab="",
       plot=TRUE)

plot3d(0, 0, 0, type="p", size=0.1, surf=NULL, # To ensure axes plot properly
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE, lwd=2)

points <- rbind(arm0, arm1, arm3, arm4, arm6)

plot3d(points$x, points$y, points$z, type="p", size=15, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

points1 <- points[1:3,]
contrasts1 <- points[0,]
for (i in seq_along(points1[,1])) {
  for (k in seq_along(points1[,1])) {
    contrasts1 <- rbind(contrasts1, points1[i,])
    contrasts1 <- rbind(contrasts1, points1[k,])
  }
}

plot3d(contrasts1$x, contrasts1$y, contrasts1$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)



points2 <- points[c(1,4,5),]
contrasts2 <- points[0,]
for (i in seq_along(points2[,1])) {
  for (k in seq_along(points2[,1])) {
    contrasts2 <- rbind(contrasts2, points2[i,])
    contrasts2 <- rbind(contrasts2, points2[k,])
  }
}

plot3d(contrasts2$x, contrasts2$y, contrasts2$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)


#
# contrasts3 <- points[c(2,4),]
# plot3d(contrasts3$x, contrasts3$y, contrasts3$z, type="l", lwd=3, surf=NULL, col="darkgrey",
#        box=FALSE, axes=FALSE,
#        plot=TRUE, add=TRUE)


axis3d('y', pos = c(0, NA, 0))
axis3d('y', pos = c(10, NA, 20), labels=FALSE, tick = FALSE)
axis3d('y', pos = c(10, NA, 0), labels=FALSE, tick = FALSE)

axis3d('z', pos = c(0, 0, NA))
axis3d('z', pos = c(10, 150, NA), labels=FALSE, tick = FALSE)
axis3d('z', pos = c(10, 0, NA), labels=FALSE, tick = FALSE)
#axis3d('z', pos = c(0, 100, NA), labels=FALSE, tick = FALSE)

axis3d('x', pos = c(NA, 0, 0), labels=FALSE, tick=FALSE)
axis3d('x', pos = c(NA, 0, 20), labels=FALSE, tick = FALSE)
axis3d('x', pos = c(NA, 150, 0), labels=FALSE, tick = FALSE)








#### Plot (version 2) ####

plot3d(A$x, A$y, A$z, col="blue", type="l", surf=NULL, lwd=2,
       box=FALSE, axes=FALSE, xlab="", ylab="", zlab="",
       plot=TRUE)

plot3d(0, 0, 0, type="p", size=0.1, surf=NULL, # To ensure axes plot properly
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE, lwd=2)

points <- rbind(arm0, arm1, arm3, arm4, arm6)

plot3d(points$x, points$y, points$z, type="p", size=15, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

points1 <- points[1:3,]
contrasts1 <- points[0,]
for (i in seq_along(points1[,1])) {
  for (k in seq_along(points1[,1])) {
    contrasts1 <- rbind(contrasts1, points1[i,])
    contrasts1 <- rbind(contrasts1, points1[k,])
  }
}

plot3d(contrasts1$x, contrasts1$y, contrasts1$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)



contrasts2 <- points[c(1,5),]
# contrasts2 <- points[0,]
# for (i in seq_along(points2[,1])) {
#   for (k in seq_along(points2[,1])) {
#     contrasts2 <- rbind(contrasts2, points2[i,])
#     contrasts2 <- rbind(contrasts2, points2[k,])
#   }
# }

plot3d(contrasts2$x, contrasts2$y, contrasts2$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

plot3d(points$x[c(2,4)], points$y[c(2,4)], points$z[c(2,4)], type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

axis3d('y', pos = c(0, NA, 0))
axis3d('y', pos = c(10, NA, 20), labels=FALSE, tick = FALSE)
axis3d('y', pos = c(10, NA, 0), labels=FALSE, tick = FALSE)

axis3d('z', pos = c(0, 0, NA))
axis3d('z', pos = c(10, 150, NA), labels=FALSE, tick = FALSE)
axis3d('z', pos = c(10, 0, NA), labels=FALSE, tick = FALSE)

axis3d('x', pos = c(NA, 0, 0), labels=FALSE, tick=FALSE)
axis3d('x', pos = c(NA, 0, 20), labels=FALSE, tick = FALSE)
axis3d('x', pos = c(NA, 150, 0), labels=FALSE, tick = FALSE)





#### Plot (version 3) ####

plot3d(A$x, A$y, A$z, col="blue", type="l", surf=NULL, lwd=2,
       box=FALSE, axes=FALSE, xlab="", ylab="", zlab="",
       plot=TRUE)

plot3d(0, 0, 0, type="p", size=0.1, surf=NULL, # To ensure axes plot properly
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

plot3d(B$x, B$y, B$z, col="red", type="l", add=TRUE, plot=TRUE, lwd=2)

points <- rbind(arm0, # placebo
                arm1, # blue 1
                #arm3, # blue 2
                arm4, # red 1
                arm6 # red 2
                )

plot3d(points$x, points$y, points$z, type="p", size=15, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

contrasts1 <- points[2:3,]
# points1 <- points[c(2,3,4),]
# contrasts1 <- points[0,]
# for (i in seq_along(points1[,1])) {
#   for (k in seq_along(points1[,1])) {
#     contrasts1 <- rbind(contrasts1, points1[i,])
#     contrasts1 <- rbind(contrasts1, points1[k,])
#   }
# }

plot3d(contrasts1$x, contrasts1$y, contrasts1$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)



contrasts2 <- points[c(1,4),]
# contrasts2 <- points[0,]
# for (i in seq_along(points2[,1])) {
#   for (k in seq_along(points2[,1])) {
#     contrasts2 <- rbind(contrasts2, points2[i,])
#     contrasts2 <- rbind(contrasts2, points2[k,])
#   }
# }

plot3d(contrasts2$x, contrasts2$y, contrasts2$z, type="l", lwd=3, surf=NULL,
       box=FALSE, axes=FALSE,
       plot=TRUE, add=TRUE)

# plot3d(points$x[c(2,4)], points$y[c(2,4)], points$z[c(2,4)], type="l", lwd=3, surf=NULL,
#        box=FALSE, axes=FALSE,
#        plot=TRUE, add=TRUE)

axis3d('y', pos = c(0, NA, 0))
axis3d('y', pos = c(10, NA, 20), labels=FALSE, tick = FALSE)
axis3d('y', pos = c(10, NA, 0), labels=FALSE, tick = FALSE)

axis3d('z', pos = c(0, 0, NA))
axis3d('z', pos = c(10, 150, NA), labels=FALSE, tick = FALSE)
axis3d('z', pos = c(10, 0, NA), labels=FALSE, tick = FALSE)

axis3d('x', pos = c(NA, 0, 0), labels=FALSE, tick=FALSE)
axis3d('x', pos = c(NA, 0, 20), labels=FALSE, tick = FALSE)
axis3d('x', pos = c(NA, 150, 0), labels=FALSE, tick = FALSE)
