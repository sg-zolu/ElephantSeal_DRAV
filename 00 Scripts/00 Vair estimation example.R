#Packages required for the analysis
library(readr)
library(MetricsWeighted)
library(Metrics)

# This is a script for estimating the DRAV from a 
# 24-hour recording from individual ml17_301, Day 2-3

# set wd
setwd("01 Data/Example_data")
# Import the dive data from Day2_3
glideList <- read.csv("I1_day2_diveOut.csv")

# Set the parameters for estimation.
# Body density values (pt) are obtained from:
# "99 Data/All Individuals/VODS_dcast_dtag.csv"

# Drag coefficient was assumed to be 0.03 based on the calculations in
# Adachi et al., 2022.
Cdf <- 0.04
# Girth: measured
Girth <- 1.85
# Frontal Area: Calculated from girth
Af <- (Girth/2/pi)^2 * pi
# Mass: measured
m <- 417
# Drag coefficient from the Hydrodynamic Performance Model
CAm <- Cdf*(Af/m)
# Check the value
CAm
# Body density: calculated
pt <- 1041.220
# gravitational acceleration
g <- 9.8
# Air density (kg/m^3)
Pair <- 1.23

# Go to "01 Vair estimation calculation" to calculate the DRAV of each dive

#After analysis
Day2_3_0.03 <- data.frame(Vair.glide)
#View(Day2_3_0.03)
setwd("01 Data")
write.csv(Day2_3_0.03, file = "01 Data/Day2_3_Vair_Cdf0.03.csv")
