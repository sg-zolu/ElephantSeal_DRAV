library(MetricsWeighted)

################################################################################
#                               Preparing data                                 #
################################################################################
# Here, the data will be separated according to the following criteria:
# 1: Dives that started gliding before 30m depth
# 2: Dives with a glide duration longer than 10 seconds

# First, remove the surfacing periods from the data as we are 
# only interested in the descent/ascent phase of the data
glideList <- subset(glideList, !is.na(glideList$dive_no))
# Factorise the dive number
glideList$dive_no <- as.factor(glideList$dive_no)
# Extract descent and ascent phases of the dive between the 
# depths of 10 to 50m
# pitch < 0 is for extracting descent phase, as the animal will be descending
# facing downwards
glideList <- subset(glideList, d >= 10 & d <= 50 & pitch < 0)
# Fill in the rms NA values with 0
glideList$heave_rms[is.na(glideList$heave_rms)] <- 0
glideList$heave_rms[is.na(glideList$sway_rms)] <- 0
glideList$heave_rms[is.na(glideList$surge_rms)] <- 0

# Now, we have a data set containing 10 to 50m of the descent phase
# from each dive.
# View(glideList)

# We will conduct an analysis whether the descent/ascent phase for each dive
# meet the criteria

# Extract the levels of the dives
dives <-unique(glideList$dive_no)
m <- length(dives)

# Create an empty list to store the extracted descents that meet the criteria
dataframe.list<-list()

# For loop to check the criteria for each dive
for(i in 1:m){
  # Subset the data into the dataframe.list
  dataframe.list[[i]] <- subset(glideList, dive_no==dives[i])
  glide <- dataframe.list[[i]] # Create a dataframe called 'glide' for analysis
  
  # Only if there are more than 10 seconds for the depth measurements
  # Given that depth was sampled at 1Hz (10 data points = 10 seconds)
  if(length(glide$d) > 9){
    # Here, we are checking the following:
    # 1. Initial gliding depth
    # 2. Whether the glide is longer than 10 seconds

    # Create a new data frame 'sum.data' that sums up the heave rms from the 
    # descent phase of the dive, which will reveal whether the seal is gliding
    # or stroking (sum = 0, meaning that the seal is gliding for that duration)
    sum.data <- data.frame(sum = (numeric(length(glide$d)-9)))
    
    # Sum the consecutive 10 data points from each point in the data
    for(p in 1:(length(sum.data$sum))){
      glide.d <- glide$d[p]
      glide.sum <- glide$heave_rms[p] + glide$heave_rms[p+1] + glide$heave_rms[p+2] + glide$heave_rms[p+3] + glide$heave_rms[p+4] + glide$heave_rms[p+5] + glide$heave_rms[p+6] + glide$heave_rms[p+7] + glide$heave_rms[p+8] + glide$heave_rms[p+9]
      sum.data$sum[p] <- glide.sum
      sum.data$d[p] <- glide.d
      
      # If the sum is 0, that means the animal is gliding for at least 10 secs
      if(sum.data$sum[p] == 0){
        # Finding the initial gliding depth
        Zeros <- sum.data[sum.data$sum == 0,]
        initial_depth <- Zeros$d[1]
        glide$initial_d <- initial_depth
      }else{
        # If the sum is not 0: glide depth is 0
        glide$initial_d <- 0
      }
    }
    
    # Filtering out dives with a glide depth greater than 30m
    if(glide$initial_d[1] > 0 & glide$initial_d[1] <= 30){
      
      # Extract data with speed measurements
      glide <- subset(glide, !is.na(glide$mvelocity))
      # Extract data where the animal is only gliding
      glide <- subset(glide, d >= glide$initial_d[1])
      # Add time from glide
      glide$t <- c(0:(length(glide$d)-1))
      # Remove columns that are unnecessary for the analysis
      glide <- glide[,c(1:4,10:14)]
      # Convert pitch into deg
      glide$pitch <- (glide$pitch*pi)/180
      
      dataframe.list[[i]] <- glide
      
      # Write a csv file for each dive
      write.csv(dataframe.list[[i]], file = 
                  paste0("dive_",i, 
                         ".csv"), row.names = FALSE)
      
    }
  }
}

################################################################################
#                               DRAV estimation                                #
################################################################################
# This section will estimate the DRAV for each dive

# Make a list to include the dives that met the criteria
fileList <- grep("dive_", list.files(getwd()), value = TRUE)
# Make a data frame that contains different values of DRAV in litres
Vair.seq <- data.frame(Vair = seq(0,0.02,by=0.0001))
# Make a data frame to contain the 'selected' DRAV
Vair.glide <- data.frame(Vair = numeric(length(fileList)))

# A for loop that will go through individual dives
for(f in 1:length(fileList)) {
  glide <- read.csv(fileList[f],sep=",", header=T) # read in the files
  
  # Estimate the velocity with different values of DRAV
  for(j in 1:length(Vair.seq$Vair)) {  
    Vair = Vair.seq$Vair[j]
    
    # Calculate the velocity using the formula
    for(i in 1:length(glide$d)) {
      if(i == 1){
        glide$tvelocity[i] <- glide$mvelocity[i]
        glide$accel[i] <- -CAm*0.5*glide$dsw[i]*glide$tvelocity[i]^2 + (glide$dsw[i]/pt-1)*g*sin(glide$pitch[i])+((Vair/m)*g*sin(glide$pitch[i])*((glide$dsw[i]+Pair*(1+0.1*glide$d[i]))/(1+0.1*glide$d[i])))
        glide$V_w[i] <- 1/(1+0.1*glide$d[i])
      }else{
        glide$tvelocity[i] <- glide$tvelocity[i-1] + glide$accel[i-1]
        glide$accel[i] <- -CAm*0.5*glide$dsw[i]*glide$tvelocity[i]^2 + (glide$dsw[i]/pt-1)*g*sin(glide$pitch[i])+((Vair/m)*g*sin(glide$pitch[i])*((glide$dsw[i]+Pair*(1+0.1*glide$d[i]))/(1+0.1*glide$d[i])))
        glide$V_w[i] <- 1/(1+0.1*glide$d[i])
      }
    }
    
    # Store the data calculated from different DRAV
    
    # Weighted MSE
    Vair.seq$MSE_w[j] <- rmse(glide$mvelocity, glide$tvelocity, w=glide$V_w)
    # R-squared value of the two velocities
    Vair.seq$LMr[j] <- summary(lm(mvelocity ~ tvelocity, data=glide))$r.squared
    # Glide depth
    Vair.seq$depth[j] <- glide$d[1]
    # Time
    Vair.seq$t[j] <- length(glide$t)
    # Pitch in deg
    Vair.seq$pitch[j] <- (mean(glide$pitch))*180/pi
    # Measured velocity
    Vair.seq$speed[j] <- mean(glide$mvelocity)
    # Dive number
    Vair.seq$dive_no[j] <- glide$dive_no[1]
    # The initial velocity
    initial_velocity <- glide$mvelocity[1]
    # The minimum velocity in the dive
    min_velocity <- min(glide$mvelocity)
  }
  
  # Extracting data from the DRAV with the least MSE
  
  # Convert the Vair value into millilitres for analysis
  Vair.glide$Vair[f] <- 1000000*Vair.seq$Vair[which.min(Vair.seq$MSE_w)]/m
  # The calculated MSE
  Vair.glide$minMSE_w[f] <- Vair.seq$MSE_w[which.min(Vair.seq$MSE_w)]
  # The initial gliding depth
  Vair.glide$initial_depth[f] <- Vair.seq$depth[which.min(Vair.seq$MSE_w)]
  # Glide duration
  Vair.glide$glide_duration[f] <- Vair.seq$t[which.min(Vair.seq$MSE_w)]
  # Pitch
  Vair.glide$avg_pitch[f] <- Vair.seq$pitch[which.min(Vair.seq$MSE_w)]
  # Velocity
  Vair.glide$avg_speed[f] <- Vair.seq$speed[which.min(Vair.seq$MSE_w)]
  # R-squared value
  Vair.glide$r_sq[f] <- Vair.seq$LMr[which.min(Vair.seq$MSE_w)]
  # Dive number
  Vair.glide$DiveNo[f] <-Vair.seq$dive_no[which.min(Vair.seq$MSE_w)]
  # Initial velocity
  Vair.glide$initial_v[f] <- initial_velocity
  # Minimum velocity
  Vair.glide$min_v[f] <- min_velocity
  # Depth with minimum velocity
  Vair.glide$d_min_v[f] <-Vair.seq$depth[which.min(Vair.seq$speed)]
}

# Now go back to "00 Vair estimation parameters" and complete the 
# "After analysis" section to save the data


